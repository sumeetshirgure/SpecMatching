`ifndef h_construct
`define h_construct

module h_construct #(
    parameter NUM_SYND        = 38,      // syndrome capacity, max defects per shot (max H nodes)
    parameter LOG_NUM_DET     = 12,      // detector id width = address width of both SRAMs
    parameter LOG_MAX_BALL    = 9,       // raw ball length width
    parameter WEIGHT_W        = 16,      // integer weight width
    parameter NUM_LANES       = 64,      // ball entries searched per cycle
    parameter LANE_W          = 28,      // one lane field, WEIGHT_W + LOG_NUM_DET
    parameter CHUNK_W         = 1792,    // one chunk of a ball-table word, NUM_LANES*LANE_W
    parameter CHUNKS_PER_BALL = 4        // chunks per ball-table word (one word per node)
) (
    input  wire                               iClk,
    input  wire                               iRstN,
    input  wire                               iStart,
    input  wire [$clog2(NUM_SYND+1)-1:0]      iM,
    input  wire [NUM_SYND*LOG_NUM_DET-1:0]    iSynFlat,
    output wire                               oBndTabRdEn,
    output wire [LOG_NUM_DET-1:0]             oBndTabAddr,
    input  wire [WEIGHT_W+LOG_MAX_BALL:0]     iBndTabRdata,
    output wire                               oBallRdEn,
    output wire [LOG_NUM_DET-1:0]             oBallAddr,
    input  wire [CHUNKS_PER_BALL*CHUNK_W-1:0] iBallRdata,
    output wire                               oEdgeValid,
    output wire [$clog2(NUM_SYND+1)-1:0]      oEdgeU,
    output wire [$clog2(NUM_SYND+1)-1:0]      oEdgeV,
    output wire [WEIGHT_W-1:0]                oEdgeW,
    output wire                               oBndValid,
    output wire [$clog2(NUM_SYND+1)-1:0]      oBndU,
    output wire [WEIGHT_W-1:0]                oBndW,
    output wire                               oDone
);
    localparam LOG_NUM_SYND_P1 = $clog2(NUM_SYND+1);                                            // local defect index width
    localparam LOG_NUM_LANES   = $clog2(NUM_LANES);                                             // lane range bound width
    localparam CIDX_W          = ($clog2(CHUNKS_PER_BALL) > 0) ? $clog2(CHUNKS_PER_BALL) : 1;
    localparam LANE_LATENCY    = 1 + $clog2(NUM_SYND+1);

    wire                   bnd_hb = iBndTabRdata[WEIGHT_W+LOG_MAX_BALL];
    wire [WEIGHT_W-1:0]    bnd_bw = iBndTabRdata[LOG_MAX_BALL +: WEIGHT_W];
    wire [LOG_MAX_BALL-1:0] bnd_size = iBndTabRdata[0 +: LOG_MAX_BALL];

    // ---------------- address generator ----------------
    wire                    go, stall, load_busy, ch_valid;
    wire [LOG_NUM_SYND_P1-1:0] rd_u;
    wire [LOG_NUM_DET-1:0]     rd_det;
    wire [CIDX_W-1:0]       rd_ch;
    wire [LOG_NUM_LANES-1:0]   rd_hi;
    wire                    hdr_valid;
    wire [LOG_NUM_SYND_P1-1:0] hdr_u;

    hc_addr_gen #(
        .NUM_SYND(NUM_SYND), .LOG_NUM_DET(LOG_NUM_DET), .LOG_MAX_BALL(LOG_MAX_BALL),
        .NUM_LANES(NUM_LANES), .CHUNKS_PER_BALL(CHUNKS_PER_BALL)
    ) u_addr_gen (
        .iClk(iClk), .iRstN(iRstN), .iGo(go), .iStall(stall),
        .iM(iM), .iSynFlat(iSynFlat),
        .oBndTabRdEn(oBndTabRdEn), .oBndTabAddr(oBndTabAddr),
        .iBndSize(bnd_size),
        .oBallRdEn(oBallRdEn), .oBallAddr(oBallAddr),
        .oChValid(ch_valid),
        .oRdU(rd_u), .oRdDet(rd_det), .oRdCh(rd_ch), .oRdHi(rd_hi),
        .oHdrValid(hdr_valid), .oHdrU(hdr_u),
        .oBusy(load_busy)
    );

    // ---------------- boundary check ----------------
    wire bnd_busy;

    hc_bnd_check #(
        .NUM_SYND(NUM_SYND), .WEIGHT_W(WEIGHT_W)
    ) u_bnd_check (
        .iClk(iClk), .iRstN(iRstN),
        .iHdrValid(hdr_valid), .iHdrU(hdr_u),
        .iHb(bnd_hb), .iBw(bnd_bw),
        .oBndValid(oBndValid), .oBndU(oBndU), .oBndW(oBndW),
        .oBusy(bnd_busy)
    );

    // ---------------- control + batch FIFO ----------------
    wire                     batch_valid;
    wire [LOG_NUM_SYND_P1-1:0]  batch_u;
    wire [LOG_NUM_DET-1:0]      batch_det;
    wire [LOG_NUM_LANES-1:0]    batch_hi;
    wire [CHUNK_W-1:0]       batch_word;
    wire                     lane_ready, lane_busy, collect_busy;

    hc_fifo #(
        .NUM_SYND(NUM_SYND), .LOG_NUM_DET(LOG_NUM_DET), .NUM_LANES(NUM_LANES),
        .CHUNK_W(CHUNK_W), .CHUNKS_PER_BALL(CHUNKS_PER_BALL)
    ) u_ctl (
        .iClk(iClk), .iRstN(iRstN), .iStart(iStart),
        .oGo(go), .oStall(stall), .iLoadBusy(load_busy),
        .iChValid(ch_valid), .iRdU(rd_u), .iRdDet(rd_det),
        .iRdCh(rd_ch), .iRdHi(rd_hi),
        .iBallRdata(iBallRdata),
        .oBatchValid(batch_valid), .oBatchU(batch_u), .oBatchDet(batch_det),
        .oBatchHi(batch_hi),
        .oBatchWord(batch_word), .iLaneReady(lane_ready),
        .iLaneBusy(lane_busy), .iCollectBusy(collect_busy), .iBndBusy(bnd_busy),
        .oDone(oDone)
    );

   wire [NUM_SYND-1:0] syn_valid;

    genvar gm;
    generate
        for (gm = 0; gm < NUM_SYND; gm = gm + 1) begin : g_mask
            localparam [LOG_NUM_SYND_P1-1:0] SLOT_NO = gm;
            assign syn_valid[gm] = (SLOT_NO < iM);
        end
    endgenerate

    wire [NUM_LANES-1:0]                hit_flat;
    wire [NUM_LANES*LOG_NUM_SYND_P1-1:0]   idx_flat;
    wire [NUM_LANES*WEIGHT_W-1:0]       w_flat;

    genvar g;
    generate
        for (g = 0; g < NUM_LANES; g = g + 1) begin : g_lane
            localparam [LOG_NUM_LANES-1:0] LANE_NO = g;
            hc_lane #(
                .NUM_SYND(NUM_SYND),
                .LOG_NUM_DET(LOG_NUM_DET), .WEIGHT_W(WEIGHT_W)
            ) u_lane (
                .iClk(iClk), .iRstN(iRstN), .iEn(lane_ready),
                .iValid(batch_valid),
                .iActive(LANE_NO <= batch_hi),
                .iLane(batch_word[g*LANE_W +: LANE_W]), .iDet(batch_det),
                .iSynFlat(iSynFlat), .iSynValid(syn_valid),
                .oHit(hit_flat[g]),
                .oIdx(idx_flat[g*LOG_NUM_SYND_P1 +: LOG_NUM_SYND_P1]),
                .oW(w_flat[g*WEIGHT_W +: WEIGHT_W])
            );
        end
    endgenerate

    reg [LANE_LATENCY-1:0]                 v_pipe;
    reg [LANE_LATENCY*LOG_NUM_SYND_P1-1:0] u_pipe;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN)        v_pipe <= {LANE_LATENCY{1'b0}};
        else if (lane_ready) v_pipe <= {v_pipe[LANE_LATENCY-2:0], batch_valid};

    always @(posedge iClk)
        if (lane_ready)
            u_pipe <= {u_pipe[(LANE_LATENCY-1)*LOG_NUM_SYND_P1-1:0], batch_u};

    assign lane_busy = |v_pipe;

    // ---------------- output collect ----------------
    hc_collect #(
        .NUM_SYND(NUM_SYND), .WEIGHT_W(WEIGHT_W), .NUM_LANES(NUM_LANES)
    ) u_collect (
        .iClk(iClk), .iRstN(iRstN),
        .iValid(v_pipe[LANE_LATENCY-1]),
        .iU(u_pipe[(LANE_LATENCY-1)*LOG_NUM_SYND_P1 +: LOG_NUM_SYND_P1]),
        .iHitFlat(hit_flat), .iIdxFlat(idx_flat), .iWFlat(w_flat),
        .oReady(lane_ready), .oBusy(collect_busy),
        .oEdgeValid(oEdgeValid), .oEdgeU(oEdgeU),
        .oEdgeV(oEdgeV), .oEdgeW(oEdgeW)
    );
endmodule
`endif
