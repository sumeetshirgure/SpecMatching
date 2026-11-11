`ifndef hc_fifo
`define hc_fifo

module hc_fifo #(
    parameter NUM_SYND        = 38,     // defects in the syndrome (sets the local index width)
    parameter LOG_NUM_DET     = 12,     // detector id width
    parameter NUM_LANES       = 64,     // entries per chunk (sets the lane bound width)
    parameter CHUNK_W         = 1792,   // one chunk, NUM_LANES*LANE_W
    parameter CHUNKS_PER_BALL = 4,      // chunks per ball-table word
    parameter FIFO_DEPTH      = 4       // batches buffered; must be a power of two and >= 2
) (
    input  wire                                                               iClk,
    input  wire                                                               iRstN,
    input  wire                                                               iStart,
    output wire                                                               oGo,
    output wire                                                               oStall,
    input  wire                                                               iLoadBusy,
    input  wire                                                               iChValid,
    input  wire [$clog2(NUM_SYND+1)-1:0]                                      iRdU,
    input  wire [LOG_NUM_DET-1:0]                                             iRdDet,
    input  wire [(($clog2(CHUNKS_PER_BALL)>0)?$clog2(CHUNKS_PER_BALL):1)-1:0] iRdCh,
    input  wire [$clog2(NUM_LANES)-1:0]                                       iRdHi,
    input  wire [CHUNKS_PER_BALL*CHUNK_W-1:0]                                 iBallRdata,
    output wire                                                               oBatchValid,
    output wire [$clog2(NUM_SYND+1)-1:0]                                      oBatchU,
    output wire [LOG_NUM_DET-1:0]                                             oBatchDet,
    output wire [$clog2(NUM_LANES)-1:0]                                       oBatchHi,
    output wire [CHUNK_W-1:0]                                                 oBatchWord,
    input  wire                                                               iLaneReady,
    input  wire                                                               iLaneBusy,
    input  wire                                                               iCollectBusy,
    input  wire                                                               iBndBusy,
    output wire                                                               oDone
);
    localparam LOG_NUM_SYND_P1 = $clog2(NUM_SYND+1);                                            // local defect index width
    localparam LOG_NUM_LANES   = $clog2(NUM_LANES);                                             // lane bound width
    localparam CIDX_W          = ($clog2(CHUNKS_PER_BALL) > 0) ? $clog2(CHUNKS_PER_BALL) : 1;
    localparam SB_W            = LOG_NUM_SYND_P1 + LOG_NUM_DET + LOG_NUM_LANES;                 // sideband {u, det_i, hi}
    localparam SLOT_W          = CHUNK_W + SB_W;                                                // one FIFO slot
    localparam PTR_W           = $clog2(FIFO_DEPTH);                                            // FIFO pointer width
    localparam CNT_W           = $clog2(FIFO_DEPTH+1);                                          // occupancy width

    // two push stages are in flight when oStall rises, so stall two slots early
    localparam [CNT_W-1:0] FIFO_HIGH = FIFO_DEPTH - 2;
    localparam [CNT_W-1:0] CNT_ZERO  = 0;

    localparam [1:0] S_IDLE = 2'd0,
                     S_RUN  = 2'd1,
                     S_DONE = 2'd2;

    reg              rd_v;
    reg [SB_W-1:0]   rd_sb;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) rd_v <= 1'b0;
        else        rd_v <= iChValid;

    always @(posedge iClk)
        rd_sb <= {iRdU, iRdDet, iRdHi};


    // ---------------- chunk cut, select replicated across slices ----------------
    localparam SEL_COPIES = 8;
    localparam SLICE_W    = (CHUNK_W + SEL_COPIES - 1) / SEL_COPIES;

    reg [SEL_COPIES*CIDX_W-1:0] rd_ch_rep;
    integer kc;
    always @(posedge iClk)
        for (kc = 0; kc < SEL_COPIES; kc = kc + 1)
            rd_ch_rep[kc*CIDX_W +: CIDX_W] <= iRdCh;

    wire [CHUNK_W-1:0] chunk;

    genvar gk;
    generate
        for (gk = 0; gk < SEL_COPIES; gk = gk + 1) begin : g_chunk_slice
            localparam integer LO  = gk * SLICE_W;
            localparam integer WID = (LO + SLICE_W <= CHUNK_W) ? SLICE_W : (CHUNK_W - LO);
            if (WID > 0) begin : g_live
                wire [CIDX_W-1:0] sel_k = rd_ch_rep[gk*CIDX_W +: CIDX_W];
                reg  [WID-1:0]    slice;
                reg  [CIDX_W:0]   ci;   // one bit wider than the index so the bound test cannot wrap
                always @* begin
                    slice = iBallRdata[LO +: WID];
                    for (ci = 1; ci < CHUNKS_PER_BALL; ci = ci + 1)
                        if (sel_k == ci[CIDX_W-1:0])
                            slice = iBallRdata[ci[CIDX_W-1:0]*CHUNK_W + LO +: WID];
                end
                assign chunk[LO +: WID] = slice;
            end
        end
    endgenerate

    // ---------------- registered chunk cut ----------------
    reg [CHUNK_W-1:0] chunk_q;
    reg [SB_W-1:0]    rd_sb2;
    reg               rd_v2;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) rd_v2 <= 1'b0;
        else        rd_v2 <= rd_v;

    always @(posedge iClk) begin
        chunk_q <= chunk;
        rd_sb2  <= rd_sb;
    end

    // ---------------- batch FIFO ----------------
    reg [SLOT_W-1:0] mem [0:FIFO_DEPTH-1];
    reg [PTR_W-1:0]  wptr, rptr;
    reg [CNT_W-1:0]  cnt;

    // --------- occupancy flags, registered from cnt_next ---------
    reg  mem_empty;
    reg  stall_q;
    wire push = rd_v2;

    // --------- registered output stage (first-word-fall-through) ---------
    reg              o_v;
    reg [SLOT_W-1:0] o_d;

    wire o_shift = !o_v || iLaneReady;        // the output register is free, or is being consumed
    wire pop     = o_shift && !mem_empty;     // ... and a slot is there to take, so one leaves

    always @(posedge iClk or negedge iRstN)
        if (!iRstN)       o_v <= 1'b0;
        else if (o_shift) o_v <= !mem_empty;

    always @(posedge iClk)
        if (pop) o_d <= mem[rptr];

    wire [CNT_W-1:0] cnt_next = (push && !pop) ? cnt + 1'b1 :
                                (pop && !push) ? cnt - 1'b1 : cnt;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) begin
            wptr      <= {PTR_W{1'b0}};
            rptr      <= {PTR_W{1'b0}};
            cnt       <= CNT_ZERO;
            mem_empty <= 1'b1;
            stall_q   <= 1'b0;
        end else begin
            if (push) wptr <= wptr + 1'b1;
            if (pop)  rptr <= rptr + 1'b1;
            cnt       <= cnt_next;
            mem_empty <= (cnt_next == CNT_ZERO);
            stall_q   <= (cnt_next >= FIFO_HIGH);
        end

    always @(posedge iClk)
        if (push) mem[wptr] <= {chunk_q, rd_sb2};

    assign oBatchValid = o_v;
    assign oBatchWord  = o_d[SB_W +: CHUNK_W];
    assign oBatchU     = o_d[LOG_NUM_DET+LOG_NUM_LANES +: LOG_NUM_SYND_P1];
    assign oBatchDet   = o_d[LOG_NUM_LANES +: LOG_NUM_DET];
    assign oBatchHi    = o_d[0 +: LOG_NUM_LANES];
    assign oStall      = stall_q;

    // ---------------- run / done ----------------
    // a batch held in the output register is not drained yet
    wire empty = mem_empty && !o_v;

    // both push stages count as busy
    wire quiet = !iLoadBusy && !rd_v && !rd_v2 && empty &&
                 !iLaneBusy && !iCollectBusy && !iBndBusy;

    reg [1:0] state;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) state <= S_IDLE;
        else
            case (state)
                S_IDLE, S_DONE: if (iStart) state <= S_RUN;
                S_RUN:          if (quiet)   state <= S_DONE;
                default:        state <= S_IDLE;
            endcase

    assign oGo   = ((state == S_IDLE) || (state == S_DONE)) && iStart;
    assign oDone = (state == S_DONE);
endmodule
`endif
