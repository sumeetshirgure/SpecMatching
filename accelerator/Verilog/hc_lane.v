`ifndef hc_lane
`define hc_lane

module hc_lane #(
    parameter NUM_SYND    = 38,      // syndrome capacity (sorted-array length)
    parameter LOG_NUM_DET = 12,      // detector id width
    parameter WEIGHT_W    = 16       // integer weight width
) (
    input  wire                            iClk,
    input  wire                            iRstN,
    input  wire                            iEn,
    input  wire                            iValid,
    input  wire                            iActive,
    input  wire [WEIGHT_W+LOG_NUM_DET-1:0] iLane,
    input  wire [LOG_NUM_DET-1:0]          iDet,
    input  wire [NUM_SYND*LOG_NUM_DET-1:0] iSynFlat,
    input  wire [NUM_SYND-1:0]             iSynValid,
    output reg                             oHit,
    output reg [$clog2(NUM_SYND+1)-1:0]    oIdx,
    output reg [WEIGHT_W-1:0]              oW
);
    localparam LOG_NUM_SYND_P1 = $clog2(NUM_SYND+1);   // local defect index width
    localparam STAGES          = LOG_NUM_SYND_P1;      // probe levels, one pipeline stage each
    localparam PIPE            = (STAGES > 1) ? STAGES - 1 : 1;             // inter-stage depth
    localparam TRI             = (STAGES > 1) ? STAGES*(STAGES-1)/2 : 1;    // packed path bits

    wire [WEIGHT_W-1:0]    e_w   = iLane[LOG_NUM_DET +: WEIGHT_W];
    wire [LOG_NUM_DET-1:0] e_tgt = iLane[0 +: LOG_NUM_DET];

    wire pass = iValid && iActive && (e_tgt > iDet);

    reg                   s1_v;
    reg [LOG_NUM_DET-1:0] s1_tgt;
    reg [WEIGHT_W-1:0]    s1_w;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN)   s1_v <= 1'b0;
        else if (iEn) s1_v <= pass;

    always @(posedge iClk)
        if (iEn) begin
            s1_tgt <= e_tgt;
            s1_w   <= e_w;
        end

    reg [PIPE-1:0]                 p_v, p_hit;
    reg [PIPE*LOG_NUM_DET-1:0]     p_key;
    reg [PIPE*WEIGHT_W-1:0]        p_w;
    reg [PIPE*LOG_NUM_SYND_P1-1:0] p_hidx;
    reg [TRI-1:0]                  p_path;   // level s holds s+1 bits at offset s*(s+1)/2

    genvar gs, gc;
    generate
        for (gs = 0; gs < STAGES; gs = gs + 1) begin : g_bs
            // ---- this level's fixed candidate slots: 2**gs of them ----
            localparam NCAND = 1 << gs;
            wire [NCAND*LOG_NUM_DET-1:0] cand_id;
            wire [NCAND-1:0]             cand_occ;
            for (gc = 0; gc < NCAND; gc = gc + 1) begin : g_cand
                localparam integer NODE =
                    gc * (1 << (STAGES - gs)) + (1 << (STAGES - 1 - gs)) - 1;
                if (NODE < NUM_SYND) begin : g_in
                    assign cand_id[gc*LOG_NUM_DET +: LOG_NUM_DET] =
                        iSynFlat[NODE*LOG_NUM_DET +: LOG_NUM_DET];
                    assign cand_occ[gc] = iSynValid[NODE];
                end else begin : g_pad
                    // past the array: +infinity, never occupied
                    assign cand_id[gc*LOG_NUM_DET +: LOG_NUM_DET] = {LOG_NUM_DET{1'b0}};
                    assign cand_occ[gc] = 1'b0;
                end
            end

            // ---- stage inputs: filter stage for gs = 0, the previous level after ----
            wire                       in_v, in_hit;
            wire [LOG_NUM_DET-1:0]     in_key;
            wire [WEIGHT_W-1:0]        in_w;
            wire [LOG_NUM_SYND_P1-1:0] in_hidx;
            wire [((gs > 0) ? gs : 1)-1:0] lpath;

            if (gs == 0) begin : g_head
                assign in_v    = s1_v;
                assign in_key  = s1_tgt;
                assign in_w    = s1_w;
                assign in_hit  = 1'b0;
                assign in_hidx = {LOG_NUM_SYND_P1{1'b0}};
                assign lpath   = 1'b0;
            end else begin : g_tail
                assign in_v    = p_v[gs-1];
                assign in_key  = p_key[(gs-1)*LOG_NUM_DET +: LOG_NUM_DET];
                assign in_w    = p_w[(gs-1)*WEIGHT_W +: WEIGHT_W];
                assign in_hit  = p_hit[gs-1];
                assign in_hidx = p_hidx[(gs-1)*LOG_NUM_SYND_P1 +: LOG_NUM_SYND_P1];
                assign lpath   = p_path[(gs-1)*gs/2 +: gs];
            end

            // ---- the probed node's array index, pure wiring off the path ----
            wire [LOG_NUM_SYND_P1-1:0] node_idx;
            if (STAGES == 1) begin : g_nidx1
                assign node_idx = 1'b0;
            end else if (gs == 0) begin : g_nidx0
                assign node_idx = {1'b0, {(STAGES-1){1'b1}}};
            end else if (gs == STAGES-1) begin : g_nidxl
                assign node_idx = {lpath, 1'b0};
            end else begin : g_nidxm
                assign node_idx = {lpath, 1'b0, {(STAGES-1-gs){1'b1}}};
            end

            // ---- the probe: one compare against the path-selected slot ----
            wire [LOG_NUM_DET-1:0] pid = cand_id[lpath*LOG_NUM_DET +: LOG_NUM_DET];
            wire                   occ = cand_occ[lpath];
            wire eq = occ && (pid == in_key);

            wire                       take   = eq && !in_hit;   // first hit freezes the index
            wire                       n_hit  = in_hit || eq;
            wire [LOG_NUM_SYND_P1-1:0] n_hidx = take ? node_idx : in_hidx;

            if (gs < STAGES-1) begin : g_mid
                // the descend direction only exists where a next level consumes it
                wire go_right = occ && (pid < in_key);
                wire [gs:0] n_path;
                if (gs == 0) begin : g_np0
                    assign n_path = go_right;
                end else begin : g_npn
                    assign n_path = {lpath, go_right};
                end

                always @(posedge iClk or negedge iRstN)
                    if (!iRstN)   p_v[gs] <= 1'b0;
                    else if (iEn) p_v[gs] <= in_v;

                always @(posedge iClk)
                    if (iEn) begin
                        p_key[gs*LOG_NUM_DET +: LOG_NUM_DET]          <= in_key;
                        p_w[gs*WEIGHT_W +: WEIGHT_W]                  <= in_w;
                        p_hit[gs]                                     <= n_hit;
                        p_hidx[gs*LOG_NUM_SYND_P1 +: LOG_NUM_SYND_P1] <= n_hidx;
                        p_path[gs*(gs+1)/2 +: gs+1]                   <= n_path;
                    end
            end else begin : g_last
                always @(posedge iClk or negedge iRstN)
                    if (!iRstN)   oHit <= 1'b0;
                    else if (iEn) oHit <= in_v && n_hit;

                always @(posedge iClk)
                    if (iEn) begin
                        oIdx <= n_hidx;
                        oW   <= in_w;
                    end
            end
        end
    endgenerate
endmodule
`endif
