`ifndef hc_lane
`define hc_lane

module hc_lane #(
    parameter NUM_SYND    = 38,      // syndrome capacity (sorted-array length)
    parameter LOG_NUM_DET = 12,      // detector id width
    parameter DELTA_W     = 12,      // signed target offset width, <= LOG_NUM_DET
    parameter WEIGHT_W    = 16       // integer weight width
) (
    input  wire                            iClk,
    input  wire                            iRstN,
    input  wire                            iEn,
    input  wire                            iValid,
    input  wire                            iActive,
    input  wire [WEIGHT_W+DELTA_W-1:0]     iLane,
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
    localparam MUX_T           = (STAGES >= 9) ? 3 : 4;                    // widest mux is 2**MUX_T : 1
    localparam PAIR_GS         = STAGES - 2;                               // last level with a pair mux
    localparam LOOK            = (PAIR_GS > MUX_T) ? 1 : 0;                // some level needs the mux cap
    localparam GO_W            = LOOK ? PAIR_GS - MUX_T : 1;               // capped levels
    localparam REP_TOP         = (LOOK ? MUX_T : PAIR_GS) - 1;             // top replicated path bank
    localparam TRI_B           = (REP_TOP >= 0) ? (REP_TOP+1)*(REP_TOP+2)/2 : 1;  // replicated path bits
    localparam CAND_W          = LOG_NUM_DET + 1;      // one candidate slot: {occ, id}

    wire [WEIGHT_W-1:0]    e_w   = iLane[DELTA_W +: WEIGHT_W];
    wire [DELTA_W-1:0]     e_dlt = iLane[0 +: DELTA_W];
    wire [LOG_NUM_DET-1:0] e_dlt_x;

    generate
        if (DELTA_W < LOG_NUM_DET) begin : g_sext
            assign e_dlt_x = {{(LOG_NUM_DET-DELTA_W){e_dlt[DELTA_W-1]}}, e_dlt};
        end else begin : g_nosext
            assign e_dlt_x = e_dlt;
        end
    endgenerate

    wire [LOG_NUM_DET-1:0] e_tgt = iDet + e_dlt_x;

    wire pass = iValid && iActive && !e_dlt[DELTA_W-1];

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

    // Each level registers the candidate PAIR the next level will probe, in parallel with its own
    // compare, so a level only needs a 2:1 pick on go_right in front of the compare.
    reg [PIPE-1:0]                 p_v, p_hit;
    reg [PIPE*LOG_NUM_DET-1:0]     p_key;
    reg [PIPE*WEIGHT_W-1:0]        p_w;
    reg [PIPE*LOG_NUM_SYND_P1-1:0] p_hidx;
    reg [TRI-1:0]                  p_path;   // level s holds s+1 bits at offset s*(s+1)/2
    reg [TRI_B-1:0]                p_path_b; // copy of p_path for the pair mux
    reg [TRI_B-1:0]                p_path_c; // the same, stored inverted
    reg [GO_W-1:0]                 p_go_b;   // go_right copies for levels past MUX_T
    reg [GO_W-1:0]                 p_go_c;   // the same, stored inverted
    reg [(1<<MUX_T)-1:0]           p_hot_b;  // one-hot select of g_look's first stage
    reg [(1<<MUX_T)-1:0]           p_hot_c;  // the same, stored inverted
    wire [GO_W*2*CAND_W-1:0]       look_pair; // the pair g_look selected for banks past MUX_T
    reg [PIPE*2*CAND_W-1:0]        p_pair;   // bank s holds level s+1's candidate pair

    genvar gs, gc;
    // ---- mux cap: levels past MUX_T narrow their pair over several stages, no added latency ----
    generate
        if (LOOK) begin : g_look
            genvar tb, ti;
            for (tb = MUX_T + 1; tb <= PAIR_GS; tb = tb + 1) begin : g_tb
                localparam L     = tb - MUX_T;
                localparam SLOTS = 1 << (tb + 1);
                // stage ti holds 2**(L+1-ti) candidates, packed one stage after another
                reg [((1 << (L + 2)) - 4)*CAND_W-1:0] q;
                wire [SLOTS*CAND_W-1:0] dcand;
                for (gc = 0; gc < SLOTS; gc = gc + 1) begin : g_dcand
                    localparam integer NODE =
                        gc * (1 << (STAGES - tb - 1)) + (1 << (STAGES - tb - 2)) - 1;
                    if (NODE < NUM_SYND) begin : g_in
                        assign dcand[gc*CAND_W +: CAND_W] =
                            {iSynValid[NODE], iSynFlat[NODE*LOG_NUM_DET +: LOG_NUM_DET]};
                    end else begin : g_pad
                        assign dcand[gc*CAND_W +: CAND_W] = {1'b0, {LOG_NUM_DET{1'b0}}};
                    end
                end
                for (ti = 0; ti < L; ti = ti + 1) begin : g_st
                    localparam W    = 1 << (L + 1 - ti);     // candidates this stage holds
                    localparam OFF  = (1 << (L + 2)) - (1 << (L + 2 - ti));  // where they start
                    if (ti == 0) begin : g_first
                        reg [(W/2)*CAND_W-1:0] pick_c, pick_b;
                        integer k;
                        always @* begin
                            pick_c = {(W/2)*CAND_W{1'b0}};
                            pick_b = {(W/2)*CAND_W{1'b0}};
                            for (k = 0; k < (1 << MUX_T); k = k + 1) begin
                                if (!p_hot_c[k])
                                    pick_c = pick_c | dcand[k*(W*CAND_W)                +: (W/2)*CAND_W];
                                if (p_hot_b[k])
                                    pick_b = pick_b | dcand[k*(W*CAND_W) + (W/2)*CAND_W +: (W/2)*CAND_W];
                            end
                        end
                        always @(posedge iClk)
                            if (iEn) begin
                                q[OFF*CAND_W +: (W/2)*CAND_W]       <= pick_c;
                                q[(OFF+W/2)*CAND_W +: (W/2)*CAND_W] <= pick_b;
                            end
                    end else begin : g_next
                        localparam PREV = OFF - 2*W;         // the stage before, 2*W candidates
                        always @(posedge iClk)
                            if (iEn) begin
                                q[OFF*CAND_W +: (W/2)*CAND_W] <= (~p_go_c[ti-1])
                                    ? q[(PREV+W)*CAND_W +: (W/2)*CAND_W]
                                    : q[PREV*CAND_W +: (W/2)*CAND_W];
                                q[(OFF+W/2)*CAND_W +: (W/2)*CAND_W] <= p_go_b[ti-1]
                                    ? q[(PREV+W+W/2)*CAND_W +: (W/2)*CAND_W]
                                    : q[(PREV+W/2)*CAND_W +: (W/2)*CAND_W];
                            end
                    end
                end
                // the last stage holds four candidates; level tb keeps two of them
                localparam LAST = (1 << (L + 2)) - 8;
                assign look_pair[(tb-MUX_T-1)*2*CAND_W +: CAND_W] = (~p_go_c[L-1])
                    ? q[(LAST+2)*CAND_W +: CAND_W] : q[LAST*CAND_W +: CAND_W];
                assign look_pair[(tb-MUX_T-1)*2*CAND_W + CAND_W +: CAND_W] = p_go_b[L-1]
                    ? q[(LAST+3)*CAND_W +: CAND_W] : q[(LAST+1)*CAND_W +: CAND_W];
            end
        end
    endgenerate
    generate
        for (gs = 0; gs < STAGES; gs = gs + 1) begin : g_bs
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

            // ---- the probed slot: a 2:1 pick out of the pair the previous level selected ----
            wire [2*CAND_W-1:0] pair;
            if (gs == 0) begin : g_root
                // level 0 probes the fixed root node, so both halves are the root
                localparam integer ROOT = (1 << (STAGES - 1)) - 1;
                if (ROOT < NUM_SYND) begin : g_in
                    assign pair =
                        {2{{iSynValid[ROOT], iSynFlat[ROOT*LOG_NUM_DET +: LOG_NUM_DET]}}};
                end else begin : g_pad
                    assign pair = {2{{1'b0, {LOG_NUM_DET{1'b0}}}}};
                end
            end else begin : g_reg
                assign pair = p_pair[(gs-1)*2*CAND_W +: 2*CAND_W];
            end

            // lpath's low bit is the previous level's go_right, which picks the half
            wire [CAND_W-1:0]      csel = lpath[0] ? pair[CAND_W +: CAND_W] : pair[0 +: CAND_W];
            wire [LOG_NUM_DET-1:0] pid  = csel[0 +: LOG_NUM_DET];
            wire                   occ  = csel[LOG_NUM_DET];

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
            wire eq = occ && (pid == in_key);

            wire                       n_hit  = in_hit || eq;

            // first hit freezes the index: hidx_sel resolves before eq, so the update is one mux
            wire [LOG_NUM_SYND_P1-1:0] hidx_sel = in_hit ? in_hidx : node_idx;
            wire [LOG_NUM_SYND_P1-1:0] n_hidx   = eq ? hidx_sel : in_hidx;

            if (gs < STAGES-1) begin : g_mid
                // ---- the next level's slot pair, selected here in parallel with the compare ----
                wire [2*CAND_W-1:0] n_pair;
                if (LOOK && gs > MUX_T) begin : g_pick
                    assign n_pair = look_pair[(gs-MUX_T-1)*2*CAND_W +: 2*CAND_W];
                end else begin : g_mux
                    localparam NNEXT = 1 << (gs + 1);
                    wire [NNEXT*CAND_W-1:0] ncand;
                    for (gc = 0; gc < NNEXT; gc = gc + 1) begin : g_ncand
                        localparam integer NODE =
                            gc * (1 << (STAGES - gs - 1)) + (1 << (STAGES - gs - 2)) - 1;
                        if (NODE < NUM_SYND) begin : g_in
                            assign ncand[gc*CAND_W +: CAND_W] =
                                {iSynValid[NODE], iSynFlat[NODE*LOG_NUM_DET +: LOG_NUM_DET]};
                        end else begin : g_pad
                            // past the array: never occupied, so neither eq nor go_right can fire
                            assign ncand[gc*CAND_W +: CAND_W] = {1'b0, {LOG_NUM_DET{1'b0}}};
                        end
                    end
                    wire [((gs > 0) ? gs : 1)-1:0] lpath_b, lpath_c;
                    if (gs == 0) begin : g_lpb0
                        assign lpath_b = 1'b0;
                        assign lpath_c = 1'b0;
                    end else begin : g_lpbn
                        assign lpath_b =  p_path_b[(gs-1)*gs/2 +: gs];
                        assign lpath_c = ~p_path_c[(gs-1)*gs/2 +: gs];
                    end
                    // slots 2*lpath and 2*lpath+1 are adjacent, so one mux over pairs picks both
                    assign n_pair[0 +: CAND_W] =
                        ncand[lpath_c*(2*CAND_W)          +: CAND_W];
                    assign n_pair[CAND_W +: CAND_W] =
                        ncand[lpath_b*(2*CAND_W) + CAND_W +: CAND_W];
                end
                // the descend direction only exists where a next level consumes it
                wire go_right = occ && (pid < in_key);
                wire [gs:0] n_path;
                if (gs == 0) begin : g_np0
                    assign n_path = go_right;
                end else begin : g_npn
                    assign n_path = {lpath, go_right};
                end

                // path copies, written only for the banks a pair mux reads
                if (gs <= REP_TOP) begin : g_path_rep
                    always @(posedge iClk)
                        if (iEn) begin
                            p_path_b[gs*(gs+1)/2 +: gs+1] <=  n_path;
                            p_path_c[gs*(gs+1)/2 +: gs+1] <= ~n_path;
                        end
                end
                if (LOOK && gs == MUX_T - 1) begin : g_hot
                    wire [(1<<MUX_T)-1:0] n_hot = {{((1<<MUX_T)-1){1'b0}}, 1'b1} << n_path;
                    always @(posedge iClk)
                        if (iEn) begin
                            p_hot_b <=  n_hot;
                            p_hot_c <= ~n_hot;
                        end
                end
                if (LOOK && gs >= MUX_T && gs < PAIR_GS) begin : g_go_rep
                    always @(posedge iClk)
                        if (iEn) begin
                            p_go_b[gs-MUX_T] <=  go_right;
                            p_go_c[gs-MUX_T] <= ~go_right;
                        end
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
                        p_pair[gs*2*CAND_W +: 2*CAND_W]               <= n_pair;
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
