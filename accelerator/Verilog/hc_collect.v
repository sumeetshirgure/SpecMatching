`ifndef hc_collect
`define hc_collect


module hc_collect #(
    parameter NUM_SYND  = 38,   // defects in the syndrome (sets the local index width)
    parameter WEIGHT_W  = 16,   // integer weight width
    parameter NUM_LANES = 64    // hits that can arrive in one cycle
) (
    input  wire                                    iClk,
    input  wire                                    iRstN,
    input  wire                                    iValid,
    input  wire [$clog2(NUM_SYND+1)-1:0]           iU,
    input  wire [NUM_LANES-1:0]                    iHitFlat,
    input  wire [NUM_LANES*$clog2(NUM_SYND+1)-1:0] iIdxFlat,
    input  wire [NUM_LANES*WEIGHT_W-1:0]           iWFlat,
    output wire                                    oReady,
    output wire                                    oBusy,
    output wire                                    oEdgeValid,
    output wire [$clog2(NUM_SYND+1)-1:0]           oEdgeU,
    output wire [$clog2(NUM_SYND+1)-1:0]           oEdgeV,
    output wire [WEIGHT_W-1:0]                     oEdgeW
);
    localparam LOG_NUM_SYND_P1 = $clog2(NUM_SYND+1);                        // local defect index width
    localparam SEL_W           = (NUM_LANES > 1) ? $clog2(NUM_LANES) : 1;   // lane selector width

    reg [NUM_LANES-1:0]          pend;       // hits of the held batch still to be emitted
    reg [LOG_NUM_SYND_P1-1:0]          pend_u;
    reg [NUM_LANES*LOG_NUM_SYND_P1-1:0]    pend_idx;
    reg [NUM_LANES*WEIGHT_W-1:0] pend_w;

    wire [NUM_LANES-1:0] sel  = pend & (~pend + {{(NUM_LANES-1){1'b0}}, 1'b1});
    wire [NUM_LANES-1:0] rem  = pend & ~sel;
    wire             take = oReady && iValid;

    assign oReady = ~|rem;
    assign oBusy  = |pend;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN)   pend <= {NUM_LANES{1'b0}};
        else if (take)  pend <= iHitFlat;
        else            pend <= rem;

    always @(posedge iClk)
        if (take) begin
            pend_u   <= iU;
            pend_idx <= iIdxFlat;
            pend_w   <= iWFlat;
        end

    wire [SEL_W-1:0]       sel_idx;
    wire [SEL_W*NUM_LANES-1:0] enc_sel;

    genvar gk, gb;
    generate
        for (gb = 0; gb < SEL_W; gb = gb + 1) begin : g_enc
            for (gk = 0; gk < NUM_LANES; gk = gk + 1) begin : g_sel
                if (((gk >> gb) & 1) != 0) begin : g_one
                    assign enc_sel[gb*NUM_LANES + gk] = sel[gk];
                end else begin : g_zero
                    assign enc_sel[gb*NUM_LANES + gk] = 1'b0;
                end
            end
            assign sel_idx[gb] = |enc_sel[gb*NUM_LANES +: NUM_LANES];
        end
    endgenerate

    assign oEdgeValid = |pend;
    assign oEdgeU     = pend_u;
    assign oEdgeV     = pend_idx[sel_idx*LOG_NUM_SYND_P1 +: LOG_NUM_SYND_P1];
    assign oEdgeW     = pend_w[sel_idx*WEIGHT_W +: WEIGHT_W];
endmodule
`endif
