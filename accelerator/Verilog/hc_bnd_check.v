`ifndef hc_bnd_check
`define hc_bnd_check

module hc_bnd_check #(
    parameter NUM_SYND = 38,   // syndrome capacity, max defects per shot (max H nodes)
    parameter WEIGHT_W = 16    // integer weight width
) (
    input  wire                          iClk,
    input  wire                          iRstN,
    input  wire                          iHdrValid,
    input  wire [$clog2(NUM_SYND+1)-1:0] iHdrU,
    input  wire                          iHb,
    input  wire [WEIGHT_W-1:0]           iBw,
    output reg                           oBndValid,
    output reg  [$clog2(NUM_SYND+1)-1:0] oBndU,
    output reg  [WEIGHT_W-1:0]           oBndW,
    output wire                          oBusy
);
    always @(posedge iClk or negedge iRstN)
        if (!iRstN) oBndValid <= 1'b0;
        else        oBndValid <= iHdrValid && iHb;

    always @(posedge iClk) begin
        oBndU <= iHdrU;
        oBndW <= iBw;
    end

    assign oBusy = iHdrValid || oBndValid;
endmodule
`endif
