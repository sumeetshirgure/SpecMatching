`ifndef ball_sram
`define ball_sram

module ball_sram #(
    parameter DEPTH = 1024,
    parameter WIDTH = 32,
    parameter INIT_FILE = "mem.hex"
) (
    input  wire                     iClk,
    input  wire                     iEn,
    input  wire [$clog2(DEPTH)-1:0] iAddr,
    output reg  [WIDTH-1:0]         oRdata
);
    reg [WIDTH-1:0] mem [0:DEPTH-1];
    initial $readmemh(INIT_FILE, mem);

    always @(posedge iClk)
        if (iEn) oRdata <= mem[iAddr];
endmodule
`endif
