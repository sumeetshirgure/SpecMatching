// Top wrapper for the h_construct cross-check under Verilator: h_construct plus the two table
// SRAMs, with nothing but the build interface exposed.
//
// (A comment whose first word is the simulator's own name is read as a pragma, so the wording
// above keeps it out of the first position.)
//
// The C++ harness (accelerator/tests/sim_main.cc) drives iStart/iM/iSynFlat and records the
// oEdge*/oBnd* stream, so the SRAM read ports are internal here and the harness never sees them.
// ball_sram holds oRdata while iEn is low, which is the behaviour hc_addr_gen's one-read-per-ball
// scheme depends on; both instances are the unmodified RTL model, and $readmemh resolves
// INIT_FILE against the simulator's working directory, so the harness runs with its CWD set to
// the vector directory.
//
// Parameters default to the generated table's geometry (params.vh), which is also what fixes this
// module's port widths -- iSynFlat is NUM_SYND*LOG_NUM_DET bits wide.
//
// Lint/sim: see accelerator/tests/run_crosscheck.sh for the verilator command line.

`ifndef hc_top
`define hc_top

`include "params.vh"

module hc_top #(
    parameter NUM_SYND        = `HC_NUM_SYND,
    parameter LOG_NUM_DET     = `HC_LOG_NUM_DET,
    parameter LOG_MAX_BALL    = `HC_LOG_MAX_BALL,
    parameter WEIGHT_W        = `HC_WEIGHT_W,
    parameter NUM_LANES       = `HC_NUM_LANES,
    parameter LANE_W          = `HC_LANE_W,
    parameter CHUNK_W         = `HC_CHUNK_W,
    parameter CHUNKS_PER_BALL = `HC_CHUNKS_PER_BALL,
    parameter SRAM_DEPTH      = `HC_SRAM_DEPTH,
    parameter BND_WIDTH       = `HC_BND_WIDTH,
    parameter BALL_WIDTH      = `HC_BALL_WIDTH
) (
    input  wire                            iClk,
    input  wire                            iRstN,
    input  wire                            iStart,
    input  wire [$clog2(NUM_SYND+1)-1:0]   iM,
    input  wire [NUM_SYND*LOG_NUM_DET-1:0] iSynFlat,
    output wire                            oEdgeValid,
    output wire [$clog2(NUM_SYND+1)-1:0]   oEdgeU,
    output wire [$clog2(NUM_SYND+1)-1:0]   oEdgeV,
    output wire [WEIGHT_W-1:0]             oEdgeW,
    output wire                            oBndValid,
    output wire [$clog2(NUM_SYND+1)-1:0]   oBndU,
    output wire [WEIGHT_W-1:0]             oBndW,
    output wire                            oDone
);

    wire                   bnd_rd_en;
    wire [LOG_NUM_DET-1:0] bnd_addr;
    wire [BND_WIDTH-1:0]   bnd_rdata;
    wire                   ball_rd_en;
    wire [LOG_NUM_DET-1:0] ball_addr;
    wire [BALL_WIDTH-1:0]  ball_rdata;

    ball_sram #(
        .DEPTH(SRAM_DEPTH), .WIDTH(BND_WIDTH), .INIT_FILE(`HC_BND_HEX)
    ) u_bnd_sram (
        .iClk(iClk), .iEn(bnd_rd_en), .iAddr(bnd_addr), .oRdata(bnd_rdata)
    );

    ball_sram #(
        .DEPTH(SRAM_DEPTH), .WIDTH(BALL_WIDTH), .INIT_FILE(`HC_BALL_HEX)
    ) u_ball_sram (
        .iClk(iClk), .iEn(ball_rd_en), .iAddr(ball_addr), .oRdata(ball_rdata)
    );

    h_construct #(
        .NUM_SYND(NUM_SYND),
        .LOG_NUM_DET(LOG_NUM_DET),
        .LOG_MAX_BALL(LOG_MAX_BALL),
        .WEIGHT_W(WEIGHT_W),
        .NUM_LANES(NUM_LANES),
        .LANE_W(LANE_W),
        .CHUNK_W(CHUNK_W),
        .CHUNKS_PER_BALL(CHUNKS_PER_BALL)
    ) u_dut (
        .iClk(iClk),
        .iRstN(iRstN),
        .iStart(iStart),
        .iM(iM),
        .iSynFlat(iSynFlat),
        .oBndTabRdEn(bnd_rd_en),
        .oBndTabAddr(bnd_addr),
        .iBndTabRdata(bnd_rdata),
        .oBallRdEn(ball_rd_en),
        .oBallAddr(ball_addr),
        .iBallRdata(ball_rdata),
        .oEdgeValid(oEdgeValid),
        .oEdgeU(oEdgeU),
        .oEdgeV(oEdgeV),
        .oEdgeW(oEdgeW),
        .oBndValid(oBndValid),
        .oBndU(oBndU),
        .oBndW(oBndW),
        .oDone(oDone)
    );

endmodule
`endif
