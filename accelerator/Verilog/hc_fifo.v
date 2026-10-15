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

    localparam [CNT_W-1:0] FIFO_HIGH = FIFO_DEPTH - 1;
    localparam [CNT_W-1:0] CNT_ZERO  = 0;

    localparam [1:0] S_IDLE = 2'd0,
                     S_RUN  = 2'd1,
                     S_DONE = 2'd2;

    // --------- align the chunk sideband with the ball word held on the SRAM output ---------
    // The loader issues one ball-table read per non-empty ball and the SRAM holds the word
    // until the next read, so the chunk itself is cut out of iBallRdata here, one cycle
    // after its sideband was issued, with no holding register of its own.
    reg              rd_v;     // iBallRdata carries a live chunk this cycle
    reg [SB_W-1:0]   rd_sb;
    reg [CIDX_W-1:0] rd_ch;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) rd_v <= 1'b0;
        else        rd_v <= iChValid;

    always @(posedge iClk) begin
        rd_sb <= {iRdU, iRdDet, iRdHi};
        rd_ch <= iRdCh;
    end


    reg [CHUNK_W-1:0] chunk;
    reg [CIDX_W:0]    ci;   // one bit wider than the index so the bound test cannot wrap

    always @* begin
        chunk = iBallRdata[0 +: CHUNK_W];
        for (ci = 1; ci < CHUNKS_PER_BALL; ci = ci + 1)
            if (rd_ch == ci[CIDX_W-1:0])
                chunk = iBallRdata[ci[CIDX_W-1:0]*CHUNK_W +: CHUNK_W];
    end

    // ---------------- batch FIFO ----------------
    reg [SLOT_W-1:0] mem [0:FIFO_DEPTH-1];
    reg [PTR_W-1:0]  wptr, rptr;
    reg [CNT_W-1:0]  cnt;

    wire empty = (cnt == CNT_ZERO);
    wire push  = rd_v;
    wire pop   = !empty && iLaneReady;

    always @(posedge iClk or negedge iRstN)
        if (!iRstN) begin
            wptr <= {PTR_W{1'b0}};
            rptr <= {PTR_W{1'b0}};
            cnt  <= CNT_ZERO;
        end else begin
            if (push) wptr <= wptr + 1'b1;
            if (pop)  rptr <= rptr + 1'b1;
            case ({push, pop})
                2'b10:   cnt <= cnt + 1'b1;
                2'b01:   cnt <= cnt - 1'b1;
                default: cnt <= cnt;
            endcase
        end

    always @(posedge iClk)
        if (push) mem[wptr] <= {chunk, rd_sb};

    wire [SLOT_W-1:0] head = mem[rptr];

    assign oBatchValid = !empty;
    assign oBatchWord  = head[SB_W +: CHUNK_W];
    assign oBatchU     = head[LOG_NUM_DET+LOG_NUM_LANES +: LOG_NUM_SYND_P1];
    assign oBatchDet   = head[LOG_NUM_LANES +: LOG_NUM_DET];
    assign oBatchHi    = head[0 +: LOG_NUM_LANES];
    assign oStall       = (cnt >= FIFO_HIGH);

    // ---------------- run / done ----------------
    wire quiet = !iLoadBusy && !rd_v && empty &&
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
