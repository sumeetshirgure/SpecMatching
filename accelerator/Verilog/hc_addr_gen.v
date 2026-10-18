`ifndef hc_addr_gen
`define hc_addr_gen

module hc_addr_gen #(
    parameter NUM_SYND        = 38,   // syndrome capacity (max H nodes)
    parameter LOG_NUM_DET     = 12,   // detector id width = ball-table address width
    parameter LOG_MAX_BALL    = 9,    // raw ball length width
    parameter NUM_LANES       = 64,   // entries per chunk
    parameter CHUNKS_PER_BALL = 4     // NUM_LANES-wide chunks per ball-table word
) (
    input  wire                                                               iClk,
    input  wire                                                               iRstN,
    input  wire                                                               iGo,
    input  wire                                                               iStall,
    input  wire [$clog2(NUM_SYND+1)-1:0]                                      iM,
    input  wire [NUM_SYND*LOG_NUM_DET-1:0]                                    iSynFlat,
    output wire                                                               oBndTabRdEn,
    output wire [LOG_NUM_DET-1:0]                                             oBndTabAddr,
    input  wire [LOG_MAX_BALL-1:0]                                            iBndSize,
    output wire                                                               oBallRdEn,
    output wire [LOG_NUM_DET-1:0]                                             oBallAddr,
    output wire                                                               oChValid,
    output wire [$clog2(NUM_SYND+1)-1:0]                                      oRdU,
    output wire [LOG_NUM_DET-1:0]                                             oRdDet,
    output wire [(($clog2(CHUNKS_PER_BALL)>0)?$clog2(CHUNKS_PER_BALL):1)-1:0] oRdCh,
    output wire [$clog2(NUM_LANES)-1:0]                                       oRdHi,
    output wire                                                               oHdrValid,
    output wire [$clog2(NUM_SYND+1)-1:0]                                      oHdrU,
    output wire                                                               oBusy
);
    localparam LOG_NUM_SYND_P1 = $clog2(NUM_SYND+1);                                            // local defect index width
    localparam LOG_NUM_LANES   = $clog2(NUM_LANES);                                             // lane bound width
    localparam CIDX_W          = ($clog2(CHUNKS_PER_BALL) > 0) ? $clog2(CHUNKS_PER_BALL) : 1;
    localparam CCNT_W          = $clog2(CHUNKS_PER_BALL + 1);                                   // holds 0 .. CHUNKS_PER_BALL
    
    localparam SZX_W = ((LOG_MAX_BALL > LOG_NUM_LANES) ? LOG_MAX_BALL : LOG_NUM_LANES) + 1;

    localparam [1:0] S_IDLE = 2'd0,
                     S_NODE = 2'd1,   // boundary SRAM read issued for syn[0]
                     S_HDR  = 2'd2,   // boundary word ready: load ball 0
                     S_SCAN = 2'd3;   // one chunk per unstalled cycle

    reg [1:0]               state;
    reg [LOG_NUM_SYND_P1-1:0]  i;       // defect whose ball is being scanned
    reg [LOG_NUM_SYND_P1-1:0]  ni;      // next defect whose boundary word is to be read
    reg [CCNT_W-1:0]        ccnt;    // chunks left in the current ball
    reg [CIDX_W-1:0]        cidx;    // chunk index into the held ball word
    reg                     first_c; // the chunk being issued is its ball's first
    reg [LOG_NUM_LANES-1:0]    bhi;     // last entry's lane offset: bound on the last chunk

    wire [LOG_NUM_DET-1:0] det_i = iSynFlat[i*LOG_NUM_DET +: LOG_NUM_DET];

    wire [LOG_NUM_SYND_P1-1:0] last_i = iM - 1'b1;
    wire                    m_zero = (iM == {LOG_NUM_SYND_P1{1'b0}});

    wire empty  = (ccnt == {CCNT_W{1'b0}});
    wire last_c = empty || (ccnt == {{(CCNT_W-1){1'b0}}, 1'b1});
    wire step   = (state == S_SCAN) && !iStall;
    wire load   = (state == S_HDR) || (step && last_c && (i != last_i));

    assign oBndTabRdEn = (state == S_NODE) || (load && (ni != iM));
    assign oBndTabAddr = iSynFlat[ni*LOG_NUM_DET +: LOG_NUM_DET];
    assign oBallRdEn   = step && !empty && first_c;
    assign oBallAddr    = det_i;
    assign oChValid     = step && !empty;
    assign oRdU         = i;
    assign oRdDet       = det_i;
    assign oRdCh        = cidx;
    assign oRdHi        = last_c ? bhi : {LOG_NUM_LANES{1'b1}};
    assign oBusy         = (state != S_IDLE);

    wire n_empty = (iBndSize == {LOG_MAX_BALL{1'b0}});

    wire [SZX_W-1:0] size_x = {{(SZX_W-LOG_MAX_BALL){1'b0}}, iBndSize};
    wire [SZX_W-1:0] ch_full = (size_x + {{(SZX_W-LOG_NUM_LANES){1'b0}}, {LOG_NUM_LANES{1'b1}}})
                               >> LOG_NUM_LANES;
    wire [SZX_W-1:0] size_m1 = size_x - 1'b1;
    wire [CCNT_W-1:0]      n_ch = n_empty ? {CCNT_W{1'b0}} : ch_full[CCNT_W-1:0];
    wire [LOG_NUM_LANES-1:0]  n_hi = size_m1[LOG_NUM_LANES-1:0];

    // ---------------- control ----------------
    always @(posedge iClk or negedge iRstN)
        if (!iRstN) begin
            state   <= S_IDLE;
            i       <= {LOG_NUM_SYND_P1{1'b0}};
            ni      <= {LOG_NUM_SYND_P1{1'b0}};
            ccnt    <= {CCNT_W{1'b0}};
            cidx    <= {CIDX_W{1'b0}};
            first_c <= 1'b0;
        end else begin
            case (state)
                S_IDLE:
                    if (iGo && !m_zero) begin
                        i     <= {LOG_NUM_SYND_P1{1'b0}};
                        ni    <= {LOG_NUM_SYND_P1{1'b0}};
                        state <= S_NODE;
                    end
                S_NODE:
                    state <= S_HDR;
                S_HDR:
                    state <= S_SCAN;
                S_SCAN:
                    if (!iStall) begin
                        ccnt    <= ccnt - 1'b1;
                        cidx    <= cidx + 1'b1;
                        first_c <= 1'b0;
                        if (last_c) begin
                            if (i == last_i) state <= S_IDLE;
                            else             i <= i + 1'b1;
                        end
                    end
                default:
                    state <= S_IDLE;
            endcase
            if (oBndTabRdEn) ni <= ni + 1'b1;
            if (load) begin
                ccnt    <= n_ch;
                cidx    <= {CIDX_W{1'b0}};
                first_c <= 1'b1;
                bhi     <= n_hi;
            end
        end

    assign oHdrValid = load;
    assign oHdrU     = ni - 1'b1;
endmodule
`endif
