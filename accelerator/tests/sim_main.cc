// Verilated harness for the h_construct cross-check.
//
//     hc_sim <vector-dir> <num-shots>
//
// For each shot it resets the DUT, drives the shot's syndrome, runs the build to oDone and writes
// the emitted H to <vector-dir>/dut_<k>.txt as "E u v w" / "B u w". Order is not canonical;
// accelerator/tests/compare.py normalises it and keeps multiplicity.
//
// Must run with its working directory set to the vector directory: ball_sram's $readmemh resolves
// INIT_FILE ("ball.hex" / "bnd.hex") against the process CWD, and it fires when the model is
// constructed. accelerator/tests/run_crosscheck.sh cds there before invoking this.
//
// Exits non-zero on a missing/!malformed vector file or on a watchdog timeout, so a DUT that never
// raises oDone fails as a timeout rather than as a short dump the comparator would blame on H.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "verilated.h"

#include "Vhc_top.h"

namespace {

// A build is ~26 chunk cycles plus the lane pipeline, and even a pathological shot in which every
// batch drains 64 hits one per cycle stays under ~1700. This is an order of magnitude above that:
// large enough never to cut a legitimate build short, small enough that a hang fails in a moment.
constexpr long WATCHDOG_CYCLES = 20000;

/// Reads `` `define <name> <integer> `` out of params.vh.
///
/// The harness takes the syndrome geometry from the same params.vh that fixed the port widths of
/// the model it is linked against, so the packing below cannot drift from the compiled geometry.
bool read_define(const std::string& path, const char* name, long* out) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (f == nullptr)
        return false;
    char line[512];
    std::string key = std::string("`define ") + name + " ";
    bool found = false;
    while (std::fgets(line, sizeof(line), f) != nullptr) {
        if (std::string(line).compare(0, key.size(), key) == 0) {
            *out = std::strtol(line + key.size(), nullptr, 10);
            found = true;
            break;
        }
    }
    std::fclose(f);
    return found;
}

/// `M` followed by `M` ascending detector ids.
bool read_syndrome(const std::string& path, std::vector<long>& dets) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (f == nullptr)
        return false;
    long m = 0;
    if (std::fscanf(f, "%ld", &m) != 1 || m < 0) {
        std::fclose(f);
        return false;
    }
    dets.clear();
    for (long i = 0; i < m; i++) {
        long d = 0;
        if (std::fscanf(f, "%ld", &d) != 1) {
            std::fclose(f);
            return false;
        }
        dets.push_back(d);
    }
    std::fclose(f);
    return true;
}

/// Packs the syndrome into iSynFlat, slot `j` at bit `j * log_num_det`.
///
/// iSynFlat is NUM_SYND * LOG_NUM_DET bits (91 for the generated table), so Verilator exposes it
/// as a wide signal and it has to be written a 32-bit word at a time rather than as one integer.
/// Slots at or above M are left zero: h_construct's syn_valid mask marks them unused, so hc_lane
/// treats them as +infinity and never matches them.
void drive_syndrome(Vhc_top* top, const std::vector<long>& dets, long log_num_det, long syn_words) {
    for (long w = 0; w < syn_words; w++)
        top->iSynFlat[w] = 0;
    for (size_t j = 0; j < dets.size(); j++) {
        uint32_t id = (uint32_t)dets[j];
        for (long b = 0; b < log_num_det; b++) {
            if (((id >> b) & 1u) == 0)
                continue;
            long bit = (long)j * log_num_det + b;
            top->iSynFlat[bit / 32] |= (uint32_t)1u << (bit % 32);
        }
    }
}

/// One clock level, then settle. Every input change is followed by an eval, and both edges are
/// evaluated, so the model sees the same waveform the event simulator did.
void clock_level(Vhc_top* top, int level) {
    top->iClk = level;
    top->eval();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: hc_sim <vector-dir> <num-shots>\n");
        return 2;
    }
    const std::string vec_dir = argv[1];
    const long num_shots = std::strtol(argv[2], nullptr, 10);
    if (num_shots < 1) {
        std::fprintf(stderr, "hc_sim: num-shots must be >= 1\n");
        return 2;
    }

    const std::string params = vec_dir + "/params.vh";
    long num_synd = 0;
    long log_num_det = 0;
    if (!read_define(params, "HC_NUM_SYND", &num_synd) || !read_define(params, "HC_LOG_NUM_DET", &log_num_det)) {
        std::fprintf(stderr, "hc_sim: cannot read HC_NUM_SYND / HC_LOG_NUM_DET from %s\n", params.c_str());
        return 2;
    }
    const long syn_words = (num_synd * log_num_det + 31) / 32;

    VerilatedContext* ctx = new VerilatedContext;
    ctx->commandArgs(argc, argv);
    // Constructing the model runs ball_sram's $readmemh, which needs CWD == vector dir.
    Vhc_top* top = new Vhc_top{ctx};

    int status = 0;
    for (long k = 0; k < num_shots && status == 0; k++) {
        std::vector<long> dets;
        const std::string syn_path = vec_dir + "/syn_" + std::to_string(k) + ".txt";
        if (!read_syndrome(syn_path, dets)) {
            std::fprintf(stderr, "hc_sim: cannot read syndrome file %s\n", syn_path.c_str());
            status = 2;
            break;
        }
        if ((long)dets.size() > num_synd) {
            std::fprintf(
                stderr, "hc_sim: shot %ld has M = %zu, above NUM_SYND = %ld\n", k, dets.size(), num_synd);
            status = 2;
            break;
        }

        // ---- reset
        top->iRstN = 0;
        top->iStart = 0;
        top->iM = 0;
        drive_syndrome(top, std::vector<long>(), log_num_det, syn_words);
        top->eval();
        for (int i = 0; i < 4; i++) {
            clock_level(top, 1);
            clock_level(top, 0);
        }
        top->iRstN = 1;
        for (int i = 0; i < 2; i++) {
            clock_level(top, 1);
            clock_level(top, 0);
        }

        // ---- drive the shot and raise iStart for exactly one posedge
        top->iM = (uint8_t)dets.size();
        drive_syndrome(top, dets, log_num_det, syn_words);
        top->iStart = 1;
        top->eval();

        std::vector<std::string> records;
        bool done = false;
        bool start_cleared = false;
        long cycles = 0;
        for (; cycles < WATCHDOG_CYCLES; cycles++) {
            // After the posedge eval the outputs hold this cycle's values, so each asserted cycle
            // contributes exactly one record -- the rate hc_collect drains hits at and the rate
            // hc_bnd_check pulses at.
            clock_level(top, 1);
            if (top->oEdgeValid) {
                unsigned u = top->oEdgeU;
                unsigned v = top->oEdgeV;
                records.push_back(
                    "E " + std::to_string(u < v ? u : v) + " " + std::to_string(u < v ? v : u) + " " +
                    std::to_string((unsigned)top->oEdgeW));
            }
            if (top->oBndValid)
                records.push_back(
                    "B " + std::to_string((unsigned)top->oBndU) + " " + std::to_string((unsigned)top->oBndW));
            done = top->oDone;
            clock_level(top, 0);
            if (!start_cleared) {
                top->iStart = 0;
                top->eval();
                start_cleared = true;
            }
            if (done)
                break;
        }

        if (!done) {
            std::fprintf(
                stderr, "hc_sim: timeout after %ld cycles on shot %ld without oDone\n", cycles, k);
            status = 1;
            break;
        }

        const std::string out_path = vec_dir + "/dut_" + std::to_string(k) + ".txt";
        FILE* out = std::fopen(out_path.c_str(), "w");
        if (out == nullptr) {
            std::fprintf(stderr, "hc_sim: cannot open %s for writing\n", out_path.c_str());
            status = 2;
            break;
        }
        size_t n_edge = 0;
        for (const std::string& r : records) {
            std::fprintf(out, "%s\n", r.c_str());
            if (r[0] == 'E')
                n_edge++;
        }
        std::fclose(out);
        std::printf(
            "hc_sim: shot %-3ld M=%-3zu edges=%-4zu boundary=%-3zu cycles=%ld\n", k, dets.size(), n_edge,
            records.size() - n_edge, cycles + 1);
    }

    top->final();
    delete top;
    delete ctx;
    return status;
}
