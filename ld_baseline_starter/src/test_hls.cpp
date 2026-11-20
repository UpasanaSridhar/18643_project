#include "ld.hpp"
#include "../hls/ld_config.hpp"
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
extern "C" void ld_top(const uint64_t *, uint32_t *, int, int);

namespace {
constexpr int kTrials = 100;
constexpr int kRunsPerTrial = 100;
}

int main(int argc, char **argv) {
    const int loci = argc > 1 ? std::atoi(argv[1]) : 10;
    const int samples = argc > 2 ? std::atoi(argv[2]) : 10;
    const unsigned seed = argc > 3 ? static_cast<unsigned>(std::strtoul(argv[3], nullptr, 10)) : 17;
    if (loci < 1 || loci > ld_cfg::MAX_LOCI || samples < 1 || samples > ld_cfg::MAX_K) {
        std::cerr << "Usage: test_hls loci<=LD_MAX_LOCI samples<=LD_MAX_K seed\n";
        return 2;
    }
    const auto h = make_input(loci, samples, seed);
    const auto reference = gemm(h, h, loci, loci, samples);
    std::vector<uint32_t> result(loci * loci, 0);
    // The HLS top-level receives an already packed fixed-stride buffer.
    std::vector<uint64_t> packed(loci * ld_cfg::WORDS, 0);
    const std::size_t input_stride = packed_stride(samples);
    for (int locus = 0; locus < loci; ++locus)
        for (std::size_t word = 0; word < input_stride; ++word)
            packed[locus * ld_cfg::WORDS + word] = h[locus * input_stride + word];

    // Validate one invocation before measuring it.
    ld_top(packed.data(), result.data(), loci, samples);
    int failures = 0;
    for (int i = 0; i < loci; ++i)
        for (int j = 0; j < loci; ++j) {
            const int at = i*loci+j;
            if (result[at] != reference[at]) {
                std::cerr << "MISMATCH (" << i << "," << j << ") ref="
                          << reference[at]
                          << " hls=" << result[at] << "\n";
                ++failures;
            }
            const double r0 = r_squared(reference[at], reference[i*loci+i], reference[j*loci+j], samples);
            const double r1 = r_squared(result[at], result[i*loci+i], result[j*loci+j], samples);
            if (!((std::isnan(r0) && std::isnan(r1)) || std::abs(r0-r1) < 1e-12)) ++failures;
        }

    // Measure 100 batches of 100 calls. Trim exactly one fastest and one
    // slowest batch, then report the average elapsed time per trial.
    double trial_total_us = 0.0;
    double min_trial_us = std::numeric_limits<double>::infinity();
    double max_trial_us = 0.0;
    volatile uint32_t timing_sink = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
        const auto start = std::chrono::steady_clock::now();
        for (int run = 0; run < kRunsPerTrial; ++run)
            ld_top(packed.data(), result.data(), loci, samples);
        const auto elapsed = std::chrono::steady_clock::now() - start;
        const double trial_us = std::chrono::duration<double, std::micro>(elapsed).count();
        trial_total_us += trial_us;
        min_trial_us = std::min(min_trial_us, trial_us);
        max_trial_us = std::max(max_trial_us, trial_us);
        timing_sink ^= result[static_cast<std::size_t>(trial) % result.size()];
    }
    const double trimmed_mean_trial_us =
        (trial_total_us - min_trial_us - max_trial_us) / (kTrials - 2);
    std::cout << (failures ? "FAIL" : "PASS")
              << " mode=popcount"
              << " MR=" << ld_cfg::MR << " NR=" << ld_cfg::NR
              << " MAX_K=" << ld_cfg::MAX_K << " POP_WIDTH=" << ld_cfg::POP_WIDTH
              << " loci=" << loci << " haplotypes=" << samples
              << " tiles=" << ((loci+ld_cfg::MR-1)/ld_cfg::MR)*((loci+ld_cfg::NR-1)/ld_cfg::NR)
              << " runtime_us_per_trial=" << std::fixed << std::setprecision(3)
              << trimmed_mean_trial_us
              << " runtime_trials=" << kTrials
              << " runs_per_trial=" << kRunsPerTrial
              << " trim=min,max"
              << " mismatches=" << failures << "\n";
    (void)timing_sink;
    if (loci <= 10) {
        std::cout << "Intersection counts:\n";
        for (int i=0;i<loci;++i) {
            for (int j=0;j<loci;++j) std::cout << std::setw(4) << result[i*loci+j];
            std::cout << '\n';
        }
    }
    return failures ? 1 : 0;
}
