#include "ld.hpp"
#include <cmath>
#include <random>
#include <stdexcept>

std::size_t packed_stride(std::size_t samples) {
    return (samples + 63) / 64;
}

std::vector<uint64_t> make_input(std::size_t loci, std::size_t samples,
                                 uint32_t seed) {
    if (!loci || !samples) throw std::invalid_argument("loci and samples must be positive");
    const std::size_t stride = packed_stride(samples);
    std::vector<uint64_t> words(loci * stride, 0);
    std::mt19937 rng(seed);
    std::bernoulli_distribution bernoulli(0.35);
    for (std::size_t i = 0; i < loci; ++i)
        for (std::size_t k = 0; k < samples; ++k)
            if (bernoulli(rng)) words[i * stride + k / 64] |= 1ULL << (k % 64);
    // Explicit edge cases for correctness; no missing calls.
    if (loci >= 1)
        for (std::size_t k = 0; k < samples; ++k)
            words[k / 64] &= ~(1ULL << (k % 64));
    if (loci >= 2)
        for (std::size_t k = 0; k < samples; ++k)
            words[stride + k / 64] |= 1ULL << (k % 64);
    if (loci >= 4)
        for (std::size_t k = 0; k < samples; ++k) {
            const uint64_t mask = 1ULL << (k % 64);
            uint64_t &destination = words[3 * stride + k / 64];
            if (words[2 * stride + k / 64] & mask) destination |= mask;
            else destination &= ~mask;
        }
    return words;
}

std::vector<uint64_t> gemm(const std::vector<uint64_t> &a,
                           const std::vector<uint64_t> &b,
                           std::size_t m, std::size_t n,
                           std::size_t samples) {
    const std::size_t stride = packed_stride(samples);
    if (!samples || a.size() != m * stride || b.size() != n * stride)
        throw std::invalid_argument("invalid packed GEMM dimensions");

    std::vector<uint64_t> result(m * n, 0);
    const uint64_t final_word_mask = samples % 64
        ? (1ULL << (samples % 64)) - 1ULL
        : ~uint64_t{0};
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t w = 0; w < stride; ++w) {
                uint64_t overlap = a[i * stride + w] & b[j * stride + w];
                if (w + 1 == stride) overlap &= final_word_mask;
                result[i * n + j] += static_cast<uint64_t>(__builtin_popcountll(overlap));
            }
    return result;
}

double r_squared(uint64_t both, uint64_t a, uint64_t b, std::size_t n) {
    const double denom = double(a) * (n - a) * double(b) * (n - b);
    if (denom == 0) return std::nan(""); // monomorphic: undefined
    const double numerator = double(n) * both - double(a) * b;
    return numerator * numerator / denom;
}
