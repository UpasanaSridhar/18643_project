#pragma once
#include <cstdint>
#include <vector>
#include <cstddef>

// Number of packed 64-bit words used by one locus.
std::size_t packed_stride(std::size_t samples);

// Returns a locus-major packed input array.  Element (locus, haplotype) is
// stored in bit (haplotype % 64) of words[locus * packed_stride(samples) +
// haplotype / 64].
std::vector<uint64_t> make_input(std::size_t loci, std::size_t samples,
                                 uint32_t seed);

// Boolean GEMM over packed haplotypes.  a has m packed rows and b has n packed
// rows, each with packed_stride(samples) words.  The returned flat m x n array
// contains unsigned 64-bit intersection counts in row-major order.
std::vector<uint64_t> gemm(const std::vector<uint64_t> &a,
                           const std::vector<uint64_t> &b,
                           std::size_t m, std::size_t n,
                           std::size_t samples);

double r_squared(uint64_t both, uint64_t a, uint64_t b, std::size_t samples);
