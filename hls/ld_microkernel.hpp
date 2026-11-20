#pragma once
#include <cstdint>
#include "ld_config.hpp"

// Fixed-stride packed panels A[MR][WORDS], B[NR][WORDS]; C[MR][NR].
// Each word holds 64 haplotypes for a single locus.
// This starter computes exact pairwise intersection counts only.
static inline uint32_t ld_chunk_popcount(uint64_t x) {
    uint32_t count = 0;
    for (int bit = 0; bit < LD_POP_WIDTH; ++bit) {
        #pragma HLS UNROLL
        count += static_cast<uint32_t>((x >> bit) & 1ULL);
    }
    return count;
}

template<int MR, int NR, int WORDS, int POP_WIDTH>
void ld_microkernel(const uint64_t *a, const uint64_t *b,
                    uint32_t *c, int m, int n, int samples) {
    uint64_t local_a[MR][WORDS];
    uint64_t local_b[NR][WORDS];
    #pragma HLS ARRAY_PARTITION variable=local_a complete dim=1
    #pragma HLS ARRAY_PARTITION variable=local_b complete dim=1
    for (int i = 0; i < MR; ++i)
        for (int w = 0; w < WORDS; ++w)
            local_a[i][w] = a[i * WORDS + w];
    for (int j = 0; j < NR; ++j)
        for (int w = 0; w < WORDS; ++w)
            local_b[j][w] = b[j * WORDS + w];

    for (int i = 0; i < MR; ++i) {
        for (int j = 0; j < NR; ++j) {
            uint32_t acc = 0;
            if (i < m && j < n) {
                for (int k = 0; k < samples; k += POP_WIDTH) {
                    #pragma HLS PIPELINE II=1
                    const int w = k / 64;
                    const int shift = k % 64;
                    uint64_t x = (local_a[i][w] & local_b[j][w]) >> shift;
                    const int remaining = samples - k;
                    if (remaining < POP_WIDTH)
                        x &= (uint64_t(1) << remaining) - 1ULL;
                    acc += ld_chunk_popcount(x);
                }
            }
            c[i * NR + j] = acc;
        }
    }
}

