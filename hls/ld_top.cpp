#include "ld_microkernel.hpp"
// Vitis HLS top-level: calculate the complete loci-by-loci LD count matrix.
// `haplotypes` is a pre-packed, locus-major [loci][WORDS] buffer, where
// WORDS = ceil(LD_MAX_K / 64); unused words and tail bits must be zero.
// `counts` is row-major [loci][loci].  Tiling stays inside the kernel.
extern "C" void ld_top(const uint64_t *haplotypes, uint32_t *counts,
                       int loci, int samples) {
    #pragma HLS INTERFACE m_axi port=haplotypes offset=slave bundle=gmem0
    #pragma HLS INTERFACE m_axi port=counts offset=slave bundle=gmem1
    #pragma HLS INTERFACE s_axilite port=loci bundle=control
    #pragma HLS INTERFACE s_axilite port=samples bundle=control
    #pragma HLS INTERFACE s_axilite port=return bundle=control

    if (loci < 1 || loci > ld_cfg::MAX_LOCI ||
        samples < 1 || samples > ld_cfg::MAX_K) return;

    uint64_t a[ld_cfg::MR][ld_cfg::WORDS];
    uint64_t b[ld_cfg::NR][ld_cfg::WORDS];
    uint32_t c[ld_cfg::MR][ld_cfg::NR];
    #pragma HLS ARRAY_PARTITION variable=a complete dim=1
    #pragma HLS ARRAY_PARTITION variable=b complete dim=1
    #pragma HLS ARRAY_PARTITION variable=c complete dim=0

    for (int i0 = 0; i0 < ld_cfg::MAX_LOCI; i0 += ld_cfg::MR) {
        if (i0 >= loci) break;
        const int m = (loci - i0 < ld_cfg::MR) ? loci - i0 : ld_cfg::MR;
        for (int j0 = 0; j0 < ld_cfg::MAX_LOCI; j0 += ld_cfg::NR) {
            if (j0 >= loci) break;
            const int n = (loci - j0 < ld_cfg::NR) ? loci - j0 : ld_cfg::NR;

            for (int i = 0; i < ld_cfg::MR; ++i)
                for (int w = 0; w < ld_cfg::WORDS; ++w)
                    a[i][w] = (i < m)
                        ? haplotypes[(i0 + i) * ld_cfg::WORDS + w] : 0;
            for (int j = 0; j < ld_cfg::NR; ++j)
                for (int w = 0; w < ld_cfg::WORDS; ++w)
                    b[j][w] = (j < n)
                        ? haplotypes[(j0 + j) * ld_cfg::WORDS + w] : 0;

            ld_microkernel<ld_cfg::MR, ld_cfg::NR, ld_cfg::WORDS, ld_cfg::POP_WIDTH>(
                &a[0][0], &b[0][0], &c[0][0], m, n, samples);
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    counts[(i0 + i) * loci + j0 + j] = c[i][j];
        }
    }
}
