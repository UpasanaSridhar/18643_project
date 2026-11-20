#include "ld.hpp"
#include "../hls/ld_config.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace {

bool bit_at(const std::vector<uint64_t>& packed, std::size_t stride,
            int locus, int haplotype) {
    return (packed[static_cast<std::size_t>(locus) * stride + haplotype / 64] >>
            (haplotype % 64)) & 1ULL;
}

} // namespace

int main(int argc, char** argv) {
    const int loci = argc > 1 ? std::atoi(argv[1]) : 10;
    const int samples = argc > 2 ? std::atoi(argv[2]) : 10;
    const unsigned seed = argc > 3
        ? static_cast<unsigned>(std::strtoul(argv[3], nullptr, 10)) : 17;
    const char* output_path = argc > 4 ? argv[4] : "hls_generated.vcf";

    if (loci < 1 || loci > ld_cfg::MAX_LOCI ||
        samples < 1 || samples > ld_cfg::MAX_K) {
        std::cerr << "Usage: generate_hls_vcf loci<=LD_MAX_LOCI "
                  << "haplotypes<=LD_MAX_K seed output.vcf\n";
        return 2;
    }

    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "Could not write " << output_path << '\n';
        return 2;
    }

    const auto haplotypes = make_input(loci, samples, seed);
    const std::size_t stride = packed_stride(samples);

    output << "##fileformat=VCFv4.2\n"
           << "##source=ld_baseline_starter_hls_fixture\n"
           << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int sample = 0; sample < samples; ++sample)
        output << "\thap" << sample;
    output << '\n';

    // Haploid calls preserve every bit of the HLS locus-by-haplotype matrix.
    // AF is intentionally non-extreme: quickLD uses it for validation, while
    // deriving the actual calls from the GT fields below.
    for (int locus = 0; locus < loci; ++locus) {
        output << "HLS\t" << (locus + 1) << "\tlocus" << locus
               << "\tA\tG\t.\tPASS\tAF=0.5;VT=SNP\tGT";
        for (int sample = 0; sample < samples; ++sample)
            output << '\t' << (bit_at(haplotypes, stride, locus, sample) ? '1' : '0');
        output << '\n';
    }

    std::cout << "Wrote " << output_path << " from make_input(loci=" << loci
              << ", haplotypes=" << samples << ", seed=" << seed << ")\n";
    return 0;
}
