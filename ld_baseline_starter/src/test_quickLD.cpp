#include "ld.hpp"
#include "vcf_loader.hpp"
#include "../hls/ld_config.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) quoted += c == '\'' ? "'\\''" : std::string(1, c);
    return quoted + "'";
}
} // namespace

void run_quickld(const vcf_matrix& input, const std::filesystem::path& output) {
    const std::filesystem::path quickld = "quickLD/bin";
    for (const char* executable : {"qLD-parse-VCF", "qLD-parse-2MDF", "qLD-compute"})
        if (!std::filesystem::is_regular_file(quickld / executable))
            throw std::runtime_error("quickLD is not built; run scripts/setup_quickld.sh first");
    if (std::filesystem::exists(output))
        throw std::runtime_error("output directory already exists: " + output.string());
    std::filesystem::create_directories(output);

    const auto vcf = output / "input.vcf";
    std::ofstream stream(vcf);
    stream << "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (std::size_t sample = 0; sample < input.samples; ++sample) stream << "\ts" << sample;
    stream << '\n';
    const std::size_t stride = packed_stride(input.samples);
    for (std::size_t locus = 0; locus < input.loci; ++locus) {
        stream << input.chrom << '\t' << locus + 1 << '\t' << input.ids[locus]
               << "\tA\tG\t.\tPASS\t.\tGT";
        for (std::size_t sample = 0; sample < input.samples; ++sample)
            stream << '\t' << ((input.packed[locus * stride + sample / 64] >> (sample % 64)) & 1ULL);
        stream << '\n';
    }
    stream.close();
    const auto parsed = output / "parsed";
    const auto mdf = output / "mdf";
    const auto report = output / "quickld_cpu_report.txt";
    const std::string prefix = shell_quote((quickld / "qLD-parse-VCF").string());
    // qLD-parse-VCF returns 1 after a successful parse in this quickLD release.
    const std::string parse = prefix + " -input " + shell_quote(vcf.string()) +
        " -output " + shell_quote(parsed.string()) + " -size 1 -chrom " + shell_quote(input.chrom);
    const int parse_status = std::system(parse.c_str());
    if (parse_status != 0 && parse_status != 256) throw std::runtime_error("qLD-parse-VCF failed");
    const std::string mdf_cmd = shell_quote((quickld / "qLD-parse-2MDF").string()) +
        " -input " + shell_quote(parsed.string()) + " -output " + shell_quote(mdf.string()) + " -ploidy haploid";
    const std::string compute = shell_quote((quickld / "qLD-compute").string()) +
        " -input " + shell_quote(mdf.string()) + " -output " + shell_quote(report.string()) +
        " -ploidy haploid -mdf -blis -threads 1 -r2limit 0";
    if (std::system(mdf_cmd.c_str()) != 0 || std::system(compute.c_str()) != 0)
        throw std::runtime_error("quickLD CPU pipeline failed");
    std::cout << "quickLD CPU report: " << report << '\n';
}
