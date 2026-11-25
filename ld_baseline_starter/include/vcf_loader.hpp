#pragma once

#include "ld.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

struct vcf_matrix {
    std::string chrom;
    std::vector<std::string> ids;
    std::size_t loci = 0;
    std::size_t samples = 0;
    std::vector<uint64_t> packed;
};

inline std::vector<std::string> vcf_split(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const std::size_t end = line.find('\t', start);
        fields.push_back(line.substr(start, end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return fields;
}

// Loads a VCF/VCF.GZ wholly into host memory, then retains a binary, haploid
// prefix that fits the compiled HLS limits.  A GT of 1 is the alternate allele;
// 0 is the reference allele. Missing and non-biallelic calls are rejected so
// the input seen by quickLD and the hardware path is unambiguous.
inline vcf_matrix load_haploid_vcf(const std::string& path,
                                   std::size_t max_loci,
                                   std::size_t max_samples) {
    gzFile file = gzopen(path.c_str(), "rb");
    if (!file) throw std::runtime_error("could not open VCF: " + path);

    std::string text;
    char buffer[65536];
    int bytes = 0;
    while ((bytes = gzread(file, buffer, sizeof(buffer))) > 0)
        text.append(buffer, static_cast<std::size_t>(bytes));
    if (bytes < 0) {
        int error_number = Z_OK;
        const char* error = gzerror(file, &error_number);
        gzclose(file);
        throw std::runtime_error(std::string("could not read VCF: ") + error);
    }
    gzclose(file);

    vcf_matrix result;
    std::size_t input_samples = 0;
    std::size_t position = 0;
    while (position < text.size()) {
        const std::size_t end = text.find('\n', position);
        std::string line = text.substr(position, end - position);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        position = end == std::string::npos ? text.size() : end + 1;
        if (line.empty() || (line[0] == '#' && line.rfind("#CHROM", 0) != 0)) continue;

        const auto fields = vcf_split(line);
        if (line.rfind("#CHROM", 0) == 0) {
            if (fields.size() < 10) throw std::runtime_error("VCF has no samples");
            input_samples = fields.size() - 9;
            result.samples = std::min(input_samples, max_samples);
            if (!result.samples) throw std::runtime_error("sample limit is zero");
            continue;
        }
        if (!input_samples) throw std::runtime_error("VCF data precedes #CHROM header");
        if (result.loci == max_loci) break;
        // quickLD's bundled data may contain records that are not representable
        // by this binary HLS kernel. Keep scanning and retain only complete,
        // biallelic haploid GT records.
        if (fields.size() != input_samples + 9 || fields[4].find(',') != std::string::npos ||
            fields[8] != "GT") continue;
        if (result.chrom.empty()) result.chrom = fields[0];
        if (fields[0] != result.chrom) continue;

        const std::size_t stride = packed_stride(result.samples);
        std::vector<uint64_t> row(stride, 0);
        bool usable = true;
        for (std::size_t sample = 0; sample < result.samples; ++sample) {
            const std::string& gt = fields[9 + sample];
            if (gt == "1") row[sample / 64] |= 1ULL << (sample % 64);
            else if (gt != "0") { usable = false; break; }
        }
        if (!usable) continue;
        result.packed.insert(result.packed.end(), row.begin(), row.end());
        result.ids.push_back(fields[2]);
        ++result.loci;
    }
    if (!result.loci) throw std::runtime_error("no VCF records were loaded");
    return result;
}
