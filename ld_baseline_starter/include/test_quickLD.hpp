#pragma once
#include "vcf_loader.hpp"
#include <filesystem>

void run_quickld(const vcf_matrix& input, const std::filesystem::path& output);
