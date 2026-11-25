// OpenCL counterpart of test_hls, using the same generated host-side input.
#include "ld.hpp"
#include "../hls/ld_config.hpp"
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2 || argc > 7) {
        std::cerr << "Usage: test_fpga_inline <kernel.xclbin> [loci] [samples] [seed] [kernel-name] [device-index]\n";
        return 2;
    }
    try {
        const int loci = argc > 2 ? std::atoi(argv[2]) : 10;
        const int samples = argc > 3 ? std::atoi(argv[3]) : 10;
        const unsigned seed = argc > 4 ? static_cast<unsigned>(std::strtoul(argv[4], nullptr, 10)) : 17;
        if (loci < 1 || loci > ld_cfg::MAX_LOCI || samples < 1 || samples > ld_cfg::MAX_K)
            throw std::runtime_error("dimensions exceed compiled HLS limits");
        const auto h = make_input(loci, samples, seed);
        const auto expected = gemm(h, h, loci, loci, samples);
        std::vector<uint64_t> input(loci * ld_cfg::WORDS, 0);
        const auto stride = packed_stride(samples);
        for (int i = 0; i < loci; ++i)
            for (std::size_t w = 0; w < stride; ++w) input[i * ld_cfg::WORDS + w] = h[i * stride + w];
        std::vector<uint32_t> counts(loci * loci);
        auto check = [](cl_int status) { if (status != CL_SUCCESS) throw std::runtime_error("OpenCL error " + std::to_string(status)); };
        cl_uint n = 0;
        check(clGetPlatformIDs(0, nullptr, &n));
        std::vector<cl_platform_id> platforms(n);
        check(clGetPlatformIDs(n, platforms.data(), nullptr));
        std::vector<cl_device_id> devices;
        for (auto platform : platforms) {
            char name[256] = {};
            check(clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(name), name, nullptr));
            if (std::string(name) != "Xilinx") continue;
            cl_uint count = 0;
            if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_ACCELERATOR, 0, nullptr, &count) == CL_SUCCESS) {
                devices.resize(count);
                check(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ACCELERATOR, count, devices.data(), nullptr));
                break;
            }
        }
        const auto index = argc > 6 ? std::strtoul(argv[6], nullptr, 10) : 0;
        if (index >= devices.size()) throw std::runtime_error("Xilinx accelerator device not found");
        const auto device = devices[index];
        std::ifstream file(argv[1], std::ios::binary);
        if (!file) throw std::runtime_error("cannot open xclbin");
        std::vector<unsigned char> binary((std::istreambuf_iterator<char>(file)), {});
        if (binary.empty()) throw std::runtime_error("empty xclbin");
        cl_int status = 0, binary_status = 0;
        cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status); check(status);
        cl_command_queue queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status); check(status);
        const unsigned char* data = binary.data();
        const size_t size = binary.size();
        cl_program program = clCreateProgramWithBinary(context, 1, &device, &size, &data, &binary_status, &status); check(status); check(binary_status);
        check(clBuildProgram(program, 1, &device, nullptr, nullptr, nullptr));
        cl_kernel kernel = clCreateKernel(program, argc > 5 ? argv[5] : "ld_top", &status); check(status);
        cl_mem in = clCreateBuffer(context, CL_MEM_READ_ONLY, input.size() * sizeof(uint64_t), nullptr, &status); check(status);
        cl_mem out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, counts.size() * sizeof(uint32_t), nullptr, &status); check(status);
        check(clSetKernelArg(kernel, 0, sizeof(in), &in));
        check(clSetKernelArg(kernel, 1, sizeof(out), &out));
        check(clSetKernelArg(kernel, 2, sizeof(loci), &loci));
        check(clSetKernelArg(kernel, 3, sizeof(samples), &samples));
        check(clEnqueueWriteBuffer(queue, in, CL_TRUE, 0, input.size() * sizeof(uint64_t), input.data(), 0, nullptr, nullptr));
        check(clEnqueueTask(queue, kernel, 0, nullptr, nullptr));
        check(clEnqueueReadBuffer(queue, out, CL_TRUE, 0, counts.size() * sizeof(uint32_t), counts.data(), 0, nullptr, nullptr));
        int failures = 0;
        for (int i = 0; i < loci; ++i) for (int j = 0; j < loci; ++j) {
            const int at = i * loci + j;
            if (counts[at] != expected[at]) {
                std::cerr << "MISMATCH (" << i << ',' << j << ") ref=" << expected[at] << " fpga=" << counts[at] << '\n';
                ++failures;
            }
            const double r0 = r_squared(expected[at], expected[i*loci+i], expected[j*loci+j], samples);
            const double r1 = r_squared(counts[at], counts[i*loci+i], counts[j*loci+j], samples);
            if (!((std::isnan(r0) && std::isnan(r1)) || std::abs(r0-r1) < 1e-12)) ++failures;
        }
        double total = 0, fastest = std::numeric_limits<double>::infinity(), slowest = 0;
        for (int trial = 0; trial < 100; ++trial) {
            const auto start = std::chrono::steady_clock::now();
            for (int run = 0; run < 100; ++run) check(clEnqueueTask(queue, kernel, 0, nullptr, nullptr));
            check(clFinish(queue));
            const double elapsed = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count();
            total += elapsed; fastest = std::min(fastest, elapsed); slowest = std::max(slowest, elapsed);
        }
        clReleaseMemObject(out); clReleaseMemObject(in); clReleaseKernel(kernel);
        clReleaseProgram(program); clReleaseCommandQueue(queue); clReleaseContext(context);
        std::cout << (failures ? "FAIL" : "PASS") << " mode=popcount device=OpenCL"
                  << " MR=" << ld_cfg::MR << " NR=" << ld_cfg::NR << " MAX_K=" << ld_cfg::MAX_K
                  << " POP_WIDTH=" << ld_cfg::POP_WIDTH << " loci=" << loci << " haplotypes=" << samples
                  << " tiles=" << ((loci+ld_cfg::MR-1)/ld_cfg::MR)*((loci+ld_cfg::NR-1)/ld_cfg::NR)
                  << " runtime_us_per_trial=" << std::fixed << std::setprecision(3) << (total-fastest-slowest)/98
                  << " runtime_trials=100 runs_per_trial=100 trim=min,max mismatches=" << failures << '\n';
        if (loci <= 10) {
            std::cout << "Intersection counts:\n";
            for (int i = 0; i < loci; ++i) { for (int j = 0; j < loci; ++j) std::cout << std::setw(4) << counts[i*loci+j]; std::cout << '\n'; }
        }
        return failures ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
