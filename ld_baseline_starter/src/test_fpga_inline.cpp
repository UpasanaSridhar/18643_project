// Standalone OpenCL smoke test for ld_top, using three inline haplotype rows.
#include "../hls/ld_config.hpp"
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: test_fpga_inline <kernel.xclbin> [kernel-name] [device-index]\n";
        return 2;
    }
    try {
        constexpr int loci = 3, samples = 8;
        if (loci > ld_cfg::MAX_LOCI || samples > ld_cfg::MAX_K)
            throw std::runtime_error("inline matrix exceeds compiled HLS limits");
        // Rows, in sample order: 10101010, 11001100, 11110000.
        const uint64_t rows[loci] = {0x55, 0x33, 0x0f};
        const uint32_t expected[loci * loci] = {4, 2, 2, 2, 4, 2, 2, 2, 4};
        std::vector<uint64_t> input(loci * ld_cfg::WORDS, 0);
        for (int i = 0; i < loci; ++i) input[i * ld_cfg::WORDS] = rows[i];
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
        const auto index = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 0;
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
        cl_kernel kernel = clCreateKernel(program, argc > 2 ? argv[2] : "ld_top", &status); check(status);
        cl_mem in = clCreateBuffer(context, CL_MEM_READ_ONLY, input.size() * sizeof(uint64_t), nullptr, &status); check(status);
        cl_mem out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, counts.size() * sizeof(uint32_t), nullptr, &status); check(status);
        check(clSetKernelArg(kernel, 0, sizeof(in), &in));
        check(clSetKernelArg(kernel, 1, sizeof(out), &out));
        check(clSetKernelArg(kernel, 2, sizeof(loci), &loci));
        check(clSetKernelArg(kernel, 3, sizeof(samples), &samples));
        check(clEnqueueWriteBuffer(queue, in, CL_TRUE, 0, input.size() * sizeof(uint64_t), input.data(), 0, nullptr, nullptr));
        check(clEnqueueTask(queue, kernel, 0, nullptr, nullptr));
        check(clEnqueueReadBuffer(queue, out, CL_TRUE, 0, counts.size() * sizeof(uint32_t), counts.data(), 0, nullptr, nullptr));
        clReleaseMemObject(out); clReleaseMemObject(in); clReleaseKernel(kernel);
        clReleaseProgram(program); clReleaseCommandQueue(queue); clReleaseContext(context);
        for (int i = 0; i < loci * loci; ++i) {
            if (counts[i] != expected[i]) throw std::runtime_error("mismatch at (" + std::to_string(i / loci) + "," + std::to_string(i % loci) + "): got " + std::to_string(counts[i]) + ", expected " + std::to_string(expected[i]));
        }
        std::cout << "PASS: inline 3 x 8 matrix, all nine intersection counts match\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
