// OpenCL host driver for the Vitis HLS ld_top kernel.
// Build with: make test_fpga OPENCL_ROOT=/path/to/opencl
// Run with:   ./test_fpga kernel.xclbin ld_top [vcf] [device-index]
#include "ld.hpp"
#include "vcf_loader.hpp"
#include "../hls/ld_config.hpp"

#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: test_fpga <kernel.xclbin> <kernel-name> [vcf] [device-index]\n";
        return 2;
    }
    const std::string xclbin = argv[1];
    const std::string kernel_name = argv[2];
    const std::string vcf = argc > 3 ? argv[3] : "quickLD/test_dataset/vcf/test_covid.vcf.gz";
    const unsigned int device_index = argc > 4 ? static_cast<unsigned int>(std::strtoul(argv[4], nullptr, 10)) : 0;

    try {
        const auto input = load_haploid_vcf(vcf, ld_cfg::MAX_LOCI, ld_cfg::MAX_K);
        std::cout << "Loaded " << input.loci << " loci and " << input.samples << " samples\n" << std::flush;
        const std::size_t stride = packed_stride(input.samples);
        std::vector<uint64_t> hls_input(input.loci * ld_cfg::WORDS, 0);
        for (std::size_t locus = 0; locus < input.loci; ++locus)
            for (std::size_t word = 0; word < stride; ++word)
                hls_input[locus * ld_cfg::WORDS + word] = input.packed[locus * stride + word];
        std::vector<uint32_t> counts(input.loci * input.loci, 0);

        auto check = [](cl_int status) { if (status != CL_SUCCESS) throw std::runtime_error("OpenCL error " + std::to_string(status)); };
        cl_uint nplatforms = 0;
        check(clGetPlatformIDs(0, nullptr, &nplatforms));
        std::vector<cl_platform_id> platforms(nplatforms);
        check(clGetPlatformIDs(nplatforms, platforms.data(), nullptr));
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
        if (device_index >= devices.size()) throw std::runtime_error("Xilinx accelerator device not found");
        const auto device = devices[device_index];
        std::ifstream binary_file(xclbin, std::ios::binary);
        if (!binary_file) throw std::runtime_error("cannot open xclbin: " + xclbin);
        std::vector<unsigned char> binary((std::istreambuf_iterator<char>(binary_file)), {});
        if (binary.empty()) throw std::runtime_error("empty xclbin");
        cl_int status = 0, binary_status = 0;
        cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status); check(status);
        cl_command_queue queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status); check(status);
        const unsigned char* binary_data = binary.data();
        const size_t binary_size = binary.size();
        cl_program program = clCreateProgramWithBinary(context, 1, &device, &binary_size, &binary_data, &binary_status, &status); check(status); check(binary_status);
        check(clBuildProgram(program, 1, &device, nullptr, nullptr, nullptr));
        cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &status); check(status);
        cl_mem in = clCreateBuffer(context, CL_MEM_READ_ONLY, hls_input.size() * sizeof(uint64_t), nullptr, &status); check(status);
        cl_mem out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, counts.size() * sizeof(uint32_t), nullptr, &status); check(status);
        int loci = static_cast<int>(input.loci), samples = static_cast<int>(input.samples);
        check(clSetKernelArg(kernel, 0, sizeof(in), &in));
        check(clSetKernelArg(kernel, 1, sizeof(out), &out));
        check(clSetKernelArg(kernel, 2, sizeof(loci), &loci));
        check(clSetKernelArg(kernel, 3, sizeof(samples), &samples));
        check(clEnqueueWriteBuffer(queue, in, CL_TRUE, 0, hls_input.size() * sizeof(uint64_t), hls_input.data(), 0, nullptr, nullptr));
        check(clEnqueueTask(queue, kernel, 0, nullptr, nullptr));
        check(clEnqueueReadBuffer(queue, out, CL_TRUE, 0, counts.size() * sizeof(uint32_t), counts.data(), 0, nullptr, nullptr));
        clReleaseMemObject(out); clReleaseMemObject(in); clReleaseKernel(kernel);
        clReleaseProgram(program); clReleaseCommandQueue(queue); clReleaseContext(context);

        const auto reference = gemm(input.packed, input.packed, input.loci, input.loci, input.samples);
        for (std::size_t index = 0; index < reference.size(); ++index)
            if (counts[index] != reference[index])
                throw std::runtime_error("FPGA count mismatch at matrix element " + std::to_string(index));
        std::cout << "PASS FPGA: " << input.loci << " loci x " << input.samples
                  << " haploid samples, xclbin=" << xclbin << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
