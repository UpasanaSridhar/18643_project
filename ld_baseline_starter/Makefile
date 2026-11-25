CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wno-unknown-pragmas
LD_PARAMS ?=
OPENCL_ROOT ?= /usr
.PHONY: all test sweep quickld-generated clean FORCE
all: test_hls

test_hls: src/ld.cpp src/test_hls.cpp include/ld.hpp hls/ld_top.cpp hls/ld_microkernel.hpp hls/ld_config.hpp FORCE
	$(CXX) $(CXXFLAGS) $(LD_PARAMS) -Iinclude src/ld.cpp src/test_hls.cpp hls/ld_top.cpp -o $@

test_fpga: src/ld.cpp src/test_fpga.cpp src/test_quickLD.cpp include/test_quickLD.hpp include/ld.hpp include/vcf_loader.hpp hls/ld_config.hpp FORCE
	$(CXX) $(CXXFLAGS) $(LD_PARAMS) -Iinclude -I$(OPENCL_ROOT)/include src/ld.cpp src/test_fpga.cpp src/test_quickLD.cpp -L$(OPENCL_ROOT)/lib -lOpenCL -lz -o $@

test_fpga_inline: src/test_fpga_inline.cpp src/ld.cpp include/ld.hpp hls/ld_config.hpp FORCE
	$(CXX) $(CXXFLAGS) $(LD_PARAMS) -Iinclude -I$(OPENCL_ROOT)/include src/test_fpga_inline.cpp src/ld.cpp -L$(OPENCL_ROOT)/lib -lOpenCL -o $@

generate_hls_vcf: src/ld.cpp src/generate_hls_vcf.cpp include/ld.hpp hls/ld_config.hpp FORCE
	$(CXX) $(CXXFLAGS) $(LD_PARAMS) -Iinclude src/ld.cpp src/generate_hls_vcf.cpp -o $@

FORCE:
test: test_hls
	./test_hls 10 10 17
	./test_hls 10 65 17
	./test_hls 11 128 23

quickld-generated:
	bash scripts/run_quickld_generated.sh

sweep:
	@set -e; for shape in '2 2' '4 4' '3 5'; do \
	  set -- $$shape; \
	  for width in 8 16 64; do \
	    $(MAKE) --no-print-directory test_hls LD_PARAMS="-DLD_MR=$$1 -DLD_NR=$$2 -DLD_MAX_K=128 -DLD_POP_WIDTH=$$width"; \
	    ./test_hls 10 10 17; ./test_hls 10 65 17; ./test_hls 11 128 23; \
	  done; \
	done
clean:
	rm -f test_hls test_quickLD test_fpga test_fpga_inline generate_hls_vcf
