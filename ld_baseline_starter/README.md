# Linkage disequilibrium baseline starter

Three pieces only: **quickLD** as the external CPU/GPU baseline, a **naive packed-popcount C++ reference**, and a **configurable complete-LD Vitis HLS kernel** for Ultra96-V2. The bit-serial implementation is deliberately postponed.

## Contents

- `scripts/setup_quickld.sh`: clone/build upstream quickLD BLIS baseline (requires network and toolchain). No quickLD source is included.
- `include/ld.hpp`, `src/ld.cpp`: class-free locus-major packed `uint64_t` arrays, a naive `gemm` (`AND + popcount`) returning a flat `uint64_t` array, and `r²`.
- `hls/ld_config.hpp`: compile-time complete-LD and microkernel parameters.
- `hls/ld_microkernel.hpp`, `hls/ld_top.cpp`: tiled popcount reduction and an HLS AXI top that calculates the complete LD matrix.
- `src/test_hls.cpp`: deterministic configurable native C++ correctness test, including edge tiles and partial 64-bit words.

### Native test (x86 or Arm)

```sh
make test
./test_hls 10 10 17   # loci, haplotypes, seed
./test_hls 10 65 17
make sweep          # MR×NR = 2×2, 4×4, 3×5; popcount width = 8,16,64
```

Change the compiled microkernel dimensions:

```sh
make test_hls LD_PARAMS='-DLD_MR=2 -DLD_NR=5 -DLD_MAX_K=128 -DLD_POP_WIDTH=16'
./test_hls 10 65 17
```

The native test compares the complete HLS LD result's exact intersection counts against a naive `__builtin_popcountll(A & B)` reference for every pair. It checks `r²` as a derived value. Monomorphic loci have undefined `r²` (NaN), but defined intersection counts.

### Microkernel parameters

| Macro | Default | Meaning |
|---|---:|---|
| `LD_MR` | 4 | Number of loci in output tile rows |
| `LD_NR` | 4 | Number of loci in output tile columns |
| `LD_MAX_K` | 128 | Maximum haplotypes per kernel call; fixes panel stride |
| `LD_MAX_LOCI` | 1024 | Maximum loci handled by one complete-LD invocation |
| `LD_POP_WIDTH` | 64 | Bits counted each reduction step; must divide 64 |

`ld_top(haplotypes, counts, loci, samples)` accepts a pre-packed, locus-major input buffer of shape `loci × WORDS`, where `WORDS = ceil(LD_MAX_K/64)`. The caller must zero unused words and unused bits in the final active word before calling the kernel. It produces a row-major `loci × loci` intersection-count matrix. The HLS top-level performs every `MR × NR` tile iteration and handles only its local final-tile padding; the host makes exactly one kernel call. The kernel calculates integer intersection counts only, not `r²`.

### Vitis HLS (Ultra96-V2)

```sh
export LD_PART=xczu3eg-sbva484-1-e
export LD_MR=4 LD_NR=4 LD_MAX_K=128 LD_POP_WIDTH=64
vitis_hls -f hls/run_hls.tcl
```

Set `LD_PART` for the actual device/toolchain. This is **not synthesized or benchmarked yet**. A `PIPELINE II=1` pragma is a request, not a measured II. Popcount is expressed as unrolled bit summation, and **DSP mapping is not guaranteed**. Use utilization and timing reports to evaluate DSP/LUT mapping.

### quickLD (separate baseline)

```sh
bash scripts/setup_quickld.sh
```

On Orin: build CPU BLIS and check the NVIDIA OpenCL runtime before using quickLD's `-gpu` option. On Ultra96-V2 Arm: run the BLIS CPU version. For phased 1000 Genomes VCF, quickLD's documented pipeline uses `qLD-parseVCF`, `qLD-parse-2MDF` with `-ploidy phased_diploid`, then `qLD-compute` with `-mdf -blis` (or `-gpu`). Check the CLI in the installed version.

### quickLD on the HLS test fixture

`generate_hls_vcf` uses the same `make_input(loci, haplotypes, seed)` call as
`test_hls`, and writes a haploid VCF where each VCF sample represents one HLS
haplotype bit. This preserves the generated binary matrix exactly, including an
odd number of haplotypes. `quickLD` reports `r²`, whereas the HLS kernel reports
intersection counts; the all-zero and all-one edge loci are correctly omitted by
quickLD because their `r²` is undefined.

After installing quickLD's BLIS build, run:

```bash
make quickld-generated
# or: bash scripts/run_quickld_generated.sh 10 65 17 output_directory
```

The script creates `hls_generated.vcf`, parses it to MDF, runs the BLIS CPU
backend with `-r2limit 0`, and writes `quickld_report.txt` in the chosen output
directory. Reuse a fresh output directory on repeated runs: upstream quickLD
renames existing parsed/MDF directories rather than overwriting them.
