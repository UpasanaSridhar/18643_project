# Usage:
# vitis_hls -f hls/run_hls.tcl
# Override via environment: LD_PART, LD_PERIOD, LD_MR, LD_NR, LD_MAX_K, LD_MAX_LOCI, LD_POP_WIDTH
# Set LD_PART to your actual Ultra96-V2 device part in your installation.
set part $::env(LD_PART)
set period 10
if {[info exists ::env(LD_PERIOD)]} {set period $::env(LD_PERIOD)}
set defs ""
foreach key {LD_MR LD_NR LD_MAX_K LD_MAX_LOCI LD_POP_WIDTH} {
    if {[info exists ::env($key)]} { append defs " -D$key=$::env($key)" }
}
open_project -reset ld_hls_project
set_top ld_top
add_files hls/ld_top.cpp -cflags "-std=c++17$defs"
open_solution -reset solution1
set_part $part
create_clock -period $period -name default
csynth_design
