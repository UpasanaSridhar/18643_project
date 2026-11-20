#pragma once
// Compile-time HLS microkernel parameters; override with -D flags.
#ifndef LD_MR
#define LD_MR 4
#endif
#ifndef LD_NR
#define LD_NR 4
#endif
#ifndef LD_MAX_K
#define LD_MAX_K 128
#endif
#ifndef LD_MAX_LOCI
#define LD_MAX_LOCI 1024
#endif
#ifndef LD_POP_WIDTH
#define LD_POP_WIDTH 64
#endif
#if LD_MR < 1 || LD_NR < 1 || LD_MAX_K < 1 || LD_MAX_LOCI < 1
#error "LD_MR, LD_NR, LD_MAX_K and LD_MAX_LOCI must be positive"
#endif
#if LD_POP_WIDTH < 1 || LD_POP_WIDTH > 64 || (64 % LD_POP_WIDTH) != 0
#error "LD_POP_WIDTH must divide 64 and be in [1,64]"
#endif
namespace ld_cfg {
constexpr int MR = LD_MR, NR = LD_NR, MAX_K = LD_MAX_K, MAX_LOCI = LD_MAX_LOCI;
constexpr int WORDS = (MAX_K + 63) / 64;
constexpr int POP_WIDTH = LD_POP_WIDTH;
}
