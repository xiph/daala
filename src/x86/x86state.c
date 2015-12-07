/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "x86int.h"
#include "cpu.h"
#if defined(OD_X86ASM)

#if defined(OD_GCC_INLINE_ASSEMBLY)
static void od_restore_fpu_mmx(void){
  __asm__ __volatile__("emms\n\t");
}
#endif

void od_state_opt_vtbl_init_x86(od_state *_state){
  od_state_opt_vtbl_init_c(_state);
  _state->cpu_flags=od_cpu_flags_get();
  if(_state->full_precision_references) {
    /*No 16 bit assembly as yet, but it will go here.*/
    if (_state->cpu_flags&OD_CPU_X86_SSE2) {
#if defined(OD_SSE2_INTRINSICS)
      _state->opt_vtbl.mc_predict1fmv = od_mc_predict1fmv16_sse2;
#endif
    }
  }
  else {
    /*8 bit assembly for those functions that work directly on 8-bit
      od_img and reference buffers.*/
    if (_state->cpu_flags&OD_CPU_X86_SSE2) {
#if defined(OD_GCC_INLINE_ASSEMBLY)
      _state->opt_vtbl.mc_blend_full = od_mc_blend_full8_sse2;
      _state->opt_vtbl.mc_blend_full_split = od_mc_blend_full_split8_sse2;
#endif
#if defined(OD_SSE2_INTRINSICS)
      _state->opt_vtbl.mc_predict1fmv = od_mc_predict1fmv8_sse2;
      _state->opt_vtbl.od_copy_nxn[4] = od_copy_16x16_8_sse2;
      _state->opt_vtbl.od_copy_nxn[5] = od_copy_32x32_8_sse2;
      _state->opt_vtbl.od_copy_nxn[6] = od_copy_64x64_8_sse2;
#endif
    }
  }
  if (_state->cpu_flags&OD_CPU_X86_SSE2) {
#if defined(OD_SSE2_INTRINSICS)
    _state->opt_vtbl.fdct_2d[0] = od_bin_fdct4x4_sse2;
    _state->opt_vtbl.idct_2d[0] = od_bin_idct4x4_sse2;
    _state->opt_vtbl.fdct_2d[1] = od_bin_fdct8x8_sse2;
    _state->opt_vtbl.idct_2d[1] = od_bin_idct8x8_sse2;
    OD_COPY(_state->opt_vtbl.filter_dering_direction,
     OD_DERING_DIRECTION_SSE2, OD_DERINGSIZES);
    OD_COPY(_state->opt_vtbl.filter_dering_orthogonal,
     OD_DERING_ORTHOGONAL_SSE2, OD_DERINGSIZES);
#endif
#if defined(OD_SSE41_INTRINSICS)
    if (_state->cpu_flags&OD_CPU_X86_SSE4_1) {
      _state->opt_vtbl.fdct_2d[0] = od_bin_fdct4x4_sse41;
      _state->opt_vtbl.idct_2d[0] = od_bin_idct4x4_sse41;
      _state->opt_vtbl.fdct_2d[1] = od_bin_fdct8x8_sse41;
      _state->opt_vtbl.idct_2d[1] = od_bin_idct8x8_sse41;
    }
#endif
#if defined(OD_AVX2_INTRINSICS)
    if (_state->cpu_flags & OD_CPU_X86_AVX2) {
      _state->opt_vtbl.fdct_2d[1] = od_bin_fdct8x8_avx2;
      _state->opt_vtbl.idct_2d[1] = od_bin_idct8x8_avx2;
    }
#endif
  }
#if defined(OD_GCC_INLINE_ASSEMBLY)
  if (_state->cpu_flags&OD_CPU_X86_MMX) {
    _state->opt_vtbl.restore_fpu=od_restore_fpu_mmx;
  }
#endif
}

#endif
