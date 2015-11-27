/*Daala video codec
Copyright (c) 2006-2015 Daala project contributors.  All rights reserved.

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

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "x86int.h"
#include "cpu.h"
#include "../mc.h"
#include "../util.h"

/*Block copy functions. Copying of overlapping regions has undefined
   behavior. Only 16x16, 32x32, 64x64 show a SIMD improvement.*/

#if defined(OD_X86ASM)
#define OD_IM_LOAD_1(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_1\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regb "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regc "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_1(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_1\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regb ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regc ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \
  "movdqu " _regd ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_16x16_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  __asm__ __volatile__(
    OD_IM_LOAD_1 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_LOAD_1 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_STORE_1("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_STORE_1("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_LOAD_1 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_LOAD_1 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_STORE_1("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_STORE_1("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    :[dst]"+r"(_dst),[src]"+r"(_src)
    :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
    :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
  );
}

#define OD_IM_LOAD_2(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_2\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "movdqu 16(%[src]), " _regb "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regc "\n\t" \
  "movdqu 16(%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_2(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_1\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "movdqu " _regb ",16(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regc ",(%[dst]) \n\t" \
  "movdqu " _regd ",16(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_32x32_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  int i;
  for (i = 0; i < 4; i++) {
    __asm__ __volatile__(
      OD_IM_LOAD_2 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_2 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_2("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_2("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_LOAD_2 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_2 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_2("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_2("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      :[dst]"+r"(_dst),[src]"+r"(_src)
      :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
      :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
    );
  }
}

#define OD_IM_LOAD_4(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_4\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "movdqu 16(%[src]), " _regb "\n\t" \
  "movdqu 32(%[src]), " _regc "\n\t" \
  "movdqu 48(%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_4(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_4\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "movdqu " _regb ",16(%[dst]) \n\t" \
  "movdqu " _regc ",32(%[dst]) \n\t" \
  "movdqu " _regd ",48(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_64x64_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  int i;
  for (i = 0; i < 16; i++) {
    __asm__ __volatile__(
      OD_IM_LOAD_4 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_4 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_4("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_4("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_LOAD_4 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_4 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_4("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_4("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      :[dst]"+r"(_dst),[src]"+r"(_src)
      :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
      :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
    );
  }
}

#endif
