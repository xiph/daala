/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "x86enc.h"
#include "x86int.h"

#if defined(OD_X86ASM)

# if defined(OD_CHECKASM)

#  include <stdio.h>

void od_mc_compute_sad8_check(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride, int w, int h, int32_t sad) {
  int32_t c_sad;
  c_sad = od_mc_compute_sad8_c(src, systride, ref, dystride, w, h);
  if (sad != c_sad) {
    fprintf(stderr, "od_mc_compute_sad %ix%i check failed: %i!=%i\n",
     w, h, sad, c_sad);
  }
  OD_ASSERT(sad == c_sad);
}
# endif

#if defined(OD_GCC_INLINE_ASSEMBLY)
/*Handle one 4x4 block with dxstride == 1.*/
int32_t od_mc_compute_sad8_4x4_sse(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride){
  ptrdiff_t srow;
  ptrdiff_t drow;
  int32_t ret;
  __asm__ __volatile__(
    "movd (%[src]), %%mm0\n"
    "movd (%[src], %[systride]), %%mm2\n"
    "punpckldq %%mm2, %%mm0\n"
    "movd (%[ref]), %%mm1\n"
    "movd (%[ref], %[dystride]), %%mm3\n"
    "punpckldq %%mm3, %%mm1\n"
    "psadbw %%mm1, %%mm0\n"
    "lea (%[src], %[systride], 2), %[srow]\n"
    "lea (%[ref], %[dystride], 2), %[drow]\n"
    : [srow]"=&r"(srow), [drow]"=r"(drow)
    : [src]"r"(src), [systride]"r"((ptrdiff_t)systride), [ref]"r"(ref),
     [dystride]"r"((ptrdiff_t)dystride)
  );
  __asm__ __volatile__(
    "movd (%[srow]), %%mm4\n"
    "movd (%[srow], %[systride]), %%mm5\n"
    "punpckldq %%mm5, %%mm4\n"
    "movd (%[drow]), %%mm6\n"
    "movd (%[drow], %[dystride]), %%mm7\n"
    "punpckldq %%mm7, %%mm6\n"
    "psadbw %%mm6, %%mm4\n"
    "paddd %%mm4, %%mm0\n"
    "movd %%mm0, %[ret]\n"
    : [ret]"=r"(ret)
    : [systride]"r"((ptrdiff_t)systride),
     [dystride]"r"((ptrdiff_t)dystride), [srow]"r"(srow), [drow]"r"(drow)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad8_check(src, systride, ref, dystride, 4, 4, ret);
# endif
  return ret;
}

/*Handle one 8x8 block with dxstride == 1.*/
int32_t od_mc_compute_sad8_8x8_sse(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride){
  const unsigned char *srow;
  const unsigned char *drow;
  int32_t ret;
  int i;
  srow = src;
  drow = ref;
  __asm__ __volatile__(
    "pxor %mm2, %mm2\n"
  );
  for (i = 0; i < 8; i++) {
    __asm__ __volatile__(
      "movq (%[srow]), %%mm0\n"
      "movq (%[drow]), %%mm1\n"
      "psadbw %%mm1, %%mm0\n"
      "paddd %%mm0, %%mm2\n"
      :
      : [srow]"r"(srow), [drow]"r"(drow)
    );
    srow += systride;
    drow += dystride;
  }
  __asm__ __volatile__(
    "movd %%mm2, %[ret]\n"
    :[ret]"=r"(ret)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad8_check(src, systride, ref, dystride, 8, 8, ret);
# endif
  return ret;
}

/*Handle one nxn block with dxstride == 1 where n is 2^ln and n >= 16.*/
OD_SIMD_INLINE int32_t od_mc_compute_sad8_nxn_sse2(const int ln,
 const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  const unsigned char *srow;
  const unsigned char *drow;
  int32_t ret;
  int n;
  int i;
  int j;
  n = 1 << ln;
  OD_ASSERT(n >= 16);
  srow = src;
  drow = ref;
  __asm__ __volatile__(
    "pxor %xmm2, %xmm2\n\t"
  );
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j += 16) {
      __asm__ __volatile__(
        "movdqu (%[drow]), %%xmm1\n"
        "movdqu (%[srow]), %%xmm0\n"
        "psadbw %%xmm1, %%xmm0\n"
        "paddq %%xmm0, %%xmm2\n"
        :
        : [srow]"r"(srow + j), [drow]"r"(drow + j)
      );
    }
    srow += systride;
    drow += dystride;
  }
  __asm__ __volatile__(
    "movdqa %%xmm2, %%xmm0\n"
    "punpckhqdq %%xmm2, %%xmm0\n"
    "paddq %%xmm2, %%xmm0\n"
    "movd %%xmm0, %[ret]\n"
    : [ret]"=r"(ret)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad8_check(src, systride, ref, dystride, n, n, ret);
# endif
  return ret;
}

/*Handle one 16x16 block with dxstride == 1.*/
int32_t od_mc_compute_sad8_16x16_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_nxn_sse2(4, src, systride, ref, dystride);
}

/*Handle one 32x32 block with dxstride == 1.*/
int32_t od_mc_compute_sad8_32x32_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_nxn_sse2(5, src, systride, ref, dystride);
}

/*Handle one 64x64 block with dxstride == 1.*/
int32_t od_mc_compute_sad8_64x64_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride) {
  return od_mc_compute_sad8_nxn_sse2(6, src, systride, ref, dystride);
}

#endif

#endif
