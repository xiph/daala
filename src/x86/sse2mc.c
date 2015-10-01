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

#if defined(OD_X86ASM)
#include <xmmintrin.h>

/*So here are the constraints we have:
  We want (as much as possible) the same code to compile on x86-32 and x86-64.
  We want the same code to be used with and without -fomit-frame-pointer and
   with and without -fPIC.
  These are, in fact, the reason we're using inline gcc assembly instead of
   separate asm files.
  Separate files using the Intel syntax might have a chance of compiling on
   win32 without requiring mingw32, but would introduce a dependency on a
   separate assembler like nasm everywhere else, would require us to maintain
   multiple versions of the code for, e.g., 32- and 64-bit processors, and
   would require horrible macro garbage to correct for things like name
   mangling, which are not consistent between platforms.
  gcc can take care of all that crap for us, _and_ I can cross-compile for
   win32 from Linux.
  That certainly made the decision easy for me.

  These are the rules:
  We cannot use either %ebp (not available without -fomit-frame-pointer) or
   %ebx (not available with -fPIC) as general purpose registers, nor can we
   rely on gcc being able to use them when it constructs inputs/outputs for us.
  This gives the register-starved x86 architecture just 5 general purpose
   registers.
  Compile without -fomit-frame-pointer and with -fPIC during testing to make
   sure you don't break register allocation.
  To make matters worse, we can't combine static variables in "m" inputs with
   an index register or scale, because with -fPIC they might be using %ebx as a
   base register, and there's no syntactically valid way to combine the two.
  Therefore we have to load them into a register first, and then index them.
  Finally, to use scratch registers that contain pointers or pointer offsets,
   we must declare a local C variable of type ptrdiff_t and let gcc allocate it
   to a register for us.
  This will automatically select a 32- or 64-bit register as appropriate for
   the target platform.
  Type long is not good enough thanks to win64, where long is still 32 bits.*/

/*Constants used for the bilinear blending weights.*/
static const unsigned short __attribute__((aligned(16),used)) OD_BIL4H[8]={
  0x0,0x1,0x2,0x3,0x0,0x1,0x2,0x3
};
static const unsigned short __attribute__((aligned(16),used)) OD_BILH[32]={
  0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
  0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,
  0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
  0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,
};
static const unsigned short __attribute__((aligned(16),used)) OD_BIL4V[64]={
  0x0,0x0,0x0,0x0,0x1,0x1,0x1,0x1,
  0x2,0x2,0x2,0x2,0x3,0x3,0x3,0x3,
  0x4,0x4,0x4,0x4,0x5,0x5,0x5,0x5,
  0x6,0x6,0x6,0x6,0x7,0x7,0x7,0x7,
  0x8,0x8,0x8,0x8,0x9,0x9,0x9,0x9,
  0xA,0xA,0xA,0xA,0xB,0xB,0xB,0xB,
  0xC,0xC,0xC,0xC,0xD,0xD,0xD,0xD,
  0xE,0xE,0xE,0xE,0xF,0xF,0xF,0xF
};
static const unsigned short __attribute__((aligned(16),used)) OD_BILV[256]={
  0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
  0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,
  0x02,0x02,0x02,0x02,0x02,0x02,0x02,0x02,
  0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
  0x04,0x04,0x04,0x04,0x04,0x04,0x04,0x04,
  0x05,0x05,0x05,0x05,0x05,0x05,0x05,0x05,
  0x06,0x06,0x06,0x06,0x06,0x06,0x06,0x06,
  0x07,0x07,0x07,0x07,0x07,0x07,0x07,0x07,
  0x08,0x08,0x08,0x08,0x08,0x08,0x08,0x08,
  0x09,0x09,0x09,0x09,0x09,0x09,0x09,0x09,
  0x0A,0x0A,0x0A,0x0A,0x0A,0x0A,0x0A,0x0A,
  0x0B,0x0B,0x0B,0x0B,0x0B,0x0B,0x0B,0x0B,
  0x0C,0x0C,0x0C,0x0C,0x0C,0x0C,0x0C,0x0C,
  0x0D,0x0D,0x0D,0x0D,0x0D,0x0D,0x0D,0x0D,
  0x0E,0x0E,0x0E,0x0E,0x0E,0x0E,0x0E,0x0E,
  0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,0x0F,
  0x10,0x10,0x10,0x10,0x10,0x10,0x10,0x10,
  0x11,0x11,0x11,0x11,0x11,0x11,0x11,0x11,
  0x12,0x12,0x12,0x12,0x12,0x12,0x12,0x12,
  0x13,0x13,0x13,0x13,0x13,0x13,0x13,0x13,
  0x14,0x14,0x14,0x14,0x14,0x14,0x14,0x14,
  0x15,0x15,0x15,0x15,0x15,0x15,0x15,0x15,
  0x16,0x16,0x16,0x16,0x16,0x16,0x16,0x16,
  0x17,0x17,0x17,0x17,0x17,0x17,0x17,0x17,
  0x18,0x18,0x18,0x18,0x18,0x18,0x18,0x18,
  0x19,0x19,0x19,0x19,0x19,0x19,0x19,0x19,
  0x1A,0x1A,0x1A,0x1A,0x1A,0x1A,0x1A,0x1A,
  0x1B,0x1B,0x1B,0x1B,0x1B,0x1B,0x1B,0x1B,
  0x1C,0x1C,0x1C,0x1C,0x1C,0x1C,0x1C,0x1C,
  0x1D,0x1D,0x1D,0x1D,0x1D,0x1D,0x1D,0x1D,
  0x1E,0x1E,0x1E,0x1E,0x1E,0x1E,0x1E,0x1E,
  0x1F,0x1F,0x1F,0x1F,0x1F,0x1F,0x1F,0x1F
};

#if defined(OD_CHECKASM)
void od_mc_predict1fmv8_check(unsigned char *_dst,const unsigned char *_src,
 int _systride,int32_t _mvx,int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz){
  unsigned char dst[OD_MVBSIZE_MAX*OD_MVBSIZE_MAX];
  int           xblk_sz;
  int           yblk_sz;
  int           failed;
  int           i;
  int           j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_predict1fmv8_c(dst,_src,_systride,_mvx,_mvy,
   _log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if(_dst[i+(j<<_log_xblk_sz)]!=dst[i+(j<<_log_xblk_sz)]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         _dst[i+(j<<_log_xblk_sz)],dst[i+(j<<_log_xblk_sz)],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_predict1fmv8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif

/*Fills 3 vectors with pairs of alternating 16 bit values for the 1D filter
   chosen for the fractional position of x or y mv.*/
OD_SIMD_INLINE void od_setup_alternating_filter_variables(
 __m128i *filter_01, __m128i *filter_23, __m128i *filter_45, int mvf) {
  uint32_t* f;
  /*Load 3 pairs of 16 bit values as 32 bit values.
    Fill each filter with these 32 values as to create a vector that
     alternates between the 2 16 bit values.*/
  f = (uint32_t *)(OD_SUBPEL_FILTER_SET[mvf]);
  *filter_01 = _mm_set1_epi32(f[0]);
  *filter_23 = _mm_set1_epi32(f[1]);
  *filter_45 = _mm_set1_epi32(f[2]);
}

OD_SIMD_INLINE __m128i od_mc_multiply_reduce_add_horizontal_4(
 __m128i src_vec, __m128i fx01, __m128i fx23, __m128i fx45) {
  __m128i src8pels;
  __m128i sums;
  __m128i madd01;
  __m128i madd23;
  __m128i madd45;
  /*Create a pattern of 0,1, 1,2, 2,3 ... 7,8*/
  src_vec = _mm_unpacklo_epi8(src_vec, _mm_srli_si128(src_vec, 1));
  /*Unpack src_vec from 8 bit unsigned integers to 16 bit integers.
    Multiply by each set of filters and then add the two horizontally adjacent
     products together.
    This results in 4 32 bit integers.
    Perform these operations for each pair of filter values.*/
  src8pels = _mm_unpacklo_epi8(src_vec, _mm_setzero_si128());
  madd01 = _mm_madd_epi16(fx01, src8pels);
  src_vec = _mm_srli_si128(src_vec, 4);
  src8pels = _mm_unpacklo_epi8(src_vec, _mm_setzero_si128());
  madd23 = _mm_madd_epi16(fx23, src8pels);
  src_vec = _mm_srli_si128(src_vec, 4);
  src8pels = _mm_unpacklo_epi8(src_vec, _mm_setzero_si128());
  madd45 = _mm_madd_epi16(fx45, src8pels);
  /*Subtract from one of the summands instead of the final value to avoid
    data hazards.*/
  madd45 = _mm_sub_epi32(madd45, _mm_set1_epi32(128 << OD_SUBPEL_COEFF_SCALE));
  /*Sum together the 3 summands.*/
  sums = _mm_add_epi32(madd01, madd23);
  sums = _mm_add_epi32(sums, madd45);
  /*Subtraction would occur here if it wasn't performed earlier.*/
  sums = _mm_packs_epi32(sums, sums);
  return sums;
}

OD_SIMD_INLINE void od_mc_predict1fmv8_horizontal_nxm(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf,
 const int xblk_sz, const int yblk_sz) {
  int i;
  int j;
  if (mvxf) {
    __m128i fx01;
    __m128i fx23;
    __m128i fx45;
    od_setup_alternating_filter_variables(&fx01, &fx23, &fx45, mvxf);
    j = -OD_SUBPEL_TOP_APRON_SZ;
    /*The mvy is of integer position*/
    if (!mvyf) {
      /*Change j such that the loop is done yblk_sz times.*/
      j = OD_SUBPEL_TOP_APRON_SZ;
      buff_p += xblk_sz*OD_SUBPEL_TOP_APRON_SZ;
      src_p += systride*OD_SUBPEL_TOP_APRON_SZ;
    }
    for (; j < yblk_sz + OD_SUBPEL_BOTTOM_APRON_SZ; j++) {
      for (i = 0; i < xblk_sz; i += 4) {
        __m128i tmp;
        __m128i sums;
        tmp = _mm_loadu_si128((__m128i *)(src_p + i - OD_SUBPEL_TOP_APRON_SZ));
        sums = od_mc_multiply_reduce_add_horizontal_4(tmp, fx01, fx23, fx45);
        /*Only store as many values as xblk_sz.*/
        if(xblk_sz >= 4) {
          OD_ASSERT(i + 4 <= xblk_sz);
          _mm_storel_epi64((__m128i *) (buff_p + i), sums);
        }
        else {
          OD_ASSERT(i + 2 <= xblk_sz);
          *((uint32_t *)(buff_p + i)) = (uint32_t)_mm_cvtsi128_si32(sums);
        }
      }
      src_p += systride;
      buff_p += xblk_sz;
    }
  }
  /*The mvx is of integer position.*/
  else {
    __m128i normalize_128;
    normalize_128 = _mm_set1_epi16(128);
    for (j = -OD_SUBPEL_TOP_APRON_SZ;
     j < yblk_sz + OD_SUBPEL_BOTTOM_APRON_SZ; j++) {
      for (i = 0; i < xblk_sz; i += 8) {
        __m128i tmp;
        __m128i src8pels;
        tmp = _mm_loadl_epi64((__m128i *)(src_p + i));
        src8pels = _mm_unpacklo_epi8(tmp, _mm_setzero_si128());
        src8pels = _mm_slli_epi16(
         _mm_subs_epi16(src8pels, normalize_128),
         OD_SUBPEL_COEFF_SCALE);
        /*Only store as many values as xblk_sz.*/
        if (xblk_sz >= 8)  {
          _mm_store_si128((__m128i *)(buff_p + i), src8pels);
        }
        else if (xblk_sz >= 4) {
          OD_ASSERT(i + 4 <= xblk_sz);
          _mm_storel_epi64((__m128i *)(buff_p + i), src8pels);
        }
        else {
          OD_ASSERT(i + 2 <= xblk_sz);
          *((uint32_t *)(buff_p + i)) = (uint32_t)_mm_cvtsi128_si32(src8pels);
        }
      }
      src_p += systride;
      buff_p += xblk_sz;
    }
  }
}

void od_mc_predict1fmv8_horizontal_2x2(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf) {
  od_mc_predict1fmv8_horizontal_nxm(buff_p, src_p, systride, mvxf, mvyf, 2, 2);
}

void od_mc_predict1fmv8_horizontal_4x4(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf) {
  od_mc_predict1fmv8_horizontal_nxm(buff_p, src_p, systride, mvxf, mvyf, 4, 4);
}

void od_mc_predict1fmv8_horizontal_8x8(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf) {
  od_mc_predict1fmv8_horizontal_nxm(buff_p, src_p, systride, mvxf, mvyf, 8, 8);
}

void od_mc_predict1fmv8_horizontal_16x16(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf) {
  od_mc_predict1fmv8_horizontal_nxm(buff_p, src_p, systride, mvxf, mvyf,
   16, 16);
}

void od_mc_predict1fmv8_horizontal_32x32(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf) {
  od_mc_predict1fmv8_horizontal_nxm(buff_p, src_p, systride, mvxf, mvyf,
   32, 32);
}

typedef void (*od_mc_predict1fmv8_horizontal_fixed_func)(int16_t *buff_p,
 const unsigned char *src_p, int systride, int mvxf, int mvyf);

#if defined(OD_SSE2_INTRINSICS)
void od_mc_predict1fmv8_sse2(unsigned char *dst,const unsigned char *src,
 int systride, int32_t mvx, int32_t mvy,
 int log_xblk_sz, int log_yblk_sz) {
  static const od_mc_predict1fmv8_horizontal_fixed_func VTBL_HORIZONTAL[5] = {
    od_mc_predict1fmv8_horizontal_2x2, od_mc_predict1fmv8_horizontal_4x4,
    od_mc_predict1fmv8_horizontal_8x8, od_mc_predict1fmv8_horizontal_16x16,
    od_mc_predict1fmv8_horizontal_32x32
  };
  int mvxf;
  int mvyf;
  int xblk_sz;
  int yblk_sz;
  int i;
  int j;
  /*Pointer to the start of an image block in local buffer (defined
     below, buff[]), where the buffer contains the top and bottom apron
     area of the image block.
    Used as output for 1st stage horizontal filtering then as input for
     2nd stage vertical filtering.*/
  int16_t *buff_p;
  /*A pointer to input row for both 1st and 2nd stage filtering.*/
  const unsigned char *src_p;
  unsigned char *dst_p;
  /*2D buffer to store the result of 1st stage (i.e. horizontal) 1D filtering
     of a block. The 1st stage filtering requires to output results for
     top and bottom aprons of input image block, because the 2nd stage
     filtering (i.e vertical) requires support region on those apron pixels.
    The size of the buffer is :
     wxh = OD_MVBSIZE_MAX x (OD_MVBSIZE_MAX + BUFF_APRON_SZ).*/
  int16_t buff[(OD_MVBSIZE_MAX + OD_SUBPEL_BUFF_APRON_SZ)
   *OD_MVBSIZE_MAX + 16];
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  src_p = src + (mvx >> 3) + (mvy >> 3)*systride;
  dst_p = dst;
  /*Fetch LSB 3 bits, i.e. fractional MV.*/
  mvxf = mvx & 0x07;
  mvyf = mvy & 0x07;
  /*Check whether mvxf and mvyf are in the range [0...7],
     i.e. downto 1/8 precision.*/
  OD_ASSERT(mvxf <= 7);
  OD_ASSERT(mvyf <= 7);
  /*MC with subpel MV?*/
  if (mvxf || mvyf) {
    /*1st stage 1D filtering, Horizontal.*/
    buff_p = buff;
    src_p -= systride*OD_SUBPEL_TOP_APRON_SZ;
    OD_ASSERT(log_xblk_sz == log_yblk_sz);
    (*VTBL_HORIZONTAL[log_xblk_sz - 1])(buff_p, src_p, systride, mvxf, mvyf);
    /*2nd stage 1D filtering, Vertical.*/
    buff_p = buff + xblk_sz*OD_SUBPEL_TOP_APRON_SZ;
    if (mvyf)
    {
      __m128i fy01;
      __m128i fy23;
      __m128i fy45;
      __m128i rounding_offset;
      od_setup_alternating_filter_variables(&fy01, &fy23, &fy45, mvyf);
      rounding_offset = _mm_set1_epi32(OD_SUBPEL_RND_OFFSET3);
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i += 8) {
          /*Keeps 8 shorts from each row.*/
          __m128i row0_src;
          __m128i row1_src;
          __m128i row2_src;
          __m128i row3_src;
          __m128i row4_src;
          __m128i row5_src;
          __m128i sums32_0to3;
          __m128i sums32_4to7;
          __m128i row01_lo;
          __m128i row01_hi;
          __m128i out;
          OD_ASSERT((buff_p + i + ((0 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          OD_ASSERT((buff_p + i + ((1 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          OD_ASSERT((buff_p + i + ((2 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          OD_ASSERT((buff_p + i + ((3 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          OD_ASSERT((buff_p + i + ((4 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          OD_ASSERT((buff_p + i + ((5 - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz) + 15)
           < buff + sizeof(buff));
          /*Load input coeffs from each row, 8 shorts at one time.*/
          row0_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((0 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          row1_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((1 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          row2_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((2 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          row3_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((3 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          row4_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((4 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          row5_src = _mm_loadu_si128((__m128i *)
           (buff_p + i + ((5 - OD_SUBPEL_TOP_APRON_SZ) << log_xblk_sz)));
          /*Row 0 and 1 together.*/
          row01_lo = _mm_unpacklo_epi16(row0_src, row1_src);
          row01_hi = _mm_unpackhi_epi16(row0_src, row1_src);
          sums32_0to3 = _mm_madd_epi16(row01_lo, fy01);
          sums32_4to7 = _mm_madd_epi16(row01_hi, fy01);
          /*Row 2 and 3 together.*/
          row01_lo = _mm_unpacklo_epi16(row2_src, row3_src);
          row01_hi = _mm_unpackhi_epi16(row2_src, row3_src);
          sums32_0to3 = _mm_add_epi32(sums32_0to3,
           _mm_madd_epi16(row01_lo, fy23));
          sums32_4to7 = _mm_add_epi32(sums32_4to7,
           _mm_madd_epi16(row01_hi, fy23));
          /*Row 4 and 5 together.*/
          row01_lo = _mm_unpacklo_epi16(row4_src, row5_src);
          row01_hi = _mm_unpackhi_epi16(row4_src, row5_src);
          sums32_0to3 = _mm_add_epi32(sums32_0to3,
           _mm_madd_epi16(row01_lo, fy45));
          sums32_4to7 = _mm_add_epi32(sums32_4to7,
           _mm_madd_epi16(row01_hi, fy45));
          /*Add rounding offset, then scale down.*/
          sums32_0to3 = _mm_add_epi32(sums32_0to3, rounding_offset);
          sums32_0to3 = _mm_srai_epi32(sums32_0to3, OD_SUBPEL_COEFF_SCALE2);
          /*Write four filter output values to destination.*/
          /*if (i + 8 <= xblk_sz) {*/
          if (xblk_sz >= 8) {
            sums32_4to7 = _mm_add_epi32(sums32_4to7, rounding_offset);
            sums32_4to7 = _mm_srai_epi32(sums32_4to7, OD_SUBPEL_COEFF_SCALE2);
            out = _mm_packs_epi32(sums32_0to3, sums32_4to7);
            out = _mm_packus_epi16(out, out);
            _mm_storel_epi64((__m128i *)(dst_p + i), out);
          }
          else if (xblk_sz >= 4) {
            OD_ASSERT(i + 4 <= xblk_sz);
            out = _mm_packs_epi32(sums32_0to3, sums32_0to3);
            out = _mm_packus_epi16(out, out);
            *((uint32_t *)(dst_p + i)) = _mm_cvtsi128_si32(out);
          }
          else {
            OD_ASSERT(i + 2 <= xblk_sz);
            out = _mm_packs_epi32(sums32_0to3, sums32_0to3);
            out = _mm_packus_epi16(out, out);
            *((uint16_t *)(dst_p + i)) =
             (uint16_t)_mm_cvtsi128_si32(out);
          }
        }
        buff_p += xblk_sz;
        dst_p += xblk_sz;
      }
    }
    /*The mvy is of integer position.*/
    else {
      __m128i rounding_offset4;
      /*Note the sign!*/
      rounding_offset4 = _mm_set1_epi16(-OD_SUBPEL_RND_OFFSET4);
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i += 8) {
          __m128i p;
          p = _mm_loadu_si128((__m128i *)(buff_p + i));
          /*p is in the range [-26584, 26456].*/
          p = _mm_max_epi16(p, rounding_offset4);
          /*p is in the range [-16448, 26456].*/
          p = _mm_sub_epi16(p, rounding_offset4);
          /*p is in the range [0, 42904].*/
          p = _mm_srli_epi16(p, OD_SUBPEL_COEFF_SCALE);
          /*p is in the range [0, 335].*/
          p = _mm_packus_epi16(p, p);
          if (xblk_sz >= 8) {
            _mm_storel_epi64((__m128i *)(dst_p + i), p);
          }
          else if (xblk_sz >= 4) {
            OD_ASSERT(i + 4 <= xblk_sz);
            *((uint32_t *)(dst_p + i)) = _mm_cvtsi128_si32(p);
          }
          else {
            OD_ASSERT(i + 2 <= xblk_sz);
            *((uint16_t *)(dst_p + i)) =
             (uint16_t)_mm_cvtsi128_si32(p);
          }
        }
        buff_p += xblk_sz;
        dst_p += xblk_sz;
      }
    }
  }
  /*MC with full-pel MV, i.e. integer position.*/
  else {
    for (j = 0; j < yblk_sz; j++) {
      OD_COPY(dst_p, src_p, xblk_sz);
      src_p += systride;
      dst_p += xblk_sz;
    }
  }
#if defined(OD_CHECKASM)
  od_mc_predict1fmv8_check(dst, src, systride, mvx, mvy,
   log_xblk_sz, log_yblk_sz);
  /*fprintf(stderr,"od_mc_predict1fmv8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}

#endif

#if defined(OD_GCC_INLINE_ASSEMBLY)

#if defined(OD_CHECKASM)
static void od_mc_blend_full8_check(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  unsigned char  dst[OD_MVBSIZE_MAX*OD_MVBSIZE_MAX];
  int            xblk_sz;
  int            yblk_sz;
  int            failed;
  int            i;
  int            j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_blend_full8_c(dst,xblk_sz,_src,_log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if(dst[i+(j<<_log_xblk_sz)]!=(_dst+j*_dystride)[i]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         dst[i+(j<<_log_xblk_sz)],(_dst+j*_dystride)[i],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_predict1fmv8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif

/*Loads a block of 16 bytes from each of the 4 images into xmm0...xmm3.
  We swap images 2 and 3 here, so that the order more closely follows the
   natural rectilinear indexing, instead of the circular indexing the rest of
   the code uses.*/
#define OD_IM_LOAD16 \
  "#OD_IM_LOAD16\n\t" \
  "mov (%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm0\n\t" \
  /*The "c" prefix here means to leave off the leading $ normally used for \
     immediates (since we want to move from an address, not an immediate). \
    Currently, this is not documented in the gcc manual, and my editor doesn't \
     even syntax highlight it correctly; who knows what versions of gcc \
     support it. \
    It _is_ documented in the internals manual: \
     http://gcc.gnu.org/onlinedocs/gccint/Output-Template.html#Output-Template \
    I found it on this thread: \
     http://gcc.gnu.org/ml/gcc-help/2006-09/msg00301.html*/ \
  "mov %c[pstride](%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm1\n\t" \
  "mov %c[pstride]*3(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm2\n\t" \
  "mov %c[pstride]*2(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm3\n\t" \

/*Unpacks an 8-bit register into two 16-bit registers.
  _rega: The input register, which will contain the low-order output.
  _regb: Will contain the high-order output.
  _zero: A register containing the value zero.*/
#define OD_IM_UNPACK(_rega,_regb,_zero) \
  "#OD_IM_UNPACK\n\t" \
  "movdqa " _rega "," _regb "\n\t" \
  "punpcklbw " _zero "," _rega "\n\t" \
  "punpckhbw " _zero "," _regb "\n\t" \

/*Bilinearly blends 2 pairs of 16-bit registers.
  _reg0a:  The low-order register of the first pair, which will contain the
            output.
  _reg0b:  The high-order register of the first pair, which will contain the
            output.
  _reg1a:  The low-order register of the second pair.
  _reg1b:  The high-order register of the second pair.
  _shift:  The number of bits of precision the blending adds (e.g., the
            precision of the weights).
  _scalea: The weights to apply to the low-order register in the second pair.
           The weights applied to the low-order register in the first pair are
            (1<<_shift)-_scalea.
  _scaleb: The weights to apply to the high-order register in the second pair.
           The weights applied to the high-order register in the first pair are
            (1<<_shift)-_scaleb.*/
#define OD_IM_BLEND(_reg0a,_reg0b,_reg1a,_reg1b,_shift,_scalea,_scaleb) \
  "#OD_IM_BLEND\n\t" \
  "psubw " _reg0a "," _reg1a "\n\t" \
  "psubw " _reg0b "," _reg1b "\n\t" \
  "psllw " _shift "," _reg0a "\n\t" \
  "pmullw " _scalea "," _reg1a "\n\t" \
  "psllw " _shift "," _reg0b "\n\t" \
  "pmullw " _scaleb "," _reg1b "\n\t" \
  "paddw " _reg1a "," _reg0a "\n\t" \
  "paddw " _reg1b "," _reg0b "\n\t" \

/*Bilinearly blends 2 pairs of 16-bit registers and rounds, shifts, and packs
  the result in an 8-bit register. This is a special case version of
  OD_IM_BLEND that handles 16-bit overflow.
  It implements the following expression:
   ((a << x) + (b - a) * scale + (1 << (y - 1))) >> y
  Expressions (a << x) and ((b - a)*scale) can potentially overflow 16-bit
   integers, so we need to rewrite them as
    (a + (((b - a) * (scale << (16 - x))) >> 16) + (1 << (y - 1 - x)))
     >> (y - x)
   in order to use pmulhw
    (a + (pmulhw((b - a), (scale << (16 - x)))) + (1 << (y - 1 - x)))
     >> (y - x)
  The (scale << (16 - x)) expression can also overflow a signed 16-bit
   integer, and there is no signed-unsigned multiply instruction, so we
   further rewrite it as
    (a + (pmulhw((b - a) << w, (scale << (16 - x - w)))) + (1 << (y - 1 - x)))
     >> (y - x)

  _reg0a:  The low-order register of the first pair, which will contain the
            output.
  _reg0b:  The high-order register of the first pair.
  _reg1a:  The low-order register of the second pair.
  _reg1b:  The high-order register of the second pair.
  _scale:  The weights to apply to the low-order and high-order registers.
  _x:      Precision of the weights.
  _y:      Pack shift amount.
  _z:      Used to prevent overflow. Callers of this macro must choose a value
           for z such that ((b - a) << w) and (scale << (16 - x - w)) don't
           overflow 16-bit registers.
  */
#define OD_IM_BLEND_AND_PACK(_reg0a,_reg0b,_reg1a,_reg1b,_scale,_x,_y,_z) \
  "#OD_IM_BLEND_AND_PACK\n\t" \
  "psubw " _reg0a "," _reg1a "\n\t" \
  "psubw " _reg0b "," _reg1b "\n\t" \
  "psllw $" _z "," _reg1a "\n\t" \
  "psllw $" _z "," _reg1b "\n\t" \
  "psllw $16-" _x "-" _z "," _scale "\n\t" \
  "pmulhw " _scale "," _reg1a "\n\t" \
  "pmulhw " _scale "," _reg1b "\n\t" \
  "paddw " _reg1a "," _reg0a "\n\t" \
  "paddw " _reg1b "," _reg0b "\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "psllw $" _y "-1-" _x ",%%xmm7\n\t" \
  OD_IM_PACK(_reg0a,_reg0b,"%%xmm7","$" _y "-" _x) \

/*Rounds, shifts, and packs two 16-bit registers into one 8 bit register.
  _rega:  The low-order register, which will contain the output.
  _regb:  The high-order register.
  _round: The register containing the rounding offset.
  _shift: The immediate that is the amount to shift.*/
#define OD_IM_PACK(_rega,_regb,_round,_shift) \
  "#OD_IM_PACK\n\t" \
  "paddw " _round "," _rega "\n\t" \
  "paddw " _round "," _regb "\n\t" \
  "psrlw " _shift "," _rega "\n\t" \
  "psrlw " _shift "," _regb "\n\t" \
  "packuswb " _regb "," _rega "\n\t" \

/*Note: We lea the address of our blending weights into a register before
   invoking this macro because of PIC.
  We are CAPABLE of referencing the memory directly in all cases, but
    a) There's no way to tell gcc to combine the base register %ebx it uses for
     PIC with our constant offset, index register, and scale values (because
     gcc just concatenates strings; it does no assembly parsing itself),
    b) gcc is the only one who knows if we're using PIC or not,
    c) gcc is the only one who knows the real name of the (static) symbol.
  The resulting performance difference is unmeasurable.*/

/*Blends 4 rows of a 4xN block (N up to 64).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_4x4(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BIL4H],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa (%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$2","%%xmm6","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$2","%%xmm6","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "lea %[OD_BIL4V],%[a]\n\t" \
  "psllw $" #_log_yblk_sz "+1,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],4)","0x10(%[a],%[row],4)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+2") \
  /*Get it back out to memory. \
    We have to do this 4 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "lea (%[dst],%[dystride]),%[a]\n\t" \
  "psrldq $4,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst])\n\t" \
  "movd %%xmm1,(%[a])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "psrldq $8,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst],%[dystride],2)\n\t" \
  "movd %%xmm1,(%[a],%[dystride],2)\n\t" \

/*Blends 2 rows of an 8xN block (N up to 32).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_8x2(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa (%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$3","%%xmm6","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$3","%%xmm6","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  "psllw $" #_log_yblk_sz "+2,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],8)","0x10(%[a],%[row],8)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+3") \
  /*Get it back out to memory. \
    We have to do this 8 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movq %%xmm0,(%[dst])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "movq %%xmm0,(%[dst],%[dystride])\n\t" \

/*Blends 1 row of a 16xN block (N up to 16).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_16x1(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa 0x10(%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$4","(%[a])","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$4","(%[a])","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a],%[row],4),%%xmm6\n\t" \
  "psllw $" #_log_yblk_sz "+3,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "%%xmm6","%%xmm6") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+4") \
  /*Get it back out to memory.*/ \
  "movdqa %%xmm0,(%[dst])\n\t" \

/*Blends half a row of a 32x32 block. We can't blend an entire row at a time
  in SSE2 so we need to split it two havles.
 _bilh_offset: Offset in the BILH table where to read the weights from.
 _bilv_offset: Offset in the BILV table where to read the weights from.
 _dst_offset: Offset in the dst pointer where to write the blended row. */
#define OD_MC_BLEND_FULL8_32_HALF(_bilh_offset, _bilv_offset, _dst_offset) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa " _bilh_offset "+0x10(%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$5", \
   _bilh_offset "(%[a])","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$5", \
   _bilh_offset "(%[a])","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "lea %[OD_BILV],%[a]\n\t" \
  "movdqa " _bilv_offset "(%[a],%[row],2),%%xmm6\n\t" \
  OD_IM_BLEND_AND_PACK("%%xmm0","%%xmm4","%%xmm2","%%xmm5","%%xmm6", \
   "5","10","1") \
  /*Get it back out to memory.*/ \
  "movdqa %%xmm0," _dst_offset "(%[dst])\n\t" \

#if 0
/*Defines a pure-C implementation with hard-coded loop limits for block sizes
   we don't want to implement manually (e.g., that have fewer than 16 bytes,
   require byte-by-byte unaligned loads, etc.).
  This should let the compiler aggressively unroll loops, etc.
  It can't vectorize it itself because of the difference in operand sizes.*/
#define OD_MC_BLEND_FULL8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full8_##_n##x##_m(unsigned char *_dst,int _dystride, \
 const unsigned char *_src[4]){ \
  int      o; \
  unsigned a; \
  unsigned b; \
  int      i; \
  int      j; \
  o=0; \
  for(j=0;j<(_m);j++){ \
    for(i=0;i<(_n);i++){ \
      a=(_src[0][o+i]<<(_log_xblk_sz))+(_src[1][o+i]-_src[0][o+i])*i; \
      b=(_src[3][o+i]<<(_log_xblk_sz))+(_src[2][o+i]-_src[3][o+i])*i; \
      _dst[i]=(unsigned char)((a<<(_log_yblk_sz))+(b-a)*j+ \
       (1<<(_log_xblk_sz)+(_log_yblk_sz))/2>>(_log_xblk_sz)+(_log_yblk_sz)); \
    } \
    o+=(_m); \
    _dst+=_dystride; \
  } \
} \

#else
/*With -O3 and inter-module optimization, gcc inlines these anyway.
  I'd rather leave the choice to the compiler.*/
#define OD_MC_BLEND_FULL8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full8_##_n##x##_m(unsigned char *_dst,int _dystride, \
 const unsigned char *_src[4]){ \
  od_mc_blend_full8_c(_dst,_dystride,_src,_log_xblk_sz,_log_yblk_sz); \
} \

#endif

OD_MC_BLEND_FULL8_C(1,1,0,0)
OD_MC_BLEND_FULL8_C(1,2,0,1)
OD_MC_BLEND_FULL8_C(1,4,0,2)
OD_MC_BLEND_FULL8_C(1,8,0,3)
OD_MC_BLEND_FULL8_C(1,16,0,4)

OD_MC_BLEND_FULL8_C(2,1,1,0)
OD_MC_BLEND_FULL8_C(2,2,1,1)
OD_MC_BLEND_FULL8_C(2,4,1,2)
OD_MC_BLEND_FULL8_C(2,8,1,3)
OD_MC_BLEND_FULL8_C(2,16,1,4)

OD_MC_BLEND_FULL8_C(4,1,2,0)
OD_MC_BLEND_FULL8_C(4,2,2,1)

static void od_mc_blend_full8_4x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_4x4(2)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    /*Note that we pass the constant 0 for [row] here.
      We'll still use it in indexing expression in the asm, but the overhead is
       negligible, and it's easier than writing a special case of
       OD_MC_BLEND_FULL8_4x4 for it.*/
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_4x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  /*We use loops like these so gcc can decide to unroll them if it wants (and
     so that it can do the jumps, etc., for us, instead of trying to figure out
     how to put that in a macro).
    row can't count by 1, because (in the 8x2 versions below) we will need to
     scale it by either 16 or 32, and an index register can only be scaled by
     16, and in the 32x32 and 64x64 versions we will want to scale it by 1/2
     and 1/4 (respectively) to reduce the size of our tables.
    Therefore we pre-scale it by 4, and do so everywhere (even for the 4x4 and
     8x2 versions) so that things are consistent.*/
  for(row=0;row<8;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_4x4(3)
      "lea (%[dst],%[dystride],4),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_4x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_4x4(4)
      "lea (%[dst],%[dystride],4),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

OD_MC_BLEND_FULL8_C(8,1,3,0)

static void od_mc_blend_full8_8x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_8x2(1)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_8x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(2)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_8x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8*8/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_8x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8*16/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(4)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x1(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_16x1(0)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_16x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16*2/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(1)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(2)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16*8/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(3)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16*16/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(4)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_32x32(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<32*32/4;row+=4){
    __asm__ __volatile__(
      /*First 16 bytes.*/ \
      OD_MC_BLEND_FULL8_32_HALF("0x00","0","0x00") \
      /*Second 16 bytes.*/ \
      "lea 4(%[row]),%[row]\t\n" \
      OD_MC_BLEND_FULL8_32_HALF("0x20","-8","0x10") \
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

typedef void (*od_mc_blend_full8_fixed_func)(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]);

/*Perform normal bilinear blending.*/
void od_mc_blend_full8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  static const od_mc_blend_full8_fixed_func
   VTBL[OD_LOG_MVBSIZE_MAX + 1][OD_LOG_MVBSIZE_MAX + 1]={
    {
      od_mc_blend_full8_1x1,od_mc_blend_full8_1x2,
      od_mc_blend_full8_1x4,od_mc_blend_full8_1x8,
      od_mc_blend_full8_1x16
    },
    {
      od_mc_blend_full8_2x1,od_mc_blend_full8_2x2,
      od_mc_blend_full8_2x4,od_mc_blend_full8_2x8,
      od_mc_blend_full8_2x16
    },
    {
      od_mc_blend_full8_4x1,od_mc_blend_full8_4x2,
      od_mc_blend_full8_4x4,od_mc_blend_full8_4x8,
      od_mc_blend_full8_4x16
    },
    {
      od_mc_blend_full8_8x1,od_mc_blend_full8_8x2,
      od_mc_blend_full8_8x4,od_mc_blend_full8_8x8,
      od_mc_blend_full8_8x16
    },
    {
      od_mc_blend_full8_16x1,od_mc_blend_full8_16x2,
      od_mc_blend_full8_16x4,od_mc_blend_full8_16x8,
      od_mc_blend_full8_16x16
    },
    /*These NULLs are placeholders because we do not currently call this
       function for these sizes.
      If you need one of them, add another OD_MC_BLEND_FULL8_C() line above.*/
    {
      NULL, NULL, NULL, NULL, NULL, od_mc_blend_full8_32x32
    }
  };
  (*VTBL[_log_xblk_sz][_log_yblk_sz])(_dst,_dystride,_src);
#if defined(OD_CHECKASM)
  od_mc_blend_full8_check(_dst,_dystride,_src,_log_xblk_sz,_log_yblk_sz);
  /*fprintf(stderr,"od_mc_blend_full8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}

#if defined(OD_CHECKASM)
void od_mc_blend_full_split8_check(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz){
  unsigned char dst[OD_MVBSIZE_MAX][OD_MVBSIZE_MAX];
  int           xblk_sz;
  int           yblk_sz;
  int           failed;
  int           i;
  int           j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_blend_full_split8_c(dst[0],sizeof(dst[0]),_src,_c,_s,
   _log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if((_dst+j*_dystride)[i]!=dst[j][i]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         (_dst+j*_dystride)[i],dst[j][i],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_blend_full_split8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif

/*Loads a block of 16 bytes from each the first 2 images into xmm0...xmm3.
  xmm2 and xmm3 contain duplicate copies of xmm0 and xmm1, or not, depending on
   whether the block edges are split or not.*/
#define OD_IM_LOAD16A \
  "#OD_IM_LOAD16A\n\t" \
  "mov (%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm0\n\t" \
  "mov %c[pstride](%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm1\n\t" \
  "mov %c[pstride]*4(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm2\n\t" \
  "mov %c[pstride]*5(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm3\n\t" \

/*Loads a block of 16 bytes from the third image into xmm2 and xmm1.
  xmm1 contains a duplicate copy of xmm2, or not, depending on whether the
   block edge is split or not.*/
#define OD_IM_LOAD16B \
  "#OD_IM_LOAD16B\n\t" \
  "mov %c[pstride]*3(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm2\n\t" \
  "mov %c[pstride]*7(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm1\n\t" \

/*Loads a block of 16 bytes from the fourth image into xmm3 and xmm1.
  xmm1 contains a duplicate copy of xmm3, or not, depending on whether the
   block edge is split or not.*/
#define OD_IM_LOAD16C \
  "#OD_IM_LOAD16C\n\t" \
  "mov %c[pstride]*2(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm3\n\t" \
  "mov %c[pstride]*6(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm1\n\t" \

/*Blends 4 rows of a 4xN block with split edges (N up to 32).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL_SPLIT8_4x4(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the first two images to blend.*/ \
  OD_IM_LOAD16A \
  /*Unpack and merge the 0 image.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  OD_IM_UNPACK("%%xmm2","%%xmm6","%%xmm7") \
  "paddw %%xmm2,%%xmm0\n\t" \
  "paddw %%xmm6,%%xmm4\n\t" \
  /*Unpack and merge the 1 image.*/ \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  "lea %[OD_BIL4H],%[a]\n\t" \
  "paddw %%xmm3,%%xmm1\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  "movdqa (%[a]),%%xmm6\n\t" \
  /*Blend the 0 and 1 images.*/ \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$2","%%xmm6","%%xmm6") \
  /*Load, unpack, and merge the 2 image.*/ \
  OD_IM_LOAD16B \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm1","%%xmm6","%%xmm7") \
  "paddw %%xmm1,%%xmm2\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  /*Load, unpack, and merge the 3 image.*/ \
  OD_IM_LOAD16C \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  /*We run out of registers here, and have to overwrite our zero. \
    Fortunately, the first row of OD_BILV is also all zeros. \
  "lea %[OD_BILV],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm7","(%[a])") \
  "paddw %%xmm1,%%xmm3\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "paddw %%xmm7,%%xmm6\n\t" \
  "pxor %%xmm7,%%xmm7\n\t"*/ \
  /*Alternate version: Saves 1 memory reference, but has a longer dependency \
     chain.*/ \
  "movdqa %%xmm1,%%xmm7\n\t" \
  "punpcklbw %[OD_BILV],%%xmm7\n\t" \
  "paddw %%xmm7,%%xmm3\n\t" \
  "pxor %%xmm7,%%xmm7\n\t" \
  "punpckhbw %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm6\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  /*End alternate version.*/ \
  "lea %[OD_BIL4H],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a]),%%xmm1\n\t" \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm6","$2","%%xmm1","%%xmm1") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "psllw $" #_log_yblk_sz "+2,%%xmm7\n\t" \
  "lea %[OD_BIL4V],%[a]\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],4)","0x10(%[a],%[row],4)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+3") \
  /*Get it back out to memory. \
    We have to do this 4 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "lea (%[dst],%[dystride]),%[a]\n\t" \
  "psrldq $4,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst])\n\t" \
  "movd %%xmm1,(%[a])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "psrldq $8,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst],%[dystride],2)\n\t" \
  "movd %%xmm1,(%[a],%[dystride],2)\n\t" \

/*Blends 2 rows of an 8xN block with split edges (N up to 16).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL_SPLIT8_8x2(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the first two images to blend.*/ \
  OD_IM_LOAD16A \
  /*Unpack and merge the 0 image.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  OD_IM_UNPACK("%%xmm2","%%xmm6","%%xmm7") \
  "paddw %%xmm2,%%xmm0\n\t" \
  "paddw %%xmm6,%%xmm4\n\t" \
  /*Unpack and merge the 1 image.*/ \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  "paddw %%xmm3,%%xmm1\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  "movdqa (%[a]),%%xmm6\n\t" \
  /*Blend the 0 and 1 images.*/ \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$3","%%xmm6","%%xmm6") \
  /*Load, unpack, and merge the 2 image.*/ \
  OD_IM_LOAD16B \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm1","%%xmm6","%%xmm7") \
  "paddw %%xmm1,%%xmm2\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  /*Load, unpack, and merge the 3 image.*/ \
  OD_IM_LOAD16C \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  /*"lea %[OD_BILV],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm7","(%[a])") \
  "paddw %%xmm1,%%xmm3\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "paddw %%xmm7,%%xmm6\n\t" \
  "pxor %%xmm7,%%xmm7\n\t"*/ \
  /*Alternate version: Saves 1 memory reference, but has a longer dependency \
     chain.*/ \
  "movdqa %%xmm1,%%xmm7\n\t" \
  "punpcklbw %[OD_BILV],%%xmm7\n\t" \
  "paddw %%xmm7,%%xmm3\n\t" \
  "pxor %%xmm7,%%xmm7\n\t" \
  "punpckhbw %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm6\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  /*End alternate version.*/ \
  "lea %[OD_BILH],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a]),%%xmm1\n\t" \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm6","$3","%%xmm1","%%xmm1") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "psllw $" #_log_yblk_sz "+3,%%xmm7\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],8)","0x10(%[a],%[row],8)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+4") \
  /*Get it back out to memory. \
    We have to do this 8 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movq %%xmm0,(%[dst])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "movq %%xmm0,(%[dst],%[dystride])\n\t" \

/*Blends 1 row of an 16xN block with split edges (N up to 16).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL_SPLIT8_16x1(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the first two images to blend.*/ \
  OD_IM_LOAD16A \
  /*Unpack and merge the 0 image.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  OD_IM_UNPACK("%%xmm2","%%xmm6","%%xmm7") \
  "paddw %%xmm2,%%xmm0\n\t" \
  "paddw %%xmm6,%%xmm4\n\t" \
  /*Unpack and merge the 1 image.*/ \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  "paddw %%xmm3,%%xmm1\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  "movdqa 0x10(%[a]),%%xmm6\n\t" \
  /*Blend the 0 and 1 images.*/ \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$4","(%[a])","%%xmm6") \
  /*Load, unpack, and merge the 2 image.*/ \
  OD_IM_LOAD16B \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm1","%%xmm6","%%xmm7") \
  "paddw %%xmm1,%%xmm2\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  /*Load, unpack, and merge the 3 image.*/ \
  OD_IM_LOAD16C \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  /*"lea %[OD_BILV],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm7","(%[a])") \
  "paddw %%xmm1,%%xmm3\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "paddw %%xmm7,%%xmm6\n\t" \
  "pxor %%xmm7,%%xmm7\n\t"*/ \
  /*Alternate version: Saves 1 memory reference, but has a longer dependency \
     chain.*/ \
  "movdqa %%xmm1,%%xmm7\n\t" \
  "punpcklbw %[OD_BILV],%%xmm7\n\t" \
  "paddw %%xmm7,%%xmm3\n\t" \
  "pxor %%xmm7,%%xmm7\n\t" \
  "punpckhbw %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm6\n\t" \
  /*End alternate version.*/ \
  "lea %[OD_BILH],%[a]\n\t" \
  "movdqa 0x10(%[a]),%%xmm1\n\t" \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm6","$4","(%[a])","%%xmm1") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "lea %[OD_BILV],%[a]\n\t" \
  "movdqa (%[a],%[row],4),%%xmm6\n\t" \
  OD_IM_BLEND_AND_PACK("%%xmm0","%%xmm4","%%xmm2","%%xmm5", "%%xmm6", \
   "4","9","1") \
  /*Get it back out to memory.*/ \
  "movdqa %%xmm0,(%[dst])\n\t" \

#if 0
/*Defines a pure-C implementation with hard-coded loop limits for block sizes
   we don't want to implement manually (e.g., that have fewer than 16 bytes,
   require byte-by-byte unaligned loads, etc.).
  This should let the compiler aggressively unroll loops, etc.
  It can't vectorize it itself because of the difference in operand sizes.*/
/*TODO: This approach (using pointer aliasing to allow us to use normal
   bilinear weights) might actually be faster than the pure-C routine we're
   currently using, which adjusts the weights.
  This should be investigated.*/
#define OD_MC_BLEND_FULL_SPLIT8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full_split8_##_n##x##_m(unsigned char *_dst, \
 int _dystride,const unsigned char *_src[8]){ \
  int      o; \
  unsigned a; \
  unsigned b; \
  int      i; \
  int      j; \
  o=0; \
  for(j=0;j<(_m);j++){ \
    for(i=0;i<(_n);i++){ \
      a=(_src[0][o+i]+_src[4+0][o+i]<<(_log_xblk_sz))+ \
       (_src[1][o+i]-_src[0][o+i]+_src[4+1][o+i]-_src[4+0][o+i])*i; \
      b=(_src[3][o+i]+_src[4+3][o+i]<<(_log_xblk_sz))+ \
       (_src[2][o+i]-_src[3][o+i]+_src[4+2][o+i]-_src[4+3][o+i])*i; \
      _dst[i]=(unsigned char)((a<<(_log_yblk_sz))+(b-a)*j+ \
       (1<<(_log_xblk_sz)+(_log_yblk_sz))>>(_log_xblk_sz)+(_log_yblk_sz)+1); \
    } \
    o+=(_m); \
    _dst+=_dystride; \
  } \
} \

#else
/*TODO: This approach (using pointer aliasing to allow us to use normal
   bilinear weights) might actually be faster than the pure-C routine we're
   currently using, which adjusts the weights.
  This should be investigated.*/
static void od_mc_blend_full_split8_bil_c(unsigned char *_dst,
 int _dystride,const unsigned char *_src[8],int _log_xblk_sz,int _log_yblk_sz){
  int      xblk_sz;
  int      yblk_sz;
  int      round;
  int      o;
  unsigned a;
  unsigned b;
  int      i;
  int      j;
  o=0;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  round=1<<(_log_xblk_sz+_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      a=((_src[0][o+i]+_src[4+0][o+i])<<_log_xblk_sz)+
       (_src[1][o+i]-_src[0][o+i]+_src[4+1][o+i]-_src[4+0][o+i])*i;
      b=((_src[3][o+i]+_src[4+3][o+i])<<_log_xblk_sz)+
       (_src[2][o+i]-_src[3][o+i]+_src[4+2][o+i]-_src[4+3][o+i])*i;
      _dst[i]=(unsigned char)(((a<<_log_yblk_sz)+(b-a)*j+
       round)>>(_log_xblk_sz+_log_yblk_sz+1));
    }
    o+=xblk_sz;
    _dst+=_dystride;
  }
}

/*With -O3 and inter-module optimization, gcc inlines these anyway.
  I'd rather leave the choice to the compiler.*/
#define OD_MC_BLEND_FULL_SPLIT8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full_split8_##_n##x##_m(unsigned char *_dst, \
 int _dystride,const unsigned char *_src[8]){ \
  od_mc_blend_full_split8_bil_c(_dst,_dystride,_src, \
  _log_xblk_sz,_log_yblk_sz); \
} \

#endif

OD_MC_BLEND_FULL_SPLIT8_C(1,1,0,0)
OD_MC_BLEND_FULL_SPLIT8_C(1,2,0,1)
OD_MC_BLEND_FULL_SPLIT8_C(1,4,0,2)
OD_MC_BLEND_FULL_SPLIT8_C(1,8,0,3)

OD_MC_BLEND_FULL_SPLIT8_C(2,1,1,0)
OD_MC_BLEND_FULL_SPLIT8_C(2,2,1,1)
OD_MC_BLEND_FULL_SPLIT8_C(2,4,1,2)
OD_MC_BLEND_FULL_SPLIT8_C(2,8,1,3)

OD_MC_BLEND_FULL_SPLIT8_C(4,1,2,0)
OD_MC_BLEND_FULL_SPLIT8_C(4,2,2,1)

static void od_mc_blend_full_split8_4x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL_SPLIT8_4x4(2)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full_split8_4x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_4x4(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

OD_MC_BLEND_FULL_SPLIT8_C(8,1,3,0)

static void od_mc_blend_full_split8_8x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL_SPLIT8_8x2(1)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full_split8_8x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_8x2(2)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full_split8_8x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<8*8/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_8x2(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full_split8_16x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<16*16/4;row+=4){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_16x1(4)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

typedef void (*od_mc_blend_full_split8_fixed_func)(unsigned char *_dst,
 int _dystride,const unsigned char *_src[8]);

/*Sets up a second set of image pointers based on the given split state to
   properly shift weight from one image to another.*/
static void od_mc_setup_split_ptrs(const unsigned char *_drc[4],
 const unsigned char *_src[4],int _c,int _s){
  int j;
  int k;
  _drc[_c]=_src[_c];
  j=(_c+(_s&1))&3;
  k=(_c+1)&3;
  _drc[k]=_src[j];
  j=(_c+(_s&2)+((_s&2)>>1))&3;
  k=(_c+3)&3;
  _drc[k]=_src[j];
  k=_c^2;
  _drc[k]=_src[k];
}

/*Perform normal bilinear blending.*/
void od_mc_blend_full_split8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz){
  static const od_mc_blend_full_split8_fixed_func
   VTBL[OD_LOG_MVBSIZE_MAX][OD_LOG_MVBSIZE_MAX]={
    {
      od_mc_blend_full_split8_1x1,od_mc_blend_full_split8_1x2,
      od_mc_blend_full_split8_1x4,od_mc_blend_full_split8_1x8,
    },
    {
      od_mc_blend_full_split8_2x1,od_mc_blend_full_split8_2x2,
      od_mc_blend_full_split8_2x4,od_mc_blend_full_split8_2x8,
    },
    {
      od_mc_blend_full_split8_4x1,od_mc_blend_full_split8_4x2,
      od_mc_blend_full_split8_4x4,od_mc_blend_full_split8_4x8,
    },
    {
      od_mc_blend_full_split8_8x1,od_mc_blend_full_split8_8x2,
      od_mc_blend_full_split8_8x4,od_mc_blend_full_split8_8x8,
    },
    /*These NULLs are placeholders because we do not currently call this
       function for these sizes.
      If you need one of them, add another OD_MC_BLEND_FULL_SPLIT8_C() line
       above.*/
    {
      NULL, NULL, NULL, NULL, od_mc_blend_full_split8_16x16
    }
  };
  /*We pack all the image pointers in one array to save a register.*/
  const unsigned char *drc[8];
  memcpy(drc,_src,sizeof(*drc)*4);
  od_mc_setup_split_ptrs(drc+4,drc,_c,_s);
  (*VTBL[_log_xblk_sz][_log_yblk_sz])(_dst,_dystride,drc);
#if defined(OD_CHECKASM)
  od_mc_blend_full_split8_check(_dst,_dystride,_src,_c,_s,
   _log_xblk_sz,_log_yblk_sz);
  /*fprintf(stderr,"od_mc_blend_full_split8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}
#endif

#endif
