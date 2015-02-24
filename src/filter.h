/*Daala video codec
Copyright (c) 2003-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_filter_H)
# define _filter_H (1)
# include "internal.h"

typedef ogg_int32_t od_coeff;
# define OD_COEFF_BITS (32)

/*There are 4 filter sizes total (4-point, 8-point, 16-point and 32-point).*/
# define OD_NFILTER_SIZES (4)

/*This is the strength reduced version of ((_a)/(1 << (_b))).
  This will not work for _b == 0, however currently this is only used for
   b == 1 anyway.*/
# define OD_UNBIASED_RSHIFT32(_a, _b) \
 (((ogg_int32_t)(((ogg_uint32_t)(_a) >> (32 - (_b))) + (_a))) >> (_b))

# define OD_DCT_RSHIFT(_a, _b) OD_UNBIASED_RSHIFT32(_a, _b)

typedef void (*od_filter_func)(od_coeff _out[], const od_coeff _in[]);

extern const od_filter_func OD_PRE_FILTER[OD_NBSIZES];
extern const od_filter_func OD_POST_FILTER[OD_NBSIZES];

/*These are the pre/post filtering functions used by Daala.
  The idea is to pre/post filter in the spatial domain (the time domain in
   signal processing terms) to improve the energy compaction as well as reduce
   or eliminate blocking artifacts.*/

extern const int OD_FILTER_PARAMS4[4];
extern const int OD_FILTER_PARAMS8[10];
extern const int OD_FILTER_PARAMS16[22];
extern const int OD_FILTER_PARAMS32[46];

void od_pre_filter4(od_coeff _y[4], const od_coeff _x[4]);
void od_post_filter4(od_coeff _x[4], const od_coeff _y[4]);
void od_pre_filter8(od_coeff _y[8], const od_coeff _x[8]);
void od_post_filter8(od_coeff _x[8], const od_coeff _y[8]);
void od_pre_filter16(od_coeff _y[16], const od_coeff _x[16]);
void od_post_filter16(od_coeff _x[16], const od_coeff _y[16]);
void od_pre_filter32(od_coeff _y[32], const od_coeff _x[32]);
void od_post_filter32(od_coeff _x[32], const od_coeff _y[32]);

void od_apply_prefilter_frame(od_coeff *c, int w, int nhsb, int nvsb,
 const unsigned char *bsize, int bstride, int dec);
void od_apply_postfilter_frame(od_coeff *c, int w, int nhsb, int nvsb,
 const unsigned char *bsize, int bstride, int dec);

extern const int OD_FILT_SIZE[OD_NBSIZES];
void od_apply_filter_sb_rows(od_coeff *c, int stride, int nhsb, int nvsb,
 int xdec, int ydec, int inv, int ln);
void od_apply_filter_sb_cols(od_coeff *c, int stride, int nhsb, int nvsb,
 int xdec, int ydec, int inv, int ln);
void od_apply_filter_hsplit(od_coeff *c0, int stride, int inv, int ln, int f);
void od_apply_filter_vsplit(od_coeff *c0, int stride, int inv, int ln, int f);

# if defined(OD_DCT_TEST) && defined(OD_DCT_CHECK_OVERFLOW)
#  include <stdio.h>

extern int od_dct_check_min[];
extern int od_dct_check_max[];

#  define OD_DCT_OVERFLOW_CHECK(val, scale, offset, idx) \
  do { \
    od_dct_check_min[(idx)] = OD_MINI(od_dct_check_min[(idx)], val); \
    od_dct_check_max[(idx)] = OD_MAXI(od_dct_check_max[(idx)], val); \
    OD_ASSERT((offset) >= 0); \
    if ((scale) > 0) { \
      if ((val) < INT_MIN/(scale)) { \
        printf("Overflow %2i: 0x%08X*0x%08X < INT_MIN\n", \
         (idx), (val), (scale)); \
      } \
      if ((val) > (INT_MAX - (offset))/(scale)) { \
        printf("Overflow %2i: 0x%08X*0x%04X + 0x%04X > INT_MAX\n", \
         (idx), (val), (scale), (offset)); \
      } \
    } \
    else if ((scale) < 0) { \
      if ((val) > (INT_MIN)/(scale)) { \
        printf("Overflow %2i: 0x%08X*-0x%04X < INT_MIN\n", \
         (idx), (val), (scale)); \
      } \
      if ((val) < (INT_MAX - (offset))/(scale)) { \
        printf("Overflow %2i: 0x%08X*-0x%08X + 0x%04X > INT_MAX\n", \
         (idx), (val), (scale), (offset)); \
      } \
    } \
  } \
  while(0)
# else
#  define OD_DCT_OVERFLOW_CHECK(val, scale, offset, idx)
# endif

#endif
