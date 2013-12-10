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

/*This should translate directly to 3 or 4 instructions for a constant _b:
#define OD_UNBIASED_RSHIFT(_a,_b) ((_a)+(((1<<(_b))-1)&-((_a)<0))>>(_b))*/
/*This version relies on a smart compiler:*/
# define OD_UNBIASED_RSHIFT(_a, _b) ((_a)/(1<<(_b)))

# if 0
#  define OD_DCT_RSHIFT(_a, _b) OD_DIV_POW2_RE(_a, _b)
# elif 1
#  define OD_DCT_RSHIFT(_a, _b) OD_UNBIASED_RSHIFT(_a, _b)
# else
#  define OD_DCT_RSHIFT(_a, _b) ((_a)>>(_b))
# endif

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

void od_pre_filter4(od_coeff _y[4], const od_coeff _x[4]);
void od_post_filter4(od_coeff _x[4], const od_coeff _y[4]);
void od_pre_filter8(od_coeff _y[8], const od_coeff _x[8]);
void od_post_filter8(od_coeff _x[8], const od_coeff _y[8]);
void od_pre_filter16(od_coeff _y[16], const od_coeff _x[16]);
void od_post_filter16(od_coeff _x[16], const od_coeff _y[16]);

# define OD_TOP_EDGE    (1<<3)
# define OD_RIGHT_EDGE  (1<<2)
# define OD_BOTTOM_EDGE (1<<1)
# define OD_LEFT_EDGE   (1<<0)

void od_apply_prefilter(od_coeff *c, int w, int bx, int by, unsigned char l,
 const unsigned char *bsize, int bstride, int xdec, int ydec, int edge);
void od_apply_postfilter(od_coeff *c, int w, int bx, int by, unsigned char l,
 const unsigned char *bsize, int bstride, int xdec, int ydec, int edge);

#endif
