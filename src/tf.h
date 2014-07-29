/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_tf_H)
# define _tf_H (1)
# include "filter.h"

/*This is an in-place, reversible, orthonormal Haar transform in 7 adds,
   1 shift (2 operations per sample).
  It is its own inverse (but requires swapping lh and hl on one side for
   bit-exact reversibility).
  It is defined in a macro here so it can be reused in various places.*/
#define OD_HAAR_KERNEL(ll, lh, hl, hh) \
  do { \
    od_coeff llmhh_2__; \
    (ll) += (hl); \
    (hh) -= (lh); \
    llmhh_2__ = ((ll) - (hh)) >> 1; \
    (lh) = llmhh_2__ - (lh); \
    (hl) = llmhh_2__ - (hl); \
    (ll) -= (lh); \
    (hh) += (hl); \
  } \
  while(0)

void od_tf_up_h_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dx, int n);

void od_tf_up_v_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dy, int n);

void od_tf_up_hv_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dx, int dy, int n);

void od_convert_block_down(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int curr_size, int dest_size, int filter);

void od_tf_filter_2d(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int n);
void od_tf_filter_inv_2d(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int n);

#endif
