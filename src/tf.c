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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "tf.h"
#include "block_size.h"
#if defined(OD_VALGRIND)
# include <valgrind/memcheck.h>
#endif

/*Increase horizontal frequency resolution of an entire block and return the LF
   half.*/
void od_tf_up_h_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dx, int n) {
  int x;
  int y;
  for (y = 0; y < n; y++) {
    for (x = 0; x < n >> 1; x++) {
      od_coeff ll;
      od_coeff lh;
      int hswap;
      ll = src[y*sstride + x];
      lh = src[y*sstride + x + dx];
      lh = ll - lh;
      ll -= OD_DCT_RSHIFT(lh, 1);
      hswap = x & 1;
      dst[y*dstride + 2*x + hswap] = ll;
      dst[y*dstride + 2*x + 1 - hswap] = lh;
    }
  }
}

/*Increase vertical frequency resolution of an entire block and return the LF
   half.*/
void od_tf_up_v_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dy, int n) {
  int x;
  int y;
  for (y = 0; y < n >> 1; y++) {
    int vswap;
    vswap = y & 1;
    for (x = 0; x < n; x++) {
      od_coeff ll;
      od_coeff hl;
      ll = src[y*sstride + x];
      hl = src[(y + dy)*sstride + x];
      hl = ll - hl;
      ll -= OD_DCT_RSHIFT(hl, 1);
      dst[(2*y + vswap)*dstride + x] = ll;
      dst[(2*y + 1 - vswap)*dstride + x] = hl;
    }
  }
}

/*Increase horizontal and vertical frequency resolution of an entire block and
   return the LF quarter.*/
void od_tf_up_hv_lp(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int dx, int dy, int n) {
  int x;
  int y;
  for (y = 0; y < n >> 1; y++) {
    int vswap;
    vswap = y & 1;
    for (x = 0; x < n >> 1; x++) {
      od_coeff ll;
      od_coeff lh;
      od_coeff hl;
      od_coeff hh;
      od_coeff lhmhl_2;
      int hswap;
      ll = src[y*sstride + x];
      lh = src[y*sstride + x + dx];
      hl = src[(y + dy)*sstride + x];
      hh = src[(y + dy)*sstride + x + dx];
      hl = ll - hl;
      lh += hh;
      lhmhl_2 = OD_DCT_RSHIFT(lh - hl, 1);
      ll += lhmhl_2;
      hh -= lhmhl_2;
      lh = ll - lh;
      hl -= hh;
      hswap = x & 1;
      dst[(2*y + vswap)*dstride + 2*x + hswap] = ll;
      dst[(2*y + vswap)*dstride + 2*x + 1 - hswap] = lh;
      dst[(2*y + 1 - vswap)*dstride + 2*x + hswap] = hl;
      dst[(2*y + 1 - vswap)*dstride + 2*x + 1 - hswap] = hh;
    }
  }
}

/*Increase horizontal and vertical frequency resolution of a 2x2 group of
   blocks, combining them into a single block.*/
void od_tf_up_hv(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int n) {
  int x;
  int y;
  for (y = 0; y < n; y++) {
    int vswap;
    vswap = y & 1;
    for (x = 0; x < n; x++) {
      od_coeff ll;
      od_coeff lh;
      od_coeff hl;
      od_coeff hh;
      od_coeff lhmhl_2;
      int hswap;
      ll = src[y*sstride + x];
      lh = src[y*sstride + x + n];
      hl = src[(y + n)*sstride + x];
      hh = src[(y + n)*sstride + x + n];
      /*This kernel is identical to that of od_tf_down_hv with the roles of
         hl and lh swapped.*/
      hl = ll - hl;
      lh += hh;
      lhmhl_2 = OD_DCT_RSHIFT(lh - hl, 1);
      ll += lhmhl_2;
      hh -= lhmhl_2;
      lh = ll - lh;
      hl -= hh;
      hswap = x & 1;
      dst[(2*y + vswap)*dstride + 2*x + hswap] = ll;
      dst[(2*y + vswap)*dstride + 2*x + 1 - hswap] = lh;
      dst[(2*y + 1 - vswap)*dstride + 2*x + hswap] = hl;
      dst[(2*y + 1 - vswap)*dstride + 2*x + 1 - hswap] = hh;
    }
  }
}

/*Increase horizontal and vertical time resolution of a block, splitting it
   into a 2x2 group of blocks.*/
void od_tf_down_hv(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, int n) {
  int x;
  int y;
  OD_ASSERT(!(n & 1));
  n >>= 1;
  for (y = 0; y < n; y++) {
    int vswap;
    vswap = y & 1;
    for (x = 0; x < n; x++) {
      od_coeff ll;
      od_coeff lh;
      od_coeff hl;
      od_coeff hh;
      od_coeff lhmhl_2;
      int hswap;
      hswap = x & 1;
      ll = src[(2*y + vswap)*sstride + 2*x + hswap];
      lh = src[(2*y + vswap)*sstride + 2*x + 1 - hswap];
      hl = src[(2*y + 1 - vswap)*sstride+2*x + hswap];
      hh = src[(2*y + 1 - vswap)*sstride+2*x + 1 - hswap];
      /*This kernel is identical to that of od_tf_up_hv with the roles of
         hl and lh swapped.*/
      lh = ll - lh;
      hl += hh;
      lhmhl_2 = OD_DCT_RSHIFT(lh - hl, 1);
      ll -= lhmhl_2;
      hh += lhmhl_2;
      hl = ll - hl;
      lh -= hh;
      dst[y*dstride + x] = ll;
      dst[y*dstride + x + n] = lh;
      dst[(y + n)*dstride + x] = hl;
      dst[(y + n)*dstride + x + n] = hh;
    }
  }
}

static void od_tf_filter(od_coeff *dst, int dstride, int n) {
  od_coeff *v;
  od_coeff *u;
  int i;
  int m;
  m = (n >> 1) - 1;
  v = dst + dstride;
  for (i = 0; i < m; i++) {
    u = v;
    v += dstride << 1;
    *u += *v >> 1;
    *v -= *u >> 1;
  }
}

static void od_tf_filter_inv(od_coeff *dst, int dstride, int n) {
  od_coeff *v;
  od_coeff *u;
  int i;
  int m;
  u = dst + dstride*(n - 1);
  m = (n >> 1) - 1;
  for (i = 0; i < m; i++) {
    v = u;
    u -= dstride << 1;
    *v += *u >> 1;
    *u -= *v >> 1;
  }
}

void od_tf_filter_2d(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int n) {
  int i,j;
  if (dst!=src) {
    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        dst[dstride*j + i] = src[sstride*j + i];
      }
    }
  }
  for (j = 0; j < n; j++) {
    od_tf_filter(&dst[dstride*j], 1, n);
  }
  for (i = 0; i < n; i++) {
    od_tf_filter(&dst[i*1], dstride, n);
  }
}

void od_tf_filter_inv_2d(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int n) {
  int i,j;
  if (dst!=src) {
    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        dst[dstride*j + i] = src[sstride*j + i];
      }
    }
  }
  for (i = 0; i < n; i++) {
    od_tf_filter_inv(&dst[i*1], dstride, n);
  }
  for (j = 0; j < n; j++) {
    od_tf_filter_inv(&dst[dstride*j], 1, n);
  }
}

void od_convert_block_down(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int curr_size, int dest_size, int filter) {
  int n;
  int j;
  int i;
  od_coeff scratch[OD_BSIZE_MAX * OD_BSIZE_MAX];
  n = 1 << (curr_size + OD_LOG_BSIZE0);
  if (curr_size == dest_size) {
    if (dst != src) {
      for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
          dst[dstride*j + i] = src[sstride*j + i];
        }
      }
    }
    return;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      scratch[OD_BSIZE_MAX*j + i] = src[sstride*j + i];
    }
  }
  if (filter) {
    od_tf_filter_inv_2d(scratch, OD_BSIZE_MAX, scratch, OD_BSIZE_MAX, n);
  }
  if (curr_size - 1 == dest_size) {
    od_tf_down_hv(dst, dstride, scratch, OD_BSIZE_MAX, n);
  }
  else {
    od_coeff scratch2[OD_BSIZE_MAX * OD_BSIZE_MAX];
    od_tf_down_hv(scratch2, OD_BSIZE_MAX, scratch, OD_BSIZE_MAX, n);
    n >>= 1;
    for (j = 0; j < 2; j++) {
      for (i = 0; i < 2; i++) {
        od_convert_block_down(&dst[dstride*n*j + n*i], dstride,
         &scratch2[OD_BSIZE_MAX*n*j + n*i], OD_BSIZE_MAX, curr_size - 1,
         dest_size, filter);
      }
    }
  }
}
