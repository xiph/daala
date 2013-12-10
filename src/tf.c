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

static void od_convert_block_up(od_coeff *dst, int dstride,
 const od_coeff *src, int sstride, const unsigned char *bsize,
 int bstride, int nx, int ny, int dest_size) {
  int i;
  int j;
  int n;
  od_coeff scratch[16][16];
  n = 1 << (dest_size + 2);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) scratch[i][j] = src[i*sstride + j];
  }
  if (dest_size > 1) {
    int sub;
    sub = 1 << (dest_size - 1);
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        if (OD_BLOCK_SIZE4x4(bsize, bstride, nx + i*sub, ny + j*sub)
         < dest_size - 1) {
          od_convert_block_up(&scratch[4*i*sub][4*j*sub], 16,
           &scratch[4*i*sub][4*j*sub], 16, bsize, bstride,
           nx + i*sub, ny + j*sub, dest_size - 1);
        }
      }
    }
  }
  /* At this point, we assume that scratch has size dest_size-1. */
  od_tf_up_hv(dst, dstride, &scratch[0][0], 16, 1 << (dest_size - 1 + 2));
}

void od_convert_block_down(od_coeff *dst, int dstride, const od_coeff *src,
 int sstride, int curr_size, int dest_size) {
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
  if (curr_size - 1 == dest_size) {
    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        scratch[OD_BSIZE_MAX*j + i] = src[sstride*j + i];
      }
    }
    od_tf_down_hv(dst, dstride, scratch, OD_BSIZE_MAX, n);
  }
  else {
    od_tf_down_hv(scratch, OD_BSIZE_MAX, src, sstride, n);
    n >>= 1;
    for (j = 0; j < 2; j++) {
      for (i = 0; i < 2; i++) {
        od_convert_block_down(&dst[dstride*n*j + n*i], dstride,
         &scratch[OD_BSIZE_MAX*n*j + n*i], OD_BSIZE_MAX, curr_size - 1,
         dest_size);
      }
    }
  }
}

void od_convert_intra_coeffs(od_coeff *(dst[4]), int dstrides[4],
 od_coeff *src, int sstride, int bx, int by,
 const unsigned char *bsize, int bstride, int has_ur) {
  /* Relative position of neighbors: up-left     up   "up-right"  left */
  static const int offsets[4][2] =
  { { -1, -1 }, { 0, -1 }, { 1, -1 }, { -1, 0 } };
  int csize;
  int n;
  csize = OD_BLOCK_SIZE4x4(bsize, bstride, bx, by);
#if defined(OD_VALGRIND)
  for (n = 0; n < 4; n++) {
    int z;
    for (z = 0; z < (4 << csize); z++) {
      VALGRIND_MAKE_MEM_UNDEFINED(dst[n] + z*dstrides[n],
       dstrides[n]*sizeof(*dst[n]));
    }
  }
#endif
  /* Loop over neighbours. */
  for (n = 0; n < 4; n++) {
    int nx;
    int ny;
    int nsize;
    if (n == 2 && !has_ur) nx = bx;
    else nx = bx + (offsets[n][0] << csize);
    ny = by + (offsets[n][1] << csize);
    nsize = OD_BLOCK_SIZE4x4(bsize, bstride, nx, ny);
    if (nsize == csize) {
      /* simply override the pointer and stride */
      dst[n] = src + 4*ny*sstride + 4*nx;
      dstrides[n] = sstride;
    }
    else if (nsize > csize) {
      /* We need to TF down. */
      int size;
      int i;
      int j;
      int off_x;
      int off_y;
      int nxa;
      int nya;
      od_coeff scratch[16][16];
      /* Aligns nx and ny to the size of the neighbour. */
      nxa = nx >> nsize << nsize;
      nya = ny >> nsize << nsize;
      size = 1 << csize << 2;
      od_convert_block_down(&scratch[0][0], 16, &src[4*nya*sstride + 4*nxa],
       sstride, nsize, csize);
      /* Find there offset in the TF'ed block that has the useful data */
      off_x = (nx-nxa) << 2;
      off_y = (ny-nya) << 2;
      /* Copy only the part we need */
      for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
          dst[n][i*dstrides[n] + j] = scratch[i + off_y][j + off_x];
        }
      }
    }
    else {
      /* We need to TF up. */
      od_convert_block_up(dst[n], dstrides[n], &src[4*ny*sstride + 4*nx],
       sstride, bsize, bstride, nx, ny, csize);
    }
#if defined(OD_VALGRIND)
    {
      int z;
      for (z = 0; z < (4 << csize); z++) {
        VALGRIND_CHECK_MEM_IS_DEFINED(dst[n] + z*dstrides[n],
         (4 << csize)*sizeof(*dst[n]));
      }
    }
#endif
  }
}
