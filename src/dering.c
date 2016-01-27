/*Daala video codec
Copyright (c) 2014-2016 Daala project contributors.  All rights reserved.

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

#include <stdlib.h>
#include "dering.h"

const od_filter_dering_direction_func
 OD_DERING_DIRECTION_C[OD_DERINGSIZES] = {
  od_filter_dering_direction_4x4_c,
  od_filter_dering_direction_8x8_c
};

const od_filter_dering_orthogonal_func
 OD_DERING_ORTHOGONAL_C[OD_DERINGSIZES] = {
  od_filter_dering_orthogonal_4x4_c,
  od_filter_dering_orthogonal_8x8_c
};

/* Generated from gen_filter_tables.c. */
const int OD_DIRECTION_OFFSETS_TABLE[8][3] = {
  {-1*OD_FILT_BSTRIDE + 1, -2*OD_FILT_BSTRIDE + 2, -3*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1, -1*OD_FILT_BSTRIDE + 2, -1*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1,  0*OD_FILT_BSTRIDE + 2,  0*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1,  1*OD_FILT_BSTRIDE + 2,  1*OD_FILT_BSTRIDE + 3  },
  { 1*OD_FILT_BSTRIDE + 1,  2*OD_FILT_BSTRIDE + 2,  3*OD_FILT_BSTRIDE + 3  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE + 1,  3*OD_FILT_BSTRIDE + 1  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE + 0,  3*OD_FILT_BSTRIDE + 0  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE - 1,  3*OD_FILT_BSTRIDE - 1  },
};

const double OD_DERING_GAIN_TABLE[OD_DERING_LEVELS] = {
  0, 0.5, 0.707, 1, 1.41, 2
};

/* Detect direction. 0 means 45-degree up-right, 2 is horizontal, and so on.
   The search minimizes the weighted variance along all the lines in a
   particular direction, i.e. the squared error between the input and a
   "predicted" block where each pixel is replaced by the average along a line
   in a particular direction. Since each direction have the same sum(x^2) term,
   that term is never computed. See Section 2, step 2, of:
   http://jmvalin.ca/notes/intra_paint.pdf */
static int od_dir_find8(const int16_t *img, int stride, int32_t *var) {
  int i;
  int cost[8] = {0};
  int partial[8][15] = {{0}};
  int best_cost = 0;
  int best_dir = 0;
  for (i = 0; i < 8; i++) {
    int j;
    for (j = 0; j < 8; j++) {
      int x;
      x = img[i*stride + j] >> OD_COEFF_SHIFT;
      partial[0][i + j] += x;
      partial[1][i + j/2] += x;
      partial[2][i] += x;
      partial[3][3 + i - j/2] += x;
      partial[4][7 + i - j] += x;
      partial[5][3 - i/2 + j] += x;
      partial[6][j] += x;
      partial[7][i/2 + j] += x;
    }
  }
  for (i = 0; i < 8; i++) {
    cost[2] += partial[2][i]*partial[2][i] >> 3;
    cost[6] += partial[6][i]*partial[6][i] >> 3;
  }
  for (i = 0; i < 7; i++) {
    cost[0] += OD_DIVU_SMALL(partial[0][i]*partial[0][i], i + 1)
     + OD_DIVU_SMALL(partial[0][14 - i]*partial[0][14 - i], i + 1);
    cost[4] += OD_DIVU_SMALL(partial[4][i]*partial[4][i], i + 1)
     + OD_DIVU_SMALL(partial[4][14 - i]*partial[4][14 - i], i + 1);
  }
  cost[0] += partial[0][7]*partial[0][8 - 1] >> 3;
  cost[4] += partial[4][7]*partial[4][8 - 1] >> 3;
  for (i = 1; i < 8; i += 2) {
    int j;
    for (j = 0; j < 4 + 1; j++) {
      cost[i] += partial[i][3 + j]*partial[i][3 + j] >> 3;
    }
    for (j = 0; j < 4 - 1; j++) {
      cost[i] += OD_DIVU_SMALL(partial[i][j]*partial[i][j], 2*j + 2)
       + OD_DIVU_SMALL(partial[i][10 - j]*partial[i][10 - j], 2*j + 2);
    }
  }
  for (i = 0; i < 8; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
  /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
  *var = best_cost - cost[(best_dir + 4) & 7];
  return best_dir;
}

#define OD_DERING_VERY_LARGE (30000)
#define OD_DERING_INBUF_SIZE ((OD_BSIZE_MAX + 2*OD_FILT_BORDER)*\
 (OD_BSIZE_MAX + 2*OD_FILT_BORDER))

/* Smooth in the direction detected. */
void od_filter_dering_direction_c(int16_t *y, int ystride, int16_t *in,
 int ln, int threshold, int dir) {
  int i;
  int j;
  int k;
  static const int taps[3] = {3, 2, 2};
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      int16_t sum;
      int16_t xx;
      int16_t yy;
      xx = in[i*OD_FILT_BSTRIDE + j];
      sum= 0;
      for (k = 0; k < 3; k++) {
        int16_t p0;
        int16_t p1;
        p0 = in[i*OD_FILT_BSTRIDE + j + OD_DIRECTION_OFFSETS_TABLE[dir][k]]
         - xx;
        p1 = in[i*OD_FILT_BSTRIDE + j - OD_DIRECTION_OFFSETS_TABLE[dir][k]]
         - xx;
        if (abs(p0) < threshold) sum += taps[k]*p0;
        if (abs(p1) < threshold) sum += taps[k]*p1;
      }
      yy = xx + ((sum + 8) >> 4);
      y[i*ystride + j] = yy;
    }
  }
}

void od_filter_dering_direction_4x4_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 2, threshold, dir);
}

void od_filter_dering_direction_8x8_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 3, threshold, dir);
}

/* Smooth in the direction orthogonal to what was detected. */
void od_filter_dering_orthogonal_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int ln, int threshold, int dir) {
  int i;
  int j;
  int offset;
  if (dir <= 4) offset = OD_FILT_BSTRIDE;
  else offset = 1;
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      int16_t athresh;
      int16_t yy;
      int16_t sum;
      int16_t p;
      /* Deringing orthogonal to the direction uses a tighter threshold
         because we want to be conservative. We've presumably already
         achieved some deringing, so the amount of change is expected
         to be low. Also, since we might be filtering across an edge, we
         want to make sure not to blur it. That being said, we might want
         to be a little bit more aggressive on pure horizontal/vertical
         since the ringing there tends to be directional, so it doesn't
         get removed by the directional filtering. */
      athresh = OD_MINI(threshold, threshold/3
       + abs(in[i*OD_FILT_BSTRIDE + j] - x[i*xstride + j]));
      yy = in[i*OD_FILT_BSTRIDE + j];
      sum = 0;
      p = in[i*OD_FILT_BSTRIDE + j + offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j + 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      y[i*ystride + j] = yy + ((3*sum + 8) >> 4);
    }
  }
}

void od_filter_dering_orthogonal_4x4_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 2, threshold, dir);
}

void od_filter_dering_orthogonal_8x8_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 3, threshold, dir);
}

/* This table approximates x^0.16 with the index being log2(x). It is clamped
   to [-.5, 3]. The table is computed as:
   round(256*min(3, max(.5, 1.08*(sqrt(2)*2.^([0:17]+8)/256/256).^.16))) */
static const int16_t OD_THRESH_TABLE_Q8[18] = {
  128, 134, 150, 168, 188, 210, 234, 262,
  292, 327, 365, 408, 455, 509, 569, 635,
  710, 768,
};

/* Compute deringing filter threshold for each 8x8 block based on the
   directional variance difference. A high variance difference means that we
   have a highly directional pattern (e.g. a high contrast edge), so we can
   apply more deringing. A low variance means that we either have a low
   contrast edge, or a non-directional texture, so we want to be careful not
   to blur. */
static void od_compute_thresh(int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int threshold, int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int32_t sb_var, int nhb, int nvb) {
  int bx;
  int by;
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      int v1;
      int v2;
      /* We use both the variance of 8x8 blocks and the variance of the
         entire superblock to determine the threshold. */
      v1 = OD_MINI(32767, var[by][bx] >> 6);
      v2 = OD_MINI(32767, sb_var/(OD_BSIZE_MAX*OD_BSIZE_MAX));
      thresh[by][bx] = threshold*OD_THRESH_TABLE_Q8[OD_CLAMPI(0,
       OD_ILOG(v1*v2) - 9, 17)] >> 8;
    }
  }
}

void od_dering(od_dering_opt_vtbl *vtbl, int16_t *y, int ystride, int16_t *x,
 int xstride, int ln, int sbx, int sby, int nhsb, int nvsb, int q, int xdec,
 int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int pli, unsigned char *bskip, int skip_stride, double gain) {
  int i;
  int j;
  int n;
  int threshold;
  int bx;
  int by;
  int16_t inbuf[OD_DERING_INBUF_SIZE];
  int16_t *in;
  int nhb;
  int nvb;
  int bsize;
  int varsum = 0;
  int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
  int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
  n = 1 << ln;
  bsize = 3 - xdec;
  nhb = nvb = n >> bsize;
  in = inbuf + OD_FILT_BORDER*OD_FILT_BSTRIDE + OD_FILT_BORDER;
  /* We avoid filtering the pixels for which some of the pixels to average
     are outside the frame. We could change the filter instead, but it would
     add special cases for any future vectorization. */
  for (i = 0; i < OD_DERING_INBUF_SIZE; i++) inbuf[i] = OD_DERING_VERY_LARGE;
  for (i = -OD_FILT_BORDER*(sby != 0); i < n
   + OD_FILT_BORDER*(sby != nvsb - 1); i++) {
    for (j = -OD_FILT_BORDER*(sbx != 0); j < n
     + OD_FILT_BORDER*(sbx != nhsb - 1); j++) {
      in[i*OD_FILT_BSTRIDE + j] = x[i*xstride + j];
    }
  }
  /* The threshold is meant to be the estimated amount of ringing for a given
     quantizer. Ringing is mostly proportional to the quantizer, but we
     use an exponent slightly smaller than unity because as quantization
     becomes coarser, the relative effect of quantization becomes slightly
     smaller as many unquantized coefficients are already close to zero. The
     value here comes from observing that on ntt-short, the best threshold for
     -v 5 appeared to be around 0.5*q, while the best threshold for -v 400
     was 0.25*q, i.e. 1-log(.5/.25)/log(400/5) = 0.84182 */
  threshold = gain*pow(q, 0.84182);
  if (pli == 0) {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        dir[by][bx] = od_dir_find8(&x[8*by*xstride + 8*bx], xstride,
         &var[by][bx]);
        varsum += var[by][bx];
      }
    }
    od_compute_thresh(thresh, threshold, var, varsum, nhb, nvb);
  }
  else {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        thresh[by][bx] = threshold;
      }
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      int xstart;
      int ystart;
      int xend;
      int yend;
      int skip;
      xstart = (sbx == 0) ? 0 : -1;
      ystart = (sby == 0) ? 0 : -1;
      xend = (2 >> xdec) + (sbx != nhsb - 1);
      yend = (2 >> xdec) + (sby != nvsb - 1);
      skip = 1;
      /* We look at whether the current block and its 4x4 surrounding (due to
         lapping) are skipped to avoid filtering the same content multiple
         times. */
      for (i = ystart; i < yend; i++) {
        for (j = xstart; j < xend; j++) {
          skip = skip && bskip[((by << 1 >> xdec) + i)*skip_stride
           + (bx << 1 >> xdec) + j];
        }
      }
      if (skip) thresh[by][bx] = 0;
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      (vtbl->filter_dering_direction[bsize - OD_LOG_BSIZE0])(
       &y[(by*ystride << bsize) + (bx << bsize)], ystride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       thresh[by][bx], dir[by][bx]);
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      in[i*OD_FILT_BSTRIDE + j] = y[i*ystride + j];
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      (vtbl->filter_dering_orthogonal[bsize - OD_LOG_BSIZE0])(
       &y[(by*ystride << bsize) + (bx << bsize)], ystride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       &x[(by*xstride << bsize) + (bx << bsize)], xstride,
       thresh[by][bx], dir[by][bx]);
    }
  }
}
