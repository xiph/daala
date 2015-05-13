/*Daala video codec
Copyright (c) 2014 Daala project contributors.  All rights reserved.

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

#include "dct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "odintrin.h"
#include "state.h"
#include "intra_paint.h"

/* All modes are numbered clockwise starting from mode 0 oriented 45 degrees
   up-right. For an NxN block, mode N is horizontal, mode 2*N is 45-degrees
   down-right, mode 3*N is vertical, and mode 4*N-1 is just next to mode 0.
   Mode 4*N means DC/gradient.  */

/* This function computes the position of the four points used to interpolate
   pixel (i,j) within an block of size 2^ln. The weights w[] are in Q7. We
   always use two edges when interpolating so we can make a smooth transition
   across edges. */
static void pixel_interp(int pi[4], int pj[4], int w[4], int m, int i, int j,
 int ln) {
  int k;
  int n;
  int rev;
  int r;
  int y0;
  int y1;
  int x0;
  int x1;
  int f0;
  int f1;
  int d0;
  int d1;
  int dir;
  n = 1 << ln;
  /* The the top and left edges, we reuse the points as-is. */
  if (i == n - 1 || j == n - 1) {
    pi[0] = i;
    pj[0] = j;
    w[0] = 128;
    for (k = 1; k < 4; k++) pi[k] = pj[k] = w[k] = 0;
    if (m == 3*n && i < n - 1) {
      pi[1] = n - 1;
      pj[1] = j;
      pi[2] = -1;
      pj[2] = j;
      w[0] = 32;
      w[1] = 96*(i + 1)/n;
      w[2] = 96 - w[1];
    }
    if (m == n && j < n - 1) {
      pi[1] = i;
      pj[1] = n - 1;
      pi[2] = i;
      pj[2] = -1;
      w[0] = 32;
      w[1] = 96*(j + 1)/n;
      w[2] = 96 - w[1];
    }
    return;
  }
  i++;
  j++;
  /* DC/Gradient mode, weights are proportional to 1/distance. The 255 alias
     is a temporary fix because mode is sometimes unsigned. */
  if (m == 4*n) {
    pi[0] = -1;
    pj[0] = j - 1;
    pi[1] = n - 1;
    pj[1] = j - 1;
    pi[2] = i - 1;
    pj[2] = -1;
    pi[3] = i - 1;
    pj[3] = n - 1;
    if (0) {
    w[0] = (n - i) << 6 >> ln;
    w[1] = i << 6 >> ln;
    w[2] = (n - j) << 6 >> ln;
    w[3] = j << 6 >> ln;
    }
    else {
      int sum;
      w[0] = 1024/i;
      w[1] = 1024/(n - i);
      w[2] = 1024/j;
      w[3] = 1024/(n - j);
      sum = w[0]+w[1]+w[2]+w[3];
      w[0] = w[0]*128/sum;
      w[1] = w[1]*128/sum;
      w[2] = w[2]*128/sum;
      w[3] = 128-w[0]-w[1]-w[2];
    }
    return;
  }
  if (m > 2*n) {
    int tmp;
    tmp = i;
    i = j;
    j = tmp;
    m = 4*n - m;
    rev = 1;
  }
  else {
    rev = 0;
  }
  dir = n - m;
  r = dir << 7 >> ln;
  y0 = (i << 7) + j*r;
  if (y0 >= 0 && y0 < (n << 7)) {
    pi[0] = y0 >> 7;
    f0 = y0 & 0x7f;
    pi[1] = OD_MINI(n, pi[0] + 1);
    pj[0] = pj[1] = 0;
    d0 = j*sqrt(128*128 + r*r);
  }
  else {
    int r_1;
    r_1 = (1 << 7 << ln) / dir;
    if (dir > 0) {
      x0 = (j << 7) - (n - i)*r_1;
      pi[0] = pi[1] = n;
      d0 = n - i;
    }
    else {
      x0 = (j << 7) + i*r_1;
      pi[0] = pi[1] = 0;
      d0 = i;
    }
    d0 = d0 * sqrt(128*128 + r_1*r_1);
    pj[0] = x0 >> 7;
    f0 = x0 & 0x7f;
    pj[1] = OD_MINI(n, pj[0] + 1);
  }

  y1 = (i << 7) - (n - j)*r;
  if (y1 >= 0 && y1 < (n << 7)) {
    pi[2] = y1 >> 7;
    f1 = y1 & 0x7f;
    pi[3] = OD_MINI(n, pi[2] + 1);
    pj[2] = pj[3] = n;
    d1 = (n - j)*sqrt(128*128 + r*r);
  }
  else {
    int r_1;
    r_1 = (1 << 7 << ln) / dir;
    if (dir > 0) {
      x1 = (j << 7) + i*r_1;
      pi[2] = pi[3] = 0;
      d1 = i;
    }
    else {
      x1 = (j << 7) - (n - i)*r_1;
      pi[2] = pi[3] = n;
      d1 = n - i;
    }
    d1 = d1 * sqrt(128*128 + r_1*r_1);
    pj[2] = x1 >> 7;
    f1 = x1 & 0x7f;
    pj[3] = OD_MINI(n, pj[2] + 1);
  }
  if (1) {
    w[0] = (128-f0)*d1/(d0+d1);
    w[1] = (f0)*d1/(d0+d1);
    w[2] = (128-f1)*d0/(d0+d1);
    w[3] = 128-w[0]-w[1]-w[2];
  }
  else {
    /* Pseudo-hanning blending -- doesn't seem to help. */
    double h = (double)d1/(double)(d0+d1);
    h = .5 - .5*cos(2*M_PI*h);
    w[0] = (128-f0)*h;
    w[1] = (f0)*h;
    w[2] = (128-f1)*(1-h);
    w[3] = 128-w[0]-w[1]-w[2];
  }
  if (rev) {
    for (k = 0; k < 4; k++) {
      int tmp = pi[k];
      pi[k] = pj[k];
      pj[k] = tmp;
    }
  }
  for (k = 0; k < 4; k++) pi[k]--;
  for (k = 0; k < 4; k++) pj[k]--;
}

/* Compute the actual interpolation (reconstruction) for a block. */
static void interp_block(unsigned char *img, const unsigned char *edge_accum, int n,
 int stride, int m) {
  int i;
  int j;
  int pi[4];
  int pj[4];
  int w[4];
  int ln;
  ln = 0;
  while (1 << ln < n) ln++;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      int k;
      int sum;
      pixel_interp(pi, pj, w, m, i, j, ln);
      sum = 0;
      for (k = 0; k < 4; k++) sum += edge_accum[pi[k]*stride + pj[k]]*w[k];
      img[i*stride + j] = OD_CLAMPI(0, (sum + 64) >> 7, 255);
    }
  }
}

/* This function finds the optimal paint direction by looking at the direction
   that minimizes the sum of the squared error caused by replacing each pixel
   of a line by the line average. For each line we would normally compute
   this as sum(x^2) - sum(x)^2/N, but since over the entire block the sum(x^2)
   terms summed over all lines will be the same regardless of the direction,
   then we don't need to compute that term and we can simply maximize
   sum(x^2)/N summed over all lines for a given direction. */
static int mode_select8b(const unsigned char *img, int n, int stride) {
  int i;
  int cost[9] = {0};
  int partial[9][2*MAXN + 1] = {{0}};
  int best_cost = 0;
  int best_dir = 0;
  for (i = 0; i < n; i++) {
    int j;
    for (j = 0; j < n; j++) {
      int x;
      x = img[i*stride + j] - 128;
      partial[0][i + j] += x;
      partial[1][i + j/2] += x;
      partial[2][i] += x;
      partial[3][n/2 - 1 + i - j/2] += x;
      partial[4][n - 1 + i - j] += x;
      partial[5][n/2 - 1 - i/2 + j] += x;
      partial[6][j] += x;
      partial[7][i/2 + j] += x;
      partial[8][0] += x;
    }
  }
  for (i = 0; i < n; i++) {
    cost[2] += partial[2][i]*partial[2][i]/n;
    cost[6] += partial[6][i]*partial[6][i]/n;
  }
  for (i = 0; i < n - 1; i++) {
    cost[0] += partial[0][i]*partial[0][i]/(i + 1)
     + partial[0][2*n - 2 - i]*partial[0][2*n - 2 - i]/(i + 1);
    cost[4] += partial[4][i]*partial[4][i]/(i + 1)
     + partial[4][2*n - 2 - i]*partial[4][2*n - 2 - i]/(i + 1);
  }
  cost[0] += partial[0][n - 1]*partial[0][n - 1]/n;
  cost[4] += partial[4][n - 1]*partial[4][n - 1]/n;
  for (i = 1; i < 8; i+=2) {
    int j;
    for (j = 0; j < n/2 + 1; j++) {
      cost[i] += partial[i][n/2 - 1 + j]*partial[i][n/2 - 1 + j]/n;
    }
    for (j = 0; j < n/2 - 1; j++) {
      cost[i] += partial[i][j]*partial[i][j]/(2*j+2);
      cost[i] += partial[i][3*n/2 - j]*partial[i][3*n/2 - j]/(2*j+2);
    }
  }
  cost[8] = (partial[8][0]/n)*(partial[8][0]/n) + 2*n*n;
  for (i = 0; i < 9; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
  /* Convert to a max of 4*N. */
  return best_dir*(4*n/8);
}

/* Accumulate sums on block edges so that we can later compute pixel values.
   There's usually two blocks used for each edge, but there can be up to 4
   in the corners. */
static void compute_edges(const unsigned char *img, int stride,
 int *edge_accum, int *edge_accum2, int *edge_count, int edge_stride, int n, int mode) {
  int i;
  int j;
  int pi[4];
  int pj[4];
  int w[4];
  int ln;
  ln = 0;
  while (1 << ln < n) ln++;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      int k;
      pixel_interp(pi, pj, w, mode, i, j, ln);
      for (k = 0; k < 4; k++) {
        edge_accum[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*w[k];
        edge_accum2[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*
          img[i*stride+j]*w[k];
        edge_count[pi[k]*edge_stride+pj[k]] += w[k];
      }
    }
  }
}

void od_intra_paint_analysis(const unsigned char *paint, int stride,
 const unsigned char *dec8, int bstride, unsigned char *mode, int mstride,
 int *edge_sum, int *edge_sum2, int *edge_count, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_intra_paint_analysis(paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_sum2, edge_count, bx, by, level);
    od_intra_paint_analysis(paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_sum2, edge_count, bx + 1, by, level);
    od_intra_paint_analysis(paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_sum2, edge_count, bx, by + 1, level);
    od_intra_paint_analysis(paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_sum2, edge_count, bx + 1, by + 1, level);
  }
  else {
    int curr_mode;
    int ln;
    int n;
    int k;
    int m;
    ln = 2 + bs;
    n = 1 << ln;
    curr_mode = mode_select8b(&paint[stride*n*by + n*bx], n, stride);
    for (k=0;k<1<<bs;k++) {
      for (m=0;m<1<<bs;m++) {
        mode[((by << bs) + k)*mstride + (bx << bs) + m] = curr_mode;
      }
    }
    compute_edges(&paint[stride*n*by + n*bx], stride,
     &edge_sum[stride*n*by + n*bx], &edge_sum2[stride*n*by + n*bx],
     &edge_count[stride*n*by + n*bx], stride, n, curr_mode);
  }
}

/* Computes the Wiener filter gain in Q8 considering (1/64)*q^2/12 as the
   noise. The 1/64 factor factor takes into consideration the fact that
   q^2/12 is really the worst case noise estimate and the fact that we'll be
   multiplying the filter gain by up to 16 when coding the strength of the
   deringing. */
#define VAR2(q) \
  do {int yy;\
  yy = ((double)edge_sum2[idx] \
   - (double)edge_sum1[idx]*edge_sum1[idx]/edge_count[idx])/edge_count[idx]; \
  paint_gain[idx] = OD_CLAMPI(0, (int)(256.*((q)*(q)/12./64)/(10+yy)), 255);} \
  while(0)

void od_paint_compute_edge_mask(od_adapt_ctx *adapt,
 unsigned char *paint, const unsigned char *img, unsigned char *paint_gain,
 int stride, const unsigned char *dec8, int bstride, unsigned char *mode,
 int mstride, int *edge_sum1, int *edge_sum2, int *edge_count, int q, int bx,
 int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_compute_edge_mask(adapt, paint, img, paint_gain, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, q, bx, by, level);
    od_paint_compute_edge_mask(adapt, paint, img, paint_gain, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, q, bx + 1, by, level);
    od_paint_compute_edge_mask(adapt, paint, img, paint_gain, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, q, bx, by + 1, level);
    od_paint_compute_edge_mask(adapt, paint, img, paint_gain, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, q, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int k;
    int idx;
    ln = 2 + bs;
    n = 1 << ln;
    if (bx == 0 && by == 0) {
      idx = -stride - 1;
      if (edge_count[idx] > 0) {
        paint[idx] = edge_sum1[idx]/edge_count[idx];
        VAR2(q);
      }
      else {
        paint[idx] = img[idx];
      }
    }
    /* Compute left edge (left column only). */
    if (bx == 0) {
      for (k = 0; k < n; k++) {
        idx = stride*(n*by + k) - 1;
        if (edge_count[idx] > 0) {
          paint[idx] = edge_sum1[idx]/edge_count[idx];
          VAR2(q);
        }
        else {
          paint[idx] = img[idx];
        }
      }
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      for (k = 0; k < n; k++) {
        idx = -stride + n*bx + k;
        if (edge_count[idx] > 0) {
          paint[idx] = edge_sum1[idx]/edge_count[idx];
          VAR2(q);
        }
        else {
          paint[idx] = img[idx];
        }
      }
    }
    /* Compute right edge stats. */
    for (k = 0; k < n - 1; k++) {
      idx = stride*(n*by + k) + n*(bx + 1) - 1;
      if (edge_count[idx] > 0) {
        paint[idx] = edge_sum1[idx]/edge_count[idx];
        VAR2(q);
      }
      else {
        paint[idx] = img[idx];
      }
    }
    /* Compute bottom edge stats. */
    for (k = 0; k < n; k++) {
      idx = stride*(n*(by + 1) - 1) + n*bx + k;
      if (edge_count[idx] > 0) {
        paint[idx] = edge_sum1[idx]/edge_count[idx];
        VAR2(q);
      }
      else {
        paint[idx] = img[idx];
      }
    }
    interp_block(&paint[stride*n*by + n*bx], &paint[stride*n*by + n*bx],
     n, stride, mode[(by*mstride + bx) << ln >> 2]);
    interp_block(&paint_gain[stride*n*by + n*bx], &paint_gain[stride*n*by + n*bx],
     n, stride, mode[(by*mstride + bx) << ln >> 2]);
  }
}

static unsigned char mask[1<<24];
static unsigned char paint_buf[1<<24];
unsigned char *paint_mask=mask+4096;
unsigned char *paint_out=paint_buf+4096;

ogg_uint16_t gain_cdf[] = {128, 256, 384, 512, 640, 768, 896, 1024, 1152};
