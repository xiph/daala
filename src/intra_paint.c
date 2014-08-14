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
void pixel_interp(int pi[4], int pj[4], int w[4], int m, int i, int j,
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
  if (i == 0 || j == 0) {
    pi[0] = i;
    pj[0] = j;
    w[0] = 128;
    for (k = 1; k < 4; k++) pi[k] = pj[k] = w[k] = 0;
    return;
  }
  /* DC/Gradient mode, weights are proportional to 1/distance. The 255 alias
     is a temporary fix because mode is sometimes unsigned. */
  if (m == 4*n) {
    pi[0] = 0;
    pj[0] = j;
    pi[1] = n;
    pj[1] = j;
    pi[2] = i;
    pj[2] = 0;
    pi[3] = i;
    pj[3] = n;
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
}

/* Compute the actual interpolation (reconstruction) for a block. */
void interp_block(unsigned char *img, unsigned char *edge_accum, int n,
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

/* Predict the bottom edge using the top, left and (optionally, if already
   coded) right edge, using the mode signaled */
void predict_bottom_edge(int *p, unsigned char *edge_accum, int n, int stride,
 int m, int has_right) {
  int i;
  if (m == 2*n) {
    for(i = 0; i < n; i++) p[i] = edge_accum[(n - i - 1)*stride];
  }
  else if (m == 3*n) {
    for(i = 0; i < n; i++) p[i] = edge_accum[i + 1];
  }
  else if (has_right && m == 0) {
    for(i = 0; i < n; i++) p[i] = edge_accum[(i + 1)*stride + n];
  }
  else if (m > n && m < 2*n) {
    int slope;
    slope = m - n;
    for (i = 0; i < n; i++) {
      p[i] = edge_accum[(n - ((i + 1)*slope+n/2)/n)*stride];
    }
  }
  else if (m > 2*n && m < 3*n) {
    int slope;
    slope = 3*n - m;
    for (i = 0; i < slope; i++) {
      p[i] = edge_accum[(n - ((i + 1)*n+slope/2)/slope)*stride];
    }
    for (; i < n; i++) {
      p[i] = edge_accum[(i + 1 - slope)];
    }
  }
  else if (m > 0 && m < n) {
    int slope;
    slope = n - m;
    if (has_right) {
      for (i = 0; i < n; i++) {
        p[i] = edge_accum[n + (n - ((n - i - 1)*slope+n/2)/n)*stride];
      }
    }
    else {
      for (i=0; i < n; i++) {
        p[i] = edge_accum[n*stride];
      }
    }
  }
  else if (m > 3*n && m < 4*n) {
    int dir;
    dir = m - 3*n;
    for (i = 0; i < n - dir; i++) p[i] = edge_accum[i + 1 + dir];
    if (has_right) {
      for (; i < n; i++)
        p[i] = edge_accum[n + (n - ((n - i - 1)*n+dir/2)/dir)*stride];
    }
    else {
      OD_ASSERT(i != 0);
      for(; i < n; i++) {
        p[i] = p[i - 1];
      }
    }
  }
  else {
    if (has_right) {
      for(i = 0; i < n; i++) p[i] = edge_accum[n*stride] +
       (double)(i+1)/n*(edge_accum[n*stride+n]-edge_accum[n*stride]);
    }
    else {
      for(i = 0; i < n; i++) p[i] = edge_accum[n*stride] +
       .25*(double)(i+1)/n*(edge_accum[n]-edge_accum[n*stride]);
    }
  }
  if (has_right) p[n - 1] = edge_accum[n*stride+n];
}

/* Predict the right edge using the top, left and (optionally, if already
   coded) bottom edge, using the mode signaled */
void predict_right_edge(int *p, unsigned char *edge_accum, int n, int stride,
 int m, int has_bottom) {
  int i;
  int dir;
  dir = n - m;
  if (m > 3*n && m < 4*n) dir = 5*n - m;
  if (0 && m == 3*n) {
    for(i = 0; i < n; i++) p[i] = edge_accum[n];
  }
  else if (m == n) {
    for(i = 0; i < n; i++) p[i] = edge_accum[(i + 1)*stride];
  }
  else if (m == 2*n) {
    for(i = 0; i < n; i++) p[i] = edge_accum[n - i - 1];
  }
  else if (has_bottom && m == 0) {
    for(i = 0; i < n; i++) p[i] = edge_accum[n*stride + i + 1];
  }
  else if (m > n && m < 2*n) {
    int slope;
    int from_top;
    slope = m - n;
    from_top = m - n;
    for (i = 0; i < from_top; i++) {
      p[i] = edge_accum[n - ((i + 1)*n+slope/2)/slope];
    }
    for (; i < n; i++) {
      p[i] = edge_accum[(i + 1 - slope)*stride];
    }
  }
  else if (m > 2*n && m < 3*n) {
    for (i = 0; i < n; i++) {
      p[i] = edge_accum[n - ((i + 1)*(3*n-m)+n/2)/n];
    }
  }
  else if (m > 0 && m < n) {
    for (i = 0; i < n - dir; i++) p[i] = edge_accum[(i + 1 + dir)*stride];
    if (has_bottom) {
      for (; i < n; i++) {
        p[i] = edge_accum[n*stride + n - ((n - i - 1)*n+dir/2)/dir];
      }
    }
    else {
      OD_ASSERT(i != 0);
      for (; i < n; i++) {
        p[i] = p[i - 1];
      }
    }
  }
  else if (m > 3*n && m < 4*n) {
    int slope;
    slope = m-3*n;
    if (has_bottom) {
      for (i=0; i < n; i++) {
        p[i] = edge_accum[n*stride + n - ((n - i - 1)*slope)/n];
      }
    }
    else {
      for (i=0; i < n; i++) {
        p[i] = edge_accum[n];
      }
    }
  }
  else {
    if (has_bottom) {
      for(i = 0; i < n; i++) p[i] = edge_accum[n]
       + (double)(i+1)/n*(edge_accum[n*stride+n]-edge_accum[n]);
    }
    else {
      for(i = 0; i < n; i++) p[i] = edge_accum[n]
       + .25*(double)(i+1)/n*(edge_accum[n*stride]-edge_accum[n]);
    }
  }
  if (has_bottom) p[n - 1] = edge_accum[n*stride+n];
}

/* This is for simulating 32-point and 64-point DCTs since we don't have those
   at the moment. */
static void od_fdct32_approx(od_coeff *y, const od_coeff *x, int xstride) {
  int i;
  for (i = 0; i < 16; i++) {
    y[i] = M_SQRT1_2*(x[2*i*xstride] + x[(2*i + 1)*xstride]);
    y[16 + i] = M_SQRT1_2*(x[2*i*xstride] - x[(2*i + 1)*xstride]);
  }
  od_bin_fdct16(y, y, 1);
  od_bin_fdct16(y + 16, y + 16, 1);
}

static void od_idct32_approx(od_coeff *x, int xstride, const od_coeff *y) {
  int i;
  od_coeff tmp[32];
  od_bin_idct16(tmp, 1, y);
  od_bin_idct16(tmp + 16, 1, y + 16);
  for (i = 0; i < 16; i++) {
    x[2*i*xstride] = M_SQRT1_2*(tmp[i] + tmp[i + 16]);
    x[(2*i + 1)*xstride] = M_SQRT1_2*(tmp[i] - tmp[i + 16]);
  }
}

static void od_fdct64_approx(od_coeff *y, const od_coeff *x, int xstride) {
  int i;
  for (i = 0; i < 32; i++) {
    y[i] = M_SQRT1_2*(x[2*i*xstride] + x[(2*i + 1)*xstride]);
    y[32 + i] = M_SQRT1_2*(x[2*i*xstride] - x[(2*i + 1)*xstride]);
  }
  od_fdct32_approx(y, y, 1);
  od_fdct32_approx(y + 32, y + 32, 1);
}

static void od_idct64_approx(od_coeff *x, int xstride, const od_coeff *y) {
  int i;
  od_coeff tmp[64];
  od_idct32_approx(tmp, 1, y);
  od_idct32_approx(tmp + 32, 1, y + 32);
  for (i = 0; i < 32; i++) {
    x[2*i*xstride] = M_SQRT1_2*(tmp[i] + tmp[i + 32]);
    x[(2*i + 1)*xstride] = M_SQRT1_2*(tmp[i] - tmp[i + 32]);
  }
}

const od_fdct_func_1d my_fdct_table[5] = {
  od_bin_fdct4,
  od_bin_fdct8,
  od_bin_fdct16,
  od_fdct32_approx,
  od_fdct64_approx,
};

const od_idct_func_1d my_idct_table[5] = {
  od_bin_idct4,
  od_bin_idct8,
  od_bin_idct16,
  od_idct32_approx,
  od_idct64_approx,
};


int od_intra_paint_mode_cdf(ogg_uint16_t *cdf, int *dir_list, int *prob_list,
 int *ctx_list, int *dc_ctx, const unsigned char *mode, int bx, int by,
 int ln, int mstride, const unsigned char *dec8, int bstride, int res) {
  int i;
  int left;
  int top;
  int topleft;
  int idx;
  int nb;
  int bs;
  int cnt;
  int prob_sum;
  int norm;
  bs = ln - 2;
  nb = 4 << ln >> res;
  idx = (by*mstride + bx) << bs;
  top = left = topleft = nb;
  if (by > 0) {
    top = mode[idx - mstride] >> res;
    top = top << bs >> dec8[(((by<<bs)-1)>>1)*bstride + (bx<<bs>>1)];
  }
  if (bx > 0) {
    left = mode[idx - 1] >> res;
    left = left << bs >> dec8[(by<<bs>>1)*bstride + (((bx<<bs)-1)>>1)];
  }
  if (bx > 0 && by > 0) {
    topleft = mode[idx - mstride - 1] >> res;
    topleft = topleft << bs >> dec8[(((by<<bs)-1)>>1)*bstride + (((bx<<bs)-1)>>1)];
  }
  /* Compensate for mixed block size. */
  OD_ASSERT(topleft <= nb);
  OD_ASSERT(left <= nb);
  OD_ASSERT(top <= nb);
  *dc_ctx = 4*(top == nb) + 2*(left == nb) + (topleft == nb);

  cnt = 0;
  prob_sum = dir_prob[0];
  for (i = 0; i< nb; i++) {
    int ctx;
    ctx = 0;
    if (i == top) ctx+=2;
    if (abs(i - top) == 1) ctx++;
    if (i == left) ctx+=2;
    if (abs(i - left) == 1) ctx++;
    if (i == topleft) ctx+=2;
    if (abs(i - topleft) == 1) ctx++;
    if (ctx > 0) {
      dir_list[cnt] = i;
      ctx_list[cnt] = ctx;
      prob_list[cnt] = dir_prob[ctx];
      /*prob_list[cnt] = prob_list[cnt]*(256-dc_prob[dc_ctx]) >> 8;*/
      prob_sum += prob_list[cnt];
      cnt++;
    }
  }
  norm = 256*(256-dc_prob[*dc_ctx])/prob_sum;
  for (i = 0; i < cnt; i++) {
    prob_list[i] = prob_list[i]*norm >> 8;

  }
  cdf[0] = dir_prob[0]*norm >> 8;
  cdf[1] = cdf[0] + dc_prob[*dc_ctx];
  for (i = 0; i < cnt; i++) {
    cdf[i + 2] = cdf[i + 1] + prob_list[i];
  }
  prob_sum = cdf[cnt + 1];
  return cnt;
}

