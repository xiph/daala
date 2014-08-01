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

#define QUANTIZE (1)

/* Warning, this will fail for images larger than 2048 x 2048 */
#define MAX_VAR_BLOCKS 1024
#define MAXN 64

static int edge_sum[1<<24] = {0};
static int edge_count[1<<24] = {0};

/* Throughout this code, mode -1 means DC/gradient. All other modes are
   numbered clockwise starting from mode 0 oriented 45 degrees up-right.
   For an NxN block, mode N is horizontal, mode 2*N is 45-degrees down-right,
   mode 3*N is vertical, and mode 4*N-1 is just next to mode 0. */

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
  if (i == 0 || j == 0) {
    pi[0] = i;
    pj[0] = j;
    w[0] = 128;
    for (k = 1; k < 4; k++) pi[k] = pj[k] = w[k] = 0;
    return;
  }
  /* DC/Gradient mode, weights are proportional to 1/distance. */
  if (m == -1) {
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

static void compare_mode(unsigned char block[MAXN][MAXN],
 unsigned char best_block[MAXN][MAXN], int *dist, int *best_dist, int id,
 int *best_id, int n, const unsigned char *img, int stride) {
  int i;
  int j;
  int curr_dist;
  curr_dist = 0;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      int e;
      e = (int)block[i][j] - (int)img[i*stride + j];
      curr_dist += e*e;
    }
  }
#if 1
  /* Give a slight bias to DC/gradient mode otherwise it doesn't get used.
     This needs to be improved. */
  if (id==-1) curr_dist -= n*n*2 + 0*curr_dist/4;
#endif
  *dist = curr_dist;
  if (curr_dist < *best_dist) {
    *best_dist = curr_dist;
    *best_id = id;
    for (i=0;i<n;i++) for (j=0;j<n;j++) best_block[i][j] = block[i][j];
  }
}

/* Select which mode to use for a block by making the (false) assumption that
   the edge is coded only based on that mode. */
static int mode_select(const unsigned char *img, int *dist, int n, int stride) {
  int i;
  int j;
  int m;
  int best_dist;
  int best_id;
  int edge_accum[MAXN+1][MAXN+1];
  int edge_count[MAXN+1][MAXN+1];
  unsigned char block[MAXN][MAXN];
  unsigned char best_block[MAXN][MAXN];
  int pi[4];
  int pj[4];
  int w[4];
  int ln;
  best_dist = 1<<30;
  best_id = 0;
  ln = 0;
  while (1 << ln < n) ln++;
  for (m = -1; m < 4*n; m++) {
    /* Only consider half the directions. */
    if (m>=0 && (m&1)) continue;
    int dist;
    for (i = 0; i <= n; i++) for (j = 0; j <= n; j++) edge_accum[i][j] = 0;
    for (i = 0; i <= n; i++) for (j = 0; j <= n; j++) edge_count[i][j] = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        int k;
        pixel_interp(pi, pj, w, m, i, j, ln);
        for (k = 0; k < 4; k++) {
          edge_accum[pi[k]][pj[k]] += (int)img[i*stride+j]*w[k];
          edge_count[pi[k]][pj[k]] += w[k];
        }
      }
    }
    for (i = 0; i <= n; i++) {
      for (j = 0; j <= n; j++) {
        if (edge_count[i][j] > 0) {
          edge_accum[i][j] = edge_accum[i][j]/edge_count[i][j];
        }
      }
    }

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        int k;
        int sum;
        sum = 0;
        pixel_interp(pi, pj, w, m, i, j, ln);
        for (k = 0; k < 4; k++) sum += edge_accum[pi[k]][pj[k]]*w[k];
        block[i][j] = OD_CLAMPI(0, (sum + 64) >> 7, 255);
      }
    }
    compare_mode(block, best_block, &dist, &best_dist, m, &best_id, n, img,
     stride);
  }
  if (dist) *dist = best_dist;
  return best_id;
}

/* Compute the final edges once the contribution of all blocks are counted.
   There's usually two blocks used for each edge, but there can be up to 4
   in the corners. */
static void compute_edges(const unsigned char *img, int *edge_accum,
 int *edge_count, int n, int stride, int m) {
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
      pixel_interp(pi, pj, w, m, i, j, ln);
      for (k = 0; k < 4; k++) {
        edge_accum[pi[k]*stride+pj[k]] += (int)img[i*stride+j]*w[k];
        edge_count[pi[k]*stride+pj[k]] += w[k];
      }
    }
  }
}

/* Compute the actual interpolation (reconstruction) for a block. */
static void interp_block(unsigned char *img, unsigned char *edge_accum, int n,
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
static void predict_bottom_edge(int *p, unsigned char *edge_accum, int n, int stride,
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
    int from_left;
    int dir;
    dir = m - 3*n;
    if (m > 2*n) from_left = 3*n - m;
    else from_left = n;
    for (i = 0; i < from_left; i++) {
      p[i] = edge_accum[(n + ((i + 1)*n+dir/2)/dir)*stride];
    }
    for (; i < n; i++) {
      p[i] = edge_accum[(i + 1 + dir)];
    }
  }
  else if (m > 0 && m < n) {
    int slope;
    slope = n - m;
    if (has_right)
    {
      for (i = 0; i < n; i++)
        p[i] = edge_accum[n + (n - ((n - i - 1)*slope+n/2)/n)*stride];
    } else {
      for (i = 0; i < n; i++) {
        for(i = 0; i < n; i++) p[i] = edge_accum[n*stride] +
         .25*(double)(i+1)/n*(edge_accum[n]-edge_accum[n*stride]);
      }
    }
  }
  else if (m > 3*n && m < 4*n) {
    int dir;
    dir = m - 3*n;
    for (i = 0; i < n - dir; i++) p[i] = edge_accum[i + 1 + dir];
    for (; i < n; i++)
      p[i] = edge_accum[n + (n - ((n - i - 1)*n+dir/2)/dir)*stride];
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
static void predict_right_edge(int *p, unsigned char *edge_accum, int n, int stride,
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

static const od_fdct_func_1d my_fdct_table[5] = {
  od_bin_fdct4,
  od_bin_fdct8,
  od_bin_fdct16,
  od_fdct32_approx,
  od_fdct64_approx,
};

static const od_idct_func_1d my_idct_table[5] = {
  od_bin_idct4,
  od_bin_idct8,
  od_bin_idct16,
  od_idct32_approx,
  od_idct64_approx,
};

/* Quantize the bottom edge using prediction. */
static void quantize_bottom_edge(unsigned char *edge_accum, int n, int stride, int q,
 int m, int has_right) {
  int x[MAXN];
  int r[MAXN];
  int p[MAXN] = {0};
  int lsize;
  int i;
  if (n == 4) lsize = 0;
  else if (n == 8) lsize = 1;
  else if (n == 16) lsize = 2;
  else lsize = 3;

  /* Quantize bottom edge. */
  predict_bottom_edge(p, edge_accum, n, stride, m, has_right);
#if !QUANTIZE
  /*printf("%d ", m);*/
  for (i = 0; i < n; i++) printf("%d ", edge_accum[n*stride + (i + 1)] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[n*stride + i + 1] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) {
    x[i] = (int)(q*floor(.5+x[i]/q));
    /*printf("%d ", (int)floor(.5+x[i]/q));*/
  }
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[n*stride + i + 1] = OD_MAXI(0, OD_MINI(255, r[i]));
}

/* Quantize the right edge using prediction. */
static void quantize_right_edge(unsigned char *edge_accum, int n, int stride, int q,
 int m, int has_bottom) {
  int x[MAXN];
  int r[MAXN];
  int p[MAXN] = {0};
  int lsize;
  int i;
  if (n == 4) lsize = 0;
  else if (n == 8) lsize = 1;
  else if (n == 16) lsize = 2;
  else lsize = 3;

  /* Quantize right edge. */
  predict_right_edge(p, edge_accum, n, stride, m, has_bottom);
#if !QUANTIZE
  for (i = 0; i < n; i++) printf("%d ", edge_accum[(i + 1)*stride + n] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[(i + 1)*stride + n] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) {
    x[i] = (int)(q*floor(.5+x[i]/q));
    /*printf("%d ", (int)floor(.5+x[i]/q));*/
  }
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[(i+1)*stride+n] = OD_MAXI(0, OD_MINI(255, r[i]));
}

/* Quantize both the right and bottom edge, changing the order to maximize
   the number of pixels we can predict. */
static void quantize_edge(unsigned char *edge_accum, int n, int stride, int q, int m) {
  /*printf("q\%d ", n);*/
  if (m > 0 && m < n) {
    quantize_right_edge(edge_accum, n, stride, q, m, 0);
    quantize_bottom_edge(edge_accum, n, stride, q, m, 1);
  }
  else {
    quantize_bottom_edge(edge_accum, n, stride, q, m, 0);
    quantize_right_edge(edge_accum, n, stride, q, m, 1);
  }
  /*printf("\n");*/
}

void od_intra_paint_encode(unsigned char *paint, const unsigned char *img,
 int w8, int h8, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride) {
  int i, j;
  for(i = 8; i < 2*h8-8; i++) {
    for(j = 8; j < 2*w8-8; j++) {
      int bs;
      bs = dec8[(i>>1)*bstride + (j>>1)];
      if (i>>bs<<bs == i && j>>bs<<bs == j) {
        int k, m;
        mode[i*mstride + j] = mode_select(&img[4*stride*i + 4*j], NULL, 4<<bs, stride);
        compute_edges(&img[4*stride*i + 4*j], &edge_sum[4*stride*i + 4*j],
         &edge_count[4*stride*i + 4*j], 4<<bs, stride, mode[i*mstride + j]);
        for (k=0;k<1<<bs;k++) {
          for (m=0;m<1<<bs;m++) {
            mode[(i+k)*mstride + j + m] = mode[i*mstride + j];
          }
        }
      }
    }
  }

  for (i = 0; i < 1<<24; i++) {
    if (edge_count[i] > 0) {
      paint[i] = edge_sum[i]/edge_count[i];
    }
  }

  for(i = 8; i < 2*h8-8; i++) {
    for(j = 8; j < 2*w8-8; j++) {
      int bs;
      bs = dec8[(i>>1)*bstride + (j>>1)];
      if (i>>bs<<bs == i && j>>bs<<bs == j) {
        quantize_edge(&paint[4*stride*i + 4*j], 4<<bs, stride, 30,
         mode[i*mstride + j]);
        interp_block(&paint[4*stride*i + 4*j], &paint[4*stride*i + 4*j],
         4<<bs, stride, mode[i*mstride + j]);
      }
    }
  }
#if 0
  for(i = 0; i < 2*h8; i++) {
    for(j = 0; j < 2*w8; j++) {
      printf("%d ", mode[i][j]);
    }
    printf("\n");
  }
#endif
}

void od_intra_paint_choose_block_size(const unsigned char *img, int stride,
 int bsize[2][2]) {
  int i;
  int j;
  int cost16;
  int cost32;
  cost16 = 0;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int k;
      int m;
      int dist;
      int cost8;
      cost8 = 0;
      for (k = 0; k < 2; k++) {
        for (m = 0; m < 2; m++) {
          mode_select(&img[(16*i + 8*k)*stride + 16*j + 8*m], &dist, 8,
           stride);
          cost8 += dist;
        }
      }
      cost8 = cost8 + 16000;
      mode_select(&img[16*stride*i + 16*j], &dist, 16, stride);
      if (dist < cost8) {
        bsize[i][j] = 2;
      }
      else {
        dist = cost8;
        bsize[i][j] = 1;
      }
      cost16 += dist;
    }
  }
  mode_select(img, &cost32, 32, stride);
  if (cost32 < cost16+32000) {
    bsize[0][0] = bsize[0][1] = bsize[1][0] = bsize[1][1] = 3;
  }
}

int switch_decision(unsigned char *img, int w, int h, int stride, int ow, int oh)
{
  int i,j;
  int h8,w8,h32,w32;
  unsigned char *dec8;
  unsigned char *mode;
  int bstride;
  int mstride;
  unsigned char *paint_out;

  (void)ow;
  (void)oh;

  paint_out = (unsigned char*)malloc((w+1)*(h+1)*sizeof(*paint_out));
  w>>=1;
  h>>=1;
  w8 = w>>2;
  h8 = h>>2;
  w32 = w>>4;
  h32 = h>>4;
  dec8 = (unsigned char*)malloc(w8*h8*sizeof(*dec8));
  mode = (unsigned char*)malloc(w8*h8*sizeof(*mode)<<2);
  bstride = w8;
  mstride = w8<<1;
  /* Replace decision with the one from process_block_size32() */
  for(i=1;i<h32-1;i++){
    for(j=1;j<w32-1;j++){
      int k,m;
      int dec[2][2];
#if 0
      od_intra_paint_choose_block_size(img+32*stride*i+32*j, stride, dec);
#else
      dec[0][0] = dec[0][1] = dec[1][0] = dec[1][1] = 2;
#endif
      for(k=0;k<4;k++)
        for(m=0;m<4;m++)
          dec8[(4*i + k)*bstride + (4*j + m)]=dec[k>>1][m>>1];
    }
  }
  od_intra_paint_encode(paint_out, img, w8, h8, stride, dec8, bstride, mode,
   mstride);
  for(i=32;i<32*(h32-1);i++) {
    for(j=32;j<32*(w32-1);j++) {
      img[i*stride + j] = paint_out[i*stride + j];
    }
  }
#if 0
  for(i=4;i<h8-4;i++){
    for(j=4;j<w8-4;j++){
      if ((i&3)==0 && (j&3)==0){
        int k;
        for(k=0;k<32;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<32;k++)
          img[(8*i+k)*stride+j*8] = 0;
      }
      if ((i&1)==0 && (j&1)==0 && dec8[i*bstride + j]==2){
        int k;
        for(k=0;k<16;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<16;k++)
          img[(8*i+k)*stride+j*8] = 0;
      }
      if (dec8[i*bstride + j]<=1){
        int k;
        for(k=0;k<8;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<8;k++)
          img[(8*i+k)*stride+j*8] = 0;
        if (dec8[i*bstride + j]==0){
          img[(8*i+4)*stride+j*8+3] = 0;
          img[(8*i+4)*stride+j*8+4] = 0;
          img[(8*i+4)*stride+j*8+5] = 0;
          img[(8*i+3)*stride+j*8+4] = 0;
          img[(8*i+5)*stride+j*8+4] = 0;
        }
      }
    }
  }
  for (i=32;i<(w32-1)*32;i++)
    img[(h32-1)*32*stride+i]=0;
  for (i=32;i<(h32-1)*32;i++)
    img[i*stride+(w32-1)*32]=0;
#endif
  free(paint_out);
  free(dec8);
  free(mode);
  return 0;
}


