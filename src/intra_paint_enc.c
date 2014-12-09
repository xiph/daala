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

#include <stdlib.h>

#include "entenc.h"
#include "generic_code.h"
#include "laplace_code.h"
#include "intra_paint.h"
#include "state.h"

static void compare_mode(unsigned char block[MAXN + 1][MAXN + 1],
 unsigned char best_block[MAXN][MAXN], int *dist, int *best_dist, int id,
 int *best_id, int n, const unsigned char *img, int stride) {
  int i;
  int j;
  int curr_dist;
  curr_dist = 0;
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      int e;
      e = (int)block[i+1][j+1] - (int)img[i*stride + j];
      curr_dist += e*e;
    }
  }
#if 1
  /* Give a slight bias to DC/gradient mode otherwise it doesn't get used.
     This needs to be improved. */
  if (id==4*n) curr_dist -= n*n*2 + 0*curr_dist/4;
#endif
  *dist = curr_dist;
  if (curr_dist < *best_dist) {
    *best_dist = curr_dist;
    *best_id = id;
    for (i=0;i<n;i++) for (j=0;j<n;j++) best_block[i][j] = block[i+1][j+1];
  }
}

/* Select which mode to use for a block by making the (false) assumption that
   the edge is coded only based on that mode. */
static int mode_select(const unsigned char *img, int *dist, int n, int stride,
 int res) {
  int i;
  int j;
  int m;
  int best_dist;
  int best_id;
  int edge_accum[MAXN+1][MAXN+1];
  int edge_count[MAXN+1][MAXN+1];
  unsigned char block[MAXN + 1][MAXN + 1];
  unsigned char best_block[MAXN][MAXN];
  int pi[4];
  int pj[4];
  int w[4];
  int ln;
  best_dist = 1<<30;
  best_id = 0;
  ln = 0;
  while (1 << ln < n) ln++;
  for (m = 0; m <= 4*n; m += 1 << res) {
    int dist;
    for (i = 0; i <= n; i++) for (j = 0; j <= n; j++) edge_accum[i][j] = 0;
    for (i = 0; i <= n; i++) for (j = 0; j <= n; j++) edge_count[i][j] = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        int k;
        pixel_interp(pi, pj, w, m, i, j, ln);
        for (k = 0; k < 4; k++) {
          /* Avoids having edges with no data. */
          w[k] = OD_MAXI(w[k], 1);
          edge_accum[1+pi[k]][1+pj[k]] += (int)img[i*stride+j]*w[k];
          edge_count[1+pi[k]][1+pj[k]] += w[k];
        }
      }
    }
    for (i = 0; i <= n; i++) {
      for (j = 0; j <= n; j++) {
        if (edge_count[i][j] > 0) {
          block[i][j] = edge_accum[i][j]/edge_count[i][j];
        }
      }
    }
    interp_block(&block[1][1], &block[1][1], n, MAXN + 1, m);
    compare_mode(block, best_block, &dist, &best_dist, m, &best_id, n, img,
     stride);
  }
  if (dist) *dist = best_dist;
  return best_id;
}

static int mode_select8b(const unsigned char *img, int *dist, int n, int stride,
 int res) {
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

/* Compute the final edges once the contribution of all blocks are counted.
   There's usually two blocks used for each edge, but there can be up to 4
   in the corners. */
static void compute_edges(const unsigned char *img, int stride,
 int *edge_accum, int *edge_count, int edge_stride, int n, int mode) {
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
        edge_count[pi[k]*edge_stride+pj[k]] += w[k];
      }
    }
  }
}

/* Compute edge variance. */
static void compute_edge_variance(const unsigned char *img, const unsigned char *paint, int stride,
 int *edge_orig, int *edge_accum2, int *edge_corr, int edge_stride, int n,
 int mode) {
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
        edge_orig[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*w[k];
        edge_accum2[pi[k]*edge_stride+pj[k]] += (int)paint[i*stride+j]*
          paint[i*stride+j]*w[k];
        edge_corr[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*paint[i*stride+j]*w[k];
      }
    }
  }
}



void od_paint_mode_select(const unsigned char *img, const unsigned char *paint,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int res,
 int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)] + 0;

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_mode_select(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by, level);
    od_paint_mode_select(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by, level);
    od_paint_mode_select(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by + 1, level);
    od_paint_mode_select(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int k;
    int m;
    int curr_mode;
    ln = 2 + bs;
    n = 1 << ln;
    curr_mode = mode_select8b(&img[stride*n*by + n*bx], NULL, n, stride, res);
    curr_mode >>= 0;
    for (k=0;k<1<<bs;k++) {
      for (m=0;m<1<<bs;m++) {
        mode[((by << bs) + k)*mstride + (bx << bs) + m] = curr_mode;
      }
    }
  }
}

void od_intra_paint_analysis(const unsigned char *img, const unsigned char *paint,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int res,
 int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_intra_paint_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by, level);
    od_intra_paint_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by, level);
    od_intra_paint_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by + 1, level);
    od_intra_paint_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    ln = 2 + bs;
    n = 1 << ln;
    compute_edges(&paint[stride*n*by + n*bx], stride,
     &edge_sum[stride*n*by + n*bx], &edge_count[stride*n*by + n*bx], stride,
     n, mode[(by<<bs)*mstride + (bx<<bs)]);
  }
}

void od_paint_var_analysis(const unsigned char *img, const unsigned char *paint,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum1, int *edge_sum2, int *edge_count, int res,
 int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_var_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by, level);
    od_paint_var_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by, level);
    od_paint_var_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by + 1, level);
    od_paint_var_analysis(img, paint, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int curr_mode;
    ln = 2 + bs;
    n = 1 << ln;
    curr_mode = mode[(by<<bs)*mstride + (bx<<bs)];
    compute_edge_variance(&img[stride*n*by + n*bx],
     &paint[stride*n*by + n*bx], stride,
     &edge_sum1[stride*n*by + n*bx], &edge_sum2[stride*n*by + n*bx],
     &edge_count[stride*n*by + n*bx], stride, n, curr_mode);
  }
}


void od_intra_paint_compute_edges(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count,
 int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by + 1, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by + 1, level);
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
      if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
      else paint[idx] = img[idx];
    }
    /* Compute left edge (left column only). */
    if (bx == 0) {
      for (k = 0; k < n; k++) {
        idx = stride*(n*by + k) - 1;
        if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
        else paint[idx] = img[idx];
      }
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      for (k = 0; k < n; k++) {
        idx = -stride + n*bx + k;
        if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
        else paint[idx] = img[idx];
      }
    }
    /* Compute right edge stats. */
    for (k = 0; k < n - 1; k++) {
      idx = stride*(n*by + k) + n*(bx + 1) - 1;
      if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
      else paint[idx] = img[idx];
    }
    /* Compute bottom edge stats. */
    for (k = 0; k < n; k++) {
      idx = stride*(n*(by + 1) - 1) + n*bx + k;
      if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
      else paint[idx] = img[idx];
    }
#if 0
    /* Refinement: one last chance to pick DC. */
    if (mode[(by*mstride + bx) << ln >> 2] != 4*n) {
      int i;
      int j;
      unsigned char block[MAXN + 1][MAXN + 1];
      unsigned char best_block[MAXN][MAXN];
      int dist;
      int best_dist;
      int best_id = 0;
      int m;
      m = mode[(by*mstride + bx) << ln >> 2];
      best_dist = 1 << 30;
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
          block[i][j] = paint[stride*(n*by + i - 1) + n*bx + j - 1];
        }
      }
      interp_block(&block[1][1], &block[1][1], n, MAXN + 1, m);
      compare_mode(block, best_block, &dist, &best_dist, m, &best_id, n,
       &img[stride*n*by + n*bx], stride);
      interp_block(&block[1][1], &block[1][1], n, MAXN + 1, 4*n);
      compare_mode(block, best_block, &dist, &best_dist, 4*n, &best_id, n,
       &img[stride*n*by + n*bx], stride);
      mode[(by*mstride + bx) << ln >> 2] = best_id;
    }
#endif
  }
}

#define SQUARE(x) ((x)*(x))
#define COVAR(xy,x,y,c) ((double)(xy)/(c) - ((double)(x)*(y)/((c)*(c))))

#define VAR(q) do {paint[idx] = OD_MINI(255, (int)(256./12/64.*(q)*(q)/(1+(double)edge_sum2[idx]/edge_count[idx] \
  - SQUARE((double)edge_sum1[idx]/edge_count[idx]))));} while(0)

#if 1
#define VAR2(q) do {int yy;\
                yy = COVAR(edge_sum2[idx], edge_sum1[idx], edge_sum1[idx], edge_count[idx]); \
                paint[idx] = OD_CLAMPI(0, (int)(256.*((q)*(q)/12./64)/(10+yy)), 255);} while(0)
#else
#define VAR2(q) do {int yy, xy;\
                yy = COVAR(edge_sum2[idx], edge_sum1[idx], edge_sum1[idx], edge_count[idx]); \
                xy = COVAR(edge_corr[idx], orig_edge_sum[idx], edge_sum1[idx], edge_count[idx]); \
                paint[idx] = OD_CLAMPI(0, (int)(256.*(yy-xy)/(10+yy)), 255);} while(0)
#endif

void od_paint_compute_edge_mask(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum1, int *edge_sum2, int *orig_edge_sum, int *edge_corr, int *edge_count,
 int q, int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, orig_edge_sum, edge_corr, edge_count, q, res, bx, by, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, orig_edge_sum, edge_corr, edge_count, q, res, bx + 1, by, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, orig_edge_sum, edge_corr, edge_count, q, res, bx, by + 1, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, orig_edge_sum, edge_corr, edge_count, q, res, bx + 1, by + 1, level);
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
      if (edge_count[idx] > 0) VAR2(q);
      else paint[idx] = img[idx];
    }
    /* Compute left edge (left column only). */
    if (bx == 0) {
      for (k = 0; k < n; k++) {
        idx = stride*(n*by + k) - 1;
        if (edge_count[idx] > 0) VAR2(q);
        else paint[idx] = img[idx];
      }
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      for (k = 0; k < n; k++) {
        idx = -stride + n*bx + k;
        if (edge_count[idx] > 0) VAR2(q);
        else paint[idx] = img[idx];
      }
    }
    /* Compute right edge stats. */
    for (k = 0; k < n - 1; k++) {
      idx = stride*(n*by + k) + n*(bx + 1) - 1;
      if (edge_count[idx] > 0) VAR2(q);
      else paint[idx] = img[idx];
    }
    /* Compute bottom edge stats. */
    for (k = 0; k < n; k++) {
      idx = stride*(n*(by + 1) - 1) + n*bx + k;
      if (edge_count[idx] > 0) VAR2(q);
      else paint[idx] = img[idx];
    }
  }
}


static void od_paint_block(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by, level);
    od_paint_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by, level);
    od_paint_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by + 1, level);
    od_paint_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    ln = 2 + bs;
    n = 1 << ln;
    interp_block(&paint[stride*n*by + n*bx], &paint[stride*n*by + n*bx],
     n, stride, mode[(by*mstride + bx) << ln >> 2]);
  }
}

static void od_paint_switch(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *img1, const unsigned char *img2, const unsigned char *ref,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_switch(adapt, enc, img1, img2, ref, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by, level);
    od_paint_switch(adapt, enc, img1, img2, ref, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by, level);
    od_paint_switch(adapt, enc, img1, img2, ref, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by + 1, level);
    od_paint_switch(adapt, enc, img1, img2, ref, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int i;
    int j;
    int dist1;
    int dist2;
    ln = 2 + bs;
    n = 1 << ln;
    dist1 = dist2 = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        int x1;
        int x2;
        int r;
        r = ref[stride*(n*by + i) + n*bx + j];
        x1 = img1[stride*(n*by + i) + n*bx + j];
        x2 = img2[stride*(n*by + i) + n*bx + j];
        dist1 += (x1-r)*(x1-r);
        dist2 += (x2-r)*(x2-r);
      }
    }
    if (dist2 < dist1) {
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          img1[stride*(n*by + i) + n*bx + j] = img2[stride*(n*by + i) + n*bx + j];
        }
      }
    }
#if 0
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) img1[stride*(n*by + i) + n*bx + j] = 255;
#endif
#if 0
    for (i = 0; i < n; i++) img1[stride*(n*by + i) + n*bx + n - 1] = 0;
    for (i = 0; i < n; i++) img1[stride*(n*by + i) + n*bx + n - 2] = 0;
    for (i = 0; i < n; i++) img1[stride*(n*by + n - 1) + n*bx + i] = 0;
    for (i = 0; i < n; i++) img1[stride*(n*by + n - 2) + n*bx + i] = 0;
#endif
#if 0
    for (i = 0; i < n; i++) img1[stride*(n*by + i) + n*bx + n - 2] = img1[stride*(n*by + i) + n*bx + n - 1];
    for (i = 0; i < n; i++) img1[stride*(n*by + n - 2) + n*bx + i] = img1[stride*(n*by + n - 1) + n*bx + i];
#endif
#if 0
    {
      int m;
      int n2;
      m = mode[(by*mstride + bx) << ln >> 2];
      n2 = n/2;
      if (m == 4*n) {
        img1[stride*(n*by + n2) + n*bx + n2] = 0;
        img1[stride*(n*by + n2) + n*bx + n2 - 1] = 0;
        img1[stride*(n*by + n2 - 1) + n*bx + n2] = 0;
        img1[stride*(n*by + n2 - 1) + n*bx + n2 - 1] = 0;
      } else if (m <= 2*n) {
        double r;
        r = (double)m/n - 1;
        for (i = 0; i < n; i++) {
          int y;
          y = floor(.5 + n2 + r*(i-n2));
          img1[stride*(n*by + y) + n*bx + i] = 0;
          img1[stride*(n*by + y - 1) + n*bx + i] = 0;
          img1[stride*(n*by + y + 1) + n*bx + i] = 0;
        }
      }
      else {
        double r;
        r = 3 - (double)m/n;
        for (i = 0; i < n; i++) {
          int x;
          x = floor(.5 + n2 + r*(i-n2));
          img1[stride*(n*by + i) + n*bx + x] = 0;
          img1[stride*(n*by + i) + n*bx + x - 1] = 0;
          img1[stride*(n*by + i) + n*bx + x + 1] = 0;
        }
      }
    }
#endif
  }
}

static int var1[1<<24];
static int var2[1<<24];
static int var_count[1<<24];
static unsigned char mask[1<<24];
static unsigned char paint_buf[1<<24];
int *orig_edge_sum = var1+4096;
int *var2_edge_sum = var2+4096;
int *edge_corr = var_count+4096;
unsigned char *paint_mask=mask+4096;
unsigned char *paint_out=paint_buf+4096;

ogg_uint16_t gain_cdf[] = {128, 256, 384, 512, 640, 768, 896, 1024, 1152};

void od_paint_dering(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res) {
  int i, j;
  for(i = 0; i < 32*w; i++) paint[-stride + i] = paint[i];
  for(i = 0; i < 32*h; i++) paint[stride*i - 1] = paint[stride*i];
  paint[-stride - 1] = paint[0];
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_paint_mode_select(paint, paint, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_count, res, j, i, 3);
      od_intra_paint_analysis(paint, paint, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_paint_var_analysis(img, paint, stride, dec8, bstride, mode,
       mstride, orig_edge_sum, var2_edge_sum, edge_corr, res, j, i, 3);
    }
  }
#if 0
  for(i = 0; i < h*8; i++) {
    for(j = 0; j < w*8; j++) {
      printf("%d ", mode[i*mstride+j]);
      /*printf("%d ", mode[i*mstride+j]<<3>>dec8[(i>>1)*bstride + (j>>1)]);*/
    }
    printf("\n");
  }
#endif

  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_intra_paint_compute_edges(adapt, enc, paint_out, paint, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_count, res, j, i, 3);
      od_paint_compute_edge_mask(adapt, enc, paint_mask, paint_mask, stride, dec8, bstride, mode,
        mstride, edge_sum, var2_edge_sum, orig_edge_sum, edge_corr, edge_count, q, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
#if 1
      od_paint_block(adapt, enc, paint_out, paint_out, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
#endif
      od_paint_block(adapt, enc, paint_mask, paint_mask, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
    }
  }
#if 1
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      int idx;
      int gi;
      int k;
      int m;
      int dist;
      int best_dist;
      int best_gain;
      best_gain = 0;
      dist = 0;
      for (k = 0; k < 32; k++) {
        for (m = 0; m < 32; m++) {
          idx = (32*i+k)*stride + 32*j + m;
          dist += (img[idx] - paint[idx])*(int)(img[idx] - paint[idx]);
        }
      }
      best_dist = dist;
      for (gi = 1; gi <= 4; gi++) {
        dist = 0;
        for (k = 0; k < 32; k++) {
          for (m = 0; m < 32; m++) {
            int x;
            int y;
            idx = (32*i+k)*stride + 32*j + m;
            x = img[idx];
            y = OD_CLAMPI(0, paint[idx] + ((OD_MINI((int)paint_mask[idx]<<gi>>1, 255)*(paint_out[idx] - paint[idx]) + 128) >> 8), 255);
            dist += (x - y)*(int)(x - y);
          }
        }
        if (dist < best_dist) {
          best_dist = dist;
          best_gain = gi;
          for (k = 0; k < 32; k++) {
            for (m = 0; m < 32; m++) {
              int y;
              idx = (32*i+k)*stride + 32*j + m;
              y = OD_CLAMPI(0, paint[idx] + ((OD_MINI((int)paint_mask[idx]<<gi>>1, 255)*(paint_out[idx] - paint[idx]) + 128) >> 8), 255);
              paint[idx] = y;
            }
          }
        }
      }
      od_encode_cdf_adapt(enc, best_gain, gain_cdf, 9, 128);
      /*printf("%d ", best_gain);*/
    }
    /*printf("\n");*/
  }
#else
  for(i = 0; i < 32*h; i++) {
    for(j = 0; j < 32*w; j++) {
      int idx;
      idx = i*stride + j;
#if 1
      paint_mask[idx] = OD_CLAMPI(0, paint[idx] + (((int)paint_mask[idx]*(paint_out[idx] - paint[idx]) + 128) >> 8), 255);
#else
      paint[idx] = paint_out[idx];
#endif
    }
  }
#if 1
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_paint_switch(adapt, enc, paint, paint_mask, img, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_count, q, res, j, i, 3);
    }
  }
#endif
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
           stride, 2);
          cost8 += dist;
        }
      }
      cost8 = cost8 + 16000;
      mode_select(&img[16*stride*i + 16*j], &dist, 16, stride, 2);
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
  mode_select(img, &cost32, 32, stride, 2);
  if (cost32 < cost16+32000) {
    bsize[0][0] = bsize[0][1] = bsize[1][0] = bsize[1][1] = 3;
  }
}
