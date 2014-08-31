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

int ex_dc = 16384;
int ex_k[4] = {16384, 16384, 16384, 16384};
ogg_int32_t paint_adapt[OD_NSB_ADAPT_CTXS] = {384, 256, 104, 128};
int dc_prob[8] = {128, 128, 128, 128, 128, 128, 128, 128};
#define DC_PROB_SPEED (4)

int dir_prob[9] = {128, 128, 128, 128, 128, 128, 128, 128, 128};
#define DIR_PROB_SPEED (4)

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
static void compute_edge_variance(const unsigned char *img, int stride,
 int *edge_accum1, int *edge_accum2, int *edge_count, int edge_stride, int n,
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
        edge_accum1[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*w[k];
        edge_accum2[pi[k]*edge_stride+pj[k]] += (int)img[i*stride+j]*
         img[i*stride+j]*w[k];
        edge_count[pi[k]*edge_stride+pj[k]] += w[k];
      }
    }
  }
}


static void encode_edge_coeffs(od_adapt_ctx *adapt, od_ec_enc *enc, od_coeff *x, int n, int ctx) {
  int i;
  int k;
  ogg_int32_t adapt_curr[OD_NSB_ADAPT_CTXS];
  int speed;
  speed = 4;
  k = 0;
  for (i = 0; i < n; i++) k += abs(x[i]);
  generic_encode(enc, &adapt->paint_edge_k_model, k, -1, &ex_k[ctx], 6);
  laplace_encode_vector(enc, x, n, k, adapt_curr, paint_adapt);

  if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
    paint_adapt[OD_ADAPT_K_Q8] += 256*adapt_curr[OD_ADAPT_K_Q8] -
      paint_adapt[OD_ADAPT_K_Q8]>>speed;
    paint_adapt[OD_ADAPT_SUM_EX_Q8] += adapt_curr[OD_ADAPT_SUM_EX_Q8] -
      paint_adapt[OD_ADAPT_SUM_EX_Q8]>>speed;
  }
  if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
    paint_adapt[OD_ADAPT_COUNT_Q8] += adapt_curr[OD_ADAPT_COUNT_Q8]-
      paint_adapt[OD_ADAPT_COUNT_Q8]>>speed;
    paint_adapt[OD_ADAPT_COUNT_EX_Q8] += adapt_curr[OD_ADAPT_COUNT_EX_Q8]-
      paint_adapt[OD_ADAPT_COUNT_EX_Q8]>>speed;
  }

}

/* Quantize the bottom edge using prediction. */
static void quantize_bottom_edge(od_adapt_ctx *adapt,od_ec_enc *enc, unsigned char *edge_accum, int n, int stride, int q,
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
  predict_bottom_edge(p, edge_accum - stride - 1, n, stride, m, has_right);
#if !QUANTIZE
  /*printf("%d ", m);*/
  for (i = 0; i < n; i++) printf("%d ", edge_accum[(n - 1)*stride + i] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[(n - 1)*stride + i] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) x[i] = (int)(floor(.5+x[i]/q));
  encode_edge_coeffs(adapt, enc, x, n, has_right || (m >= n && m <= 3*n));
  for (i = 0; i < n; i++) x[i] = q*x[i];
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[(n - 1)*stride + i] = OD_MAXI(0, OD_MINI(255, r[i]));
}

/* Quantize the right edge using prediction. */
static void quantize_right_edge(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *edge_accum, int n, int stride, int q,
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
  predict_right_edge(p, edge_accum - stride - 1, n, stride, m, has_bottom);
#if !QUANTIZE
  for (i = 0; i < n; i++) printf("%d ", edge_accum[i*stride + n - 1] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[i*stride + n - 1] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) x[i] = (int)(floor(.5+x[i]/q));
  encode_edge_coeffs(adapt, enc, x, n, (has_bottom || (m >= n && m <= 3*n)) + 2);
  for (i = 0; i < n; i++) x[i] = q*x[i];
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[i*stride + n - 1] = OD_MAXI(0, OD_MINI(255, r[i]));
}

/* Quantize both the right and bottom edge, changing the order to maximize
   the number of pixels we can predict. */
void quantize_edge(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *edge_accum, int n, int stride, int q,
 int m, int dc_quant) {
  /*printf("%d ", m);*/
  if (dc_quant) {
    int pred;
    int res;
    int qdc;
    int i;
    pred = floor(.5 + .73*(edge_accum[n - 1] + edge_accum[(n - 1)*stride]) - .46 *edge_accum[-stride - 1]);
    res = edge_accum[(n - 1)*stride + n - 1]-pred;
    qdc = OD_MAXI(1, q/8);
    res = (int)floor(.5+res/qdc);
    generic_encode(enc, &adapt->paint_dc_model, abs(res), -1, &ex_dc, 6);
    if (res != 0) od_ec_enc_bits(enc, res > 0, 1);
    /*printf("DC %d\n", res);*/
    res = res*qdc;
    edge_accum[(n - 1)*stride + n - 1] = res + pred;
    for (i = 0; i < n - 1; i++) {
      edge_accum[(n - 1)*stride + i] = edge_accum[(n - 1)*stride - 1]
       + i*(edge_accum[(n - 1)*stride + n - 1]-edge_accum[(n - 1)*stride - 1])/n;
      edge_accum[i*stride + n - 1] = edge_accum[-stride + n - 1]
       + i*(edge_accum[(n - 1)*stride + n - 1]-edge_accum[-stride + n - 1])/n;
    }
  }
  else if (m >= 0 && m < 2*n) {
    quantize_right_edge(adapt, enc, edge_accum, n, stride, q, m, 0);
    quantize_bottom_edge(adapt, enc, edge_accum, n, stride, q, m, 1);
  }
  else {
    quantize_bottom_edge(adapt, enc, edge_accum, n, stride, q, m, 0);
    quantize_right_edge(adapt, enc, edge_accum, n, stride, q, m, 1);
  }
  /*printf("\n");*/
}

void od_intra_paint_analysis(const unsigned char *img,
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
    od_intra_paint_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by, level);
    od_intra_paint_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx + 1, by, level);
    od_intra_paint_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, res, bx, by + 1, level);
    od_intra_paint_analysis(img, stride, dec8, bstride,
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
    curr_mode = mode_select(&img[stride*n*by + n*bx], NULL, n, stride, res);
    compute_edges(&img[stride*n*by + n*bx], stride,
     &edge_sum[stride*n*by + n*bx], &edge_count[stride*n*by + n*bx], stride,
     n, curr_mode);
    mode[(by<<bs)*mstride + (bx<<bs)] = curr_mode;
    for (k=0;k<1<<bs;k++) {
      for (m=0;m<1<<bs;m++) {
        mode[((by << bs) + k)*mstride + (bx << bs) + m] = curr_mode;
      }
    }
  }
}

void od_paint_var_analysis(const unsigned char *img,
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
    od_paint_var_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by, level);
    od_paint_var_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by, level);
    od_paint_var_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by + 1, level);
    od_paint_var_analysis(img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int curr_mode;
    ln = 2 + bs;
    n = 1 << ln;
    curr_mode = mode[(by<<bs)*mstride + (bx<<bs)];
    compute_edge_variance(&img[stride*n*by + n*bx], stride,
     &edge_sum1[stride*n*by + n*bx], &edge_sum2[stride*n*by + n*bx],
     &edge_count[stride*n*by + n*bx], stride, n, curr_mode);
  }
}


void od_intra_paint_mode_encode(od_ec_enc *enc, const unsigned char *mode, int bx, int by,
 int ln, int mstride, const unsigned char *dec8, int bstride, int res) {
  int i;
  int m;
  int idx;
  int dc_ctx;
  int nb;
  int bs;
  int prob_list[10];
  int dir_list[10];
  int ctx_list[10];
  int cnt;
  int in_list;
  ogg_uint16_t cdf[16];
  bs = ln - 2;
  nb = 4 << ln >> res;
  idx = (by*mstride + bx) << bs;
  m = mode[idx] >> res;
  OD_ASSERT(m <= nb);
  cnt = od_intra_paint_mode_cdf(cdf, dir_list, prob_list, ctx_list, &dc_ctx,
   mode, bx, by, ln, mstride, dec8, bstride, res);

  if (m == nb) {
    od_ec_encode_cdf_unscaled(enc, 1, cdf, cnt + 2);
    dc_prob[dc_ctx] += (256 - dc_prob[dc_ctx]) >> DC_PROB_SPEED;
  } else {
    dc_prob[dc_ctx] -= dc_prob[dc_ctx] >> DC_PROB_SPEED;

    in_list = 0;
    for (i = 0; i < cnt; i++) {
      if (dir_list[i] == m) {
        od_ec_encode_cdf_unscaled(enc, i + 2, cdf, cnt + 2);
        in_list = 1;
        dir_prob[ctx_list[i]] += (256 - dir_prob[ctx_list[i]]) >> DIR_PROB_SPEED;
      }
      else {
        dir_prob[ctx_list[i]] -= dir_prob[ctx_list[i]] >> DIR_PROB_SPEED;
      }
    }
    if (!in_list) {
      int other_id;
      od_ec_encode_cdf_unscaled(enc, 0, cdf, cnt + 2);
      other_id = m;
      /* Don't count the directions we had in the list. */
      for (i = 0; i < cnt; i++) if (dir_list[i] < m) other_id--;
      od_ec_enc_uint(enc, other_id, nb - cnt);
      dir_prob[0] += (256 - dir_prob[0]) >> DIR_PROB_SPEED;
    } else {
      dir_prob[0] -= dir_prob[0] >> DIR_PROB_SPEED;
    }
  }
}

static void quantize_initial_edge(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, int ln,
 int stride, int q) {
  int pred;
  int i;
  int n;
  od_coeff x[MAXN];
  od_coeff r[MAXN];
  n = 1 << ln;
  pred = paint[-stride];
  for (i = 0; i < n; i++) r[i] = paint[i*stride] - pred;
  my_fdct_table[ln - 2](x, r, 1);
  for (i = 0; i < n; i++) x[i] = (int)(floor(.5+x[i]/q));
  encode_edge_coeffs(adapt, enc, x, n, 0);
  for (i = 0; i < n; i++) x[i] = q*x[i];
  my_idct_table[ln - 2](r, 1, x);
  for (i = 0; i < n; i++) paint[i*stride] = r[i] + pred;
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
#if 1
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
#define VAR do {paint[idx] = OD_MINI(255, 600/(int)(1+sqrt((double)edge_sum2[idx]/edge_count[idx] \
  - SQUARE((double)edge_sum1[idx]/edge_count[idx]))));} while(0)

void od_paint_compute_edge_mask(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum1, int *edge_sum2, int *edge_count,
 int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx, by + 1, level);
    od_paint_compute_edge_mask(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum1, edge_sum2, edge_count, res, bx + 1, by + 1, level);
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
      if (edge_count[idx] > 0) VAR;
      else paint[idx] = img[idx];
    }
    /* Compute left edge (left column only). */
    if (bx == 0) {
      for (k = 0; k < n; k++) {
        idx = stride*(n*by + k) - 1;
        if (edge_count[idx] > 0) VAR;
        else paint[idx] = img[idx];
      }
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      for (k = 0; k < n; k++) {
        idx = -stride + n*bx + k;
        if (edge_count[idx] > 0) VAR;
        else paint[idx] = img[idx];
      }
    }
    /* Compute right edge stats. */
    for (k = 0; k < n - 1; k++) {
      idx = stride*(n*by + k) + n*(bx + 1) - 1;
      if (edge_count[idx] > 0) VAR;
      else paint[idx] = img[idx];
    }
    /* Compute bottom edge stats. */
    for (k = 0; k < n; k++) {
      idx = stride*(n*(by + 1) - 1) + n*bx + k;
      if (edge_count[idx] > 0) VAR;
      else paint[idx] = img[idx];
    }
  }
}

void od_intra_paint_quant_block(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
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
    od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by, level);
    od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by, level);
    od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by + 1, level);
    od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int dc_quant;
    ln = 2 + bs;
    n = 1 << ln;
    if (bx == 0 && by == 0) {
      /* For now, we just encode the first pixel as is. */
      od_ec_enc_bits(enc, paint[-stride - 1], 8);
    }
    /* Compute left edge (left column only). */
    if (bx == 0) {
      quantize_initial_edge(adapt, enc, &paint[stride*n*by - 1], ln, stride, q);
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      quantize_initial_edge(adapt, enc, &paint[-stride + n*bx], ln, 1, q);
    }
    od_intra_paint_mode_encode(enc, mode, bx, by, ln, mstride, dec8, bstride, res);
    /* Only use special DC quantization when the two adjacent blocks are
       also DC. In the future, we could also treat each edge separately. */
    dc_quant = mode[(by*mstride + bx) << ln >> 2]==4*n
     && mode[(by*mstride + bx + 1) << ln >> 2] == 4*n
     && mode[((by + 1)*mstride + bx) << ln >> 2] == 4*n;
    quantize_edge(adapt, enc, &paint[stride*n*by + n*bx], n, stride, q,
     mode[(by*mstride + bx) << ln >> 2], dc_quant);
    interp_block(&paint[stride*n*by + n*bx], &paint[stride*n*by + n*bx],
     n, stride, mode[(by*mstride + bx) << ln >> 2]);
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

static int od_intra_paint_dist(const unsigned char *img, int stride,
 unsigned char *mode, int hb, int vb, int bx, int by, int n) {
  int edge_accum[(3*MAXN + 1)*(3*MAXN + 1)];
  int edge_count[(3*MAXN + 1)*(3*MAXN + 1)];
  unsigned char block[(MAXN + 1)*(MAXN + 1)];
  int i, j;
  int dist;
  char buf[1024];
  int off=0;
  /*fprintf(stdout,"bx=%i by=%i hb=%i vb=%i n=%i mode=%i\n",bx,by,hb,vb,n,mode[hb*(by - 1) + bx]);*/
  for (i = 0; i < (3*MAXN + 1)*(3*MAXN + 1); i++) {
    edge_accum[i] = edge_count[i] = 0;
  }
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      int tx, ty;
      tx = bx + i - 1;
      ty = by + j - 1;
      if (ty >= 0 && ty < vb && tx >= 0 && tx < hb) {
        sprintf(&buf[off],"(%2i,%2i)=%2i ",tx,ty,mode[hb*ty + tx]);
        off+=11;
        compute_edges(&img[stride*ty*n + tx*n], stride,
         &edge_accum[(3*MAXN + 1)*(n*j + 1) + n*i + 1], &edge_count[(3*MAXN + 1)*(n*j + 1) + n*i + 1],
         3*MAXN + 1, n, mode[hb*ty + tx]);
      }
    }
  }
  /*for (j = 0; j < 3*n + 1; j++) {
    for (i = 0; i < 3*n + 1; i++) {
      fprintf(stdout,"%d ",edge_accum[(3*MAXN + 1)*j + i]);
    }
    fprintf(stdout,"\n");
  }*/
  for (i = 0; i < (MAXN + 1)*(MAXN + 1); i++) {
    block[i] = 0;
  }
  for (j = 0; j <= n; j++) {
    for (i = 0; i <= n; i++) {
      if (edge_count[(3*MAXN + 1)*(n + j) + n + i] > 0) {
        block[(MAXN + 1)*j + i] =
         edge_accum[(3*MAXN + 1)*(n + j) + n + i]/edge_count[(3*MAXN + 1)*(n + j) + n + i];
      }
    }
  }
  /*fprintf(stdout,"pre interp\n");
  for (j = 0; j < n + 1; j++) {
    for (i = 0; i < n + 1; i++) {
      fprintf(stdout,"%d ",block[(MAXN + 1)*j + i]);
    }
    fprintf(stdout,"\n");
  }*/
  interp_block(&block[(MAXN + 1)*1 + 1], &block[(MAXN + 1)*1 + 1], n, MAXN + 1, mode[hb*by + bx]);
  /*fprintf(stdout,"post interp\n");
  for (j = 0; j < n + 1; j++) {
    for (i = 0; i < n + 1; i++) {
      fprintf(stdout,"%d ",block[(MAXN + 1)*j + i]);
    }
    fprintf(stdout,"\n");
  }*/
  dist = 0;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      int e;
      e = (int)block[(MAXN + 1)*(j + 1) + (i + 1)] - (int)img[stride*(by*n + j) + (bx*n + i)];
      dist += e*e;
    }
  }
  /*fprintf(stderr,"%s: %i\n",buf,dist);*/
  return dist;
}

static void od_intra_paint_find_modes(const unsigned char *img, int stride,
 int hb, int vb, unsigned char *mode, int n, int res, int max_steps) {
  int i, j, k, l;
  /* Initialize the modes with a greedy search */
  for (j = 0; j < vb; j++) {
    for (i = 0; i < hb; i++) {
      mode[hb*j + i] = mode_select(&img[stride*n*j + n*i], NULL, n, stride, res);
    }
  }

#define PRINT_PAINT_DEBUG (1)

#if 1
  {
    int max_blocks;
    int max_modes;
    int *nmodes;
    unsigned char *modes;
    unsigned char *last;
    unsigned char *prev;
    int *costs;
    int step;
    int step_pre, step_post;
    step = 0;
    max_blocks = OD_MAXI(hb, vb);
    max_modes = (4*n >> res) + 1;
    nmodes = (int *)malloc(max_blocks*sizeof(*nmodes));
    modes = (unsigned char *)malloc(max_modes*max_blocks*sizeof(*modes));
    last = (unsigned char *)malloc(max_blocks*sizeof(*last));
    prev = (unsigned char *)malloc(max_blocks*max_modes*sizeof(*prev));
    costs = (int *)malloc(max_blocks*max_modes*sizeof(*costs));
    do {
      int nh, nv;
#if PRINT_PAINT_DEBUG
      fprintf(stdout,"step %i\n",step);
      for (j = 0; j < vb; j++) {
        for (i = 0; i < hb; i++) {
          fprintf(stdout,"%2d ",mode[hb*j + i]);
        }
        fprintf(stdout,"\n");
      }
      fflush(stdout);
#endif
      step_pre = 0;
      for (j = 0; j < vb; j++) {
        for (i = 0; i < hb; i++) {
          step_pre += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n);
        }
      }
      /* Compute optimal mode assignment per row, holding neighbors constant */
      nv = 0;
      for (j = 0; j < vb; j++) {
        int pre_cost, post_cost;
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"row %i\n",j);
        fflush(stdout);
#endif
        /* Keep a copy of the current row */
        for (i = 0; i < hb; i++) {
          last[i] = mode[hb*j + i];
        }
        /* Create a table of modes to consider per block */
        for (i = 0; i < hb; i++) {
          unsigned char *b;
          b = &modes[max_modes*i];
          for (k = 0; k <= 4*n; k += 1 << res) {
            *b=k;
            b++;
          }
          nmodes[i] = (4*n >> res) + 1;
        }
        /* Compute the row cost pre optimization. */
        pre_cost = 0;
        for (i = 0; i < hb; i++) {
          pre_cost += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
           (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j - 1, n) : 0) +
           (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j + 1, n) : 0);
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"%i: pre_cost=%i\n",j,pre_cost);
#endif
        /* Initialize the first stage costs. */
        for (k = 0; k < nmodes[0]; k++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"(%i,%i) mode=%i\n",0,j,mode[hb*j+0]);
#endif
          mode[hb*j + 0] = modes[max_modes*0 + k];
          /* Compute the cost as MSE distortion. */
          prev[max_modes*0 + k] = -1;
          costs[max_modes*0 + k] = 0;
          if (0 + 1 == hb) {
            costs[max_modes*0 + k] = od_intra_paint_dist(img, stride, mode, hb, vb, 0, j, n) +
             (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, 0, j - 1, n) : 0) +
             (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, 0, j + 1, n) : 0);
          }
#if PRINT_PAINT_DEBUG
          if (j < 1) {
            fprintf(stdout,"block (0, 0), mode=%i cost=%i\n",
            modes[max_modes*0 + k],
            costs[max_modes*0 + k]);
            fflush(stdout);
          }
#endif
        }
        for (i = 1; i < hb; i++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"col %i\n",i);
          for (k = 0; k < nmodes[i - 1]; k++) {
            fprintf(stdout,"  %i: cost=%i prev=%i\n",
             modes[max_modes*(i - 1) + k],
             costs[max_modes*(i - 1) + k],
             i > 1 ? modes[max_modes*(i - 2) + prev[max_modes*(i - 1) + k]] : -1);
          }
#endif
          /* Compute the optimal mode to select for the previous block, given
             each mode for the current block. */
          for (k = 0; k < nmodes[i]; k++) {
            int cost;
            mode[hb*j + i] = modes[max_modes*i + k];
            mode[hb*j + (i - 1)] = modes[max_modes*(i - 1) + 0];
            if (i > 1) {
              mode[hb*j + (i - 2)] = modes[max_modes*(i - 2) + prev[max_modes*(i - 1) + 0]];
            }
            cost = costs[max_modes*(i - 1) + 0] +
             od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j, n) +
             (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j - 1, n) : 0) +
             (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j + 1, n) : 0);
#if PRINT_PAINT_DEBUG
            fprintf(stdout,"(%2i, %2i): modes=(%i,%i) cost=%i\n",
              modes[max_modes*(i - 1) + 0],
              modes[max_modes*i + k],
              mode[hb*j + i - 1],
              mode[hb*j + i],
              cost);
#endif
            prev[max_modes*i + k] = 0;
            costs[max_modes*i + k] = cost;
            for (l = 1; l < nmodes[i - 1]; l++) {
              mode[hb*j + (i - 1)] = modes[max_modes*(i - 1) + l];
              if (i > 1) {
                mode[hb*j + (i - 2)] = modes[max_modes*(i - 2) + prev[max_modes*(i - 1) + l]];
              }
              cost = costs[max_modes*(i - 1) + l] +
               od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j, n) +
               (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j - 1, n) : 0) +
               (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j + 1, n) : 0);
#if PRINT_PAINT_DEBUG
              fprintf(stdout,"(%2i, %2i): modes=(%i,%i) cost=%i\n",
               modes[max_modes*(i - 1) + l],
               modes[max_modes*i + k],
               mode[hb*j + i - 1],
               mode[hb*j + i],
               cost);
#endif
              if (cost < costs[max_modes*i + k]) {
                prev[max_modes*i + k] = l;
                costs[max_modes*i + k] = cost;
              }
            }
            if (i + 1 == hb) {
              mode[hb*j + (i - 1)] = modes[max_modes*(i - 1) + prev[max_modes*i + k]];
              if (i > 1) {
                mode[hb*j + (i - 2)] = modes[max_modes*(i - 2) + prev[max_modes*(i - 1) + prev[max_modes*i + k]]];
              }
              costs[max_modes*i + k] +=
               od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
               (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j - 1, n) : 0) +
               (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j + 1, n) : 0);
            }
          }
        }
#if PRINT_PAINT_DEBUG
        for (i = 0; i < hb; i++) {
          for (k = 0; k < nmodes[i]; k++) {
            fprintf(stdout,"%2d ",i==0?-1:prev[max_modes*i + k]);
          }
          fprintf(stdout,"\n");
        }
        for (k = 0; k < nmodes[hb - 1]; k++) {
          fprintf(stdout,"  %i: cost=%i prev=%i\n",
           modes[max_modes*(hb - 1) + k],
           costs[max_modes*(hb - 1) + k],
           i > 1 ? modes[max_modes*(hb - 2) + prev[max_modes*(hb - 1) + k]] : -1);
        }
#endif
        /* Find the mode with the lowest cost for the last block. */
        l = 0;
        for (k = 1; k < nmodes[hb - 1]; k++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"cost block=%i mode=%i is %i\n",hb-1,k,costs[max_modes*(hb - 1) + k]);
#endif
          if (costs[max_modes*(hb - 1) + k] < costs[max_modes*(hb - 1) + l]) {
            l = k;
          }
        }
        /* Walk backwards setting the mode for each block in the row. */
        /*if (last[hb - 1] != modes[max_modes*(hb - 1) + l]) {
          mode[hb*j + (hb - 1)] = modes[max_modes*(hb - 1) + l];
          printf("first changing block %i\n",hb -1);
          nv++;
        }*/
        for (i = hb; i-- > 0; ) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"i=%i l=%i last=%i prev=%i\n",i,l,last[i],prev[max_modes*i + l]);
#endif
          mode[hb*j + i] = modes[max_modes*i + l];
          if (last[i] != modes[max_modes*i + l]) {
#if PRINT_PAINT_DEBUG
            printf("changing block %i from %i to %i\n",i,last[i],modes[max_modes*i + l]);
#endif
            nv++;
          }
          l = prev[max_modes*i + l];
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"row %i finished (%i):\n",j,nv);
        for (k = 0; k < 3; k++) {
          if (j + k - 1 >= 0 && j + k - 1 < vb) {
            fprintf(stdout,"%2i: ",j + k - 1);
            for (i = 0; i < hb; i++) {
              fprintf(stdout,"%2d ",mode[hb*(j + k - 1) + i]);
            }
            fprintf(stdout,"\n");
          }
        }
#endif
        post_cost = 0;
        for (i = 0; i < hb; i++) {
          post_cost += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
           (j - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j - 1, n) : 0) +
           (j + 1 < vb ? od_intra_paint_dist(img, stride, mode, hb, vb, i, j + 1, n) : 0);
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"%i: post_cost %i\n",j,post_cost);
#endif
      }
      /* Compute optimal mode assignment per column, holding neighbors constant */
      nh = 0;
      for (i = 0; i < hb; i++) {
        int pre_cost, post_cost;
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"col %i\n",j);
        fflush(stdout);
#endif
        /* Keep a copy of the current column */
        for (j = 0; j < vb; j++) {
          last[j] = mode[hb*j + i];
        }
        /* Create a table of modes to consider per block */
        for (j = 0; j < vb; j++) {
          unsigned char *b;
          b = &modes[max_modes*j];
          for (k = 0; k <= 4*n; k += 1 << res) {
            *b=k;
            b++;
          }
          nmodes[j] = (4*n >> res) + 1;
        }
        /* Compute the column cost pre optimization. */
        pre_cost = 0;
        for (j = 0; j < vb; j++) {
          pre_cost += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
           (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j, n) : 0) +
           (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, j, n) : 0);
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"%i: pre_cost=%i\n",i,pre_cost);
#endif
        /* Initialize the first stage costs. */
        for (k = 0; k < nmodes[0]; k++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"(%i,%i) mode=%i\n",0, i, mode[hb*0 + i]);
#endif
          mode[hb*0 + i] = modes[max_modes*0 + k];
          /* Compute the cost as MSE distortion. */
          prev[max_modes*0 + k] = -1;
          costs[max_modes*0 + k] = 0;
          if (0 + 1 == vb) {
            costs[max_modes*0 + k] = od_intra_paint_dist(img, stride, mode, hb, vb, i, 0, n) +
             (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, 0, n) : 0) +
             (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, 0, n) : 0);
          }
#if PRINT_PAINT_DEBUG
          if (i < 1) {
            fprintf(stdout,"block (0, 0), mode=%i cost=%i\n",
            modes[max_modes*0 + k],
            costs[max_modes*0 + k]);
            fflush(stdout);
          }
#endif
        }
        for (j = 1; j < vb; j++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"row %i\n",j);
          for (k = 0; k < nmodes[j - 1]; k++) {
            fprintf(stdout,"  %i: cost=%i prev=%i\n",
             modes[max_modes*(j - 1) + k],
             costs[max_modes*(j - 1) + k],
             j > 1 ? modes[max_modes*(j - 2) + prev[max_modes*(j - 1) + k]] : -1);
          }
#endif
          /* Compute the optimal mode to select for the previous block, given
             each mode for the current block. */
          for (k = 0; k < nmodes[j]; k++) {
            int cost;
            mode[hb*j + i] = modes[max_modes*j + k];
            mode[hb*(j - 1) + i] = modes[max_modes*(j - 1) + 0];
            if (j > 1) {
              mode[hb*(j - 2) + i] = modes[max_modes*(j - 2) + prev[max_modes*(j - 1) + 0]];
            }
            cost = costs[max_modes*(j - 1) + 0] +
             od_intra_paint_dist(img, stride, mode, hb, vb, i, j - 1, n) +
             (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j - 1, n) : 0) +
             (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, j - 1, n) : 0);
#if PRINT_PAINT_DEBUG
            fprintf(stdout,"(%2i, %2i): modes=(%i,%i) cost=%i\n",
              modes[max_modes*(j - 1) + 0],
              modes[max_modes*j + k],
              mode[hb*(j - 1) + i],
              mode[hb*j + i],
              cost);
#endif
            prev[max_modes*j + k] = 0;
            costs[max_modes*j + k] = cost;
            for (l = 1; l < nmodes[j - 1]; l++) {
              mode[hb*(j - 1) + i] = modes[max_modes*(j - 1) + l];
              if (j > 1) {
                mode[hb*(j - 2) + i] = modes[max_modes*(j - 2) + prev[max_modes*(j - 1) + l]];
              }
              cost = costs[max_modes*(j - 1) + l] +
               od_intra_paint_dist(img, stride, mode, hb, vb, i, j - 1, n) +
               (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j - 1, n) : 0) +
               (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, j - 1, n) : 0);
#if PRINT_PAINT_DEBUG
              fprintf(stdout,"(%2i, %2i): modes=(%i,%i) cost=%i\n",
               modes[max_modes*(j - 1) + l],
               modes[max_modes*j + k],
               mode[hb*(j - 1) + i],
               mode[hb*j + i],
               cost);
#endif
              if (cost < costs[max_modes*j + k]) {
                prev[max_modes*j + k] = l;
                costs[max_modes*j + k] = cost;
              }
            }
            if (j + 1 == vb) {
              mode[hb*(j - 1) + i] = modes[max_modes*(j - 1) + prev[max_modes*j + k]];
              if (j > 1) {
                mode[hb*(j - 2) + i] = modes[max_modes*(j - 2) + prev[max_modes*(j - 1) + prev[max_modes*j + k]]];
              }
              costs[max_modes*j + k] +=
               od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
               (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j, n) : 0) +
               (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, j, n) : 0);
            }
          }
        }
#if PRINT_PAINT_DEBUG
        for (j = 0; j < vb; j++) {
          for (k = 0; k < nmodes[j]; k++) {
            fprintf(stdout,"%2d ",j==0?-1:prev[max_modes*j + k]);
          }
          fprintf(stdout,"\n");
        }
        for (k = 0; k < nmodes[vb - 1]; k++) {
          fprintf(stdout,"  %i: cost=%i prev=%i\n",
           modes[max_modes*(vb - 1) + k],
           costs[max_modes*(vb - 1) + k],
           i > 1 ? modes[max_modes*(vb - 2) + prev[max_modes*(vb - 1) + k]] : -1);
        }
#endif
        /* Find the mode with the lowest cost for the last block. */
        l = 0;
        for (k = 1; k < nmodes[vb - 1]; k++) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"cost block=%i mode=%i is %i\n",vb-1,k,costs[max_modes*(vb - 1) + k]);
#endif
          if (costs[max_modes*(vb - 1) + k] < costs[max_modes*(vb - 1) + l]) {
            l = k;
          }
        }
        /* Walk backwards setting the mode for each block in the row. */
        /*if (last[hb - 1] != modes[max_modes*(hb - 1) + l]) {
          mode[hb*j + (hb - 1)] = modes[max_modes*(hb - 1) + l];
          printf("first changing block %i\n",hb -1);
          nv++;
        }*/
        for (j = vb; j-- > 0; ) {
#if PRINT_PAINT_DEBUG
          fprintf(stdout,"j=%i l=%i last=%i prev=%i\n",j,l,last[j],prev[max_modes*j + l]);
#endif
          mode[hb*j + i] = modes[max_modes*j + l];
          if (last[j] != modes[max_modes*j + l]) {
#if PRINT_PAINT_DEBUG
            printf("changing block %i from %i to %i\n",j,last[j],modes[max_modes*j + l]);
#endif
            nh++;
          }
          l = prev[max_modes*j + l];
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"col %i finished (%i):\n",i,nh);
        for (k = 0; k < 3; k++) {
          if (i + k - 1 >= 0 && i + k - 1 < hb) {
            fprintf(stdout,"%2i: ",i + k - 1);
            for (j = 0; j < vb; j++) {
              fprintf(stdout,"%2d ",mode[hb*j + i + k - 1]);
            }
            fprintf(stdout,"\n");
          }
        }
#endif
        post_cost = 0;
        for (j = 0; j < vb; j++) {
          post_cost += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n) +
           (i - 1 >= 0 ? od_intra_paint_dist(img, stride, mode, hb, vb, i - 1, j, n) : 0) +
           (i + 1 < hb ? od_intra_paint_dist(img, stride, mode, hb, vb, i + 1, j, n) : 0);
        }
#if PRINT_PAINT_DEBUG
        fprintf(stdout,"%i: post_cost=%i\n",i,post_cost);
#endif
      }
      step_post = 0;
      for (j = 0; j < vb; j++) {
        for (i = 0; i < hb; i++) {
          step_post += od_intra_paint_dist(img, stride, mode, hb, vb, i, j, n);
        }
      }
      fprintf(stderr,"Step %i: nv=%i nh=%i pre=%i post=%i\n",step,nv,nh,step_pre,step_post);
      step++;
    }
    while (step_pre > step_post && step < max_steps);
    free(nmodes);
    free(modes);
    free(last);
    free(prev);
    free(costs);
  }
#endif
}

void od_intra_paint_encode2(od_adapt_ctx *adapt, od_ec_enc *enc,
 unsigned char *paint, const unsigned char *img, int stride, int nhsb, int nvsb,
 unsigned char *dec8, unsigned char *mode4, unsigned char *mode8,
 unsigned char *mode16, unsigned char *mode32, int *edge_sum, int *edge_count,
 int q, int res) {
  int i,j,k,l;
  int w4,h4,w8,h8,w16,h16;
  (void)mode8;
  (void)mode32;
  /* Compute the MSE optimal mode selection for 16x16 blocks using DP */
  w16 = nhsb << 1;
  h16 = nvsb << 1;
  od_intra_paint_find_modes(img, stride, w16, h16, mode16, 16, res, 1000);

  /* Set the paint decision to be all 16x16 blocks */
  w8 = nhsb << 2;
  h8 = nvsb << 2;
  for (j = 0; j < h8; j++) {
    for (i = 0; i < w8; i++) {
      dec8[w8*j + i] = 2;
    }
  }

  /* Set the mode decision at the 4x4 scale */
  w4 = nhsb << 3;
  h4 = nvsb << 3;
  for (j = 0; j < h16; j++) {
    for (i = 0; i < w16; i++) {
      for (l = 0; l < 4; l++) {
        for (k = 0; k < 4; k++) {
          mode4[w4*((j << 2) + l) + (i << 2) + k] = mode16[w16*j + i];
        }
      }
    }
  }

#if 1
  {
    for (j = 0; j < h4; j++) {
      for (i = 0; i < w4; i++) {
        printf("%d ",mode4[j*w4+i]);
      }
      printf("\n");
    }
  }
#endif

  for (j = 0; j < nvsb; j++) {
    for (i = 0; i < nhsb; i++) {
      od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, w8,
       mode4, w4, edge_sum, edge_count, res, i, j, 3);
    }
  }

  for (j = 0; j < nvsb; j++) {
    for (i = 0; i < nhsb; i++) {
      od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, w8,
       mode4, w4, edge_sum, edge_count, q, res, i, j, 3);
    }
  }
}

void od_intra_paint_encode(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res) {
  int i, j;
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_intra_paint_analysis(img, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, res, j, i, 3);
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
      od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_count, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
    }
  }
}

static int var1[1<<24];
static int var2[1<<24];
static int var_count[1<<24];
static unsigned char mask[1<<24];
static unsigned char paint_buf[1<<24];
int *var1_edge_sum = var1+4096;
int *var2_edge_sum = var2+4096;
int *var_edge_count = var_count+4096;
unsigned char *paint_mask=mask+4096;
unsigned char *paint_out=paint_buf+4096;

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
      od_intra_paint_analysis(img, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_paint_var_analysis(img, stride, dec8, bstride, mode,
       mstride, var1_edge_sum, var2_edge_sum, var_edge_count, res, j, i, 3);
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
        mstride, var1_edge_sum, var2_edge_sum, var_edge_count, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_paint_block(adapt, enc, paint_out, paint_out, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
      od_paint_block(adapt, enc, paint_mask, paint_mask, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
    }
  }
  for(i = 0; i < 32*h; i++) {
    for(j = 0; j < 32*w; j++) {
      int idx;
      idx = i*stride + j;
      paint[idx] = OD_CLAMPI(0, paint[idx] + ((int)paint_mask[idx]*(paint_out[idx] - paint[idx]) >> 8), 255);
    }
  }
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

#if 0

/* Eventually we shouldn't be using this function, but calling
   od_intra_paint_choose_block_size() and od_intra_paint_encode() directly. */
int compute_intra_paint(od_ec_enc *enc, unsigned char *img, int w, int h, int stride)
{
  int i,j;
  int h8,w8,h32,w32;
  unsigned char *dec8;
  unsigned char *mode;
  unsigned char *paint_out;
  int *edge_sum;
  int *edge_count;
  int bstride;
  int mstride;

  generic_model_init(&paint_dc_model);
  generic_model_init(&paint_edge_k_model);

  paint_out = (unsigned char*)calloc(stride*(h+32)*sizeof(*paint_out), 1);
  edge_sum = (int*)calloc((w+32)*(h+32), sizeof(*edge_sum));
  edge_count = (int*)calloc((w+32)*(h+32), sizeof(*edge_count));
  w32 = w>>5;
  h32 = h>>5;
  w8 = w32<<2;
  h8 = h32<<2;
  dec8 = (unsigned char*)malloc(w8*h8*sizeof(*dec8));
  mode = (unsigned char*)calloc((w8 + 1)*(h8 + 1)*sizeof(*mode)<<2, 1);
  bstride = w8;
  mstride = (w8 + 1)<<1;
  /* Replace decision with the one from process_block_size32() */
  for(i=0;i<h32;i++){
    for(j=0;j<w32;j++){
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
  od_intra_paint_encode(enc, paint_out, img, w32, h32, stride, dec8, bstride, mode,
   mstride, edge_sum, edge_count, 30, 1);
  for(i=0;i<32*h32;i++) {
    for(j=0;j<32*w32;j++) {
      img[i*stride + j] = paint_out[i*stride + j];
    }
  }
#if 0
  for(i=0;i<h8;i++){
    for(j=0;j<w8;j++){
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
  for (i=0;i<w32*32;i++)
    img[(h32*32-1)*stride+i]=0;
  for (i=0;i<h32*32;i++)
    img[i*stride+w32*32-1]=0;
#endif
  free(paint_out);
  free(edge_sum);
  free(edge_count);
  free(dec8);
  free(mode);
  return 0;
}

#endif
