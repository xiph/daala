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
  predict_bottom_edge(p, edge_accum, n, stride, m, has_right);
#if !QUANTIZE
  /*printf("%d ", m);*/
  for (i = 0; i < n; i++) printf("%d ", edge_accum[n*stride + (i + 1)] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[n*stride + i + 1] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) x[i] = (int)(floor(.5+x[i]/q));
  encode_edge_coeffs(adapt, enc, x, n, has_right || (m >= n && m <= 3*n));
  for (i = 0; i < n; i++) x[i] = q*x[i];
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[n*stride + i + 1] = OD_MAXI(0, OD_MINI(255, r[i]));
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
  predict_right_edge(p, edge_accum, n, stride, m, has_bottom);
#if !QUANTIZE
  for (i = 0; i < n; i++) printf("%d ", edge_accum[(i + 1)*stride + n] - p[i]);
#endif
  for (i = 0; i < n; i++) r[i] = edge_accum[(i + 1)*stride + n] - p[i];
  my_fdct_table[lsize](x, r, 1);
#if QUANTIZE
  for (i = 0; i < n; i++) x[i] = (int)(floor(.5+x[i]/q));
  encode_edge_coeffs(adapt, enc, x, n, (has_bottom || (m >= n && m <= 3*n)) + 2);
  for (i = 0; i < n; i++) x[i] = q*x[i];
#endif
  my_idct_table[lsize](r, 1, x);
  for (i = 0; i < n; i++) r[i] += p[i];
  for (i = 0; i < n; i++) edge_accum[(i+1)*stride+n] = OD_MAXI(0, OD_MINI(255, r[i]));
}

/* Quantize both the right and bottom edge, changing the order to maximize
   the number of pixels we can predict. */
void quantize_edge(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *edge_accum, int n, int stride, int q,
 int m, int dc_quant) {
  /*printf("q\%d ", n);*/
  if (dc_quant) {
    int pred;
    int res;
    int qdc;
    int i;
    pred = floor(.5 + .73*(edge_accum[n] + edge_accum[n*stride]) - .46 *edge_accum[0]);
    res = edge_accum[n*stride+n]-pred;
    qdc = OD_MAXI(1, q/8);
    res = (int)floor(.5+res/qdc);
    generic_encode(enc, &adapt->paint_dc_model, abs(res), -1, &ex_dc, 6);
    if (res != 0) od_ec_enc_bits(enc, res > 0, 1);
    /*printf("DC %d\n", res);*/
    res = res*qdc;
    edge_accum[n*stride+n] = res + pred;
    for (i = 1; i < n; i++) {
      edge_accum[n*stride+i] = edge_accum[n*stride] + i*(edge_accum[n*stride+n]-edge_accum[n*stride])/n;
      edge_accum[i*stride+n] = edge_accum[n] + i*(edge_accum[n*stride+n]-edge_accum[n])/n;
    }
  }
  else if (m > 0 && m < n) {
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
    compute_edges(&img[stride*n*by + n*bx], &edge_sum[stride*n*by + n*bx],
     &edge_count[stride*n*by + n*bx], n, stride, curr_mode);
    mode[(by<<bs)*mstride + (bx<<bs)] = curr_mode;
    for (k=0;k<1<<bs;k++) {
      for (m=0;m<1<<bs;m++) {
        mode[((by << bs) + k)*mstride + (bx << bs) + m] = curr_mode;
      }
    }
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
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res, int bx, int by, int level) {
  int bs;
  bs = dec8[(by<<level>>1)*bstride + (bx<<level>>1)];

  OD_ASSERT(bs <= level);
  if (bs < level) {
    level--;
    bx <<= 1;
    by <<= 1;
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx, by + 1, level);
    od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride,
     mode, mstride, edge_sum, edge_count, q, res, bx + 1, by + 1, level);
  }
  else {
    int ln;
    int n;
    int k;
    int idx;
    ln = 2 + bs;
    n = 1 << ln;
    /* Compute left edge (left column only). */
    if (bx == 0) {
      for (k = 0; k < n; k++) {
        idx = stride*(n*by + k) - 1;
        if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
        else paint[idx] = img[idx];
      }
      quantize_initial_edge(adapt, enc, &paint[stride*(n*by + 1)], ln, stride, q);
    }
    /* Compute top edge (top row only). */
    if (by == 0) {
      for (k = 0; k < n; k++) {
        idx = -stride + n*bx + k;
        if (edge_count[idx] > 0) paint[idx] = edge_sum[idx]/edge_count[idx];
        else paint[idx] = img[idx];
      }
      quantize_initial_edge(adapt, enc, &paint[n*bx + 1], ln, 1, q);
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
    od_intra_paint_mode_encode(enc, mode, bx, by, ln, mstride, dec8, bstride, res);
    /* Only use special DC quantization when the two adjacent blocks are
       also DC. In the future, we could also treat each edge separately. */
    dc_quant = mode[(by*mstride + bx) << ln >> 2]==4*n
     && mode[(by*mstride + bx + 1) << ln >> 2] == 4*n
     && mode[((by + 1)*mstride + bx) << ln >> 2] == 4*n;
    /*quantize_edge(adapt, enc, &paint[stride*n*by + n*bx], n, stride, q,
     mode[(by*mstride + bx) << ln >> 2], dc_quant);*/
    interp_block(&paint[stride*n*by + n*bx], &paint[stride*n*by + n*bx],
     n, stride, mode[(by*mstride + bx) << ln >> 2]);
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
      printf("%d ", mode[i*mstride+j]<<3>>dec8[(i>>1)*bstride + (j>>1)]);
    }
    printf("\n");
  }
#endif

  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_intra_paint_compute_edges(adapt, enc, paint, img, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_count, q, res, j, i, 3);
    }
  }
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      od_intra_paint_quant_block(adapt, enc, paint, img, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_count, q, res, j, i, 3);
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
