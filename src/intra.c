/*Daala video codec
Copyright (c) 2012-2013 Daala project contributors.  All rights reserved.

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
#include <limits.h>
#include "block_size.h"
#include "filter.h"
#include "intra.h"
#include "tf.h"
#include "state.h"

const od_intra_mult_func OD_INTRA_MULT[OD_NBSIZES+1] = {
  od_intra_pred4x4_mult,
  od_intra_pred8x8_mult,
  od_intra_pred16x16_mult
};

void od_intra_pred4x4_mult(double *pred, int pred_stride, od_coeff *blocks[4],
 int strides[4], int mode) {
  const ogg_uint16_t *index;
  const double *weights;
  int j;
  int i;
  int k;
  index = OD_PRED_INDEX_4x4 + OD_PRED_OFFSETS_4x4[mode];
  weights = OD_PRED_WEIGHTS_4x4 + OD_PRED_OFFSETS_4x4[mode];
  for (j = 0; j < 4; j++) {
    for (i = 0; i < 4; i++) {
      double sum;
      sum = 0;
      for (k = OD_PRED_MULTS_4x4[mode][j][i]; k-- > 0;) {
        int id;
        int x;
        int y;
        id = *index;
        x = id & 0x3;
        id >>= 2;
        y = id & 0x3;
        id >>= 2;
        sum += blocks[id][strides[id]*y + x]*(*weights);
        index++;
        weights++;
      }
      pred[pred_stride*j + i] = sum;
    }
  }
}

void od_intra_pred8x8_mult(double *pred, int pred_stride, od_coeff *blocks[4],
 int strides[4], int mode) {
  const ogg_uint16_t *index;
  const double *weights;
  int j;
  int i;
  int k;
  index = OD_PRED_INDEX_8x8 + OD_PRED_OFFSETS_8x8[mode];
  weights = OD_PRED_WEIGHTS_8x8 + OD_PRED_OFFSETS_8x8[mode];
  for (j = 0; j < 8; j++) {
    for (i = 0; i < 8; i++) {
      double sum;
      sum = 0;
      for (k = OD_PRED_MULTS_8x8[mode][j][i]; k-- > 0;) {
        int id;
        int x;
        int y;
        id = *index;
        x = id & 0x7;
        id >>= 3;
        y = id & 0x7;
        id >>= 3;
        sum += blocks[id][strides[id]*y + x]*(*weights);
        index++;
        weights++;
      }
      pred[pred_stride*j + i] = sum;
    }
  }
}

void od_intra_pred16x16_mult(double *pred, int pred_stride,
 od_coeff *blocks[4], int strides[4], int mode) {
  const ogg_uint16_t *index;
  const double *weights;
  int j;
  int i;
  int k;
  index = OD_PRED_INDEX_16x16 + OD_PRED_OFFSETS_16x16[mode];
  weights = OD_PRED_WEIGHTS_16x16 + OD_PRED_OFFSETS_16x16[mode];
  for (j = 0; j < 16; j++) {
    for (i = 0; i < 16; i++) {
      double sum;
      sum = 0;
      for (k = OD_PRED_MULTS_16x16[mode][j][i]; k-- > 0;) {
        int id;
        int x;
        int y;
        id = *index;
        x = id & 0xf;
        id >>= 4;
        y = id & 0xf;
        id >>= 4;
        sum += blocks[id][strides[id]*y + x]*(*weights);
        index++;
        weights++;
      }
      pred[pred_stride*j + i] = sum;
    }
  }
}

const od_intra_dist_func OD_INTRA_DIST[OD_NBSIZES+1] = {
  od_intra_pred4x4_dist,
  od_intra_pred8x8_dist,
  od_intra_pred16x16_dist
};

void od_intra_pred4x4_dist(ogg_uint32_t *dist, const od_coeff *c,
 int stride, od_coeff *neighbors[4], int neighbor_strides[4]) {
  double p[4*4];
  ogg_uint32_t satd;
  int mode;
  int i;
  int j;
  for (mode = 0; mode < OD_INTRA_NMODES; mode++) {
    od_intra_pred4x4_mult(p, 4, neighbors, neighbor_strides, mode);
    satd = 0;
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        satd += abs(floor(c[stride*i + j] - p[i*4 + j] + 0.5))*
         OD_SATD_WEIGHTS_4x4[i*4 + j];
      }
    }
    dist[mode] = satd;
  }
}

void od_intra_pred8x8_dist(ogg_uint32_t *dist, const od_coeff *c,
 int stride, od_coeff *neighbors[4], int neighbor_strides[4]) {
  double p[8*8];
  ogg_uint32_t satd;
  int mode;
  int i;
  int j;
  for (mode = 0; mode < OD_INTRA_NMODES; mode++) {
    od_intra_pred8x8_mult(p, 8, neighbors, neighbor_strides, mode);
    satd = 0;
    for (i = 0; i < 8; i++) {
      for (j = 0; j < 8; j++) {
        satd += abs(floor(c[stride*i + j] - p[i*8 + j] + 0.5))*
         OD_SATD_WEIGHTS_8x8[i*8 + j];
      }
    }
    dist[mode] = satd;
  }
}

void od_intra_pred16x16_dist(ogg_uint32_t *dist, const od_coeff *c,
 int stride, od_coeff *neighbors[4], int neighbor_strides[4]) {
  double p[16*16];
  ogg_uint32_t satd;
  int mode;
  int i;
  int j;
  for (mode = 0; mode < OD_INTRA_NMODES; mode++) {
    od_intra_pred16x16_mult(p, 16, neighbors, neighbor_strides, mode);
    satd = 0;
    for (i = 0; i < 16; i++) {
      for (j = 0; j < 16; j++) {
        satd += abs(floor(c[stride*i + j] - p[i*16 + j] + 0.5))*
         OD_SATD_WEIGHTS_16x16[i*16 + j];
      }
    }
    dist[mode] = satd;
  }
}

/*These are used to weight the samples we use to build the chroma-from-luma
   model.
  The weights are a result of minimizing the mean squared error for DC on
   subset1/3 with 10% of the outliers removed, subject to the weights
   adding to 1 and having no values less than 1/256.
  We look up the weights using the INTRA modes of the corresponding luma
   blocks.*/
const int OD_INTRA_CHROMA_WEIGHTS_Q8[OD_INTRA_NMODES][3] = {
  { 1,  89, 166},
  { 1,  45, 210},
  { 1,  19, 236},
  { 1,   1, 254},
  {27,  17, 212},
  {25,  92, 139},
  {31, 181,  44},
  { 1, 254,   1},
  { 1, 231,  24},
  { 1, 168,  87},
};

void od_chroma_pred(od_coeff *p, const od_coeff *c, const od_coeff *l,
 int stride, int bx, int by, int ln, int xdec, int ydec,
 const unsigned char *bsize, int bstride, const int weights_q8[3]) {
  static const int BLOCK_DX[3] = { -1,  0, -1 };
  static const int BLOCK_DY[3] = { -1, -1, 0 };
  static const int AC_DX[3] = { 1, 0, 1 };
  static const int AC_DY[3] = { 0, 1, 1 };
  ogg_int64_t xx;
  ogg_int64_t xy;
  ogg_int32_t alpha_q8;
  ogg_int32_t beta_q8;
  ogg_int32_t lc_sum_q8;
  ogg_int32_t cc_sum_q8;
  int oshift;
  int cshift;
  int bi;
  int i;
  int j;
  int n;
  /*Solve for a simple predictive model using the UL, U, L neighbors:
    chroma DC = luma DC * alpha + beta
    chroma AC = luma AC * alpha*/
  lc_sum_q8 = cc_sum_q8 = 0;
  xx = xy = 0;
  /*Because we may be predicting from neighbors of different sizes we must
     shift to a common scale.
    Because there is enough dynamic range we use the largest scale plus Q8.*/
  oshift = 11 - ln;
  cshift = 3 - ln;
  OD_ASSERT(xdec == ydec);
  OD_ASSERT(cshift - xdec >= 0);
  for (bi = 0; bi < 3; bi++) {
    od_coeff lc;
    od_coeff cc;
    int nx;
    int ny;
    int nsize;
    int nshift;
    int boffs;
    int ci;
    int w_q8;
    w_q8 = weights_q8[bi];
    nx = bx + BLOCK_DX[bi];
    ny = by + BLOCK_DY[bi];
    nsize = OD_BLOCK_SIZE4x4(bsize, bstride, nx << xdec, ny << ydec);
    nsize = OD_MAXI(nsize - xdec, 0);
    nx = nx >> nsize << nsize;
    ny = ny >> nsize << nsize;
    nshift = 3 - nsize;
    boffs = (ny << 2)*stride + (nx << 2);
    /*Resampled prediction is scaled by a factor of 2.*/
    lc = *(l + boffs) << (nshift - xdec);
    cc = *(c + boffs) << nshift;
    xx += lc*(ogg_int64_t)lc*w_q8;
    xy += lc*(ogg_int64_t)cc*w_q8;
    lc_sum_q8 += lc*w_q8;
    cc_sum_q8 += cc*w_q8;
    for (ci = 0; ci < 3; ci++) {
      int coffs;
      coffs = boffs + AC_DY[ci]*stride + AC_DX[ci];
      lc = *(l + coffs) << (nshift - xdec);
      cc = *(c + coffs) << nshift;
      xx += lc*(ogg_int64_t)lc*w_q8;
      xy += lc*(ogg_int64_t)cc*w_q8;
    }
  }
  l += (by << 2)*stride + (bx << 2);
  xx -= (lc_sum_q8*(ogg_int64_t)lc_sum_q8 + 128) >> 8;
  xy -= (cc_sum_q8*(ogg_int64_t)lc_sum_q8 + 128) >> 8;
  if (abs(xx) > abs(xy) >> 1) alpha_q8 = (xy << 8)/xx;
  else alpha_q8 = 0;
  if (abs(alpha_q8) > 128) alpha_q8 = 0;
  beta_q8 = cc_sum_q8 - ((alpha_q8*lc_sum_q8 + 128) >> 8);
  /*Alpha is scaled by the amount needed to bring the luma at this block to
    the working scale.*/
  alpha_q8 <<= cshift - xdec;
  p[0] = (l[0]*alpha_q8 + beta_q8 + (1 << (oshift - 1))) >> oshift;
  n = 1 << (2 + ln);
  for (i = 0; i < n; i++) {
    for (j = i == 0; j < n; j++) {
      p[i*n + j] = (l[i*stride+j]*alpha_q8 + (1 << (oshift - 1))) >> oshift;
    }
  }
}

const od_intra_get_func OD_INTRA_GET[OD_NBSIZES+1] = {
  od_intra_pred4x4_get,
  od_intra_pred8x8_get,
  od_intra_pred16x16_get
};

void od_intra_pred4x4_get(od_coeff *_out,
 od_coeff *_neighbors[4], int _neighbor_strides[4], int _mode) {
  double p[4*4];
  int    i;
  int    j;
  od_intra_pred4x4_mult(p, 4, _neighbors, _neighbor_strides, _mode);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      _out[4*i+j] = (od_coeff)floor(p[i*4+j]+0.5);
    }
  }
}

void od_intra_pred8x8_get(od_coeff *_out,
 od_coeff *_neighbors[4], int _neighbor_strides[4], int _mode) {
  double p[8*8];
  int    i;
  int    j;
  od_intra_pred8x8_mult(p, 8, _neighbors, _neighbor_strides, _mode);
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      _out[8*i+j] = (od_coeff)floor(p[i*8+j]+0.5);
    }
  }
}

void od_intra_pred16x16_get(od_coeff *_out,
 od_coeff *_neighbors[4], int _neighbor_strides[4], int _mode) {
  double p[16*16];
  int    i;
  int    j;
  od_intra_pred16x16_mult(p, 16, _neighbors, _neighbor_strides, _mode);
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      _out[16*i+j] = (od_coeff)floor(p[i*16+j]+0.5);
    }
  }
}

void od_intra_pred_cdf(ogg_uint16_t _cdf[],
 unsigned char _probs[][OD_INTRA_NCONTEXTS], int _nmodes,
 int _left, int _upleft, int _up) {
  unsigned p[OD_INTRA_NMODES+1];
  int      mi;
  int      sum;
  int      curr_cdf;
  sum = 0;
  for (mi = 0; mi < _nmodes; mi++) {
    p[mi] = _probs[mi][(_left == mi)*4+(_upleft == mi)*2+(_up == mi)];
    p[mi] = OD_MAXI(p[mi], 1);
    sum += p[mi];
  }
  curr_cdf = 0;
  for (mi = 0; mi < _nmodes; mi++) {
    /*Apply probability combination here: p[mi] *= (sum-p[mi])/(1-p[mi]).*/
    /*FIXME: Make this fixed-point.*/
    p[mi] = OD_MINI(8192, OD_MAXI(1,
     (int)(p[mi]*(sum-p[mi])/(float)(256-p[mi]))));
    curr_cdf += p[mi];
    _cdf[mi] = curr_cdf;
  }
}

int od_intra_pred_search(const ogg_uint16_t _cdf[],
 const ogg_uint32_t _dist[], int _nmodes, ogg_uint16_t _lambda) {
  int best_score;
  int best_mode;
  int mi;
  /*FIXME: Compute the log2() in fixed-point.*/
  best_score = _dist[0]-_lambda*OD_LOG2(_cdf[0]);
  best_mode = 0;
  for (mi = 1; mi < _nmodes; mi++) {
    int score;
    /*FIXME: Compute the log2() in fixed-point.*/
    score = _dist[mi]-_lambda*OD_LOG2(_cdf[mi]-_cdf[mi-1]);
    if (score < best_score) {
      best_score = score;
      best_mode = mi;
    }
  }
  return best_mode;
}

void od_intra_pred_update(unsigned char _probs[][OD_INTRA_NCONTEXTS],
 int _nmodes, int _mode, int _left, int _upleft, int _up) {
  int mi;
  for (mi = 0; mi < _nmodes; mi++) {
    int id = (_left == mi)*4+(_upleft == mi)*2+(_up == mi);
    _probs[mi][id] -= _probs[mi][id] >> OD_INTRA_ADAPT_SPEED;
    if (_mode == mi) {
      _probs[mi][id] = _probs[mi][id] + (256 >> OD_INTRA_ADAPT_SPEED);
    }
    /* Bound certainty to avoid problems with the renormalization. Also
       avoids high costs if the modeling isn't very good. */
    _probs[mi][id] = OD_MAXI(5, OD_MINI(251, _probs[mi][id]));
  }
}

void od_resample_luma_coeffs(od_coeff *l, int lstride,
 const od_coeff *c, int cstride, int xdec, int ydec, int ln, int cln) {
  int n;
  n = 4 << ln;
  if (cln == 0 && (xdec || ydec)) {
    if (xdec) {
      if (ydec) od_tf_up_hv_lp(l, lstride, c, cstride, n, n, n);
      else od_tf_up_h_lp(l, lstride, c, cstride, n, n);
    }
    else {
      OD_ASSERT(ydec);
      od_tf_up_v_lp(l, lstride, c, cstride, n, n);
    }
  }
  else {
    /*When the transform we code chroma with is smaller than the luma one,
       downsampling just requires copying the upper lower quarter coeffs.*/
    int x;
    int y;
    OD_ASSERT(xdec == ydec);
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        l[y*lstride + x] = c[y*cstride + x];
      }
    }
  }
}
