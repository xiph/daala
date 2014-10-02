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

void od_hv_intra_pred(od_coeff *pred, od_coeff *d, int w, int bx, int by,
 unsigned char *bsize, int bstride, int ln) {
  int i;
  od_coeff *t;
  double g1;
  double g2;
  int top;
  int left;
  int n;
  n = 1 << (ln + OD_LOG_BSIZE0);
  top = by > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx, by - 1) == ln;
  left = bx > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx - 1, by) == ln;
  t = &d[((by << OD_LOG_BSIZE0))*w + (bx << OD_LOG_BSIZE0)];
  g1 = g2 = 0;
  if (top) for (i = 1; i < 4; i++) g1 += t[-n*w + i]*(double)t[-n*w + i];
  if (left) for (i = 1; i < 4; i++) g2 += t[-n + i*w]*(double)t[-n + i*w];
  if (top) for (i = 4; i < n; i++) pred[i] = t[-n*w + i];
  if (left) for (i = 4; i < n; i++) pred[i*n] = t[-n + i*w];
  if (g1 > g2) {
    if (top) for (i = 1; i < 4; i++) pred[i] = t[-n*w + i];
  }
  else {
    if (left) for (i = 1; i < 4; i++) pred[i*n] = t[-n + i*w];
  }
}

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
