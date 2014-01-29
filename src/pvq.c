/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

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

#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>
#include "logging.h"
#include <math.h>
#include <string.h>
#include "filter.h"

#define MAXN 256
#define EPSILON 1e-30

static const int od_layout16_offsets[4] = { 0, 32, 64, 192 };
extern const index_pair od_zigzag16[];
const band_layout od_layout16 = {
  od_zigzag16,
  16,
  3,
  od_layout16_offsets
};

const int od_layout8_offsets[4] = { 0, 8, 16, 48 };
extern const index_pair od_zigzag8[];
const band_layout od_layout8 = {
  od_zigzag8,
  8,
  3,
  od_layout8_offsets
};

static const int od_layout4_offsets[2] = { 0, 15 };
extern const index_pair od_zigzag4[];
const band_layout od_layout4 = {
  od_zigzag4,
  4,
  1,
  od_layout4_offsets
};

/* Table of combined "use prediction" pvq flags for 8x8 trained on
   subset1. Should eventually make this adaptive. */
const ogg_uint16_t pred8_cdf[16] = {
  22313, 22461, 22993, 23050, 23418, 23468, 23553, 23617,
  29873, 30181, 31285, 31409, 32380, 32525, 32701, 32768
};

const ogg_uint16_t pred16_cdf[16][8] = {
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 }
};

void od_bands_from_raster(const band_layout *layout, od_coeff *dst,
 od_coeff *src, int stride) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    dst[i] = src[layout->dst_table[i][1]*stride + layout->dst_table[i][0]];
  }
}

void od_raster_from_bands(const band_layout *layout, od_coeff *dst,
 int stride, od_coeff *src) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    dst[layout->dst_table[i][1]*stride + layout->dst_table[i][0]] = src[i];
  }
}

void od_band_pseudo_zigzag(od_coeff *dst,  int n, od_coeff *src, int stride,
 int interleave) {
  od_coeff tmp1[1024];
  od_bands_from_raster(&od_layout4, dst+1, src, stride);
  if (n >= 8) {
    int i;
    od_bands_from_raster(&od_layout8, dst+16, src, stride);
    if (interleave) {
      for (i = 0; i < 8; i++) {
        tmp1[2*i] = dst[16+i];
        tmp1[2*i+1] = dst[24+i];
      }
      for (i = 0; i < 16; i++) {
        dst[16+i] = tmp1[i];
      }
    }
  }
  if (n >= 16) {
    int i;
    od_bands_from_raster(&od_layout16, dst+64, src, stride);
    if (interleave) {
      for (i = 0; i < 32; i++) {
        tmp1[2*i] = dst[64+i];
        tmp1[2*i+1] = dst[96+i];
      }
      for (i = 0; i < 64; i++) {
        dst[64+i] = tmp1[i];
      }
    }
  }
  dst[0] = src[0];
}

void od_band_pseudo_dezigzag(od_coeff *dst,  int stride, od_coeff *src,
 int n, int interleave) {
  od_raster_from_bands(&od_layout4, dst, stride, src+1);
  if (n >= 8) {
    if (interleave) {
      int i;
      od_coeff tmp1[1024];
      for (i = 0; i < 16; i++) {
        tmp1[i] = src[16 + i];
      }
      for (i = 0; i < 8; i++) {
        src[16+i] = tmp1[2*i];
        src[24+i] = tmp1[2*i + 1];
      }
    }
    od_raster_from_bands(&od_layout8, dst, stride, src+16);
  }
  if (n >= 16) {
    if (interleave) {
      int i;
      od_coeff tmp1[1024];
      for (i = 0; i < 64; i++) {
        tmp1[i] = src[64 + i];
      }
      for (i = 0; i < 32; i++) {
        src[64+i] = tmp1[2*i];
        src[96+i] = tmp1[2*i + 1];
      }
    }
    od_raster_from_bands(&od_layout16, dst, stride, src+64);
  }
  dst[0] = src[0];
}

/* Double-precision PVQ search just to make sure our tests aren't limited
   by numerical accuracy.
   @param [in]      x      input vector to quantize
   @param [in]      n      number of dimensions
   @param [in]      k      number of pulses
   @param [out]     y      optimal codevector found
   @return                 cosine distance between x and y (between 0 and 1)
*/
static double pvq_search_double(const double *x, int n, int k, od_coeff *y) {
  int i, j;
  double xy;
  double yy;
  double X[1024];
  double xx;
  xx = xy = yy = 0;
  for (j = 0; j < n; j++) {
    X[j] = fabs(x[j]);
    xx += X[j]*X[j];
  }
  for (j = 0; j < n; j++) y[j] = 0;
  /* Search one pulse at a time */
  for (i = 0; i < k; i++) {
    int pos;
    double best_xy;
    double best_yy;
    pos = 0;
    best_xy = -10;
    best_yy = 1;
    for (j = 0; j < n; j++) {
      double tmp_xy;
      double tmp_yy;
      tmp_xy = xy + X[j];
      tmp_yy = yy + 2*y[j] + 1;
      tmp_xy *= tmp_xy;
      if (j == 0 || tmp_xy*best_yy > best_xy*tmp_yy) {
        best_xy = tmp_xy;
        best_yy = tmp_yy;
        pos = j;
      }
    }
    xy = xy + X[pos];
    yy = yy + 2*y[pos] + 1;
    y[pos]++;
  }
  for (i = 0; i < n; i++) {
    if (x[i] < 0) y[i] = -y[i];
  }
  return xy/(1e-100 + sqrt(xx*yy));
}

#define ACTIVITY (1.)

/* Computes Householder reflection that aligns the reference r to one
   of the dimensions (return value). The reflection vector is returned
   in r and x is reflected. */
int compute_householder(double *r, int n, double gr, int *sign) {
  int m;
  int i;
  int s;
  double maxr;
  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m = 0;
  maxr = 0;
  for (i = 0; i < n; i++) {
    if (fabs(r[i]) > maxr) {
      maxr = fabs(r[i]);
      m = i;
    }
  }
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "max r: %f %f %d", maxr, r[m], m));
  s = r[m] > 0 ? 1 : -1;
  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m] += gr*s;
  *sign = s;
  return m;
}

/* Applies Householder reflection from compute_householder(). The reflection
   is its own inverse. */
static void apply_householder(double *x, const double *r, int n) {
  int i;
  double proj;
  double proj_1;
  double l2r;
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += r[i]*r[i];
  }
  /* Apply Householder reflection */
  proj = 0;
  for (i = 0; i < n; i++) {
    proj += r[i]*x[i];
  }
  proj_1 = proj*2./(1e-100 + l2r);
  for (i = 0; i < n; i++) {
    x[i] -= r[i]*proj_1;
  }
}

/* Encodes the gain in such a way that the return value increases with the
   distance |x-ref|, so that we can encode a zero when x=ref. The value x=0
   is not covered because it is only allowed in the noref case. */
static int neg_interleave(int x, int ref) {
  if (x < ref) return -2*(x - ref) - 1;
  else if (x < 2*ref) return 2*(x - ref);
  else return x-1;
}

static int neg_deinterleave(int x, int ref) {
  if (x < 2*ref-1) {
    if (x & 1) return ref - 1 - (x >> 1);
    else return ref + (x >> 1);
  }
  else return x+1;
}

void pvq_synthesis(od_coeff *x0, od_coeff *y, const double *r, int n, int noref,
 int qg, double gain_offset, double theta, int m, int s, double q) {
  int i;
  int yy;
  double qcg;
  double norm;
  double x[MAXN];
  double g;
  if (noref) {
    qcg = qg;
    yy = 0;
    for (i = 0; i < n; i++) {
      yy += y[i]*(ogg_int32_t)y[i];
    }
    norm = sqrt(1./(1e-100 + yy));
    for (i = 0; i < n; i++) {
      x[i] = y[i]*norm;
    }
  }
  else {
    qcg = qg+gain_offset;
    if (qg == 0) qcg = 0;
    yy = 0;
    for (i = 0; i < n; i++) {
      yy += y[i]*(ogg_int32_t)y[i];
    }
    norm = sqrt(1./(1e-100 + yy));
    for (i = 0; i < n; i++) {
      x[i] = y[i]*norm*sin(theta);
    }
    x[m] = -s*cos(theta);
    apply_householder(x, r, n);
  }
  g = q*pow(qcg, 1./ACTIVITY);
  for (i = 0; i < n; i++) {
    x[i] *= g;
  }
  for (i = 0; i < n; i++) {
    x0[i] = floor(.5 + x[i]);
  }
}

/* This does PVQ quantization with prediction, trying several possible gains
   and angles. See draft-valin-videocodec-pvq and
   http://jmvalin.ca/slides/pvq.pdf for more details.
   @param [in,out] x0        coefficients being quantized (before and after)
   @param [in]     r0        reference, aka predicted coefficients
   @param [in]     n         number of dimensions
   @param [in]     q0        quantization step size
   @param [out]    y         pulse vector (i.e. selected PVQ codevector)
   @param [out]    itheta    angle between input and reference (-1 if noref)
   @param [out]    max_theta maximum value of itheta that could have been
   @param [out]    vk        total number of pulses
   @return         gain      index of the quatized gain
*/
int pvq_theta(od_coeff *x0, od_coeff *r0, int n, int q0, od_coeff *y, int *itheta,
 int *max_theta, int *vk) {
  double l2x;
  double l2r;
  double g;
  double gr;
  double x[MAXN];
  double r[MAXN];
  od_coeff y_tmp[MAXN];
  int i;
  /* Number of pulses. */
  int k;
  /* Companded gain of x and reference, normalized to q. */
  double cg;
  double cgr;
  int icgr;
  int qg;
  double best_dist;
  double q;
  /* Sign of Householder reflection. */
  int s;
  /* Dimension on which Householder reflects. */
  int m;
  double theta;
  double corr;
  int best_k;
  double best_qtheta;
  double gain_offset;
  /* Quantized companded gain */
  double qcg;
  int noref;
  double lambda;
  /* Normalized lambda. At high rate, this would be log(2)/6, but we're
     making RDO a bit less aggressive for now. */
  lambda = .05;
  q = q0;
  OD_ASSERT(n > 1);
  l2x = 0;
  l2r = 0;
  corr = 0;
  for (i = 0; i < n; i++) {
    x[i] = x0[i];
    r[i] = r0[i];
    l2x += x[i]*x[i];
    l2r += r[i]*r[i];
    corr += x[i]*r[i];
  }
  g = sqrt(l2x);
  gr = sqrt(l2r);
  /* Normalize gain by quantization step size and apply companding
     (if ACTIVITY != 1). */
  cg = pow(g/q, ACTIVITY);
  cgr = pow(gr/q, ACTIVITY);
  /* gain_offset is meant to make sure one of the quantized gains has
     exactly the same gain as the reference. */
  icgr = floor(.5+cgr);
  gain_offset = cgr-icgr;
  qg = 0;
  best_dist = 1e100;
  best_k = 0;
  best_qtheta = 0;
  noref = 0;
  m = 0;
  s = 1;
  if (corr > 0) {
    /* Perform theta search only if prediction is useful. */
    corr = corr/(1e-100+g*gr);
    corr = OD_MAXF(OD_MINF(corr, 1.), -1.);
    theta = acos(corr);
    m = compute_householder(r, n, gr, &s);
    apply_householder(x, r, n);
    x[m] = 0;
    /* Search for the best gain within a reasonable range. */
    for (i = OD_MAXI(1, (int)floor(cg-gain_offset)-1);
     i <= (int)ceil(cg-gain_offset); i++) {
      int j;
      int ts;
      qcg = i+gain_offset;
      if (i == 0) qcg = 0;
      /* Set angular resolution (in ra) to match the encoded gain */
      ts = (int)floor(.5 + qcg*M_PI/2);
      /* Special case for low gains -- will need to be tuned anyway */
      if (qcg < 1.4) ts = 1;
      /* Search for the best angle within a reasonable range. */
      for (j = OD_MAXI(0, (int)floor(.5+theta*2/M_PI*ts)-1);
       j <= OD_MINI(ts-1, (int)ceil(theta*2/M_PI*ts)); j++) {
        double cos_dist;
        double dist;
        double dist_theta;
        double qtheta;
        if (ts != 0) qtheta = j*.5*M_PI/ts;
        else qtheta = 0;
        /* Sets K according to gain and theta, based on the high-rate
           PVQ distortion curves D~=N^2/(24*K^2). Low-rate will have to be
           perceptually tuned anyway.  */
        k = floor(.5 + qcg*sin(qtheta)*sqrt(n/2));
        cos_dist = pvq_search_double(x, n, k, y_tmp);
        /* See Jmspeex' Journal of Dubious Theoretical Results. */
        dist_theta = 2 - 2*cos(theta - qtheta)
         + sin(theta)*sin(qtheta)*(2 - 2*cos_dist);
        dist = (qcg - cg)*(qcg - cg) + qcg*cg*dist_theta;
        /* Do approximate RDO -- should eventually be improved. */
        dist += lambda*log2(n)*k;
        if (j == 0) dist -= lambda*2.;
        if (i == icgr) dist -= lambda*2.;
        if (dist < best_dist) {
          best_dist = dist;
          qg = i;
          best_k = k;
          best_qtheta = qtheta;
          *itheta = j;
          *max_theta = ts;
          OD_COPY(y, y_tmp, n);
        }
      }
    }
  }
  /* Don't bother with no-reference version if there's a reasonable
     correlation */
  if (corr < .5 || cg < 2.) {
    double x1[MAXN];
    for (i = 0; i < n; i++) x1[i] = x0[i];
    /* Search for the best gain (haven't determined reasonable range yet). */
    for (i = 0; i <= ceil(cg); i++) {
      double cos_dist;
      double dist;
      qcg = i;
      /* See K from above. */
      k = floor(.5 + qcg*sqrt(n/2));
      cos_dist = pvq_search_double(x1, n, k, y_tmp);
      /* See Jmspeex' Journal of Dubious Theoretical Results. */
      dist = (qcg - cg)*(qcg - cg) + qcg*cg*(2 - 2*cos_dist);
      dist += lambda*(log2(n)*k-2.);
      if (dist <= best_dist) {
        best_dist = dist;
        qg = i;
        noref = 1;
        best_k = k;
        *itheta = -1;
        *max_theta = 0;
        OD_COPY(y, y_tmp, n);
      }
    }
  }
  k = best_k;
  theta = best_qtheta;
  /* Synthesize like the decoder would. */
  pvq_synthesis(x0, y, r, n, noref, qg, gain_offset, theta, m, s, q);
  /* Remove dimension m if we're using theta. */
  if (!noref) {
    for (i = m; i < n - 1; i++) y[i] = y[i+1];
  }
  *vk = k;
  /* Encode gain differently depending on whether we use prediction or not. */
  return noref ? qg : neg_interleave(qg, icgr);
}

int od_compute_max_theta(const od_coeff *r, int n, int q, double *gr,
 double *qcg, int *qg, double *gain_offset, int noref) {
  double l2r;
  double cgr;
  int icgr;
  int ts;
  int i;
  l2r=0;
  for (i = 0; i < n; i++) l2r += r[i]*r[i];
  *gr = sqrt(l2r);
  cgr = pow(*gr/q, ACTIVITY);
  icgr = floor(.5+cgr);
  *gain_offset = cgr-icgr;
  if (!noref) *qg = neg_deinterleave(*qg, icgr);
  *qcg = *qg;
  if (!noref) *qcg += *gain_offset;
  if (*qg == 0) *qcg = 0;
  if (noref) ts = 0;
  else {
    /* Set angular resolution (in ra) to match the encoded gain */
    ts = (int)floor(.5 + *qcg*M_PI/2);
    /* Special case for low gains -- will need to be tuned anyway */
    if (*qcg < 1.4) ts = 1;
  }
  return ts;
}

double od_compute_k_theta(int *k, double qcg, int t, int max_theta, int noref,
 int n) {
  double theta;
  if (noref) {
    *k = floor(.5 + qcg*sqrt(n/2));
    theta = -1;
  }
  else {
    if (max_theta != 0) theta = t*.5*M_PI/max_theta;
    else theta = 0;
    /* Sets K according to gain and theta, based on the high-rate
       PVQ distortion curves D~=N^2/(24*K^2). Low-rate will have to be
       perceptually tuned anyway.  */
    *k = (int)floor(.5 + qcg*sin(theta)*sqrt(n/2));
  }
  return theta;
}

