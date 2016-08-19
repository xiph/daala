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

#include <stdlib.h>
#include <stdio.h>
#include "internal.h"
#include "block_size.h"
#include "entenc.h"
#include "entcode.h"
#include "filter.h"
#include "pvq_encoder.h"
#include "partition.h"

#define OD_PVQ_RATE_APPROX (0)
/*Shift to ensure that the upper bound (i.e. for the max blocksize) of the
   dot-product of the 1st band of chroma with the luma ref doesn't overflow.*/
#define OD_CFL_FLIP_SHIFT (OD_LIMIT_BSIZE_MAX + 0)

static void od_encode_pvq_codeword(od_ec_enc *ec, od_pvq_codeword_ctx *adapt,
 const od_coeff *in, int n, int k) {
  int i;
  od_encode_band_pvq_splits(ec, adapt, in, n, k, 0);
  for (i = 0; i < n; i++) if (in[i]) od_ec_enc_bits(ec, in[i] < 0, 1);
}

/* Computes 1/sqrt(i) using a table for small values. */
static double od_rsqrt_table(int i) {
  static double table[16] = {
    1.000000, 0.707107, 0.577350, 0.500000,
    0.447214, 0.408248, 0.377964, 0.353553,
    0.333333, 0.316228, 0.301511, 0.288675,
    0.277350, 0.267261, 0.258199, 0.250000};
  if (i <= 16) return table[i-1];
  else return 1./sqrt(i);
}

/*Computes 1/sqrt(start+2*i+1) using a lookup table containing the results
   where 0 <= i < table_size.*/
static double od_custom_rsqrt_dynamic_table(const double* table,
 const int table_size, const double start, const int i) {
  if (i < table_size) return table[i];
  else return od_rsqrt_table(start + 2*i + 1);
}

/*Fills tables used in od_custom_rsqrt_dynamic_table for a given start.*/
static void od_fill_dynamic_rqrt_table(double *table, const int table_size,
 const double start) {
  int i;
  for (i = 0; i < table_size; i++)
    table[i] = od_rsqrt_table(start + 2*i + 1);
}

/** Find the codepoint on the given PSphere closest to the desired
 * vector. Double-precision PVQ search just to make sure our tests
 * aren't limited by numerical accuracy.
 *
 * @param [in]      xcoeff  input vector to quantize (x in the math doc)
 * @param [in]      n       number of dimensions
 * @param [in]      k       number of pulses
 * @param [out]     ypulse  optimal codevector found (y in the math doc)
 * @param [out]     g2      multiplier for the distortion (typically squared
 *                          gain units)
 * @param [in] pvq_norm_lambda enc->pvq_norm_lambda for quantized RDO
 * @param [in]      prev_k  number of pulses already in ypulse that we should
 *                          reuse for the search (or 0 for a new search)
 * @return                  cosine distance between x and y (between 0 and 1)
 */
static double pvq_search_rdo_double(const od_val16 *xcoeff, int n, int k,
 od_coeff *ypulse, double g2, double pvq_norm_lambda, int prev_k) {
  int i, j;
  double xy;
  double yy;
  /* TODO - This blows our 8kB stack space budget and should be fixed when
   converting PVQ to fixed point. */
  double x[MAXN];
  double xx;
  double lambda;
  double norm_1;
  int rdo_pulses;
  double delta_rate;
  xx = xy = yy = 0;
  for (j = 0; j < n; j++) {
    x[j] = fabs((float)xcoeff[j]);
    xx += x[j]*x[j];
  }
  norm_1 = 1./sqrt(1e-30 + xx);
  lambda = pvq_norm_lambda/(1e-30 + g2);
  i = 0;
  if (prev_k > 0 && prev_k <= k) {
    /* We reuse pulses from a previous search so we don't have to search them
       again. */
    for (j = 0; j < n; j++) {
      ypulse[j] = abs(ypulse[j]);
      xy += x[j]*ypulse[j];
      yy += ypulse[j]*ypulse[j];
      i += ypulse[j];
    }
  }
  else if (k > 2) {
    double l1_norm;
    double l1_inv;
    l1_norm = 0;
    for (j = 0; j < n; j++) l1_norm += x[j];
    l1_inv = 1./OD_MAXF(l1_norm, 1e-100);
    for (j = 0; j < n; j++) {
      double tmp;
      tmp = k*x[j]*l1_inv;
      ypulse[j] = OD_MAXI(0, (int)floor(tmp));
      xy += x[j]*ypulse[j];
      yy += ypulse[j]*ypulse[j];
      i += ypulse[j];
    }
  }
  else OD_CLEAR(ypulse, n);

  /* Only use RDO on the last few pulses. This not only saves CPU, but using
     RDO on all pulses actually makes the results worse for reasons I don't
     fully understand. */
  rdo_pulses = 1 + k/4;
  /* Rough assumption for now, the last position costs about 3 bits more than
     the first. */
  delta_rate = 3./n;
  /* Search one pulse at a time */
  for (; i < k - rdo_pulses; i++) {
    int pos;
    double best_xy;
    double best_yy;
    pos = 0;
    best_xy = -10;
    best_yy = 1;
    for (j = 0; j < n; j++) {
      double tmp_xy;
      double tmp_yy;
      tmp_xy = xy + x[j];
      tmp_yy = yy + 2*ypulse[j] + 1;
      tmp_xy *= tmp_xy;
      if (j == 0 || tmp_xy*best_yy > best_xy*tmp_yy) {
        best_xy = tmp_xy;
        best_yy = tmp_yy;
        pos = j;
      }
    }
    xy = xy + x[pos];
    yy = yy + 2*ypulse[pos] + 1;
    ypulse[pos]++;
  }
  /* Search last pulses with RDO. Distortion is D = (x-y)^2 = x^2 - 2*x*y + y^2
     and since x^2 and y^2 are constant, we just maximize x*y, plus a
     lambda*rate term. Note that since x and y aren't normalized here,
     we need to divide by sqrt(x^2)*sqrt(y^2). */
  for (; i < k; i++) {
    double rsqrt_table[4];
    int rsqrt_table_size = 4;
    int pos;
    double best_cost;
    pos = 0;
    best_cost = -1e5;
    /*Fill the small rsqrt lookup table with inputs relative to yy.
      Specifically, the table of n values is filled with
       rsqrt(yy + 1), rsqrt(yy + 2 + 1) .. rsqrt(yy + 2*(n-1) + 1).*/
    od_fill_dynamic_rqrt_table(rsqrt_table, rsqrt_table_size, yy);
    for (j = 0; j < n; j++) {
      double tmp_xy;
      double tmp_yy;
      tmp_xy = xy + x[j];
      /*Calculate rsqrt(yy + 2*ypulse[j] + 1) using an optimized method.*/
      tmp_yy = od_custom_rsqrt_dynamic_table(rsqrt_table, rsqrt_table_size,
       yy, ypulse[j]);
      tmp_xy = 2*tmp_xy*norm_1*tmp_yy - lambda*j*delta_rate;
      if (j == 0 || tmp_xy > best_cost) {
        best_cost = tmp_xy;
        pos = j;
      }
    }
    xy = xy + x[pos];
    yy = yy + 2*ypulse[pos] + 1;
    ypulse[pos]++;
  }
  for (i = 0; i < n; i++) {
    if (xcoeff[i] < 0) ypulse[i] = -ypulse[i];
  }
  return xy/(1e-100 + sqrt(xx*yy));
}

/** Encodes the gain so that the return value increases with the
 * distance |x-ref|, so that we can encode a zero when x=ref. The
 * value x=0 is not covered because it is only allowed in the noref
 * case.
 *
 * @param [in]      x      quantized gain to encode
 * @param [in]      ref    quantized gain of the reference
 * @return                 interleave-encoded quantized gain value
 */
static int neg_interleave(int x, int ref) {
  if (x < ref) return -2*(x - ref) - 1;
  else if (x < 2*ref) return 2*(x - ref);
  else return x-1;
}

int od_vector_is_null(const od_coeff *x, int len) {
  int i;
  for (i = 0; i < len; i++) if (x[i]) return 0;
  return 1;
}

static double od_pvq_rate(int qg, int icgr, int theta, int ts,
 const od_adapt_ctx *adapt, const od_coeff *y0, int k, int n,
 int is_keyframe, int pli, int speed) {
  double rate;
  if (k == 0) rate = 0;
  else if (speed > 0) {
    int i;
    int sum;
    double f;
    /* Compute "center of mass" of the pulse vector. */
    sum = 0;
    for (i = 0; i < n - (theta != -1); i++) sum += i*abs(y0[i]);
    f = sum/(double)(k*n);
    /* Estimates the number of bits it will cost to encode K pulses in
       N dimensions based on hand-tuned fit for bitrate vs K, N and
       "center of mass". */
    rate = (1 + .4*f)*n*OD_LOG2(1 + OD_MAXF(0, log(n*2*(1*f + .025))*k/n)) + 3;
  }
  else {
    od_ec_enc ec;
    od_pvq_codeword_ctx cd;
    int tell;
    od_ec_enc_init(&ec, 1000);
    OD_COPY(&cd, &adapt->pvq.pvq_codeword_ctx, 1);
    tell = od_ec_enc_tell_frac(&ec);
    od_encode_pvq_codeword(&ec, &cd, y0, n - (theta != -1), k);
    rate = (od_ec_enc_tell_frac(&ec)-tell)/8.;
    od_ec_enc_clear(&ec);
  }
  if (qg > 0 && theta >= 0) {
    /* Approximate cost of entropy-coding theta */
    rate += .9*OD_LOG2(ts);
    /* Adding a cost to using the H/V pred because it's going to be off
       most of the time. Cost is optimized on subset1, while making
       sure we don't hurt the checkerboard image too much.
       FIXME: Do real RDO instead of this arbitrary cost. */
    if (is_keyframe && pli == 0) rate += 6;
    if (qg == icgr) rate -= .5;
  }
  return rate;
}

/* This stores the information about a PVQ search candidate, so we can sort
   based on K. */
typedef struct {
  int gain;
  int k;
  od_val32 qtheta;
  int theta;
  int ts;
  od_val32 qcg;
} pvq_search_item;

int items_compare(pvq_search_item *a, pvq_search_item *b) {
  return a->k - b->k;
}

/** Perform PVQ quantization with prediction, trying several
 * possible gains and angles. See draft-valin-videocodec-pvq and
 * http://jmvalin.ca/slides/pvq.pdf for more details.
 *
 * @param [out]    out       coefficients after quantization
 * @param [in]     x0        coefficients before quantization
 * @param [in]     r0        reference, aka predicted coefficients
 * @param [in]     n         number of dimensions
 * @param [in]     q0        quantization step size
 * @param [out]    y         pulse vector (i.e. selected PVQ codevector)
 * @param [out]    itheta    angle between input and reference (-1 if noref)
 * @param [out]    max_theta maximum value of itheta that could have been
 * @param [out]    vk        total number of pulses
 * @param [in]     beta      per-band activity masking beta param
 * @param [out]    skip_diff distortion cost of skipping this block
 *                           (accumulated)
 * @param [in]     robust    make stream robust to error in the reference
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [in]     pli       plane index
 * @param [in]     adapt     probability adaptation context
 * @param [in]     qm        QM with magnitude compensation
 * @param [in]     qm_inv    Inverse of QM with magnitude compensation
 * @param [in] pvq_norm_lambda enc->pvq_norm_lambda for quantized RDO
 * @return         gain      index of the quatized gain
*/
static int pvq_theta(od_coeff *out, const od_coeff *x0, const od_coeff *r0,
 int n, int q0, od_coeff *y, int *itheta, int *max_theta, int *vk,
 double beta, double *skip_diff, int robust, int is_keyframe, int pli,
 const od_adapt_ctx *adapt, const int16_t *qm,
 const int16_t *qm_inv, double pvq_norm_lambda, int speed) {
  od_val32 g;
  od_val32 gr;
  od_coeff y_tmp[MAXN];
  int i;
  /* Number of pulses. */
  int k;
  /* Companded gain of x and reference, normalized to q. */
  od_val32 cg;
  od_val32 cgr;
  int icgr;
  int qg;
  /* Best RDO cost (D + lamdba*R) so far. */
  double best_cost;
  double dist0;
  /* Distortion (D) that corresponds to the best RDO cost. */
  double best_dist;
  double dist;
  /* Sign of Householder reflection. */
  int s;
  /* Dimension on which Householder reflects. */
  int m;
  od_val32 theta;
  double corr;
  int best_k;
  od_val32 best_qtheta;
  od_val32 gain_offset;
  int noref;
  double skip_dist;
  int cfl_enabled;
  int skip;
  double gain_weight;
  od_val16 x16[MAXN];
  od_val16 r16[MAXN];
  int xshift;
  int rshift;
  /* Give more weight to gain error when calculating the total distortion. */
  gain_weight = 1.4;
  OD_ASSERT(n > 1);
  corr = 0;
#if !defined(OD_FLOAT_PVQ)
  /* Shift needed to make x fit in 16 bits even after rotation.
     This shift value is not normative (it can be changed without breaking
     the bitstream) */
  xshift = OD_MAXI(0, od_vector_log_mag(x0, n) - 15);
  /* Shift needed to make the reference fit in 15 bits, so that the Householder
     vector can fit in 16 bits.
     This shift value *is* normative, and has to match the decoder. */
  rshift = OD_MAXI(0, od_vector_log_mag(r0, n) - 14);
#else
  xshift = 0;
  rshift = 0;
#endif
  for (i = 0; i < n; i++) {
#if defined(OD_FLOAT_PVQ)
    /*This is slightly different from the original float PVQ code,
       where the qm was applied in the accumulation in od_pvq_compute_gain and
       the vectors were od_coeffs, not od_val16 (i.e. double).*/
    x16[i] = x0[i]*(double)qm[i]*OD_QM_SCALE_1;
    r16[i] = r0[i]*(double)qm[i]*OD_QM_SCALE_1;
#else
    x16[i] = OD_SHR_ROUND(x0[i]*qm[i], OD_QM_SHIFT + xshift);
    r16[i] = OD_SHR_ROUND(r0[i]*qm[i], OD_QM_SHIFT + rshift);
#endif
    corr += OD_MULT16_16(x16[i], r16[i]);
  }
  cfl_enabled = is_keyframe && pli != 0 && !OD_DISABLE_CFL;
  cg  = od_pvq_compute_gain(x16, n, q0, &g, beta, xshift);
  cgr = od_pvq_compute_gain(r16, n, q0, &gr, beta, rshift);
  if (cfl_enabled) cgr = OD_CGAIN_SCALE;
  /* gain_offset is meant to make sure one of the quantized gains has
     exactly the same gain as the reference. */
#if defined(OD_FLOAT_PVQ)
  icgr = (int)floor(.5 + cgr);
#else
  icgr = OD_SHR_ROUND(cgr, OD_CGAIN_SHIFT);
#endif
  gain_offset = cgr - OD_SHL(icgr, OD_CGAIN_SHIFT);
  /* Start search with null case: gain=0, no pulse. */
  qg = 0;
  dist = gain_weight*cg*cg*OD_CGAIN_SCALE_2;
  best_dist = dist;
  best_cost = dist + pvq_norm_lambda*od_pvq_rate(0, 0, -1, 0, adapt, NULL, 0,
   n, is_keyframe, pli, speed);
  noref = 1;
  best_k = 0;
  *itheta = -1;
  *max_theta = 0;
  OD_CLEAR(y, n);
  best_qtheta = 0;
  m = 0;
  s = 1;
  corr = corr/(1e-100 + g*(double)gr/OD_SHL(1, xshift + rshift));
  corr = OD_MAXF(OD_MINF(corr, 1.), -1.);
  if (is_keyframe) skip_dist = gain_weight*cg*cg*OD_CGAIN_SCALE_2;
  else {
    skip_dist = gain_weight*(cg - cgr)*(cg - cgr)
     + cgr*(double)cg*(2 - 2*corr);
    skip_dist *= OD_CGAIN_SCALE_2;
  }
  if (!is_keyframe) {
    /* noref, gain=0 isn't allowed, but skip is allowed. */
    od_val32 scgr;
    scgr = OD_MAXF(0,gain_offset);
    if (icgr == 0) {
      best_dist = gain_weight*(cg - scgr)*(cg - scgr)
       + scgr*(double)cg*(2 - 2*corr);
      best_dist *= OD_CGAIN_SCALE_2;
    }
    best_cost = best_dist + pvq_norm_lambda*od_pvq_rate(0, icgr, 0, 0, adapt,
     NULL, 0, n, is_keyframe, pli, speed);
    best_qtheta = 0;
    *itheta = 0;
    *max_theta = 0;
    noref = 0;
  }
  dist0 = best_dist;
  if (n <= OD_MAX_PVQ_SIZE && !od_vector_is_null(r0, n) && corr > 0) {
    od_val16 xr[MAXN];
    int gain_bound;
    int prev_k;
    pvq_search_item items[20];
    int idx;
    int nitems;
    double cos_dist;
    idx = 0;
    gain_bound = OD_SHR(cg - gain_offset, OD_CGAIN_SHIFT);
    /* Perform theta search only if prediction is useful. */
    theta = OD_ROUND32(OD_THETA_SCALE*acos(corr));
    m = od_compute_householder(r16, n, gr, &s, rshift);
    od_apply_householder(xr, x16, r16, n);
    prev_k = 0;
    for (i = m; i < n - 1; i++) xr[i] = xr[i + 1];
    /* Compute all candidate PVQ searches within a reasonable range of gain
       and theta. */
    for (i = OD_MAXI(1, gain_bound - 1); i <= gain_bound + 1; i++) {
      int j;
      od_val32 qcg;
      int ts;
      int theta_lower;
      int theta_upper;
      /* Quantized companded gain */
      qcg = OD_SHL(i, OD_CGAIN_SHIFT) + gain_offset;
      /* Set angular resolution (in ra) to match the encoded gain */
      ts = od_pvq_compute_max_theta(qcg, beta);
      theta_lower = OD_MAXI(0, (int)floor(.5 +
       theta*OD_THETA_SCALE_1*2/M_PI*ts) - 2);
      theta_upper = OD_MINI(ts - 1, (int)ceil(theta*OD_THETA_SCALE_1*2/M_PI*ts));
      /* Include the angles within a reasonable range. */
      for (j = theta_lower; j <= theta_upper; j++) {
        od_val32 qtheta;
        qtheta = od_pvq_compute_theta(j, ts);
        k = od_pvq_compute_k(qcg, j, qtheta, 0, n, beta, robust || is_keyframe);
        items[idx].gain = i;
        items[idx].theta = j;
        items[idx].k = k;
        items[idx].qcg = qcg;
        items[idx].qtheta = qtheta;
        items[idx].ts = ts;
        idx++;
      }
    }
    nitems = idx;
    cos_dist = 0;
    /* Sort PVQ search candidates in ascending order of pulses K so that
       we can reuse all the previously searched pulses across searches. */
    qsort(items, nitems, sizeof(items[0]),
     (int (*)(const void *, const void *))items_compare);
    /* Search for the best gain/theta in order. */
    for (idx = 0; idx < nitems; idx++) {
      int j;
      od_val32 qcg;
      int ts;
      double cost;
      double dist_theta;
      double sin_prod;
      od_val32 qtheta;
      /* Quantized companded gain */
      qcg = items[idx].qcg;
      i = items[idx].gain;
      j = items[idx].theta;
      /* Set angular resolution (in ra) to match the encoded gain */
      ts = items[idx].ts;
      /* Search for the best angle within a reasonable range. */
      qtheta = items[idx].qtheta;
      k = items[idx].k;
      /* Compute the minimal possible distortion by not taking the PVQ
         cos_dist into account. */
      dist_theta = 2 - 2.*od_pvq_cos(theta - qtheta)*OD_TRIG_SCALE_1;
      dist = gain_weight*(qcg - cg)*(qcg - cg) + qcg*(double)cg*dist_theta;
      dist *= OD_CGAIN_SCALE_2;
      /* If we have no hope of beating skip (including a 1-bit worst-case
         penalty), stop now. */
      if (dist > dist0 + 1.0*pvq_norm_lambda && k != 0) continue;
      sin_prod = od_pvq_sin(theta)*OD_TRIG_SCALE_1*od_pvq_sin(qtheta)*
       OD_TRIG_SCALE_1;
      /* PVQ search, using a gain of qcg*cg*sin(theta)*sin(qtheta) since
         that's the factor by which cos_dist is multiplied to get the
         distortion metric. */
      if (k == 0) {
        cos_dist = 0;
        OD_CLEAR(y_tmp, n-1);
      }
      else if (k != prev_k) {
        cos_dist = pvq_search_rdo_double(xr, n - 1, k, y_tmp,
         qcg*(double)cg*sin_prod*OD_CGAIN_SCALE_2, pvq_norm_lambda, prev_k);
      }
      prev_k = k;
      /* See Jmspeex' Journal of Dubious Theoretical Results. */
      dist_theta = 2 - 2.*od_pvq_cos(theta - qtheta)*OD_TRIG_SCALE_1
       + sin_prod*(2 - 2*cos_dist);
      dist = gain_weight*(qcg - cg)*(qcg - cg) + qcg*(double)cg*dist_theta;
      dist *= OD_CGAIN_SCALE_2;
      /* Do approximate RDO. */
      cost = dist + pvq_norm_lambda*od_pvq_rate(i, icgr, j, ts, adapt, y_tmp,
       k, n, is_keyframe, pli, speed);
      if (cost < best_cost) {
        best_cost = cost;
        best_dist = dist;
        qg = i;
        best_k = k;
        best_qtheta = qtheta;
        *itheta = j;
        *max_theta = ts;
        noref = 0;
        OD_COPY(y, y_tmp, n - 1);
      }
    }
  }
  /* Don't bother with no-reference version if there's a reasonable
     correlation. The only exception is luma on a keyframe because
     H/V prediction is unreliable. */
  if (n <= OD_MAX_PVQ_SIZE &&
   ((is_keyframe && pli == 0) || corr < .5
   || cg < (od_val32)(OD_SHL(2, OD_CGAIN_SHIFT)))) {
    int gain_bound;
    int prev_k;
    gain_bound = OD_SHR(cg, OD_CGAIN_SHIFT);
    prev_k = 0;
    /* Search for the best gain (haven't determined reasonable range yet). */
    for (i = OD_MAXI(1, gain_bound); i <= gain_bound + 1; i++) {
      double cos_dist;
      double cost;
      od_val32 qcg;
      qcg = OD_SHL(i, OD_CGAIN_SHIFT);
      k = od_pvq_compute_k(qcg, -1, -1, 1, n, beta, robust || is_keyframe);
      /* Compute the minimal possible distortion by not taking the PVQ
         cos_dist into account. */
      dist = gain_weight*(qcg - cg)*(qcg - cg);
      dist *= OD_CGAIN_SCALE_2;
      if (dist > dist0 && k != 0) continue;
      cos_dist = pvq_search_rdo_double(x16, n, k, y_tmp,
       qcg*(double)cg*OD_CGAIN_SCALE_2, pvq_norm_lambda, prev_k);
      prev_k = k;
      /* See Jmspeex' Journal of Dubious Theoretical Results. */
      dist = gain_weight*(qcg - cg)*(qcg - cg)
       + qcg*(double)cg*(2 - 2*cos_dist);
      dist *= OD_CGAIN_SCALE_2;
      /* Do approximate RDO. */
      cost = dist + pvq_norm_lambda*od_pvq_rate(i, 0, -1, 0, adapt, y_tmp, k,
       n, is_keyframe, pli, speed);
      if (cost <= best_cost) {
        best_cost = cost;
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
  skip = 0;
  if (noref) {
    if (qg == 0) skip = OD_PVQ_SKIP_ZERO;
  }
  else {
    if (!is_keyframe && qg == 0) {
      skip = (icgr ? OD_PVQ_SKIP_ZERO : OD_PVQ_SKIP_COPY);
    }
    if (qg == icgr && *itheta == 0 && !cfl_enabled) skip = OD_PVQ_SKIP_COPY;
  }
  /* Synthesize like the decoder would. */
  if (skip) {
    if (skip == OD_PVQ_SKIP_COPY) OD_COPY(out, r0, n);
    else OD_CLEAR(out, n);
  }
  else {
    if (noref) gain_offset = 0;
    g = od_gain_expand(OD_SHL(qg, OD_CGAIN_SHIFT) + gain_offset, q0, beta);
    od_pvq_synthesis_partial(out, y, r16, n, noref, g, theta, m, s,
     qm_inv);
  }
  *vk = k;
  *skip_diff += skip_dist - best_dist;
  /* Encode gain differently depending on whether we use prediction or not.
     Special encoding on inter frames where qg=0 is allowed for noref=0
     but not noref=1.*/
  if (is_keyframe) return noref ? qg : neg_interleave(qg, icgr);
  else return noref ? qg - 1 : neg_interleave(qg + 1, icgr + 1);
}

/** Encodes a single vector of integers (eg, a partition within a
 *  coefficient block) using PVQ
 *
 * @param [in,out] ec         range encoder
 * @param [in]     qg         quantized gain
 * @param [in]     theta      quantized post-prediction theta
 * @param [in]     max_theta  maximum possible quantized theta value
 * @param [in]     in         coefficient vector to code
 * @param [in]     n          number of coefficients in partition
 * @param [in]     k          number of pulses in partition
 * @param [in,out] model      entropy encoder state
 * @param [in,out] adapt      adaptation context
 * @param [in,out] exg        ExQ16 expectation of gain value
 * @param [in,out] ext        ExQ16 expectation of theta value
 * @param [in]     nodesync   do not use info that depend on the reference
 * @param [in]     cdf_ctx    selects which cdf context to use
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [in]     code_skip  whether the "skip rest" flag is allowed
 * @param [in]     skip_rest  when set, we skip all higher bands
 * @param [in]     encode_flip whether we need to encode the CfL flip flag now
 * @param [in]     flip       value of the CfL flip flag
 */
static void pvq_encode_partition(od_ec_enc *ec,
                                 int qg,
                                 int theta,
                                 int max_theta,
                                 const od_coeff *in,
                                 int n,
                                 int k,
                                 generic_encoder model[3],
                                 od_adapt_ctx *adapt,
                                 int *exg,
                                 int *ext,
                                 int nodesync,
                                 int cdf_ctx,
                                 int is_keyframe,
                                 int code_skip,
                                 int skip_rest,
                                 int encode_flip,
                                 int flip) {
  int noref;
  int id;
  noref = (theta == -1);
  id = (qg > 0) + 2*OD_MINI(theta + 1,3) + 8*code_skip*skip_rest;
  if (is_keyframe) {
    OD_ASSERT(id != 8);
    if (id >= 8) id--;
  }
  else {
    OD_ASSERT(id != 10);
    if (id >= 10) id--;
  }
  /* Jointly code gain, theta and noref for small values. Then we handle
     larger gain and theta values. For noref, theta = -1. */
  od_encode_cdf_adapt(ec, id, &adapt->pvq.pvq_gaintheta_cdf[cdf_ctx][0],
   8 + 7*code_skip, adapt->pvq.pvq_gaintheta_increment);
  if (encode_flip) {
    /* We could eventually do some smarter entropy coding here, but it would
       have to be good enough to overcome the overhead of the entropy coder.
       An early attempt using a "toogle" flag with simple adaptation wasn't
       worth the trouble. */
    od_ec_enc_bits(ec, flip, 1);
  }
  if (qg > 0) {
    int tmp;
    tmp = *exg;
    generic_encode(ec, &model[!noref], qg - 1, -1, &tmp, 2);
    OD_IIR_DIADIC(*exg, qg << 16, 2);
  }
  if (theta > 1 && (nodesync || max_theta > 3)) {
    int tmp;
    tmp = *ext;
    generic_encode(ec, &model[2], theta - 2, nodesync ? -1 : max_theta - 3,
     &tmp, 2);
    OD_IIR_DIADIC(*ext, theta << 16, 2);
  }
  od_encode_pvq_codeword(ec, &adapt->pvq.pvq_codeword_ctx, in,
   n - (theta != -1), k);
}

/** Quantizes a scalar with rate-distortion optimization (RDO)
 * @param [in] x      unquantized value
 * @param [in] q      quantization step size
 * @param [in] delta0 rate increase for encoding a 1 instead of a 0
 * @param [in] pvq_norm_lambda enc->pvq_norm_lambda for quantized RDO
 * @retval quantized value
 */
int od_rdo_quant(od_coeff x, int q, double delta0, double pvq_norm_lambda) {
  int threshold;
  /* Optimal quantization threshold is 1/2 + lambda*delta_rate/2. See
     Jmspeex' Journal of Dubious Theoretical Results for details. */
  threshold = 128 + OD_CLAMPI(0, (int)(256*pvq_norm_lambda*delta0/2), 128);
  if (abs(x) < q*threshold/256) {
    return 0;
  }
  else {
    return OD_DIV_R0(x, q);
  }
}

#if OD_SIGNAL_Q_SCALING
void od_encode_quantizer_scaling(daala_enc_ctx *enc, int q_scaling,
 int sbx, int sby, int skip) {
  int nhsb;
  OD_ASSERT(skip == !!skip);
  nhsb = enc->state.nhsb;
  OD_ASSERT(sbx < nhsb);
  OD_ASSERT(sby < enc->state.nvsb);
  OD_ASSERT(!skip || q_scaling == 0);
  enc->state.sb_q_scaling[sby*nhsb + sbx] = q_scaling;
  if (!skip) {
    int above;
    int left;
    /* use value from neighbour if possible, otherwise use 0 */
    above = sby > 0 ? enc->state.sb_q_scaling[(sby - 1)*enc->state.nhsb + sbx]
     : 0;
    left = sbx > 0 ? enc->state.sb_q_scaling[sby*enc->state.nhsb + (sbx - 1)]
     : 0;
    od_encode_cdf_adapt(&enc->ec, q_scaling,
     enc->state.adapt.q_cdf[above + left*4], 4,
     enc->state.adapt.q_increment);
  }
}
#endif

/** Encode a coefficient block (excepting DC) using PVQ
 *
 * @param [in,out] enc     daala encoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [in]     in      coefficient block to quantize and encode
 * @param [out]    out     quantized coefficient block
 * @param [in]     q0      scale/quantizer
 * @param [in]     pli     plane index
 * @param [in]     bs      log of the block size minus two
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     robust  make stream robust to error in the reference
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [in]     q_scaling scaling factor to apply to quantizer
 * @param [in]     bx      x-coordinate of this block
 * @param [in]     by      y-coordinate of this block
 * @param [in]     qm      QM with magnitude compensation
 * @param [in]     qm_inv  Inverse of QM with magnitude compensation
 * @return         Returns 1 if both DC and AC coefficients are skipped,
 *                 zero otherwise
 */
int od_pvq_encode(daala_enc_ctx *enc,
                   od_coeff *ref,
                   const od_coeff *in,
                   od_coeff *out,
                   int q0,
                   int pli,
                   int bs,
                   const double *beta,
                   int robust,
                   int is_keyframe,
                   int q_scaling,
                   int bx,
                   int by,
                   const int16_t *qm,
                   const int16_t *qm_inv,
                   int speed){
  int theta[PVQ_MAX_PARTITIONS];
  int max_theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  od_coeff y[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int *exg;
  int *ext;
  int nb_bands;
  int i;
  const int *off;
  int size[PVQ_MAX_PARTITIONS];
  generic_encoder *model;
  double skip_diff;
  int tell;
  uint16_t *skip_cdf;
  od_rollback_buffer buf;
  int dc_quant;
  int flip;
  int cfl_encoded;
  int skip_rest;
  int skip_dir;
  int skip_theta_value;
  const unsigned char *pvq_qm;
  double dc_rate;
#if !OD_SIGNAL_Q_SCALING
  OD_UNUSED(q_scaling);
  OD_UNUSED(bx);
  OD_UNUSED(by);
#endif
  pvq_qm = &enc->state.pvq_qm_q4[pli][0];
  exg = &enc->state.adapt.pvq.pvq_exg[pli][bs][0];
  ext = enc->state.adapt.pvq.pvq_ext + bs*PVQ_MAX_PARTITIONS;
  skip_cdf = enc->state.adapt.skip_cdf[2*bs + (pli != 0)];
  model = enc->state.adapt.pvq.pvq_param_model;
  nb_bands = OD_BAND_OFFSETS[bs][0];
  off = &OD_BAND_OFFSETS[bs][1];
  dc_quant = OD_MAXI(1, q0*pvq_qm[od_qm_get_index(bs, 0)] >> 4);
  tell = 0;
  for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
  skip_diff = 0;
  flip = 0;
  /*If we are coding a chroma block of a keyframe, we are doing CfL.*/
  if (pli != 0 && is_keyframe) {
    od_val32 xy;
    xy = 0;
    /*Compute the dot-product of the first band of chroma with the luma ref.*/
    for (i = off[0]; i < off[1]; i++) {
#if defined(OD_FLOAT_PVQ)
      xy += ref[i]*(double)qm[i]*OD_QM_SCALE_1*
       (double)in[i]*(double)qm[i]*OD_QM_SCALE_1;
#else
      od_val32 rq;
      od_val32 inq;
      rq = ref[i]*qm[i];
      inq = in[i]*qm[i];
      xy += OD_SHR(rq*(int64_t)inq, OD_SHL(OD_QM_SHIFT + OD_CFL_FLIP_SHIFT,
       1));
#endif
    }
    /*If cos(theta) < 0, then |theta| > pi/2 and we should negate the ref.*/
    if (xy < 0) {
      flip = 1;
      for(i = off[0]; i < off[nb_bands]; i++) ref[i] = -ref[i];
    }
  }
  for (i = 0; i < nb_bands; i++) {
    int q;
    q = OD_MAXI(1, q0*pvq_qm[od_qm_get_index(bs, i + 1)] >> 4);
    qg[i] = pvq_theta(out + off[i], in + off[i], ref + off[i], size[i],
     q, y + off[i], &theta[i], &max_theta[i],
     &k[i], beta[i], &skip_diff, robust, is_keyframe, pli, &enc->state.adapt,
     qm + off[i], qm_inv + off[i], enc->pvq_norm_lambda, speed);
  }
  od_encode_checkpoint(enc, &buf);
  if (is_keyframe) out[0] = 0;
  else {
    dc_rate = -OD_LOG2((double)(skip_cdf[3] - skip_cdf[2])/
     (double)(skip_cdf[2] - skip_cdf[1]));
    out[0] = od_rdo_quant(in[0] - ref[0], dc_quant, dc_rate,
     enc->pvq_norm_lambda);
  }
  tell = od_ec_enc_tell_frac(&enc->ec);
  /* Code as if we're not skipping. */
  od_encode_cdf_adapt(&enc->ec, 2 + (out[0] != 0), skip_cdf,
   4 + (pli == 0 && bs > 0), enc->state.adapt.skip_increment);
#if OD_SIGNAL_Q_SCALING
  if (bs == OD_NBSIZES - 1 && pli == 0) {
    od_encode_quantizer_scaling(enc, q_scaling, bx >> (OD_NBSIZES - 1),
     by >> (OD_NBSIZES - 1), 0);
  }
#endif
  cfl_encoded = 0;
  skip_rest = 1;
  skip_theta_value = is_keyframe ? -1 : 0;
  for (i = 1; i < nb_bands; i++) {
    if (theta[i] != skip_theta_value || qg[i]) skip_rest = 0;
  }
  skip_dir = 0;
  if (nb_bands > 1) {
    for (i = 0; i < 3; i++) {
      int j;
      int tmp;
      tmp = 1;
      for (j = i + 1; j < nb_bands; j += 3) {
        if (theta[j] != skip_theta_value || qg[j]) tmp = 0;
      }
      skip_dir |= tmp << i;
    }
  }
  if (theta[0] == skip_theta_value && qg[0] == 0 && skip_rest) nb_bands = 0;
  for (i = 0; i < nb_bands; i++) {
    int encode_flip;
    /* Encode CFL flip bit just after the first time it's used. */
    encode_flip = pli != 0 && is_keyframe && theta[i] != -1 && !cfl_encoded;
    if (i == 0 || (!skip_rest && !(skip_dir & (1 << ((i - 1)%3))))) {
      pvq_encode_partition(&enc->ec, qg[i], theta[i], max_theta[i], y + off[i],
       size[i], k[i], model, &enc->state.adapt, exg + i, ext + i,
       robust || is_keyframe, (pli != 0)*OD_NBSIZES*PVQ_MAX_PARTITIONS
       + bs*PVQ_MAX_PARTITIONS + i, is_keyframe, i == 0 && (i < nb_bands - 1),
       skip_rest, encode_flip, flip);
    }
    if (i == 0 && !skip_rest && bs > 0) {
      od_encode_cdf_adapt(&enc->ec, skip_dir,
       &enc->state.adapt.pvq.pvq_skip_dir_cdf[(pli != 0) + 2*(bs - 1)][0], 7,
       enc->state.adapt.pvq.pvq_skip_dir_increment);
    }
    if (encode_flip) cfl_encoded = 1;
  }
  tell = od_ec_enc_tell_frac(&enc->ec) - tell;
  /* Account for the rate of skipping the AC, based on the same DC decision
     we made when trying to not skip AC. */
  {
    double skip_rate;
    if (out[0] != 0) {
      skip_rate = -OD_LOG2((skip_cdf[1] - skip_cdf[0])/
     (double)skip_cdf[3 + (pli == 0 && bs > 0)]);
    }
    else {
      skip_rate = -OD_LOG2(skip_cdf[0]/
     (double)skip_cdf[3 + (pli == 0 && bs > 0)]);
    }
    tell -= (int)floor(.5+8*skip_rate);
  }
  if (nb_bands == 0 || skip_diff <= enc->pvq_norm_lambda/8*tell) {
    if (is_keyframe) out[0] = 0;
    else {
      dc_rate = -OD_LOG2((double)(skip_cdf[1] - skip_cdf[0])/
       (double)skip_cdf[0]);
      out[0] = od_rdo_quant(in[0] - ref[0], dc_quant, dc_rate,
       enc->pvq_norm_lambda);
    }
    /* We decide to skip, roll back everything as it was before. */
    od_encode_rollback(enc, &buf);
    od_encode_cdf_adapt(&enc->ec, out[0] != 0, skip_cdf,
     4 + (pli == 0 && bs > 0), enc->state.adapt.skip_increment);
#if OD_SIGNAL_Q_SCALING
    if (bs == OD_NBSIZES - 1 && pli == 0) {
      int skip;
      skip = out[0] == 0;
      if (skip) {
        q_scaling = 0;
      }
      od_encode_quantizer_scaling(enc, q_scaling, bx >> (OD_NBSIZES - 1),
       by >> (OD_NBSIZES - 1), skip);
    }
#endif
    if (is_keyframe) for (i = 1; i < 1 << (2*bs + 4); i++) out[i] = 0;
    else for (i = 1; i < 1 << (2*bs + 4); i++) out[i] = ref[i];
    if (out[0] == 0) return 1;
  }
  return 0;
}
