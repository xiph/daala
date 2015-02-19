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
#include "logging.h"
#include "entenc.h"
#include "entcode.h"
#include "laplace_code.h"
#include "pvq_code.h"
#include "adapt.h"
#include "filter.h"

#define OD_PVQ_RATE_APPROX (0)

static void od_encode_pvq_codeword(od_ec_enc *ec, od_adapt_ctx *adapt,
 const od_coeff *in, int n, int k, int noref, int ln) {
  if (k == 1 && n < 16) {
    int cdf_id;
    int i;
    int pos;
    cdf_id = 2*(n == 15) + !noref;
    pos = 32;
    for (i = 0; i < n - !noref; i++) {
      if (in[i]) {
        pos = i;
        break;
      }
    }
    OD_ASSERT(pos < n - !noref);
    od_encode_cdf_adapt(ec, pos, adapt->pvq_k1_cdf[cdf_id], n - !noref,
     adapt->pvq_k1_increment);
    od_ec_enc_bits(ec, in[pos] < 0, 1);
  }
  else {
    int speed = 5;
    int *pvq_adapt;
    int adapt_curr[OD_NSB_ADAPT_CTXS] = { 0 };
    pvq_adapt = adapt->pvq_adapt + 4*(2*ln + noref);
    laplace_encode_vector(ec, in, n - !noref, k, adapt_curr,
     pvq_adapt);
    if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
      pvq_adapt[OD_ADAPT_K_Q8] += (256*adapt_curr[OD_ADAPT_K_Q8]
       - pvq_adapt[OD_ADAPT_K_Q8]) >> speed;
      pvq_adapt[OD_ADAPT_SUM_EX_Q8] += (adapt_curr[OD_ADAPT_SUM_EX_Q8]
       - pvq_adapt[OD_ADAPT_SUM_EX_Q8]) >> speed;
    }
    if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
      pvq_adapt[OD_ADAPT_COUNT_Q8] += (adapt_curr[OD_ADAPT_COUNT_Q8]
       - pvq_adapt[OD_ADAPT_COUNT_Q8]) >> speed;
      pvq_adapt[OD_ADAPT_COUNT_EX_Q8] += (adapt_curr[OD_ADAPT_COUNT_EX_Q8]
       - pvq_adapt[OD_ADAPT_COUNT_EX_Q8]) >> speed;
    }
  }
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
 * @return                  cosine distance between x and y (between 0 and 1)
 */
static double pvq_search_rdo_double(const double *xcoeff, int n, int k,
 od_coeff *ypulse, double g2) {
  int i, j;
  double xy;
  double yy;
  double x[1024];
  double xx;
  double lambda;
  double norm_1;
  int rdo_pulses;
  double delta_rate;
  xx = xy = yy = 0;
  for (j = 0; j < n; j++) {
    x[j] = fabs(xcoeff[j]);
    xx += x[j]*x[j];
  }
  norm_1 = 1./sqrt(1e-30 + xx);
  lambda = OD_PVQ_LAMBDA/(1e-30 + g2);
  i = 0;
  if (k > 2) {
    double l1_norm;
    double l1_inv;
    l1_norm = 0;
    for (j = 0; j < n; j++) l1_norm += x[j];
    l1_inv = 1./OD_MAXF(l1_norm, 1e-100);
    for (j = 0; j < n; j++) {
      ypulse[j] = OD_MAXI(0, (int)floor(k*x[j]*l1_inv));
      xy += x[j]*ypulse[j];
      yy += ypulse[j]*ypulse[j];
      i += ypulse[j];
    }
  }
  else {
    for (j = 0; j < n; j++) ypulse[j] = 0;
  }
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
  /* Search last pulses with RDO. Distortion is D = (x-y)^2 = x^2 - x*y + y^2
     and since x^2 and y^2 are constant, we just maximize x*y, plus a
     lambda*rate term. Note that since x and y aren't normalized here,
     we need to divide by sqrt(x^2)*sqrt(y^2). */
  for (; i < k; i++) {
    int pos;
    double best_cost;
    pos = 0;
    best_cost = -1e5;
    for (j = 0; j < n; j++) {
      double tmp_xy;
      double tmp_yy;
      tmp_xy = xy + x[j];
      tmp_yy = yy + 2*ypulse[j] + 1;
      tmp_xy = 2*tmp_xy*norm_1*od_rsqrt_table(tmp_yy) - lambda*j*delta_rate;
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
 int is_keyframe, int pli, int ln) {
  double rate;
#if OD_PVQ_RATE_APPROX
  /* Estimates the number of bits it will cost to encode K pulses in
     N dimensions based on experimental data for bitrate vs K. */
  rate = n*OD_LOG2(1+log(n*2)*k/n);
  (void)adapt;
  (void)m;
  (void)y0;
#else
  if (k > 0){
    od_ec_enc ec;
    od_adapt_ctx ad;
    int tell;
    od_ec_enc_init(&ec, 1000);
    OD_COPY(&ad, adapt, 1);
    tell = od_ec_enc_tell_frac(&ec);
    od_encode_pvq_codeword(&ec, &ad, y0, n, k, theta == -1, ln);
    rate = (od_ec_enc_tell_frac(&ec)-tell)/8.;
    od_ec_enc_clear(&ec);
  }
  else rate = 0;
#endif
  if (qg > 0 && theta >= 0) {
    /* Approximate cost of entropy-coding theta */
    rate += .9*OD_LOG2(ts);
    /* Adding a cost to using the H/V pred because it's going to be off
       most of the time. Cost is optimized on subset1. */
    if (is_keyframe && pli == 0) rate += 2.5;
    if (qg == icgr) rate -= .5;
  }
  return rate;
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
 * @param [in]     ln        log of the block size minus two
 * @return         gain      index of the quatized gain
*/
static int pvq_theta(od_coeff *out, od_coeff *x0, od_coeff *r0, int n, int q0,
 od_coeff *y, int *itheta, int *max_theta, int *vk,
 double beta, double *skip_diff, int robust, int is_keyframe, int pli,
 const od_adapt_ctx *adapt, int ln) {
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
  /* Best RDO cost (D + lamdba*R) so far. */
  double best_cost;
  /* Distortion (D) that corresponds to the best RDO cost. */
  double best_dist;
  double dist;
  /* Sign of Householder reflection. */
  int s;
  /* Dimension on which Householder reflects. */
  int m;
  double theta;
  double corr;
  int best_k;
  double best_qtheta;
  double gain_offset;
  int noref;
  double lambda;
  double skip_dist;
  int cfl_enabled;
  int skip;
  double gain_weight;
  lambda = OD_PVQ_LAMBDA;
  /* Give more weight to gain error when calculating the total distortion. */
  gain_weight = 1.4;
  OD_ASSERT(n > 1);
  corr = 0;
  for (i = 0; i < n; i++) {
    x[i] = x0[i];
    r[i] = r0[i];
    corr += x[i]*r[i];
  }
  cfl_enabled = is_keyframe && pli != 0 && !OD_DISABLE_CFL;
  cg  = od_pvq_compute_gain(x0, n, q0, &g, beta);
  cgr = od_pvq_compute_gain(r0, n, q0, &gr, beta);
  if (pli != 0 && is_keyframe && !OD_DISABLE_CFL) cgr = 1;
  /* gain_offset is meant to make sure one of the quantized gains has
     exactly the same gain as the reference. */
  icgr = (int)floor(.5+cgr);
  gain_offset = cgr-icgr;
  /* Start search with null case: gain=0, no pulse. */
  qg = 0;
  dist = gain_weight*cg*cg;
  best_dist = dist;
  best_cost = dist + lambda*od_pvq_rate(0, 0, -1, 0, adapt, NULL, 0, n,
   is_keyframe, pli, ln);
  noref = 1;
  best_k = 0;
  *itheta = -1;
  *max_theta = 0;
  OD_CLEAR(y, n);
  best_qtheta = 0;
  m = 0;
  s = 1;
  corr = corr/(1e-100 + g*gr);
  corr = OD_MAXF(OD_MINF(corr, 1.), -1.);
  skip_dist = gain_weight*(cg - cgr)*(cg - cgr) + cgr*cg*(2 - 2*corr);
  if (!is_keyframe) {
    /* noref, gain=0 isn't allowed, but skip is allowed. */
    double scgr;
    scgr = OD_MAXF(0,gain_offset);
    if (icgr == 0) {
      best_dist = gain_weight*(cg - scgr)*(cg - scgr) + scgr*cg*(2 - 2*corr);
    }
    best_cost = best_dist + lambda*od_pvq_rate(0, icgr, 0, 0, adapt, NULL,
     0, n, is_keyframe, pli, ln);
    best_qtheta = 0;
    *itheta = 0;
    *max_theta = 0;
    noref = 0;
  }
  if (!od_vector_is_null(r0, n) && corr > 0) {
    /* Perform theta search only if prediction is useful. */
    theta = acos(corr);
    m = od_compute_householder(r, n, gr, &s);
    od_apply_householder(x, r, n);
    for (i = m; i < n - 1; i++) x[i] = x[i + 1];
    /* Search for the best gain within a reasonable range. */
    for (i = OD_MAXI(1, (int)floor(cg-gain_offset));
     i <= (int)ceil(cg-gain_offset); i++) {
      int j;
      double qcg;
      int ts;
      /* Quantized companded gain */
      qcg = i+gain_offset;
      /* Set angular resolution (in ra) to match the encoded gain */
      ts = od_pvq_compute_max_theta(qcg, beta);
      /* Search for the best angle within a reasonable range. */
      for (j = OD_MAXI(0, (int)floor(.5+theta*2/M_PI*ts)-1);
       j <= OD_MINI(ts-1, (int)ceil(theta*2/M_PI*ts)); j++) {
        double cos_dist;
        double cost;
        double dist_theta;
        double qtheta = od_pvq_compute_theta(j, ts);
        k = od_pvq_compute_k(qcg, j, qtheta, 0, n, beta, robust || is_keyframe);
        /* PVQ search, using a gain of qcg*cg*sin(theta)*sin(qtheta) since
           that's the factor by which cos_dist is multiplied to get the
           distortion metric. */
        cos_dist = pvq_search_rdo_double(x, n - 1, k, y_tmp,
         qcg*cg*sin(theta)*sin(qtheta));
        /* See Jmspeex' Journal of Dubious Theoretical Results. */
        dist_theta = 2 - 2*cos(theta - qtheta)
         + sin(theta)*sin(qtheta)*(2 - 2*cos_dist);
        dist = gain_weight*(qcg - cg)*(qcg - cg) + qcg*cg*dist_theta;
        /* Do approximate RDO. */
        cost = dist + lambda*od_pvq_rate(i, icgr, j, ts, adapt, y_tmp, k, n,
         is_keyframe, pli, ln);
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
  }
  /* Don't bother with no-reference version if there's a reasonable
     correlation */
  if (corr < .5 || cg < 2.) {
    double x1[MAXN];
    for (i = 0; i < n; i++) x1[i] = x0[i];
    /* Search for the best gain (haven't determined reasonable range yet). */
    for (i = OD_MAXI(1, (int)floor(cg)); i <= ceil(cg); i++) {
      double cos_dist;
      double cost;
      double qcg;
      qcg = i;
      k = od_pvq_compute_k(qcg, -1, -1, 1, n, beta, robust || is_keyframe);
      cos_dist = pvq_search_rdo_double(x1, n, k, y_tmp, qcg*cg);
      /* See Jmspeex' Journal of Dubious Theoretical Results. */
      dist = gain_weight*(qcg - cg)*(qcg - cg) + qcg*cg*(2 - 2*cos_dist);
      /* Do approximate RDO. */
      cost = dist + lambda*od_pvq_rate(i, 0, -1, 0, adapt, y_tmp, k, n,
       is_keyframe, pli, ln);
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
    g = od_gain_expand(qg + gain_offset, q0, beta);
    od_pvq_synthesis_partial(out, y, r, n, noref, g, theta, m, s);
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
 * @param [in,out] model      entropy encoder state
 * @param [in,out] adapt      adaptation context
 * @param [in,out] exg        ExQ16 expectation of gain value
 * @param [in,out] ext        ExQ16 expectation of theta value
 * @param [in]     nodesync   do not use info that depend on the reference
 * @param [in]     cdf_ctx    selects which cdf context to use
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @param [in]     code_skip  whether the "skip rest" flag is allowed
 * @param [in]     skip_rest  when set, we skip all higher bands
 * @param [in]     ln         log of the block size minus two
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
                                 int ln) {
  int noref;
  int id;
  noref = (theta == -1);
  id = (qg > 0) + 2*OD_MINI(theta + 1,3) + 8*code_skip*skip_rest;
  if (!is_keyframe) {
    OD_ASSERT(id != 10);
    if (id >= 10) id--;
  }
  /* Jointly code gain, theta and noref for small values. Then we handle
     larger gain and theta values. For noref, theta = -1. */
  od_encode_cdf_adapt(ec, id, &adapt->pvq_gaintheta_cdf[cdf_ctx][0],
   8 + (8 - !is_keyframe)*code_skip, adapt->pvq_gaintheta_increment);
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
  od_encode_pvq_codeword(ec, adapt, in, n, k, theta == -1, ln);
}

/** Quantizes a scalar with rate-distortion optimization (RDO)
 * @param [in] x      unquantized value
 * @param [in] q      quantization step size
 * @param [in] delta0 rate increase for encoding a 1 instead of a 0
 * @retval quantized value
 */
int od_rdo_quant(od_coeff x, int q, double delta0) {
  int threshold;
  /* Optimal quantization threshold is 1/2 + lambda*delta_rate/2. See
     Jmspeex' Journal of Dubious Theoretical Results for details. */
  threshold = 128 + OD_CLAMPI(0, (int)(256*OD_PVQ_LAMBDA*delta0/2), 128);
  if (abs(x) < q * threshold / 256) {
    return 0;
  } else {
    return OD_DIV_R0(x, q);
  }
}

/** Encode a coefficient block (excepting DC) using PVQ
 *
 * @param [in,out] enc     daala encoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [in]     in      coefficient block to quantize and encode
 * @param [out]    out     quantized coefficient block
 * @param [in]     q       scale/quantizer
 * @param [in]     pli     plane index
 * @param [in]     ln      log of the block size minus two
 * @param [in]     qm      per-band quantization matrix
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     robust  make stream robust to error in the reference
 * @param [in]     is_keyframe whether we're encoding a keyframe
 */
int od_pvq_encode(daala_enc_ctx *enc,
                   od_coeff *ref,
                   od_coeff *in,
                   od_coeff *out,
                   int q0,
                   int pli,
                   int ln,
                   const double *beta,
                   int robust,
                   int is_keyframe){
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
  unsigned tell;
  ogg_uint16_t *skip_cdf;
  od_rollback_buffer buf;
  int dc_quant;
  int flip;
  int cfl_encoded;
  int skip_rest;
  int skip_dir;
  int skip_theta_value;
  const unsigned char *qm;
  qm = &enc->state.pvq_qm_q4[pli][0];
  exg = &enc->state.adapt.pvq_exg[pli][ln][0];
  ext = enc->state.adapt.pvq_ext + ln*PVQ_MAX_PARTITIONS;
  skip_cdf = enc->state.adapt.skip_cdf[pli];
  model = enc->state.adapt.pvq_param_model;
  nb_bands = OD_BAND_OFFSETS[ln][0];
  off = &OD_BAND_OFFSETS[ln][1];
  dc_quant = OD_MAXI(1, q0*qm[od_qm_get_index(ln, 0)] >> 4);
  tell = 0;
  for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
  skip_diff = 0;
  flip = 0;
  if (pli != 0 && is_keyframe) {
    double xy;
    xy = 0;
    for (i = 1; i < 16; i++) xy += ref[i]*(double)in[i];
    if (xy < 0) {
      flip = 1;
      for(i = 1; i < off[nb_bands]; i++) ref[i] = -ref[i];
    }
  }
  for (i = 0; i < nb_bands; i++) {
    int q;
    q = OD_MAXI(1, q0*qm[od_qm_get_index(ln, i + 1)] >> 4);
    qg[i] = pvq_theta(out + off[i], in + off[i], ref + off[i], size[i],
     q, y + off[i], &theta[i], &max_theta[i],
     &k[i], beta[i], &skip_diff, robust, is_keyframe, pli, &enc->state.adapt,
     ln);
  }
  if (!is_keyframe) {
    double dc_rate;
    od_encode_checkpoint(enc, &buf);
    dc_rate = -OD_LOG2((double)(skip_cdf[1]-skip_cdf[0])/(double)skip_cdf[0]);
    out[0] = od_rdo_quant(in[0] - ref[0], dc_quant, dc_rate);
    /* Code as if we're not skipping. */
    od_encode_cdf_adapt(&enc->ec, (out[0] != 0), skip_cdf,
     4, enc->state.adapt.skip_increment);
    /* Excluding skip flag from the rate since it's minor and would be prone
       to greedy decision issues. */
    tell = od_ec_enc_tell_frac(&enc->ec);
  }
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
  if (!is_keyframe && theta[0] == 0 && qg[0] == 0 && skip_rest) nb_bands = 0;
  for (i = 0; i < nb_bands; i++) {
    if (i == 0 || (!skip_rest && !(skip_dir & (1 << ((i - 1)%3))))) {
      pvq_encode_partition(&enc->ec, qg[i], theta[i], max_theta[i], y + off[i],
       size[i], k[i], model, &enc->state.adapt, exg + i, ext + i,
       robust || is_keyframe, (pli != 0)*OD_NBSIZES*PVQ_MAX_PARTITIONS
       + ln*PVQ_MAX_PARTITIONS + i, is_keyframe, i == 0 && (i < nb_bands - 1),
       skip_rest, ln);
    }
    if (i == 0 && !skip_rest && ln > 0) {
      od_encode_cdf_adapt(&enc->ec, skip_dir,
       &enc->state.adapt.pvq_skip_dir_cdf[(pli != 0) + 2*(ln - 1)][0], 7,
       enc->state.adapt.pvq_skip_dir_increment);
    }
    /* Encode CFL flip bit just after the first time it's used. */
    if (pli!=0 && is_keyframe && theta[i] != -1 && !cfl_encoded) {
      /* We could eventually do some smarter entropy coding here, but it would
         have to be good enough to overcome the overhead of the entropy coder.
         An early attempt using a "toogle" flag with simple adaptation wasn't
         worth the trouble. */
      od_ec_enc_bits(&enc->ec, flip, 1);
      cfl_encoded = 1;
    }
  }
  if (!is_keyframe) {
    tell = od_ec_enc_tell_frac(&enc->ec) - tell;
    if (nb_bands == 0 || skip_diff <= OD_PVQ_LAMBDA/8*tell) {
      double dc_rate;
      dc_rate = -OD_LOG2((double)(skip_cdf[3]-skip_cdf[2])/
       (double)(skip_cdf[2]-skip_cdf[1]));
      out[0] = od_rdo_quant(in[0] - ref[0], dc_quant, dc_rate);
      /* We decide to skip, roll back everything as it was before. */
      od_encode_rollback(enc, &buf);
      od_encode_cdf_adapt(&enc->ec, 2 + (out[0] != 0), skip_cdf,
       4, enc->state.adapt.skip_increment);
      for (i = 1; i < 1 << (2*ln + 4); i++) out[i] = ref[i];
      if ((out[0] == 0)) return 1;
    }
  }
  return 0;
}
