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

/* This is the PVQ equivalent of a quantization matrix, except that
   the values are per-band. */

#if OD_DISABLE_MASKING

static const int od_pvq_qm4_luma[2] = {8, 16};
static const int od_pvq_qm8_luma[5] = {16, 16, 32, 32, 64};
static const int od_pvq_qm16_luma[8] = {16, 16, 16, 16, 32, 32, 32, 64};

static const double od_pvq_mask4_luma[1] = {1.};
static const double od_pvq_mask8_luma[4] = {1., 1., 1., 1.};
static const double od_pvq_mask16_luma[7] = {1., 1., 1., 1., 1., 1., 1.};

#else

static const int od_pvq_qm4_luma[2] = {8, 16};
static const int od_pvq_qm8_luma[5] = {16, 16, 44, 44, 72};
static const int od_pvq_qm16_luma[8] = {16, 13, 18, 18, 36, 40, 40, 80};

static const double od_pvq_mask4_luma[1] = {1.};
static const double od_pvq_mask8_luma[4] = {1.5, 1.5, 1.5, 1.5};
static const double od_pvq_mask16_luma[7] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

#endif


static const int od_pvq_qm4_chroma[2] = {8, 32};
static const int od_pvq_qm8_chroma[5] = {8, 24, 64, 64, 100};
static const int od_pvq_qm16_chroma[8] = {8, 16, 24, 24, 32, 64, 64, 100};

static const double od_pvq_mask4_chroma[1] = {1.};
static const double od_pvq_mask8_chroma[4] = {1., 1., 1., 1.};
static const double od_pvq_mask16_chroma[7] = {1., 1., 1., 1., 1., 1., 1.};

const int * const od_pvq_qm[3][OD_NBSIZES] = {
  {od_pvq_qm4_luma, od_pvq_qm8_luma, od_pvq_qm16_luma},
  {od_pvq_qm4_chroma, od_pvq_qm8_chroma, od_pvq_qm16_chroma},
  {od_pvq_qm4_chroma, od_pvq_qm8_chroma, od_pvq_qm16_chroma}
};

const double * const od_pvq_mask[3][OD_NBSIZES] = {
  {od_pvq_mask4_luma, od_pvq_mask8_luma, od_pvq_mask16_luma},
  {od_pvq_mask4_chroma, od_pvq_mask8_chroma, od_pvq_mask16_chroma},
  {od_pvq_mask4_chroma, od_pvq_mask8_chroma, od_pvq_mask16_chroma}
};

/** Find the codepoint on the given PSphere closest to the desired
 * vector. Double-precision PVQ search just to make sure our tests
 * aren't limited by numerical accuracy.
 *
 * @param [in]      xcoeff  input vector to quantize (x in the math doc)
 * @param [in]      n       number of dimensions
 * @param [in]      k       number of pulses
 * @param [out]     ypulse  optimal codevector found (y in the math doc)
 * @return                  cosine distance between x and y (between 0 and 1)
 */
static double pvq_search_double(const double *xcoeff, int n, int k,
                                od_coeff *ypulse) {
  int i, j;
  double xy;
  double yy;
  double x[1024];
  double xx;
  xx = xy = yy = 0;
  for (j = 0; j < n; j++) {
    x[j] = fabs(xcoeff[j]);
    xx += x[j]*x[j];
  }
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
  /* Search one pulse at a time */
  for (; i < k; i++) {
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
  for (i = 0; i < n; i++) {
    if (xcoeff[i] < 0) ypulse[i] = -ypulse[i];
  }
  return xy/(1e-100 + sqrt(xx*yy));
}

/** Computes Householder reflection that aligns the reference r to the
 *  dimension in r with the greatest absolute value. The reflection
 *  vector is returned in r.
 *
 * @param [in,out]  r      reference vector to be reflected, reflection
 *                         also returned in r
 * @param [in]      n      number of dimensions in r
 * @param [in]      gr     gain of reference vector
 * @param [out]     sign   sign of reflection
 * @return                 dimension number to which reflection aligns
 **/
static int compute_householder(double *r, int n, double gr, int *sign) {
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

/** Applies Householder reflection from compute_householder(). The
 * reflection is its own inverse.
 *
 * @param [in,out]  x      vector to be reflected
 * @param [in]      r      reflection
 * @param [in]      n      number of dimensions in x,r
 */
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

/** Inverse of neg_interleave; decodes the interleaved gain.
 *
 * @param [in]      x      quantized/interleaved gain to decode
 * @param [in]      ref    quantized gain of the reference
 * @return                 original quantized gain value
 */
int neg_deinterleave(int x, int ref) {
  if (x < 2*ref-1) {
    if (x & 1) return ref - 1 - (x >> 1);
    else return ref + (x >> 1);
  }
  else return x+1;
}

/** Computes the raw and quantized/companded gain of a given input
 * vector
 *
 * @param [in]      x      vector of input data
 * @param [in]      n      number of elements in vector x
 * @param [in]      q      quantizer
 * @param [out]     g      raw gain
 * @return                 quantized/companded gain
 */
double pvq_compute_gain(od_coeff *x, int n, double q, double *g, double mask){
  int i;
  double acc=0;
  for (i = 0; i < n; i++) acc += x[i]*(double)x[i];
  *g = sqrt(acc);
  /* Normalize gain by quantization step size and apply companding
     (if ACTIVITY != 1). */
  return pow(*g, 1./mask)/q;
}

/** Compute theta quantization range from quantized/companded gain
 *
 * @param [in]      qcg    quantized companded gain value
 * @return                 max theta value
 */
int pvq_compute_max_theta(double qcg, double mask){
  /* Set angular resolution (in ra) to match the encoded gain */
  int ts = (int)floor(.5 + qcg*M_PI/(2*mask));
  /* Special case for low gains -- will need to be tuned anyway */
  if (qcg < 1.4) ts = 1;
  return ts;
}

/** Decode quantized theta value from coded value
 *
 * @param [in]      t          quantized companded gain value
 * @param [in]      max_theta  maximum theta value
 * @return                     decoded theta value
 */
double pvq_compute_theta(int t, int max_theta) {
  if (max_theta != 0) return t*.5*M_PI/max_theta;
  return 0;
}

/** Compute the number of pulses used for PVQ encoding a vector from
 * available metrics (encode and decode side)
 *
 * @param [in]      qcg        quantized companded gain value
 * @param [in]      theta      PVQ error angle theta
 * @param [in]      noref      indicates present or lack of reference
 *                             (prediction)
 * @param [in]      n          number of elements to be coded
 * @return                     number of pulses to use for coding
 */
int pvq_compute_k(double qcg, double theta, int noref, int n, double mask) {
  if (noref) {
    return OD_MAXI(1, (int)floor(.5 + qcg*sqrt(n/2)/mask));
  }
  else {
    if (theta == 0) return 0;
    /* Sets K according to gain and theta, based on the high-rate
       PVQ distortion curves D~=N^2/(24*K^2). Low-rate will have to be
       perceptually tuned anyway.  */
    return OD_MAXI(1, (int)floor(.5 + qcg*sin(theta)*sqrt(n/2)/mask));
  }
}

/** Synthesizes one parition of coefficient values from a PVQ-encoded
 * vector.  This 'partial' version is called by the encode loop where
 * the Householder reflection has already been computed and there's no
 * need to recompute it.
 *
 * @param [out]     xcoeff  output coefficient partition (x in math doc)
 * @param [in]      ypulse  PVQ-encoded values (y in the math doc); in
 *                          the noref case, this vector has n entries,
 *                          in the reference case it contains n-1 entries
 *                          (the m-th entry is not included)
 * @param [in]      r       reference vector (prediction)
 * @param [in]      n       number of elements in this partition
 * @param [in]      noref   indicates presence or lack of prediction
 * @param [in]      qg      decoded quantized vector gain
 * @param [in]      go      gain offset for predicted case
 * @param [in]      theta   decoded theta (prediction error)
 * @param [in]      m       alignment dimension of Householder reflection
 * @param [in]      s       sign of Householder reflection
 * @param [in]      q       gain quantizer
 */
static double pvq_synthesis_partial(od_coeff *xcoeff, od_coeff *ypulse,
                                  const double *r, int n,
                                  int noref, int qg, double go,
                                  double theta, int m, int s, double q,
                                  double mask) {
  int i;
  int yy;
  double qcg;
  double norm;
  double g;
  double x[MAXN];
  int nn = n-(!noref); /* when noref==0, vector in is sized n-1 */

  yy = 0;
  for (i = 0; i < nn; i++)
    yy += ypulse[i]*(ogg_int32_t)ypulse[i];
  norm = sqrt(1./(1e-100 + yy));

  if (noref) {
    qcg = qg;
    for (i = 0; i < n; i++)
      x[i] = ypulse[i]*norm;
  }
  else{
    qcg = qg==0 ? 0 : qg+go;
    norm *= sin(theta);
    for (i = 0; i < m; i++)
      x[i] = ypulse[i]*norm;
    x[m] = -s*cos(theta);
    for (i = m; i < nn; i++)
      x[i+1] = ypulse[i]*norm;
    apply_householder(x, r, n);
  }

  g = pow(q*qcg, mask);
  for (i = 0; i < n; i++) {
    xcoeff[i] = floor(.5 + x[i]*g);
  }
  return g;
}

/** Synthesizes one parition of coefficient values from a PVQ-encoded
 * vector.
 *
 * @param [out]     xcoeff  output coefficient partition (x in math doc)
 * @param [in]      ypulse  PVQ-encoded values (y in math doc); in the noref
 *                          case, this vector has n entries, in the
 *                          reference case it contains n-1 entries
 *                          (the m-th entry is not included)
 * @param [in]      r       reference vector (prediction)
 * @param [in]      n       number of elements in this partition
 * @param [in]      gr      gain of the reference vector (prediction)
 * @param [in]      noref   indicates presence or lack of prediction
 * @param [in]      qg      decoded quantized vector gain
 * @param [in]      go      gain offset for predicted case
 * @param [in]      theta   decoded theta (prediction error)
 * @param [in]      m       alignment dimension of Householder reflection
 * @param [in]      s       sign of Householder reflection
 * @param [in]      q       gain quantizer
 */
double pvq_synthesis(od_coeff *xcoeff, od_coeff *ypulse, double *r, int n,
                   double gr, int noref, int qg, double go,
                   double theta, double q, double mask) {
  int s = 0;
  int m = noref ? 0 : compute_householder(r, n, gr, &s);

  return pvq_synthesis_partial(xcoeff, ypulse, r, n, noref, qg,
                        go, theta, m, s, q, mask);
}

/* Estimates the number of bits it will cost to encode K pulses in
   N dimensions. This could be improved by using a table, but in the
   short-term, the approximation is more general. It is based on experimental
   data for bitrate vs K. */
static double pvq_rate_approx(int n, int k)
{
  return n*OD_LOG2(1+log(n*2)*k/n);
}

/** Computes the effect of masking from other bands.
 *
 * @param [in]     inter     Combined gain from other bands
 * @param [in]     curr      Quantized gain from current band
 * @param [in]     mask      Activity masking exponent (beta)
 * @return         ratio     Reduction of K/theta allocation
 */
double pvq_interband_masking(double inter, double curr, double mask) {
  return pow(curr*curr/(curr*curr+inter*inter), .5-.5/mask);
}

/** Perform PVQ quantization with prediction, trying several
 * possible gains and angles. See draft-valin-videocodec-pvq and
 * http://jmvalin.ca/slides/pvq.pdf for more details.
 *
 * @param [in,out] x0        coefficients being quantized (before and after)
 * @param [in]     r0        reference, aka predicted coefficients
 * @param [in]     n         number of dimensions
 * @param [in]     q0        quantization step size
 * @param [out]    y         pulse vector (i.e. selected PVQ codevector)
 * @param [out]    itheta    angle between input and reference (-1 if noref)
 * @param [out]    max_theta maximum value of itheta that could have been
 * @param [out]    vk        total number of pulses
 * @return         gain      index of the quatized gain
*/
int pvq_theta(od_coeff *x0, od_coeff *r0, int n, int q0, od_coeff *y, int *itheta,
 int *max_theta, int *vk, double *mask_gain, double mask) {
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
  int noref;
  double lambda;
  double mask_ratio;
  /* Normalized lambda. Since we normalize the gain by q, the distortion is
     normalized by q^2 and lambda does not need the q^2 factor. At high rate,
     this would be log(2)/6, but we're making RDO a bit less aggressive for
     now. */
  lambda = .025;
  /* Quantization step calibration to account for the activity masking. */
  q = q0*pow(256<<OD_COEFF_SHIFT, 1./mask - 1);
  OD_ASSERT(n > 1);
  corr = 0;
  for (i = 0; i < n; i++) {
    x[i] = x0[i];
    r[i] = r0[i];
    corr += x[i]*r[i];
  }

  cg  = pvq_compute_gain(x0, n, q, &g, mask);
  cgr = pvq_compute_gain(r0, n, q, &gr, mask);
  /* gain_offset is meant to make sure one of the quantized gains has
     exactly the same gain as the reference. */
  icgr = floor(.5+cgr);
  gain_offset = cgr-icgr;
  /* Start search with null case: gain=0, no pulse. */
  qg = 0;
  best_dist = cg*cg;
  noref = 1;
  best_k = 0;
  *itheta = -1;
  *max_theta = 0;
  OD_CLEAR(y, n);
  best_qtheta = 0;
  noref = 1;
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
      double qcg;
      int ts;
      /* Quantized companded gain */
      qcg = i+gain_offset;
      mask_ratio = pvq_interband_masking(*mask_gain, pow(q*qcg, mask), mask);
      /* Set angular resolution (in ra) to match the encoded gain */
      ts = pvq_compute_max_theta(mask_ratio*qcg, mask);
      /* Search for the best angle within a reasonable range. */
      for (j = OD_MAXI(0, (int)floor(.5+theta*2/M_PI*ts)-1);
       j <= OD_MINI(ts-1, (int)ceil(theta*2/M_PI*ts)); j++) {
        double cos_dist;
        double dist;
        double dist_theta;
        double qtheta = pvq_compute_theta(j, ts);
        k = pvq_compute_k(mask_ratio*qcg, qtheta, 0, n, mask);
        cos_dist = pvq_search_double(x, n, k, y_tmp);
        /* See Jmspeex' Journal of Dubious Theoretical Results. */
        dist_theta = 2 - 2*cos(theta - qtheta)
         + sin(theta)*sin(qtheta)*(2 - 2*cos_dist);
        dist = (qcg - cg)*(qcg - cg) + qcg*cg*dist_theta;
        /* Do approximate RDO. */
        dist += lambda*pvq_rate_approx(n, k);
        dist += lambda*(1+log2(ts));
        if (j == 0) dist -= lambda*2.;
        if (i == icgr) dist -= lambda*2.;
        if (dist < best_dist) {
          best_dist = dist;
          qg = i;
          best_k = k;
          best_qtheta = qtheta;
          *itheta = j;
          *max_theta = ts;
          noref = 0;
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
    for (i = OD_MAXI(1, (int)floor(cg) - 1); i <= ceil(cg); i++) {
      double cos_dist;
      double dist;
      double qcg;
      qcg = i;
      mask_ratio = pvq_interband_masking(*mask_gain, pow(q*qcg, mask), mask);
      k = pvq_compute_k(mask_ratio*qcg, -1, 1, n, mask);
      cos_dist = pvq_search_double(x1, n, k, y_tmp);
      /* See Jmspeex' Journal of Dubious Theoretical Results. */
      dist = (qcg - cg)*(qcg - cg) + qcg*cg*(2 - 2*cos_dist);
      /* Do approximate RDO. */
      dist += lambda*pvq_rate_approx(n, k);
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
  /* Remove dimension m if we're using theta. */
  if (!noref) {
    for (i = m; i < n - 1; i++) y[i] = y[i+1];
  }
  /* Synthesize like the decoder would. */
  *mask_gain = pvq_synthesis_partial(x0, y, r, n, noref, qg, gain_offset, theta, m, s, q,
   mask);
  *vk = k;
  /* Encode gain differently depending on whether we use prediction or not. */
  return noref ? qg : neg_interleave(qg, icgr);
}

