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

#define EPSILON 1e-30

/* These are the PVQ equivalent of quantization matrices, except that
   the values are per-band. */

#if OD_DISABLE_QM

static const unsigned char od_flat_qm_q4[OD_QM_SIZE] = {
  16, 16,
  16, 16, 16, 16,
  15, 15, 15, 15, 15, 15,
  13, 13, 13, 13, 13, 13, 13, 13
};

const od_qm_entry OD_DEFAULT_QMS[][OD_NPLANES_MAX] = {
  {{15, 256, od_flat_qm_q4},
   {15, 448, od_flat_qm_q4},
   {15, 320, od_flat_qm_q4}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}
};

#else

static const unsigned char od_low_rate_qm_q4[OD_QM_SIZE] = {
  16, 19,
  16, 16, 24, 32,
  16, 15, 17, 19, 23, 31,
  16, 14, 14, 16, 16, 18, 21, 29
};

static const unsigned char od_high_rate_qm_q4[OD_QM_SIZE] = {
  19, 26,
  15, 17, 52, 83,
  11, 11, 13, 16, 39, 65,
  7, 7, 8, 9, 11, 16, 35, 54
};

const od_qm_entry OD_DEFAULT_QMS[][OD_NPLANES_MAX] = {
  {{20, 256, od_high_rate_qm_q4},
   {20, 448, od_high_rate_qm_q4},
   {20, 320, od_high_rate_qm_q4}},
  {{70, 256, od_low_rate_qm_q4},
   {70, 256, od_low_rate_qm_q4},
   {70, 256, od_low_rate_qm_q4}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}
};
#endif

#if OD_DISABLE_MASKING

static const double OD_PVQ_BETA4_LUMA[1] = {1.};
static const double OD_PVQ_BETA8_LUMA[4] = {1., 1., 1., 1.};
static const double OD_PVQ_BETA16_LUMA[7] = {1., 1., 1., 1., 1., 1., 1.};
static const double OD_PVQ_BETA32_LUMA[10] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1.};

#else

static const double OD_PVQ_BETA4_LUMA[1] = {1.};
static const double OD_PVQ_BETA8_LUMA[4] = {1.5, 1.5, 1.5, 1.5};
static const double OD_PVQ_BETA16_LUMA[7] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
static const double OD_PVQ_BETA32_LUMA[10] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5};

#endif

static const double OD_PVQ_BETA4_CHROMA[1] = {1.};
static const double OD_PVQ_BETA8_CHROMA[4] = {1., 1., 1., 1.};
static const double OD_PVQ_BETA16_CHROMA[7] = {1., 1., 1., 1., 1., 1., 1.};
static const double OD_PVQ_BETA32_CHROMA[10] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1.};

const double *const OD_PVQ_BETA[OD_NPLANES_MAX][OD_NBSIZES] = {
  {OD_PVQ_BETA4_LUMA, OD_PVQ_BETA8_LUMA,
   OD_PVQ_BETA16_LUMA, OD_PVQ_BETA32_LUMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA}
};

/* Indexing for the packed quantization matrices. */
int od_qm_get_index(int ln, int band) {
  static const int offsets[OD_NPLANES_MAX] = {0, 2, 6, 12};
  /* The -band/3 term is due to the fact that we force corresponding horizontal
     and vertical bands to have the same quantization. */
  return offsets[ln] + band - band/3;
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
int od_compute_householder(double *r, int n, double gr, int *sign) {
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
void od_apply_householder(double *x, const double *r, int n) {
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

/** Gain companding: raises gain to the power 1/beta for activity masking.
 *
 * @param [in]  g     real (uncompanded) gain
 * @param [in]  q0    uncompanded quality parameter
 * @param [in]  beta  activity masking beta param (exponent)
 * @return            g^(1/beta)
 */
static double od_gain_compand(double g, int q0, double beta) {
  if (beta == 1) return g/q0;
  else return OD_COMPAND_SCALE*pow(g*OD_COMPAND_SCALE_1, 1./beta)/q0;
}

/** Gain expanding: raises gain to the power beta for activity masking.
 *
 * @param [in]  cg    companded gain
 * @param [in]  q0    uncompanded quality parameter
 * @param [in]  beta  activity masking beta param (exponent)
 * @return            g^beta
 */
double od_gain_expand(double cg, int q0, double beta) {
  if (beta == 1) return cg*q0;
  else if (beta == 1.5) {
    cg *= q0*OD_COMPAND_SCALE_1;
    return OD_COMPAND_SCALE*cg*sqrt(cg);
  }
  else {
    return OD_COMPAND_SCALE*pow(cg*q0*OD_COMPAND_SCALE_1, beta);
  }
}

/** Computes the raw and quantized/companded gain of a given input
 * vector
 *
 * @param [in]      x      vector of input data
 * @param [in]      n      number of elements in vector x
 * @param [in]      q0     quantizer
 * @param [out]     g      raw gain
 * @param [in]      beta   activity masking beta param
 * @return                 quantized/companded gain
 */
double od_pvq_compute_gain(od_coeff *x, int n, int q0, double *g, double beta){
  int i;
  double acc=0;
  for (i = 0; i < n; i++) acc += x[i]*(double)x[i];
  *g = sqrt(acc);
  /* Normalize gain by quantization step size and apply companding
     (if ACTIVITY != 1). */
  return od_gain_compand(*g, q0, beta);
}

/** Compute theta quantization range from quantized/companded gain
 *
 * @param [in]      qcg    quantized companded gain value
 * @return                 max theta value
 */
int od_pvq_compute_max_theta(double qcg, double beta){
  /* Set angular resolution (in ra) to match the encoded gain */
  int ts = (int)floor(.5 + qcg*M_PI/(2*beta));
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
double od_pvq_compute_theta(int t, int max_theta) {
  if (max_theta != 0) return OD_MINI(t, max_theta - 1)*.5*M_PI/max_theta;
  return 0;
}

/** Compute the number of pulses used for PVQ encoding a vector from
 * available metrics (encode and decode side)
 *
 * @param [in]      qcg        quantized companded gain value
 * @param [in]      itheta     quantizized PVQ error angle theta
 * @param [in]      theta      PVQ error angle theta
 * @param [in]      noref      indicates present or lack of reference
 *                             (prediction)
 * @param [in]      n          number of elements to be coded
 * @param [in]      beta       activity masking beta param
 * @param [in]      nodesync   do not use info that depend on the reference
 * @return                     number of pulses to use for coding
 */
int od_pvq_compute_k(double qcg, int itheta, double theta, int noref, int n,
 double beta, int nodesync) {
  if (noref) {
    if (qcg == 0) return 0;
    if (n == 15 && qcg == 1 && beta > 1.25) return 1;
    else return OD_MAXI(1, (int)floor(.5 + (qcg - .2)*sqrt((n+3)/2)/beta));
  }
  else {
    if (itheta == 0) return 0;
    /* Sets K according to gain and theta, based on the high-rate
       PVQ distortion curves (see PVQ document). Low-rate will have to be
       perceptually tuned anyway. We subtract 0.2 from the radius as an
       approximation for the fact that the coefficients aren't identically
       distributed within a band so at low gain the number of dimensions that
       are likely to have a pulse is less than n. */
    if (nodesync) {
      return OD_MAXI(1, (int)floor(.5 + (itheta - .2)*sqrt((n + 2)/2)));
    }
    else {
      return OD_MAXI(1, (int)floor(.5 + (qcg*sin(theta) - .2)*
       sqrt((n + 2)/2)/beta));
    }
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
 * @param [in]      theta   decoded theta (prediction error)
 * @param [in]      m       alignment dimension of Householder reflection
 * @param [in]      s       sign of Householder reflection
 */
void od_pvq_synthesis_partial(od_coeff *xcoeff, const od_coeff *ypulse,
 const double *r, int n, int noref, double g, double theta, int m, int s) {
  int i;
  int yy;
  double scale;
  int nn;
  OD_ASSERT(g != 0);
  nn = n-(!noref); /* when noref==0, vector in is sized n-1 */
  yy = 0;
  for (i = 0; i < nn; i++)
    yy += ypulse[i]*(ogg_int32_t)ypulse[i];
  if (yy == 0) scale = 0;
  else scale = g/sqrt(yy);
  if (noref) {
    for (i = 0; i < n; i++) xcoeff[i] = (od_coeff)floor(.5 + ypulse[i]*scale);
  }
  else{
    double x[MAXN];
    scale *= sin(theta);
    for (i = 0; i < m; i++)
      x[i] = ypulse[i]*scale;
    x[m] = -s*g*cos(theta);
    for (i = m; i < nn; i++)
      x[i+1] = ypulse[i]*scale;
    od_apply_householder(x, r, n);
    for (i = 0; i < n; i++) xcoeff[i] = (od_coeff)floor(.5 + x[i]);
  }
}

