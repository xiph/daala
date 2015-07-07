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

#define OD_MASKING_DISABLED 0
#define OD_MASKING_ENABLED 1

/*These tables were generated using compute_basis.c. If OD_FILT_SIZE is
   changed, they have to be regenerated.*/
static const double MAG4[] = {
  0.774125, 0.877780, 0.925934, 0.951682
};
static const double MAG8[] = {
  0.836776, 0.844316, 0.917307, 0.924980,
  0.948172, 0.936507, 0.968913, 0.967917
};
static const double MAG16[] = {
  0.921737, 0.868401, 0.925373, 0.958481,
  0.959319, 0.954073, 0.962690, 0.975782,
  0.974046, 0.967441, 0.968526, 0.979529,
  0.985361, 0.982844, 0.983440, 0.993243
};
static const double MAG32[] = {
  0.961865, 0.926229, 0.935907, 0.950836,
  0.962498, 0.972889, 0.979745, 0.979867,
  0.980251, 0.978192, 0.976537, 0.978706,
  0.981138, 0.984588, 0.987381, 0.987904,
  0.987045, 0.985931, 0.983917, 0.983186,
  0.983692, 0.987112, 0.989474, 0.992827,
  0.992394, 0.991791, 0.991204, 0.990484,
  0.992098, 0.994740, 0.995867, 1.000695
};
static const double MAG64[] = {
  0.980898, 0.962312, 0.964001, 0.965698,
  0.968751, 0.971755, 0.975340, 0.978820,
  0.984460, 0.984691, 0.986874, 0.988789,
  0.989932, 0.990482, 0.990709, 0.990440,
  0.988047, 0.989613, 0.988654, 0.988571,
  0.988610, 0.988864, 0.989477, 0.989908,
  0.993607, 0.991850, 0.992486, 0.993251,
  0.993588, 0.994034, 0.994137, 0.993913,
  0.993533, 0.993267, 0.992192, 0.992293,
  0.992051, 0.991815, 0.991758, 0.992005,
  0.995101, 0.992828, 0.993454, 0.994170,
  0.994912, 0.995532, 0.995988, 0.995961,
  0.994525, 0.996416, 0.995954, 0.995933,
  0.995813, 0.995513, 0.995332, 0.995681,
  0.998732, 0.996535, 0.997001, 0.997798,
  0.998515, 0.999156, 0.999503, 0.999964
};

static const double MAG4_CHROMA_420[] = {
  0.870774, 0.872037, 0.949493, 0.947936
};
static const double MAG8_CHROMA_420[] = {
  0.936496, 0.892830, 0.938452, 0.970087,
  0.974272, 0.967954, 0.974035, 0.990480
};
static const double MAG16_CHROMA_420[] = {
  0.968807, 0.940969, 0.947977, 0.957741,
  0.969762, 0.978644, 0.984885, 0.988009,
  0.987424, 0.985569, 0.984215, 0.984462,
  0.987205, 0.991415, 0.994985, 0.998237
};
static const double MAG32_CHROMA_420[] = {
  0.985068, 0.970006, 0.969893, 0.973192,
  0.973444, 0.975881, 0.979601, 0.981070,
  0.984989, 0.987520, 0.988830, 0.990983,
  0.992376, 0.992884, 0.993447, 0.993381,
  0.993712, 0.994060, 0.993294, 0.992392,
  0.991338, 0.992410, 0.992051, 0.993874,
  0.993488, 0.994162, 0.995318, 0.995925,
  0.997475, 0.999027, 0.998303, 1.001413
};
static const double MAG64_CHROMA_420[] = {
  0.992453, 0.984930, 0.985137, 0.985029,
  0.985514, 0.985784, 0.986269, 0.986854,
  0.989932, 0.987780, 0.988269, 0.989175,
  0.989951, 0.990466, 0.991145, 0.991839,
  0.990773, 0.993191, 0.993618, 0.994221,
  0.994662, 0.995259, 0.995826, 0.995996,
  0.999070, 0.996624, 0.996835, 0.996948,
  0.997022, 0.996973, 0.996993, 0.996996,
  0.996871, 0.996828, 0.996598, 0.996688,
  0.996845, 0.996407, 0.996327, 0.996435,
  0.999173, 0.996216, 0.995981, 0.996173,
  0.996595, 0.996334, 0.996512, 0.996627,
  0.994976, 0.997113, 0.997248, 0.997548,
  0.997943, 0.998121, 0.998291, 0.998687,
  1.001696, 0.999133, 0.999315, 0.999621,
  0.999745, 0.999905, 0.999936, 1.000075
};

const double *OD_BASIS_MAG[2][OD_NBSIZES + 1] = {
  {
    MAG4, MAG8, MAG16, MAG32, MAG64
  },
  {
    MAG4_CHROMA_420, MAG8_CHROMA_420, MAG16_CHROMA_420, MAG32_CHROMA_420,
    MAG64_CHROMA_420
  }
};

/* Quantization matrices for 8x8. For other block sizes, we currently just do
   resampling. */
/* Flat quantization, i.e. optimize for PSNR. */
const int OD_QM8_Q4_QM_FLAT[] = {
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
};
# if 0
/* M1: MPEG2 matrix for inter (which has a dead zone). */
const int OD_QM8_Q4[] = {
  16, 17, 18, 19, 20, 21, 22, 23,
  17, 18, 19, 20, 21, 22, 23, 24,
  18, 19, 20, 21, 22, 23, 24, 25,
  19, 20, 21, 22, 23, 24, 26, 27,
  20, 21, 22, 23, 25, 26, 27, 28,
  21, 22, 23, 24, 26, 27, 28, 30,
  22, 23, 24, 26, 27, 28, 30, 31,
  23, 24, 25, 27, 28, 30, 31, 33};
# endif
# if 0
/* M2: MPEG2 matrix for intra (no dead zone). */
const int OD_QM8_Q4[] = {
  16, 16, 19, 22, 22, 26, 26, 27,
  16, 16, 22, 22, 26, 27, 27, 29,
  19, 22, 26, 26, 27, 29, 29, 35,
  22, 24, 27, 27, 29, 32, 34, 38,
  26, 27, 29, 29, 32, 35, 38, 46,
  27, 29, 34, 34, 35, 40, 46, 56,
  29, 34, 34, 37, 40, 48, 56, 69,
  34, 37, 38, 40, 48, 58, 69, 83
};
# endif
# if 0
/* M3: Taken from dump_psnrhvs. */
const int OD_QM8_Q4[] = {
  16, 16, 17, 20, 24, 29, 36, 42,
  16, 17, 17, 19, 22, 26, 31, 37,
  17, 17, 21, 23, 26, 30, 34, 40,
  20, 19, 23, 28, 31, 35, 39, 45,
  24, 22, 26, 31, 36, 41, 46, 51,
  29, 26, 30, 35, 41, 47, 52, 58,
  36, 31, 34, 39, 46, 52, 59, 66,
  42, 37, 40, 45, 51, 58, 66, 73
};
# endif
# if 1
/* M4: a compromise equal to .5*(M3 + .5*(M2+transpose(M2))) */
const int OD_QM8_Q4_QM_HVS[] = {
  16, 16, 18, 21, 24, 28, 32, 36,
  16, 17, 20, 21, 24, 27, 31, 35,
  18, 20, 24, 25, 27, 31, 33, 38,
  21, 21, 25, 28, 30, 34, 37, 42,
  24, 24, 27, 30, 34, 38, 43, 49,
  28, 27, 31, 34, 38, 44, 50, 58,
  32, 31, 33, 37, 43, 50, 58, 68,
  36, 35, 38, 42, 49, 58, 68, 78
};
#endif

/* These are the PVQ equivalent of quantization matrices, except that
   the values are per-band. */

static const unsigned char OD_FLAT_LUMA_QM_Q4[2][OD_QM_SIZE] = {
/* Flat quantization for PSNR. The DC component isn't 16 because the DC
   magnitude compensation is done here for inter (Haar DC doesn't need it).
   Masking disabled: */
 {
  27, 16,
  23, 16, 16, 16,
  19, 16, 16, 16, 16, 16,
  17, 16, 16, 16, 16, 16, 16, 16
 },
/* The non-flat AC coefficients compensate for the non-linear scaling caused
   by activity masking. The values are currently hand-tuned so that the rate
   of each band remains roughly constant when enabling activity masking
   on intra.
   Masking enabled: */
 {
  27, 16,
  23, 18, 28, 32,
  19, 14, 20, 20, 28, 32,
  17, 11, 16, 14, 16, 16, 23, 28
 }
};

static const unsigned char OD_FLAT_CHROMA_QM_Q4[2][OD_QM_SIZE] = {
/* Chroma quantization is different because of the reduced lapping.
   FIXME: Use the same matrix as luma for 4:4:4.
   Masking disabled: */
 {
  21, 16,
  18, 16, 16, 16,
  17, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
 },
/* The AC part is flat for chroma because it has no activity masking.
   Masking enabled: */
 {
  21, 16,
  18, 16, 16, 16,
  17, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
 }
};

/* No interpolation, always use od_flat_qm_q4, but use a different scale for
   each plane.
   FIXME: Add interpolation and properly tune chroma. */
const od_qm_entry OD_DEFAULT_QMS[2][2][OD_NPLANES_MAX] = {
 /* Masking disabled */
 {{{15, 256, OD_FLAT_LUMA_QM_Q4[OD_MASKING_DISABLED]},
   {15, 448, OD_FLAT_CHROMA_QM_Q4[OD_MASKING_DISABLED]},
   {15, 320, OD_FLAT_CHROMA_QM_Q4[OD_MASKING_DISABLED]}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}},
 /* Masking enabled */
 {{{15, 256, OD_FLAT_LUMA_QM_Q4[OD_MASKING_ENABLED]},
   {15, 448, OD_FLAT_CHROMA_QM_Q4[OD_MASKING_ENABLED]},
   {15, 320, OD_FLAT_CHROMA_QM_Q4[OD_MASKING_ENABLED]}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}}
};

/* Constants for the beta parameter, which controls how activity masking is
   used.
   beta = 1 / (1 - alpha), so when beta is 1, alpha is 0 and activity
   masking is disabled. When beta is 1.5, activity masking is used. Note that
   activity masking is neither used for 4x4 blocks nor for chroma. */

static const double OD_PVQ_BETA4_LUMA[1] = {1.};
static const double OD_PVQ_BETA8_LUMA[4] = {1., 1., 1., 1.};
static const double OD_PVQ_BETA16_LUMA[7] = {1., 1., 1., 1., 1., 1., 1.};
static const double OD_PVQ_BETA32_LUMA[10] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1.};
static const double OD_PVQ_BETA64_LUMA[13] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1., 1., 1., 1.};

static const double OD_PVQ_BETA4_LUMA_MASKING[1] = {1.};
static const double OD_PVQ_BETA8_LUMA_MASKING[4] = {1.5, 1.5, 1.5, 1.5};
static const double OD_PVQ_BETA16_LUMA_MASKING[7] = {1.5, 1.5, 1.5, 1.5, 1.5,
 1.5, 1.5};
static const double OD_PVQ_BETA32_LUMA_MASKING[10] = {1.5, 1.5, 1.5, 1.5,
 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
static const double OD_PVQ_BETA64_LUMA_MASKING[13] = {1.5, 1.5, 1.5, 1.5,
 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

static const double OD_PVQ_BETA4_CHROMA[1] = {1.};
static const double OD_PVQ_BETA8_CHROMA[4] = {1., 1., 1., 1.};
static const double OD_PVQ_BETA16_CHROMA[7] = {1., 1., 1., 1., 1., 1., 1.};
static const double OD_PVQ_BETA32_CHROMA[10] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1.};
static const double OD_PVQ_BETA64_CHROMA[13] = {1., 1., 1., 1., 1., 1., 1.,
 1., 1., 1., 1., 1., 1.};

const double *const OD_PVQ_BETA[2][OD_NPLANES_MAX][OD_NBSIZES + 1] = {
 {{OD_PVQ_BETA4_LUMA, OD_PVQ_BETA8_LUMA,
   OD_PVQ_BETA16_LUMA, OD_PVQ_BETA32_LUMA,
   OD_PVQ_BETA64_LUMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA}},
 {{OD_PVQ_BETA4_LUMA_MASKING, OD_PVQ_BETA8_LUMA_MASKING,
   OD_PVQ_BETA16_LUMA_MASKING, OD_PVQ_BETA32_LUMA_MASKING,
   OD_PVQ_BETA64_LUMA_MASKING},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA}}
};

/* Apply the quantization matrix and the magnitude compensation. We need to
   compensate for the magnitude because lapping causes some basis functions to
   be smaller, so they would end up being quantized too finely (the same
   error in the quantized domain would result in a smaller pixel domain
   error). */
void od_apply_qm(od_coeff *out, int out_stride, od_coeff *in, int in_stride,
 int bs, int dec, int inverse, const int *qm) {
  int i;
  int j;
  for (i = 0; i < 4 << bs; i++) {
    for (j = 0; j < 4 << bs; j++) {
      double mag;
      mag = OD_BASIS_MAG[dec][bs][i]*OD_BASIS_MAG[dec][bs][j];
      if (i == 0 && j == 0) {
        mag = 1;
      }
      else {
        mag /= 0.0625*qm[(i << 1 >> bs)*8 + (j << 1 >> bs)];
      }
      if (inverse) {
        out[i*out_stride + j] = (od_coeff)floor(.5 + in[i*in_stride + j]/mag);
      }
      else {
        out[i*out_stride + j] = (od_coeff)floor(.5 + in[i*in_stride + j]*mag);
      }
    }
  }
}

/* Indexing for the packed quantization matrices. */
int od_qm_get_index(int bs, int band) {
  static const int offsets[OD_NPLANES_MAX] = {0, 2, 6, 12};
  /* The -band/3 term is due to the fact that we force corresponding horizontal
     and vertical bands to have the same quantization. */
  return offsets[bs] + band - band/3;
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
 * @param [in]      beta   activity masking beta param
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
 * @param [in]      g       decoded quantized vector gain
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
    yy += ypulse[i]*(int32_t)ypulse[i];
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

