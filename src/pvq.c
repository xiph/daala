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
#include "partition.h"
#include <stdlib.h>
#include <stdio.h>
#include "logging.h"
#include <math.h>
#include <string.h>
#include "filter.h"

/*These tables were generated using compute_basis.c. If OD_FILT_SIZE is
   changed, they have to be regenerated.*/
static const double MAG4[] = {
  0.870774, 0.872037, 0.949493, 0.947936
};
static const double MAG8[] = {
  0.936496, 0.892830, 0.938452, 0.970087,
  0.974272, 0.967954, 0.974035, 0.990480
};
static const double MAG16[] = {
  0.968807, 0.940969, 0.947977, 0.957741,
  0.969762, 0.978644, 0.984885, 0.988009,
  0.987424, 0.985569, 0.984215, 0.984462,
  0.987205, 0.991415, 0.994985, 0.998237
};
static const double MAG32[] = {
  0.985068, 0.970006, 0.969893, 0.973192,
  0.973444, 0.975881, 0.979601, 0.981070,
  0.984989, 0.987520, 0.988830, 0.990983,
  0.992376, 0.992884, 0.993447, 0.993381,
  0.993712, 0.994060, 0.993294, 0.992392,
  0.991338, 0.992410, 0.992051, 0.993874,
  0.993488, 0.994162, 0.995318, 0.995925,
  0.997475, 0.999027, 0.998303, 1.001413,
};
static const double MAG64[] = {
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
const int OD_QM8_Q4_FLAT[] = {
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
const int OD_QM8_Q4_HVS[] = {
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

void od_adapt_pvq_ctx_reset(od_pvq_adapt_ctx *state, int is_keyframe) {
  od_pvq_codeword_ctx *ctx;
  int i;
  int pli;
  int bs;
  ctx = &state->pvq_codeword_ctx;
  generic_model_init(&state->pvq_param_model[0]);
  generic_model_init(&state->pvq_param_model[1]);
  generic_model_init(&state->pvq_param_model[2]);
  for (i = 0; i < 2*OD_NBSIZES; i++) {
    ctx->pvq_adapt[4*i + OD_ADAPT_K_Q8] = 384;
    ctx->pvq_adapt[4*i + OD_ADAPT_SUM_EX_Q8] = 256;
    ctx->pvq_adapt[4*i + OD_ADAPT_COUNT_Q8] = 104;
    ctx->pvq_adapt[4*i + OD_ADAPT_COUNT_EX_Q8] = 128;
  }
  ctx->pvq_k1_increment = 128;
  OD_CDFS_INIT(ctx->pvq_k1_cdf, ctx->pvq_k1_increment);
  for (pli = 0; pli < OD_NPLANES_MAX; pli++) {
    for (bs = 0; bs < OD_NBSIZES; bs++)
    for (i = 0; i < PVQ_MAX_PARTITIONS; i++) {
      state->pvq_exg[pli][bs][i] = 2 << 16;
    }
  }
  for (i = 0; i < OD_NBSIZES*PVQ_MAX_PARTITIONS; i++) {
    state->pvq_ext[i] = is_keyframe ? 24576 : 2 << 16;
  }
  state->pvq_gaintheta_increment = 128;
  OD_CDFS_INIT(state->pvq_gaintheta_cdf, state->pvq_gaintheta_increment >> 2);
  state->pvq_skip_dir_increment = 128;
  OD_CDFS_INIT(state->pvq_skip_dir_cdf, state->pvq_skip_dir_increment >> 2);

}

/* QMs are arranged from smallest to largest blocksizes, for decimation=0,
   and again for decimation=1.*/
int od_qm_offset(int bs, int xydec)
{
    return bs*2*OD_QM_BSIZE + xydec*OD_QM_BSIZE;
}

/* Initialize the quantization matrix with the magnitude compensation applied.
   We need to compensate for the magnitude because lapping causes some basis
   functions to be smaller, so they would end up being quantized too finely
   (the same error in the quantized domain would result in a smaller pixel
   domain error). */
void od_init_qm(int16_t *x, int16_t *x_inv, const int *qm) {
  int i;
  int j;
  int16_t y[OD_QM_BSIZE];
  int16_t y_inv[OD_QM_BSIZE];
  int16_t *x1;
  int16_t *x1_inv;
  int off;
  int bs;
  int xydec;
  for (bs = 0; bs < OD_NBSIZES; bs++) {
    for (xydec = 0; xydec < 2; xydec++) {
      off = od_qm_offset(bs, xydec);
      x1 = x + off;
      x1_inv = x_inv + off;
      for (i = 0; i < 4 << bs; i++) {
        for (j = 0; j < 4 << bs; j++) {
          double mag;
#if OD_DEBLOCKING || OD_DISABLE_FILTER
          mag = 1.0;
#else
          mag = OD_BASIS_MAG[xydec][bs][i]*OD_BASIS_MAG[xydec][bs][j];
#endif
          if (i == 0 && j == 0) {
            mag = 1.0;
          }
          else {
            mag /= 0.0625*qm[(i << 1 >> bs)*8 + (j << 1 >> bs)];
            OD_ASSERT(mag > 0.0);
          }
          /*Convert to fit in 16 bits.*/
          y[i*(4 << bs) + j] = (int16_t)OD_MINI(OD_QM_SCALE_MAX,
           (int32_t)floor(.5 + mag*OD_QM_SCALE));
          y_inv[i*(4 << bs) + j] = (int16_t)floor(.5
           + OD_QM_SCALE*OD_QM_INV_SCALE/(double)y[i*(4 << bs) + j]);
        }
      }
      od_raster_to_coding_order_16(x1, 4 << bs, y, 4 << bs);
      od_raster_to_coding_order_16(x1_inv, 4 << bs, y_inv, 4 << bs);
    }
  }
}

/* Indexing for the packed quantization matrices. */
int od_qm_get_index(int bs, int band) {
  /* The -band/3 term is due to the fact that we force corresponding horizontal
     and vertical bands to have the same quantization. */
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  return bs*(bs + 1) + band - band/3;
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
 * @param [in]      qm     QM with magnitude compensation
 * @return                 quantized/companded gain
 */
double od_pvq_compute_gain(od_coeff *x, int n, int q0, double *g, double beta,
 const int16_t *qm){
  int i;
  double acc=0;
  for (i = 0; i < n; i++) {
    acc += x[i]*(double)x[i]*qm[i]*OD_QM_SCALE_1*
     qm[i]*OD_QM_SCALE_1;
  }
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
 * @param [in]      qm_inv  inverse of the QM with magnitude compensation
 */
void od_pvq_synthesis_partial(od_coeff *xcoeff, const od_coeff *ypulse,
 const double *r, int n, int noref, double g, double theta, int m, int s,
 const int16_t *qm_inv) {
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
    for (i = 0; i < n; i++) {
      xcoeff[i] = (od_coeff)floor(.5
       + (ypulse[i]*scale)*(qm_inv[i]*OD_QM_INV_SCALE_1));
    }
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
    for (i = 0; i < n; i++) {
      xcoeff[i] = (od_coeff)floor(.5 + (x[i]*(qm_inv[i]*OD_QM_INV_SCALE_1)));
    }
  }
}

