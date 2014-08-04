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
 * @param [in]     is_keyframe whether we're encoding a keyframe
 */
static void pvq_encode_partition(od_ec_enc *ec,
                                 int qg,
                                 int theta,
                                 int max_theta,
                                 const od_coeff *in,
                                 int n,
                                 int k,
                                 generic_encoder model[3],
                                 int *adapt,
                                 int *exg,
                                 int *ext,
                                 int is_keyframe) {

  int adapt_curr[OD_NSB_ADAPT_CTXS] = { 0 };
  int speed = 5;
  int noref;
  noref = (theta == -1);
  generic_encode(ec, &model[!noref], qg, -1, exg, 2);
  if (!noref && max_theta > 1) {
    if (is_keyframe) {
      int tmp;
      tmp = max_theta**ext;
      generic_encode(ec, &model[2], theta, max_theta-1, &tmp, 2);
      /* Adapt expectation as fraction of max_theta */
      *ext += (theta*65536/max_theta - *ext) >> 5;
    }
    else generic_encode(ec, &model[2], theta, max_theta-1, ext, 2);
  }
  laplace_encode_vector(ec, in, n - (theta >= 0), k, adapt_curr, adapt);

  if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
    adapt[OD_ADAPT_K_Q8] += 256*adapt_curr[OD_ADAPT_K_Q8] -
     adapt[OD_ADAPT_K_Q8]>>speed;
    adapt[OD_ADAPT_SUM_EX_Q8] += adapt_curr[OD_ADAPT_SUM_EX_Q8] -
     adapt[OD_ADAPT_SUM_EX_Q8]>>speed;
  }
  if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
    adapt[OD_ADAPT_COUNT_Q8] += adapt_curr[OD_ADAPT_COUNT_Q8]-
     adapt[OD_ADAPT_COUNT_Q8]>>speed;
    adapt[OD_ADAPT_COUNT_EX_Q8] += adapt_curr[OD_ADAPT_COUNT_EX_Q8]-
     adapt[OD_ADAPT_COUNT_EX_Q8]>>speed;
  }
}

void code_flag(od_ec_enc *ec, int val, unsigned *prob0)
{
  od_ec_encode_bool_q15(ec, val, *prob0);
  if (val) {
    *prob0 = *prob0 - (*prob0 >> OD_NOREF_ADAPT_SPEED);
  }
  else {
    *prob0 = *prob0 + ((32768 - *prob0) >> OD_NOREF_ADAPT_SPEED);
  }
}

/** Encode a coefficient block (excepting DC) using PVQ
 *
 * @param [in,out] enc     daala encoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [in]     in      coefficient block to quantize and encode
 * @param [out]    out     quantized coefficient block
 * @param [in]     q       scale/quantizer
 * @param [in]     ln      log of the block size minus two
 * @param [in]     qm      per-band quantization matrix
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     is_keyframe whether we're encoding a keyframe
 */
void pvq_encode(daala_enc_ctx *enc,
                od_coeff *ref,
                od_coeff *in,
                od_coeff *out,
                int q,
                int pli,
                int ln,
                const int *qm,
                const double *beta,
                const double *inter_band,
                int is_keyframe){
  int theta[PVQ_MAX_PARTITIONS];
  int max_theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  od_coeff y[16*16];
  int *adapt;
  int *exg;
  int *ext;
  int nb_bands;
  int i;
  const int *off;
  int size[PVQ_MAX_PARTITIONS];
  double g[PVQ_MAX_PARTITIONS] = {0};
  generic_encoder *model;
  unsigned *noref_prob;
  double skip_diff;
  unsigned tell;
  ogg_uint16_t *skip_cdf;
  od_rollback_buffer buf;
  adapt = enc->adapt.pvq_adapt;
  exg = &enc->adapt.pvq_exg[pli][ln][0];
  ext = enc->adapt.pvq_ext + ln*PVQ_MAX_PARTITIONS;
  noref_prob = enc->adapt.pvq_noref_prob + ln*PVQ_MAX_PARTITIONS;
  skip_cdf = enc->adapt.skip_cdf[pli];
  model = enc->adapt.pvq_param_model;
  nb_bands = od_band_offsets[ln][0];
  off = &od_band_offsets[ln][1];
  tell = 0;
  for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
  skip_diff = 0;
  for (i = 0; i < nb_bands; i++) {
    int j;
    double mask;
    mask = 0;
    for (j = 0; j < i; j++) mask += *inter_band++*g[j];
    g[i] = mask;
    qg[i] = pvq_theta(out + off[i], in + off[i], ref + off[i], size[i],
     OD_MAXI(1, q*qm[i + 1] >> 4), y + off[i], &theta[i], &max_theta[i],
     &k[i], &g[i], beta[i], &skip_diff);
  }
  if (!is_keyframe) {
    od_encode_checkpoint(enc, &buf);
    /* Code as if we're not skipping. */
    od_encode_cdf_adapt(&enc->ec, (out[0] != 0), skip_cdf,
     4, enc->adapt.skip_increment);
    /* Excluding skip flag from the rate since it's minor and would be prone
       to greedy decision issues. */
    tell = od_ec_enc_tell_frac(&enc->ec);
  }
  if (!is_keyframe && ln > 0) {
    int id;
    id = 0;
    /* Jointly code the noref flags for the first 4 bands. */
    for (i = 0; i < 4; i++) id = (id << 1) + (theta[i] == -1);
    od_encode_cdf_adapt(&enc->ec, id, enc->adapt.pvq_noref_joint_cdf[ln - 1],
     16, enc->adapt.pvq_noref_joint_increment);
    if (ln >= 2) {
      int nb_norefs;
      id = 0;
      nb_norefs = 0;
      /* Context for the last 3 bands is how many of the first 4 bands are
         noref. */
      for (i = 0; i < 4; i++) nb_norefs += (theta[i] == -1);
      for (i = 0; i < 3; i++) id = (id << 1) + (theta[i + 4] == -1);
      od_encode_cdf_adapt(&enc->ec, id,
       enc->adapt.pvq_noref2_joint_cdf[nb_norefs], 8,
       enc->adapt.pvq_noref_joint_increment);
    }
  }
  else {
    for (i = 0; i < nb_bands; i++) {
      if (!(is_keyframe && vector_is_null(ref + off[i], size[i])))
        code_flag(&enc->ec, theta[i] != -1, &noref_prob[i]);
    }
  }
  for (i = 0; i < nb_bands; i++) {
    pvq_encode_partition(&enc->ec, qg[i], theta[i], max_theta[i], y + off[i],
      size[i], k[i], model, adapt, exg + i, ext + i, is_keyframe);
  }
  if (!is_keyframe) {
    tell = od_ec_enc_tell_frac(&enc->ec) - tell;
    if (skip_diff < OD_PVQ_LAMBDA/8*tell) {
      /* We decide to skip, roll back everything as it was before. */
      od_encode_rollback(enc, &buf);
      od_encode_cdf_adapt(&enc->ec, 2 + (out[0] != 0), skip_cdf,
       4, enc->adapt.skip_increment);
      for (i = 1; i < 1 << (2*ln + 4); i++) out[i] = ref[i];
    }
  }
}
