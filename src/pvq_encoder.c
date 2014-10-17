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
 * @param [in]     nodesync   do not use info that depend on the reference
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
                                 int nodesync) {

  int adapt_curr[OD_NSB_ADAPT_CTXS] = { 0 };
  int speed = 5;
  int noref;
  noref = (theta == -1);
  generic_encode(ec, &model[!noref], qg, -1, exg, 2);
  if (!noref && (max_theta > 1 || nodesync)) {
    generic_encode(ec, &model[2], theta, nodesync ? -1 : max_theta - 1, ext,
     2);
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
void pvq_encode(daala_enc_ctx *enc,
                od_coeff *ref,
                od_coeff *in,
                od_coeff *out,
                int q,
                int pli,
                int ln,
                const int *qm,
                const double *beta,
                int robust,
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
  generic_encoder *model;
  unsigned *noref_prob;
  double skip_diff;
  unsigned tell;
  ogg_uint16_t *skip_cdf;
  od_rollback_buffer buf;
  int dc_quant;
  int use_cfl;
  int flip;
  adapt = enc->state.adapt.pvq_adapt;
  exg = &enc->state.adapt.pvq_exg[pli][ln][0];
  ext = enc->state.adapt.pvq_ext + ln*PVQ_MAX_PARTITIONS;
  noref_prob = enc->state.adapt.pvq_noref_prob + ln*PVQ_MAX_PARTITIONS;
  skip_cdf = enc->state.adapt.skip_cdf[pli];
  model = enc->state.adapt.pvq_param_model;
  nb_bands = od_band_offsets[ln][0];
  off = &od_band_offsets[ln][1];
  dc_quant = OD_MAXI(1, q*qm[0] >> 4);
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
  use_cfl = 0;
  for (i = 0; i < nb_bands; i++) {
    qg[i] = pvq_theta(out + off[i], in + off[i], ref + off[i], size[i],
     OD_MAXI(1, q*qm[i + 1] >> 4), y + off[i], &theta[i], &max_theta[i],
     &k[i], beta[i], &skip_diff, robust, is_keyframe, pli);
    if (pli!=0 && is_keyframe && theta[i] != -1) use_cfl = 1;
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
  if (!is_keyframe && ln > 0) {
    int id;
    id = 0;
    /* Jointly code the noref flags for the first 4 bands. */
    for (i = 0; i < 4; i++) id = (id << 1) + (theta[i] == -1);
    od_encode_cdf_adapt(&enc->ec, id, enc->state.adapt.pvq_noref_joint_cdf[ln - 1],
     16, enc->state.adapt.pvq_noref_joint_increment);
    if (ln >= 2) {
      int nb_norefs;
      id = 0;
      nb_norefs = 0;
      /* Context for the last 3 bands is how many of the first 4 bands are
         noref. */
      for (i = 0; i < 4; i++) nb_norefs += (theta[i] == -1);
      for (i = 0; i < 3; i++) id = (id << 1) + (theta[i + 4] == -1);
      od_encode_cdf_adapt(&enc->ec, id,
       enc->state.adapt.pvq_noref2_joint_cdf[nb_norefs], 8,
       enc->state.adapt.pvq_noref_joint_increment);
    }
  }
  else {
    for (i = 0; i < nb_bands; i++) {
      code_flag(&enc->ec, theta[i] != -1, &noref_prob[i]);
    }
  }
  if (use_cfl) {
    /* We could eventually do some smarter entropy coding here, but it would
       have to be good enough to overcome the overhead of the entropy coder.
       An early attempt using a "toogle" flag with simple adaptation wasn't
       worth the trouble. */
    od_ec_enc_bits(&enc->ec, flip, 1);
  }
  for (i = 0; i < nb_bands; i++) {
    pvq_encode_partition(&enc->ec, qg[i], theta[i], max_theta[i], y + off[i],
     size[i], k[i], model, adapt, exg + i, ext + i, robust || is_keyframe);
  }
  if (!is_keyframe) {
    tell = od_ec_enc_tell_frac(&enc->ec) - tell;
    if (skip_diff < OD_PVQ_LAMBDA/8*tell) {
      double dc_rate;
      dc_rate = -OD_LOG2((double)(skip_cdf[3]-skip_cdf[2])/
       (double)(skip_cdf[2]-skip_cdf[1]));
      out[0] = od_rdo_quant(in[0] - ref[0], dc_quant, dc_rate);
      /* We decide to skip, roll back everything as it was before. */
      od_encode_rollback(enc, &buf);
      od_encode_cdf_adapt(&enc->ec, 2 + (out[0] != 0), skip_cdf,
       4, enc->state.adapt.skip_increment);
      for (i = 1; i < 1 << (2*ln + 4); i++) out[i] = ref[i];
    }
  }
}
