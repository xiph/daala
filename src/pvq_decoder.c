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
#include "decint.h"
#include "entdec.h"
#include "entcode.h"
#include "logging.h"
#include "laplace_code.h"
#include "pvq_code.h"

/** Decodes a single vector of integers (eg, a partition within a
 *  coefficient block) encoded using PVQ
 *
 * @param [in,out] ec      range encoder
 * @param [in]     q       scale/quantizer
 * @param [in]     n       number of coefficients in partition
 * @param [in,out] model   entropy decoder state
 * @param [in,out] adapt   adaptation context
 * @param [in,out] exg     ExQ16 expectation of decoded gain value
 * @param [in,out] ext     ExQ16 expectation of decoded theta value
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [in]     noref   boolean indicating absence of reference
 * @param [in,out] mask_gain input masking from other bands, output masking for
 *                           other bands
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     is_keyframe whether we're encoding a keyframe
 */
static void pvq_decode_partition(od_ec_dec *ec,
                                 int q0,
                                 int n,
                                 generic_encoder model[3],
                                 int *adapt,
                                 int *exg,
                                 int *ext,
                                 od_coeff *ref,
                                 od_coeff *out,
                                 int noref,
                                 double *mask_gain,
                                 double beta,
                                 int is_keyframe) {
  int adapt_curr[OD_NSB_ADAPT_CTXS] = {0};
  int speed;
  int k;
  double qcg;
  int max_theta;
  int itheta;
  double theta;
  double gr;
  double gain_offset;
  od_coeff y[1024];
  double r[1024];
  int qg;
  double q;
  double mask_ratio;
  /* Quantization step calibration to account for the activity masking. */
  q = q0*pow(256<<OD_COEFF_SHIFT, 1./beta - 1);
  speed = 5;
  theta = 0;
  gr = 0;
  gain_offset = 0;

  /* read quantized gain */
  qg = generic_decode(ec, &model[!noref], -1, exg, 2);

  if(!noref){
    /* we have a reference; compute its gain */
    double cgr;
    int icgr;
    int i;
    cgr = pvq_compute_gain(ref, n, q, &gr, beta);
    icgr = floor(.5+cgr);
    /* quantized gain is interleave encoded when there's a reference;
       deinterleave it now */
    qg = neg_deinterleave(qg, icgr);
    gain_offset = cgr-icgr;
    qcg = qg + gain_offset;
    mask_ratio = pvq_interband_masking(*mask_gain, pow(q*qcg, 2*beta), beta);
    /* read and decode first-stage PVQ error theta */
    max_theta = pvq_compute_max_theta(mask_ratio*qcg, beta);
    if (max_theta > 1) {
      if (is_keyframe) {
        int tmp;
        tmp = max_theta**ext;
        itheta = generic_decode(ec, &model[2], max_theta-1, &tmp, 2);
        /* Adapt expectation as fraction of max_theta */
        *ext += (itheta*65536/max_theta - *ext) >> 5;
      }
      else itheta = generic_decode(ec, &model[2], max_theta - 1, ext, 2);
    }
    else itheta = 0;
    theta = pvq_compute_theta(itheta, max_theta);
    for (i = 0; i < n; i++) r[i] = ref[i];
  }
  else{
    qcg=qg;
    mask_ratio = pvq_interband_masking(*mask_gain, pow(q*qcg, 2*beta), beta);
  }

  if (qg != 0) {
    k = pvq_compute_k(mask_ratio*qcg, theta, noref, n, beta);
    /* when noref==0, y is actually size n-1 */
    laplace_decode_vector(ec, y, n-(!noref), k, adapt_curr, adapt);
  } else {
    OD_CLEAR(y, n);
  }
  *mask_gain = pvq_synthesis(out, y, r, n, gr, noref, qg, gain_offset, theta,
   q, beta);

  if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
    adapt[OD_ADAPT_K_Q8]
      += 256*adapt_curr[OD_ADAPT_K_Q8]-adapt[OD_ADAPT_K_Q8]>>speed;
    adapt[OD_ADAPT_SUM_EX_Q8]
      += adapt_curr[OD_ADAPT_SUM_EX_Q8]-adapt[OD_ADAPT_SUM_EX_Q8]>>speed;
  }
  if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
    adapt[OD_ADAPT_COUNT_Q8]
      += adapt_curr[OD_ADAPT_COUNT_Q8]-adapt[OD_ADAPT_COUNT_Q8]>>speed;
    adapt[OD_ADAPT_COUNT_EX_Q8]
      += adapt_curr[OD_ADAPT_COUNT_EX_Q8]-adapt[OD_ADAPT_COUNT_EX_Q8]>>speed;
  }
}

static int decode_flag(od_ec_dec *ec, unsigned *prob0)
{
  int val;
  val = od_ec_decode_bool_q15(ec, *prob0);
  if (val) {
    *prob0 = *prob0 - (*prob0 >> OD_NOREF_ADAPT_SPEED);
  }
  else {
    *prob0 = *prob0 + ((32768 - *prob0) >> OD_NOREF_ADAPT_SPEED);
  }
  return val;
}

/** Decodes a coefficient block (except for DC) encoded using PVQ
 *
 * @param [in,out] dec     daala decoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [in]     q       quantizer
 * @param [in]     ln      log of the block size minus two
 * @param [in]     qm      per-band quantization matrix
 * @param [in]     beta    per-band activity masking beta param
 * @param [in]     is_keyframe whether we're encoding a keyframe
 */
void pvq_decode(daala_dec_ctx *dec,
                od_coeff *ref,
                od_coeff *out,
                int q,
                int pli,
                int ln,
                const int *qm,
                const double *beta,
                const double *inter_band,
                int is_keyframe){

  int noref[PVQ_MAX_PARTITIONS];
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
  int skip;
  adapt = dec->adapt.pvq_adapt;
  exg = &dec->adapt.pvq_exg[pli][ln][0];
  ext = dec->adapt.pvq_ext + ln*PVQ_MAX_PARTITIONS;
  noref_prob = dec->adapt.pvq_noref_prob + ln*PVQ_MAX_PARTITIONS;
  model = dec->adapt.pvq_param_model;
  nb_bands = od_band_offsets[ln][0];
  off = &od_band_offsets[ln][1];
  if (is_keyframe) skip = 0;
  else {
    skip = od_decode_cdf_adapt(&dec->ec, dec->adapt.skip_cdf[pli], 4,
     dec->adapt.skip_increment);
    out[0] = skip&1;
    skip >>= 1;
  }
  if (skip) {
    for (i = 1; i < 1 << (2*ln + 4); i++) out[i] = ref[i];
  }
  else {
    for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
    if (!is_keyframe && ln > 0) {
      int id;
      id = od_decode_cdf_adapt(&dec->ec,
       dec->adapt.pvq_noref_joint_cdf[ln - 1], 16,
       dec->adapt.pvq_noref_joint_increment);
      for (i = 0; i < 4; i++) noref[i] = (id >> (3 - i)) & 1;
      if (ln >= 2) {
        int nb_norefs;
        nb_norefs = 0;
        for (i = 0; i < 4; i++) nb_norefs += noref[i];
        id = od_decode_cdf_adapt(&dec->ec,
         dec->adapt.pvq_noref2_joint_cdf[nb_norefs], 8,
         dec->adapt.pvq_noref_joint_increment);
        for (i = 0; i < 3; i++) noref[i + 4] = (id >> (2 - i)) & 1;
      }
    }
    else {
      for (i = 0; i < nb_bands; i++) {
        if (is_keyframe && vector_is_null(ref + off[i], size[i])) noref[i] = 1;
        else noref[i] = !decode_flag(&dec->ec, &noref_prob[i]);
      }
    }
    for (i = 0; i < nb_bands; i++) {
      int j;
      double mask;
      mask = 0;
      for (j = 0; j < i; j++) mask += *inter_band++*g[j];
      g[i] = mask;
      pvq_decode_partition(&dec->ec, OD_MAXI(1, q*qm[i + 1] >> 4), size[i],
       model, adapt, exg + i, ext + i, ref + off[i], out + off[i], noref[i],
       &g[i], beta[i], is_keyframe);
    }
  }
}
