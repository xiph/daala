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
 */
static void pvq_encode_partition(od_ec_enc *ec,
                                 int qg,
                                 int theta,
                                 int max_theta,
                                 const od_coeff *in,
                                 int n,
                                 int k,
                                 generic_encoder *model,
                                 int *adapt,
                                 int *exg,
                                 int *ext) {

  int adapt_curr[OD_NSB_ADAPT_CTXS] = { 0 };
  int speed = 5;

  generic_encode(ec, &model[theta!=-1], qg, -1, exg, 2);
  if (theta >= 0 && max_theta > 0)
    generic_encode(ec, &model[2], theta, max_theta-1, ext, 2);

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

/** Encode a coefficient block (excepting DC) using PVQ
 *
 * @param [in,out] enc     daala encoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [in]     in      coefficient block to quantize and encode
 * @param [out]    out     quantized coefficient block
 * @param [in]     q       scale/quantizer
 * @param [in]     n       number of coefficients on one side of block
 */
void pvq_encode(daala_enc_ctx *enc,
                od_coeff *ref,
                od_coeff *in,
                od_coeff *out,
                int q,
                int ln,
                const int *qm,
                const double *mask){
  int theta[PVQ_MAX_PARTITIONS];
  int max_theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  int *adapt;
  int *exg;
  int *ext;
  int predflags8;
  int predflags16;
  int nb_bands;
  int i;
  const int *off;
  int size[PVQ_MAX_PARTITIONS];
  double g[PVQ_MAX_PARTITIONS] = {0};
  generic_encoder *model;
  adapt = enc->state.pvq_adapt;
  exg = enc->state.pvq_exg;
  ext = enc->state.pvq_ext;
  model = enc->state.pvq_gain_model;
  nb_bands = od_band_offsets[ln][0];
  off = &od_band_offsets[ln][1];
  for (i = 0; i < nb_bands; i++) size[i] = off[i+1] - off[i];
  for (i = 0; i < nb_bands; i++) {
    if (i == 1) g[1] = g[2] = g[3] = INTER_MASKING*g[0];
    qg[i] = pvq_theta(in+off[i], ref+off[i], size[i], q*qm[i+1] >> 4, out+off[i],
     &theta[i], &max_theta[i], &k[i], &g[i], mask[i]);
  }
  if (ln == 0) {
    od_ec_encode_bool_q15(&enc->ec, theta[0] != -1, PRED4_PROB);
  } else {
    predflags8 = 8*(theta[0] != -1) + 4*(theta[1] != -1) + 2*(theta[2] != -1)
     + (theta[3] != -1);
    od_ec_encode_cdf_q15(&enc->ec, predflags8, pred8_cdf, 16);
    if (ln >= 2) {
      predflags16 = 4*(theta[4] != -1) + 2*(theta[5] != -1)
       + (theta[6] != -1);
      od_ec_encode_cdf_q15(&enc->ec, predflags16, pred16_cdf[predflags8], 8);
    }
  }
  for (i = 0; i < nb_bands; i++) {
    pvq_encode_partition(&enc->ec, qg[i], theta[i], max_theta[i], out+off[i],
     size[i], k[i], model, adapt, exg+i, ext+i);
  }
}
