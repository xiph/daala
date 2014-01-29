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

  generic_encode(ec, model, qg, exg, 2);
  if (theta >= 0 && max_theta > 0)
    generic_encode(ec, model, theta, ext, 2);

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
                int n){
  int theta[PVQ_MAX_PARTITIONS];
  int max_theta[PVQ_MAX_PARTITIONS];
  int qg[PVQ_MAX_PARTITIONS];
  int k[PVQ_MAX_PARTITIONS];
  int *adapt;
  int *exg;
  int *ext;
  int predflags8;
  int predflags16;
  generic_encoder *model;
  adapt = enc->state.pvq_adapt;
  exg = enc->state.pvq_exg;
  ext = enc->state.pvq_ext;
  model = &enc->state.pvq_gain_model;

  qg[0] = pvq_theta(in+1, ref+1, 15, q, out+1,
                    &theta[0], &max_theta[0], &k[0]);

  if (n==4){

    od_ec_encode_bool_q15(&enc->ec, theta[0] != -1, PRED4_PROB);
    pvq_encode_partition(&enc->ec, qg[0], theta[0], max_theta[0], out+1,
                        15, k[0], model, adapt, exg, ext);

  }
  else{

    qg[1] = pvq_theta(in+16, ref+16, 8, q, out+16,
                      &theta[1], &max_theta[1], &k[1]);
    qg[2] = pvq_theta(in+24, ref+24, 8, q, out+24,
                      &theta[2], &max_theta[2], &k[2]);
    qg[3] = pvq_theta(in+32, ref+32, 32, q, out+32,
                      &theta[3], &max_theta[3], &k[3]);
    predflags8 = 8*(theta[0] != -1) + 4*(theta[1] != -1) + 2*(theta[2] != -1)
      + (theta[3] != -1);
    od_ec_encode_cdf_q15(&enc->ec, predflags8, pred8_cdf, 16);

    if (n >= 16) {
      qg[4] = pvq_theta(in+64, ref+64, 32, q, out+64,
                        &theta[4], &max_theta[4], &k[4]);
      qg[5] = pvq_theta(in+96, ref+96, 32, q, out+96,
                        &theta[5], &max_theta[5], &k[5]);
      qg[6] = pvq_theta(in+128, ref+128, 128, q, out+128,
                          &theta[6], &max_theta[6], &k[6]);

      predflags16 = 4*(theta[4] != -1) + 2*(theta[5] != -1)
        + (theta[6] != -1);
      od_ec_encode_cdf_q15(&enc->ec, predflags16, pred16_cdf[predflags8], 8);
    }

    pvq_encode_partition(&enc->ec, qg[0], theta[0], max_theta[0], out+1,
                         15, k[0], model, adapt, exg, ext);
    pvq_encode_partition(&enc->ec, qg[1], theta[1], max_theta[1], out+16,
                         8, k[1], model, adapt, exg+1, ext+1);
    pvq_encode_partition(&enc->ec, qg[2], theta[2], max_theta[2], out+24,
                         8, k[2], model, adapt, exg+2, ext+2);
    pvq_encode_partition(&enc->ec, qg[3], theta[3], max_theta[3], out+32,
                         32, k[3], model, adapt, exg+3, ext+3);

    if (n >= 16) {
      pvq_encode_partition(&enc->ec, qg[4], theta[4], max_theta[4], out+64,
                           32, k[4], model, adapt, exg+4, ext+4);
      pvq_encode_partition(&enc->ec, qg[5], theta[5], max_theta[5], out+96,
                           32, k[5], model, adapt, exg+5, ext+5);
      pvq_encode_partition(&enc->ec, qg[6], theta[6], max_theta[6], out+128,
                           128, k[6], model, adapt, exg+6, ext+6);
    }
  }
}
