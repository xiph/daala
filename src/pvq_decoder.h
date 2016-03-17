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

#if !defined(_pvq_decoder_H)
# define _pvq_decoder_H (1)
# include "internal.h"
# include "filter.h"
# include "pvq.h"
# include "entdec.h"
# include "decint.h"

void od_decode_band_pvq_splits(od_ec_dec *ec, od_pvq_codeword_ctx *adapt,
 od_coeff *y, int n, int k, int level);

#if OD_ACCOUNTING
# define laplace_decode_special(dec, decay, max, str) laplace_decode_special_(dec, decay, max, str)
# define laplace_decode(dec, ex_q8, k, str) laplace_decode_(dec, ex_q8, k, str)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) laplace_decode_vector_(dec, y, n, k, curr, means, str)
#else
# define laplace_decode_special(dec, decay, max, str) laplace_decode_special_(dec, decay, max)
# define laplace_decode(dec, ex_q8, k, str) laplace_decode_(dec, ex_q8, k)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) laplace_decode_vector_(dec, y, n, k, curr, means)
#endif

int laplace_decode_special_(od_ec_dec *dec, unsigned decay, int max OD_ACC_STR);
int laplace_decode_(od_ec_dec *dec, unsigned ex_q8, int k OD_ACC_STR);
void laplace_decode_vector_(od_ec_dec *dec, od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means
                                  OD_ACC_STR);


void od_pvq_decode(daala_dec_ctx *dec, od_coeff *ref, od_coeff *out, int q0,
 int pli, int bs, const double *beta, int robust, int is_keyframe,
 unsigned int *flags, int block_skip, const int16_t *qm,
 const int16_t *qm_inv);

#endif
