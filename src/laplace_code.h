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

#if !defined(_laplace_code_H)
# define _laplace_code_H

# include "entenc.h"
# include "entdec.h"
# include "filter.h"
# include "adapt.h"

#if OD_ACCOUNTING
# define laplace_decode_special(dec, decay, max, str) laplace_decode_special_(dec, decay, max, str)
# define laplace_decode(dec, ex_q8, k, str) laplace_decode_(dec, ex_q8, k, str)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) laplace_decode_vector_(dec, y, n, k, curr, means, str)
#else
# define laplace_decode_special(dec, decay, max, str) laplace_decode_special_(dec, decay, max)
# define laplace_decode(dec, ex_q8, k, str) laplace_decode_(dec, ex_q8, k)
#define laplace_decode_vector(dec, y, n, k, curr, means, str) laplace_decode_vector_(dec, y, n, k, curr, means)
#endif

extern const uint16_t EXP_CDF_TABLE[][16];
extern const uint16_t LAPLACE_OFFSET[];

void laplace_encode_special(od_ec_enc *enc, int x, unsigned decay, int max);
void laplace_encode(od_ec_enc *enc, int x, int ex_q8, int k);
void laplace_encode_vector(od_ec_enc *enc, const od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means);

int laplace_decode_special_(od_ec_dec *dec, unsigned decay, int max OD_ACC_STR);
int laplace_decode_(od_ec_dec *dec, int ex_q8, int k OD_ACC_STR);
void laplace_decode_vector_(od_ec_dec *dec, od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means
                                  OD_ACC_STR);

#endif
