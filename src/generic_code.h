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

#if !defined(_generic_code_H)
# define _generic_code_H

# include "entenc.h"
# include "entdec.h"

# define GENERIC_TABLES 12

typedef struct {
  /** cdf for multiple expectations of x */
  uint16_t cdf[GENERIC_TABLES][16];
  /** Frequency increment for learning the cdfs */
  int increment;
} generic_encoder;

#define OD_IIR_DIADIC(y, x, shift) ((y) += ((x) - (y)) >> (shift))

void generic_model_init(generic_encoder *model);

#define OD_CDFS_INIT(cdf, val) od_cdf_init(&cdf[0][0],\
 sizeof(cdf)/sizeof(cdf[0]), sizeof(cdf[0])/sizeof(cdf[0][0]), val, val)

#define OD_CDFS_INIT_FIRST(cdf, val, first) od_cdf_init(&cdf[0][0],\
 sizeof(cdf)/sizeof(cdf[0]), sizeof(cdf[0])/sizeof(cdf[0][0]), val, first)

#define OD_SINGLE_CDF_INIT(cdf, val) od_cdf_init(cdf,\
 1, sizeof(cdf)/sizeof(cdf[0]), val, val)

#define OD_SINGLE_CDF_INIT_FIRST(cdf, val, first) od_cdf_init(cdf,\
 1, sizeof(cdf)/sizeof(cdf[0]), val, first)

void od_cdf_init(uint16_t *cdf, int ncdfs, int nsyms, int val, int first);

void od_encode_cdf_adapt(od_ec_enc *ec, int val, uint16_t *cdf, int n,
 int increment);

int od_decode_cdf_adapt(od_ec_dec *ec, uint16_t *cdf, int n,
 int increment);

void generic_encode(od_ec_enc *enc, generic_encoder *model, int x, int max,
 int *ex_q16, int integration);
double generic_encode_cost(generic_encoder *model, int x, int max,
 int *ex_q16);

int generic_decode(od_ec_dec *dec, generic_encoder *model, int max,
 int *ex_q16, int integration);

int log_ex(int ex_q16);

void generic_model_update(generic_encoder *model, int *ex_q16, int x, int xs,
 int id, int integration);

#endif
