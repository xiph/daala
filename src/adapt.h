/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_adapt_H)
# define _adapt_H (1)

# include "../include/daala/daala_integer.h"

# define OD_NSB_ADAPT_CTXS (4)

# define OD_ADAPT_K_Q8        0
# define OD_ADAPT_SUM_EX_Q8   1
# define OD_ADAPT_COUNT_Q8    2
# define OD_ADAPT_COUNT_EX_Q8 3

# define OD_ADAPT_NO_VALUE (-2147483647-1)
# define OD_NADAPT_CTXS_MAX (OD_NSB_ADAPT_CTXS)

typedef struct {
  int32_t curr;
  int32_t mean;
} od_adapt2d_data;

typedef struct {
  od_adapt2d_data *data;
  /*Adaptation rates and initial values.*/
  const int32_t *params;
  /*Number of values in a row.*/
  int nhv;
  /*Number of contexts being modeled.*/
  int nctx;
} od_adapt2d_ctx;

void od_adapt2d_init(od_adapt2d_ctx *ctx, int nhv,
 int nctx, const int32_t *params);
void od_adapt2d_clear(od_adapt2d_ctx *ctx);

void od_adapt2d_row_init(od_adapt2d_ctx *ctx);

void od_adapt2d_hmean_init(const od_adapt2d_ctx *ctx, int32_t *hmean);

void od_adapt2d_get_stats(const od_adapt2d_ctx *ctx, int xpos,
 const int32_t *hmean, int32_t *means);

void od_adapt2d_forward(od_adapt2d_ctx *ctx, int xpos, int32_t *hmean,
 const int32_t *curr);

void od_adapt2d_row_backward(od_adapt2d_ctx *ctx);

#endif
