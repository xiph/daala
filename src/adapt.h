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

#ifndef _ADAPT_H
#define _ADAPT_H

#include <ogg/ogg.h>

/* One line per adapted parameter, first column is the init value,
   second column is the adaptation speed */
static const ogg_int32_t od_adapt_params[][2] = {
#define OD_ADAPT_K_Q8        0
    {2031, 2},
#define OD_ADAPT_SUM_EX_Q8   1
    {216, 2},
#define OD_ADAPT_COUNT_Q8    2
    {104, 1},
#define OD_ADAPT_COUNT_EX_Q8 3
    {128, 1},
};

#define OD_ADAPT_NO_VALUE (-2147483648)

#define NB_ADAPT_CTX ((int)sizeof(od_adapt_params)/(int)sizeof(od_adapt_params[0]))

typedef struct{
  ogg_int32_t curr[NB_ADAPT_CTX];
  ogg_int32_t mean[NB_ADAPT_CTX];
} od_adapt_ctx;

typedef struct{
  od_adapt_ctx *ctx;
  int nhmbs;
} od_adapt_row_ctx;

/* Once per frame */
void od_adapt_row_init(od_adapt_row_ctx *row);

/* Once per macroblock line */
void od_adapt_hmean_init(od_adapt_ctx *hmean);

/* Compute stats for current mb */
void od_adapt_update_stats(const od_adapt_row_ctx *curr_row, int id, const od_adapt_ctx *hmean, od_adapt_ctx *out);

void od_adapt_mb(od_adapt_row_ctx *row, int id, od_adapt_ctx *hmean, const od_adapt_ctx *curr);

void od_adapt_row(od_adapt_row_ctx *row, od_adapt_ctx *hmean);

#endif
