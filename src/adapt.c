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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>

#include "internal.h"
#include "adapt.h"

/*This implements a 2D moving average filter.
  Every point is a combination of a forward horizontal moving average, a
   vertical moving average, and a backwards moving average on the prior row.
  This gives a computationally efficient casual 2d filter with a
   symmetrical response.
  The filter supports adapting multiple distinct values at once for better
   memory locality during the backwards pass.*/

/*Initializes a 2D moving average filter.*/
void od_adapt_init(od_adapt_ctx *ctx, int nhv,
 int nctx, const ogg_int32_t *params) {
  OD_ASSERT(nctx <= OD_NADAPT_CTXS_MAX);
  ctx->data = (od_adapt_data *)_ogg_malloc(sizeof(*ctx->data)*nhv*nctx);
  ctx->nhv = nhv;
  ctx->nctx = nctx;
  ctx->params = params;
}

/*Frees a 2D moving average filter.*/
void od_adapt_clear(od_adapt_ctx *ctx) {
  _ogg_free(ctx->data);
  ctx->data = NULL;
}

/*Reinitializes a row context.
  This must be called once per-frame per filter.*/
void od_adapt_row_init(od_adapt_ctx *ctx) {
  const ogg_int32_t *params;
  od_adapt_data *data;
  int nhv;
  int nctx;
  int r;
  nhv = ctx->nhv;
  nctx = ctx->nctx;
  data = ctx->data;
  params = ctx->params;
  for (r = 0; r < nhv; r++) {
    int i;
    for (i = 0; i < nctx; i++) {
      data[r*nctx + i].mean = 2*params[nctx + i]
       - (params[nctx + i] >> params[i]);
    }
  }
}

/*Initializes a set of running average means for the 2D moving average.*/
void od_adapt_hmean_init(const od_adapt_ctx *ctx, ogg_int32_t *hmean) {
  const ogg_int32_t *inits;
  int nctx;
  int i;
  nctx = ctx->nctx;
  inits = &ctx->params[nctx];
  for (i = 0; i < nctx; i++) {
    hmean[i] = inits[i];
  }
}

/*Compute expected values for the current position.
  xpos: current x offset.
  hmeans: horizontal running averages.
  means: current values.*/
void od_adapt_get_stats(const od_adapt_ctx *ctx, int xpos,
 const ogg_int32_t *hmean, ogg_int32_t *means) {
  const ogg_int32_t *adapts;
  od_adapt_data *data;
  int i;
  int nctx;
  nctx = ctx->nctx;
  data = &ctx->data[xpos*nctx];
  adapts = ctx->params;
  for (i = 0; i < nctx; i++) {
    means[i] = data[i].mean + (hmean[i] >> adapts[i]);
  }
}

/*Update the filter for the current position.
  xpos: current x offset.
  hmeans: horizontal running averages.
  curr: current values.*/
void od_adapt_forward(od_adapt_ctx *ctx, int xpos,
 ogg_int32_t *hmean, const ogg_int32_t *curr) {
  od_adapt_data *data;
  const ogg_int32_t *adapts;
  int i;
  int nctx;
  nctx = ctx->nctx;
  data = &ctx->data[xpos*nctx];
  adapts = ctx->params;
  for (i = 0; i < nctx; i++) {
    if (curr[i] != OD_ADAPT_NO_VALUE) {
      ogg_int32_t adapt;
      adapt = adapts[i];
      data[i].curr = curr[i] >> 1;
      hmean[i] += (data[i].curr - hmean[i]) >> adapt;
      data[i].mean += (hmean[i] - data[i].mean) >> adapt;
    }
    else {
      data[i].curr = OD_ADAPT_NO_VALUE;
    }
  }
}

/*Update the 2D filter for a finished row.
  This must be run at the end of every horizontal scan.*/
void od_adapt_row_backward(od_adapt_ctx *ctx) {
  const ogg_int32_t *adapts;
  od_adapt_data *data;
  ogg_int32_t hmean[OD_NADAPT_CTXS_MAX];
  int r;
  int nctx;
  int nhv;
  od_adapt_hmean_init(ctx, hmean);
  nctx = ctx->nctx;
  nhv = ctx->nhv;
  data = ctx->data;
  adapts = ctx->params;
  for (r = nhv; r-- > 0;) {
    int i;
    int roff;
    roff = r*nctx;
    for (i = 0; i < nctx; i++) {
      if (data[roff + i].curr != OD_ADAPT_NO_VALUE) {
        ogg_int32_t adapt;
        adapt = adapts[i];
        data[roff + i].mean += (hmean[i] >> adapt) - (hmean[i] >> 2*adapt);
        hmean[i] += (data[roff + i].curr - hmean[i]) >> adapt;
      }
    }
  }
}
