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

#include "adapt.h"

/* Once per frame */
void od_adapt_row_init(od_adapt_row_ctx *row)
{
  int i,r;
  for(r=0;r<row->nhsb;r++)
    for(i=0;i<NB_ADAPT_CTX;i++)
      row->ctx[r].mean[i] = 2*od_adapt_params[i][0]-(od_adapt_params[i][0]>>od_adapt_params[i][1]);
}

/* Once per macroblock line */
void od_adapt_hmean_init(od_adapt_ctx *hmean)
{
  int i;
  for(i=0;i<NB_ADAPT_CTX;i++)
    hmean->mean[i] = od_adapt_params[i][0];
}

/* Compute stats for current mb */
void od_adapt_update_stats(const od_adapt_row_ctx *row, int sbx,
    const od_adapt_ctx *hmean, od_adapt_ctx *curr)
{
  int i;
  for(i=0;i<NB_ADAPT_CTX;i++)
    curr->mean[i] = row->ctx[sbx].mean[i] + (hmean->mean[i]>>od_adapt_params[i][1]);
}

void od_adapt_sb(od_adapt_row_ctx *row, int sbx, od_adapt_ctx *hmean, const od_adapt_ctx *curr)
{
  int i;
  for (i=0;i<NB_ADAPT_CTX;i++)
  {
    if (curr->curr[i]!=OD_ADAPT_NO_VALUE)
    {
      row->ctx[sbx].curr[i] = curr->curr[i]>>1;
      hmean->mean[i] += (row->ctx[sbx].curr[i] - hmean->mean[i])>>od_adapt_params[i][1];
      row->ctx[sbx].mean[i] += (hmean->mean[i]-row->ctx[sbx].mean[i])>>od_adapt_params[i][1];
    } else {
      row->ctx[sbx].curr[i] = OD_ADAPT_NO_VALUE;
    }
  }
}

void od_adapt_row(od_adapt_row_ctx *row, od_adapt_ctx *hmean)
{
  int r;
  od_adapt_hmean_init(hmean);
  for(r=row->nhsb;r-->0;)
  {
    int i;
    for(i=0;i<NB_ADAPT_CTX;i++)
    {
      if (row->ctx[r].curr[i]!=OD_ADAPT_NO_VALUE)
      {
        row->ctx[r].mean[i] += (hmean->mean[i]>>od_adapt_params[i][1])
                             -(hmean->mean[i]>>2*od_adapt_params[i][1]);
        hmean->mean[i] += (row->ctx[r].curr[i]-hmean->mean[i])>>od_adapt_params[i][1];
      }
    }
  }
}
