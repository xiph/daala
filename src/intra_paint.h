/*Daala video codec
Copyright (c) 2014 Daala project contributors. All rights reserved.

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

#if !defined(_intra_paint_H)
# define _intra_paint_H (1)
# include "internal.h"

#include "dct.h"
#include "odintrin.h"
#include "entenc.h"
#include "state.h"

extern int dir_prob[9];
extern int dc_prob[8];


/* Warning, this will fail for images larger than 2048 x 2048 */
#define MAX_VAR_BLOCKS 1024

void od_paint_dering(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_sum2, int *edge_count, int q);

void od_intra_paint_analysis(const unsigned char *paint, int stride,
 const unsigned char *dec8, int bstride, unsigned char *mode, int mstride,
 int *edge_sum, int *edge_sum2, int *edge_count, int bx, int by, int level);

void od_paint_compute_edge_mask(od_adapt_ctx *adapt,
 unsigned char *paint, const unsigned char *img, unsigned char *paint_mask,
 int stride, const unsigned char *dec8, int bstride, unsigned char *mode,
 int mstride, int *edge_sum1, int *edge_sum2, int *edge_count, int q, int bx,
 int by, int level);

void od_paint_dering_decode(od_adapt_ctx *adapt, od_ec_dec *dec, unsigned char *paint,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_sum2, int *edge_count, int q);

extern unsigned char *paint_mask;
extern unsigned char *paint_out;

#endif
