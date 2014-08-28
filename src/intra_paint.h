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


#define QUANTIZE (1)

/* Warning, this will fail for images larger than 2048 x 2048 */
#define MAX_VAR_BLOCKS 1024
#define MAXN 64

void pixel_interp(int pi[4], int pj[4], int w[4], int m, int i, int j, int ln);

void interp_block(unsigned char *img, const unsigned char *edge_accum, int n,
 int stride, int m);

void predict_bottom_edge(int *p, unsigned char *edge_accum, int n, int stride,
 int m, int has_right);

void predict_right_edge(int *p, unsigned char *edge_accum, int n, int stride,
 int m, int has_bottom);

extern const od_fdct_func_1d my_fdct_table[5];
extern const od_idct_func_1d my_idct_table[5];

int od_intra_paint_mode_cdf(ogg_uint16_t *cdf, int *dir_list, int *prob_list,
 int *ctx_list, int *dc_ctx, const unsigned char *mode, int bx, int by,
 int ln, int mstride, const unsigned char *dec8, int bstride, int res);

int compute_intra_paint(od_ec_enc *enc, unsigned char *img, int w, int h,
 int stride);

void od_intra_paint_encode(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res);

void od_intra_paint_encode2(od_adapt_ctx *adapt, od_ec_enc *enc,
 unsigned char *paint, const unsigned char *img, int stride, int nhsb, int nvsb,
 unsigned char *dec8, unsigned char *mode4, unsigned char *mode8,
 unsigned char *mode16, unsigned char *mode32, int *edge_sum, int *edge_count,
 int q, int res);

void od_paint_dering(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_count, int q,
 int res);

#endif
