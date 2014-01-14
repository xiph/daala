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

#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"

# include "filter.h"
# include "pvq_code.h"

# define USE_PSEUDO_ZIGZAG (1)

typedef unsigned char index_pair[2];

typedef struct {
  const index_pair *const dst_table;
  int size;
  int nb_bands;
  const int *const band_offsets;
} band_layout;

extern const band_layout od_layout4;
extern const band_layout od_layout8;
extern const band_layout od_layout16;

#define PRED4_PROB (26376)
extern const ogg_uint16_t pred8_cdf[16];

void od_bands_from_raster(const band_layout *layout, od_coeff *dst,
 od_coeff *src, int stride);
void od_raster_from_bands(const band_layout *layout, od_coeff *dst,
 int stride, od_coeff *src);
void od_band_pseudo_zigzag(od_coeff *dst,  int n, od_coeff *src, int stride,
 int interleave);
void od_band_pseudo_dezigzag(od_coeff *dst,  int stride, od_coeff *src,
 int n, int interleave);

int quant_pvq_theta(ogg_int32_t *x0, const ogg_int32_t *r0,
 ogg_int16_t *scale0, int *y, int n, int q0, int *qg, int shift, int intra);

int pvq_unquant_k(const ogg_int32_t *_r, int _n, int _qg, int _scale,
 int shift, int intra);

int quant_pvq(od_coeff *_x, const od_coeff *_r, ogg_int16_t *_scale,
 int *y, int N, int Q, int *qg, int shift, int intra);

void dequant_pvq(od_coeff *_x, const od_coeff *_r, ogg_int16_t *_scale,
 int N, int _Q, int qg, int shift, int intra);

int quant_scalar(ogg_int32_t *_x, const ogg_int32_t *_r,
 ogg_int16_t *_scale, int *y, int N, int Q, ogg_int32_t *_adapt);

int quant_pvq_noref(ogg_int32_t *_x, float gr,
 ogg_int16_t *_scale, int *y, int N, int Q);


/* New PVQ implementation */

int pvq_theta(od_coeff *x0, od_coeff *r0, int n, int q0, od_coeff *y, int *itheta,
 int *max_theta, int *vk);

int compute_householder(double *r, int n, double gr, int *sign);

void pvq_synthesis(od_coeff *x0, od_coeff *y, const double *r, int n, int noref,
 int qg, double gain_offset, double theta, int m, int s, double q);

int od_compute_max_theta(const od_coeff *r, int n, int q, double *gr,
 double *qcg, int *qg, double *gain_offset, int noref);

double od_compute_k_theta(int *k, double qcg, int t, int max_theta, int noref,
 int n);

#endif
