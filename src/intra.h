/*Daala video codec
Copyright (c) 2012-2013 Daala project contributors.  All rights reserved.

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

#if !defined(_intra_H)
# define _intra_H (1)
# include "filter.h"

# define OD_INTRA_NMODES (10)

# define OD_INTRA_NCONTEXTS (8)

typedef void (*od_intra_mult_func)(double *_p,int _pred_stride,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);

extern const od_intra_mult_func OD_INTRA_MULT[OD_NBSIZES];

extern const double OD_INTRA_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][2*4][2*4];
extern const unsigned char OD_INTRA_PRED_PROB_4x4[3]
 [OD_INTRA_NMODES][OD_INTRA_NCONTEXTS];

void od_intra_pred4x4_mult(double *_p,int _pred_stride,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);
void od_intra_pred8x8_mult(double *_p,int _pred_stride,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);
void od_intra_pred16x16_mult(double *_p,int _pred_stride,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);

/*Fetches intra prediction to a 4x4 block of coefficients at _c, using
   UR, UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    original coefficients (before quantization).
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
   The _mode parameter contains the target mode.
  On output:
   {{_out[0],_out[1],_out[2],_out[3]},{_out[4],_out[4+1],...}} contains
    the input coefficients with the prediction subtracted.*/
void od_intra_pred4x4_get(od_coeff *_out,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);
void od_intra_pred8x8_get(od_coeff *_out,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);
void od_intra_pred16x16_get(od_coeff *_out,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _mode);

void od_intra_pred4x4_dist(ogg_uint32_t *_dist,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _pli);
void od_intra_pred8x8_dist(ogg_uint32_t *_dist,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _pli);
void od_intra_pred16x16_dist(ogg_uint32_t *_dist,
 const od_coeff *_c,int _stride,
 const od_coeff *_ur,int _ur_stride,int _pli);

extern const signed char OD_INTRA_CHROMA_WEIGHTS_Q6[OD_INTRA_NMODES][3];

void od_chroma_pred4x4(od_coeff *_p,const od_coeff *_c,
 const od_coeff *_l,int _stride,const int _weights_q8[3]);

ogg_uint32_t od_chroma_pred4x4_dist(const od_coeff *_c,
 const od_coeff *_l,int _stride,const int _weights_q8[3],int _pli);

void od_intra_pred_cdf(ogg_uint16_t _cdf[],
 const unsigned char _probs[][OD_INTRA_NCONTEXTS],const ogg_uint16_t _p0[],
 int _nmodes,int _left,int _upleft,int _up);

int od_intra_pred_search(const ogg_uint16_t _cdf[],
 const ogg_uint32_t _dist[],int _nmodes,ogg_uint16_t _lambda);

void od_intra_pred_update(ogg_uint16_t _p0[],int _nmodes,int _mode,
 int _left,int _upleft,int _up);

void od_resample_luma_coeffs(od_coeff *l, int lstride,
 const od_coeff *c, int cstride, int xdec, int ydec, int n);

extern double OD_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][4*4*4];
extern int OD_PRED_MULTS_4x4[OD_INTRA_NMODES][4][4];
extern int OD_PRED_PARAMX_4x4[OD_INTRA_NMODES][4][4][4*4*4];
extern int OD_PRED_PARAMY_4x4[OD_INTRA_NMODES][4][4][4*4*4];
extern double OD_PRED_WEIGHTS_8x8[OD_INTRA_NMODES][8][8][4*8*8];
extern int OD_PRED_MULTS_8x8[OD_INTRA_NMODES][8][8];
extern int OD_PRED_PARAMX_8x8[OD_INTRA_NMODES][8][8][4*8*8];
extern int OD_PRED_PARAMY_8x8[OD_INTRA_NMODES][8][8][4*8*8];
extern double OD_PRED_WEIGHTS_16x16[OD_INTRA_NMODES][16][16][4*16*16];
extern int OD_PRED_MULTS_16x16[OD_INTRA_NMODES][16][16];
extern int OD_PRED_PARAMX_16x16[OD_INTRA_NMODES][16][16][4*16*16];
extern int OD_PRED_PARAMY_16x16[OD_INTRA_NMODES][16][16][4*16*16];


#endif
