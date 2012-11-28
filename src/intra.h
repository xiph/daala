/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
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

extern const double OD_INTRA_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][2*4][2*4];
extern const unsigned char OD_INTRA_PRED_PROB_4x4[3][OD_INTRA_NMODES][OD_INTRA_NCONTEXTS];

void od_intra_pred4x4_mult(const od_coeff *_c,int _stride,int _mode,double *_p);

/*Applies intra prediction to a 4x4 block of coefficients at _c, using
   UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    original coefficients (before quantization).
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
  On output:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    the input coefficients with the prediction subtracted.
  Return: The intra prediction mode used (0...OD_INTRA_NMODES-1).*/
int od_intra_pred4x4_apply(od_coeff *_c,int _stride);

/*Fetches intra prediction to a 4x4 block of coefficients at _c, using
   UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    original coefficients (before quantization).
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
   The _mode parameter contains the target mode.
  On output:
   {{_out[0],_out[1],_out[2],_out[3]},{_out[4],_out[4+1],...}} contains
    the input coefficients with the prediction subtracted.*/
void od_intra_pred4x4_get(od_coeff *_out,od_coeff *_c,int _stride, int _mode);

void od_intra_pred4x4_dist(od_coeff *_dist, od_coeff *_c,int _stride, int _pli);

/*Unapplies intra prediction to a 4x4 block of coefficients at _c, using
   UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    unquantized coefficients.
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
  On output:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    the input coefficients with the prediction added.
  Return: The intra prediction mode used (0...OD_INTRA_NMODES-1).*/
void od_intra_pred4x4_unapply(od_coeff *_c,int _stride,int _mode);

void od_intra_pred_cdf(ogg_uint16_t cdf[OD_INTRA_NMODES],
    const unsigned char probs[OD_INTRA_NMODES][OD_INTRA_NCONTEXTS],
    const ogg_uint16_t p0[OD_INTRA_NMODES],int left,int upleft,int up);

int od_intra_pred_search(ogg_uint16_t p0[OD_INTRA_NMODES],
    const ogg_uint16_t cdf[OD_INTRA_NMODES],
    const od_coeff dist[OD_INTRA_NMODES], ogg_uint16_t lambda, int left,
    int upleft, int up);

#endif
