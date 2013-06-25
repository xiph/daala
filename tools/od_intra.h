/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_tools_od_intra_H)
# define _tools_od_intra_H (0)

#include <stdio.h>
#include "od_defs.h"

typedef void (*ne_intra_mult_func)(double *_pred,int _stride,
 const od_coeff *_coeff,int _mode);

extern const ne_intra_mult_func NE_INTRA_MULT[OD_NBSIZES];

void od_intra_init();
void od_intra_clear();

extern double NE_PRED_OFFSETS_4x4[OD_INTRA_NMODES][4][4];
extern int NE_PRED_MULTS_4x4[OD_INTRA_NMODES][4][4];
extern double *NE_PRED_WEIGHTS_4x4[OD_INTRA_NMODES];
extern int *NE_PRED_INDEX_4x4[OD_INTRA_NMODES];

extern double NE_PRED_OFFSETS_8x8[OD_INTRA_NMODES][8][8];
extern int NE_PRED_MULTS_8x8[OD_INTRA_NMODES][8][8];
extern double *NE_PRED_WEIGHTS_8x8[OD_INTRA_NMODES];
extern int *NE_PRED_INDEX_8x8[OD_INTRA_NMODES];

extern double NE_PRED_OFFSETS_16x16[OD_INTRA_NMODES][16][16];
extern int NE_PRED_MULTS_16x16[OD_INTRA_NMODES][16][16];
extern double *NE_PRED_WEIGHTS_16x16[OD_INTRA_NMODES];
extern int *NE_PRED_INDEX_16x16[OD_INTRA_NMODES];

void update_predictors(int _mode,double *_beta_0,double *_beta_1,int *_mask);
void print_predictors(FILE *_fp);

#endif
