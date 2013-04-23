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

#if !defined(_tools_tran_tools_H)
# define _tools_tran_tools_H (0)

#include <stdio.h>
#include "intra_fit_tools.h"
#include "../src/internal.h"

#define NE_DISABLE_FILTER (0)

typedef struct image_ctx image_ctx;

struct image_ctx{
  const char *name;
  int         nxblocks;
  int         nyblocks;
};

void image_ctx_init(image_ctx *_this,const char *_name,int _nxblocks,
 int _nyblocks);

typedef struct trans_data trans_data;

struct trans_data{
  int     sz;
  int     n;
  double *mean;
  double *cov;
  double *work;
};

void trans_data_init(trans_data *_this,int _sz);
void trans_data_clear(trans_data *_this);
void trans_data_add(trans_data *_this,const unsigned char *_data);
void trans_data_combine(trans_data *_a,const trans_data *_b);
void trans_data_correct(trans_data *_this);
void trans_data_normalize(trans_data *_this);
void trans_data_collapse(trans_data *_this,int _n,double *_r);
void trans_data_expand(trans_data *_this,int _n,const double *_r);
void trans_data_print(trans_data *_this,FILE *_fp);

void covariance_collapse(const double *_cov,int _sz,int _n,double *_r,
 double *_work);
void covariance_expand(double *_cov,int _sz,int _n,const double *_r);

typedef struct trans_ctx trans_ctx;

struct trans_ctx{
  image_ctx  img;
  trans_data td;
};

typedef struct param_data param_data;

struct param_data{
  int    p[B_SZ/2-1];
  int    q[B_SZ/2-1];
  int    s[B_SZ/2];
  double grgt[B_SZ*B_SZ];
  double hi[2*B_SZ*B_SZ];
  double cg;
  double sba;
  double w;
  double bo;
};

void auto_regressive_collapsed(double *_out,int _sz,int _n,double _r);

void analysis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n,const int *_f);
void synthesis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n,const int *_f);

double coding_gain_1d(const double _r[2*B_SZ*2*B_SZ],const int *_f);
double coding_gain_1d_collapsed(const double _r[2*B_SZ],const int *_f);
double coding_gain_2d(const double _r[2*B_SZ*2*B_SZ*2*B_SZ*2*B_SZ],const int *_f);
double coding_gain_2d_collapsed(const double _r[2*B_SZ*2*B_SZ],const int *_f);

extern const double *SUBSET1_1D[OD_NBSIZES];
extern const double *SUBSET3_1D[OD_NBSIZES];

extern const double *SUBSET1_2D[OD_NBSIZES];
extern const double *SUBSET3_2D[OD_NBSIZES];

#endif
