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

#if !defined(_od_covmat_H)
# define _od_covmat_H (0)

#include <stdio.h>

typedef struct od_covmat od_covmat;

struct od_covmat{
  int     sz;
  double  w;
  double *mean;
  double *cov;
  double *work;
};

void od_covmat_init(od_covmat *_this,int _sz);
void od_covmat_clear(od_covmat *_this);
void od_covmat_reset(od_covmat *_this);
void od_covmat_add(od_covmat *_this,const double *_data,double _w);
void od_covmat_combine(od_covmat *_a,const od_covmat *_b);
void od_covmat_update(od_covmat *_this,const double *_cov,const double *_mean,
 double _w);
void od_covmat_correct(od_covmat *_this);
void od_covmat_normalize(od_covmat *_this);
void od_covmat_collapse(od_covmat *_this,int _n,double *_r);
void od_covmat_expand(od_covmat *_this,int _n,const double *_r);
void od_covmat_print(od_covmat *_this,FILE *_fp);

#endif
