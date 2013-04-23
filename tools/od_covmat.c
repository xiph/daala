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

#include <math.h>
#include <stdlib.h>
#include "od_defs.h"
#include "od_covmat.h"

#define PRINT_COLLAPSE (0)

void od_covmat_init(od_covmat *_this,int _sz){
  _this->sz=_sz;
  _this->mean=(double *)malloc(sizeof(*_this->mean)*_this->sz);
  _this->cov=(double *)malloc(sizeof(*_this->cov)*_this->sz*_this->sz);
  _this->work=(double *)malloc(sizeof(*_this->work)*_this->sz);
  od_covmat_reset(_this);
}

void od_covmat_clear(od_covmat *_this){
  free(_this->mean);
  free(_this->cov);
  free(_this->work);
}

void od_covmat_reset(od_covmat *_this){
  int j;
  int i;
  _this->w=0;
  for(j=0;j<_this->sz;j++){
    _this->mean[j]=0;
    for(i=0;i<_this->sz;i++){
      _this->cov[_this->sz*j+i]=0;
    }
  }
}

void cov_update(double *_mean,double *_cov,double *_work,double *_weight,int _n,
 const double *_data,double _w){
  double s;
  int j;
  int i;
  (*_weight)+=_w;
  s=1/(*_weight);
  for(i=0;i<_n;i++){
    _work[i]=_data[i]-_mean[i];
#if FAST_MATH
    _mean[i]+=_work[i]*s;
#else
    _mean[i]+=_work[i]/_w;
#endif
  }
  s=((*_weight)-_w)/(*_weight);
  for(j=0;j<_n;j++){
    for(i=0;i<_n;i++){
#if FAST_MATH
      _cov[_n*j+i]+=_work[j]*_work[i]*s;
#else
      _cov[_n*j+i]+=_work[j]*_work[i]*((*_weight)-_w)/(*_weight);
#endif
    }
  }
}

void od_covmat_add(od_covmat *_this,const double *_data,double _w){
/* why does this cause performance penalty?? 15.8% at -03 and 8.7% at -O2 */
#if 0
  cov_update(_this->mean,_this->cov,_this->work,&_this->w,_this->sz,_data,_w);
#else
  double s;
  int    i;
  int    j;
  _this->w+=_w;
  s=1/_this->w;
  for(i=0;i<_this->sz;i++){
    _this->work[i]=_data[i]-_this->mean[i];
#if FAST_MATH
    _this->mean[i]+=_this->work[i]*s;
#else
    _this->mean[i]+=_this->work[i]/_this->w;
#endif
  }
  s=(_this->w-_w)/_this->w;
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
#if FAST_MATH
      _this->cov[_this->sz*j+i]+=_this->work[j]*_this->work[i]*s;
#else
      _this->cov[_this->sz*j+i]+=
       _this->work[j]*_this->work[i]*(_this->w-_w)/_this->w;
#endif
    }
  }
#endif
}

void od_covmat_combine(od_covmat *_a,const od_covmat *_b){
  double s;
  int    i;
  int    j;
  if(_b->w==0){
    return;
  }
  s=((double)_b->w)/(_a->w+_b->w);
  for(i=0;i<_a->sz;i++){
    _a->work[i]=_b->mean[i]-_a->mean[i];
#if FAST_MATH
    _a->mean[i]+=_a->work[i]*s;
#else
    _a->mean[i]+=_a->work[i]*_b->w/(_a->w+_b->w);
#endif
  }
  s*=_a->w;
  for(i=0;i<_a->sz;i++){
    for(j=0;j<_a->sz;j++){
#if FAST_MATH
      _a->cov[_a->sz*i+j]+=_b->cov[_a->sz*i+j]+_a->work[i]*_a->work[j]*s;
#else
      _a->cov[_a->sz*i+j]+=
       _b->cov[_a->sz*i+j]+_a->work[i]*_a->work[j]*_a->w*_b->w/(_a->w+_b->w);
#endif
    }
  }
  _a->w+=_b->w;
}

void od_covmat_correct(od_covmat *_this){
  double s;
  int    j;
  int    i;
  s=1/_this->w;
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
#if FAST_MATH
      _this->cov[_this->sz*j+i]*=s;
#else
      _this->cov[_this->sz*j+i]/=_this->w;
#endif
    }
  }
}

void od_covmat_normalize(od_covmat *_this){
  int i;
  int j;
  for(i=0;i<_this->sz;i++){
    _this->work[i]=sqrt(_this->cov[_this->sz*i+i]);
  }
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
      _this->cov[_this->sz*j+i]/=_this->work[j]*_this->work[i];
    }
  }
}

static void covariance_collapse(const double *_cov,int _sz,int _n,double *_r,
 double *_work){
  int m;
  int i;
  int j;
  int u;
  int v;
  m=_sz/_n;
  for(i=0;i<_sz;i++){
    _r[i]=0;
    _work[i]=0;
  }
  for(v=0;v<_n;v++){
    for(u=v;u<_n;u++){
      for(j=0;j<m;j++){
        for(i=j;i<m;i++){
          _r[(u-v)*m+(i-j)]+=_cov[(v*m+j)*_sz+(u*m+i)];
          _r[(u-v)*m+(i-j)]+=_cov[(v*m+i)*_sz+(u*m+j)];
          _r[(u-v)*m+(i-j)]+=_cov[(u*m+j)*_sz+(v*m+i)];
          _r[(u-v)*m+(i-j)]+=_cov[(u*m+i)*_sz+(v*m+j)];
          _work[(u-v)*m+(i-j)]+=4;
        }
      }
    }
  }
  for(i=0;i<_sz;i++){
    _r[i]/=_work[i];
  }
}

void od_covmat_collapse(od_covmat *_this,int _n,double *_r){
  covariance_collapse(_this->cov,_this->sz,_n,_r,_this->work);
#if PRINT_COLLAPSE
  for(i=0;i<_this->sz;i++){
    fprintf(stderr,"%s  %- 24.18G",i>0?",":"",_r[i]);
  }
  fprintf(stderr,"\n");
#endif
}

static void covariance_expand(double *_cov,int _sz,int _n,const double *_r){
  int i;
  int j;
  int u;
  int v;
  int m;
  m=_sz/_n;
  for(v=0;v<_n;v++){
    for(u=0;u<_n;u++){
      for(j=0;j<m;j++){
        for(i=0;i<m;i++){
          _cov[(v*m+j)*_sz+u*m+i]=_r[abs(v-u)*m+abs(i-j)];
        }
      }
    }
  }
}

void od_covmat_expand(od_covmat *_this,int _n,const double *_r){
  covariance_expand(_this->cov,_this->sz,_n,_r);
}

void od_covmat_print(od_covmat *_this,FILE *_fp){
  int i;
  int j;
  for(i=0;i<_this->sz;i++){
    fprintf(_fp,"%s  %- 24.18G",i>0?",":"",_this->mean[i]);
  }
  fprintf(_fp,"\n");
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
      fprintf(_fp,"%s  %- 24.18G",i>0?",":"",_this->cov[_this->sz*j+i]);
    }
    fprintf(_fp,"\n");
  }
}
