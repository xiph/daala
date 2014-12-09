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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include "od_defs.h"
#include "od_filter.h"
#include "trans_tools.h"

#define PRINT_COLLAPSE (0)
#define PRINT_CG_MATH (0)

void image_ctx_init(image_ctx *_this,const char *_name,int _nxblocks,
 int _nyblocks){
  _this->name=_name;
  _this->nxblocks=_nxblocks;
  _this->nyblocks=_nyblocks;
}

void trans_data_init(trans_data *_this,int _sz){
  int i;
  int j;
  _this->sz=_sz;
  _this->n=0;
  _this->mean=(double *)malloc(sizeof(*_this->mean)*_this->sz);
  _this->cov=(double *)malloc(sizeof(*_this->cov)*_this->sz*_this->sz);
  _this->work=(double *)malloc(sizeof(*_this->work)*_this->sz);
  for(j=0;j<_this->sz;j++){
    _this->mean[j]=0;
    for(i=0;i<_this->sz;i++){
      _this->cov[_this->sz*j+i]=0;
    }
  }
}

void trans_data_clear(trans_data *_this){
  free(_this->mean);
  free(_this->cov);
  free(_this->work);
}

void trans_data_add(trans_data *_this,const unsigned char *_data){
  double s;
  int    i;
  int    j;
  _this->n++;
  s=1.0/_this->n;
  for(i=0;i<_this->sz;i++){
    _this->work[i]=_data[i]-_this->mean[i];
#if FAST_MATH
    _this->mean[i]+=_this->work[i]*s;
#else
    _this->mean[i]+=_this->work[i]/_this->n;
#endif
  }
  s=(_this->n-1.0)/_this->n;
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
#if FAST_MATH
      _this->cov[_this->sz*j+i]+=_this->work[j]*_this->work[i]*s;
#else
      _this->cov[_this->sz*j+i]+=
       _this->work[j]*_this->work[i]*(_this->n-1)/_this->n;
#endif
    }
  }
}

void trans_data_combine(trans_data *_a,const trans_data *_b){
  double s;
  int    i;
  int    j;
  if(_b->n==0){
    return;
  }
  s=((double)_b->n)/(_a->n+_b->n);
  for(i=0;i<_a->sz;i++){
    _a->work[i]=_b->mean[i]-_a->mean[i];
#if FAST_MATH
    _a->mean[i]+=_a->work[i]*s;
#else
    _a->mean[i]+=_a->work[i]*_b->n/(_a->n+_b->n);
#endif
  }
  s*=_a->n;
  for(i=0;i<_a->sz;i++){
    for(j=0;j<_a->sz;j++){
#if FAST_MATH
      _a->cov[_a->sz*i+j]+=_b->cov[_a->sz*i+j]+_a->work[i]*_a->work[j]*s;
#else
      _a->cov[_a->sz*i+j]+=
       _b->cov[_a->sz*i+j]+_a->work[i]*_a->work[j]*_a->n*_b->n/(_a->n+_b->n);
#endif
    }
  }
  _a->n+=_b->n;
}

void trans_data_correct(trans_data *_this){
  double s;
  int    j;
  int    i;
  s=1.0/_this->n;
  for(j=0;j<_this->sz;j++){
    for(i=0;i<_this->sz;i++){
#if FAST_MATH
      _this->cov[_this->sz*j+i]*=s;
#else
      _this->cov[_this->sz*j+i]/=_this->n;
#endif
    }
  }
}

void trans_data_normalize(trans_data *_this){
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

void trans_data_collapse(trans_data *_this,int _n,double *_r){
  covariance_collapse(_this->cov,_this->sz,_n,_r,_this->work);
}

void trans_data_expand(trans_data *_this,int _n,const double *_r){
  covariance_expand(_this->cov,_this->sz,_n,_r);
}

void trans_data_print(trans_data *_this,FILE *_fp){
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

void covariance_collapse(const double *_cov,int _sz,int _n,double *_r,
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
#if PRINT_COLLAPSE
  for(i=0;i<_sz;i++){
    fprintf(stderr,"%s  %- 24.18G",i>0?",":"",_r[i]);
  }
  fprintf(stderr,"\n");
#endif
}

void covariance_expand(double *_cov,int _sz,int _n,const double *_r){
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

/*The true forward 4-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT4_BASIS[4][4]={
  {
     0.5,                                 0.5,
     0.5,                                 0.5
  },
  {
     0.65328148243818826392832158671359,  0.27059805007309849219986160268319,
    -0.27059805007309849219986160268319, -0.65328148243818826392832158671359
  },
  {
     0.5,                                -0.5,
    -0.5,                                 0.5
  },
  {
     0.27059805007309849219986160268319, -0.65328148243818826392832158671359,
     0.65328148243818826392832158671359, -0.27059805007309849219986160268319
  },
};

/*The true forward 8-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT8_BASIS[8][8]={
  {
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.49039264020161522456309111806712,  0.41573480615127261853939418880895,
     0.27778511650980111237141540697427,  0.097545161008064133924142434238511,
    -0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
    -0.41573480615127261853939418880895, -0.49039264020161522456309111806712
  },
  {
     0.46193976625574337806409159469839,  0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
     0.19134171618254488586422999201520,  0.46193976625574337806409159469839
  },
  {
     0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
    -0.49039264020161522456309111806712, -0.27778511650980111237141540697427,
     0.27778511650980111237141540697427,  0.49039264020161522456309111806712,
     0.097545161008064133924142434238511,-0.41573480615127261853939418880895
  },
  {
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.27778511650980111237141540697427, -0.49039264020161522456309111806712,
     0.097545161008064133924142434238511, 0.41573480615127261853939418880895,
    -0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
     0.49039264020161522456309111806712, -0.27778511650980111237141540697427
  },
  {
     0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
     0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520,  0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839,  0.19134171618254488586422999201520
  },
  {
     0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
     0.41573480615127261853939418880895, -0.49039264020161522456309111806712,
     0.49039264020161522456309111806712, -0.41573480615127261853939418880895,
     0.27778511650980111237141540697427, -0.097545161008064133924142434238511
  }
};

/*The true forward 8-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT16_BASIS[16][16]={
  {
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25
  },{
     0.35185093438159561475651955574599,  0.33832950029358816956728612239141,
     0.31180625324666780814299756569942,  0.27330046675043937205975812924953,
     0.22429189658565907106266034730521,  0.16666391461943662432490137779072,
     0.10263113188058934528909230943878,  0.034654292299772865648933749133244,
    -0.034654292299772865648933749133244,-0.10263113188058934528909230943878,
    -0.16666391461943662432490137779072, -0.22429189658565907106266034730521,
    -0.27330046675043937205975812924953, -0.31180625324666780814299756569942,
    -0.33832950029358816956728612239141, -0.35185093438159561475651955574599
  },{
     0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
     0.19642373959677554531947434191430,  0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
     0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
     0.29396890060483967924361677615282,  0.34675996133053686545540479789161
  },{
     0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
     0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
    -0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
    -0.27330046675043937205975812924953, -0.10263113188058934528909230943878,
     0.10263113188058934528909230943878,  0.27330046675043937205975812924953,
     0.35185093438159561475651955574599,  0.31180625324666780814299756569942,
     0.16666391461943662432490137779072, -0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521, -0.33832950029358816956728612239141
  },{
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680
  },{
     0.31180625324666780814299756569942,  0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953, -0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878,  0.22429189658565907106266034730521,
     0.35185093438159561475651955574599,  0.16666391461943662432490137779072,
    -0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
    -0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
     0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
    -0.034654292299772865648933749133244,-0.31180625324666780814299756569942
  },{
     0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
     0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
     0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
     0.34675996133053686545540479789161,  0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343, 0.29396890060483967924361677615282
  },{
     0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
    -0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
     0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
    -0.31180625324666780814299756569942, -0.22429189658565907106266034730521,
     0.22429189658565907106266034730521,  0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878, -0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244, 0.33832950029358816956728612239141,
     0.16666391461943662432490137779072, -0.27330046675043937205975812924953
  },{
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25
  },{
     0.22429189658565907106266034730521, -0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878,  0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244,-0.33832950029358816956728612239141,
     0.16666391461943662432490137779072,  0.27330046675043937205975812924953,
    -0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
     0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
    -0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
     0.31180625324666780814299756569942, -0.22429189658565907106266034730521
  },{
     0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
     0.068974844820735753083989390917343, 0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
     0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
     0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161,  0.19642373959677554531947434191430
  },{
     0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
     0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
    -0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
     0.034654292299772865648933749133244,-0.31180625324666780814299756569942,
     0.31180625324666780814299756569942, -0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953,  0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878, -0.22429189658565907106266034730521,
     0.35185093438159561475651955574599, -0.16666391461943662432490137779072
  },{
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160
  },{
     0.10263113188058934528909230943878, -0.27330046675043937205975812924953,
     0.35185093438159561475651955574599, -0.31180625324666780814299756569942,
     0.16666391461943662432490137779072,  0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521,  0.33832950029358816956728612239141,
    -0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
    -0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
     0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
     0.27330046675043937205975812924953, -0.10263113188058934528909230943878
  },{
     0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
     0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
     0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
     0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282,  0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430,  0.068974844820735753083989390917343
  },{
     0.034654292299772865648933749133244,-0.10263113188058934528909230943878,
     0.16666391461943662432490137779072, -0.22429189658565907106266034730521,
     0.27330046675043937205975812924953, -0.31180625324666780814299756569942,
     0.33832950029358816956728612239141, -0.35185093438159561475651955574599,
     0.35185093438159561475651955574599, -0.33832950029358816956728612239141,
     0.31180625324666780814299756569942, -0.27330046675043937205975812924953,
     0.22429189658565907106266034730521, -0.16666391461943662432490137779072,
     0.10263113188058934528909230943878, -0.034654292299772865648933749133244
  }
};

void ne_fdct4_double(double _x[],const double _y[],int _in_stride){
  double t[4];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++){
      t[j]+=DCT4_BASIS[j][i]*_y[i*_in_stride];
    }
  }
  for(j=0;j<4;j++){
    _x[j]=t[j];
  }
}

void ne_idct4_double(double _x[],int _out_stride,const double _y[]){
  double t[4];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++){
      t[j]+=DCT4_BASIS[i][j]*_y[i];
    }
  }
  for(j=0;j<4;j++){
    _x[j*_out_stride]=t[j];
  }
}

void ne_fdct8_double(double _x[],const double _y[],int _in_stride){
  double t[8];
  int    i;
  int    j;
  for(j=0;j<8;j++){
    t[j]=0;
    for(i=0;i<8;i++){
      t[j]+=DCT8_BASIS[j][i]*_y[i*_in_stride];
    }
  }
  for(j=0;j<8;j++){
    _x[j]=t[j];
  }
}

void ne_idct8_double(double _x[],int _out_stride,const double _y[]){
  double t[8];
  int    i;
  int    j;
  for(j=0;j<8;j++){
    t[j]=0;
    for(i=0;i<8;i++){
      t[j]+=DCT8_BASIS[i][j]*_y[i];
    }
  }
  for(j=0;j<8;j++){
    _x[j*_out_stride]=t[j];
  }
}

void ne_fdct16_double(double _x[],const double _y[],int _in_stride){
  double t[16];
  int    i;
  int    j;
  for(j=0;j<16;j++){
    t[j]=0;
    for(i=0;i<16;i++){
      t[j]+=DCT16_BASIS[j][i]*_y[i*_in_stride];
    }
  }
  for(j=0;j<16;j++){
    _x[j]=t[j];
  }
}

void ne_idct16_double(double _x[],int _out_stride,const double _y[]){
  double t[16];
  int    i;
  int    j;
  for(j=0;j<16;j++){
    t[j]=0;
    for(i=0;i<16;i++){
      t[j]+=DCT16_BASIS[i][j]*_y[i];
    }
  }
  for(j=0;j<16;j++){
    _x[j*_out_stride]=t[j];
  }
}

typedef void (*ne_fdct_func_1d)(double *_out,const double *_in,int _in_stride);
typedef void (*ne_idct_func_1d)(double *_out,int _out_stride,const double *_in);

const ne_fdct_func_1d OD_FDCT_1D_DOUBLE[OD_NBSIZES]={
  ne_fdct4_double,
  ne_fdct8_double,
  ne_fdct16_double
};

const ne_idct_func_1d OD_IDCT_1D_DOUBLE[OD_NBSIZES]={
  ne_idct4_double,
  ne_idct8_double,
  ne_idct16_double
};

/* 2-dimensional AR(1) model using parameter _r */
void auto_regressive_collapsed(double *_out,int _sz,int _n,double _r){
  int m;
  int j;
  int i;
  m=_sz/_n;
  for(j=0;j<_n;j++){
    _out[j*m]=j==0?1:_out[(j-1)*m]*_r;
    for(i=1;i<m;i++){
      _out[j*m+i]=_out[j*m+i-1]*_r;
    }
  }
}

void analysis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n,const int *_f){
  int    i;
  int    j;
  double t[2*B_SZ];
  for(i=0;i<_n;i++){
    for(j=0;j<2*B_SZ;j++){
      t[j]=_in[j*_in_stride+i];
    }
#if !NE_DISABLE_FILTER
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[0],&t[0],_f);
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ],&t[B_SZ],_f);
#else
# error "Need a prefilter implementation for this block size."
#endif
#endif
#if !NE_DISABLE_TRANSFORM
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_FDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ/2],&t[B_SZ/2],1);
#else
# error "Need an fDCT implementation for this block size."
#endif
#endif
    for(j=0;j<B_SZ;j++){
      _out[j*_out_stride+i]=t[B_SZ/2+j];
    }
  }
}

void synthesis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n,const int *_f){
  int    i;
  int    j;
  double t[2*B_SZ];
  for(i=0;i<_n;i++){
    for(j=0;j<2*B_SZ;j++){
      t[j]=0;
    }
    for(j=0;j<B_SZ;j++){
      t[B_SZ/2+j]=_in[j*_in_stride+i];
    }
#if !NE_DISABLE_TRANSFORM
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ/2],1,&t[B_SZ/2]);
#else
# error "Need an iDCT implementation for this block size."
#endif
#endif
#if !NE_DISABLE_FILTER
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[0],&t[0],_f);
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ],&t[B_SZ],_f);
#else
# error "Need a postfilter implementation for this block size."
#endif
#endif
    for(j=0;j<2*B_SZ;j++){
      _out[j*_out_stride+i]=t[j];
    }
  }
}

#if PRINT_CG_MATH
static void print_matrix(FILE *_fp,const char *_label,const double *_mat,
 int _stride,int _n,int _m){
  int i;
  int j;
  fprintf(_fp,"%s=\n",_label);
  for(j=0;j<_n;j++){
    for(i=0;i<_m;i++){
      fprintf(_fp,"%s  %- 12.6G",i>0?",":"",_mat[j*_m+i]);
    }
    fprintf(_fp,"\n");
  }
}
#endif

double coding_gain_1d(const double _r[2*B_SZ*2*B_SZ],const int *_f){
  int    j;
  int    i;
  double r[B_SZ];
  double rgt[2*B_SZ*B_SZ];
  double grgt[B_SZ*B_SZ];
  double cg;
  /* R*G^T */ 
  for(j=0;j<2*B_SZ;j++){
    analysis(&rgt[j*B_SZ],1,&_r[2*B_SZ*j],1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"rgt",rgt,B_SZ,2*B_SZ,B_SZ);
#endif
  /* G*R*G^T */
  analysis(grgt,B_SZ,rgt,B_SZ,B_SZ,_f);
#if PRINT_CG_MATH
  print_matrix(stdout,"grgt",grgt,B_SZ,B_SZ,B_SZ);
#endif
  /* H */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      r[i]=i==j?1:0;
    }
    synthesis(&rgt[j],B_SZ,r,1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hi",rgt,B_SZ,2*B_SZ,B_SZ);
  fprintf(stdout,"G,H=\n");
#endif
  /* (G*R*G^T)_ii * (H^T*H)_ii */
  cg=0;
  for(j=0;j<B_SZ;j++){
    double h;
    h=0;
    for(i=0;i<2*B_SZ;i++){
      h+=rgt[i*B_SZ+j]*rgt[i*B_SZ+j];
    }
#if PRINT_CG_MATH
    fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",grgt[j*B_SZ+j],h);
#endif
    cg-=10*log10(grgt[j*B_SZ+j]*h);
  }
  return cg/B_SZ;
}

double coding_gain_1d_collapsed(const double _r[2*B_SZ],const int *_f){
  int    j;
  int    i;
  double r[2*B_SZ];
  double rgt[2*B_SZ*B_SZ];
  double grgt[B_SZ*B_SZ];
  double cg;
  /* R*G^T */ 
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      r[i]=_r[abs(i-j)];
    }
    analysis(&rgt[j*B_SZ],1,r,1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"rgt",rgt,B_SZ,2*B_SZ,B_SZ);
#endif
  /* G*R*G^T */
  analysis(grgt,B_SZ,rgt,B_SZ,B_SZ,_f);
#if PRINT_CG_MATH
  print_matrix(stdout,"grgt",grgt,B_SZ,B_SZ,B_SZ);
#endif
  /* H */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      r[i]=i==j?1:0;
    }
    synthesis(&rgt[j],B_SZ,r,1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hi",rgt,B_SZ,2*B_SZ,B_SZ);
  fprintf(stdout,"G,H=\n");
#endif
  /* (G*R*G^T)_ii * (H^T*H)_ii */
  cg=0;
  for(j=0;j<B_SZ;j++){
    double h;
    h=0;
    for(i=0;i<2*B_SZ;i++){
      h+=rgt[i*B_SZ+j]*rgt[i*B_SZ+j];
    }
#if PRINT_CG_MATH
    fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",grgt[j*B_SZ+j],h);
#endif
    cg-=10*log10(grgt[j*B_SZ+j]*h);
  }
  return cg/B_SZ;
}

double coding_gain_2d(const double _r[2*B_SZ*2*B_SZ*2*B_SZ*2*B_SZ],const int *_f){
  int    v;
  int    u;
  int    j;
  int    i;
  double r[2*B_SZ];
  double s[2*B_SZ*B_SZ];
  double rggt[2*B_SZ*2*B_SZ*B_SZ*B_SZ];
  double ggrggt[B_SZ*B_SZ*B_SZ*B_SZ];
  double cg;
  /* R*(G2*P*G1)^T */
  for(v=0;v<2*B_SZ;v++){
    for(j=0;j<2*B_SZ;j++){
      for(u=0;u<2*B_SZ;u++){
        analysis(&s[u*B_SZ],1,&_r[2*B_SZ*2*B_SZ*2*B_SZ*u+2*B_SZ*v+j],2*B_SZ*2*B_SZ,1,_f);
      }
      analysis(&rggt[(v*2*B_SZ+j)*B_SZ*B_SZ],B_SZ,s,B_SZ,B_SZ,_f);
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"rggt",rggt,B_SZ*B_SZ,4*B_SZ*B_SZ,B_SZ*B_SZ);
#endif
  /* G1*P*G2*R*(G2*P*G1)^T */
  for(v=0;v<B_SZ;v++){
    for(j=0;j<B_SZ;j++){
      for(u=0;u<2*B_SZ;u++){
        analysis(&s[u*B_SZ],1,&rggt[2*B_SZ*B_SZ*B_SZ*u+v*B_SZ+j],B_SZ*B_SZ,1,_f);
      }
      analysis(&ggrggt[(v*B_SZ+j)*B_SZ*B_SZ],B_SZ,s,B_SZ,B_SZ,_f);
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"ggrggt",ggrggt,B_SZ*B_SZ,B_SZ*B_SZ,B_SZ*B_SZ);
#endif
  /* H */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      r[i]=i==j?1:0;
    }
    synthesis(&s[j],B_SZ,r,1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hi",s,B_SZ,2*B_SZ,B_SZ);
#endif
  /* H1*P*H2 */
  for(i=0;i<2*B_SZ;i++){
    r[i]=0;
  }
  for(u=0;u<2*B_SZ;u++){
    for(j=0;j<B_SZ;j++){
      r[B_SZ]=s[B_SZ*u+j];
      for(i=0;i<B_SZ;i++){
        synthesis(&rggt[B_SZ*B_SZ*u*2*B_SZ+j+i*B_SZ],B_SZ*B_SZ,&r[B_SZ-i],1,1,_f);
      }
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hhi",rggt,B_SZ*B_SZ,4*B_SZ*B_SZ,B_SZ*B_SZ);
  fprintf(stdout,"G,H=\n");
#endif
  /* ((H1*P*H2)^T*H1*P*H2)_ii */
  for(j=0;j<B_SZ*B_SZ;j++){
    s[j]=0;
    for(i=0;i<4*B_SZ*B_SZ;i++){
      s[j]+=rggt[i*B_SZ*B_SZ+j]*rggt[i*B_SZ*B_SZ+j];
    }
  }
  /* (G1*P*G2*R*(G1*P*G2)^T)_ii * ((H1*P*H2)^T*H1*P*H2)_ii */
  cg=0;
  for(j=0;j<B_SZ*B_SZ;j++){
    fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",ggrggt[B_SZ*B_SZ*j+j],s[j]);
    cg-=10*log10(ggrggt[B_SZ*B_SZ*j+j]*s[j]);
  }
  return cg/(B_SZ*B_SZ);
}

double coding_gain_2d_collapsed(const double _r[2*B_SZ*2*B_SZ],const int *_f){
  int    v;
  int    u;
  int    j;
  int    i;
  double r[2*B_SZ];
  double s[2*B_SZ*B_SZ];
  double rggt[2*B_SZ*2*B_SZ*B_SZ*B_SZ];
  double ggrggt[B_SZ*B_SZ*B_SZ*B_SZ];
  double cg;
#if NN_SEARCH
  double basis[B_SZ];
  for(v=0;v<B_SZ;v++)basis[v]=v<(B_SZ>>1)?0:255;
  (NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(basis,basis,_f);
  for(v=0;v<B_SZ;v++)if(basis[v]<0)return -100;
#endif
  /* R*(G2*P*G1)^T */
  for(v=0;v<2*B_SZ;v++){
    for(j=0;j<2*B_SZ;j++){
      for(u=0;u<2*B_SZ;u++){
        for(i=0;i<2*B_SZ;i++){
          r[i]=_r[abs(u-v)*2*B_SZ+abs(i-j)];
        }
        analysis(&s[u*B_SZ],1,r,1,1,_f);
      }
      analysis(&rggt[(v*2*B_SZ+j)*B_SZ*B_SZ],B_SZ,s,B_SZ,B_SZ,_f);
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"rggt",rggt,B_SZ*B_SZ,4*B_SZ*B_SZ,B_SZ*B_SZ);
#endif
  /* G1*P*G2*R*(G2*P*G1)^T */
  for(v=0;v<B_SZ;v++){
    for(j=0;j<B_SZ;j++){
      for(u=0;u<2*B_SZ;u++){
        analysis(&s[u*B_SZ],1,&rggt[2*B_SZ*B_SZ*B_SZ*u+v*B_SZ+j],B_SZ*B_SZ,1,_f);
      }
      analysis(&ggrggt[(v*B_SZ+j)*B_SZ*B_SZ],B_SZ,s,B_SZ,B_SZ,_f);
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"ggrggt",ggrggt,B_SZ*B_SZ,B_SZ*B_SZ,B_SZ*B_SZ);
#endif
  /* H */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      r[i]=i==j?1:0;
    }
    synthesis(&s[j],B_SZ,r,1,1,_f);
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hi",s,B_SZ,2*B_SZ,B_SZ);
#endif
  /* H1*P*H2 */
  for(i=0;i<2*B_SZ;i++){
    r[i]=0;
  }
  for(u=0;u<2*B_SZ;u++){
    for(j=0;j<B_SZ;j++){
      r[B_SZ]=s[B_SZ*u+j];
      for(i=0;i<B_SZ;i++){
        synthesis(&rggt[B_SZ*B_SZ*u*2*B_SZ+j+i*B_SZ],B_SZ*B_SZ,&r[B_SZ-i],1,1,_f);
      }
    }
  }
#if PRINT_CG_MATH
  print_matrix(stdout,"hhi",rggt,B_SZ*B_SZ,4*B_SZ*B_SZ,B_SZ*B_SZ);
  fprintf(stdout,"G,H=\n");
#endif
  /* ((H1*P*H2)^T*H1*P*H2)_ii */
  for(j=0;j<B_SZ*B_SZ;j++){
    s[j]=0;
    for(i=0;i<4*B_SZ*B_SZ;i++){
      s[j]+=rggt[i*B_SZ*B_SZ+j]*rggt[i*B_SZ*B_SZ+j];
    }
  }
  /* (G1*P*G2*R*(G1*P*G2)^T)_ii * ((H1*P*H2)^T*H1*P*H2)_ii */
  cg=0;
  for(j=0;j<B_SZ*B_SZ;j++){
#if PRINT_CG_MATH
    fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",ggrggt[B_SZ*B_SZ*j+j],s[j]);
#endif
    cg-=10*log10(ggrggt[B_SZ*B_SZ*j+j]*s[j]);
  }
  return cg/(B_SZ*B_SZ);
}
