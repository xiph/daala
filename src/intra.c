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

#include <stdlib.h>
#include <limits.h>
#include "filter.h"
#include "intra.h"
#include "tf.h"

const od_intra_mult_func OD_INTRA_MULT[OD_NBSIZES+1]={
  od_intra_pred4x4_mult,
  od_intra_pred8x8_mult,
  od_intra_pred16x16_mult
};

void od_intra_pred4x4_mult(double *_pred,int _pred_stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  const int *paramx;
  const int *paramy;
  const double *weights;
  paramx=OD_PRED_PARAMX_4x4_DATA+OD_PRED_OFFSETS_4x4[_mode];
  paramy=OD_PRED_PARAMY_4x4_DATA+OD_PRED_OFFSETS_4x4[_mode];
  weights=OD_PRED_WEIGHTS_4x4_DATA+OD_PRED_OFFSETS_4x4[_mode];
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      double sum=0;
      for(k=0;k<OD_PRED_MULTS_4x4[_mode][j][i];k++){
        od_coeff *neighbor;
        int neighbor_stride;
        int neighbori;
        /*The values in the arrays are relative to the x & y of the upper-left
           block.*/
        x=*paramx++;
        y=*paramy++;
        OD_ASSERT(0<=x&&x<12&&0<=y&&y<8);
        neighbori = (x>>2)+3*(y>>2);
        x&=3;
        y&=3;
        neighbor=_neighbors[neighbori];
        neighbor_stride=_neighbor_strides[neighbori];
        sum+=
         neighbor[neighbor_stride*y+x]**weights++;
      }
      _pred[_pred_stride*j+i]=sum;
    }
  }
}

void od_intra_pred8x8_mult(double *_pred,int _pred_stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
      _pred[_pred_stride*j+i]=0;
      for(k=0;k<OD_PRED_MULTS_8x8[_mode][j][i];k++){
        od_coeff *neighbor;
        int neighbor_stride;
        int neighbori;
        x=OD_PRED_PARAMX_8x8[_mode][j][i][k];
        y=OD_PRED_PARAMY_8x8[_mode][j][i][k];
        OD_ASSERT(0<=x&&x<24&&0<=y&&y<16);
        if(y<8){
          if(x<8)neighbori=0;
          else if(x<16){
            neighbori=1;
            x-=8;
          }
          else{
            neighbori=2;
            x-=16;
          }
        }
        else{
          neighbori=3;
          y-=8;
        }
        neighbor=_neighbors[neighbori];
        neighbor_stride=_neighbor_strides[neighbori];
        _pred[_pred_stride*j+i]+=
         neighbor[neighbor_stride*y+x]*OD_PRED_WEIGHTS_8x8[_mode][j][i][k];
      }
    }
  }
}

void od_intra_pred16x16_mult(double *_pred,int _pred_stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
      _pred[_pred_stride*j+i]=0;
      for(k=0;k<OD_PRED_MULTS_16x16[_mode][j][i];k++){
        od_coeff *neighbor;
        int neighbor_stride;
        int neighbori;
        x=OD_PRED_PARAMX_16x16[_mode][j][i][k];
        y=OD_PRED_PARAMY_16x16[_mode][j][i][k];
        OD_ASSERT(0<=x&&x<48&&0<=y&&y<32);
        if(y<16){
          if(x<16)neighbori=0;
          else if(x<32){
            neighbori=1;
            x-=16;
          }
          else{
            neighbori=2;
            x-=32;
          }
        }
        else{
          neighbori=3;
          y-=16;
        }
        neighbor=_neighbors[neighbori];
        neighbor_stride=_neighbor_strides[neighbori];
        _pred[_pred_stride*j+i]+=
         neighbor[neighbor_stride*y+x]*OD_PRED_WEIGHTS_16x16[_mode][j][i][k];
      }
    }
  }
}

static const float OD_SATD_WEIGHTS_4x4[3][4*4]={
  {
    0.230317,0.329547,0.457088,0.517193,0.315611,0.389662,0.525445,0.577099,
    0.424841,0.506366,0.648566,0.696878,0.473835,0.543352,0.681457,0.720876
  },
  {
    0.458235,0.617704,0.811118,0.890949,0.601811,0.711183,0.906126,0.962359,
    0.776575,0.882826,1.038644,1.073067,0.851901,0.929541,1.065617,1.084039
  },
  {
    0.481204,0.631404,0.829253,0.905577,0.615889,0.736232,0.930215,0.983703,
    0.802367,0.908049,1.060383,1.091684,0.872241,0.953467,1.087123,1.106003
  }
};

static const float OD_SATD_WEIGHTS_8x8[3][8*8]={
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  },
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  },
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  }
};

static const float OD_SATD_WEIGHTS_16x16[3][16*16]={
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  },
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  },
  {
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000
  }
};

const od_intra_dist_func OD_INTRA_DIST[OD_NBSIZES+1]={
  od_intra_pred4x4_dist,
  od_intra_pred8x8_dist,
  od_intra_pred16x16_dist
};

void od_intra_pred4x4_dist(ogg_uint32_t *_dist,const od_coeff *_c,int _stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4], int _pli){
  double p[4*4];
  float  satd;
  int    mode;
  int    i;
  int    j;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    od_intra_pred4x4_mult(p,4,_neighbors,_neighbor_strides,mode);
    satd=0;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        satd+=
         fabs((_c[_stride*i+j]-p[i*4+j])*OD_SATD_WEIGHTS_4x4[_pli][i*4+j]);
      }
    }
    _dist[mode]=satd;
  }
}

void od_intra_pred8x8_dist(ogg_uint32_t *_dist,const od_coeff *_c,int _stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _pli){
  double p[8*8];
  float  satd;
  int    mode;
  int    i;
  int    j;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    od_intra_pred8x8_mult(p,8,_neighbors,_neighbor_strides,mode);
    satd=0;
    for(i=0;i<8;i++){
      for(j=0;j<8;j++){
        satd+=
         fabs((_c[_stride*i+j]-p[i*8+j])*OD_SATD_WEIGHTS_8x8[_pli][i*8+j]);
      }
    }
    _dist[mode]=satd;
  }
}

void od_intra_pred16x16_dist(ogg_uint32_t *_dist,
 const od_coeff *_c,int _stride,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _pli){
  double p[16*16];
  float  satd;
  int    mode;
  int    i;
  int    j;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    od_intra_pred16x16_mult(p,16,_neighbors,_neighbor_strides,mode);
    satd=0;
    for(i=0;i<16;i++){
      for(j=0;j<16;j++){
        satd+=
         fabs((_c[_stride*i+j]-p[i*16+j])*OD_SATD_WEIGHTS_16x16[_pli][i*16+j]);
      }
    }
    _dist[mode]=satd;
  }
}

/*These are used to weight the samples we use to build the chroma-from-luma
   model.
  They're constructed from the DC weights used to predict the DC coefficient in
   each INTRA mode.
  We look up the weights using the INTRA modes of the corresponding luma
   blocks.*/
const signed char OD_INTRA_CHROMA_WEIGHTS_Q6[OD_INTRA_NMODES][3]={
  { -3, 35, 32},
  {-23, 41, 46},
  {  5,-29, 88},
  { 14,-13, 63},
  { 34,-10, 40},
  { 43,  9, 12},
  { 32, 39, -8},
  { 11, 65,-12},
  {  5, 80,-21},
  { 64,  0, -1}
};

void od_chroma_pred4x4(od_coeff *_p,const od_coeff *_c,
 const od_coeff *_l,int _stride,const int _weights_q8[3]){
  static const int BLOCK_DX[3]={-4,0,-4};
  static const int BLOCK_DY[3]={-4,-4,0};
  static const int AC_DX[3]={1,0,1};
  static const int AC_DY[3]={0,1,1};
  ogg_int64_t xx;
  ogg_int64_t xy;
  ogg_int32_t alpha_q8;
  ogg_int32_t beta_q8;
  ogg_int32_t lc_sum_q8;
  ogg_int32_t cc_sum_q8;
  int         bi;
  int         i;
  int         j;
  /*Solve for a simple predictive model using the UL, U, L neighbors:
    chroma DC = luma DC * alpha + beta
    chroma AC = luma AC * alpha*/
  lc_sum_q8=cc_sum_q8=0;
  xx=xy=0;
  for(bi=0;bi<3;bi++){
    od_coeff lc;
    od_coeff cc;
    int      boffs;
    int      ci;
    int      w_q8;
    boffs=BLOCK_DY[bi]*_stride+BLOCK_DX[bi];
    lc=*(_l+boffs);
    cc=*(_c+boffs);
    w_q8=_weights_q8[bi];
    xx+=lc*(ogg_int64_t)lc*w_q8;
    xy+=lc*(ogg_int64_t)cc*w_q8;
    lc_sum_q8+=lc*w_q8;
    cc_sum_q8+=cc*w_q8;
    for(ci=0;ci<3;ci++){
      int coffs;
      coffs=boffs+AC_DY[ci]*_stride+AC_DX[ci];
      lc=*(_l+coffs);
      cc=*(_c+coffs);
      xx+=lc*(ogg_int64_t)lc*w_q8;
      xy+=cc*(ogg_int64_t)cc*w_q8;
    }
  }
  xx-=lc_sum_q8*(ogg_int64_t)lc_sum_q8+128>>8;
  xy-=cc_sum_q8*(ogg_int64_t)lc_sum_q8+128>>8;
  if(abs(xx)>abs(xy)>>1)alpha_q8=(xy<<8)/xx;
  else alpha_q8=0;
  beta_q8=cc_sum_q8-(alpha_q8*lc_sum_q8+128>>8);
  _p[0]=_l[0]*alpha_q8+beta_q8+128>>8;
  for(i=0;i<4;i++){
    for(j=i==0;j<4;j++){
      _p[i*4+j]=_l[i*_stride+j]*alpha_q8+128>>8;
    }
  }
}

ogg_uint32_t od_chroma_pred4x4_dist(const od_coeff *_c,
 const od_coeff *_l,int _stride,const int _weights_q8[3],int _pli){
  od_coeff p[4*4];
  float    satd;
  int      i;
  int      j;
  od_chroma_pred4x4(p,_c,_l,_stride,_weights_q8);
  satd=0;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      satd+=abs(_c[_stride*i+j]-p[i*4+j])*OD_SATD_WEIGHTS_4x4[_pli][i*4+j];
    }
  }
  return (ogg_uint32_t)satd;
}

const od_intra_get_func OD_INTRA_GET[OD_NBSIZES+1]={
  od_intra_pred4x4_get,
  od_intra_pred8x8_get,
  od_intra_pred16x16_get
};

void od_intra_pred4x4_get(od_coeff *_out,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  double p[4*4];
  int    i;
  int    j;
  od_intra_pred4x4_mult(p,4,_neighbors,_neighbor_strides,_mode);
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _out[4*i+j]=(od_coeff)floor(p[i*4+j]+0.5);
    }
  }
}

void od_intra_pred8x8_get(od_coeff *_out,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  double p[8*8];
  int    i;
  int    j;
  od_intra_pred8x8_mult(p,8,_neighbors,_neighbor_strides,_mode);
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      _out[8*i+j]=(od_coeff)floor(p[i*8+j]+0.5);
    }
  }
}

void od_intra_pred16x16_get(od_coeff *_out,
 od_coeff *_neighbors[4],int _neighbor_strides[4],int _mode){
  double p[16*16];
  int    i;
  int    j;
  od_intra_pred16x16_mult(p,16,_neighbors,_neighbor_strides,_mode);
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      _out[16*i+j]=(od_coeff)floor(p[i*16+j]+0.5);
    }
  }
}

void od_intra_pred_cdf(ogg_uint16_t _cdf[],
 const unsigned char _probs[][OD_INTRA_NCONTEXTS],const ogg_uint16_t _p0[],
 int _nmodes,int _left,int _upleft,int _up){
  unsigned p[OD_INTRA_NMODES+1];
  int      mi;
  int      sum;
  int      curr_cdf;
  sum=0;
  for(mi=0;mi<_nmodes;mi++){
    p[mi]=_probs[mi][(_left==mi)*4+(_upleft==mi)*2+(_up==mi)];
    p[mi]=OD_MAXI(OD_MAXI(p[mi],_p0[mi]>>8),1);
    sum+=p[mi];
  }
  curr_cdf=0;
  for(mi=0;mi<_nmodes;mi++){
    /*Apply probability combination here: p[mi] *= (sum-p[mi])/(1-p[mi]).*/
    /*FIXME: Make this fixed-point.*/
    p[mi]=p[mi]*(sum-p[mi])/(float)(256-p[mi]);
    curr_cdf+=p[mi];
    _cdf[mi]=curr_cdf;
  }
}

int od_intra_pred_search(const ogg_uint16_t _cdf[],
 const ogg_uint32_t _dist[],int _nmodes,ogg_uint16_t _lambda){
  int best_score;
  int best_mode;
  int mi;
  /*FIXME: Compute the log2() in fixed-point.*/
  best_score=_dist[0]-_lambda*(M_LOG2E/128)*log(_cdf[0]);
  best_mode=0;
  for(mi=1;mi<_nmodes;mi++){
    int score;
    /*FIXME: Compute the log2() in fixed-point.*/
    score=_dist[mi]-_lambda*(M_LOG2E/128)*log(_cdf[mi]-_cdf[mi-1]);
    if(score<best_score){
      best_score=score;
      best_mode=mi;
    }
  }
  return best_mode;
}

void od_intra_pred_update(ogg_uint16_t _p0[],int _nmodes,int _mode,
 int _left,int _upleft,int _up){
  int mi;
  for(mi=1;mi<_nmodes;mi++){
    /*p0 is the probability of choosing mode mi when none of its neighbors
       are in mode mi.
      So decay p0 if this block satisfies that condition.*/
    if(_left!=mi&&_up!=mi&&_upleft!=mi)_p0[mi]-=_p0[mi]+256>>9;
  }
  /*And bump up the probability in the mode we actually chose.*/
  if(_left!=_mode&&_up!=_mode&&_upleft!=_mode){
    /*Arbitrary ceiling at 0.75 to prevent insane p0 values.*/
    if(_p0[_mode]<24576)_p0[_mode]+=64;
  }
}

void od_resample_luma_coeffs(od_coeff *l, int lstride,
 const od_coeff *c, int cstride, int xdec, int ydec, int n) {
  if (xdec) {
    if (ydec) od_tf_up_hv_lp(l, lstride, c, cstride, n, n, n);
    else od_tf_up_h_lp(l, lstride, c, cstride, n, n);
  }
  else{
    OD_ASSERT(ydec);
    if (ydec) od_tf_up_v_lp(l, lstride, c, cstride, n, n);
  }
}
