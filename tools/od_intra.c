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
#include "od_intra.h"
#include "../src/intra.h"

static void ne_intra_pred4x4_mult(double *_pred,int _stride,
 const od_coeff *_coeff,int _mode){
  double *weight;
  int    *index;
  int    j;
  int    i;
  int    k;
  weight=NE_PRED_WEIGHTS_4x4[_mode];
  index=NE_PRED_INDEX_4x4[_mode];
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      double sum;
#if ZERO_MEAN
      sum=0;
#else
      sum=NE_PRED_OFFSETS_4x4[_mode][j][i];
#endif
      for(k=NE_PRED_MULTS_4x4[_mode][j][i];k-->0;){
        sum+=_coeff[*index]*(*weight);
        index++;
        weight++;
      }
      _pred[_stride*j+i]=sum;
    }
  }
}

static void ne_intra_pred8x8_mult(double *_pred,int _stride,
 const od_coeff *_coeff,int _mode){
  double *weight;
  int    *index;
  int     j;
  int     i;
  int     k;
  weight=NE_PRED_WEIGHTS_8x8[_mode];
  index=NE_PRED_INDEX_8x8[_mode];
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
      double  sum;
#if ZERO_MEAN
      sum=0;
#else
      sum=NE_PRED_OFFSETS_8x8[_mode][j][i];
#endif
      for(k=NE_PRED_MULTS_8x8[_mode][j][i];k-->0;){
        sum+=_coeff[*index]*(*weight);
        index++;
        weight++;
      }
      _pred[_stride*j+i]=sum;
    }
  }
}

static void ne_intra_pred16x16_mult(double *_pred,int _stride,
 const od_coeff *_coeff,int _mode){
  double *weight;
  int    *index;
  int     j;
  int     i;
  int     k;
  weight=NE_PRED_WEIGHTS_16x16[_mode];
  index=NE_PRED_INDEX_16x16[_mode];
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
      double  sum;
#if ZERO_MEAN
      sum=0;
#else
      sum=NE_PRED_OFFSETS_16x16[_mode][j][i];
#endif
      for(k=NE_PRED_MULTS_16x16[_mode][j][i];k-->0;){
        sum+=_coeff[*index]*(*weight);
        index++;
        weight++;
      }
      _pred[_stride*j+i]=sum;
    }
  }
}

const ne_intra_mult_func NE_INTRA_MULT[OD_NBSIZES]={
  ne_intra_pred4x4_mult,
  ne_intra_pred8x8_mult,
  ne_intra_pred16x16_mult
};

void od_intra_init(){
  int m;
  for(m=0;m<OD_INTRA_NMODES;m++){
    const double *weights;
    const ogg_uint16_t *index;
    int w;
    int b_sz;
    int j;
    int i;
    int k;
    w=0;
    b_sz=4;
    NE_PRED_WEIGHTS_4x4[m]=(double *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(double));
    NE_PRED_INDEX_4x4[m]=(int *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(int));
    index=OD_PRED_INDEX_4x4+OD_PRED_OFFSETS_4x4[m];
    weights=OD_PRED_WEIGHTS_4x4+OD_PRED_OFFSETS_4x4[m];
    for(j=0;j<b_sz;j++){
      for(i=0;i<b_sz;i++){
        NE_PRED_MULTS_4x4[m][j][i]=OD_PRED_MULTS_4x4[m][j][i];
        for(k=OD_PRED_MULTS_4x4[m][j][i];k-->0;w++){
          NE_PRED_WEIGHTS_4x4[m][w]=*weights;
          NE_PRED_INDEX_4x4[m][w]=*index;
          weights++;
          index++;
        }
      }
    }
    w=0;
    b_sz=8;
    NE_PRED_WEIGHTS_8x8[m]=(double *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(double));
    NE_PRED_INDEX_8x8[m]=(int *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(int));
    index=OD_PRED_INDEX_8x8+OD_PRED_OFFSETS_8x8[m];
    weights=OD_PRED_WEIGHTS_8x8+OD_PRED_OFFSETS_8x8[m];
    for(j=0;j<b_sz;j++){
      for(i=0;i<b_sz;i++){
        NE_PRED_MULTS_8x8[m][j][i]=OD_PRED_MULTS_8x8[m][j][i];
        for(k=OD_PRED_MULTS_8x8[m][j][i];k-->0;w++){
          NE_PRED_WEIGHTS_8x8[m][w]=*weights;
          NE_PRED_INDEX_8x8[m][w]=*index;
          weights++;
          index++;
        }
      }
    }
    w=0;
    b_sz=16;
    NE_PRED_WEIGHTS_16x16[m]=(double *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(double));
    NE_PRED_INDEX_16x16[m]=(int *)malloc(b_sz*b_sz*4*b_sz*b_sz*sizeof(int));
    index=OD_PRED_INDEX_16x16+OD_PRED_OFFSETS_16x16[m];
    weights=OD_PRED_WEIGHTS_16x16+OD_PRED_OFFSETS_16x16[m];
    for(j=0;j<b_sz;j++){
      for(i=0;i<b_sz;i++){
        NE_PRED_MULTS_16x16[m][j][i]=OD_PRED_MULTS_16x16[m][j][i];
        for(k=OD_PRED_MULTS_16x16[m][j][i];k-->0;w++){
          NE_PRED_WEIGHTS_16x16[m][w]=*weights;
          NE_PRED_INDEX_16x16[m][w]=*index;
          weights++;
          index++;
        }
      }
    }
  }
}

void od_intra_clear(){
  int m;
  for(m=0;m<OD_INTRA_NMODES;m++){
    free(NE_PRED_WEIGHTS_4x4[m]);
    free(NE_PRED_INDEX_4x4[m]);
    free(NE_PRED_WEIGHTS_8x8[m]);
    free(NE_PRED_INDEX_8x8[m]);
    free(NE_PRED_WEIGHTS_16x16[m]);
    free(NE_PRED_INDEX_16x16[m]);
  }
}

void update_predictors(int _mode,double *_beta_0,double *_beta_1,int *_mask){
  double *weights;
  int    *index;
  int     i;
  int     j;
#if B_SZ==4
  weights=NE_PRED_WEIGHTS_4x4[_mode];
  index=NE_PRED_INDEX_4x4[_mode];
#elif B_SZ==8
  weights=NE_PRED_WEIGHTS_8x8[_mode];
  index=NE_PRED_INDEX_8x8[_mode];
#elif B_SZ==16
  weights=NE_PRED_WEIGHTS_16x16[_mode];
  index=NE_PRED_INDEX_16x16[_mode];
#else
# error "Need predictors for this block size."
#endif
  for(i=0;i<B_SZ*B_SZ;i++){
    int y;
    int x;
    int by;
    int bx;
    int n;
    y=i/B_SZ;
    x=i%B_SZ;
#if B_SZ==4
    NE_PRED_OFFSETS_4x4[_mode][y][x]=_beta_0[i];
#elif B_SZ==8
    NE_PRED_OFFSETS_8x8[_mode][y][x]=_beta_0[i];
#elif B_SZ==16
    NE_PRED_OFFSETS_16x16[_mode][y][x]=_beta_0[i];
#else
# error "Need predictors for this block size."
#endif
    n=0;
    for(by=0;by<=1;by++){
      for(bx=0;bx<=(1-by)<<1;bx++){
        for(j=0;j<B_SZ*B_SZ;j++){
          int v;
          int u;
          int ij;
          v=j/B_SZ;
          u=j%B_SZ;
          ij=4*B_SZ*B_SZ*i+(3*by+bx)*B_SZ*B_SZ+j;
          if(_mask[ij]){
            *weights=_beta_1[ij];
            *index=((((3*by+bx)<<B_SZ_LOG)|v)<<B_SZ_LOG)|u;
            weights++;
            index++;
            n++;
          }
        }
      }
    }
#if B_SZ==4
    NE_PRED_MULTS_4x4[_mode][y][x]=n;
#elif B_SZ==8
    NE_PRED_MULTS_8x8[_mode][y][x]=n;
#elif B_SZ==16
    NE_PRED_MULTS_16x16[_mode][y][x]=n;
#else
# error "Need predictors for this block size."
#endif
  }
}

void print_predictors(FILE *_fp){
  int m;
  int i;
  int j;
  int k;
  int w;
  int offset[OD_INTRA_NMODES];
  fprintf(_fp,"double NE_PRED_OFFSETS_%ix%i[OD_INTRA_NMODES][%i][%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ);
  for (m=0;m<OD_INTRA_NMODES;m++) {
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for (j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for (i=0;i<B_SZ;i++){
#if B_SZ==4
        fprintf(_fp,"%s  %- 24.18G",i>0?",":"",NE_PRED_OFFSETS_4x4[m][j][i]);
#elif B_SZ==8
        fprintf(_fp,"%s  %- 24.18G",i>0?",":"",NE_PRED_OFFSETS_8x8[m][j][i]);
#elif B_SZ==16
        fprintf(_fp,"%s  %- 24.18G",i>0?",":"",NE_PRED_OFFSETS_16x16[m][j][i]);
#else
# error "Unsupported block size."
#endif
      }
      fprintf(_fp,"  }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fprintf(_fp,"const int OD_PRED_MULTS_%ix%i[OD_INTRA_NMODES][%i][%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for(i=0;i<B_SZ;i++){
#if B_SZ==4
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_4x4[m][j][i]);
#elif B_SZ==8
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_8x8[m][j][i]);
#elif B_SZ==16
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_16x16[m][j][i]);
#else
# error "Unsupported block size."
#endif
      }
      fprintf(_fp,"   }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fprintf(_fp,"const double OD_PRED_WEIGHTS_%ix%i[]={",B_SZ,B_SZ);
  offset[0]=0;
  for(m=0;m<OD_INTRA_NMODES;m++){
    w=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
      fprintf(_fp,"\n/* Mode %i (%i %i) */",m,i,j);
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++,w++){
          fprintf(_fp,"%s%- 24.18G,",k&0x7?"":"\n  ",
           NE_PRED_WEIGHTS_4x4[m][w]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++,w++){
          fprintf(_fp,"%s%- 24.18G,",k&0x7?"":"\n  ",
           NE_PRED_WEIGHTS_8x8[m][w]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++,w++){
          fprintf(_fp,"%s%- 24.18G,",k&0x7?"":"\n  ",
           NE_PRED_WEIGHTS_16x16[m][w]);
        }
#else
# error "Unsupported block size."
#endif
      }
    }
    if(m>0){
      offset[m]=offset[m-1]+w;
    }
  }
  fprintf(_fp,"\n};\n");
  fprintf(_fp,"const ogg_uint16_t OD_PRED_INDEX_%ix%i[]={",B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    w=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"\n/* Mode %i (%i %i) */",m,i,j);
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++,w++){
          fprintf(_fp,"%s0x%02x,",k&0xf?"":"\n  ",NE_PRED_INDEX_4x4[m][w]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++,w++){
          fprintf(_fp,"%s0x%02x,",k&0xf?"":"\n  ",NE_PRED_INDEX_8x8[m][w]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++,w++){
          fprintf(_fp,"%s0x%02x,",k&0xf?"":"\n  ",NE_PRED_INDEX_16x16[m][w]);
        }
#else
# error "Unsupported block size."
#endif
      }
    }
  }
  fprintf(_fp,"\n};\n");
  fprintf(_fp,"const int OD_PRED_OFFSETS_%ix%i[OD_INTRA_NMODES] = {\n  ",
   B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"%s%i",m>0?", ":"",offset[m]);
  }
  fprintf(_fp,"\n};\n");
}
