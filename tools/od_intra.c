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

#include "od_intra.h"

static void ne_intra_pred4x4_mult(double *_pred,int _pred_stride,
 const od_coeff *_coeff,int _coeff_stride,int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
#if ZERO_MEAN
      _pred[_pred_stride*j+i]=0;
#else
      _pred[_pred_stride*j+i]=NE_PRED_OFFSETS_4x4[_mode][j][i];
#endif
      for(k=0;k<NE_PRED_MULTS_4x4[_mode][j][i];k++){
        x=NE_PRED_PARAMX_4x4[_mode][j][i][k];
        y=NE_PRED_PARAMY_4x4[_mode][j][i][k];
        _pred[_pred_stride*j+i]+=
         _coeff[_coeff_stride*(y-4)+(x-4)]*NE_PRED_WEIGHTS_4x4[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred8x8_mult(double *_pred,int _pred_stride,
 const od_coeff *_coeff,int _coeff_stride,int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
#if ZERO_MEAN
      _pred[_pred_stride*j+i]=0;
#else
      _pred[_pred_stride*j+i]=NE_PRED_OFFSETS_8x8[_mode][j][i];
#endif
      for(k=0;k<NE_PRED_MULTS_8x8[_mode][j][i];k++){
        x=NE_PRED_PARAMX_8x8[_mode][j][i][k];
        y=NE_PRED_PARAMY_8x8[_mode][j][i][k];
        _pred[_pred_stride*j+i]+=
         _coeff[_coeff_stride*(y-8)+(x-8)]*NE_PRED_WEIGHTS_8x8[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred16x16_mult(double *_pred,int _pred_stride,
 const od_coeff *_coeff,int _coeff_stride,int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
#if ZERO_MEAN
      _pred[_pred_stride*j+i]=0;
#else
      _pred[_pred_stride*j+i]=NE_PRED_OFFSETS_16x16[_mode][j][i];
#endif
      for(k=0;k<NE_PRED_MULTS_16x16[_mode][j][i];k++){
        x=NE_PRED_PARAMX_16x16[_mode][j][i][k];
        y=NE_PRED_PARAMY_16x16[_mode][j][i][k];
        _pred[_pred_stride*j+i]+=_coeff[_coeff_stride*(y-16)+(x-16)]*
         NE_PRED_WEIGHTS_16x16[_mode][j][i][k];
      }
    }
  }
}

const ne_intra_mult_func NE_INTRA_MULT[OD_NBSIZES]={
  ne_intra_pred4x4_mult,
  ne_intra_pred8x8_mult,
  ne_intra_pred16x16_mult
};

void update_predictors(int _mode,double *_beta_0,double *_beta_1,
 int _mask[B_SZ*B_SZ*4*B_SZ*B_SZ]){
  int i;
  int j;
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
#if B_SZ==4
            NE_PRED_WEIGHTS_4x4[_mode][y][x][n]=_beta_1[ij];
            NE_PRED_PARAMX_4x4[_mode][y][x][n]=B_SZ*bx+u;
            NE_PRED_PARAMY_4x4[_mode][y][x][n]=B_SZ*by+v;
#elif B_SZ==8
            NE_PRED_WEIGHTS_8x8[_mode][y][x][n]=_beta_1[ij];
            NE_PRED_PARAMX_8x8[_mode][y][x][n]=B_SZ*bx+u;
            NE_PRED_PARAMY_8x8[_mode][y][x][n]=B_SZ*by+v;
#elif B_SZ==16
            NE_PRED_WEIGHTS_16x16[_mode][y][x][n]=_beta_1[ij];
            NE_PRED_PARAMX_16x16[_mode][y][x][n]=B_SZ*bx+u;
            NE_PRED_PARAMY_16x16[_mode][y][x][n]=B_SZ*by+v;
#else
# error "Need predictors for this block size."
#endif
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
  int j;
  int i;
  int k;
  fprintf(_fp,"double NE_PRED_OFFSETS_%ix%i[OD_INTRA_NMODES][%i][%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for(i=0;i<B_SZ;i++){
#if B_SZ==4
        fprintf(_fp,"%s%- 24.18G",i>0?",":" ",NE_PRED_OFFSETS_4x4[m][j][i]);
#elif B_SZ==8
        fprintf(_fp,"%s%- 24.18G",i>0?",":"",NE_PRED_OFFSETS_8x8[m][j][i]);
#elif B_SZ==16
        fprintf(_fp,"%s%- 24.18G",i>0?",":"",NE_PRED_OFFSETS_16x16[m][j][i]);
#else
# error "Unsupported block size."
#endif
      }
      fprintf(_fp," }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fprintf(_fp,"int NE_PRED_MULTS_%ix%i[OD_INTRA_NMODES][%i][%i]={\n",
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
  fprintf(_fp,
   "double NE_PRED_WEIGHTS_%ix%i[OD_INTRA_NMODES][%i][%i][4*%i*%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s}%s\n",k==0?"  0  ":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fprintf(_fp,"int NE_PRED_PARAMX_%ix%i[OD_INTRA_NMODES][%i][%i][4*%i*%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s   }%s\n",k==0?"  0":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fprintf(_fp,"int NE_PRED_PARAMY_%ix%i[OD_INTRA_NMODES][%i][%i][4*%i*%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s   }%s\n",k==0?"  0":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fflush(_fp);
}

void print_predictors_nonsparse(FILE *_fp){
  int m;
  int i;
  int j;
  int k;
  int l;
  int bx;
  int by;
  fprintf(_fp,"double NE_PRED_OFFSETS_%ix%i[OD_INTRA_NMODES][%i][%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for(i=0;i<B_SZ;i++){
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
  fprintf(_fp,
   "double NE_PRED_WEIGHTS_%ix%i[OD_INTRA_NMODES][%i][%i][2*%i][3*%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ,B_SZ,B_SZ);
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"  {\n");
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        double w[2*B_SZ][3*B_SZ];
        fprintf(_fp,"/* Mode %i (%i,%i) */\n      {\n",m,j,i);
        for(k=0;k<2*B_SZ;k++){
          for(l=0;l<3*B_SZ;l++){
            w[k][l]=0;
          }
        }
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          int u=NE_PRED_PARAMX_4x4[m][j][i][k];
          int v=NE_PRED_PARAMY_4x4[m][j][i][k];
          w[v][u]=NE_PRED_WEIGHTS_4x4[m][j][i][k];
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          int u=NE_PRED_PARAMX_8x8[m][j][i][k];
          int v=NE_PRED_PARAMY_8x8[m][j][i][k];
          w[v][u]=NE_PRED_WEIGHTS_8x8[m][j][i][k];
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          int u=NE_PRED_PARAMX_16x16[m][j][i][k];
          int v=NE_PRED_PARAMY_16x16[m][j][i][k];
          w[v][u]=NE_PRED_WEIGHTS_16x16[m][j][i][k];
        }
#else
# error "Unsupported block size."
#endif
        for(by=0;by<=1;by++){
          for(k=0;k<B_SZ;k++){
            fprintf(_fp,"        {");
            for(bx=0;bx<=(1-by)<<1;bx++){
              for(l=0;l<B_SZ;l++){
                fprintf(_fp,"%s  %- 24.18G",bx>0||l>0?",":"",
                 w[B_SZ*by+k][B_SZ*bx+l]);
              }
            }
            fprintf(_fp,"  }%s\n",by<1||k<B_SZ-1?",":"");
          }
        }
        fprintf(_fp,"      }%s\n",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
}
