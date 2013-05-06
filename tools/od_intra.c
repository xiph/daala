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

#include "od_intra.h"

static void ne_intra_pred4x4_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      _p[4*j+i]=NE_PRED_OFFSETS_4x4[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_4x4[_mode][j][i];k++){
        x=NE_PRED_PARAMX_4x4[_mode][j][i][k];
        y=NE_PRED_PARAMY_4x4[_mode][j][i][k];
        _p[4*j+i]+=_c[_stride*(y-4)+(x-4)]*NE_PRED_WEIGHTS_4x4[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred8x8_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
      _p[8*j+i]=NE_PRED_OFFSETS_8x8[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_8x8[_mode][j][i];k++){
        x=NE_PRED_PARAMX_8x8[_mode][j][i][k];
        y=NE_PRED_PARAMY_8x8[_mode][j][i][k];
        _p[8*j+i]+=_c[_stride*(y-8)+(x-8)]*NE_PRED_WEIGHTS_8x8[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred16x16_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
      _p[16*j+i]=NE_PRED_OFFSETS_16x16[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_16x16[_mode][j][i];k++){
        x=NE_PRED_PARAMX_16x16[_mode][j][i][k];
        y=NE_PRED_PARAMY_16x16[_mode][j][i][k];
        _p[16*j+i]+=
         _c[_stride*(y-16)+(x-16)]*NE_PRED_WEIGHTS_16x16[_mode][j][i][k];
      }
    }
  }
}

const od_intra_mult_func NE_INTRA_MULT[OD_NBSIZES]={
  ne_intra_pred4x4_mult,
  ne_intra_pred8x8_mult,
  ne_intra_pred16x16_mult
};

void print_betas(FILE *_fp){
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
