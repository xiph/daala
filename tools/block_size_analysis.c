/*Daala video codec
Copyright (c) 2002-2013 Daala project contributors.  All rights reserved.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/odintrin.h"
#include "vidinput.h"
#include "../src/filter.h"
#include "../src/dct.h"
#include "../src/pvq.h"
#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#endif
#include "getopt.h"
#include "../src/block_size.h"
#include "../src/block_size_enc.h"
#include <math.h>

static void usage(char **_argv){
  fprintf(stderr,"Usage: %s [options] <input> <output>\n"
   "    <reference> and <input> may be either YUV4MPEG or Ogg Theora files.\n\n"
   "    Options:\n"
   "      --intra    Intraframes only.\n"
   "      --limit N  Only read N frames from input.\n"
   "      --ext   N  Encode the final frame N more times.\n"
   "      --fps   N  Override output fps.\n"
   "      --ref   N  Reference stream.\n"
   "      --pvqk  N  PVQ K parameter.\n",_argv[0]);
}

static const char *CHROMA_TAGS[4]={" C420jpeg",""," C422jpeg"," C444"};


/* Warning, this will fail for images larger than 2024 x 2024 */
#define MAX_VAR_BLOCKS 1024
#define SQUARE(x) ((int)(x)*(int)(x))

/* Actual 2D coding gains of lapped transforms (the 32x32 one is made-up). We divide by 6 to get bits. */
#define CG4 (15.943/6)
#define CG8 (16.7836/6)
#define CG16 (16.9986/6)
#define CG32 (17.1/6)


#define OFF8   (1)
#define OFF16  (2)

#define OFF16_8  (1)
#define OFF32_8  (1)

#define COUNT16_8  (3+2*OFF16_8)
#define COUNT32_8  (7+2*OFF32_8)


#define PSY_LAMBDA .65

int switch_decision(unsigned char *img, int w, int h, int stride, int ow, int oh)
{
  int i,j;
  int h8,w8,h32,w32;
  static unsigned char dec8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
#if 0
  int h4,w4,h16,w16;
  static int Sx[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int Sxx[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int Sx4[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];
  static int Sxx4[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];

  static int var[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int var_1[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int var8[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];
  static int var8_1[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];
  static float dummy[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static float dummy8[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];


  static float nmr4[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];
  static float nmr8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
  static float cg8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
  static float nmr16[MAX_VAR_BLOCKS>>3][MAX_VAR_BLOCKS>>3];
  static float cg16[MAX_VAR_BLOCKS>>3][MAX_VAR_BLOCKS>>3];
  static float nmr32[MAX_VAR_BLOCKS>>4][MAX_VAR_BLOCKS>>4];
  static float cg32[MAX_VAR_BLOCKS>>4][MAX_VAR_BLOCKS>>4];

  const unsigned char *x;
#endif

  (void)ow;
  (void)oh;

  w>>=1;
  h>>=1;
  w8 = w>>2;
  h8 = h>>2;
  w32 = w>>4;
  h32 = h>>4;
#if 0
  w4 = w>>1;
  h4 = h>>1;
  w16 = w>>3;
  h16 = h>>3;
  x = img;
  for(i=0;i<h;i++){
    for(j=0;j<w;j++){
      Sx[i][j]=x[2*j]+x[2*j+1]+x[stride+2*j]+x[stride+2*j+1];
      Sxx[i][j]=SQUARE(x[2*j])+SQUARE(x[2*j+1])+SQUARE(x[stride+2*j])+SQUARE(x[stride+2*j+1]);
    }
    x+=2*stride;
  }
  for(i=0;i<h4;i++){
    for(j=0;j<w4;j++){
      Sxx4[i][j] = Sxx[2*i][2*j] + Sxx[2*i][2*j+1] + Sxx[2*i+1][2*j] + Sxx[2*i+1][2*j+1];
      Sx4[i][j]  = Sx [2*i][2*j] + Sx [2*i][2*j+1] + Sx [2*i+1][2*j] + Sx [2*i+1][2*j+1];
    }
  }

  for(i=0;i<h-1;i++){
    for(j=0;j<w-1;j++){
      int sum_x;
      int sum_xx;
      int var_floor;
      sum_x=Sx[i][j]+Sx[i][j+1]+Sx[i+1][j]+Sx[i+1][j+1];
      sum_xx=Sxx[i][j]+Sxx[i][j+1]+Sxx[i+1][j]+Sxx[i+1][j+1];
      var[i][j]=(sum_xx-(SQUARE(sum_x)>>4))>>5;
      var_floor = 4+(sum_x>>8);
      if (var[i][j]<var_floor)var[i][j]=var_floor;
      /*printf("%d ", var[i][j]);*/
      var_1[i][j] = 16384/var[i][j];
    }
    /*printf("\n");*/
  }

  for(i=0;i<h4-1;i++){
    for(j=0;j<w4-1;j++){
      int sum_x;
      int sum_xx;
      int var_floor;
      sum_x =Sx4 [i][j]+Sx4 [i][j+1]+Sx4 [i+1][j]+Sx4 [i+1][j+1];
      sum_xx=Sxx4[i][j]+Sxx4[i][j+1]+Sxx4[i+1][j]+Sxx4[i+1][j+1];
      var8[i][j]=(sum_xx-(SQUARE(sum_x)>>6))>>5;
      var_floor = 4+(sum_x>>8);
      if (var8[i][j]<var_floor)var8[i][j]=var_floor;
      /*printf("%d ", var8[i][j]);*/
      var8_1[i][j] = 16384/var8[i][j];
    }
    /*printf("\n");*/
  }

  for(i=1;i<h4-1;i++){
    for(j=1;j<w4-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var_1=0;
      float psy=0;
      float noise;
      for(k=0;k<3;k++){
        for(m=0;m<3;m++){
          sum_var+=var[2*i-1+k][2*j-1+m];
        }
      }
      noise = sum_var/(3*3);
      for(k=0;k<3;k++){
        for(m=0;m<3;m++){
          psy += OD_LOG2(1+noise*var_1[2*i-1+k][2*j-1+m]/16384.);
        }
      }
      psy /= (3*3);
      psy -= 1;
      nmr4[i][j] = psy;
      /*printf("%f ", nmr4[i][j]);*/
    }
    /*printf("\n");*/
  }

  for(i=1;i<h8-1;i++){
    for(j=1;j<w8-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var_1=0;
      float nmr4_avg;
      float cgl, cgs;
      float noise;
      float psy;
      for(k=0;k<COUNT8;k++){
        for(m=0;m<COUNT8;m++){
          sum_var  +=var[4*i-OFF8+k][4*j-OFF8+m];
        }
      }
      noise = sum_var/(COUNT8*COUNT8);
      psy=0;
      for(k=0;k<COUNT8;k++){
        for(m=0;m<COUNT8;m++){
          psy += OD_LOG2(1+noise*var_1[4*i-OFF8+k][4*j-OFF8+m]/16384.);
        }
      }
      psy /= (COUNT8*COUNT8);
      psy -= 1;
      nmr8[i][j] = psy;
      nmr4_avg = .25f*(nmr4[2*i][2*j]+nmr4[2*i][2*j+1]+nmr4[2*i+1][2*j]+nmr4[2*i+1][2*j+1]);
      cgs = CG4 - PSY_LAMBDA*(nmr4_avg);
      cgl = CG8 - PSY_LAMBDA*(nmr8[i][j]);
      if (cgl>=cgs)
      {
        dec8[i][j] = 1;
        cg8[i][j] = CG8;
      } else {
        nmr8[i][j] = nmr4_avg;
        dec8[i][j] = 0;
        cg8[i][j] = CG4;
      }
      /*printf("%d ", dec8[i][j]);*/
    }
    /*printf("\n");*/
  }

  for(i=1;i<h16-1;i++){
    for(j=1;j<w16-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var8=0;
      int sum_var_1=0;
      float nmr8_avg;
      float cgl,cgs;
      float noise;
      float noise8;
      float psy;
      float psy8;
      for(k=0;k<COUNT16;k++){
        for(m=0;m<COUNT16;m++){
          sum_var+=var[8*i-OFF16+k][8*j-OFF16+m];
        }
      }
      noise = sum_var/(float)(COUNT16*COUNT16);
      for(k=0;k<COUNT16_8;k++){
        for(m=0;m<COUNT16_8;m++){
          sum_var8+=var8[4*i-OFF16_8+k][4*j-OFF16_8+m];
        }
      }
      noise8 = sum_var8/(float)(COUNT16_8*COUNT16_8);
      psy=0;
      for(k=0;k<COUNT16;k++){
        for(m=0;m<COUNT16;m++){
          psy += OD_LOG2(1+noise*var_1[8*i-OFF16+k][8*j-OFF16+m]/16384.);
        }
      }
      psy /= (COUNT16*COUNT16);
      psy -= 1;
      psy8=0;
      for(k=0;k<COUNT16_8;k++){
        for(m=0;m<COUNT16_8;m++){
          psy8 += OD_LOG2(1+noise8*var8_1[4*i-OFF16_8+k][4*j-OFF16_8+m]/16384.);
        }
      }
      psy8 /= (COUNT16_8*COUNT16_8);
      psy8 -= 1;
      psy = OD_MAXF(psy, .25*psy8);
      /*psy = .5*(psy+psy8);*/
      nmr16[i][j] = psy;
      nmr8_avg = .25f*(nmr8[2*i][2*j]+nmr8[2*i][2*j+1]+nmr8[2*i+1][2*j]+nmr8[2*i+1][2*j+1]);
      cg16[i][j] = .25*(cg8[2*i][2*j] + cg8[2*i][2*j+1] + cg8[2*i+1][2*j] + cg8[2*i+1][2*j+1]);
      cgs = cg16[i][j] - PSY_LAMBDA*(nmr8_avg);
      cgl = CG16 - PSY_LAMBDA*(nmr16[i][j]);
      /*printf("%f ", psy);*/
      if (cgl>=cgs)
      {
        dec8[2*i][2*j] = 2;
        dec8[2*i][2*j+1] = 2;
        dec8[2*i+1][2*j] = 2;
        dec8[2*i+1][2*j+1] = 2;
        cg16[i][j] = CG16;
      } else {
        nmr16[i][j] = nmr8_avg;
      }
    }
    /*printf("\n");*/
  }

#if 1
  for(i=1;i<h32-1;i++){
    for(j=1;j<w32-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var_1=0;
      int sum_var8=0;
      float nmr16_avg;
      float cgl,cgs;
      float noise, psy;
      float noise8, psy8;
      for(k=0;k<COUNT32;k++){
        for(m=0;m<COUNT32;m++){
          sum_var  +=var[16*i-OFF32+k][16*j-OFF32+m];
        }
      }
      noise = sum_var/(float)(COUNT32*COUNT32);
      for(k=0;k<COUNT32_8;k++){
        for(m=0;m<COUNT32_8;m++){
          sum_var8+=var8[8*i-OFF32_8+k][8*j-OFF32_8+m];
        }
      }
      noise8 = sum_var8/(float)(COUNT32_8*COUNT32_8);
      psy=0;
      for(k=0;k<COUNT32;k++){
        for(m=0;m<COUNT32;m++){
          psy += OD_LOG2(1.+noise*var_1[16*i-OFF32+k][16*j-OFF32+m]/16384.);
        }
      }
      psy /= (COUNT32*COUNT32);
      psy -= 1;
      psy8=0;
      for(k=0;k<COUNT32_8;k++){
        for(m=0;m<COUNT32_8;m++){
          psy8 += OD_LOG2(1+noise8*var8_1[8*i-OFF32_8+k][8*j-OFF32_8+m]/16384.);
        }
      }
      psy8 /= (COUNT32_8*COUNT32_8);
      psy8 -= 1;
      psy = OD_MAXF(psy, .25*psy8);
      /*psy = .5*(psy+psy8);*/
     /*psy += psy8;*/
      nmr32[i][j] = psy;
      nmr16_avg = .25f*(nmr16[2*i][2*j]+nmr16[2*i][2*j+1]+nmr16[2*i+1][2*j]+nmr16[2*i+1][2*j+1]);
      cg32[i][j] = .25*(cg16[2*i][2*j] + cg16[2*i][2*j+1] + cg16[2*i+1][2*j] + cg16[2*i+1][2*j+1]);
      cgs = cg32[i][j] - PSY_LAMBDA*(nmr16_avg);
      cgl = CG32 - PSY_LAMBDA*(nmr32[i][j]);
      /*printf("%f ", psy8);*/
      if (cgl>=cgs)
      {
        for(k=0;k<4;k++){
          for(m=0;m<4;m++){
            dec8[4*i+k][4*j+m]=3;
          }
        }
        cg32[i][j] = CG32;
      } else {
        nmr32[i][j] = nmr16_avg;
      }
    }
    /*printf("\n");*/
  }
#endif
#endif
  /* Replace decision with the one from process_block_size32() */
  if (1)
  {
    BlockSizeComp bs;

    for(i=1;i<h32-1;i++){
      for(j=1;j<w32-1;j++){
        int k,m;
        int dec[4][4];
        process_block_size32(&bs, img+32*stride*i+32*j, stride, NULL, 0, dec);
        for(k=0;k<4;k++)
          for(m=0;m<4;m++)
            dec8[4*i+k][4*j+m]=dec[k][m];
#if 0
        for(k=0;k<16;k++)
        {
          for(m=0;m<16;m++)
          {
            var[16*i+k][16*j+m]=bs.img_stats.Var4[k+3][m+3];
            var_1[16*i+k][16*j+m]=bs.img_stats.invVar4[k+3][m+3];
          }
        }
        for(k=0;k<8;k++)
        {
          for(m=0;m<8;m++)
          {
            dummy[8*i+k][8*j+m]=bs.psy4[k][m];
          }
        }
        for(k=0;k<4;k++)
        {
          for(m=0;m<4;m++)
          {
            dummy8[4*i+k][4*j+m]=bs.psy8[k][m];
          }
        }
        for(k=0;k<8;k++)
        {
          for(m=0;m<8;m++)
          {
            var8[8*i+k][8*j+m]=bs.img_stats.Var8[k+2][m+2];
            var8_1[8*i+k][8*j+m]=bs.img_stats.invVar8[k+2][m+2];
          }
        }
#endif
      }
    }
    if (0) {
      int count32[28] = {0};
      int count16[144] = {0};
      int stats32[28] = {0};
      int stats16[144] = {0};
      int stats8[144][16] = {{0}};
      for (i = 2; i < h32 - 2; i++) {
        for (j = 2; j < w32 - 2; j++) {
          int k;
          int m;
          unsigned char *bsize;
          int stride;
          int id32;
          bsize = &dec8[4*i][4*j];
          stride = sizeof(dec8[0])/sizeof(dec8[0][0]);
          id32 = od_block_size_prob32(bsize, stride);
          count32[id32]++;
          if (bsize[0] == 3) stats32[id32]++;
          else for (k = 0; k < 2; k++) for (m = 0; m < 2; m++) {
            int id16;
            id16 = od_block_size_cdf16_id(&bsize[2*k*stride + 2*m], stride);
            count16[id16]++;
            OD_ASSERT(bsize[2*k*stride + 2*m]<=2);
            stats16[id16] += bsize[2*k*stride + 2*m] == 2;
            if (bsize[2*k*stride + 2*m] != 2) {
              int id8;
              OD_ASSERT(bsize[2*k*stride + 2*m]<=1);
              OD_ASSERT(bsize[2*k*stride + 2*m + 1]<=1);
              OD_ASSERT(bsize[2*k*stride + stride + 2*m]<=1);
              OD_ASSERT(bsize[2*k*stride + stride + 2*m + 1]<=1);
              id8 = 8*bsize[2*k*stride + 2*m] + 4*bsize[2*k*stride + 2*m + 1]
               + 2*bsize[2*k*stride + stride + 2*m]
               + bsize[2*k*stride + stride + 2*m + 1];
              stats8[id16][id8]++;
            }
          }
        }
      }
      for (i = 0; i < 28; i++) printf("%d ", count32[i]);
      for (i = 0; i < 28; i++) printf("%d ", stats32[i]);
      for (i = 0; i < 144; i++) printf("%d ", count16[i]);
      for (i = 0; i < 144; i++) printf("%d ", stats16[i]);
      for (i = 0; i < 144; i++) for (j = 0; j < 16; j++) printf("%d ", stats8[i][j]);
      printf("\n");
    }
  }

#if 0
  for(i=0;i<h;i++){
    for(j=0;j<w;j++){
      printf("%d ", var[i][j]);
    }
    printf("\n");
  }
#endif
#if 0
  for(i=8;i<h4-8;i++){
    for(j=8;j<w4-8;j++){
      printf("%f ", dummy[i][j]);
    }
    printf("\n");
  }
#endif
#if 0
  for(i=4;i<h8-4;i++){
    for(j=4;j<w8-4;j++){
      printf("%f ", dummy8[i][j]);
    }
    printf("\n");
  }
#endif

#if 0
  for(i=4;i<h8-4;i++){
    for(j=4;j<w8-4;j++){
      printf("%d ", dec8[i][j]);
    }
    printf("\n");
  }
#endif

#if 0
  fprintf(stderr, "size : %dx%d\n", (w<<1), (h<<1));
  for(i=0;i<(h<<1);i++){
    for(j=0;j<1296;j++){
      putc(dec8[i>>3][j>>3], stdout);
    }
  }
#endif
#if 0
  {
    /*Raw mode data with offsets to match the 4x8 training tool's padding.*/
    int wn=(ow-16)>>2;
    int hn=(oh-16)>>2;
    fprintf(stderr, "size : %dx%d\n", wn, hn);
    for(i=0;i<hn;i++){
      for(j=0;j<wn;j++){
        int posi=(i>>1)+1;
        int posj=(j>>1)+1;
        if(posi>=4 && posi<(h8-4) &&posj>=4 && posj<(w8-4))putc(dec8[posi][posj], stdout);
        else putc(0, stdout);
      }
    }
  }
#endif
#if 1
  for(i=4;i<h8-4;i++){
    for(j=4;j<w8-4;j++){
      if ((i&3)==0 && (j&3)==0){
        int k;
        for(k=0;k<32;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<32;k++)
          img[(8*i+k)*stride+j*8] = 0;
      }
      if ((i&1)==0 && (j&1)==0 && dec8[i][j]==2){
        int k;
        for(k=0;k<16;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<16;k++)
          img[(8*i+k)*stride+j*8] = 0;
      }
      if (dec8[i][j]<=1){
        int k;
        for(k=0;k<8;k++)
          img[i*stride*8+j*8+k] = 0;
        for(k=0;k<8;k++)
          img[(8*i+k)*stride+j*8] = 0;
        if (dec8[i][j]==0){
          img[(8*i+4)*stride+j*8+3] = 0;
          img[(8*i+4)*stride+j*8+4] = 0;
          img[(8*i+4)*stride+j*8+5] = 0;
          img[(8*i+3)*stride+j*8+4] = 0;
          img[(8*i+5)*stride+j*8+4] = 0;
        }
      }
    }
  }
  for (i=32;i<(w32-1)*32;i++)
    img[(h32-1)*32*stride+i]=0;
  for (i=32;i<(h32-1)*32;i++)
    img[i*stride+(w32-1)*32]=0;
#endif

#if 0 /* 32x32 decision data */
  for(i=2;i<h32-2;i++){
    for(j=2;j<w32-2;j++){
      int i8, j8;
      int k;
      i8 = 4*i-1;
      j8 = 4*j-1;
      for(k=0;k<5;k++)
        printf("%d ", dec8[i8+k][j8]);
      for(k=1;k<5;k++)
        printf("%d ", dec8[i8][j8+k]);
      printf("%d\n", dec8[i8+1][j8+1]);
    }
  }
#endif

#if 0 /* 16x16 decision data */
  for(i=4;i<h16-4;i++){
    for(j=4;j<w16-4;j++){
      int i8, j8;
      int k;
      i8 = 2*i-1;
      j8 = 2*j-1;
      if (dec8[i8+1][j8+1]==3)
        continue;
      if (1)
      {
        int sum=0;
        /*for(k=0;k<3;k++)
          printf("%d ", dec8[i8+k][j8]);
        for(k=1;k<3;k++)
          printf("%d ", dec8[i8][j8+k]);*/
        for(k=0;k<3;k++)
          sum += dec8[i8+k][j8];
        for(k=1;k<3;k++)
          sum += dec8[i8][j8+k];
        printf("%d ", sum);
      } else {
        int up, left;
        up = (dec8[i8][j8+1]>=2) ? 2+dec8[i8][j8+1] : dec8[i8][j8+1]*2+dec8[i8][j8+2];
        left = (dec8[i8+1][j8]>=2) ? 2+dec8[i8+1][j8] : dec8[i8+1][j8]*2+dec8[i8+2][j8];
        printf("%d ", dec8[i8][j8]*36+up*6+left);
      }
      if (dec8[i8+1][j8+1]==2)
        printf("16\n");
      else
        printf("%d\n", 8*dec8[i8+1][j8+1]+4*dec8[i8+1][j8+2]+2*dec8[i8+2][j8+1]+dec8[i8+2][j8+2]);
    }
    printf("\n");
  }
#endif

#if 0 /* 8x8 decision data */
  for(i=8;i<h8-8;i++){
    for(j=8;j<w8-8;j++){
      if (dec8[i][j]<=1)
      {
        printf("%d %d %d %d\n", dec8[i-1][j-1], dec8[i-1][j], dec8[i][j-1], dec8[i][j]);
      }
    }
  }
#endif

  return 0;
}

/*Applies vert then horiz prefilters of size _n.*/
void prefilter_image(od_coeff *_img, int _w, int _h, int _n){
  int x,y,j;
  /*Pre-filter*/
  switch(_n){
    case 4:
      for(y=0;y<_h;y++){
        for(x=2;x<_w-2;x+=4){
          od_pre_filter4(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      for(y=2;y<_h-2;y+=4){
        for(x=0;x<_w;x++){
          od_coeff tmp[4];
          for(j=0;j<4;j++)tmp[j]=_img[(y+j)*_w+x];
          od_pre_filter4(tmp,tmp);
          for(j=0;j<4;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      break;
    case 8:
      for(y=0;y<_h;y++){
        for(x=4;x<_w-4;x+=8){
          od_pre_filter8(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      for(y=4;y<_h-4;y+=8){
        for(x=0;x<_w;x++){
          od_coeff tmp[16];
          for(j=0;j<8;j++)tmp[j]=_img[(y+j)*_w+x];
          od_pre_filter8(tmp,tmp);
          for(j=0;j<8;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      break;
    case 16:
      for(y=0;y<_h;y++){
        for(x=8;x<_w-8;x+=16){
          od_pre_filter16(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      for(y=8;y<_h-8;y+=16){
        for(x=0;x<_w;x++){
          od_coeff tmp[16];
          for(j=0;j<16;j++)tmp[j]=_img[(y+j)*_w+x];
          od_pre_filter16(tmp,tmp);
          for(j=0;j<16;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      break;
    default:
      break;
  }
}

/*Applies horiz then vert postfilters of size _n.*/
void postfilter_image(od_coeff *_img, int _w, int _h, int _n){
  int x,y,j;
  /*Pre-filter*/
  switch(_n){
    case 4:
      for(y=2;y<_h-2;y+=4){
        for(x=0;x<_w;x++){
          od_coeff tmp[4];
          for(j=0;j<4;j++)tmp[j]=_img[(y+j)*_w+x];
          od_post_filter4(tmp,tmp);
          for(j=0;j<4;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      for(y=0;y<_h;y++){
        for(x=2;x<_w-2;x+=4){
          od_post_filter4(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      break;
    case 8:
      for(y=4;y<_h-4;y+=8){
        for(x=0;x<_w;x++){
          od_coeff tmp[16];
          for(j=0;j<8;j++)tmp[j]=_img[(y+j)*_w+x];
          od_post_filter8(tmp,tmp);
          for(j=0;j<8;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      for(y=0;y<_h;y++){
        for(x=4;x<_w-4;x+=8){
          od_post_filter8(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      break;
    case 16:
      for(y=8;y<_h-8;y+=16){
        for(x=0;x<_w;x++){
          od_coeff tmp[16];
          for(j=0;j<16;j++)tmp[j]=_img[(y+j)*_w+x];
          od_post_filter16(tmp,tmp);
          for(j=0;j<16;j++)_img[(y+j)*_w+x]=tmp[j];
        }
      }
      for(y=0;y<_h;y++){
        for(x=8;x<_w-8;x+=16){
          od_post_filter16(&_img[y*_w+x],&_img[y*_w+x]);
        }
      }
      break;
    default:
      break;
  }
}

/*static const int fzig16[256] = {0,16,1,2,17,32,48,33,18,3,4,19,34,49,64,80,65,50,35,20,5,6,21,36,51,66,81,96,112,97,82,67,52,37,22,7,8,23,38,53,68,83,98,113,128,144,129,114,99,84,69,54,39,24,9,10,25,40,55,70,85,100,115,130,145,160,176,161,146,131,116,101,86,71,56,41,26,11,12,27,42,57,72,87,102,117,132,147,162,177,192,208,193,178,163,148,133,118,103,88,73,58,43,28,13,14,29,44,59,74,89,104,119,134,149,164,179,194,209,224,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,31,46,61,76,91,106,121,136,151,166,181,196,211,226,241,242,227,212,197,182,167,152,137,122,107,92,77,62,47,63,78,93,108,123,138,153,168,183,198,213,228,243,244,229,214,199,184,169,154,139,124,109,94,79,95,110,125,140,155,170,185,200,215,230,245,246,231,216,201,186,171,156,141,126,111,127,142,157,172,187,202,217,232,247,248,233,218,203,188,173,158,143,159,174,189,204,219,234,249,250,235,220,205,190,175,191,206,221,236,251,252,237,222,207,223,238,253,254,239,255};
static const int izig16[256] = {0,2,3,9,10,20,21,35,36,54,55,77,78,104,105,135,1,4,8,11,19,22,34,37,53,56,76,79,103,106,134,136,5,7,12,18,23,33,38,52,57,75,80,102,107,133,137,164,6,13,17,24,32,39,51,58,74,81,101,108,132,138,163,165,14,16,25,31,40,50,59,73,82,100,109,131,139,162,166,189,15,26,30,41,49,60,72,83,99,110,130,140,161,167,188,190,27,29,42,48,61,71,84,98,111,129,141,160,168,187,191,210,28,43,47,62,70,85,97,112,128,142,159,169,186,192,209,211,44,46,63,69,86,96,113,127,143,158,170,185,193,208,212,227,45,64,68,87,95,114,126,144,157,171,184,194,207,213,226,228,65,67,88,94,115,125,145,156,172,183,195,206,214,225,229,240,66,89,93,116,124,146,155,173,182,196,205,215,224,230,239,241,90,92,117,123,147,154,174,181,197,204,216,223,231,238,242,249,91,118,122,148,153,175,180,198,203,217,222,232,237,243,248,250,119,121,149,152,176,179,199,202,218,221,233,236,244,247,251,254,120,150,151,177,178,200,201,219,220,234,235,245,246,252,253,255};
*/

#define ROUNDUP_32(x) (((x)+31)&~31)



#define MAXB 16
#define SQUARE(x) ((int)(x)*(int)(x))

int oc_ilog32(unsigned _v){
  int ret;
  static const unsigned char OC_DEBRUIJN_IDX32[32]={
     0, 1,28, 2,29,14,24, 3,30,22,20,15,25,17, 4, 8,
    31,27,13,23,21,19,16, 7,26,12,18, 6,11, 5,10, 9};
  _v|=_v>>1;
  _v|=_v>>2;
  _v|=_v>>4;
  _v|=_v>>8;
  _v|=_v>>16;
  ret=_v&1;
  _v=(_v>>1)+1;
  ret+=OC_DEBRUIJN_IDX32[_v*0x77CB531U>>27&0x1F];
  return ret;

}



void quant_scalar_gain(ogg_int32_t *_x,ogg_int16_t *_scale,int *y,int N,int Q){
  float gain0, gain1;
  float Q_1;
  int i;

  (void)_scale;

  Q*=15;
  Q_1 = 1.f/Q;
  gain0=0;
  gain1=0;
  for (i=0;i<N;i++)
  {
    int qx;
    float s = _x[i]*Q_1;
    float bias = s>0?-.49:.49;
    gain0 += s*s;
    qx = (int)floor(.5+s+bias);
    y[i] = qx;
    gain1 += qx*qx;
    _x[i] = Q*qx;
  }
  gain0 = sqrt(gain0/(1e-15+gain1));
  for (i=0;i<N;i++)
    _x[i] *= gain0;
}


static void process_plane(od_coeff *_img, od_coeff *_refi, int _w, int _h, int _pli, int _pvq_k){
  int x;
  int y;
  int j;
  int i;
  int free_ref;
  static int count=0;
  (void)_pvq_k;
  _w = ROUNDUP_32(_w);
  _h = ROUNDUP_32(_h);
  if(!_refi){
    _refi=calloc(ROUNDUP_32(_w)*ROUNDUP_32(_h),sizeof(od_coeff));
    free_ref=1;
  }else free_ref=0;

  prefilter_image(_img,_w,_h,16);


  /*for (i=0;i<1000000;i++){
    int tmp[16];
    for(j=0;j<16;j++)
      tmp[j] = rand()%255-127;
    od_bin_idct16(tmp, tmp);
    for(j=0;j<16;j++)printf("%d ", tmp[j]);
    printf("\n");
  }
  exit(0);*/
  /*Block processing.*/
  for(y=0;y<_h;y+=16){
    for(x=0;x<_w;x+=16){
      od_coeff coeffs[256];
      ogg_int32_t zi[256];
/*      int         out[256];*/
/*      ogg_int16_t scale[256];*/

/*      for(j=0;j<256;j++)scale[j]=1;*/
      od_bin_fdct16x16(coeffs,16,&_img[(y)*_w+x],_w);

      for(i=0;i<16;i++){
        for(j=0;j<16;j++){
          zi[16*i+j]=floor(.5+coeffs[16*i+j]);
        }
      }
      if (_pli==-1){
#if 0
        int foo[256];
        ogg_int32_t x[256];
        /*quant_scalar_gain(&zi[1],NULL,foo,255,200);*/
        extract(&zi[4], x, 2, 4, 16);
        quant_scalar_gain(x,NULL,foo,8,200);
        interleave(x, &zi[4], 2, 4, 16);

        extract(&zi[64], x, 4, 2, 16);
        quant_scalar_gain(x,NULL,foo,8,200);
        interleave(x, &zi[64], 4, 2, 16);

        extract(&zi[64+2], x, 4, 6, 16);
        extract(&zi[32+4], x+24, 2, 4, 16);
        quant_scalar_gain(x,NULL,foo,32,600);
        interleave(x, &zi[64+2], 4, 6, 16);
        interleave(x+24, &zi[32+4], 2, 4, 16);
#endif
#if 0
        extract(&zi[8], x, 4, 8, 16);
        /*quant_pvq(x, r, scale, out, 32, 1200./f, &qg);*/
        quant_scalar_gain(x,NULL,foo,32,600);
        interleave(x, &zi[8], 4, 8, 16);

        extract(&zi[128], x, 8, 4, 16);
        /*quant_pvq(x, r, scale, out, 32, 1200./f, &qg);*/
        quant_scalar_gain(x,NULL,foo,32,600);
        interleave(x, &zi[128], 8, 4, 16);

        extract(&zi[128+4], x, 8, 12, 16);
        extract(&zi[64+8], x+96, 4, 8, 16);
        /*quant_pvq(x, r, scale, out, 128, 1200./f, &qg);*/
        quant_scalar_gain(x,NULL,foo,128,1200);
        interleave(x, &zi[128+4], 8, 12, 16);
        interleave(x+96, &zi[64+8], 4, 8, 16);
#endif
      }
      /*for(j=0;j<256;j++)coeffs[j]=zi[j];*/
      for(i=0;i<16;i++){
        for(j=0;j<16;j++){
          coeffs[16*i+j]=floor(.5+zi[16*i+j]);
        }
      }
      od_bin_idct16x16(&_img[(y)*_w+x],_w,coeffs,16);
    }
    /*printf("\n");*/
  }
  postfilter_image(_img,_w,_h,16);
  if(free_ref)free(_refi);
  count++;
}

int main(int _argc,char **_argv){
  const char *optstring = "hv?";
  const struct option long_options[]={
    {"ref",required_argument,NULL,0},
    {"limit",required_argument,NULL,0},
    {"fps",required_argument,NULL,0},
    {"ext",required_argument,NULL,0},
    {"pvqk",required_argument,NULL,0},
    {"intra",no_argument,NULL,0},
    {NULL,0,NULL,0}
  };
  FILE *fin;
  FILE *fout;
  video_input vid1;
  video_input_info info1;
  video_input vid2;
  video_input_info info2;
  int frameno;
  int pli;
  unsigned char *outline;
  od_coeff *refi[3];
  od_coeff *iimg[3];
  int xdec[3];
  int ydec[3];
  int w[3];
  int h[3];
  int pvq_k;
  int fps;
  int extend;
  int limit;
  int intra;
  char refname[1024];
  int ref_in;
  int long_option_index;
  int c;
  pvq_k=32;
  fps=-1;
  ref_in=0;
  limit=0;
  extend=0;
  intra=0;
  while((c=getopt_long(_argc,_argv,optstring,long_options,&long_option_index))!=EOF){
    switch(c){
      case 0:
        if(strcmp(long_options[long_option_index].name,"ref")==0){
          ref_in=1;
          strncpy(refname,optarg,1023);
        } else if (strcmp(long_options[long_option_index].name,"pvqk")==0){
          pvq_k=atoi(optarg);
        } else if (strcmp(long_options[long_option_index].name,"limit")==0){
          limit=atoi(optarg);
        } else if (strcmp(long_options[long_option_index].name,"ext")==0){
          extend=atoi(optarg);
        } else if (strcmp(long_options[long_option_index].name,"intra")==0){
          intra=1;
        } else if (strcmp(long_options[long_option_index].name,"fps")==0){
          fps=atoi(optarg);
        }
      break;
      case 'v':
      case '?':
      case 'h':
      default:{
        usage(_argv);
        exit(EXIT_FAILURE);
      }break;
    }
  }
  if(optind+2!=_argc){
    usage(_argv);
    exit(EXIT_FAILURE);
  }
  fin=strcmp(_argv[optind],"-")==0?stdin:fopen(_argv[optind],"rb");
  if(fin==NULL){
    fprintf(stderr,"Unable to open '%s' for extraction.\n",_argv[optind]);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr,"Opening %s as input%s...\n",_argv[optind],ref_in?"":" and reference");
  if(video_input_open(&vid1,fin)<0)exit(EXIT_FAILURE);
  video_input_get_info(&vid1,&info1);
  if(ref_in){
    fin=fopen(refname,"rb");
    if(fin==NULL){
      fprintf(stderr,"Unable to open '%s' for extraction.\n",refname);
      exit(EXIT_FAILURE);
    }
    fprintf(stderr,"Opening %s as reference...\n",refname);
    if(video_input_open(&vid2,fin)<0)exit(EXIT_FAILURE);
    video_input_get_info(&vid2,&info2);
    /*Check to make sure these videos are compatible.*/
    if(info1.pic_w!=info2.pic_w||info1.pic_h!=info2.pic_h){
      fprintf(stderr,"Video resolution does not match.\n");
      exit(EXIT_FAILURE);
    }
    if(info1.pixel_fmt!=info2.pixel_fmt){
      fprintf(stderr,"Pixel formats do not match.\n");
      exit(EXIT_FAILURE);
    }
    if((info1.pic_x&!(info1.pixel_fmt&1))!=(info2.pic_x&!(info2.pixel_fmt&1))||
     (info1.pic_y&!(info1.pixel_fmt&2))!=(info2.pic_y&!(info2.pixel_fmt&2))){
      fprintf(stderr,"Chroma subsampling offsets do not match.\n");
      exit(EXIT_FAILURE);
    }
    if(info1.fps_n*(ogg_int64_t)info2.fps_d!=
     info2.fps_n*(ogg_int64_t)info1.fps_d){
      fprintf(stderr,"Warning: framerates do not match.\n");
    }
    if(info1.par_n*(ogg_int64_t)info2.par_d!=
     info2.par_n*(ogg_int64_t)info1.par_d){
     fprintf(stderr,"Warning: aspect ratios do not match.\n");
    }
  }
  for(pli=0;pli<3;pli++){
    /*Planes padded up to a multiple of 32px*/
    xdec[pli]=pli&&!(info1.pixel_fmt&1);
    ydec[pli]=pli&&!(info1.pixel_fmt&2);
    h[pli]=ROUNDUP_32(info1.pic_h>>ydec[pli]);
    w[pli]=ROUNDUP_32(info1.pic_w>>xdec[pli]);
    refi[pli]=malloc(w[pli]*h[pli]*sizeof(od_coeff));
    iimg[pli]=malloc(w[pli]*h[pli]*sizeof(od_coeff));
  }
  outline=malloc(sizeof(*outline)*info1.pic_w);
  fout=strcmp(_argv[optind+1],"-")==0?stdout:fopen(_argv[optind+1],"wb");
  if(fout==NULL){
    fprintf(stderr,"Error opening output file \"%s\".\n",_argv[optind+1]);
    return 1;
  }
  fprintf(fout,"YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   info1.pic_w,info1.pic_h, fps > 0 ? (unsigned) fps : (unsigned) info1.fps_n,
   fps > 0 ? 1U : (unsigned) info1.fps_d, info1.par_n, info1.par_d,
   CHROMA_TAGS[ydec[1] ? xdec[1] ? 0 : 2 : 3]);

  for(frameno=0;;frameno++){
    video_input_ycbcr in;
    video_input_ycbcr ref;
    int             ret1=0;
    int             ret2=0;
    char            tag1[5];
    char            tag2[5];
    if(!limit||frameno<limit){
      ret1=video_input_fetch_frame(&vid1,in,tag1);
      if(ref_in)ret2=video_input_fetch_frame(&vid2,ref,tag2);
    }
    if(ret1==0){
      if(extend<1)break;
      extend--;
      /*If we're extending, keep feeding back the output to the reference input.*/
      for(pli=0;pli<3;pli++){
        int x;
        int y;
        for(y=0;y<h[pli];y++){
          for(x=0;x<w[pli];x++){
            refi[pli][y*w[pli]+x]=iimg[pli][y*w[pli]+x];
          }
        }
      }
    }
    if(ref_in&&ret1!=0&&ret2==0){
      fprintf(stderr,"Warning: Reference ended before input!\n");
      break;
    }
    for(pli=0;pli<3;pli++){
      int x;
      int y;
      if (pli==0)
        switch_decision(in[pli].data, w[pli], h[pli], in[pli].stride, info1.pic_w, info1.pic_h);

      for(y=0;y<h[pli];y++){
        for(x=0;x<w[pli];x++){
          int cy=OD_MINI(y+(int)(info1.pic_y>>ydec[pli]),(int)info1.pic_h>>ydec[pli]);
          int cx=OD_MINI(x+(int)(info1.pic_x>>xdec[pli]),(int)info1.pic_w>>xdec[pli]);
          iimg[pli][y*w[pli]+x]=128*(in[pli].data[cy*in[pli].stride+cx]-128);
        }
      }
      if(ref_in&&ret2!=0){
        for(y=0;y<h[pli];y++){
          for(x=0;x<w[pli];x++){
            int cy=OD_MINI(y+(int)(info1.pic_y>>ydec[pli]),(int)info1.pic_h>>ydec[pli]);
            int cx=OD_MINI(x+(int)(info1.pic_x>>xdec[pli]),(int)info1.pic_w>>xdec[pli]);
            refi[pli][y*w[pli]+x]=128*(ref[pli].data[cy*in[pli].stride+cx]-128);
          }
        }
      }
      process_plane(iimg[pli],(ref_in||frameno>0)&&!intra?refi[pli]:NULL,info1.pic_w>>xdec[pli],info1.pic_h>>ydec[pli],pli,pvq_k);
      if(!ref_in){
        for(y=0;y<h[pli];y++){
          for(x=0;x<w[pli];x++){
            refi[pli][y*w[pli]+x]=iimg[pli][y*w[pli]+x];
          }
        }
      }
    }
    fprintf(fout,"FRAME\n");
    for(pli=0;pli<3;pli++){
      int x;
      int y;
      for(y=0;y<(int)info1.pic_h>>ydec[pli];y++){
        for(x=0;x<(int)info1.pic_w>>xdec[pli];x++)outline[x]=OD_CLAMP255((int)floor(.5+(1./128)*iimg[pli][y*w[pli]+x])+128);
        if(fwrite(outline,
         (info1.pic_w>>xdec[pli]),1,fout)<1){
          fprintf(stderr,"Error writing to output.\n");
          return EXIT_FAILURE;
        }
      }
    }
    fprintf(stderr, "Completed frame %d.\n",frameno);
  }
  video_input_close(&vid1);
  if(ref_in)video_input_close(&vid2);
  if(fout!=stdout)fclose(fout);
  free(outline);
  for(pli=0;pli<3;pli++)free(refi[pli]);
  for(pli=0;pli<3;pli++)free(iimg[pli]);
  return EXIT_SUCCESS;
}
