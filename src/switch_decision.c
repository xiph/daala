
#include <stdio.h>
#include <math.h>
#include "switch_decision.h"

/* Warning, this will fail for images larger than 2024 x 2024 */
#define MAX_VAR_BLOCKS 1024
#define SQUARE(x) ((int)(x)*(int)(x))

#define CG4 8.6
#define CG8 9.57
#define CG16 9.81
#define CG32 9.93

/* Fudge factor giving more importance to NMR compared to coding gain.
   Partially compensates for the fact that edges have a lower coding gain */
#define fudge 1.5

int switch_decision(const unsigned char *img, int w, int h, int stride)
{
  int i,j;
  int h4, w4,h8,w8,h16,w16,h32,w32;
  static int Sx[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int Sxx[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int var[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static int var_1[MAX_VAR_BLOCKS][MAX_VAR_BLOCKS];
  static float nmr4[MAX_VAR_BLOCKS>>1][MAX_VAR_BLOCKS>>1];
  static float nmr8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
  static float cg8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
  static float nmr16[MAX_VAR_BLOCKS>>3][MAX_VAR_BLOCKS>>3];
  static float cg16[MAX_VAR_BLOCKS>>3][MAX_VAR_BLOCKS>>3];
  static float nmr32[MAX_VAR_BLOCKS>>4][MAX_VAR_BLOCKS>>4];
  static float cg32[MAX_VAR_BLOCKS>>4][MAX_VAR_BLOCKS>>4];

  static int dec8[MAX_VAR_BLOCKS>>2][MAX_VAR_BLOCKS>>2];
  const unsigned char *x;
  w>>=1;
  h>>=1;
  w4 = w>>1;
  h4 = h>>1;
  w8 = w>>2;
  h8 = h>>2;
  w16 = w>>3;
  h16 = h>>3;
  w32 = w>>4;
  h32 = h>>4;
  x = img;
  for(i=0;i<h;i++){
    for(j=0;j<w;j++){
      Sx[i][j]=x[2*j]+x[2*j+1]+x[stride+2*j]+x[stride+2*j+1];
      Sxx[i][j]=SQUARE(x[2*j])+SQUARE(x[2*j+1])+SQUARE(x[stride+2*j])+SQUARE(x[stride+2*j+1]);
    }
    x+=2*stride;
  }

  for(i=0;i<h-1;i++){
    for(j=0;j<w-1;j++){
      int sum_x;
      int sum_xx;
      int var_floor;
      sum_x=Sx[i][j]+Sx[i][j+1]+Sx[i+1][j]+Sx[i+1][j+1];
      sum_xx=Sxx[i][j]+Sxx[i][j+1]+Sxx[i+1][j]+Sxx[i+1][j+1];
      var[i][j]=(sum_xx-(SQUARE(sum_x)>>4))>>5;
      var_floor = 4+(sum_x>>6);
      if (var[i][j]<var_floor)var[i][j]=var_floor;
      /*printf("%d ", var[i][j]);*/
      var_1[i][j] = 16384/var[i][j];
    }
    /*printf("\n");*/
  }

  for(i=1;i<h4-1;i++){
    for(j=1;j<w4-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var_1=0;
      float maxvar = 0;
      for(k=0;k<3;k++){
        for(m=0;m<3;m++){
          if (var[2*i-1+k][2*j-1+m]>maxvar)maxvar=var[2*i-1+k][2*j-1+m];
          sum_var+=var[2*i-1+k][2*j-1+m];
          sum_var_1+=var_1[2*i-1+k][2*j-1+m];
        }
      }
      //sum_var=sqrt(sum_var*maxvar*9);
      nmr4[i][j] = (float)sum_var*(float)sum_var_1/81.;
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
      float maxvar = 0;
      for(k=0;k<7;k++){
        for(m=0;m<7;m++){
          if (var[4*i-2+k][4*j-2+m]>maxvar)maxvar=var[4*i-2+k][4*j-2+m];
          sum_var+=var[4*i-2+k][4*j-2+m];
          sum_var_1+=var_1[4*i-2+k][4*j-2+m];
        }
      }
      //sum_var=sqrt(sum_var*maxvar*7*7);
      nmr8[i][j] = (float)sum_var*(float)sum_var_1/(7*7*7*7);
      nmr4_avg = .25f*(nmr4[2*i][2*j]+nmr4[2*i][2*j+1]+nmr4[2*i+1][2*j]+nmr4[2*i+1][2*j+1]);
      //if (nmr4_avg>nmr8[i][j])nmr4_avg=nmr8[i][j];
      cgs = CG4/6. - fudge*.5*log2(nmr4_avg);
      cgl = CG8/6. - fudge*.5*log2(nmr8[i][j]);
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
      int sum_var_1=0;
      float nmr8_avg;
      float cgl,cgs;
      float bias=0;
      float maxvar = 0;

      /* Bias against 16x16 is there's a block that's better with 4x4 than 8x8 */
      if (!(dec8[2*i][2*j]&&dec8[2*i][2*j+1]&&dec8[2*i+1][2*j]&&dec8[2*i+1][2*j+1]))
        bias=-.00;
      for(k=0;k<15;k++){
        for(m=0;m<15;m++){
          if (var[8*i-4+k][8*j-4+m]>maxvar)maxvar=var[8*i-4+k][8*j-4+m];
          sum_var+=var[8*i-4+k][8*j-4+m];
          sum_var_1+=var_1[8*i-4+k][8*j-4+m];
        }
      }
      //sum_var=sqrt(sum_var*maxvar*15*15);
      nmr16[i][j] = (float)sum_var*(float)sum_var_1/(15*15*15*15);
      nmr8_avg = .25f*(nmr8[2*i][2*j]+nmr8[2*i][2*j+1]+nmr8[2*i+1][2*j]+nmr8[2*i+1][2*j+1]);
      //if (nmr8_avg>nmr16[i][j])nmr8_avg=nmr16[i][j];
      cg16[i][j] = .25*(cg8[2*i][2*j] + cg8[2*i][2*j+1] + cg8[2*i+1][2*j] + cg8[2*i+1][2*j+1]);
      cgs = cg16[i][j]/6 - fudge*.5*log2(nmr8_avg);
      cgl = CG16/6 - fudge*.5*log2(nmr16[i][j]);
      if (cgl+bias>=cgs)
      {
        dec8[2*i][2*j] = 2;
        dec8[2*i][2*j+1] = 2;
        dec8[2*i+1][2*j] = 2;
        dec8[2*i+1][2*j+1] = 2;
        cg16[i][j] = CG16;
      } else {
        nmr16[i][j] = nmr8_avg;
      }
      /*printf("%f ", nmr16[i][j]);*/
    }
    /*printf("\n");*/
  }

#if 1
  for(i=1;i<h32-1;i++){
    for(j=1;j<w32-1;j++){
      int k,m;
      int sum_var=0;
      int sum_var_1=0;
      float nmr16_avg;
      float cgl,cgs;
      float bias = 0;
      float maxvar = 0;

      /* Bias against 32x32 is there's a block that's better with 4x4 or 8x8 than 16x16 */
      if (!(dec8[4*i][4*j]==2&&dec8[4*i][4*j+2]==2&&dec8[4*i+2][4*j]==2&&dec8[4*i+2][4*j+2]==2))
        bias=-.0;
      for(k=0;k<31;k++){
        for(m=0;m<31;m++){
          if (var[16*i-8+k][16*j-8+m]>maxvar)maxvar=var[16*i-8+k][16*j-8+m];
          sum_var+=var[16*i-8+k][16*j-8+m];
          sum_var_1+=var_1[16*i-8+k][16*j-8+m];
        }
      }
      //sum_var=sqrt(sum_var*maxvar*31*31);
      nmr32[i][j] = (float)sum_var*(float)sum_var_1/(31*31*31*31);
      nmr16_avg = .25f*(nmr16[2*i][2*j]+nmr16[2*i][2*j+1]+nmr16[2*i+1][2*j]+nmr16[2*i+1][2*j+1]);
      //if (nmr16_avg>nmr32[i][j])nmr16_avg=nmr32[i][j];
      cg32[i][j] = .25*(cg16[2*i][2*j] + cg16[2*i][2*j+1] + cg16[2*i+1][2*j] + cg16[2*i+1][2*j+1]);
      cgs = cg32[i][j]/6 - fudge*.5*log2(nmr16_avg);
      cgl = CG32/6 - fudge*.5*log2(nmr32[i][j]);
      if (cgl+bias>=cgs)
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
      /*printf("%f ", nmr16[i][j]);*/
    }
    /*printf("\n");*/
  }
#endif
#if 1
  for(i=4;i<h8-4;i++){
    for(j=4;j<w8-4;j++){
      printf("%d ", dec8[i][j]);
    }
    printf("\n");
  }
#endif
  return 0;
}
