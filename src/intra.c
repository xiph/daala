#include <stdlib.h>
#include <limits.h>
#include "filter.h"
#include "intra.h"


int od_intra_pred4x4_apply(od_coeff *_c,int _stride){
  od_coeff c[2*4][2*4];
  od_coeff phat[OD_INTRA_NMODES][4][4];
  unsigned satd;
  unsigned best_satd;
  int      mode;
  int      best_mode;
  int      i;
  int      j;
  int      k;
  int      l;
  best_satd=UINT_MAX;
  best_mode=0;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    for(k=0;k<8;k++){
      for(l=0;l<8;l++){
        c[k][l]=*(_c+_stride*(k-4)+l-4);
      }
    }
    satd=0;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        double p;
        p=0;
        for(k=0;k<8;k++){
          for(l=0;l<8;l++){
            p+=OD_INTRA_PRED_WEIGHTS_4x4[mode][i][j][k][l]*c[k][l];
          }
        }
        phat[mode][i][j]=(od_coeff)floor(p+0.5);
        satd+=abs(c[i+4][j+4]-phat[mode][i][j]);
      }
    }
    if(satd<best_satd){
      best_satd=satd;
      best_mode=mode;
    }
  }
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      *(_c+_stride*i+j)-=phat[best_mode][i][j];
    }
  }
  return mode;
}

void od_intra_pred4x4_unapply(od_coeff *_c,int _stride,int _mode){
  int i;
  int j;
  int k;
  int l;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      double p;
      p=0;
      for(k=0;k<8;k++){
        for(l=0;l<8;l++){
          p+=OD_INTRA_PRED_WEIGHTS_4x4[_mode][i][j][k][l]*
           *(_c+_stride*(k-4)+l-4);
        }
      }
      *(_c+_stride*i+j)+=floor(p+0.5);
    }
  }
}
