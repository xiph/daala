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
  return best_mode;
}

void od_intra_pred4x4_dist(od_coeff *_dist, od_coeff *_c,int _stride, int _pli){
  od_coeff c[2*4][2*4];
  float    satd;
  int      mode;
  int      i;
  int      j;
  int      k;
  int      l;
  const float satd_weights2[3][4*4] =
  {{0.230317, 0.329547, 0.457088, 0.517193, 0.315611, 0.389662, 0.525445, 0.577099,
    0.424841, 0.506366, 0.648566, 0.696878, 0.473835, 0.543352, 0.681457, 0.720876,},
   {0.458235, 0.617704, 0.811118, 0.890949, 0.601811, 0.711183, 0.906126, 0.962359,
    0.776575, 0.882826, 1.038644, 1.073067, 0.851901, 0.929541, 1.065617, 1.084039,},
   {0.481204, 0.631404, 0.829253, 0.905577, 0.615889, 0.736232, 0.930215, 0.983703,
    0.802367, 0.908049, 1.060383, 1.091684, 0.872241, 0.953467, 1.087123, 1.106003,}
  };
  for(k=0;k<8;k++){
    for(l=0;l<8;l++){
      c[k][l]=*(_c+_stride*(k-4)+l-4);
    }
  }
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    satd=0;
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        float p;
        p=0;
        for(k=0;k<8;k++){
          for(l=0;l<8;l++){
            p+=OD_INTRA_PRED_WEIGHTS_4x4[mode][i][j][k][l]*c[k][l];
          }
        }
        satd+=fabs((c[i+4][j+4]-p)*satd_weights2[_pli][i*4+j]);
      }
    }
    _dist[mode]=satd;
  }
}

void od_intra_pred4x4_get(od_coeff *_out,od_coeff *_c,int _stride, int _mode){
  od_coeff c[2*4][2*4];
  int      i;
  int      j;
  int      k;
  int      l;
  for(k=0;k<8;k++){
    for(l=0;l<8;l++){
      c[k][l]=*(_c+_stride*(k-4)+l-4);
    }
  }
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      double p;
      p=0;
      for(k=0;k<8;k++){
        for(l=0;l<8;l++){
          p+=OD_INTRA_PRED_WEIGHTS_4x4[_mode][i][j][k][l]*c[k][l];
        }
      }
      *(_out+4*i+j)=(od_coeff)floor(p+0.5);
    }
  }
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

void od_intra_pred_cdf(ogg_uint16_t cdf[OD_INTRA_NMODES],
    const unsigned char probs[OD_INTRA_NMODES][OD_INTRA_NCONTEXTS],
    const ogg_uint16_t p0[OD_INTRA_NMODES],int left,int upleft,int up)
{
  int m;
  unsigned p[OD_INTRA_NMODES];
  int sum=0;
  int curr_cdf;
  for(m=0;m<OD_INTRA_NMODES;m++)
  {
    p[m]=probs[m][(left==m)*4+(upleft==m)*2+(up==m)];
    p[m]=OD_MAXI(OD_MAXI(p[m], p0[m]>>8), 1);
    sum+=p[m];
  }

  curr_cdf=0;
  for(m=0;m<OD_INTRA_NMODES;m++)
  {
    /* Apply probability combination here: p[m] *= (sum-p[m])/(1-p[m])*/
    /* FIXME: Make this fixed-point */
    p[m] = p[m]*(sum-p[m])/(float)(256-p[m]);

    curr_cdf+=p[m];
    cdf[m]=curr_cdf;
  }
}

int od_intra_pred_search(ogg_uint16_t p0[OD_INTRA_NMODES],
    const ogg_uint16_t cdf[OD_INTRA_NMODES],
    const od_coeff dist[OD_INTRA_NMODES], ogg_uint16_t lambda, int left,
    int upleft, int up)
{
  int m;
  int best_mode;
  unsigned best_score;

  /* FIXME: Compute the log2() in fixed-point */
  best_score=dist[0]+(lambda*-log(cdf[0]/(float)cdf[OD_INTRA_NMODES-1])/log(2))/128;
  /*printf("[%f+l*%f= %f] ", dist[0]/64., -log((cdf[0])/(float)cdf[OD_INTRA_NMODES-1])/log(2), best_score/64.);*/
  best_mode=0;

  for(m=1;m<OD_INTRA_NMODES;m++){
    unsigned score;
    /* FIXME: Compute the log2() in fixed-point */
    /*printf("[%f+l*%f= ", dist[m]/64., -log((cdf[m]-cdf[m-1])/(float)cdf[OD_INTRA_NMODES-1])/log(2));*/
    score=dist[m]+(lambda*-log((cdf[m]-cdf[m-1])/(float)cdf[OD_INTRA_NMODES-1])/log(2))/128;
    /*printf("%f] ", score/64.);*/
    if (score<best_score){
      best_score=score;
      best_mode=m;
    }
    if (left!=m && up!=m && upleft!=m)
      p0[m] -= (p0[m]+256)>>9;
  }
  /*printf("\n");*/
  if (left!=best_mode && up!=best_mode && upleft!=best_mode)
  {
    /* Arbitrary ceiling at 0.75 to prevent insane p0 values */
    if (p0[best_mode]<24576)
      p0[best_mode]+=64;
  }
  return best_mode;
}
