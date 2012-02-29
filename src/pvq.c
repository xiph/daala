
#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXN 256
#define EPSILON 1e-15

static void pvq_search(float *x,int N2,int K,int *out){
  float L1;
  float g;
  int   i;
  int   j;
  int   left;
  int   s[MAXN];
  float xy;
  float yy;

  L1=0;
  for (i=0;i<N2;i++){
    s[i]=x[i]>0?1:-1;
    x[i]=fabs(x[i]);
    L1+=x[i];
  }
  g=K/(EPSILON+L1);
  left=K;
  xy=0;
  yy=0;
  for (i=0;i<N2;i++){
    out[i]=floor(g*x[i]);
    left-=out[i];
    xy+=x[i]*out[i];
    yy+=out[i]*out[i];
  }

  for(i=0;i<left;i++){
    float best_num;
    float best_den;
    int   best_id;
    best_num=-1;
    best_den=0;
    yy+=1;
    best_id = 0;
    for(j=0;j<N2;j++){
      float tmp_xy;
      float tmp_yy;
      tmp_xy=xy+x[j];
      tmp_yy=yy+2*out[j];
      tmp_xy*=tmp_xy;
      if (tmp_xy*best_den > best_num*tmp_yy){
        best_num=tmp_xy;
        best_den=tmp_yy;
        best_id=j;
      }
    }
    xy+=x[best_id];
    yy+=2*out[best_id];
    out[best_id]++;
  }

  g=1/(EPSILON+sqrt(yy));
  for(i=0;i<N2;i++){
    x[i]=s[i]*g*out[i];
  }
}


int quant_pvq(ogg_int16_t *_x,const ogg_int16_t *_r,const int *_q,
     int *out,int N,int K){
  float L2x, L2r;
  float g;
  float gr_1;
  float x[MAXN];
  float r[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float v_norm;
  float proj;

  L2x=L2r=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    r[i]=_r[i];
    L2x+=x[i]*x[i];
    L2r+=r[i]*r[i];
  }

  g=sqrt(L2x);
  gr_1=1./(EPSILON+sqrt(L2r));

  for(i=0;i<N;i++){
    r[i]*=gr_1;
  }

  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }
  /*printf("max r: %f %f %d\n", maxr, r[m], m);*/
  s=r[m]>0?1:-1;

  r[m]+=s;

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }
  v_norm=1./(EPSILON+sqrt(L2r));
  for(i=0;i<N;i++){
    r[i]*=v_norm;
  }

  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
    /*printf("%f ", x[i]);*/
  }
  /*printf("\n");*/

  for(i=0;i<N;i++){
    x[i]-=2*r[i]*proj;
    /*printf("%f ", x[i]);*/
  }
  /*printf("\n");*/

  L2x=0;
  for(i=0;i<N;i++){
    L2x+=x[i]*x[i];
  }
  /*printf("L2 = %f\n", L2x);*/

  pvq_search(x, N, K, out);

  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  for(i=0;i<N;i++){
    x[i]=g*(x[i]-2*r[i]*proj);
  }

  for(i=0;i<N;i++){
    _x[i]=floor(.5+x[i]);
  }

  return m;
}

/*#define MAIN*/
#ifdef MAIN
int main()
{
  int i;
  ogg_int16_t x[256], r[256];
  int out[256];
  printf("original x:\n");
  for(i=0;i<64;i++){
    r[i]=rand()%255-127;
    x[i]=r[i]+(rand()%65-31);
    /*r[i]=0;*/
    printf("%d ", x[i]);
  }
  printf("\n");
  printf("original r:\n");
  for(i=0;i<64;i++){
    printf("%d ", r[i]);
  }
  printf("\n");
  quant_pvq(x, r, NULL, out, 64, 64);
  printf("quantized x:\n");
  for(i=0;i<64;i++){
    printf("%d ", x[i]);
  }
  printf("\n");
  return 0;
}
#endif


