
#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXN 256
#define EPSILON 1e-30

/* This is a "standard" pyramid vector quantizer search */
static void pvq_search(float *x,float *scale,float *scale_1,float g,int N,int K,int *y){
  float L1;
  float L2;
  float L1_proj;
  int   i;
  int   j;
  int   left;
  int   s[MAXN];
  float xy; /* sum(x*y) */
  float yy; /* sum(y*y) */

  /* Apply inverse scaling */
  if (scale!=NULL){
    for (i=0;i<N;i++){
      x[i]*=scale_1[i];
    }
  }
  /* Remove the sign and compute projection on the pyramid */
  L1=0;
  for (i=0;i<N;i++){
    s[i]=x[i]>0?1:-1;
    x[i]=fabs(x[i]);
    L1+=x[i];
  }
  L1_proj=K/(EPSILON+L1);
  left=K;
  xy=0;
  yy=0;
  /* Find the first pulses based on the projection */
  for (i=0;i<N;i++){
    y[i]=floor(L1_proj*x[i]);
    left-=y[i];
    xy+=x[i]*y[i];
    yy+=y[i]*y[i];
  }

  /* Find the remaining pulses "the long way" by maximizing xy^2/yy */
  for(i=0;i<left;i++){
    float best_num;
    float best_den;
    int   best_id;
    best_num=-1;
    best_den=0;
    yy+=1;
    best_id = 0;
    for(j=0;j<N;j++){
      float tmp_xy;
      float tmp_yy;
      tmp_xy=xy+x[j];
      tmp_yy=yy+2*y[j];
      tmp_xy*=tmp_xy;
      /* Trick to avoid having to divide by the denominators */
      if (tmp_xy*best_den > best_num*tmp_yy){
        best_num=tmp_xy;
        best_den=tmp_yy;
        best_id=j;
      }
    }
    xy+=x[best_id];
    yy+=2*y[best_id];
    y[best_id]++;
  }

  if (scale!=NULL){
    L2=0;
    for(i=0;i<N;i++){
      float tmp;
      tmp=scale[i]*y[i];
      L2+=tmp*tmp;
    }
    /* Scale x to unit norm (with scaling) */
    g/=(EPSILON+sqrt(L2));
    for(i=0;i<N;i++){
      x[i]=s[i]*g*scale[i]*y[i];
    }
  } else {
    /* Scale x to unit norm (without scaling) */
    g/=(EPSILON+sqrt(yy));
    for(i=0;i<N;i++){
      x[i]=s[i]*g*y[i];
    }
  }
}



#define GAIN_EXP (4./3.)
#define GAIN_EXP_1 (1./GAIN_EXP)

int quant_pvq_theta(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q){
  float L2x,L2r;
  float g;
  float gr;
  float x[MAXN];
  float r[MAXN];
  float scale[MAXN];
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  float xm;
  float L2_m;
  float theta;
  int qg, qt;
  int theta_res;
  float Qtheta;
  int K;

  for(i=0;i<N;i++){
    scale[i]=_scale[i];
    scale_1[i]=1./scale[i];
  }

  L2x=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    r[i]=_r[i];
    L2x+=x[i]*x[i];
  }

  g=sqrt(L2x);
  /*printf("%f\n", g);*/
  /* Round towards zero as a slight bias */
  qg = floor(pow(g/Q,GAIN_EXP_1));
  g = Q*pow(qg, GAIN_EXP);
  theta_res = floor(.5+ (M_PI/2)*qg/GAIN_EXP );
  /*if(N==16)printf("%d ", qg);*/
  /*printf("%d %d ", qg, theta_res);*/

  /*if (g>100000 && g0>100000)
    printf("%f %f\n", g, g0);*/
  /*for(i=0;i<N;i++){
    x[i]*=scale_1[i];
    r[i]*=scale_1[i];
  }*/

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }
  gr=sqrt(L2r);

  /* This is where we can skip */
  /*
  if (K<=0)
  {
    for(i=0;i<N;i++){
      _x[i]=_r[i];
    }
    return 0;
  }
*/

  /*printf("%f ", xc0);*/
  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  /*printf("max r: %f %f %d\n", maxr, r[m], m);*/
  s=r[m]>0?1:-1;

  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m]+=gr*s;

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }

  /* Apply Householder reflection */
  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+L2r);
  for(i=0;i<N;i++){
    x[i]-=r[i]*proj;
  }

  xm=-x[m]*s;
  x[m]=0;
  L2_m=0;
  for(i=0;i<N;i++){
    L2_m+=x[i]*x[i];
  }
  theta=atan2(sqrt(L2_m),xm);

  /* Quantize theta and compute K */
  if (theta_res>0){
    Qtheta=(M_PI/2.)/theta_res;
    /* Round towards zero as a slight bias */
    qt = floor(theta/Qtheta);
    theta = Qtheta*qt;
    if (qt==0){
      K=0;
    }else{
      int K_large;
      float tau=Qtheta/sin(theta);
      K=floor(.5+2./(tau*tau));
      K_large = floor(.5+sqrt(2*(N-1))/tau);
      if (K>K_large){
        K=K_large;
      }
      if(K<1)
        K=1;
    }
  }else{
    theta=0;
    K=0;
  }
  /*printf("%d %d\n", K, N);*/

  pvq_search(x,NULL,NULL,1,N,K,y);
  /*if (N==32)for(i=0;i<N;i++)printf("%d ", y[i]);printf("\n");*/

  for(i=0;i<N;i++){
    x[i]*=sin(theta);
  }
  x[m]=-s*cos(theta);

  /* Apply Householder reflection again to get the quantized coefficients */
  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+L2r);
  for(i=0;i<N;i++){
    x[i]-=r[i]*proj;
  }

  L2x=0;
  for(i=0;i<N;i++){
    float tmp=x[i]/* *scale[i]*/;
    L2x+=tmp*tmp;
  }
  g/=EPSILON+sqrt(L2x);
  for(i=0;i<N;i++){
    x[i]*=g/* *scale[i]*/;
  }

  for(i=0;i<N;i++){
    _x[i]=floor(.5+x[i]);
  }

  /*printf("xc1=%f\n", xc1);*/
  /*printf("y[m]=%d\n", y[m]);*/
  return m;
}



int quant_pvq(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q){
  float L2x,L2r;
  float g;
  float gr;
  float x[MAXN];
  float r[MAXN];
  float scale[MAXN];
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  int qg;
  int K;

  for(i=0;i<N;i++){
    scale[i]=_scale[i];
    scale_1[i]=1./scale[i];
  }

  L2x=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    r[i]=_r[i];
    L2x+=x[i]*x[i];
  }

  g=sqrt(L2x);

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }
  gr=sqrt(L2r);

  /*printf("%f\n", g);*/
  float cg, cgr;
  cg = pow(g/Q,GAIN_EXP_1);
  cgr = pow(gr/Q,GAIN_EXP_1);

  /* Round towards zero as a slight bias */
  qg = floor(.5+cg-cgr);
  cg = cgr+qg;
  if (cg<0)cg=0;
  g = Q*pow(cg, GAIN_EXP);
  K = floor(.5+ 2.3*(M_PI/2)*(cg)/GAIN_EXP );
  /*if(N==16)printf("%d ", qg);*/
  /*printf("%d %d ", qg, theta_res);*/

  /*if (g>100000 && g0>100000)
    printf("%f %f\n", g, g0);*/
  /*for(i=0;i<N;i++){
    x[i]*=scale_1[i];
    r[i]*=scale_1[i];
  }*/

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }
  gr=sqrt(L2r);

  /* This is where we can skip */
  /*
  if (K<=0)
  {
    for(i=0;i<N;i++){
      _x[i]=_r[i];
    }
    return 0;
  }
*/

  /*printf("%f ", xc0);*/
  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  /*printf("max r: %f %f %d\n", maxr, r[m], m);*/
  s=r[m]>0?1:-1;

  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m]+=gr*s;

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }

  /* Apply Householder reflection */
  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+L2r);
  for(i=0;i<N;i++){
    x[i]-=r[i]*proj;
  }

  /*printf("%d ", qt);*/
  /*printf("%d %d\n", K, N);*/

  pvq_search(x,NULL,NULL,1,N,K,y);
  /*if (N==32)for(i=0;i<N;i++)printf("%d ", y[i]);printf("\n");*/

  /*printf("%d %d\n", K-y[m], N);*/

  /* Apply Householder reflection again to get the quantized coefficients */
  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+L2r);
  for(i=0;i<N;i++){
    x[i]-=r[i]*proj;
  }

  L2x=0;
  for(i=0;i<N;i++){
    float tmp=x[i]/* *scale[i]*/;
    L2x+=tmp*tmp;
  }
  g/=EPSILON+sqrt(L2x);
  for(i=0;i<N;i++){
    x[i]*=g/* *scale[i]*/;
  }

  for(i=0;i<N;i++){
    _x[i]=floor(.5+x[i]);
  }

  /*printf("xc1=%f\n", xc1);*/
  /*printf("y[m]=%d\n", y[m]);*/
  return m;
}

int quant_scalar(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q){
  int i;
  float Q2, Q2_1;
  int K;

  K=0;
  Q2=1.5*Q;
  if (N==64)
    Q2*=.8;
  Q2_1=1./Q2;
  for (i=0;i<N;i++)
  {
    int qi;
    float tmp;
    tmp=Q2_1*(_x[i]-_r[i]);
    if (tmp<.8&&tmp>-.8)
      qi=0;
    else
      qi=floor(.5+tmp);
    K+=abs(qi);
    _x[i]=_r[i]+Q2*qi;
  }
  printf("%d %d\n", K, N);
  return 0;
}
/*#define MAIN*/
#ifdef MAIN
int main()
{
  int i;
  ogg_int16_t x[256], r[256];
  ogg_int16_t scale[256];
  int out[256];
  for(i=0;i<64;i++){
    scale[i]=32+64-i;
  }
  printf("original x:\n");
  for(i=0;i<64;i++){
    r[i]=rand()%255-127;
    x[i]=r[i]+(rand()%65-31);
    /*r[i]=0;*/
    printf("%d ",x[i]);
  }
  printf("\n");
  printf("original r:\n");
  for(i=0;i<64;i++){
    printf("%d ",r[i]);
  }
  printf("\n");
  quant_pvq(x,r,NULL,scale,out,64,64, 1);
  printf("quantized x:\n");
  for(i=0;i<64;i++){
    printf("%d ",x[i]);
  }
  printf("\n");
  return 0;
}
#endif


