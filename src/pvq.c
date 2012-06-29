/*
    Daala video codec
    Copyright (C) 2012 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXN 256
#define EPSILON 1e-30

/* This is a "standard" pyramid vector quantizer search */
static void pvq_search(float *x,float *scale,float *scale_1,float g,int N,int K,int *y,int m,float lambda){
  float L1;
  float L2;
  float L1_proj;
  int   i;
  int   j;
  int   left;
  int   s[MAXN];
  float xy; /* sum(x*y) */
  float yy; /* sum(y*y) */
  float xx;

  /* Apply inverse scaling */
  if (scale!=NULL){
    for (i=0;i<N;i++){
      x[i]*=scale_1[i];
    }
  }
  /* Remove the sign and compute projection on the pyramid */
  L1=0;
  xx=0;
  for (i=0;i<N;i++){
    s[i]=x[i]>0?1:-1;
    x[i]=fabs(x[i]);
    xx += x[i]*x[i];
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
    float best_cost;
    best_num=-1;
    best_den=1e-15;
    best_cost=0;
    yy+=1;
    best_id = 0;
    for(j=0;j<N;j++){
      float tmp_xy;
      float tmp_yy;
      float cost;
      tmp_xy=xy+x[j];
      tmp_yy=yy+2*y[j];
      tmp_xy*=tmp_xy;
      cost=(j==m)?0:lambda;
      /* Trick to avoid having to divide by the denominators */
      if (tmp_xy*best_den > best_num*tmp_yy){
      /*if (tmp_xy/sqrt(xx*tmp_yy)+best_cost > best_num/sqrt(xx*best_den)+cost){*/
        best_num=tmp_xy;
        best_den=tmp_yy;
        best_cost = cost;
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

  for(i=0;i<N;i++){
    y[i]*=s[i];
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
  int K;
  OD_ASSERT(N>1);

  Q *= .42;

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

#if 1
  {
    float cg, cgr;
    L2r=0;
    for(i=0;i<N;i++){
      L2r+=r[i]*r[i];
    }
    gr=sqrt(L2r);

    cg = pow(g/Q,GAIN_EXP_1)-1.;
    if (cg<0)
      cg=0;
    cgr = pow(gr/Q,GAIN_EXP_1);

    /* Round towards zero as a slight bias */
    qg = floor(.5+cg-cgr);
    /*printf("%d ", qg);*/
    /*g = Q*pow(cg, GAIN_EXP);*/
    cg = cgr+qg;
    if (cg<0)cg=0;
    g = Q*pow(cg, GAIN_EXP);
    qg = floor(.5+cg);
  }
#else
  /*printf("%f\n", g);*/
  /* Round towards zero as a slight bias */
  qg = floor(pow(g/Q,GAIN_EXP_1)-.5);
  if (qg<0)
    qg=0;
  g = Q*pow(qg, GAIN_EXP);
#endif
  /*if(N==16)printf("%d ", qg);*/

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
  if (qg>0){
    float beta;
    float lambda;
    float theta_mod;
    int flip=0;

    if (theta>M_PI/2){
      flip=1;
      theta=M_PI-theta;
    }

#if 1
    lambda = 1./(qg*qg);
    beta = 2*sin(theta)*sin(theta)-lambda;
    if (beta<=0 /*|| theta<.5*M_PI/qg*/){
      theta=0;
    }else{
      theta = theta-lambda*cos(theta)*sin(theta)/beta;
      if (theta<0){
        theta=0;
      }
    }

#else
    lambda=2;
    beta = 1-lambda*cos(theta)/(qg*sin(theta));
    if (beta>0)
      beta = .5+.5*sqrt(beta);
    else
      beta = 0;
#endif
    theta_mod = theta/1.2389 - (theta/1.2389)*(theta/1.2389)/6.;
    qt = floor(qg*theta_mod);
    theta_mod = qt/(float)qg;
    theta = 1.2389 * 3*(1-sqrt(1-theta_mod/1.5));
    if (flip)
      theta=M_PI-theta;
    if (qt==0){
      K=0;
    }else{
      int K_large;
      K = qt*qt;
      K_large = sqrt(qt*N);
      if (K>K_large){
        K=K_large;
      }
    }
  }else{
    theta=0;
    K=0;
  }
  /*printf("%d %d\n", K, N);*/

  pvq_search(x,NULL,NULL,1,N,K,y,m,0);
  /*for(i=0;i<N;i++)printf("%d ", y[i]);*/
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
    ogg_int16_t *_scale,int *y,int N,int Q,int *qg){
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
  int K,ym;
  float cg, cgr;
  OD_ASSERT(N>1);
  Q*=.50;
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
  cg = pow(g/Q,GAIN_EXP_1)-1.;
  if (cg<0)
    cg=0;
  cgr = pow(gr/Q,GAIN_EXP_1);

  /* Round towards zero as a slight bias */
  *qg = floor(.5+cg-cgr);
  /*printf("%d ", qg);*/
  /*g = Q*pow(cg, GAIN_EXP);*/
  cg = cgr+*qg;
  if (cg<0)cg=0;
  g = Q*pow(cg, GAIN_EXP);
#if 0
  K = floor(.5+ 1.3*(M_PI/2)*(cg)/GAIN_EXP );
#else
  if (cg==0){
    K=0;
  }else{
    int K_large;
    K = cg*cg;
    K_large = sqrt(cg*N);
    if (K>K_large){
      K=K_large;
    }
  }
#endif
  /*if(N==16)printf("%d ", qg);*/

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

  pvq_search(x,NULL,NULL,1,N,K,y,m,.4/(cg*cg));
  /*printf("%d ", K-abs(y[m]));*/
  /*for(i=0;i<N;i++)printf("%d ", (m==i)?0:y[i]);*/

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

  /* Move y[m] to the front */
  ym = y[m];
  for (i=m;i>=1;i--)
    y[i] = y[i-1];
  y[0] = ym;

  /*printf("%d ", *qg);*/
  /*printf("xc1=%f\n", xc1);*/
  /*printf("y[m]=%d\n", y[m]);*/
  return m;
}


int quant_pvq_noref(ogg_int32_t *_x,float gr,
    ogg_int16_t *_scale,int *y,int N,int Q){
  float L2x;
  float g;
  float x[MAXN];
  float scale[MAXN];
  float scale_1[MAXN];
  int   i;
  int K;
  float cg, cgr;
  int qg;
  OD_ASSERT(N>1);

  Q *= 1.52;

  for(i=0;i<N;i++){
    scale[i]=_scale[i];
    scale_1[i]=1./scale[i];
  }

  L2x=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    L2x+=x[i]*x[i];
  }

  g=sqrt(L2x);


  cg = pow(g/Q,GAIN_EXP_1)-1.;
  if (cg<0)
    cg=0;
  cgr = pow(gr/Q,GAIN_EXP_1);

  /* Round towards zero as a slight bias */
  qg = floor(.5+cg-cgr);
  cg = cgr+qg;
  if (cg<0)cg=0;
  g = Q*pow(cg, GAIN_EXP);
  qg = floor(.5+cg);

  K = floor(.5+cg*cg);
  pvq_search(x,NULL,NULL,1,N,K,y,0,0);

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

  return qg;
}


int quant_scalar(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q){
  int i;
  float Q2, Q2_1;
  int K;
  OD_ASSERT(N>0);

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
  return 0;
}
