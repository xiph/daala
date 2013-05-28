/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>
#include "logging.h"
#include <math.h>

#define MAXN 256
#define EPSILON 1e-30

typedef struct {
  int i;
  float rd;
} RDOEntry;

#if 0
static int compare(const RDOEntry *a, const RDOEntry *b)
{
  if (a->rd<b->rd)
    return -1;
  else if (a->rd>b->rd)
    return 1;
  else
    return 0;
}
#endif

#define SWAP(x,a,b) do{RDOEntry tmp=x[b];x[b]=x[a];x[a]=tmp;}while(0);
static void find_nbest(RDOEntry *x, int n, int len)
{
  int begin, end;
  begin=0;
  end=len;
  if (n<=0)
    return;
  while(1) {
    int i;
    int index;
    int pivot;
    float pval;
    pivot=(end+begin)/2;
    pval = x[pivot].rd;
    index = begin;
    SWAP(x,pivot,end-1);
    for (i=begin;i<end-1;i++)
    {
      if (x[i].rd<pval)
      {
        SWAP(x,i,index);
        index++;
      }
    }
    SWAP(x,index,end-1);
    if (index<n-1)
    {
      begin=index+1;
    } else if (index>n)
    {
      end=index;
    } else {
      break;
    }
  };
}

/* This is a "standard" pyramid vector quantizer search */
static void pvq_search_rdo(float *x,float *scale,float *scale_1,float g,int N,int K,int *y,int m,float lambda){
  float L1;
  float L2;
  float L1_proj;
  float dist_scale;
  int   i;
  int   j;
  int   left;
  int   s[MAXN];
  float xy; /* sum(x*y) */
  float yy; /* sum(y*y) */
  float xx;
  float p[MAXN];
  RDOEntry rd[MAXN];
  /* Some simple RDO constants that should eventually depend on the state of the pvq encoders */
  const float rate_ym = .3; /* Cost advantage of y[m] compared to y[0] */
  const float rate_lin = .1; /* Cost penalty of y[n+1] compared to y[n] */

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
  {
    float rescale = 1./sqrt(1e-15+xx);
    for (i=0;i<N;i++){
      x[i] *= rescale;
    }
    L1 *= rescale;
    xx=1;
  }
  dist_scale = 1.f/(1e-15+L1*L1);
  L1_proj=K/(EPSILON+L1);
  left=K;
  xy=0;
  yy=0;
  /* Find the first pulses based on the projection */
  for (i=0;i<N;i++){
    p[i]=L1_proj*x[i];
    y[i]=floor(p[i]);
    p[i] -= y[i];
    rd[i].rd = dist_scale*(1-2*p[i]) + rate_lin*lambda*i ;
    rd[i].i = i;
    left-=y[i];
    xy+=x[i]*y[i];
    yy+=y[i]*y[i];
  }
  /* Revert the rate cost of m and replace by the special case cost */
  rd[m].rd -= rate_ym*lambda + rate_lin*lambda*m;
#if 0
  qsort(rd, N-1, sizeof(RDOEntry), compare);
#else
  find_nbest(rd,left-1,N);
#endif
#if 1
  i=0;
  while(left>1)
  {
    int ii=rd[i].i;
    y[ii]++;
    left--;
    xy+=x[ii];
    yy+=2*y[ii]-1;
    i++;
  }
#endif
  /* Find the remaining pulses "the long way" by minimizing the RDO cost function */
  for(i=0;i<left;i++){
    int   best_id;
    float best_cost; /* cost has reversed sign */
    best_cost=0;
    yy+=1;
    best_id = 0;
    for(j=0;j<N;j++){
      float tmp_xy;
      float tmp_yy;
      float cost; /* cost has reversed sign */
      tmp_xy=xy+x[j];
      tmp_yy=yy+2*y[j];
      /* RDO cost function (reversed sign) */
      cost = 2*tmp_xy/sqrt(tmp_yy)-rate_lin*lambda*j;
      /* Revert the rate cost of m and replace by the special case cost */
      if (j==m)
        cost+=rate_ym*lambda+rate_lin*lambda*m;
      if (cost > best_cost){
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
    best_num=-1;
    best_den=1e-15;
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
      /*if (tmp_xy/sqrt(xx*tmp_yy)+best_cost > best_num/sqrt(xx*best_den)+cost){*/
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

  for(i=0;i<N;i++){
    y[i]*=s[i];
  }
}



#define GAIN_EXP (4./3.)
#define GAIN_EXP_1 (1./GAIN_EXP)

int quant_pvq_theta(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int _Q, int *qg){
  float L2x,L2r;
  float g;
  float gr;
  float x[MAXN];
  float r[MAXN];
  float scale[MAXN];
  float Q;
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  float xm;
  float L2_m;
  float theta;
  float cg;              /* Companded gain of x*/
  float cgq;
  float cgr;             /* Companded gain of r*/
  int qt;
  int K;
  float lambda;
  OD_ASSERT(N>1);

  /* Just some calibration -- should eventually go away */
  Q=pow(_Q*1.3,GAIN_EXP_1); /* Converts Q to the "companded domain" */

  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */
  lambda = 0.10*Q*Q;

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

  /* compand gains */
  cg = pow(g,GAIN_EXP_1)/Q;
  cgr = pow(gr,GAIN_EXP_1)/Q;

  /* Doing some RDO on the gain, start by rounding down */
  *qg = floor(cg-cgr);
  cgq = cgr+*qg;
  if (cgq<1e-15) cgq=1e-15;
  /* Cost difference between rounding up or down */
  if ( 2*(cgq-cg)+1 + (lambda/(Q*Q))*(2. + (N-1)*log2(1+1./(cgq)))  < 0)
  {
    (*qg)++;
    cgq = cgr+*qg;
  }

  cg = cgr+*qg;
  if (cg<0)cg=0;
  /* This is the actual gain the decoder will apply */
  g = pow(Q*cg, GAIN_EXP);

  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "max r: %f %f %d", maxr, r[m], m));
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
  if (cg>0){
    float beta;
    float lambda;
    float theta_mod;
    int flip=0;

    if (theta>M_PI/2){
      flip=1;
      theta=M_PI-theta;
    }

#if 1
    lambda = 1./(cg*cg);
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
    qt = floor(cg*theta_mod);
    theta_mod = qt/(float)cg;
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
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", K, N));

  /*pvq_search(x,NULL,NULL,1,N,K,y,m,0);*/
  pvq_search_rdo(x,NULL,NULL,1,N,K,y,m,.0*lambda/(cg*cg));

  for(i=0;i<N;i++) {
    OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", y[i]));
  }
  if (N==32) {
    for(i=0;i<N;i++) {
      OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", y[i]));
    }
    OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, " "));
  }

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

  
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "y[m]=%d", y[m]));
  return m;
}

/** PVQ quantizer based on a reference
 *
 * This function computes a Householder reflection that causes the reference
 * to have a single non-zero value. The reflection is then applied to x, and
 * the result is PVQ-quantized. If the prediction is good, then most of the
 * PVQ "pulses" end up at the same position as the non-zero value in the
 * reflected reference. This is probably a good way to quantize intra blocks
 * with intra prediction.
 *
 * @param [in,out] _x     coefficients being quantized (output is the quantized values)
 * @param [in]     _r     reference
 * @param [in]     _scale quantization matrix (unused for now)
 * @param [out]    _y     quantization output (to be encoded)
 * @param [in]     N      length of vectors _x, _y, and _scale
 * @param [in]     Q      quantization resolution (lower means higher quality)
 * @param [out]    qg     quantized gain (to be encoded)
 *
 * @retval position that should have the most pulses in _y
 */
int quant_pvq(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int _Q,int *qg){
  float L2x,L2r;
  float g;               /* L2-norm of x */
  float gr;              /* L2-norm of r */
  float x[MAXN];
  float r[MAXN];
  float scale[MAXN];
  float Q;
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  int   K,ym;
  float cg;              /* Companded gain of x*/
  float cgq;
  float cgr;             /* Companded gain of r*/
  float lambda;
  OD_ASSERT(N>1);

  /* Just some calibration -- should eventually go away */
  Q=pow(_Q*1.3,GAIN_EXP_1); /* Converts Q to the "companded domain" */
  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */
  lambda = 0.10*Q*Q;

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

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%f", g));
  /* compand gain of x and subtract a constant for "pseudo-RDO" purposes */
  cg = pow(g,GAIN_EXP_1)/Q;
  if (cg<0)
    cg=0;
  /* FIXME: Make that 0.2 adaptive */
  cgr = pow(gr,GAIN_EXP_1)/Q+.2;

  /* Gain quantization. Round to nearest because we've already reduced cg.
     Maybe we should have a dead zone */
#if 0
  *qg = floor(.5+cg-cgr);
#else
  /* Doing some RDO on the gain, start by rounding down */
  *qg = floor(cg-cgr);
  cgq = cgr+*qg;
  if (cgq<1e-15) cgq=1e-15;
  /* Cost difference between rounding up or down */
  if ( 2*(cgq-cg)+1 + (lambda/(Q*Q))*(2. + (N-1)*log2(1+1./(cgq)))  < 0)
  {
    (*qg)++;
    cgq = cgr+*qg;
  }
#endif
  cg = cgr+*qg;
  if (cg<0)cg=0;
  /* This is the actual gain the decoder will apply */
  g = pow(Q*cg, GAIN_EXP);

  /* Compute the number of pulses K based on the quantized gain -- still work
     to do here */
#if 0
  K = floor(.5+ 1.*(M_PI/2)*(cg)/GAIN_EXP );
#else
  if (cg==0){
    K=0;
  }else{
    int K_large;
    K = floor(.5+0.6*cg*cg);
    K_large = floor(.5+1.5*cg*sqrt(N/2));
    if (K>K_large){
      K=K_large;
    }
  }
#endif
  if (K==0)
  {
    g=0;
    cg=0;
  }

  /* TODO: actually add a reasonable log statement here.
    
    if(N==16) OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", *qg));*/

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

  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "max r: %f %f %d", maxr, r[m], m));
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

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", K, N));

  /* Normalize lambda for quantizing on the unit circle */
  /* FIXME: See if we can avoid setting lambda to zero! */
  pvq_search_rdo(x,NULL,NULL,1,N,K,y,m,.0*lambda/(cg*cg));
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", K-abs(y[m])));
  for(i=0;i<N;i++) {
    OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", (m==i)?0:y[i]));
  }

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", K-y[m], N));

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
  /* Make y[0] positive when prediction is good  */
  y[0] = -ym*s;

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", *qg));
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "y[m]=%d", y[m]));
  return m;
}

int pvq_unquant_k(const ogg_int32_t *_r,int _n,int _qg, int _scale){
  int    i;
  int    vk;
  double Q;
  double cgr;
  Q=pow(_scale*1.3,1./(4./3.));
  vk=0;
  for(i=0;i<_n;i++)vk+=_r[i]*_r[i];
  cgr=pow(sqrt(vk),1./(4./3.))/Q+.2;
  cgr=cgr+_qg;
  if(cgr<0)cgr=0;
  if(cgr==0){
    vk=0;
  }else{
    int K_large;
    vk=floor(.5+0.6*cgr*cgr);
    K_large=floor(.5+1.5*cgr*sqrt(_n/2));
    if(vk>K_large){
      vk=K_large;
    }
  }
  return vk;
}

/** PVQ dequantizer based on a reference
 *
 * @param [in,out] _x     coefficients being dequantized (output is the dequantized values)
 * @param [in]     _r     reference
 * @param [in]     _scale quantization matrix (unused for now)
 * @param [in]     N      length of vectors _x, _y, and _scale
 * @param [in]     Q      quantization resolution (lower means higher quality)
 * @param [in]    qg     quantized gain
 *
 */
void dequant_pvq(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int N,int _Q,int qg){
  float L2x,L2r;
  float g;               /* L2-norm of x */
  float gr;              /* L2-norm of r */
  float x[MAXN];
  float r[MAXN];
  float scale[MAXN];
  float Q;
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  int   xm;
  float cg;              /* Companded gain of x*/
  float cgr;             /* Companded gain of r*/
  OD_ASSERT(N>1);

  /* Just some calibration -- should eventually go away */
  Q=pow(_Q*1.3,GAIN_EXP_1); /* Converts Q to the "companded domain" */
  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */

  for(i=0;i<N;i++){
    scale[i]=_scale[i];
    scale_1[i]=1./scale[i];
  }

  L2r=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    r[i]=_r[i];
    L2r+=r[i]*r[i];
  }
  gr=sqrt(L2r);
  cgr = pow(gr,GAIN_EXP_1)/Q+.2;
  cg = cgr+qg;
  if (cg<0)cg=0;
  g = pow(Q*cg, GAIN_EXP);

  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }
  s=r[m]>0?1:-1;

  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m]+=gr*s;

  L2r=0;
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }

  /* Move x[m] back */
  xm = x[0];
  for (i=0;i<m;i++)
    x[i] = x[i+1];
  /* x[0] is positive when prediction is good  */
  x[m] = -xm*s;

  /* Apply Householder reflection to the quantized coefficients */
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
}

/** PVQ quantizer with no reference
 *
 * This quantizer just applies companding on the gain and does the PVQ search
 *
 * @param [in,out] _x     coefficients being quantized (output is the quantized values)
 * @param [in]     gr     reference gain (i.e. predicted valud of the gain)
 * @param [in]     _scale quantization matrix (unused for now)
 * @param [out]    _y     quantization output (to be encoded)
 * @param [in]     N      length of vectors _x, _y, and _scale
 * @param [in]     Q      quantization resolution (lower means higher quality)
 *
 * @retval quantized gain (to be encoded)
 */
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

  qg = floor(.5+cg-cgr);
  cg = cgr+qg;
  if (cg<0)cg=0;
  g = Q*pow(cg, GAIN_EXP);
  qg = floor(.5+cg);

  K = floor(.5+cg*cg);
  pvq_search(x,NULL,NULL,1,N,K,y);

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
    ogg_int16_t *_scale,int *y,int N,int Q, od_adapt_ctx *_adapt){
  int i;
  float Q2, Q2_1;
  int K;
  float lambda;
  int Kn;
  OD_ASSERT(N>0);

  K=0;
  Q2=1.4*Q;
  lambda = .115;
  Q2_1=1./Q2;
  for (i=0;i<N;i++)
  {
    float tmp;
    _x[i] -= _r[i];
    tmp=Q2_1*_x[i];
    y[i]=floor(.5+tmp);
    K+=abs(y[i]);
  }


#if 0
  /* Attempt at modelling the rate more accurately -- doesn't work for now */
  mean_k_q8=_adapt->mean_k_q8;
  mean_sum_ex_q8=_adapt->mean_sum_ex_q8;
  if(mean_k_q8<1<<23)expQ8=256*mean_k_q8/(1+mean_sum_ex_q8);
  else expQ8=mean_k_q8/(1+(mean_sum_ex_q8>>8));

  Kn=0;
  for (i=N-1;i>=0;i--)
  {
    float e;
    int Ex;
    int decay;
    float rate0, rate_1;
    float cost0,cost_1;
    float ay;

    ay=abs(y[i]);
    if (y[i]==0)
      continue;
    e = Q2_1*_x[i]-y[i];

    Kn += abs(y[i]);
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if(Ex>Kn*128)
      Ex=Kn*128;

    decay = OD_MINI(255,(int)((256*Ex/(Ex+256) + 8*Ex*Ex/(256*(Kn+2)*(Kn)*(Kn)))));
    OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "(%d %d %d %d ", expQ8, Kn, Ex, decay));
    rate0 = -log2(pow(decay/256.,ay)*(1-decay/256.)/(1-pow(decay/256.,Kn+1)));
    if(Kn==1){
      rate_1=0;
    } else {
      Ex=(2*expQ8*(Kn-1)+(N-i))/(2*(N-i));
      if(Ex>(Kn-1)*256)
        Ex=(Kn-1)*256;

      decay = OD_MINI(255,(int)((256*Ex/(Ex+256) + 8*Ex*Ex/(256*(Kn+1)*(Kn-1)*(Kn-1)))));
      OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", decay));
      rate_1 = -log2(pow(decay/256.,ay-1)*(1-decay/256.)/(1-pow(decay/256.,Kn)));
    }

    cost0 = e*e + lambda*rate0;
    cost_1 = (e+1)*(e+1) + lambda*rate_1;
    OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%f %f %f %f) ",
                    e*e, (e+1)*(e+1), rate0, rate_1));
    if (cost_1<cost0)
    {
      if (y[i]>0)
        y[i]--;
      else
        y[i]++;
      Kn--;
    }
  }
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, " "));
#else
  if (K!=0) {
    float alpha;
    float r, r_1;
    float E0;
    float bias;
    alpha = (float)_adapt->mean[OD_ADAPT_K_Q8]/(float)_adapt->mean[OD_ADAPT_SUM_EX_Q8];
    if(alpha<1)
      alpha=1;
    r = 1-alpha/N;
    if (r>.1)
      r = .1;
    E0 = alpha*K/N;
    r_1 = 1./r;
    bias = -lambda/E0;
    if (bias<-.49)
      bias=-.49;
    K=0;
    for (i=0;i<N;i++)
    {
      float tmp;
      tmp=Q2_1*_x[i];
      if (tmp>0)
        tmp+=bias;
      else
        tmp-=bias;
      y[i]=floor(.5+tmp);
      K+=abs(y[i]);
      bias *= r_1;
      if (bias<-.49)
        bias=-.49;
    }
  }
  Kn=K;
#endif

  for (i=0;i<N;i++)
    _x[i]=_r[i]+Q2*y[i];
  return Kn;
}
