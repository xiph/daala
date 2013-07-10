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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>
#include "logging.h"
#include <math.h>
#include "filter.h"

#define MAXN 256
#define EPSILON 1e-30

static const int od_layout16_offsets[4] = {0, 32, 64, 192};
extern const index_pair od_zigzag16[];
const band_layout od_layout16 = {
  od_zigzag16,
  16,
  3,
  od_layout16_offsets
};

const int od_layout8_offsets[4] = {0, 8, 16, 48};
extern const index_pair od_zigzag8[];
const band_layout od_layout8 = {
  od_zigzag8,
  8,
  3,
  od_layout8_offsets
};

static const int od_layout4_offsets[2] = {0, 15};
extern const index_pair od_zigzag4[];
const band_layout od_layout4 = {
  od_zigzag4,
  4,
  1,
  od_layout4_offsets
};

void od_bands_from_raster(const band_layout *layout, od_coeff *dst,
  od_coeff *src, int stride) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    dst[i] = src[layout->dst_table[i][1]*stride + layout->dst_table[i][0]];
  }
}

void od_raster_from_bands(const band_layout *layout, od_coeff *src,
 int stride, od_coeff *dst) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    src[layout->dst_table[i][1]*stride + layout->dst_table[i][0]] = dst[i];
  }
}

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

#if 0
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
#endif

/* This is a "standard" pyramid vector quantizer search */
static void pvq_search_rdo(int *x,float *scale,float *scale_1,float g,int N,int K,int *y,int m,float lambda){
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
  float xn[MAXN];
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
      xn[i] = x[i]*rescale;
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
    p[i]=L1_proj*xn[i];
    y[i]=floor(p[i]);
    p[i] -= y[i];
    rd[i].rd = dist_scale*(1-2*p[i]) + rate_lin*lambda*i ;
    rd[i].i = i;
    left-=y[i];
    xy+=xn[i]*y[i];
    yy+=y[i]*y[i];
  }
  /* Revert the rate cost of m and replace by the special case cost */
  rd[m].rd -= rate_ym*lambda + rate_lin*lambda*m;
#if 0
  find_nbest(rd,left-1,N);
  i=0;
  while(left>1)
  {
    int ii=rd[i].i;
    y[ii]++;
    left--;
    xy+=xn[ii];
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
      tmp_xy=xy+xn[j];
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
    xy+=xn[best_id];
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
      xn[i]=s[i]*g*scale[i]*y[i];
    }
  } else {
    /* Scale x to unit norm (without scaling) */
    g/=(EPSILON+sqrt(yy));
    for(i=0;i<N;i++){
      xn[i]=s[i]*g*y[i];
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

#define CSCALE (64)
#define CSCALE_1 (0.015625f)
#define CSHIFT (6)

/* FIXME: Replace with an actual fixed-point sqrt() */
int od_sqrt(int x) {
  return floor(.5+32768*sqrt(x));
}
/* Raises X to the power 3/4, Q15 input, Q3 output.
   We use a lookup-based domain reduction to the [0.5,1[ range. */
int od_gain_compander(int x) {
  return floor(.5 + CSCALE*pow(x*(1.f/32768), .75));
}

/* Raises X to the power 4/3, Q3 input, Q15 output.
   We use a lookup-based domain reduction to the [0.5,1[ range. */
int od_gain_expander(int x) {
  ogg_int32_t base;
  ogg_int16_t correction;
  ogg_int16_t xn;
  int ilog2x;
  /*const int C[3] = {2236, 25953, 18092};*/
  const int C[3] = {2243, 25982, 18059};
  /* Generated using:
     printf("%d, ", round(32768*(2.^([0:15]-3).^(4/3))));printf("\n"); */
  static const int expand_table[19] = {
   128, 323, 813, 2048, 5161, 13004, 32768, 82570, 208064, 524288, 1321123, 3329021,
   8388608, 21137968, 53264341, 134217728, 338207482, 852229450, 2147483647
  };
  if (x == 0) return 0;
  ilog2x = OD_ILOG(x);
  base = expand_table[ilog2x];
  /*printf("(%d)\n", base);*/
  /* Normalize x to [16384,32768[ */
  xn = x << (18 - ilog2x) >> 3;
  correction = (xn*(C[1] + (C[2]*xn >> 16)) >> 15) - C[0];
  /* FIXME: Use a 32x16 macro here. We'll have to finish with a << 1 too */
  return (long long)base*correction >> 15;
}

static int compute_k_from_gain(int cg, int N) {
  /* Compute the number of pulses K based on the quantized gain -- still work
     to do here */
  int k;
  if (cg == 0) {
    k = 0;
  }
  else {
    int K_large;
    k = floor(.5 + 0.6*CSCALE_1*CSCALE_1*cg*cg);
    K_large = floor(.5 + 1.5*CSCALE_1*cg*sqrt(N/2));
    if (k > K_large) {
      k = K_large;
    }
  }
  return k;
}

void pvq_synth(od_coeff *x, int *xn, od_coeff *r, int l2r, int cg,
  int q, int n, int syn_shift) {
  int i;
  int proj;
  int proj_1;
  int l2x;
  int g;
  int shift;
  int maxval;
  int round;
  l2x = 0;
  for (i = 0; i < n; i++) {
    l2x += xn[i]*xn[i];
  }
  /* Apply Householder reflection to the quantized coefficients */
  proj = 0;
  for (i = 0; i < n; i++) {
    proj += r[i]*xn[i];
  }
  g = od_gain_expander(q*CSCALE_1*cg)/OD_MAXI(1, od_sqrt(l2x) + 16384 >> 15);
  /* FIXME: Do this without the 64-bit math. */
  proj_1 = (ogg_int64_t)g*2*proj/OD_MAXI(1,l2r);
  maxval = OD_MAXI(g, abs(proj_1));
  shift = OD_MAXI(0, OD_ILOG(maxval)-16);
  g >>= shift;
  proj_1 >>= shift;
  shift = 15 - shift;
  round = (1 << syn_shift >> 1) << shift;
  for (i = 0; i < n; i++) {
    x[i] = (g*xn[i] - r[i]*proj_1 + round) >> (shift + syn_shift);
  }
}

#if 0 /* Disabled until it gets integerized */
int quant_pvq_theta(ogg_int32_t *x0,const ogg_int32_t *r0, ogg_int16_t *scale0,
  int *y, int n, int q0, int *qg, int shift, int intra)
/*(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int _Q, int *qg)*/{
  int l2x;
  float L2x;
  int l2r;
  float g;
  int gr;
  int x[MAXN];
  float xn[MAXN];
  int r[MAXN];
  float scale[MAXN];
  int q;
  float scale_1[MAXN];
  int   i;
  int   m;
  float s;
  float maxr=-1;
  float proj;
  float xm;
  float L2_m;
  float theta;
  int cg;              /* Companded gain of x*/
  int cgq;
  int cgr;             /* Companded gain of r*/
  int qt=0;
  int k;
  float lambda;
  OD_ASSERT(n>1);

  /* Just some calibration -- should eventually go away */
  q = od_gain_compander((1 << shift)*q0*1.3*32768);

  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */
  lambda = 0.10*q*q*CSCALE_1*CSCALE_1;

  for(i=0;i<n;i++){
    scale[i]=scale0[i];
    scale_1[i]=1./scale[i];
  }
  (void)scale_1[0];

  l2x=0;
  for(i=0;i<n;i++){
    x[i]=x0[i] << shift;
    r[i]=r0[i] << shift;
    l2x+=x[i]*x[i];
  }

  g=od_sqrt(l2x);

  l2r=0;
  for(i=0;i<n;i++){
    l2r+=r[i]*r[i];
  }
  gr=od_sqrt(l2r);

  /* compand gains */
  cg = (CSCALE*od_gain_compander(g) + q/2)/q;
  (void)intra;
  cgr = (CSCALE*od_gain_compander(gr) + q/2)/q;

  /* Doing some RDO on the gain, start by rounding down */
  *qg = floor(CSCALE_1*(cg-cgr));
  cgq = cgr+CSCALE**qg;
  if (cgq<1e-15) cgq=1e-15;
  /* Cost difference between rounding up or down */
  if ( 2*CSCALE_1*(cgq-cg)+1 + (lambda/(q*q*CSCALE_1*CSCALE_1))*(2. + (n-1)*log2(1+(float)CSCALE/(cgq)))  < 0)
  {
    (*qg)++;
    cgq = cgr+CSCALE**qg;
  }

  cg = cgr+CSCALE**qg;
  if (cg<0)cg=0;
  /* This is the actual gain the decoder will apply */
  g = od_gain_expander(q*CSCALE_1*cg);

  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<n;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "max r: %f %d %d", maxr, r[m], m));
  s=r[m]>0?1:-1;

  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m]+=(1.f/32768)*gr*s;

  l2r=0;
  for(i=0;i<n;i++){
    l2r+=r[i]*r[i];
  }

  /* Apply Householder reflection */
  proj=0;
  for(i=0;i<n;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+l2r);
  for(i=0;i<n;i++){
    x[i]-=r[i]*proj;
  }

  xm=-x[m]*s;
  x[m]=0;
  L2_m=0;
  for(i=0;i<n;i++){
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
    lambda = 1./(CSCALE_1*CSCALE_1*cg*cg);
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
    qt = floor(CSCALE_1*cg*theta_mod);
    theta_mod = qt/(float)(CSCALE_1*cg);
    theta = 1.2389 * 3*(1-sqrt(1-theta_mod/1.5));
    if (flip)
      theta=M_PI-theta;
    if (qt==0){
      k=0;
    }else{
      int K_large;
      k = qt*qt;
      K_large = sqrt(qt*n);
      if (k>K_large){
        k=K_large;
      }
    }
  }else{
    theta=0;
    k=0;
  }
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", k, n));

#if 0
  x[m] = 0;
  pvq_search (x,NULL,NULL,1,n,k,y);
#else
  {
    x[m] = 0;
    pvq_search_rdo(x,NULL,NULL,1,n,k,y,m,.0*lambda/(CSCALE_1*CSCALE_1*cg*cg));
    {
      float sum=1e-15;
      for (i=0;i<n;i++) sum += y[i]*y[i];
      sum = 1/sqrt(sum);
      for (i=0;i<n;i++) xn[i] = y[i]*sum;
    }
  }
#endif

  for(i=0;i<n;i++) {
    OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", y[i]));
  }
  if (n==32) {
    for(i=0;i<n;i++) {
      OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", y[i]));
    }
    OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "\n"));
  }

  for(i=0;i<n;i++){
    xn[i]*=sin(theta);
  }
  xn[m]=-s*cos(theta);

  /* Apply Householder reflection again to get the quantized coefficients */
  proj=0;
  for(i=0;i<n;i++){
    proj+=r[i]*xn[i];
  }
  proj*=2.F/(EPSILON+l2r);
  for(i=0;i<n;i++){
    xn[i]-=r[i]*proj;
  }

  L2x=0;
  for(i=0;i<n;i++){
    float tmp=xn[i]/* *scale[i]*/;
    L2x+=tmp*tmp;
  }
  g/=EPSILON+32768*sqrt(L2x);
  for(i=0;i<n;i++){
    xn[i]*=g/* *scale[i]*/;
  }

  for(i=0;i<n;i++){
    x0[i]=floor(.5+xn[i]*(1./(1<<shift)));
  }

  for (i = m; i >= 1; i--)
    y[i] = y[i-1];
  y[0] = qt;
  
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "y[m]=%d", y[m]));
  return m;
}
#endif

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
int quant_pvq(ogg_int32_t *x0,const ogg_int32_t *r0, ogg_int16_t *scale0,
 int *y, int n, int q0, int *qg, int shift, int intra) {
  int l2x;
  int l2r;
  int g;               /* L2-norm of x */
  int gr;              /* L2-norm of r */
  int x[MAXN];
  int r[MAXN];
  int scale[MAXN];
  int q;
  int scale_1[MAXN];
  int i;
  int m;
  int s;
  int maxr;
  int proj;
  float proj_1;
  int k;
  int ym;
  int cg;              /* Companded gain of x*/
  float cgq;
  int cgr;             /* Companded gain of r*/
  float lambda;
  int gain_offset;
  OD_ASSERT(n > 1);
  maxr = -1;
  /* Just some calibration -- should eventually go away */
  /* Converts Q to the "companded domain" */
  q = od_gain_compander((1 << shift)*q0*1.3*32768);
  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */
  lambda = 0.10*q*q*CSCALE_1*CSCALE_1;
  for (i = 0; i < n; i++) {
    scale[i] = scale0[i];
    scale_1[i] = 1./scale[i];
  }
  (void)scale_1[0];
  l2x = 0;
  for (i = 0; i < n; i++) {
    x[i] = x0[i] << shift;
    r[i] = r0[i] << shift;
    l2x += x[i]*x[i];
  }
  g = od_sqrt(l2x);
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += r[i]*r[i];
  }
  gr = od_sqrt(l2r);
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d", g));
  /* compand gain of x and subtract a constant for "pseudo-RDO" purposes */
  cg = (CSCALE*od_gain_compander(g) + q/2)/q;
  if (cg < 0) cg = 0;
  /* FIXME: Make that offset adaptive */
  gain_offset = intra ? 13*q : 0;
  cgr = (CSCALE*od_gain_compander(gr) + q/2 + gain_offset)/q;
  /* Gain quantization. Round to nearest because we've already reduced cg.
     Maybe we should have a dead zone */
  /* Doing some RDO on the gain, start by rounding down */
  *qg = (cg - cgr) >> CSHIFT;
  cgq = CSCALE_1*cgr+*qg;
  if (cgq < 1e-15) cgq = 1e-15;
  /* Cost difference between rounding up or down */
  if (2*(cgq - CSCALE_1*cg) + 1
   + (CSCALE*CSCALE*lambda/(q*q))*(2. + (n - 1)*log2(1 + 1./(cgq))) < 0) {
    (*qg)++;
    cgq = CSCALE_1*cgr + *qg;
  }
  cg = cgr + CSCALE**qg;
  if (cg < 0) cg=0;
  /* This is the actual gain the decoder will apply */
  k = compute_k_from_gain(cg, n);
  if (k == 0) {
    cg = 0;
  }
  /* TODO: actually add a reasonable log statement here.
    
    if(N==16) OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", *qg));*/

  /*for(i=0;i<N;i++){
    x[i]*=scale_1[i];
    r[i]*=scale_1[i];
  }*/
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
  m = 0;
  for (i = 0; i < n; i++) {
    if (fabs(r[i]) > maxr) {
      maxr = fabs(r[i]);
      m = i;
    }
  }
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "max r: %d %d %d", maxr, r[m], m));
  s = r[m] > 0 ? 1 : -1;
  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m] += (gr*s + 16384) >> 15;
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += r[i]*r[i];
  }
  /* Apply Householder reflection */
  proj = 0;
  for (i = 0; i < n; i++) {
    proj += r[i]*x[i];
  }
  proj_1 = proj*2.F/(EPSILON + l2r);
  for (i = 0; i < n; i++) {
    x[i] -= r[i]*proj_1;
  }
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", k, n));
  /* Normalize lambda for quantizing on the unit circle */
  /* FIXME: See if we can avoid setting lambda to zero! */
  pvq_search_rdo(x, NULL, NULL, 1, n, k, y, m, .0*lambda*CSCALE*CSCALE/(cg*cg));
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", k-abs(y[m])));
  for (i = 0; i < n; i++) {
    OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", (m == i) ? 0 : y[i]));
  }
  OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "\n"));
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d", k - y[m], n));
  pvq_synth(x0, y, r, l2r, cg, q, n, shift);
  /* Move y[m] to the front */
  ym = y[m];
  for (i = m; i >= 1; i--)
    y[i] = y[i-1];
  /* Make y[0] positive when prediction is good  */
  y[0] = -ym*s;
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", *qg));
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "y[m]=%d", y[m]));
  return m;
}

int pvq_unquant_k(const ogg_int32_t *r, int n, int qg, int q0,
 int shift, int intra) {
  int i;
  int q;
  int l2r;
  int gr;
  int cg;
  int cgr;
  int gain_offset;
  q = od_gain_compander((1 << shift)*q0*1.3*32768);
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += r[i]*r[i] << 2*shift;
  }
  gr = od_sqrt(l2r);
  /* FIXME: Make that offset adaptive */
  gain_offset = intra ? 13*q : 0;
  cgr = (CSCALE*od_gain_compander(gr) + q/2 + gain_offset)/q;
  cg = cgr + CSCALE*qg;
  if (cg < 0) cg=0;
  return compute_k_from_gain(cg, n);
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
void dequant_pvq(ogg_int32_t *x0, const ogg_int32_t *r0, ogg_int16_t *scale0,
 int n, int q0, int qg, int shift, int intra) {
  int l2r;
  int gr;              /* L2-norm of r */
  int x[MAXN];
  int r[MAXN];
  int scale[MAXN];
  int q;
  int scale_1[MAXN];
  int i;
  int m;
  int s;
  int maxr;
  int xm;
  int cg;              /* Companded gain of x*/
  int cgr;             /* Companded gain of r*/
  int gain_offset;
  OD_ASSERT(n>1);
  maxr = -1;
  /* Just some calibration -- should eventually go away */
  /* Converts Q to the "companded domain" */
  q = od_gain_compander((1 << shift)*q0*1.3*32768);
  /* High rate predicts that the constant should be log(2)/6 = 0.115, but in
     practice, it should be lower. */
  for (i = 0; i < n; i++) {
    scale[i] = scale0[i];
    scale_1[i] = 1./scale[i];
  }
  (void)scale_1[0];
  l2r = 0;
  for (i = 0; i < n; i++) {
    x[i] = x0[i];
    r[i] = r0[i] << shift;
    l2r += r[i]*r[i];
  }
  gr = od_sqrt(l2r);
  gain_offset = intra ? 13*q : 0;
  cgr = (CSCALE*od_gain_compander(gr) + q/2 + gain_offset)/q;
  cg = cgr + CSCALE*qg;
  if (cg < 0) cg = 0;
  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m = 0;
  for (i = 0; i < n; i++) {
    if (abs(r[i]) > maxr) {
      maxr = abs(r[i]);
      m = i;
    }
  }
  s = r[m] > 0 ? 1 : -1;
  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m] += (gr*s + 16384) >> 15;
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += r[i]*r[i];
  }
  /* Move x[m] back */
  xm = x[0];
  for (i = 0; i < m; i++) x[i] = x[i+1];
  /* x[0] is positive when prediction is good  */
  x[m] = -xm*s;
  pvq_synth(x0, x, r, l2r, cg, q, n, shift);
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
    ogg_int16_t *_scale,int *y,int N,int Q) {
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
  (void)scale_1[0];

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
    ogg_int16_t *_scale,int *y,int N,int Q, ogg_int32_t *_adapt){
  int i;
  float Q2, Q2_1;
  int K;
  float lambda;
  int Kn;
  OD_ASSERT(N>0);
  (void)_scale[0];

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
  OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "\n"));
#else
  if (K!=0) {
    float alpha;
    float r, r_1;
    float E0;
    float bias;
    alpha = (float)_adapt[OD_ADAPT_K_Q8]/(float)_adapt[OD_ADAPT_SUM_EX_Q8];
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
