
#include "pvq.h"
#include <stdlib.h>
#include <stdio.h>

#define MAXN 256
#define EPSILON 1e-15

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


int quant_pvq(ogg_int16_t *_x,const ogg_int16_t *_r,const int *_q,
    ogg_int16_t *_scale,int *y,int N,int K,int paper_scaling){
  float L2x, L2r;
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
  float min_scale=1e15;

  for(i=0;i<N;i++){
    scale[i]=_scale[i];
    scale_1[i]=1./scale[i];
    if(scale[i]<min_scale)min_scale=scale[i];
  }

  L2x=L2r=0;
  for(i=0;i<N;i++){
    x[i]=_x[i];
    r[i]=_r[i];
    L2x+=x[i]*x[i];
  }

  if(paper_scaling){
    for(i=0;i<N;i++){
      x[i]*=scale_1[i];
      r[i]*=scale_1[i];
    }
  }
  for(i=0;i<N;i++){
    L2r+=r[i]*r[i];
  }
  g=sqrt(L2x);
  gr=sqrt(L2r);

  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m=0;
  for(i=0;i<N;i++){
    if(fabs(r[i])>maxr){
      maxr=fabs(r[i]);
      m=i;
    }
  }

  if(!paper_scaling){
    /* FIXME: It's not clear whether we should use the smallest scaling factor, the largest,
     * something in-between, something else, or just do the scaling before the Householder
     * reflection. By reducing the scaling for m, we increase the resolution in the
     * area of the prediction. Counter-intuitively, this *increases* the error, but at the
     * same time it makes y[m] larger, which is easier to code. */
    scale[m]=min_scale;
    scale_1[m]=1./scale[m];
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

  if(paper_scaling){
    pvq_search(x,NULL,NULL,g,N,K,y);
  } else {
    pvq_search(x,scale,scale_1,g,N,K,y);
  }
  /* Apply Householder reflection again to get the quantized coefficients */
  proj=0;
  for(i=0;i<N;i++){
    proj+=r[i]*x[i];
  }
  proj*=2.F/(EPSILON+L2r);
  for(i=0;i<N;i++){
    x[i]-=r[i]*proj;
  }

  if(paper_scaling){
    L2x=0;
    for(i=0;i<N;i++){
      float tmp=x[i]*scale[i];
      L2x+=tmp*tmp;
    }
    g/=EPSILON+sqrt(L2x);
    for(i=0;i<N;i++){
      x[i]*=g*scale[i];
    }
  }

  for(i=0;i<N;i++){
    _x[i]=floor(.5+x[i]);
  }

  printf("y[m]=%d\n", y[m]);
  return m;
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


