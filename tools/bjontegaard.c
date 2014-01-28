/*Daala video codec
Copyright (c) 2014 Daala project contributors.  All rights reserved.

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

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cholesky.h"
#include "qr.h"
#include "svd.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define PRINT_COMP (0)
#define USE_SVD    (0)
#define USE_CHOL   (0)
#define USE_QR     (1)

#if !(USE_SVD + USE_CHOL + USE_QR == 1)
# error "One of USE_SVD, USE_CHOL or USE_QR must be enabled"
#endif

#if PRINT_COMP
void printMatrix(const char *_name,const double *_m,int _r,int _c) {
  int j;
  int i;
  printf("%s=[\n",_name);
  for (j=0;j<_r;j++) {
    for (i=0;i<_c;i++) {
      printf("%s%G",i>0?",":"",(double)_m[j*_c+i]);
    }
    printf("\n");
  }
  printf("]\n");
}
#endif

static void polyfit(const double *_x,const double *_y,int _n,double *_c,
 int _deg) {
  int     m;
  double *x;
#if USE_SVD || USE_CHOL
  double *xtx;
  double *xty;
  int     k;
#endif
  int     j;
  int     i;
  m=_deg+1;
  /* Compute the Vandermonde matrix X */
  x=(double *)malloc(sizeof(*x)*_n*m);
  for (j=0;j<_n;j++) {
    x[j*m+0]=1;
    for (i=1;i<m;i++) {
      x[j*m+i]=x[j*m+i-1]*_x[j];
    }
  }
#if PRINT_COMP
  printMatrix("x",x,_n,m);
#endif
#if USE_SVD || USE_CHOL
  /* Compute X^T*X */
  xtx=(double *)malloc(sizeof(*xtx)*m*m);
  for (j=0;j<m;j++) {
    for (i=0;i<m;i++) {
      xtx[j*m+i]=0;
      for (k=0;k<_n;k++) {
        xtx[j*m+i]+=x[k*m+j]*x[k*m+i];
      }
    }
  }
#if PRINT_COMP
  printMatrix("xtx",xtx,m,m);
#endif
  /* Compute X^T*Y */
  xty=(double *)malloc(sizeof(*xty)*m);
  for (j=0;j<m;j++) {
    xty[j]=0;
    for (k=0;k<_n;k++) {
      xty[j]+=x[k*m+j]*_y[k];
    }
  }
#if PRINT_COMP
  printMatrix("xty",xty,m,1);
#endif
#endif
#if USE_SVD
  {
    double  *w;
    double **wp;
    double  *s;
    /* Compute (X^T*X)^-1 */
    w=(double *)malloc(sizeof(*w)*m*2*m);
    wp=(double **)malloc(sizeof(*wp)*2*m);
    s=(double *)malloc(sizeof(*s)*m);
    for (j=0;j<m;j++) {
      for (i=0;i<m;i++) {
        w[j*m+i]=xtx[j*m+i];
      }
    }
    for (j=0;j<2*m;j++) {
      wp[j]=&w[j*m];
    }
    svd_pseudoinverse(wp,s,m,m);
#if PRINT_COMP
    printMatrix("xtx^-1",w,m,m);
#endif
    /* Compute (X^T*X)^-1*(X^T*Y) */
    for (j=0;j<m;j++) {
      _c[j]=0;
      for (k=0;k<m;k++) {
        _c[j]+=w[j*m+k]*xty[k];
      }
    }
    free(w);
    free(wp);
    free(s);
  }
#elif USE_CHOL
  {
    double *r;
    int    *pivot;
    int     rank;
    double *tau;
    double *work;
    r=(double *)malloc(sizeof(*r)*UT_SZ(m,m));
    pivot=(int *)malloc(sizeof(*pivot)*m);
    tau=(double *)malloc(sizeof(*tau)*m);
    work=(double *)malloc(sizeof(*work)*m);
    for (j=0;j<m;j++) {
      for (i=j;i<m;i++) {
        r[UT_IDX(j,i,m)]=xtx[j*m+i];
      }
    }
    rank=cholesky(r,pivot,DBL_EPSILON,m);
    chdecomp(r,tau,rank,m);
    chsolve(r,pivot,tau,_c,xty,work,rank,m);
    free(r);
    free(pivot);
    free(tau);
    free(work);
  }
#elif USE_QR
  {
    double *qrt;
    double *d;
    double *b;
    qrt=(double *)malloc(sizeof(*qrt)*_n*m);
    d=(double *)malloc(sizeof(*d)*MIN(m,_n));
    /* Compute x^T as input */
    for (j=0;j<_n;j++) {
      for (i=0;i<m;i++) {
        qrt[i*_n+j]=x[j*m+i];
      }
    }
    qrdecomp_hh(qrt,_n,d,NULL,0,NULL,m,_n);
    b=(double *)malloc(sizeof(*b)*_n);
    for (j=0;j<_n;j++) {
      b[j]=_y[j];
    }
    /* Solve _x*b=_y using QR decomposition */
    qrsolve_hh(qrt,_n,d,b,_n,m,_n,1);
    for (j=0;j<m;j++) {
      _c[j]=b[j];
    }
    free(qrt);
    free(d);
    free(b);
  }
#endif
  free(x);
#if USE_SVD || USE_CHOL
  free(xtx);
  free(xty);
#endif
#if PRINT_COMP
  printMatrix("c",_c,m,1);
#endif
}

static void polyint(const double *_c,int _deg,double *_ci) {
  int i;
  _ci[0]=0;
  for (i=1;i<_deg+2;i++) {
    _ci[i]=_c[i-1]/i;
  }
}

static double polyval(const double *_c,int _deg,double _x) {
  double x;
  double val;
  int    i;
  x=1;
  val=0;
  for (i=0;i<_deg+1;i++) {
    val+=_c[i]*x;
    x*=_x;
  }
  return(val);
}

static double min(double *_v,int _n) {
  int    i;
  double min;
  min=_v[0];
  for (i=1;i<_n;i++) {
    min=MIN(min,_v[i]);
  }
  return(min);
}

static double max(double *_v,int _n) {
  int    i;
  double max;
  max=_v[0];
  for (i=1;i<_n;i++) {
    max=MAX(max,_v[i]);
  }
  return(max);
}

#define DELIM (",")

static int parseInts(int *_p,int _n,char *_s) {
  int   i;
  char *s;
  i=0;
  s=strtok(_s,DELIM);
  while (s!=NULL&&i<_n) {
    _p[i]=atoi(s);
    i++;
    s=strtok(NULL,DELIM);
  }
  return s==NULL&&i==_n;
}

static int parseDoubles(double *_p,int _n,char *_s) {
  int   i;
  char *s;
  i=0;
  s=strtok(_s,DELIM);
  while (s!=NULL&&i<_n) {
    _p[i]=atof(s);
    i++;
    s=strtok(NULL,DELIM);
  }
  return s==NULL&&i==_n;
}

#define RATE (0)
#define PSNR (1)

#define DEG (3)

/* Compute the Bjontegaard metric between two RD-curves */
int main(int _argc,char *_argv[]) {
  int     type;
  int     n1;
  int    *area1;
  int    *size1;
  double *rate1;
  double *psnr1;
  int     n2;
  int    *area2;
  int    *size2;
  double *rate2;
  double *psnr2;
  int     i;
  double  p1[DEG+1];
  double  p1i[DEG+2];
  double  p2[DEG+1];
  double  p2i[DEG+2];
  double  min_int;
  double  max_int;
  double  int1;
  double  int2;
  double  avg_diff;
  if (_argc!=10) {
    printf("usage: bjontegaard <type> <n1> <area1> <size1> <psnr1> <n2> <area2> <size2> <psnr2>\n");
    return EXIT_FAILURE;
  }
  type=atoi(_argv[1]);
  n1=atoi(_argv[2]);
  n2=atoi(_argv[6]);
  area1=(int *)malloc(sizeof(*area1)*n1);
  size1=(int *)malloc(sizeof(*size1)*n1);
  rate1=(double *)malloc(sizeof(*rate1)*n1);
  psnr1=(double *)malloc(sizeof(*psnr1)*n1);
  area2=(int *)malloc(sizeof(*area2)*n2);
  size2=(int *)malloc(sizeof(*size2)*n2);
  rate2=(double *)malloc(sizeof(*rate2)*n2);
  psnr2=(double *)malloc(sizeof(*psnr2)*n2);
  if (!parseInts(area1,n1,_argv[3])) {
    printf("error parsing %i ints from '%s'\n",n1,_argv[3]);
    return EXIT_FAILURE;
  }
  if (!parseInts(size1,n1,_argv[4])) {
    printf("error parsing %i ints from '%s'\n",n1,_argv[4]);
    return EXIT_FAILURE;
  }
  if (!parseDoubles(psnr1,n1,_argv[5])) {
    printf("error parsing %i doubles from '%s'\n",n1,_argv[5]);
    return EXIT_FAILURE;
  }
  if (!parseInts(area2,n2,_argv[7])) {
    printf("error parsing %i ints from '%s'\n",n1,_argv[7]);
    return EXIT_FAILURE;
  }
  if (!parseInts(size2,n2,_argv[8])) {
    printf("error parsing %i ints from '%s'\n",n1,_argv[8]);
    return EXIT_FAILURE;
  }
  if (!parseDoubles(psnr2,n2,_argv[9])) {
    printf("error parsing %i doubles from '%s'\n",n1,_argv[9]);
    return EXIT_FAILURE;
  }
  /* Compute the rate in log scale */
  for (i=0;i<n1;i++) {
    rate1[i]=log(((double)size1[i])/area1[i]);
  }
  for (i=0;i<n2;i++) {
    rate2[i]=log(((double)size2[i])/area2[i]);
  }
  if (type==RATE) {
    /* find best fit least squares polynomial for each curve */
    polyfit(psnr1,rate1,n1,p1,DEG);
    polyfit(psnr2,rate2,n2,p2,DEG);
    /* find the overlapping range */
    min_int=MAX(min(psnr1,n1),min(psnr2,n2));
    max_int=MIN(max(psnr1,n1),max(psnr2,n2));
  }
  else {
    /* find best fit least squares polynomial for each curve */
    polyfit(rate1,psnr1,n1,p1,DEG);
    polyfit(rate2,psnr2,n2,p2,DEG);
    /* find the overlapping range */
    min_int=MAX(min(rate1,n1),min(rate2,n2));
    max_int=MIN(max(rate1,n1),max(rate2,n2));
  }
  /* compute the integral */
  polyint(p1,DEG,p1i);
  polyint(p2,DEG,p2i);
  int1=polyval(p1i,DEG+1,max_int)-polyval(p1i,DEG+1,min_int);
  int2=polyval(p2i,DEG+1,max_int)-polyval(p2i,DEG+1,min_int);
  /* find average difference */
  avg_diff=(int2-int1)/(max_int-min_int);
  if (type==RATE) {
    avg_diff=(exp(avg_diff)-1)*100;
  }
  printf("%G\n",avg_diff);
  free(size1);
  free(area1);
  free(rate1);
  free(psnr1);
  free(size2);
  free(area2);
  free(rate2);
  free(psnr2);
  return EXIT_SUCCESS;
}
