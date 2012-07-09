#include <stdio.h>
#include <float.h>
#include <math.h>
#include "cholesky.h"
#include "pythag.h"

/*Implements several Cholesky factorization routines, including an unpivoted
   routine (for well-conditioned positive definite matrices), a basic diagonal
   pivoting routine, and a strong rank-revealing routine \cite{GM04}, a
   complete orthogonal decomposition for the rank deficient case, and a minimum
   norm least squares solver.
  @ARTICLE{GM04,
    author="Ming Gu and Luiza Miranian",
    title="Strong Rank Revealing {Cholesky} Factorization",
    journal="Electronic Transactions on Numerical Analysis",
    volume=17,
    pages="76--92",
    month=Feb,
    year=2004
  }*/

/*Expand the factorization to encompass the next row of the input matrix.*/
static void ch_update(double *_rr,double _alpha,int _k,int _n){
  int i;
  int j;
  _rr[UT_IDX(_k,_k,_n)]=_alpha;
  for(i=_k+1;i<_n;i++)_rr[UT_IDX(_k,i,_n)]/=_alpha;
  for(i=_k+1;i<_n;i++){
    double t;
    t=_rr[UT_IDX(_k,i,_n)];
    for(j=i;j<_n;j++)_rr[UT_IDX(i,j,_n)]-=t*_rr[UT_IDX(_k,j,_n)];
  }
}

/*Shrink the factorization, converting the last row back to input form.*/
static void ch_downdate(double *_rr,int _k,int _n){
  double alpha;
  int    i;
  int    j;
  alpha=_rr[UT_IDX(_k,_k,_n)];
  for(i=_k+1;i<_n;i++){
    double t;
    t=_rr[UT_IDX(_k,i,_n)];
    for(j=i;j<_n;j++)_rr[UT_IDX(i,j,_n)]+=t*_rr[UT_IDX(_k,j,_n)];
  }
  for(i=_k;i<_n;i++)_rr[UT_IDX(_k,i,_n)]*=alpha;
}

static int cholesky_unpivoted(double *_rr,double _tol,int _n){
  double akk;
  int    k;
  /*We derive the tolerance from \cite{High90}.
    Higham reported that akk*50*DBL_EPSILON was always sufficient in his
     numerical experiments on matrices up to 50x50.
    We use the empirical bound of 10 on ||W||_2 observed when pivoting, even
     though we do no pivoting here, so this is optimistic.
  @INPROCEEDINGS{High90,
    author="Nicholas J. Higham",
    title="Chapter 9: Analysis of the {Cholesky} Decomposition of a
     Semi-Definite Matrix",
    editor="Maurice G. Cox and Sven J. Hammarling",
    booktitle="Reliable Numerical Computation",
    publisher="Oxford University Press",
    pages="161--185",
    year=1990
  }*/
  akk=0;
  for(k=0;k<_n;k++)if(_rr[UT_IDX(k,k,_n)]>akk)akk=_rr[UT_IDX(k,k,_n)];
  _tol*=40*_n*(_n+1)*akk;
  for(k=0;k<_n&&_rr[UT_IDX(k,k,_n)]>_tol;k++){
    ch_update(_rr,sqrt(_rr[UT_IDX(k,k,_n)]),k,_n);
  }
  return k;
}

/*Pivot: swap row and column _i of C^(_k) with row and column _k.
  Note _rr[UT_IDX(_k,_k,_n)] is not set: it is assumed this will be set by the
   caller, and the appropriate value must have already been saved.*/
static void ch_pivot(double *_rr,double *_wwt,int *_pivot,
 int _k,int _i,int _n){
  double t;
  int    j;
  for(j=0;j<_k;j++)CP_SWAP(_rr[UT_IDX(j,_k,_n)],_rr[UT_IDX(j,_i,_n)],t);
  for(j=_k+1;j<_i;j++)CP_SWAP(_rr[UT_IDX(_k,j,_n)],_rr[UT_IDX(j,_i,_n)],t);
  for(j=_i+1;j<_n;j++)CP_SWAP(_rr[UT_IDX(_k,j,_n)],_rr[UT_IDX(_i,j,_n)],t);
  _rr[UT_IDX(_i,_i,_n)]=_rr[UT_IDX(_k,_k,_n)];
  if(_wwt!=NULL){
    for(j=0;j<_k;j++)CP_SWAP(_wwt[SLT_IDX(_i,j)],_wwt[SLT_IDX(_k,j)],t);
  }
  CP_SWAP(_pivot[_k],_pivot[_i],j);
}

int cholesky(double *_rr,int *_pivot,double _tol,int _n){
  double akk;
  int    pi;
  int    j;
  int    k;
  if(_pivot==NULL)return cholesky_unpivoted(_rr,_tol,_n);
  /*Find the first pivot element.*/
  akk=0;
  pi=-1;
  for(j=0;j<_n;j++)if(_rr[UT_IDX(j,j,_n)]>akk){
    pi=j;
    akk=_rr[UT_IDX(j,j,_n)];
  }
  _tol*=40*_n*(_n+1)*akk;
  /*Initialize the pivot list.*/
  for(k=0;k<_n;k++)_pivot[k]=k;
  for(k=0;pi>=0;k++){
    double t;
    if(pi!=k)ch_pivot(_rr,NULL,_pivot,k,pi,_n);
    ch_update(_rr,sqrt(akk),k,_n);
    /*Find the next pivot element.*/
    akk=_tol;
    pi=-1;
    for(j=k+1;j<_n;j++)if(_rr[UT_IDX(j,j,_n)]>akk){
      akk=_rr[UT_IDX(j,j,_n)];
      pi=j;
    }
  }
  return k;
}

/*Forward substitution.*/
static void ch_fwd_sub(const double *_rr,double *_x,int _k,int _n){
  int i;
  for(i=0;i<_k;i++){
    int j;
    _x[i]/=_rr[UT_IDX(i,i,_n)];
    for(j=i+1;j<_k;j++)_x[j]-=_rr[UT_IDX(i,j,_n)]*_x[i];
  }
}

/*Back substitution.*/
static void ch_back_sub(const double *_rr,double *_x,int _k,int _n){
  int i;
  for(i=_k;i-->0;){
    int j;
    for(j=i+1;j<_k;j++)_x[i]-=_rr[UT_IDX(i,j,_n)]*_x[j];
    _x[i]/=_rr[UT_IDX(i,i,_n)];
  }
}

void chdecomp(double *_rr,double *_tau,int _r,int _n){
  int k;
  int i;
  int j;
  /*See Section 4 of \cite{HL69} for a derivation for a general matrix.
    We ignore the orthogonal matrix Q on the left, since we're already upper
     trapezoidal (and it would cancel with its transpose in the product R^T.R).
    @ARTICLE{HL69,
      author="Richard J. Hanson and Charles L. Lawson",
      title="Extensions and Applications of the Householder Algorithm for
       Solving Linear Least Squares",
      journal="Mathematics of Computation",
      volume=23,
      number=108,
      pages="787--812",
      month=Oct,
      year=1969
    }*/
  for(k=_r;k-->0;){
    double alpha;
    double beta;
    double s;
    double d2;
    /*Apply the Householder reflections from the previous rows.*/
    for(i=_r;--i>k;){
      s=_rr[UT_IDX(k,i,_n)];
      for(j=_r;j<_n;j++)s+=_rr[UT_IDX(k,j,_n)]*_rr[UT_IDX(i,j,_n)];
      s*=_tau[i];
      /*Note the negative here: we add an extra scale by -1 to the i'th column
         so that the diagonal entry remains positive.*/
      _rr[UT_IDX(k,i,_n)]=s-_rr[UT_IDX(k,i,_n)];
      for(j=_r;j<_n;j++)_rr[UT_IDX(k,j,_n)]-=s*_rr[UT_IDX(i,j,_n)];
    }
    /*Compute the reflection which zeros the right part of this row.*/
    alpha=_rr[UT_IDX(k,k,_n)];
    beta=alpha;
    for(j=_r;j<_n;j++){
      if(fabs(_rr[UT_IDX(k,j,_n)])>beta)beta=fabs(_rr[UT_IDX(k,j,_n)]);
    }
    s=1/beta;
    d2=(alpha*s)*(alpha*s);
    for(j=_r;j<_n;j++)d2+=(_rr[UT_IDX(k,j,_n)]*s)*(_rr[UT_IDX(k,j,_n)]*s);
    beta*=sqrt(d2);
    _tau[k]=alpha/beta+1;
    s=1/(alpha+beta);
    _rr[UT_IDX(k,k,_n)]=beta;
    for(j=_r;j<_n;j++)_rr[UT_IDX(k,j,_n)]*=s;
  }
}

int chsolve_worksz(const int *_pivot,int _n){
  return _pivot!=NULL?_n:0;
}

void chsolve(const double *_rr,const int *_pivot,const double *_tau,double *_x,
 const double *_b,double *_work,int _r,int _n){
  double *y;
  double  s;
  int     i;
  int     j;
  int     k;
  if(_pivot!=NULL){
    y=_work!=NULL?_work:(double *)malloc(chsolve_worksz(_pivot,_n)*sizeof(*y));
    for(i=0;i<_n;i++)y[i]=_b[_pivot[i]];
  }
  else{
    memmove(_x,_b,_n*sizeof(*_x));
    y=_x;
  }
  if(_r<_n){
    for(k=_r;k-->0;){
      s=y[k];
      for(j=_r;j<_n;j++)s+=y[j]*_rr[UT_IDX(k,j,_n)];
      s*=_tau[k];
      y[k]=s-y[k];
      for(j=_r;j<_n;j++)y[j]-=s*_rr[UT_IDX(k,j,_n)];
    }
  }
  ch_fwd_sub(_rr,y,_r,_n);
  ch_back_sub(_rr,y,_r,_n);
  if(_r<_n){
    memset(y+_r,0,(_n-_r)*sizeof(*y));
    for(k=0;k<_r;k++){
      s=-y[k];
      for(j=_r;j<_n;j++)s+=y[j]*_rr[UT_IDX(k,j,_n)];
      s*=_tau[k];
      y[k]=-(s+y[k]);
      for(j=_r;j<_n;j++)y[j]-=s*_rr[UT_IDX(k,j,_n)];
    }
  }
  if(_pivot!=NULL){
    for(i=0;i<_n;i++)_x[_pivot[i]]=y[i];
    if(_work==NULL)free(y);
  }
}
