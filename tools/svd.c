#include <float.h>
#include <math.h>
#include <string.h>
#include "svd.h"
#include "pythag.h"

/*Computes the singular value decomposition of a matrix into matrices
   U\Sigma V^T, where U and V are orthonormal and \Sigma is diagonal.
  _w:     An _nrows by (_nrows+_ncols) matrix of storage for input and output.
          On input, the first _nrows rows contain the input matrix.
          On output, the first _nrows rows contain the columns of U, each
           scaled by the corresponding singular value.
          The caller must divide by the appropriate value to obtain a unit
           vector, if one is desired.
          The next _ncols rows contain the rows of V (_not_ V^T).
  _s:     On output contains the squares of the _ncols singular values.
  _nrows: The number of rows of the matrix.
  _ncols: The number of columns of the matrix.
  Return: The estimated column rank of the matrix.*/
int svd(double **_w,double *_s,int _nrows,int _ncols){
  double e2;
  double tol;
  int    i;
  int    j;
  int    col_rank;
  int    slimit;
  int    sweepi;
  int    nrots;
  e2=8*_nrows*DBL_EPSILON*DBL_EPSILON;
  tol=0.125*DBL_EPSILON;
  for(i=0;i<_ncols;i++)for(j=0;j<_ncols;j++)_w[i+_nrows][j]=i==j;
  col_rank=_ncols;
  slimit=_ncols>24?_ncols>>2:6;
  sweepi=0;
  do{
    int k;
    nrots=0;
    for(j=0;j<col_rank-1;j++)for(k=j+1;k<col_rank;k++){
      double vt;
      double ee;
      double ff;
      double gg;
      double c;
      double s;
      double x;
      double y;
      ff=ee=gg=0;
      for(i=0;i<_nrows;i++){
        ee+=_w[i][j]*_w[i][j];
        ff+=_w[i][j]*_w[i][k];
        gg+=_w[i][k]*_w[i][k];
      }
      _s[j]=ee;
      _s[k]=gg;
      if(ee>=gg){
        if(ee<=e2*_s[0]||fabs(ff)<=tol*ee)continue;
        ff/=ee;
        gg=1-gg/ee;
        vt=pythag(2*ff,gg);
        c=sqrt(fabs(0.5*(1+gg/vt)));
        s=ff/(vt*c);
      }
      else{
        ff/=gg;
        ee=ee/gg-1;
        vt=pythag(2*ff,ee);
        s=sqrt(fabs(0.5*(1-ee/vt)));
        if(ff<0)s=-s;
        c=ff/(vt*s);
      }
      for(i=0;i<_nrows+_ncols;i++){
        x=_w[i][j];
        y=_w[i][k];
        _w[i][j]=x*c+y*s;
        _w[i][k]=y*c-x*s;
      }
      nrots++;
    }
    while(col_rank>1&&_s[col_rank-1]<=_s[0]*tol+tol*tol)col_rank--;
  }
  while(nrots>0&&++sweepi<=slimit);
  return col_rank-(col_rank>0&&_s[0]<=0);
}

/*Computes the Moore-Penrose pseudoinverse of a matrix using an SVD.
  _w:     An _nrows by (_nrows+_ncols) matrix of storage for input and output.
          On input, the first _nrows rows contain the input matrix.
          On output, the first _nrows rows contain the transpose of the
           pseudoinverse.
          The next _ncols rows are temporary storage.
          On output the contents are undefined.
  _s:     _ncols temporary values.
          On output the contents are undefined.
  _nrows: The number of rows of the matrix.
  _ncols: The number of columns of the matrix.
  Return: The estimated column rank of the matrix.*/
int svd_pseudoinverse(double **_w,double *_s,int _nrows,int _ncols){
  int rank;
  int i;
  int j;
  int k;
  rank=svd(_w,_s,_nrows,_ncols);
  /*There are generally far fewer columns than rows.
    Fold the singular value inverses into V.*/
  for(j=0;j<_ncols;j++)for(i=0;i<rank;i++)_w[_nrows+j][i]/=_s[i];
  /*Recompose the remaining matrices.*/
  for(i=0;i<_nrows;i++){
    memset(_s,0,_ncols*sizeof(*_s));
    for(j=0;j<_ncols;j++)for(k=0;k<rank;k++)_s[j]+=_w[i][k]*_w[_nrows+j][k];
    memcpy(_w[i],_s,_ncols*sizeof(*_s));
  }
  return rank;
}
