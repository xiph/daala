#if !defined(_cpack_linalg_cholesky_H)
# define _cpack_linalg_cholesky_H (1)
# include "matidx.h"

/*Computes the Cholesky decomposition of the symmetric positive semidefinite
   matrix A=R^T.R, R upper-triangular.
  @param _rr    On input, contains the upper-triangular portion of A packed
                 linearly in row-wise fashion, i.e.,
                 _rr[UT_IDX(i,j,_n)]=A_ij, i<=j.
                On output, returns the upper-trapezoidal matrix R packed
                 linearly in row-wise fashion, i.e.,
                 _rr[UT_IDX(i,j,_n)]=R_ij, i<=j.
                Only the first r rows of R are computed, where r is the rank
                 of A; the remaining rows are taken to be zero, but their
                 contents in _rr are undefined.
  @param _pivot Returns the pivot list, so that column i of R corresponds to
                 row and column _pivot[i] of A.
                That is, if A'=R^T.R, then A_{_pivot[i],_pivot[j]}=A'_ij.
                This may be NULL to skip pivoting, but will result in less
                 accuracy and a less reliable positive-definite test.
  @param _tol   The expected relative error in the elements of A (e.g.,
                 DBL_EPSILON).
  @param _n     The number of rows and columns in A.
  @return The number of positive eigenvalues of A encountered, r.
          If A is positive semidefinite and pivoting was enabled, this is also
           the estimated rank of A.*/
int cholesky(double *_rr,int *_pivot,double _tol,int _n);

/*Computes the complete orthogonal decomposition of the upper-trapezoidal _r by
   _n matrix R, in the case that A was positive semidefinite.
  That is, computes R' and V such that R={{R',0}}.V, R' is an _r by _r
   upper-triangular matrix, and V is orthogonal.
  This allows computation of the pseudoinverse and the minimum norm solution of
   rank-deficient least squares problems (e.g., by chsolve()).
  @param _rr   On input, contains the upper-trapezoidal matrix R packed linearly
                in row-wise fashion, i.e., _rr[UT_IDX(i,j,_n)]=R_ij, i<=j.
               On output, returns the new upper-triangular matrix packed into
                the left portion of _rr, and the Householder vectors used to
                annihilate the last _n-_r columns in the right portion.
  @param _tau  Returns the _r scale factors of the Householder reflections.
  @param _r    The number of rows (the rank) of R, _r<=_n.
  @param _n    The number of columns of R.*/
void chdecomp(double *_rr,double *_tau,int _r,int _n);

/*Returns the number of elements needed in the work array for chsolve().
  The contents of _pivot need not be initialized, but it must be non-NULL if
   pivoting will be used.*/
int chsolve_worksz(const int *_pivot,int _n);

/*Solves the system R^T.R.x=b given a full-rank upper-triangular matrix R or a
   complete orthogonal decomposition (for rank-deficient R).
  @param _rr    In the case that _r==_n, the full-rank upper-triangular matrix
                 R as computed by cholesky() or rrcholesky().
                Otherwise, the complete orthogonal decomposition as computed by
                 chdecomp().
  @param _pivot The pivot list as computed by cholesky() or rrcholesky(), so
                 that column i of R corresponds to element _pivot[i] of _x and
                 _b.
                This may be NULL if no pivoting was done.
  @param _tau   The scale values of the Householder reflectors for a
                 rank-deficient solution.
                This paramter is ignored if _r==_n.
  @param _x     Returns the solution vector.
                This may be aliased to _b to compute the solution in-place.
  @param _b     Contains the right hand side.
  @param _work  Work array with chsolve_worksz() elements, or NULL to allocate
                 one internally.
  @param _r     The number of rows in R.
  @param _n     The number of columns in R.*/
void chsolve(const double *_rr,const int *_pivot,const double *_tau,double *_x,
 const double *_b,double *_work,int _r,int _n);

#endif
