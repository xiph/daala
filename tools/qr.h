/*Daala video codec
Copyright (c) 2007-2014 Daala project contributors.  All rights reserved.

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

#if !defined(_qr_H)
# define _qr_H (1)
# include "matidx.h"

/*Computes the (skinny) QR-decomposition A=Q.R of an m by n matrix A using
   Householder reflections.
  The m by MIN(m, n) matrix Q is column-orthonormal (so that Q^T.Q=I), and
   the MIN(m, n) by n matrix R is upper-trapezoidal (so that R_{i, j} = 0 for
   all j < i).
  Because Householder reflections operate on columns, the decomposition is done
   on the transpose of the matrix (e.g., column-major order) for cache
   efficiency.
  The result of the computation is stored in place, and not in direct form, but
   individual components can be extracted with qrdecomp_hh_get_r() and
   qrdecomp_hh_get_q().
  aat: The n by m transpose of the matrix to decompose.
       The entries to the left of the diagonal are replaced with the entries of
        R^T, and the entries to the right of the diagonal are replaced by the
        Householder vectors used to compute Q.
  aat_stride: The row stride of aat (column stride of A).
  d: Scratch space for the diagonal of R.
     This must have storage for MIN(m, n) entries.
  qqt: Outputs the MIN(n, m) by m orthogonal matrix Q^T.
       This may be NULL to skip computation of this matrix if it is not needed.
  qqt_stride: The row stride of qqt (column stride of Q).
  rr: Outputs the upper-trapezoidal matrix R.
      This may be NULL to skip the output of this matrix in upper-triangular
       form if it is not needed, otherwise it must have storage for
       UT_SZ(MIN(m, n), n) entries and rr[UT_IDX(i, j, n)] will contain
       R_{i, j} for i <= j on exit.
  Return: The rank of the matrix aat, computed as the number of non-zero
           entries on the diagonal of R.*/
int qrdecomp_hh(double *aat,int aat_stride, double *d,
 double *qqt, int qqt_stride, double *rr, int n, int m);

/*Finds the least-squares solution of an over-determined linear system
   AX ~= B using a previously computed full-rank QR-decomposition of A.
  The solution for l different column vectors of B is found simultaneously.
  qrt: The QR-decomposition, as computed by qrdecomp_hh().
       The rank of the decomposition (as returned by qrdecomp_hh()) must be n.
  qrt_stride: The row stride of qrt.
  d: The diagonal of R, as computed by qrdecomp_hh().
  bbt: The l by m transpose of the l right hand sides B.
       These are stored as row vectors to improve cache efficiency and to make
        the common case of a single solution vector simpler to handle.
       They are replaced by the solution vectors on output.*/
void qrsolve_hh(double *qrt, int qrt_stride, double *d,
 double *bbt, int bbt_stride, int n, int m, int l);

#endif
