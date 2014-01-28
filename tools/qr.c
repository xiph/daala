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

#include <errno.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "qr.h"
#include "matidx.h"

/*Compute the maximum magnitude of any element in a vector.*/
static double vinfnorm(const double *v, int n) {
  double d;
  int i;
  d = 0;
  for (i = 0; i <  n; i++) {
    double avi;
    avi = fabs(v[i]);
    d = avi > d ? avi : d;
  }
  return d;
}

/*Compute the norm of a vector without destructive underflow.*/
static double v2norm(const double *v, int n) {
  double s;
  int sb;
  s = frexp(vinfnorm(v, n), &sb);
  if (s > 0 && s <= 1) {
    double vis;
    double d2;
    int i;
    d2 = 0;
    s = ldexp(1.0, -sb);
    for (i = 0; i < n; i++) {
      vis = v[i]*s;
      d2 += vis*vis;
    }
    return ldexp(sqrt(d2), sb);
  }
  else return s;
}

/*Compute the dot product of two vectors.*/
static double vdot(const double *u, const double *v, int n) {
  double d;
  int i;
  d = 0;
  for (i = 0; i < n; i++) d += u[i]*v[i];
  return d;
}

int qrdecomp_hh(double *aat, int aat_stride, double *d,
 double *qqt, int qqt_stride, double *rr, int n, int m) {
  int rank;
  int i;
  int j;
  int k;
  int l;
  rank = 0;
  l = m < n ? m : n;
  for (k = 0; k < l; k++) {
    double *aatk;
    double d2;
    aatk = aat + k*aat_stride;
    d2 = v2norm(aatk + k, m - k);
    if (d2 != 0) {
      double e;
      double s;
      if (aatk[k] < 0) d2 = -d2;
      for (i = k; i < m; i++) aatk[i] /= d2;
      e = ++aatk[k];
      for (j = k + 1; j < n; j++) {
        double *aatj;
        aatj = aat + j*aat_stride;
        s = -vdot(aatk + k, aatj + k, m - k)/e;
        for (i = k; i < m; i++) aatj[i] += s*aatk[i];
        if (rr != NULL) rr[UT_IDX(k, j, n)] = aatj[k];
      }
      rank++;
    }
    d[k] = -d2;
    if (rr != NULL) rr[UT_IDX(k, k, n)] = d[k];
  }
  /*Uncomment (along with code below for Q) to compute the _unique_
     factorization with the diagonal of R strictly non-negative.
    Unfortunately, this will not match the encoded Q and R in qrt, preventing
     the user from mixing and matching the explicit and implicit
     decompositions.*/
  /*if(rr != NULL) {
    for (k = 0; k < l; k++) {
      if (d[i] < 0) {
        for(j = k; j < n; j++) rr[UT_IDX(k, j, n)] = -rr[UT_IDX(k, j, n)];
      }
    }
  }*/
  if(qqt != NULL) {
    for (k = l; k-- > 0;) {
      double *aatk;
      double *qqtj;
      double e;
      aatk = aat + k*aat_stride;
      qqtj = qqt + k*qqt_stride;
      memset(qqtj, 0, k*sizeof(*qqtj));
      for (i = k; i < m; i++) qqtj[i] = -aatk[i];
      qqtj[k]++;
      e = aatk[k];
      if(e != 0)for(j = k + 1; j < l; j++) {
        double s;
        qqtj = qqt + j*qqt_stride;
        s = -vdot(aatk + k, qqtj + k, m - k)/e;
        for (i = k; i < m; i++) qqtj[i] += s*aatk[i];
      }
    }
    /*Uncomment (along with code above for R) to compute the _unique_
       factorization with the diagonal of R strictly non-negative.
      Unfortunately, this will not match the encoded Q and R in qrt, preventing
       the user from mixing and matching the explicit and implicit
       decompositions.*/
    /*for (k = 0; k < l; k++) if(d[k] < 0) {
      double *qqtk;
      qqtk = qqt + k*qqt_stride;
      for (i = 0; i < m; i++) qqtk[i] = -qqtk[i];
    }*/
  }
  return rank;
}

void qrsolve_hh(double *qrt, int qrt_stride, double *d,
 double *bbt, int bbt_stride, int n,int m,int l) {
  int i;
  int j;
  int k;
  /*Compose with Q^T.*/
  for (k = 0; k < n; k++) {
    double *qrtk;
    double e;
    qrtk = qrt + k*qrt_stride;
    e = qrtk[k];
    for (j = 0; j < l; j++) {
      double *bbtj;
      double s;
      bbtj = bbt + j*bbt_stride;
      s = -vdot(qrtk + k, bbtj + k, m - k)/e;
      for (i = k; i < m; i++) bbtj[i] += s*qrtk[i];
    }
  }
  /*Back-substitute through R.*/
  for (k = n; k-- > 0;) {
    double *qrtk;
    qrtk = qrt + k*qrt_stride;
    for (j = 0; j < l; j++) {
      double *bbtj;
      bbtj = bbt + j*bbt_stride;
      bbtj[k] /= d[k];
      for (i = 0; i < k; i++) bbtj[i] -= bbtj[k]*qrtk[i];
    }
  }
}
