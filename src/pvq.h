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

#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"
# include "filter.h"

extern const double *od_basis_mag[2][OD_NBSIZES];
extern const int OD_QM8[];

# define PVQ_MAX_PARTITIONS (1 + 3*(OD_NBSIZES-1))

# define OD_NOREF_ADAPT_SPEED (4)
/* FIXME: this comment is no longer accurate, we're using a smaller lambda */
/* Normalized lambda. Since we normalize the gain by q, the distortion is
   normalized by q^2 and lambda does not need the q^2 factor. At high rate,
   this would be log(2)/6, but we're using a slightly more aggressive value
   taken from:
   Li, Xiang, et al. "Laplace distribution based Lagrangian rate distortion
   optimization for hybrid video coding." Circuits and Systems for Video
   Technology, IEEE Transactions on 19.2 (2009): 193-205.
   */
# define OD_PVQ_LAMBDA (.106)

#define OD_PVQ_SKIP_ZERO 1
#define OD_PVQ_SKIP_COPY 2

#define MAXN 512

#define OD_COMPAND_SCALE (256 << OD_COEFF_SHIFT)
#define OD_COMPAND_SCALE_1 (1./OD_COMPAND_SCALE)

#define OD_QM_SIZE (20)

int od_qm_get_index(int ln, int band);

typedef struct od_qm_entry {
  int interp_q;
  int scale_q8;
  const unsigned char *qm_q4;
} od_qm_entry;

extern const od_qm_entry OD_DEFAULT_QMS[][OD_NPLANES_MAX];

extern const double *const OD_PVQ_BETA[OD_NPLANES_MAX][OD_NBSIZES];

void od_apply_qm(od_coeff *out, int out_stride, od_coeff *in, int in_stride,
 int ln, int dec, int inverse);
int od_compute_householder(double *r, int n, double gr, int *sign);
void od_apply_householder(double *x, const double *r, int n);
void od_pvq_synthesis_partial(od_coeff *xcoeff, const od_coeff *ypulse,
                                  const double *r, int n,
                                  int noref, double g,
                                  double theta, int m, int s);

double od_gain_expand(double cg, int q0, double beta);

double od_pvq_compute_gain(od_coeff *x, int n, int q0, double *g, double beta);
int od_pvq_compute_max_theta(double qcg, double beta);
double od_pvq_compute_theta(int t, int max_theta);
int od_pvq_compute_k(double qcg, int itheta, double theta, int noref, int n,
 double beta, int nodesync);

int od_vector_is_null(const od_coeff *x, int len);

#endif
