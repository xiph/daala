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

# define PVQ_MAX_PARTITIONS (7)

# define INTER_MASKING .25

extern const int * const od_pvq_qm[3][OD_NBSIZES];
extern const double * const od_pvq_mask[3][OD_NBSIZES];

int neg_deinterleave(int x, int ref);

int pvq_theta(od_coeff *x0, od_coeff *r0, int n, int q0, od_coeff *y, int *itheta,
 int *max_theta, int *vk, double *mask_gain, double mask);

double pvq_compute_gain(od_coeff *x, int n, double q, double *g, double mask);
int pvq_compute_max_theta(double qcg, double mask);
double pvq_compute_theta(int t, int max_theta);
int pvq_compute_k(double qcg, double theta, int noref, int n, double mask);

double pvq_synthesis(od_coeff *x0, od_coeff *y, double *r, int n, double gr,
 int noref, int qg, double gain_offset, double theta, double q, double mask);

double pvq_interband_masking(double inter, double curr, double mask);

#endif
