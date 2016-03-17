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
# include "generic_code.h"

extern const double *OD_BASIS_MAG[2][OD_NBSIZES + 1];
extern const int OD_QM8_Q4_FLAT[];
extern const int OD_QM8_Q4_HVS[];

extern const uint16_t EXP_CDF_TABLE[][16];
extern const uint16_t LAPLACE_OFFSET[];

# define PVQ_MAX_PARTITIONS (9)

# define OD_NOREF_ADAPT_SPEED (4)
/* Normalized lambda for PVQ quantizer. Since we normalize the gain by q, the
   distortion is normalized by q^2 and lambda does not need the q^2 factor.
   At high rate, this would be log(2)/6, but we're using a slightly more
   aggressive value, closer to:
   Li, Xiang, et al. "Laplace distribution based Lagrangian rate distortion
   optimization for hybrid video coding." Circuits and Systems for Video
   Technology, IEEE Transactions on 19.2 (2009): 193-205.
   */
# define OD_PVQ_LAMBDA (.147)

#define OD_PVQ_SKIP_ZERO 1
#define OD_PVQ_SKIP_COPY 2

/* Maximum size for coding a PVQ band. */
#define OD_MAX_PVQ_SIZE (128)

#if defined(OD_FLOAT_PVQ)
#define OD_QM_SHIFT (15)
#else
#define OD_QM_SHIFT (11)
#endif
#define OD_QM_SCALE (1 << OD_QM_SHIFT)
#if defined(OD_FLOAT_PVQ)
#define OD_QM_SCALE_1 (1./OD_QM_SCALE)
#endif
#define OD_QM_SCALE_MAX 32767
#define OD_QM_INV_SHIFT (12)
#define OD_QM_INV_SCALE (1 << OD_QM_INV_SHIFT)
#if defined(OD_FLOAT_PVQ)
#define OD_QM_INV_SCALE_1 (1./OD_QM_INV_SCALE)
#endif
#define OD_QM_OFFSET(bs) ((((1 << 2*bs) - 1) << 2*OD_LOG_BSIZE0)/3)
#define OD_QM_STRIDE (OD_QM_OFFSET(OD_NBSIZES))
#define OD_QM_BUFFER_SIZE (2*OD_QM_STRIDE)

#if !defined(OD_FLOAT_PVQ)
#define OD_THETA_SCALE (32768*2./M_PI)
#define OD_TRIG_SCALE (32768)
#define OD_CGAIN_SHIFT (8)
#define OD_CGAIN_SCALE (1 << OD_CGAIN_SHIFT)
#else
#define OD_THETA_SCALE (1)
#define OD_TRIG_SCALE (1)
#define OD_CGAIN_SCALE (1)
#endif
#define OD_THETA_SCALE_1 (1./OD_THETA_SCALE)
#define OD_TRIG_SCALE_1 (1./OD_TRIG_SCALE)
#define OD_CGAIN_SCALE_1 (1./OD_CGAIN_SCALE)
#define OD_CGAIN_SCALE_2 (OD_CGAIN_SCALE_1*OD_CGAIN_SCALE_1)

/* Largest PVQ partition is half the coefficients of largest block size. */
#define MAXN (OD_BSIZE_MAX*OD_BSIZE_MAX/2)

#define OD_COMPAND_SHIFT (8 + OD_COEFF_SHIFT)
#define OD_COMPAND_SCALE (1 << OD_COMPAND_SHIFT)
#define OD_COMPAND_SCALE_1 (1./OD_COMPAND_SCALE)

#define OD_QM_SIZE (OD_NBSIZES*(OD_NBSIZES + 1))

#define OD_FLAT_QM 0
#define OD_HVS_QM  1

# define OD_NSB_ADAPT_CTXS (4)

# define OD_ADAPT_K_Q8        0
# define OD_ADAPT_SUM_EX_Q8   1
# define OD_ADAPT_COUNT_Q8    2
# define OD_ADAPT_COUNT_EX_Q8 3

# define OD_ADAPT_NO_VALUE (-2147483647-1)

typedef struct od_pvq_adapt_ctx  od_pvq_adapt_ctx;
typedef struct od_pvq_codeword_ctx od_pvq_codeword_ctx;

struct od_pvq_codeword_ctx {
  int                 pvq_adapt[2*OD_NBSIZES*OD_NSB_ADAPT_CTXS];
  int                 pvq_k1_increment;
  /* CDFs are size 16 despite the fact that we're using less than that. */
  uint16_t            pvq_k1_cdf[12][16];
  uint16_t            pvq_split_cdf[14*7][8];
  int                 pvq_split_increment;
};

struct od_pvq_adapt_ctx {
  od_pvq_codeword_ctx pvq_codeword_ctx;
  generic_encoder     pvq_param_model[3];
  int                 pvq_ext[OD_NBSIZES*PVQ_MAX_PARTITIONS];
  int                 pvq_exg[OD_NPLANES_MAX][OD_NBSIZES][PVQ_MAX_PARTITIONS];
  int                 pvq_gaintheta_increment;
  uint16_t        pvq_gaintheta_cdf[2*OD_NBSIZES*PVQ_MAX_PARTITIONS][16];
  int                 pvq_skip_dir_increment;
  uint16_t        pvq_skip_dir_cdf[2*(OD_NBSIZES-1)][7];
};

void od_adapt_pvq_ctx_reset(od_pvq_adapt_ctx *state, int is_keyframe);
int od_pvq_size_ctx(int n);
int od_pvq_k1_ctx(int n, int orig_size);

od_val16 od_pvq_sin(od_val32 x);
od_val16 od_pvq_cos(od_val32 x);
#if !defined(OD_FLOAT_PVQ)
int od_vector_log_mag(const od_coeff *x, int n);
#endif

int od_qm_get_index(int bs, int band);

extern const double *const OD_PVQ_BETA[2][OD_NPLANES_MAX][OD_NBSIZES + 1];

void od_init_qm(int16_t *x, int16_t *x_inv, const int *qm);
int od_compute_householder(od_val16 *r, int n, od_val32 gr, int *sign,
 int shift);
void od_apply_householder(od_val16 *out, const od_val16 *x, const od_val16 *r,
 int n);
void od_pvq_synthesis_partial(od_coeff *xcoeff, const od_coeff *ypulse,
                                  const od_val16 *r, int n,
                                  int noref, od_val32 g,
                                  od_val32 theta, int m, int s,
                                  const int16_t *qm_inv);
od_val32 od_gain_expand(od_val32 cg, int q0, double beta);
od_val32 od_pvq_compute_gain(const od_val16 *x, int n, int q0, od_val32 *g,
 double beta, int bshift);
int od_pvq_compute_max_theta(od_val32 qcg, double beta);
od_val32 od_pvq_compute_theta(int t, int max_theta);
int od_pvq_compute_k(od_val32 qcg, int itheta, od_val32 theta, int noref,
 int n, double beta, int nodesync);

int od_vector_is_null(const od_coeff *x, int len);
int od_qm_offset(int bs, int xydec);

#endif
