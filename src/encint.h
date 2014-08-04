/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

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

#if !defined(_encint_H)
# define _encint_H (1)
# include "../include/daala/daaladec.h"
# include "../include/daala/daalaenc.h"
# include "state.h"
# include "entenc.h"
#if defined(OD_ACCOUNTING)
# include "accounting.h"
#endif

typedef struct daala_enc_ctx od_enc_ctx;
typedef struct od_mv_est_ctx od_mv_est_ctx;
typedef struct od_enc_opt_vtbl od_enc_opt_vtbl;
typedef struct od_rollback_buffer od_rollback_buffer;

/*Constants for the packet state machine specific to the encoder.*/
/*No packet currently ready to output.*/
# define OD_PACKET_EMPTY       (0)
/*A packet ready to output.*/
# define OD_PACKET_READY       (1)
/*The number of fractional bits of precision in our \lambda values.*/
# define OD_LAMBDA_SCALE       (2)
/*The number of bits of precision to add to distortion values to match
   \lambda*R.*/
# define OD_ERROR_SCALE        (OD_LAMBDA_SCALE + OD_BITRES)

struct od_enc_opt_vtbl {
  int (*mc_compute_sad_4x4_xstride_1)(const unsigned char *src,
   int systride, const unsigned char *ref, int dystride);
  int (*mc_compute_sad_8x8_xstride_1)(const unsigned char *src,
   int systride, const unsigned char *ref, int dystride);
  int (*mc_compute_sad_16x16_xstride_1)(const unsigned char *src,
   int systride, const unsigned char *ref, int dystride);
};

struct daala_enc_ctx{
  od_state state;
  od_adapt_ctx adapt;
  od_enc_opt_vtbl opt_vtbl;
  oggbyte_buffer obb;
  od_ec_enc ec;
  int packet_state;
  int quantizer[OD_NPLANES_MAX];
  od_mv_est_ctx *mvest;
#if defined(OD_ENCODER_CHECK)
  struct daala_dec_ctx *dec;
#endif
#if defined(OD_ACCOUNTING)
  /*Account for where bits are spent in encoding.*/
  od_acct acct;
#endif
};

/** Holds important encoder information so we can roll back decisions */
struct od_rollback_buffer {
  od_ec_enc ec;
  od_adapt_ctx adapt;
};

void od_encode_checkpoint(const daala_enc_ctx *enc, od_rollback_buffer *rbuf);
void od_encode_rollback(daala_enc_ctx *enc, const od_rollback_buffer *rbuf);

od_mv_est_ctx *od_mv_est_alloc(od_enc_ctx *enc);
void od_mv_est_free(od_mv_est_ctx *est);
void od_mv_est(od_mv_est_ctx *est, int ref, int lambda);

int od_mc_compute_sad_4x4_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int od_mc_compute_sad_8x8_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int od_mc_compute_sad_16x16_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int od_mc_compute_sad_c(const unsigned char *_src, int _systride,
 const unsigned char *_ref, int _dystride, int _dxstride, int _w, int _h);

void od_enc_opt_vtbl_init_c(od_enc_ctx *enc);

# if defined(OD_X86ASM)
void od_enc_opt_vtbl_init_x86(od_enc_ctx *enc);
# endif

#endif
