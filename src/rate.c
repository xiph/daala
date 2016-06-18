/*Daala video codec
Copyright (c) 2002-2016 Daala project contributors.  All rights reserved.

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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "encint.h"

#include "quantizer.h"

#define OD_Q57(v) ((int64_t)(v) << 57)
#define OD_F_Q45(v) ((int64_t)(((v)*((int64_t)1 << 45))))
#define OD_F_Q12(v) ((int32_t)(((v)*((int32_t)1 << 12))))

static const int64_t OD_ATANH_LOG2[32]={
  0x32B803473F7AD0F4LL,0x2F2A71BD4E25E916LL,0x2E68B244BB93BA06LL,
  0x2E39FB9198CE62E4LL,0x2E2E683F68565C8FLL,0x2E2B850BE2077FC1LL,
  0x2E2ACC58FE7B78DBLL,0x2E2A9E2DE52FD5F2LL,0x2E2A92A338D53EECLL,
  0x2E2A8FC08F5E19B6LL,0x2E2A8F07E51A485ELL,0x2E2A8ED9BA8AF388LL,
  0x2E2A8ECE2FE7384ALL,0x2E2A8ECB4D3E4B1ALL,0x2E2A8ECA94940FE8LL,
  0x2E2A8ECA6669811DLL,0x2E2A8ECA5ADEDD6ALL,0x2E2A8ECA57FC347ELL,
  0x2E2A8ECA57438A43LL,0x2E2A8ECA57155FB4LL,0x2E2A8ECA5709D510LL,
  0x2E2A8ECA5706F267LL,0x2E2A8ECA570639BDLL,0x2E2A8ECA57060B92LL,
  0x2E2A8ECA57060008LL,0x2E2A8ECA5705FD25LL,0x2E2A8ECA5705FC6CLL,
  0x2E2A8ECA5705FC3ELL,0x2E2A8ECA5705FC33LL,0x2E2A8ECA5705FC30LL,
  0x2E2A8ECA5705FC2FLL,0x2E2A8ECA5705FC2FLL
};

static int od_ilog64(int64_t v){
  static const unsigned char OD_DEBRUIJN_IDX64[64]={
   0, 1, 2, 7, 3,13, 8,19, 4,25,14,28, 9,34,20,40,
   5,17,26,38,15,46,29,48,10,31,35,54,21,50,41,57,
   63, 6,12,18,24,27,33,39,16,37,45,47,30,53,49,56,
   62,11,23,32,36,44,52,55,61,22,43,51,60,42,59,58
  };
  int ret;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  ret = (int)v & 1;
  v = (v >> 1) + 1;
  ret += OD_DEBRUIJN_IDX64[v * 0x218A392CD3D5DBF >> 58 & 0x3F];
  return ret;
}

/*Computes the binary exponential of logq57.
  input: a log base 2 in Q57 format
  output: a 64 bit integer in Q0 (no fraction) */
static int64_t od_bexp64(int64_t logq57){
  int64_t w;
  int64_t z;
  int ipart;
  ipart = (int)(logq57 >> 57);
  if (ipart < 0) return 0;
  if (ipart >= 63) return 0x7FFFFFFFFFFFFFFFLL;
  z = logq57 - OD_Q57(ipart);
  if (z) {
    int64_t mask;
    long wlo;
    int i;
    /*C doesn't give us 64x64->128 muls, so we use CORDIC.
      This is not particularly fast, but it's not being used in time-critical
       code; it is very accurate.*/
    /*z is the fractional part of the log in Q62 format.
      We need 1 bit of headroom since the magnitude can get larger than 1
       during the iteration, and a sign bit.*/
    z <<= 5;
    /*w is the exponential in Q61 format (since it also needs headroom and can
       get as large as 2.0); we could get another bit if we dropped the sign,
       but we'll recover that bit later anyway.
      Ideally this should start out as
        \lim_{n->\infty} 2^{61}/\product_{i=1}^n \sqrt{1-2^{-2i}}
       but in order to guarantee convergence we have to repeat iterations 4,
        13 (=3*4+1), and 40 (=3*13+1, etc.), so it winds up somewhat larger.*/
    w = 0x26A3D0E401DD846DLL;
    for (i = 0;; i++) {
      mask = -(z < 0);
      w += ((w >> (i + 1)) + mask) ^ mask;
      z -= (OD_ATANH_LOG2[i] + mask) ^ mask;
      /*Repeat iteration 4.*/
      if (i >= 3) break;
      z <<= 1;
    }
    for (;; i++) {
      mask = -(z < 0);
      w += ((w >> (i + 1)) + mask) ^ mask;
      z -= (OD_ATANH_LOG2[i] + mask) ^ mask;
      /*Repeat iteration 13.*/
      if (i >= 12) break;
      z <<= 1;
    }
    for(; i < 32; i++){
      mask = -(z < 0);
      w += ((w >> (i + 1)) + mask) ^ mask;
      z = (z - ((OD_ATANH_LOG2[i] + mask) ^ mask)) << 1;
    }
    wlo = 0;
    /*Skip the remaining iterations unless we really require that much
       precision.
      We could have bailed out earlier for smaller iparts, but that would
       require initializing w from a table, as the limit doesn't converge to
       61-bit precision until n=30.*/
    if (ipart > 30) {
      /*For these iterations, we just update the low bits, as the high bits
         can't possibly be affected.
        OD_ATANH_LOG2 has also converged (it actually did so one iteration
         earlier, but that's no reason for an extra special case).*/
      for (;; i++) {
        mask = -(z < 0);
        wlo += ((w >> i) + mask) ^ mask;
        z -= (OD_ATANH_LOG2[31] + mask) ^ mask;
        /*Repeat iteration 40.*/
        if (i >= 39) break;
        z <<= 1;
      }
      for (; i<61; i++) {
        mask = -(z < 0);
        wlo += ((w >> i) + mask) ^ mask;
        z = (z - ((OD_ATANH_LOG2[31] + mask) ^ mask)) << 1;
      }
    }
    w = (w << 1) + wlo;
  }
  else {
    w = (int64_t)1 << 62;
  }
  if (ipart < 62) {
    w = ((w >> (61 - ipart)) + 1) >> 1;
  }
  return w;
}


/*Computes the binary log of w
  input: a 64-bit integer in Q0 (no fraction)
  output: a 64-bit log in Q57 */
static int64_t od_blog64(int64_t w){
  int64_t z;
  int ipart;
  if (w <= 0) return -1;
  ipart = od_ilog64(w) - 1;
  if (ipart > 61) {
    w >>= ipart - 61;
  }
  else {
    w <<= 61 - ipart;
  }
  z = 0;
  if (w & (w - 1)) {
    int64_t x;
    int64_t y;
    int64_t u;
    int64_t mask;
    int i;
    /*C doesn't give us 64x64->128 muls, so we use CORDIC.
      This is not particularly fast, but it's not being used in time-critical
       code; it is very accurate.*/
    /*z is the fractional part of the log in Q61 format.*/
    /*x and y are the cosh() and sinh(), respectively, in Q61 format.
      We are computing z = 2*atanh(y/x) = 2*atanh((w - 1)/(w + 1)).*/
    x = w + ((int64_t)1 << 61);
    y = w - ((int64_t)1 << 61);
    for (i = 0; i < 4; i++) {
      mask = -(y < 0);
      z += ((OD_ATANH_LOG2[i] >> i) + mask) ^ mask;
      u = x >> (i + 1);
      x -= ((y >> (i + 1)) + mask) ^ mask;
      y -= (u + mask) ^ mask;
    }
    /*Repeat iteration 4.*/
    for (i--; i < 13; i++) {
      mask = -(y < 0);
      z += ((OD_ATANH_LOG2[i] >> i) + mask) ^ mask;
      u = x >> (i + 1);
      x -= ((y >> (i + 1)) + mask) ^ mask;
      y -= (u + mask) ^ mask;
    }
    /*Repeat iteration 13.*/
    for (i--; i < 32; i++) {
      mask = -(y < 0);
      z += ((OD_ATANH_LOG2[i] >> i) + mask) ^ mask;
      u = x >> (i + 1);
      x -= ((y >> (i + 1)) + mask) ^ mask;
      y -= (u + mask) ^ mask;
    }
    /*OD_ATANH_LOG2 has converged.*/
    for (; i < 40; i++) {
      mask = -(y < 0);
      z += ((OD_ATANH_LOG2[31] >> i) + mask) ^ mask;
      u = x >> (i + 1);
      x -= ((y >> (i + 1)) + mask) ^ mask;
      y -= (u + mask) ^ mask;
    }
    /*Repeat iteration 40.*/
    for (i--; i < 62; i++) {
      mask = -(y < 0);
      z += ((OD_ATANH_LOG2[31] >> i) + mask) ^ mask;
      u = x >> (i + 1);
      x -= ((y >> (i + 1)) + mask) ^ mask;
      y -= (u + mask) ^ mask;
    }
    z = (z + 8) >> 4;
  }
  return OD_Q57(ipart) + z;
}

static void od_enc_rc_reset(od_enc_ctx *enc) {
  /* Stub to establish API */
  OD_UNUSED(enc);
}

int od_enc_rc_resize(od_enc_ctx *enc) {
  /* Stub to establish API */
  OD_UNUSED(enc);
  return OD_EIMPL;
}

int od_enc_rc_init(od_enc_ctx *enc, long bitrate) {
  od_rc_state *rc;
  OD_UNUSED(enc);
  rc = &enc->rc;
  if(rc->target_bitrate > 0){
    /*State has already been initialized; rather than reinitialize,
      adjust the buffering for the new target rate. */
    rc->target_bitrate = bitrate;
    return od_enc_rc_resize(enc);
  }
  rc->target_bitrate = bitrate;
  if(bitrate>0){
    /*The buffer size is set equal to the keyframe interval, clamped to the
       range [12,256] frames.
      The 12 frame minimum gives us some chance to distribute bit estimation
       errors.
      The 256 frame maximum means we'll require 8-10 seconds of pre-buffering
       at 24-30 fps, which is not unreasonable.*/
    rc->reservoir_frame_delay = enc->state.info.keyframe_rate > 256 ? 256 :
     enc->state.info.keyframe_rate;
    /*By default, enforce all buffer constraints.*/
    rc->drop_frames=1;
    rc->cap_overflow=1;
    rc->cap_underflow=0;
    rc->twopass_state=0;
    od_enc_rc_reset(enc);
  }
  /* Stub to establish API */
  return OD_EIMPL;
}

void od_enc_rc_clear(od_enc_ctx *enc) {
  /*No-op until we get to two-pass support.*/
  OD_UNUSED(enc);
}

static int quality_to_quantizer(int quality) {
  /*A quality of -1 indicates the value is unset.
    A quality of 0 indicates a request for lossless.
    Above zero, there's a simple mapping from quality to baseline quantizer to
    produce a relatively smooth rd curve.*/
  return quality <= 0 ? quality :
   (quality << OD_COEFF_SHIFT >> OD_QUALITY_SHIFT) +
   (1 << OD_COEFF_SHIFT >> 1);
}

void od_enc_rc_select_quantizers_and_lambdas(od_enc_ctx *enc,
 int is_golden_frame, int frame_type) {
  /*Can't use the OD_LOSSLESS macro, as it uses state.quantizer to intuit,
     and we've not set it yet.*/
  if (enc->quality == 0) {
    /*Lossless coding requested.*/
    enc->rc.base_quantizer = 0;
    enc->target_quantizer = 0;
    enc->state.coded_quantizer = 0;
    enc->state.quantizer = 0;
  }
  else {
    int quantizer;
    int64_t log_quantizer;
    int frame_subtype;
    int lossy_quantizer_min;
    int lossy_quantizer_max;
    int32_t mqp_Q12[OD_FRAME_NSUBTYPES];
    int64_t dqp_Q45[OD_FRAME_NSUBTYPES];
    /*Quantizer selection sticks to the codable, lossy portion of the quantizer
       range.*/
    lossy_quantizer_min =
     od_codedquantizer_to_quantizer(1);
    lossy_quantizer_max =
     od_codedquantizer_to_quantizer(OD_N_CODED_QUANTIZERS - 1);
    /*P-frames can be golden, and thus boosted.
      Boosted and un-boosted P-frames are treated as different subtypes for
       convenience. */
    frame_subtype = is_golden_frame && frame_type == OD_P_FRAME ?
     OD_GOLDEN_P_FRAME : frame_type;
    /*Stash quantizer modulation by frame type.*/
    mqp_Q12[OD_I_FRAME] = OD_F_Q12(OD_MQP_I);
    mqp_Q12[OD_P_FRAME] = OD_F_Q12(OD_MQP_P);
    mqp_Q12[OD_B_FRAME] = OD_F_Q12(OD_MQP_B);
    mqp_Q12[OD_GOLDEN_P_FRAME] = OD_F_Q12(OD_MQP_I);
    dqp_Q45[OD_I_FRAME] = OD_F_Q45(OD_DQP_I);
    dqp_Q45[OD_P_FRAME] = OD_F_Q45(OD_DQP_P);
    dqp_Q45[OD_B_FRAME] = OD_F_Q45(OD_DQP_B);
    dqp_Q45[OD_GOLDEN_P_FRAME] = OD_F_Q45(OD_DQP_I);
    if (enc->quality == -1) {
      /*A quality of -1 means quality was unset; use a default.*/
      enc->rc.base_quantizer = quality_to_quantizer(10);
    }
    else {
      enc->rc.base_quantizer = quality_to_quantizer(enc->quality);
    }
    /*As originally written, qp modulation is applied to the coded quantizer.
      Because we now have and use a more precise target quantizer for various
       calculation, that needs to be modulated as well.
      Calculate what is, effectively, a fractional coded quantizer. */
    /*Get the log2 quantizer in Q57 (normalized for coefficient shift).*/
    log_quantizer = od_blog64(enc->rc.base_quantizer) - OD_Q57(OD_COEFF_SHIFT);
    /*log_quantizer to Q21.*/
    log_quantizer >>= 36;
    /*scale log quantizer, result is Q33.*/
    log_quantizer *= OD_LOG_QUANTIZER_BASE_Q12;
    /*Add Q33 offset to Q33 log_quantizer.*/
    log_quantizer += OD_LOG_QUANTIZER_OFFSET_Q45 >> 12;
    /*Modulate quantizer according to frame type; result is Q45.*/
    log_quantizer *= mqp_Q12[frame_subtype];
    /*Add Q45 boost/cut to Q45 fractional coded quantizer.*/
    log_quantizer += dqp_Q45[frame_subtype];
    /*Back to log2 quantizer in Q57.*/
    log_quantizer = (log_quantizer - OD_LOG_QUANTIZER_OFFSET_Q45)*
     OD_LOG_QUANTIZER_EXP_Q12 + OD_Q57(OD_COEFF_SHIFT);
    /*Convert Q57 log2 quantizer to unclamped linear target quantizer value.*/
    quantizer = od_bexp64(log_quantizer);
    /*Clamp and save target quantizer.*/
    enc->target_quantizer =
     OD_CLAMPI(lossy_quantizer_min, quantizer, lossy_quantizer_max);
    /*Coded quantizer is modulated independently to preserve the
       desired, coarser integer rouding behavior.
      Specifically, we want to make sure that an integer coded quantizer boost
       produces exactly that coded quantizer boost, and doesn't narrowly miss
       high or low.
      That narrow miss is a relatively much larger/smaller change than would be
       made to the target quantizer.*/
    log_quantizer =
      (int64_t)od_quantizer_to_codedquantizer(enc->rc.base_quantizer) << 33;
    log_quantizer *= mqp_Q12[frame_subtype];
    log_quantizer += dqp_Q45[frame_subtype];
    enc->state.coded_quantizer =
     OD_CLAMPI(1, log_quantizer >> 45, OD_N_CODED_QUANTIZERS - 1);
    enc->state.quantizer =
     od_codedquantizer_to_quantizer(enc->state.coded_quantizer);
  }
  /*Generate encoding lambdas from target and actual quantizers.*/
  /*Motion estimation normalized lambda is 2320000 ~= 0.55313
     (or sqrt(0.30595)) in Q22.
    The lower bound of 40 is there because we do not yet consider PVQ noref
     flags during the motion search, so we waste far too many bits trying to
     predict unpredictable areas when lambda is too small.
    Hopefully when we fix that, we can remove the limit.*/
  enc->mv_rdo_lambda =
   OD_MAXI(((2320000 + (((1 << OD_COEFF_SHIFT) - 1) >> 1)) >> OD_COEFF_SHIFT)*
   enc->target_quantizer >> (22 - OD_LAMBDA_SCALE), 40);
  /*We need a normalized PVQ lambda based on the target (not actual) quantizer
     for use within PVQ after we've already performed quantization.*/
  enc->pvq_norm_lambda = OD_PVQ_LAMBDA;
  /*The PVQ RDO lambda is used for RDO calculations involving unquantized
     data.*/
  enc->pvq_rdo_lambda = OD_PVQ_LAMBDA*
   enc->target_quantizer*enc->target_quantizer;
  /*The blocksize RDO lambda is calculated from the PVQ RDO lambda.*/
  enc->bs_rdo_lambda = OD_PVQ_LAMBDA*(1./(1 << OD_BITRES))*
    enc->target_quantizer*enc->target_quantizer;
  /*The deringing filter uses yet another adjusted lambda*/
  enc->dering_lambda = 0.67*OD_PVQ_LAMBDA*
    enc->target_quantizer*enc->target_quantizer;

  OD_UNUSED(is_golden_frame);
  OD_UNUSED(frame_type);
}

int od_enc_rc_update_state(od_enc_ctx *enc, long bits,
 int is_golden_frame, int frame_type, int droppable) {
  OD_UNUSED(enc);
  OD_UNUSED(bits);
  OD_UNUSED(is_golden_frame);
  OD_UNUSED(frame_type);
  OD_UNUSED(droppable);
  /* Stub to establish API */
  return OD_EIMPL;
}

int od_enc_rc_2pass_out(od_enc_ctx *enc, unsigned char **buf) {
  OD_UNUSED(enc);
  OD_UNUSED(buf);
  if (enc->rc.target_bitrate <= 0 ||
   (enc->state.cur_time >=0 && enc->rc.twopass_state != 1)) {
    return OD_EINVAL;
  }
  /* Stub to establish API */
  return OD_EIMPL;
}

int od_enc_rc_2pass_in(od_enc_ctx *enc, unsigned char *buf, size_t bytes) {
  OD_UNUSED(enc);
  OD_UNUSED(buf);
  OD_UNUSED(bytes);
  if(enc->rc.target_bitrate <=0 ||
   (enc->state.cur_time >= 0 && enc->rc.twopass_state != 2)) {
    return OD_EINVAL;
  }
  /* Stub to establish API */
  return OD_EIMPL;
}
