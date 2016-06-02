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

void od_enc_rc_select_quantizers_and_lambdas(od_enc_ctx *enc,
 int is_golden_frame, int frame_type) {
  /*For now, use the quantizer selected by legacy code.*/
  enc->target_quantizer = enc->state.quantizer;
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
