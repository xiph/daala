/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

#include <string.h>
#include "encint.h"

int daala_encode_flush_header(daala_enc_ctx *_enc, daala_comment *_dc,
 ogg_packet *_op) {
  daala_info *info = &_enc->state.info;
  if (_enc == NULL || _op == NULL) return OD_EFAULT;
  switch (_enc->packet_state) {
    case OD_PACKET_INFO_HDR:
    {
      int pli;
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb, 0x80);
      oggbyte_writecopy(&_enc->obb, "daala", 5);
      oggbyte_write1(&_enc->obb, info->version_major);
      oggbyte_write1(&_enc->obb, info->version_minor);
      oggbyte_write1(&_enc->obb, info->version_sub);
      OD_ASSERT(info->pic_width > 0);
      oggbyte_write4(&_enc->obb, info->pic_width);
      OD_ASSERT(info->pic_height > 0);
      oggbyte_write4(&_enc->obb, info->pic_height);
      oggbyte_write4(&_enc->obb, info->pixel_aspect_numerator);
      oggbyte_write4(&_enc->obb, info->pixel_aspect_denominator);
      oggbyte_write4(&_enc->obb, info->timebase_numerator);
      oggbyte_write4(&_enc->obb, info->timebase_denominator);
      oggbyte_write4(&_enc->obb, info->frame_duration);
      OD_ASSERT(info->keyframe_granule_shift < 32);
      oggbyte_write1(&_enc->obb, info->keyframe_granule_shift);
      OD_ASSERT((info->nplanes >= 1) && (info->nplanes <= OD_NPLANES_MAX));
      oggbyte_write1(&_enc->obb, info->nplanes);
      for (pli = 0; pli < info->nplanes; ++pli) {
        oggbyte_write1(&_enc->obb, info->plane_info[pli].xdec);
        oggbyte_write1(&_enc->obb, info->plane_info[pli].ydec);
      }
      _op->b_o_s = 1;
    }
    break;
    case OD_PACKET_COMMENT_HDR:
    {
      const char *vendor;
      ogg_uint32_t  vendor_len;
      int           ci;
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb, 0x81);
      oggbyte_writecopy(&_enc->obb, "daala", 5);
      vendor = daala_version_string();
      vendor_len = strlen(vendor);
      oggbyte_write4(&_enc->obb, vendor_len);
      oggbyte_writecopy(&_enc->obb, vendor, vendor_len);
      oggbyte_write4(&_enc->obb, _dc->comments);
      for (ci = 0; ci < _dc->comments; ci++) {
        if (_dc->user_comments[ci] != NULL) {
          oggbyte_write4(&_enc->obb, _dc->comment_lengths[ci]);
          oggbyte_writecopy(&_enc->obb, _dc->user_comments[ci],
           _dc->comment_lengths[ci]);
        }
        else oggbyte_write4(&_enc->obb, 0);
      }
      _op->b_o_s = 0;
    }
    break;
    case OD_PACKET_SETUP_HDR:
    {
      /*TODO: Produce header contents.*/
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb, 0x82);
      oggbyte_writecopy(&_enc->obb, "daala", 5);
      _op->b_o_s = 0;
    }
    break;
    /*No more headers to emit.*/
    default: return 0;
  }
  /*This is kind of fugly: we hand the user a buffer which they do not own.
    We will overwrite it when the next packet is output, so the user better be
     done with it by then.
    Vorbis is little better: it hands back buffers that it will free the next
     time the headers are requested, or when the encoder is cleared.
    Hopefully libogg2 will make this much cleaner.*/
  _op->packet = oggbyte_get_buffer(&_enc->obb);
  _op->bytes = oggbyte_bytes(&_enc->obb);
  _op->e_o_s = 0;
  _op->granulepos = 0;
  /*Is this smart? Vorbis does not even set this field.*/
  _op->packetno = 0;
  return ++_enc->packet_state+3;
}
