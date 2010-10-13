/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#include <string.h>
#include "encint.h"

int daala_encode_flush_header(daala_enc_ctx *_enc,daala_comment *_dc,
 ogg_packet *_op){
  if(_enc==NULL||_op==NULL)return OD_EFAULT;
  /*TODO: Produce header contents.*/
  switch(_enc->packet_state){
    case OD_PACKET_INFO_HDR:{
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb,0x80);
      oggbyte_writecopy(&_enc->obb,"daala",5);
      _op->b_o_s=1;
    }break;
    case OD_PACKET_COMMENT_HDR:{
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb,0x81);
      oggbyte_writecopy(&_enc->obb,"daala",5);
      _op->b_o_s=0;
    }break;
    case OD_PACKET_SETUP_HDR:{
      oggbyte_reset(&_enc->obb);
      oggbyte_write1(&_enc->obb,0x82);
      oggbyte_writecopy(&_enc->obb,"daala",5);
      _op->b_o_s=0;
    }break;
    /*No more headers to emit.*/
    default:return 0;
  }
  /*This is kind of fugly: we hand the user a buffer which they do not own.
    We will overwrite it when the next packet is output, so the user better be
     done with it by then.
    Vorbis is little better: it hands back buffers that it will free the next
     time the headers are requested, or when the encoder is cleared.
    Hopefully libogg2 will make this much cleaner.*/
  _op->packet=oggbyte_get_buffer(&_enc->obb);
  _op->bytes=oggbyte_bytes(&_enc->obb);
  _op->e_o_s=0;
  _op->granulepos=0;
  /*Is this smart? Vorbis does not even set this field.*/
  _op->packetno=0;
  return ++_enc->packet_state+3;
}
