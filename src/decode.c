/*
    Daala video codec
    Copyright (C) 2010 Timothy B. Terriberry and Daala project contributors

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


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "decint.h"

static int od_dec_init(od_dec_ctx *_dec,const daala_info *_info,
 const daala_setup_info *_setup){
  int ret;
  ret=od_state_init(&_dec->state,_info);
  if(ret<0)return ret;
  _dec->packet_state=OD_PACKET_INFO_HDR;
  return 0;
}

static void od_dec_clear(od_dec_ctx *_dec){
  od_state_clear(&_dec->state);
}

daala_dec_ctx *daala_decode_alloc(const daala_info *_info,
 const daala_setup_info *_setup){
  od_dec_ctx *dec;
  if(_info==NULL)return NULL;
  dec=(od_dec_ctx *)_ogg_malloc(sizeof(*dec));
  if(od_dec_init(dec,_info,_setup)<0){
    _ogg_free(dec);
    return NULL;
  }
  return dec;
}

void daala_decode_free(daala_dec_ctx *_dec){
  if(_dec!=NULL){
    od_dec_clear(_dec);
    _ogg_free(_dec);
  }
}

int daala_decode_ctl(daala_dec_ctx *_dec,int _req,void *_buf,
 size_t _buf_sz){
  switch(_req){
    default:return OD_EIMPL;
  }
}
