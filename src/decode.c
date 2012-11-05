/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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
