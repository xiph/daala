/* Copyright (c) 2001-2011 Timothy B. Terriberry
   Copyright (c) 2008-2009 Xiph.Org Foundation */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#if !defined(_entdec_H)
# define _entdec_H (1)
# include <limits.h>
# include "entcode.h"



typedef struct od_ec_ctx od_ec_dec;



/*See entdec.c for further documentation.*/


void od_ec_dec_init(od_ec_dec *_this,
 unsigned char *_buf,ogg_uint32_t _storage);

unsigned od_ec_decode(od_ec_dec *_this,unsigned _ft);
unsigned od_ec_decode_bin(od_ec_dec *_this,unsigned _ftb);

void od_ec_dec_update(od_ec_dec *_this,unsigned _fl,unsigned _fh,unsigned _ft);

int od_ec_dec_bool(od_ec_dec *_this,unsigned _fz);
int od_ec_dec_icdf_ft(od_ec_dec *_this,
 const unsigned char *_icdf,unsigned _ft);
int od_ec_dec_icdf16_ft(od_ec_dec *_this,
 const ogg_uint16_t *_icdf,unsigned _ft);
int od_ec_dec_icdf(od_ec_dec *_this,const unsigned char *_icdf,unsigned _ftb);
int od_ec_dec_icdf16(od_ec_dec *_this,const ogg_uint16_t *_icdf,unsigned _ftb);

ogg_uint32_t od_ec_dec_uint(od_ec_dec *_this,ogg_uint32_t _ft);

ogg_uint32_t od_ec_dec_bits(od_ec_dec *_this,unsigned _ftb);

#endif
