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
#if !defined(_entenc_H)
# define _entenc_H (1)
# include <stddef.h>
# include "entcode.h"



typedef struct od_ec_enc od_ec_enc;



struct od_ec_enc{
  /*The common encoder/decoder state.*/
  od_ec_ctx     base;
  /*A buffer for output bytes with their associated carry flags.*/
  ogg_uint16_t *precarry_buf;
  /*The size of the pre-carry buffer.*/
  ogg_uint32_t  precarry_storage;
};


/*See entenc.c for further documentation.*/


void od_ec_enc_init(od_ec_enc *_this,ogg_uint32_t _size);
void od_ec_enc_reset(od_ec_enc *_this);
void od_ec_enc_clear(od_ec_enc *_this);

void od_ec_encode(od_ec_enc *_this,unsigned _fl,unsigned _fh,unsigned _ft);
void od_ec_encode_bin(od_ec_enc *_this,
 unsigned _fl,unsigned _fh,unsigned _ftb);
void od_ec_enc_bool(od_ec_enc *_this,int _val,unsigned _fz);
void od_ec_enc_icdf_ft(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ft);
void od_ec_enc_icdf(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ftb);
void od_ec_enc_icdf16(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_icdf,unsigned _ftb);
void od_ec_enc_icdf16_ft(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_icdf,unsigned _ft);

void od_ec_enc_uint(od_ec_enc *_this,ogg_uint32_t _fl,ogg_uint32_t _ft);

void od_ec_enc_bits(od_ec_enc *_this,ogg_uint32_t _fl,unsigned _ftb);

void od_ec_enc_patch_initial_bits(od_ec_enc *_this,
 unsigned _val,unsigned _nbits);
unsigned char *od_ec_enc_done(od_ec_enc *_this,ogg_uint32_t *_nbytes);

#endif
