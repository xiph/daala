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



typedef struct od_ec_dec od_ec_dec;



/*The entropy decoder context.*/
struct od_ec_dec{
  /*The start of the current input buffer.*/
  const unsigned char *buf;
  /*The read pointer for the raw bits.*/
  const unsigned char *eptr;
  /*Bits that will be read from/written at the end.*/
  od_ec_window         end_window;
  /*Number of valid bits in end_window.*/
  int                  nend_bits;
  /*An offset used to keep track of tell after reaching the end of the stream.
    This is constant throughout most of the decoding process, but becomes
     important once we hit the end of the buffer and stop incrementing pointers
     (and instead pretend cnt/nend_bits have lots of bits).*/
  ogg_int32_t          tell_offs;
  /*The end of the current input buffer.*/
  const unsigned char *end;
  /*The read pointer for the entropy-coded bits.*/
  const unsigned char *bptr;
  /*The difference between the coded value and the low end of the current
     range.*/
  od_ec_window         dif;
  /*The number of values in the current range.*/
  ogg_uint16_t         rng;
  /*The number of bits of data in the current value.*/
  ogg_int16_t          cnt;
  /*The saved normalization factor from od_ec_decode().*/
  ogg_uint32_t         ext;
  /*Nonzero if an error occurred.*/
  int                  error;
};




/*See entdec.c for further documentation.*/


void od_ec_dec_init(od_ec_dec *_this,
 const unsigned char *_buf,ogg_uint32_t _storage);

unsigned od_ec_decode(od_ec_dec *_this,unsigned _ft);
unsigned od_ec_decode_q15(od_ec_dec *_this);
unsigned od_ec_decode_unscaled(od_ec_dec *_this,unsigned _ft);
unsigned od_ec_decode_unscaled_dyadic(od_ec_dec *_this,unsigned _ftb);

void od_ec_dec_update(od_ec_dec *_this,unsigned _fl,unsigned _fh,unsigned _ft);

int od_ec_decode_bool(od_ec_dec *_this,unsigned _fz,unsigned _ft);
int od_ec_decode_bool_q15(od_ec_dec *_this,unsigned _fz);
int od_ec_decode_cdf(od_ec_dec *_this,const ogg_uint16_t *_cdf,int _nsyms);
int od_ec_decode_cdf_q15(od_ec_dec *_this,const ogg_uint16_t *_cdf,int _nsyms);
int od_ec_decode_cdf_unscaled(od_ec_dec *_this,
 const ogg_uint16_t *_cdf,int _nsyms);
int od_ec_decode_cdf_unscaled_dyadic(od_ec_dec *_this,
 const ogg_uint16_t *_cdf,int _nsyms,unsigned _ftb);

ogg_uint32_t od_ec_dec_uint(od_ec_dec *_this,ogg_uint32_t _ft);

ogg_uint32_t od_ec_dec_bits(od_ec_dec *_this,unsigned _ftb);

int od_ec_dec_tell(od_ec_dec *_this);
ogg_uint32_t od_ec_dec_tell_frac(od_ec_dec *_this);

#endif
