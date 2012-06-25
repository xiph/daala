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
#if !defined(_entcode_H)
# define _entcode_H (1)
# include <limits.h>
# include <stddef.h>
# include "internal.h"



/*OPT: od_ec_window must be at least 32 bits, but if you have fast arithmetic
   on a larger type, you can speed up the decoder by using it here.*/
typedef ogg_uint32_t     od_ec_window;
typedef struct od_ec_ctx od_ec_ctx;



# define OD_EC_WINDOW_SIZE ((int)sizeof(od_ec_window)*CHAR_BIT)

/*The number of bits to use for the range-coded part of unsigned integers.*/
# define OD_EC_UINT_BITS   (4)

/*The resolution of fractional-precision bit usage measurements, i.e.,
   3 => 1/8th bits.*/
# define OD_BITRES         (3)



/*The entropy encoder/decoder context.
  This serves both as a decoder context and as the base for an encoder context,
   so that common functions like od_ec_tell() can be used on either one.*/
struct od_ec_ctx{
   /*Buffered input/output.*/
   unsigned char *buf;
   /*The size of the buffer.*/
   ogg_uint32_t   storage;
   /*The offset at which the last byte containing raw bits was read/written.*/
   ogg_uint32_t   end_offs;
   /*Bits that will be read from/written at the end.*/
   od_ec_window   end_window;
   /*Number of valid bits in end_window.*/
   int            nend_bits;
   /*The total number of whole bits read/written.
     This does not include partial bits currently in the range coder.*/
   int            nbits_total;
   /*The offset at which the next range coder byte will be read/written.*/
   ogg_uint32_t   offs;
   /*In the decoder: the difference between the top of the current range and
      the input value, minus one.
     In the encoder: the low end of the current range.*/
   od_ec_window   val;
   /*The number of values in the current range.*/
   ogg_uint16_t   rng;
   /*The number of bits of data in the current value.*/
   ogg_int16_t    cnt;
   /*Decoder-only: the saved normalization factor from od_ec_decode().*/
   ogg_uint32_t   ext;
   /*Nonzero if an error occurred.*/
   int            error;
};

#define od_ec_range_bytes(/*od_ec_ctx **/_this) \
 (((od_ec_ctx *)(_this))->offs)

#define od_ec_get_error(/*od_ec_ctx **/_this) \
 (((od_ec_ctx *)(_this))->error)

/*Returns the number of bits "used" by the encoded or decoded symbols so far.
  This same number can be computed in either the encoder or the decoder, and is
   suitable for making coding decisions.
  Return: The number of bits.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
#define od_ec_tell(/*od_ec_ctx **/_this) \
  (((od_ec_ctx *)(_this))->nbits_total)


/*See entcode.c for further documentation.*/


ogg_uint32_t od_ec_tell_frac(od_ec_ctx *_this);

#endif
