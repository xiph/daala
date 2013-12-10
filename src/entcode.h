/*Daala video codec
Copyright (c) 2001-2013 Daala project contributors.  All rights reserved.

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

#if !defined(_entcode_H)
# define _entcode_H (1)
# include <limits.h>
# include <stddef.h>
# include "internal.h"

/*OPT: od_ec_window must be at least 32 bits, but if you have fast arithmetic
   on a larger type, you can speed up the decoder by using it here.*/
typedef ogg_uint32_t od_ec_window;

# define OD_EC_WINDOW_SIZE ((int)sizeof(od_ec_window)*CHAR_BIT)

/*The number of bits to use for the range-coded part of unsigned integers.*/
# define OD_EC_UINT_BITS (4)

/*The resolution of fractional-precision bit usage measurements, i.e.,
   3 => 1/8th bits.*/
# define OD_BITRES (3)

extern const ogg_uint16_t OD_UNIFORM_CDFS_Q15[135];

/*Returns a Q15 CDF for a uniform probability distribution of the given size.
  n: The size of the distribution.
     This must be at least 2, and no more than 16.*/
# define OD_UNIFORM_CDF_Q15(n) \
   (OD_UNIFORM_CDFS_Q15 + ((n)*((n) - 1) >> 1) - 1)

/*See entcode.c for further documentation.*/

OD_WARN_UNUSED_RESULT ogg_uint32_t od_ec_tell_frac(ogg_uint32_t nbits_total,
 ogg_uint32_t rng);

#endif
