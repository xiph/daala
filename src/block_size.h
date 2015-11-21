/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_block_size_h)
# define _block_size_h

# include "odintrin.h"
# include "../include/daala/daala_integer.h"
# include "entenc.h"

# define OD_BLOCK_SIZE4x4(bsize, bstride, bx, by) \
   OD_BLOCK_SIZE8x8(bsize, bstride, (bx) >> 1, (by) >> 1)
# define OD_BLOCK_SIZE8x8(bsize, bstride, bx, by) \
   ((bsize)[(by)*(bstride) + (bx)])

/*Possible block sizes, note that OD_BLOCK_NXN = log2(N) - 2.*/
#define OD_BLOCK_4X4 (0)
#define OD_BLOCK_8X8 (1)
#define OD_BLOCK_16X16 (2)
#define OD_BLOCK_32X32 (3)
#define OD_BLOCK_SIZES (OD_BLOCK_32X32+1)

#endif
