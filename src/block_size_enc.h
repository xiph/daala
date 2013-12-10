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

#if !defined(_block_size_enc_h)
# define _block_size_enc_h

# include "entenc.h"

/* None of these values should be larger than OFF32 or else the sun will
   explode */
# define OFF8  (1)
# define OFF16 (2)
# define OFF32 (3)

# define COUNT8  (3 + 2*OFF8)
# define COUNT16 (7 + 2*OFF16)
# define COUNT32 (15 + 2*OFF32)


/* None of these values should be larger than OFF32/2 of else the universe will
   collapse */
# define OFF8_16 (1)
# define OFF8_32 (1)

# define COUNT8_16  (3 + 2*OFF8_16)
# define COUNT8_32  (7 + 2*OFF8_32)

# define SIZE2_SUMS (16 + 2*OFF32)
# define SIZE4_SUMS (15 + 2*OFF32)
# define SIZE8_SUMS (7 + 2*OFF8_32)

/*#define STRIDE (32+4*OFF32)*/
# define BLOCK_OFFSET(stride) ((2*OFF32)*(stride) + (2*OFF32))

/* This struct can be made a lot smaller by using temporary values,
   but as it is, it's much easier to debug and modify */
typedef struct {
  ogg_int32_t Sx2[SIZE2_SUMS][SIZE2_SUMS];
  ogg_int32_t Sxx2[SIZE2_SUMS][SIZE2_SUMS];
  ogg_int32_t Sx4[SIZE4_SUMS][SIZE4_SUMS];
  ogg_int32_t Sxx4[SIZE4_SUMS][SIZE4_SUMS];
  ogg_int32_t Sx8[SIZE8_SUMS][SIZE8_SUMS];
  ogg_int32_t Sxx8[SIZE8_SUMS][SIZE8_SUMS];
  ogg_int32_t Var4[SIZE4_SUMS][SIZE4_SUMS];
  ogg_int32_t invVar4[SIZE4_SUMS][SIZE4_SUMS];
  ogg_int32_t Var8[SIZE8_SUMS][SIZE8_SUMS];
  ogg_int32_t invVar8[SIZE8_SUMS][SIZE8_SUMS];
} BlockStats;

typedef struct {
  BlockStats img_stats;
  BlockStats psy_stats;

  unsigned char res[2*SIZE2_SUMS][2*SIZE2_SUMS];

  /* 4x4 metrics */
  ogg_int32_t noise4_4[8][8];
  ogg_int32_t noise4_8[4][4];
  ogg_int32_t noise4_16[2][2];
  ogg_int32_t noise4_32;

  /* 8x8 metrics */
  ogg_int32_t noise8_16[2][2];
  ogg_int32_t noise8_32;

  float psy4[8][8];
  float psy8[4][4];
  float psy16[2][2];
  float psy32;

  float cg8[4][4];
  float cg16[2][2];
  float cg32;

  float dec_gain8[4][4];
  float dec_gain16[2][2];
} BlockSizeComp;

void process_block_size32(BlockSizeComp *bs, const unsigned char *psy_img,
 int stride, const unsigned char *img, int pred_stride, int dec[4][4]);

void od_block_size_encode(od_ec_enc *enc,
 const unsigned char *bsize, int stride);

#endif
