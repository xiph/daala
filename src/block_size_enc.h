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
# include "state.h"

/*Maximum overlap returned by od_overlap_var4x4 in units of 2x2 blocks.*/
# define OD_MAX_OVERLAP (3)
/*Maximum overlap returned by od_overlap_var8x8 in units of 4x4 blocks.*/
# define OD_MAX_OVERLAP_8 (1)

/*Number of 2x2 blocks along one direction of a superblock, plus OD_MAX_OVERLAP
   extra blocks on each side to account for lapping.*/
# define OD_SIZE2_SUMS (16 + 2*OD_MAX_OVERLAP)
/*Number of overlapping 4x4 blocks along one direction of a superblock, plus
   OD_MAX_OVERLAP extra blocks on each side to account for lapping.
  The distance between two overlapping blocks is 2 pixels.*/
# define OD_SIZE4_SUMS (15 + 2*OD_MAX_OVERLAP)
/*Number of overlapping 8x8 blocks along one direction of a superblock, plus
   OD_MAX_OVERLAP_8 extra blocks on each side to account for lapping.
  The distance between two overlapping blocks is 4 pixels.*/
# define OD_SIZE8_SUMS (7 + 2*OD_MAX_OVERLAP_8)

/*Offset between the first pixel used in `od_compute_stats` and the pixel at
   offset (0,0) in the current superblock.*/
# define OD_BLOCK_OFFSET(stride)\
  ((2*OD_MAX_OVERLAP)*(stride) + (2*OD_MAX_OVERLAP))

/*This struct can be made a lot smaller by using temporary values,
  but as it is, it's much easier to debug and modify.*/
typedef struct {
  /*Sx2[OD_MAX_OVERLAP + y][OD_MAX_OVERLAP + x] =
     Sum of the values in the 2x2 block at offset (2*x, 2*y).*/
  int32_t Sx2[OD_SIZE2_SUMS][OD_SIZE2_SUMS];
  /*Sxx2[OD_MAX_OVERLAP + y][OD_MAX_OVERLAP + x] =
     Sum of the squared values in the 2x2 block at offset (2*x, 2*y).*/
  int32_t Sxx2[OD_SIZE2_SUMS][OD_SIZE2_SUMS];
  /*Sx4[OD_MAX_OVERLAP + y][OD_MAX_OVERLAP + x] =
     Sum of the values in the 4x4 block at offset (2*x, 2*y).*/
  int32_t Sx4[OD_SIZE4_SUMS][OD_SIZE4_SUMS];
  /*Sxx4[OD_MAX_OVERLAP + y][OD_MAX_OVERLAP + x] =
     Sum of the squared values in the 2x2 block at offset (2*x, 2*y).*/
  int32_t Sxx4[OD_SIZE4_SUMS][OD_SIZE4_SUMS];
  /*Sx8[OD_MAX_OVERLAP_8 + y][OD_MAX_OVERLAP_8 + x] =
     Sum of the values in the 8x8 block at offset (4*x, 4*y).*/
  int32_t Sx8[OD_SIZE8_SUMS][OD_SIZE8_SUMS];
  /*Sxx8[OD_MAX_OVERLAP_8 + y][OD_MAX_OVERLAP_8 + x] =
     Sum of the squared values in the 8x8 block at offset (4*x, 4*y).*/
  int32_t Sxx8[OD_SIZE8_SUMS][OD_SIZE8_SUMS];
  /*Var4[OD_MAX_OVERLAP + y][OD_MAX_OVERLAP + x] =
     Variance of the 4x4 block at offset (2*x, 2*y).*/
  int32_t Var4[OD_SIZE4_SUMS][OD_SIZE4_SUMS];
  /*invVar4[y][x] = 16384/Var4[y][x]*/
  int32_t invVar4[OD_SIZE4_SUMS][OD_SIZE4_SUMS];
  /*Var8[OD_MAX_OVERLAP_8 + y][OD_MAX_OVERLAP_8 + x] =
     Variance of the 8x8 block at offset (4*x, 4*y).*/
  int32_t Var8[OD_SIZE8_SUMS][OD_SIZE8_SUMS];
  /*invVar8[y][x] = 16384/Var8[y][x]*/
  int32_t invVar8[OD_SIZE8_SUMS][OD_SIZE8_SUMS];
} od_superblock_stats;

typedef struct {
  od_superblock_stats img_stats;
  od_superblock_stats psy_stats;

  signed char res[2*OD_SIZE2_SUMS][2*OD_SIZE2_SUMS];

  /* 4x4 metrics */
  int32_t noise4_4[8][8];
  int32_t noise4_8[4][4];
  int32_t noise4_16[2][2];
  int32_t noise4_32;

  /* 8x8 metrics */
  int32_t noise8_16[2][2];
  int32_t noise8_32;

  float psy4[8][8];
  float psy8[4][4];
  float psy16[2][2];
  float psy32;

  float cg8[4][4];
  float cg16[2][2];
  float cg32;

  float dec_gain8[4][4];
  float dec_gain16[2][2];
} od_block_size_comp;

void od_split_superblock(od_block_size_comp *bs,
 const unsigned char *psy_img, int stride,
 const unsigned char *pred, int pred_stride,
  int dec[OD_BSIZE_GRID][OD_BSIZE_GRID], int q);

#endif
