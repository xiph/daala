/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "filter.h"
#include "partition.h"

/* The tables below specify how coefficient blocks are translated to
   and from PVQ partition coding scan order for 4x4, 8x8 and 16x16 */
static const int od_layout16_offsets[4] = { 0, 32, 64, 192 };
extern const index_pair od_zigzag16[];
const band_layout od_layout16 = {
  od_zigzag16,
  16,
  3,
  od_layout16_offsets
};

const int od_layout8_offsets[4] = { 0, 8, 16, 48 };
extern const index_pair od_zigzag8[];
const band_layout od_layout8 = {
  od_zigzag8,
  8,
  3,
  od_layout8_offsets
};

static const int od_layout4_offsets[2] = { 0, 15 };
extern const index_pair od_zigzag4[];
const band_layout od_layout4 = {
  od_zigzag4,
  4,
  1,
  od_layout4_offsets
};

/** Perform a single stage of conversion from a coefficient block in
 * raster order into coding scan order
 *
 * @param [in]     layout  scan order specification
 * @param [out]    dst     destination vector
 * @param [in]     src     source coefficient block
 * @param [int]    int     source vector row stride
 */
static void od_band_from_raster(const band_layout *layout, od_coeff *dst,
 od_coeff *src, int stride) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    dst[i] = src[layout->dst_table[i][1]*stride + layout->dst_table[i][0]];
  }
}

/** Perform a single stage of conversion from a vector in coding scan
    order back into a coefficient block in raster order
 *
 * @param [in]     layout  scan order specification
 * @param [out]    dst     destination coefficient block
 * @param [in]     src     source vector
 * @param [int]    stride  destination vector row stride
 */
static void od_raster_from_band(const band_layout *layout, od_coeff *dst,
 int stride, od_coeff *src) {
  int i;
  int len;
  len = layout->band_offsets[layout->nb_bands];
  for (i = 0; i < len; i++) {
    dst[layout->dst_table[i][1]*stride + layout->dst_table[i][0]] = src[i];
  }
}

/** Converts a coefficient block in raster order into a vector in
 * coding scan order with the PVQ partitions laid out one after
 * another.  This works in stages; the 4x4 conversion is applied to
 * the coefficients nearest DC, then the 8x8 applied to the 8x8 block
 * nearest DC that was not already coded by 4x4, then 16x16 following
 * the same pattern.
 *
 * @param [out]    dst        destination vector
 * @param [in]     n          block size (along one side)
 * @param [in]     src        source coefficient block
 * @param [in]     stride     source vector row stride
 * @param [in]     interleave interleaves entries for the scalar
                              (non-pvq) case
 */
void od_raster_to_coding_order(od_coeff *dst,  int n, od_coeff *src, int stride,
 int interleave) {
  od_coeff tmp1[1024];
  od_band_from_raster(&od_layout4, dst+1, src, stride);
  if (n >= 8) {
    int i;
    od_band_from_raster(&od_layout8, dst+16, src, stride);
    if (interleave) {
      for (i = 0; i < 8; i++) {
        tmp1[2*i] = dst[16+i];
        tmp1[2*i+1] = dst[24+i];
      }
      for (i = 0; i < 16; i++) {
        dst[16+i] = tmp1[i];
      }
    }
  }
  if (n >= 16) {
    int i;
    od_band_from_raster(&od_layout16, dst+64, src, stride);
    if (interleave) {
      for (i = 0; i < 32; i++) {
        tmp1[2*i] = dst[64+i];
        tmp1[2*i+1] = dst[96+i];
      }
      for (i = 0; i < 64; i++) {
        dst[64+i] = tmp1[i];
      }
    }
  }
  dst[0] = src[0];
}

/** Converts a vector in coding scan order witht he PVQ partitions
 * laid out one after another into a coefficient block in raster
 * order. This works in stages in the reverse order of raster->scan
 * order; the 16x16 conversion is applied to the coefficients that
 * don't appear in an 8x8 block, then the 8x8 applied to the 8x8 block
 * sans the 4x4 block it contains, then 4x4 is converted sans DC.
 *
 * @param [out]    dst        destination coefficient block
 * @param [in]     stride     destination vector row stride
 * @param [in]     src        source vector
 * @param [in]     n          block size (along one side)
 * @param [in]     interleave de-interleaves entries for
                              the scalar (non-pvq) case
 */
void od_coding_order_to_raster(od_coeff *dst,  int stride, od_coeff *src,
 int n, int interleave) {
  od_raster_from_band(&od_layout4, dst, stride, src+1);
  if (n >= 8) {
    if (interleave) {
      int i;
      od_coeff tmp1[1024];
      for (i = 0; i < 16; i++) {
        tmp1[i] = src[16 + i];
      }
      for (i = 0; i < 8; i++) {
        src[16+i] = tmp1[2*i];
        src[24+i] = tmp1[2*i + 1];
      }
    }
    od_raster_from_band(&od_layout8, dst, stride, src+16);
  }
  if (n >= 16) {
    if (interleave) {
      int i;
      od_coeff tmp1[1024];
      for (i = 0; i < 64; i++) {
        tmp1[i] = src[64 + i];
      }
      for (i = 0; i < 32; i++) {
        src[64+i] = tmp1[2*i];
        src[96+i] = tmp1[2*i + 1];
      }
    }
    od_raster_from_band(&od_layout16, dst, stride, src+64);
  }
  dst[0] = src[0];
}
