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

