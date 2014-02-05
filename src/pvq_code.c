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

#include "pvq_code.h"

/* Table of combined "use prediction" pvq flags for 8x8 trained on
   subset1. Should eventually make this adaptive. */
const ogg_uint16_t pred8_cdf[16] = {
  22313, 22461, 22993, 23050, 23418, 23468, 23553, 23617,
  29873, 30181, 31285, 31409, 32380, 32525, 32701, 32768
};

const ogg_uint16_t pred16_cdf[16][8] = {
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 },
  { 4096,  8192, 12288, 16384, 20480, 24576, 28672, 32768 }
};
