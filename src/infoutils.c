/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

#include <string.h>
#include "encint.h"

#define DAALA_PIX_FMT_NUM 20

struct daala_pix_fmt {
  const char *name;
  int bitdepth;
  int nplanes;
  int xdec[4];
  int ydec[4];
};

static const struct daala_pix_fmt pix_fmt_list[DAALA_PIX_FMT_NUM] = {
  { "444",        8, 3, {0, 0, 0},    {0, 0, 0}    },
  { "444p10",    10, 3, {0, 0, 0},    {0, 0, 0}    },
  { "444p12",    12, 3, {0, 0, 0},    {0, 0, 0}    },
  { "444p14",    14, 3, {0, 0, 0},    {0, 0, 0}    },
  { "444p16",    16, 3, {0, 0, 0},    {0, 0, 0}    },
  { "444alpha",   8, 4, {0, 0, 0, 0}, {0, 0, 0, 0} },
  { "422",        8, 3, {0, 1, 0},    {0, 1, 0}    },
  { "422p10",    10, 3, {0, 1, 0},    {0, 1, 0}    },
  { "422p12",    12, 3, {0, 1, 0},    {0, 1, 0}    },
  { "422p14",    14, 3, {0, 1, 0},    {0, 1, 0}    },
  { "422p16",    16, 3, {0, 1, 0},    {0, 1, 0}    },
  { "411",        8, 3, {0, 2, 0},    {0, 2, 0}    },
  { "420",        8, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420jpeg",    8, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420mpeg2",   8, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420paldv",   8, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420p10",    10, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420p12",    12, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420p14",    14, 3, {0, 1, 1},    {0, 1, 1}    },
  { "420p16",    16, 3, {0, 1, 1},    {0, 1, 1}    },
};

static int daala_lookup_bitdepth_code(int bitdepth)
{
  if (bitdepth == 8)
    return OD_BITDEPTH_MODE_8;
  else if (bitdepth == 10)
    return OD_BITDEPTH_MODE_10;
  else
    return OD_BITDEPTH_MODE_12;
}

/* Set format from a code */
int daala_set_pix_info(daala_info *info, unsigned char code)
{
  int i;
  if (code < 0 || code > DAALA_PIX_FMT_NUM)
    return 1;
  info->bitdepth_mode = daala_lookup_bitdepth_code(pix_fmt_list[code].bitdepth);
  info->nplanes = pix_fmt_list[code].nplanes;
  for (i = 0; i < info->nplanes; ++i) {
    info->plane_info[i].xdec = pix_fmt_list[code].xdec[i];
    info->plane_info[i].ydec = pix_fmt_list[code].ydec[i];
  }
  return 0;
}

/* Look up name and set info, used by the encoder */
unsigned char daala_lookup_pix_fmt(const char *name)
{
  int i;
  for (i = 0; i < DAALA_PIX_FMT_NUM; i++) {
    if (!strncmp(pix_fmt_list[i].name, name, strlen(name))) {
      return i;
    }
  }
  return 255;
}
