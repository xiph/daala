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


#if !defined(_metrics_H)
# define _metrics_H (1)

# include <stdio.h>
# include "internal.h"

enum od_metrics_category {
  OD_MET_CAT_TECHNIQUE = 0,
  OD_MET_CAT_PLANE = 1,
  OD_MET_NCATS
};

typedef enum od_metrics_category od_metrics_category;

enum od_metrics_technique {
  OD_MET_TECH_UNKNOWN = 0,
  OD_MET_TECH_FRAME = 1,
  OD_MET_TECH_BLOCK_SIZE = 2,
  OD_MET_TECH_INTRA_MODE = 3,
  OD_MET_TECH_DC_COEFF = 4,
  OD_MET_TECH_AC_COEFFS = 5,
  OD_MET_TECH_MOTION_VECTORS = 6,
  OD_MET_NTECHS
};

typedef enum od_metrics_technique od_metrics_technique;

enum od_metrics_plane {
  OD_MET_PLANE_UNKNOWN = 0,
  OD_MET_PLANE_FRAME = 1,
  OD_MET_PLANE_LUMA = 2,
  OD_MET_PLANE_CB = 3,
  OD_MET_PLANE_CR = 4,
  OD_MET_PLANE_ALPHA = 5,
  OD_MET_NPLANES
};

typedef enum od_metrics_plane od_metrics_plane;

/*typedef enum od_metrics_band od_metrics_band;*/

#if defined(OD_METRICS)
# define OD_METRICS_UPDATE(metrics, frac_bits, cat, value) \
 od_metrics_update(metrics, frac_bits, cat, value)
#else
# define OD_METRICS_UPDATE(metrics, frac_bits, cat, value)
#endif

#define OD_METRICS_SIZE (OD_MET_NTECHS*OD_MET_NPLANES)

typedef struct od_metrics od_metrics;

struct od_metrics {
  FILE *fp;
  ogg_uint32_t last_frac_bits;
  unsigned int state[OD_MET_NCATS];
  ogg_uint32_t frac_bits[OD_METRICS_SIZE];
};

void od_metrics_init(od_metrics *metrics);
void od_metrics_clear(od_metrics *metrics);
void od_metrics_reset(od_metrics *metrics);
void od_metrics_update_frac_bits(od_metrics *metrics, ogg_uint32_t frac_bits);
void od_metrics_set_category(od_metrics *metrics, od_metrics_category cat,
 unsigned int value);
void od_metrics_update(od_metrics *metrics, ogg_uint32_t frac_bits,
 od_metrics_category cat, unsigned int value);
void od_metrics_print(od_metrics *metrics, FILE *_fp);
void od_metrics_write(od_metrics *metrics, ogg_int64_t cur_time);

#endif
