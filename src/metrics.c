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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include "logging.h"
#include "metrics.h"

static const char *OD_MET_CATEGORY_NAMES[OD_MET_NCATS] = {
  "technique",
  "plane"
};

static const char *OD_MET_TECHNIQUE_NAMES[OD_MET_NTECHS] = {
  "unknown",
  "frame",
  "block-size",
  "intra-mode",
  "dc-coeff",
  "ac-coeffs",
  "motion-vectors"
};

static const char *OD_MET_PLANE_NAMES[OD_MET_NPLANES] = {
  "unknown",
  "frame",
  "luma",
  "cb",
  "cr",
  "alpha"
};

static const char * const *OD_MET_CATEGORY_VALUE_NAMES[OD_MET_NCATS] = {
  OD_MET_TECHNIQUE_NAMES,
  OD_MET_PLANE_NAMES
};

static const unsigned int OD_MET_INDICES[OD_MET_NCATS] = {
  OD_MET_NTECHS,
  OD_MET_NPLANES
};

void od_metrics_init(od_metrics *metrics) {
  char  fname[1024];
  char *pre;
  char *suf;
  int   rv;
  pre = "metrics-";
  suf = getenv("OD_METRICS_SUFFIX");
  if (!suf) {
    pre = "metrics";
    suf = "";
  }
  rv = snprintf(fname, sizeof(fname), "%s%s.json", pre, suf);
  OD_ASSERT(rv >= 0 && ((size_t)rv) < sizeof(fname));
  metrics->fp = fopen(fname, "w");
  OD_ASSERT(metrics->fp);
  od_metrics_reset(metrics);
}

void od_metrics_clear(od_metrics *metrics) {
  OD_ASSERT(metrics->fp);
  fclose(metrics->fp);
}

void od_metrics_reset(od_metrics *metrics) {
  int i;
  /* Calling od_ec_enc_tell_frac() on a reset od_ec_enc struct returns 8. */
  metrics->last_frac_bits = 8;
  /* Set the initial state for each category to unknown. */
  for (i = 0; i < OD_MET_NCATS; i++) {
    metrics->state[i] = 0;
  }
  for (i = 0; i < OD_METRICS_SIZE; i++) {
    metrics->frac_bits[i] = 0;
  }
}

static int od_metrics_index(unsigned int state[OD_MET_NCATS]) {
  int index;
  od_metrics_category cat;
  index = state[0];
  for (cat = 1; cat < OD_MET_NCATS; cat++) {
    index = (index*OD_MET_INDICES[OD_MET_NCATS - cat]) + state[cat];
  }
  return index;
}

void od_metrics_update_frac_bits(od_metrics *metrics, ogg_uint32_t frac_bits) {
  ogg_uint32_t frac_bits_diff;
  frac_bits_diff = frac_bits - metrics->last_frac_bits;
  metrics->frac_bits[od_metrics_index(metrics->state)] += frac_bits_diff;
  metrics->last_frac_bits = frac_bits;
}

void od_metrics_set_category(od_metrics *metrics, od_metrics_category cat,
 unsigned int value) {
  OD_ASSERT(cat < OD_MET_NCATS);
  OD_ASSERT(value < OD_MET_INDICES[cat]);
  metrics->state[cat] = value;
}

void od_metrics_update(od_metrics *metrics, ogg_uint32_t frac_bits,
 od_metrics_category cat, unsigned int value) {
  od_metrics_update_frac_bits(metrics, frac_bits);
  od_metrics_set_category(metrics, cat, value);
}

static int od_metrics_next_state(unsigned int state[OD_MET_NCATS],
 od_metrics_category skip) {
  od_metrics_category i;
  i = skip == 0;
  while (i < OD_MET_NCATS) {
    state[i]++;
    if (state[i] < OD_MET_INDICES[i]) {
      return 1;
    }
    state[i] = 0;
    i++;
    if (skip != OD_MET_NCATS && i == skip) {
      i++;
    }
  }
  return 0;
}

static ogg_uint32_t od_metrics_get_total(od_metrics *metrics,
 od_metrics_category cat, unsigned int value) {
  unsigned int state[OD_MET_NCATS];
  ogg_uint32_t total;
  od_metrics_category i;
  for (i = 0; i < OD_MET_NCATS; i++) {
    state[i] = 0;
  }
  state[cat] = value;
  total = 0;
  do {
    total += metrics->frac_bits[od_metrics_index(state)];
  }
  while (od_metrics_next_state(state, cat));
  return total;
}

void od_metrics_print_state(FILE *_fp, unsigned int state[OD_MET_NCATS]) {
  int cat;
  for (cat = 0; cat < OD_MET_NCATS; cat++) {
    fprintf(_fp, "%s%s", cat > 0 ? "," : "",
     OD_MET_CATEGORY_VALUE_NAMES[cat][state[cat]]);
  }
}

void od_metrics_print(od_metrics *metrics, FILE *_fp) {
  unsigned int state[OD_MET_NCATS];
  int i;
  for (i = 0; i < OD_MET_NCATS; i++) {
    state[i] = 0;
  }
  do {
    int index;
    index = od_metrics_index(state);
    od_metrics_print_state(_fp, state);
    fprintf(_fp, " (%i): %i\n", index, metrics->frac_bits[index]);
  }
  while (od_metrics_next_state(state, OD_MET_NCATS));
}

void od_metrics_write(od_metrics *metrics, ogg_int64_t cur_time) {
  od_metrics_category cat;
  unsigned int value;
  long fsize;
  OD_ASSERT(metrics->fp);
  fsize = ftell(metrics->fp);
  if (fsize == 0) {
    fprintf(metrics->fp, "[");
  }
  else {
    fseek(metrics->fp, fsize - 1, SEEK_SET);
    fprintf(metrics->fp, ",\n");
  }
  fprintf(metrics->fp, "{\n");
  fprintf(metrics->fp, "  \"frame\": %" OD_I64FMT ",\n", (long long)cur_time);
  fprintf(metrics->fp, "  \"total\": %u,\n", metrics->last_frac_bits-8);
  for (cat = 0; cat < OD_MET_NCATS; cat++) {
    fprintf(metrics->fp, "%s  \"%s\": {\n", cat > 0 ? ",\n" : "",
     OD_MET_CATEGORY_NAMES[cat]);
    for (value = 0; value < OD_MET_INDICES[cat]; value++) {
      fprintf(metrics->fp, "%s    \"%s\": %i", value > 0 ? ",\n" : "",
       OD_MET_CATEGORY_VALUE_NAMES[cat][value],
       od_metrics_get_total(metrics, cat, value));
    }
    fprintf(metrics->fp, "\n  }");
  }
  fprintf(metrics->fp, "\n}]");
}
