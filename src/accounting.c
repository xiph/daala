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
#include "accounting.h"

static const char *OD_ACCT_CATEGORY_NAMES[OD_ACCT_NCATS] = {
  "technique",
  "plane"
};

static const char *OD_ACCT_TECHNIQUE_NAMES[OD_ACCT_NTECHS] = {
  "unknown",
  "frame",
  "block-size",
  "intra-mode",
  "dc-coeff",
  "ac-coeffs",
  "motion-vectors"
};

static const char *OD_ACCT_PLANE_NAMES[OD_ACCT_NPLANES] = {
  "unknown",
  "frame",
  "luma",
  "cb",
  "cr",
  "alpha"
};

static const char * const *OD_ACCT_CATEGORY_VALUE_NAMES[OD_ACCT_NCATS] = {
  OD_ACCT_TECHNIQUE_NAMES,
  OD_ACCT_PLANE_NAMES
};

static const unsigned int OD_ACCT_INDICES[OD_ACCT_NCATS] = {
  OD_ACCT_NTECHS,
  OD_ACCT_NPLANES
};

void od_acct_init(od_acct *acct) {
  char  fname[1024];
  char *pre;
  char *suf;
  int   rv;
  pre = "acct-";
  suf = getenv("OD_ACCT_SUFFIX");
  if (!suf) {
    pre = "acct";
    suf = "";
  }
  rv = snprintf(fname, sizeof(fname), "%s%s.json", pre, suf);
  OD_ASSERT(rv >= 0 && ((size_t)rv) < sizeof(fname));
  acct->fp = fopen(fname, "w");
  OD_ASSERT(acct->fp);
  od_acct_reset(acct);
}

void od_acct_clear(od_acct *acct) {
  OD_ASSERT(acct->fp);
  fclose(acct->fp);
}

void od_acct_reset(od_acct *acct) {
  int i;
  /* Calling od_ec_enc_tell_frac() on a reset od_ec_enc struct returns 8. */
  acct->last_frac_bits = 8;
  /* Set the initial state for each category to unknown. */
  for (i = 0; i < OD_ACCT_NCATS; i++) {
    acct->state[i] = 0;
  }
  for (i = 0; i < OD_ACCT_SIZE; i++) {
    acct->frac_bits[i] = 0;
  }
}

static int od_acct_index(unsigned int state[OD_ACCT_NCATS]) {
  int index;
  od_acct_category cat;
  index = state[0];
  for (cat = 1; cat < OD_ACCT_NCATS; cat++) {
    index = (index*OD_ACCT_INDICES[OD_ACCT_NCATS - cat]) + state[cat];
  }
  return index;
}

void od_acct_update_frac_bits(od_acct *acct, ogg_uint32_t frac_bits) {
  ogg_uint32_t frac_bits_diff;
  frac_bits_diff = frac_bits - acct->last_frac_bits;
  acct->frac_bits[od_acct_index(acct->state)] += frac_bits_diff;
  acct->last_frac_bits = frac_bits;
}

void od_acct_set_category(od_acct *acct, od_acct_category cat,
 unsigned int value) {
  OD_ASSERT(cat < OD_ACCT_NCATS);
  OD_ASSERT(value < OD_ACCT_INDICES[cat]);
  acct->state[cat] = value;
}

void od_acct_update(od_acct *acct, ogg_uint32_t frac_bits,
 od_acct_category cat, unsigned int value) {
  od_acct_update_frac_bits(acct, frac_bits);
  od_acct_set_category(acct, cat, value);
}

static int od_acct_next_state(unsigned int state[OD_ACCT_NCATS],
 od_acct_category skip) {
  od_acct_category i;
  i = skip == 0;
  while (i < OD_ACCT_NCATS) {
    state[i]++;
    if (state[i] < OD_ACCT_INDICES[i]) {
      return 1;
    }
    state[i] = 0;
    i++;
    if (skip != OD_ACCT_NCATS && i == skip) {
      i++;
    }
  }
  return 0;
}

static ogg_uint32_t od_acct_get_total(od_acct *acct,
 od_acct_category cat, unsigned int value) {
  unsigned int state[OD_ACCT_NCATS];
  ogg_uint32_t total;
  od_acct_category i;
  for (i = 0; i < OD_ACCT_NCATS; i++) {
    state[i] = 0;
  }
  state[cat] = value;
  total = 0;
  do {
    total += acct->frac_bits[od_acct_index(state)];
  }
  while (od_acct_next_state(state, cat));
  return total;
}

void od_acct_print_state(FILE *_fp, unsigned int state[OD_ACCT_NCATS]) {
  int cat;
  for (cat = 0; cat < OD_ACCT_NCATS; cat++) {
    fprintf(_fp, "%s%s", cat > 0 ? "," : "",
     OD_ACCT_CATEGORY_VALUE_NAMES[cat][state[cat]]);
  }
}

void od_acct_print(od_acct *acct, FILE *_fp) {
  unsigned int state[OD_ACCT_NCATS];
  int i;
  for (i = 0; i < OD_ACCT_NCATS; i++) {
    state[i] = 0;
  }
  do {
    int index;
    index = od_acct_index(state);
    od_acct_print_state(_fp, state);
    fprintf(_fp, " (%i): %i\n", index, acct->frac_bits[index]);
  }
  while (od_acct_next_state(state, OD_ACCT_NCATS));
}

void od_acct_write(od_acct *acct, ogg_int64_t cur_time) {
  od_acct_category cat;
  unsigned int value;
  long fsize;
  OD_ASSERT(acct->fp);
  fsize = ftell(acct->fp);
  if (fsize == 0) {
    fprintf(acct->fp, "[");
  }
  else {
    fseek(acct->fp, fsize - 1, SEEK_SET);
    fprintf(acct->fp, ",\n");
  }
  fprintf(acct->fp, "{\n");
  fprintf(acct->fp, "  \"frame\": %" OD_I64FMT ",\n", (long long)cur_time);
  fprintf(acct->fp, "  \"total\": %u,\n", acct->last_frac_bits-8);
  for (cat = 0; cat < OD_ACCT_NCATS; cat++) {
    fprintf(acct->fp, "%s  \"%s\": {\n", cat > 0 ? ",\n" : "",
     OD_ACCT_CATEGORY_NAMES[cat]);
    for (value = 0; value < OD_ACCT_INDICES[cat]; value++) {
      fprintf(acct->fp, "%s    \"%s\": %i", value > 0 ? ",\n" : "",
       OD_ACCT_CATEGORY_VALUE_NAMES[cat][value],
       od_acct_get_total(acct, cat, value));
    }
    fprintf(acct->fp, "\n  }");
  }
  fprintf(acct->fp, "\n}]");
}
