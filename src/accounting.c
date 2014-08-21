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
#include <stdarg.h>
#include "logging.h"
#include "accounting.h"

/*Entropy accounting:
  The main function you will use is od_ec_acct_record(), which takes the
   symbol, the number of possible symbols, and a variable number of arguments
   for the context. Use this right before you call od_ec_encode_*() for the
   symbol.
  Before you call od_ec_acct_record(), you must create a label to hold the
   events. Use od_ec_acct_add_label() to add a new label and tell the system
   how many pieces of context you will use.
  For simply recording the symbol and having no context it would look like:
    od_ec_acct_add_label(&enc->ec.acct, "motion-flags-level-2", 0);
    ...
    od_ec_acct_record(&enc->ec.acct, "motion-flags-level-2", mvp->valid, 2);
    od_ec_encode_bool_q15(&enc->ec, mvp->valid, 16384);
  Here's an example record cal of the same thing but with 4 pieces of
   context:
   od_ec_acct_record(&enc->ec.acct, "motion-flags-level-4", mvp->valid, 2,
    vx > 0 && vy > 0 ? grid[vy - 1][vx - 1].valid : 0,
    vy > 0 ? grid[vy - 1][vx + 1].valid : 0,
    grid[vy - 1][vx].mv[0] == grid[vy][vx + 1].mv[0] &&
    grid[vy - 1][vx].mv[1] == grid[vy][vx + 1].mv[1],
    grid[vy + 1][vx].mv[0] == grid[vy][vx + 1].mv[0] &&
    grid[vy + 1][vx].mv[1] == grid[vy][vx + 1].mv[1]);
  By default these results will be written out per frame into
   ec-acct.json. You can override the filename with OD_EC_ACCT_SUFFIX in the
   environment.*/

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
  const char *pre;
  const char *suf;
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
  int cat;
  index = state[0];
  for (cat = (int)OD_ACCT_CAT_PLANE; cat < OD_ACCT_NCATS; cat++) {
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
 int skip) {
  int i;
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
 int cat, unsigned int value) {
  unsigned int state[OD_ACCT_NCATS];
  ogg_uint32_t total;
  int i;
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
  for (cat = OD_ACCT_CAT_TECHNIQUE; cat < OD_ACCT_NCATS; cat++) {
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
  int cat;
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

void od_ec_acct_init(od_ec_acct *acct) {
  char  fname[1024];
  const char *pre;
  const char *suf;
  int   rv;
  pre = "ec-acct-";
  suf = getenv("OD_EC_ACCT_SUFFIX");
  if (!suf) {
    pre = "ec-acct";
    suf = "";
  }
  rv = snprintf(fname, sizeof(fname), "%s%s.json", pre, suf);
  OD_ASSERT(rv >= 0 && ((size_t)rv) < sizeof(fname));
  acct->fp = fopen(fname, "w");
  OD_ASSERT(acct->fp);
  acct->data = NULL;
  od_ec_acct_reset(acct);
}

void od_ec_acct_clear(od_ec_acct *acct) {
  od_ec_acct_data *data;
  od_ec_acct_data *old;
  OD_ASSERT(acct->fp);
  fclose(acct->fp);
  data = acct->data;
  acct->data = NULL;
  while (data) {
    _ogg_free(data->values);
    old = data;
    data = data->next;
    _ogg_free(old);
  }
}

void od_ec_acct_reset(od_ec_acct *acct) {
  od_ec_acct_data *data;
  data = acct->data;
  while (data) {
    data->used = 0;
    data = data->next;
  }
}

void od_ec_acct_add_label(od_ec_acct *acct, const char *label, int ncontext) {
  od_ec_acct_data *data;
  od_ec_acct_data *old_data;
  int i;
  old_data = NULL;
  data = acct->data;
  while (data) {
    if (strcmp(label, data->label) == 0) {
      break;
    }
    old_data = data;
    data = data->next;
  }
  if (data == NULL) {
    data = (od_ec_acct_data *)_ogg_malloc(sizeof(od_ec_acct_data));
    OD_ASSERT(data);
    data->label = label;
    data->capacity = 128;
    data->used = 0;
    /*Records are composed of a symbol, the number of possible symbols, and
      ncontext items of context.*/
    data->reclen = ncontext + 2;
    data->values = (int **)_ogg_malloc(128 * sizeof(int *));
    for (i = 0; i < 128; i++) {
      data->values[i] = (int *)_ogg_malloc(data->reclen * sizeof(int));
    }
    OD_ASSERT(data->values);
    data->next = NULL;
    if (old_data == NULL) {
      acct->data = data;
    } else {
      old_data->next = data;
    }
  }
}

void od_ec_acct_record(od_ec_acct *acct, const char *label, int val, int n, ...) {
  va_list ap;
  od_ec_acct_data *data;
  int i;
  int old_capacity;
  data = acct->data;
  while (data) {
    if (strcmp(label, data->label) == 0) {
      break;
    }
    data = data->next;
  }
  OD_ASSERT(data);
  if (data->used >= data->capacity) {
    old_capacity = data->capacity;
    data->capacity *= 2;
    data->values = (int **)_ogg_realloc(data->values, data->capacity * sizeof(int *));
    for (i = old_capacity; i < data->capacity; i++) {
      data->values[i] = (int *)_ogg_malloc(data->reclen * sizeof(int));
    }
  }
  data->values[data->used][0] = val;
  data->values[data->used][1] = n;
  va_start(ap, n);
  for (i = 2; i < data->reclen; i++) {
    data->values[data->used][i] = va_arg(ap, int);
  }
  va_end(ap);
  data->used++;
}

void od_ec_acct_write(od_ec_acct *acct) {
  long fsize;
  od_ec_acct_data *data;
  int i;
  int j;
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
  data = acct->data;
  while (data) {
    fprintf(acct->fp, "  \"%s\": [\n", data->label);
    for (i = 0; i < data->used; i++) {
      fprintf(acct->fp, "    [");
      for (j = 0; j < data->reclen; j++) {
        fprintf(acct->fp, "%s%d", j > 0 ? "," : "", data->values[i][j]);
      }
      fprintf(acct->fp, "]%s\n", i == (data->used - 1) ? "" : ",");
    }
    fprintf(acct->fp, "  ]%s", data->next ? ",\n" : "");
    data = data->next;
  }
  fprintf(acct->fp, "\n}]");
}
