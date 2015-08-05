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

/* Simple linear-time lookup in the dictionary. This could be made much
   faster, but it's fine for now. */
int od_accounting_dict_lookup(od_accounting_dict *dict, const char *str) {
  int i;
  for (i = 0; i < dict->nb_str; i++) {
    if (strcmp(dict->str[i], str) == 0) break;
  }
  OD_ASSERT(i < MAX_SYMBOL_TYPES);
  if (i == dict->nb_str) {
    dict->str[i] = malloc(strlen(str) + 1);
    strcpy(dict->str[i], str);
    dict->nb_str++;
  }
  return i;
}

void od_accounting_init(od_accounting_internal *acct) {
  acct->nb_syms_alloc = 1000;
  acct->acct.syms = malloc(sizeof(acct->acct.syms[0])*acct->nb_syms_alloc);
  acct->acct.dict.nb_str = 0;
  od_accounting_reset(acct);
}

void od_accounting_reset(od_accounting_internal *acct) {
  acct->acct.nb_syms = 0;
  acct->curr_x = acct->curr_y = acct->curr_level = acct->curr_layer = -1;
  acct->last_tell = 0;
}

void od_accounting_clear(od_accounting_internal *acct) {
  int i;
  free(acct->acct.syms);
  for (i = 0; i < acct->acct.dict.nb_str; i++) {
    free(acct->acct.dict.str[i]);
  }
}

void od_accounting_set_location(od_accounting_internal *acct, int layer,
 int level, int x, int y) {
  acct->curr_x = x;
  acct->curr_y = y;
  acct->curr_level = level;
  acct->curr_layer = layer;

}

void od_accounting_record(od_accounting_internal *acct, char *str,
 int bits_q3) {
  od_acct_symbol curr;
  int id;
  OD_ASSERT(acct->curr_x >= 0);
  OD_ASSERT(acct->curr_y >= 0);
  OD_ASSERT(bits_q3 <= 255);
  curr.x = acct->curr_x;
  curr.y = acct->curr_y;
  curr.level = acct->curr_level;
  curr.layer = acct->curr_layer;
  curr.bits_q3 = bits_q3;
  id = od_accounting_dict_lookup(&acct->acct.dict, str);
  OD_ASSERT(id <= 255);
  curr.id = id;
  if (acct->acct.nb_syms == acct->nb_syms_alloc) {
    acct->nb_syms_alloc *= 2;
    acct->acct.syms = realloc(acct->acct.syms,
     sizeof(acct->acct.syms[0])*acct->nb_syms_alloc);
    OD_ASSERT(acct->acct.syms != NULL);
  }
  acct->acct.syms[acct->acct.nb_syms++] = curr;
}
