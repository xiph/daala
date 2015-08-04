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

static int od_acct_hash(const char *str) {
  uint32_t val;
  const unsigned char *ustr;
  val = 0;
  ustr = (const unsigned char*)str;
  /* This is about the worst hash one can design, but it should be good enough
     here. */
  while (*ustr) val += *ustr++;
  return val % OD_ACCT_HASH_SIZE;
}

/* Dictionary lookup based on an open-addressing hash table. */
int od_accounting_dict_lookup(od_accounting_internal *acct, const char *str) {
  int hash;
  od_accounting_dict *dict;
  dict = &acct->acct.dict;
  hash = od_acct_hash(str);
  while (acct->hash_dict[hash] != -1) {
    if (strcmp(dict->str[acct->hash_dict[hash]], str) == 0) {
      return acct->hash_dict[hash];
    }
    hash++;
    if (hash == OD_ACCT_HASH_SIZE) hash = 0;
  }
  /* No match found. */
  OD_ASSERT(dict->nb_str + 1 < MAX_SYMBOL_TYPES);
  acct->hash_dict[hash] = dict->nb_str;
  dict->str[dict->nb_str] = malloc(strlen(str) + 1);
  strcpy(dict->str[dict->nb_str], str);
  dict->nb_str++;
  return dict->nb_str - 1;
}

void od_accounting_init(od_accounting_internal *acct) {
  int i;
  acct->nb_syms_alloc = 1000;
  acct->acct.syms = malloc(sizeof(acct->acct.syms[0])*acct->nb_syms_alloc);
  acct->acct.dict.nb_str = 0;
  OD_ASSERT(OD_ACCT_HASH_SIZE > 2*MAX_SYMBOL_TYPES);
  for (i = 0; i < OD_ACCT_HASH_SIZE; i++) acct->hash_dict[i] = -1;
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
  id = od_accounting_dict_lookup(acct, str);
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
