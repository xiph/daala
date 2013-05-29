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

#include "../encint.h"

#include <stdlib.h>
#include <check.h>

static daala_info di;
static daala_enc_ctx *dd;
static ogg_packet op;

void setup(void) {
  daala_info_init(&di);
  di.nplanes = 2;
  dd = daala_encode_create(&di);
  ck_assert(dd != NULL);
}

void teardown(void) {
  daala_encode_free(dd);
}

START_TEST(encode_packet_info) {
  int rv;

  rv = daala_encode_flush_header(dd, NULL, &op);
  ck_assert_int_ne(OD_EFAULT, rv);
}
END_TEST


Suite *headerencode_suite() {
  Suite *s = suite_create("HeaderEncode");
  TCase *tc = tcase_create("HeaderEncode");
  tcase_add_checked_fixture (tc, setup, teardown);
  tcase_add_test(tc, encode_packet_info);
  suite_add_tcase(s, tc);

  return s;
}
