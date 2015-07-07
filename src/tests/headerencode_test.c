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
#include "config.h"
#endif

#include "../encint.h"

#include <stdlib.h>
#include <check.h>

static daala_info di;
static daala_enc_ctx *dd;
static ogg_packet op;
static daala_comment dc;
static daala_setup_info *dsi;
static daala_info di2;
static daala_comment dc2;

static char *make_string(const char *str) {
  size_t len = strlen(str);
  char *ret = (char *)malloc(len + 1);
  memcpy(ret, str, len);
  ret[len] = '\0';
  return ret;
}

void setup(void) {
  /*Zero the struct and set the versions and keyframe_granule_shift.*/
  daala_info_init(&di);
  di.pic_width = 176;
  di.pic_height = 144;
  di.pixel_aspect_numerator = 128;
  di.pixel_aspect_denominator = 117;
  di.timebase_numerator = 30000;
  di.timebase_denominator = 1001;
  di.frame_duration = 30;
  di.nplanes = 2;
  dd = daala_encode_create(&di);
  ck_assert(dd != NULL);
  daala_comment_init(&dc);
  dc.comments = 2;
  dc.user_comments = (char **)malloc(sizeof(*dc.user_comments)*2);
  ck_assert(dc.user_comments != NULL);
  dc.comment_lengths = (int *)malloc(sizeof(*dc.comment_lengths)*2);
  ck_assert(dc.comment_lengths != NULL);
  dc.user_comments[0] = make_string("COMMENT=Comment 0");
  ck_assert(dc.user_comments[0] != NULL);
  dc.comment_lengths[0] = strlen(dc.user_comments[0]);
  dc.user_comments[1] = make_string("COMMENT=Comment 1 (this one longer)");
  ck_assert(dc.user_comments[1] != NULL);
  dc.comment_lengths[1] = strlen(dc.user_comments[1]);
}

void teardown(void) {
  daala_encode_free(dd);
}

START_TEST(encode_info_header) {
  int rv;
  daala_info_init(&di2);
  rv = daala_encode_flush_header(dd, &dc, &op);
  ck_assert_int_ne(OD_EFAULT, rv);
  rv = daala_decode_header_in(&di2, &dc2, &dsi, &op);
  /*TODO: ck_assert_int_le() gives me linker errors, so hack this for now.*/
  ck_assert_int_eq(1, rv >= 0);
  ck_assert_int_eq(0, memcmp(&di, &di2, sizeof(di2)));
}
END_TEST

START_TEST(encode_comment_header) {
  int rv;
  int i;
  daala_comment_init(&dc2);
  rv = daala_encode_flush_header(dd, &dc, &op);
  ck_assert_int_ne(OD_EFAULT, rv);
  rv = daala_decode_header_in(&di2, &dc2, &dsi, &op);
  /*TODO: ck_assert_int_le() gives me linker errors, so hack this for now.*/
  ck_assert_int_eq(1, rv >= 0);
  ck_assert_int_eq(0, strcmp(daala_version_string(), dc2.vendor));
  ck_assert_int_eq(dc.comments, dc2.comments);
  for (i = 0; i < dc.comments; i++) {
    ck_assert_int_eq(dc.comment_lengths[i], dc2.comment_lengths[i]);
    ck_assert_int_eq(0, strcmp(dc.user_comments[i], dc2.user_comments[i]));
  }
}
END_TEST

Suite *headerencode_suite() {
  Suite *s = suite_create("HeaderEncode");
  TCase *tc = tcase_create("HeaderEncode");
  tcase_add_unchecked_fixture (tc, setup, teardown);
  tcase_add_test(tc, encode_info_header);
  tcase_add_test(tc, encode_comment_header);
  suite_add_tcase(s, tc);
  return s;
}
