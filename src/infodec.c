/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

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

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "../include/daala/daaladec.h"
#include "decint.h"
#include "encint.h"

#define OD_INT32_MAX (0x7FFFFFFFUL)
#define OD_UINT32_MAX (0xFFFFFFFFUL)

static daala_setup_info *daala_setup_create() {
  return _ogg_malloc(sizeof(daala_setup_info));
}

void daala_setup_free(daala_setup_info *setup) {
  if (setup != NULL) {
    _ogg_free(setup);
  }
}

static int daala_read_length_and_string(oggbyte_buffer *obb, int *lenp,
 char **strp) {
  ogg_uint32_t len;
  char *str;
  if (oggbyte_read4(obb, &len)) return OD_EBADHEADER;
  /*Check that it is in bounds for where it is going, i.e., an int.*/
  if (len > INT_MAX) return OD_EBADHEADER;
  /*Check that enough bytes are left.*/
  if ((int)len > oggbyte_bytes_left(obb)) return OD_EBADHEADER;
  /*This is safe because oggbyte_bytes_left() is less than OD_MEM_DIFF_MAX.*/
  str = _ogg_malloc(len + 1);
  if (!str) return OD_EFAULT;
  str[len] = '\0';
  /*Since we made sure we had enough space above, this function should never
     fail.*/
  OD_ALWAYS_TRUE(!oggbyte_readcopy(obb, str, len));
  if (lenp) *lenp = (int)len;
  *strp = str;
  return 0;
}

static int od_comment_unpack(daala_comment *dc, oggbyte_buffer *obb) {
  ogg_uint32_t tmp;
  long comments;
  int rv;
  int i;
  daala_comment_init(dc);
  rv = daala_read_length_and_string(obb, NULL, &dc->vendor);
  if (rv) return rv;
  if (oggbyte_read4(obb, &tmp)) return OD_EBADHEADER;
  /*Check that this will fit into an integer, since that's where it is going.*/
  if (tmp > INT_MAX) return OD_EBADHEADER;
  /*Safe because LONG_MAX >= INT_MAX.*/
  comments = tmp;
  /*Check that we have enough room in the buffer for at least this many
     0-length comments.*/
  if ((comments > (LONG_MAX >> 2))
   || ((comments << 2) > oggbyte_bytes_left(obb))) {
    return OD_EBADHEADER;
  }
  dc->comments = comments;
  dc->comment_lengths =
   (int *)_ogg_malloc(sizeof(*dc->comment_lengths)*dc->comments);
  dc->user_comments =
   (char **)_ogg_malloc(sizeof(*dc->user_comments)*dc->comments);
  if (dc->comment_lengths == NULL || dc->user_comments == NULL) {
    dc->comments = 0;
    return OD_EFAULT;
  }
  for (i = 0; i < dc->comments; i++) {
    char *comment;
    int comment_len;
    if (daala_read_length_and_string(obb, &comment_len, &comment)) {
      dc->comments = i;
      return OD_EBADHEADER;
    }
    dc->comment_lengths[i] = comment_len;
    dc->user_comments[i] = comment;
  }
  return 0;
}

int daala_decode_header_in(daala_info *info,
 daala_comment *dc, daala_setup_info **ds, const ogg_packet *op) {
  oggbyte_buffer obb;
  int packtype;
  int rv;
  char daala[5];
  if (info == NULL || dc == NULL || ds == NULL) return OD_EFAULT;
  if (op == NULL) return OD_EBADHEADER;
  if (op->bytes < 0 || op->bytes > OD_MEM_DIFF_MAX) return OD_EFAULT;
  oggbyte_readinit(&obb, op->packet, op->bytes);
  packtype = oggbyte_read1(&obb);
  if (!(packtype & 0x80) && info->pic_width > 0
   && dc->vendor != NULL && *ds != NULL) {
    /*If we're at a data packet and we have received all three headers, we're
       done.*/
    return 0;
  }
  /*Check the codec string.*/
  rv = oggbyte_readcopy(&obb, daala, sizeof(daala));
  if (rv != 0) return OD_EBADHEADER;
  if (memcmp(daala, "daala", sizeof(daala)) != 0) return OD_EBADHEADER;
  switch (packtype) {
    /*Codec info header.*/
    case 0x80:
    {
      int pli;
      ogg_uint32_t tmp;
      int tmpi;
      /*This should be the first packet, and we should not have already read
         the info header packet yet.*/
      if (!op->b_o_s || info->pic_width) return OD_EBADHEADER;
      tmpi = oggbyte_read1(&obb);
      if (tmpi < 0) return OD_EBADHEADER;
      info->version_major = tmpi;
      tmpi = oggbyte_read1(&obb);
      if (tmpi < 0) return OD_EBADHEADER;
      info->version_minor = tmpi;
      tmpi = oggbyte_read1(&obb);
      if (tmpi < 0) return OD_EBADHEADER;
      info->version_sub = tmpi;
      if (info->version_major > OD_VERSION_MAJOR
       || (info->version_major == OD_VERSION_MAJOR
       && info->version_minor > OD_VERSION_MINOR)) {
        return OD_EVERSION;
      }
      if (oggbyte_read4(&obb, &tmp)) return OD_EBADHEADER;
      if (tmp > OD_INT32_MAX) return OD_EBADHEADER;
      info->pic_width = tmp;
      if (oggbyte_read4(&obb, &tmp))
        if (tmp > OD_INT32_MAX) return OD_EBADHEADER;
      info->pic_height = tmp;
      if (oggbyte_read4(&obb, &info->pixel_aspect_numerator)) {
        return OD_EBADHEADER;
      }
      if (oggbyte_read4(&obb, &info->pixel_aspect_denominator)) {
        return OD_EBADHEADER;
      }
      if (oggbyte_read4(&obb, &info->timebase_numerator)) {
        return OD_EBADHEADER;
      }
      if (oggbyte_read4(&obb, &info->timebase_denominator)) {
        return OD_EBADHEADER;
      }
      if (oggbyte_read4(&obb, &info->frame_duration)) return OD_EBADHEADER;
      tmpi = oggbyte_read1(&obb);
      if (tmpi < 0 || tmpi >= 32) return OD_EBADHEADER;
      info->keyframe_granule_shift = tmpi;
      info->nplanes = oggbyte_read1(&obb);
      if ((info->nplanes < 1) || (info->nplanes > OD_NPLANES_MAX)) {
        return OD_EBADHEADER;
      }
      for (pli = 0; pli < info->nplanes; pli++) {
        tmpi = oggbyte_read1(&obb);
        if (tmpi < 0) return OD_EBADHEADER;
        info->plane_info[pli].xdec = !!tmpi;
        tmpi = oggbyte_read1(&obb);
        if (tmpi < 0) return OD_EBADHEADER;
        info->plane_info[pli].ydec = !!tmpi;
      }
      return 3;
    }
    case 0x81:
    {
      /*Check that we have read the info header and have not read the
         comment header.*/
      if (!info->pic_width || dc->vendor != NULL) return OD_EBADHEADER;
      if (od_comment_unpack(dc, &obb)) {
        daala_comment_clear(dc);
        return OD_EBADHEADER;
      }
      return 2;
    }
    case 0x82:
    {
      /*Check that we have read the info header and the comment header,
         and not setup header.*/
      if (!info->pic_width || dc->vendor == NULL || *ds != NULL) {
        return OD_EBADHEADER;
      }
      *ds = daala_setup_create();
      return 1;
    }
  }
  return OD_EBADHEADER;
}
