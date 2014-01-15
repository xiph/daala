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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "internal.h"

void daala_info_init(daala_info *_info) {
  OD_CLEAR(_info, 1);
  _info->version_major = OD_VERSION_MAJOR;
  _info->version_minor = OD_VERSION_MINOR;
  _info->version_sub = OD_VERSION_SUB;
  _info->keyframe_granule_shift = 31;
  /*TODO: Set other defaults.*/
}

void daala_info_clear(daala_info *_info) {
  OD_CLEAR(_info, 1);
}

void daala_comment_init(daala_comment *_dc) {
  OD_CLEAR(_dc, 1);
}

void daala_comment_clear(daala_comment *_dc) {
  if (_dc != NULL) {
    int ci;
    for (ci = 0; ci < _dc->comments; ci++) _ogg_free(_dc->user_comments[ci]);
    _ogg_free(_dc->user_comments);
    _ogg_free(_dc->comment_lengths);
    _ogg_free(_dc->vendor);
    OD_CLEAR(_dc, 1);
  }
}
