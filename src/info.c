/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#include <stdlib.h>
#include <string.h>
#include "internal.h"

void daala_info_init(daala_info *_info){
  memset(_info,0,sizeof(*_info));
  _info->version_major=OD_VERSION_MAJOR;
  _info->version_minor=OD_VERSION_MINOR;
  _info->version_sub=OD_VERSION_SUB;
  _info->keyframe_granule_shift=31;
  /*TODO: Set other defaults.*/
}

void daala_info_clear(daala_info *_info){
  memset(_info,0,sizeof(*_info));
}

void daala_comment_init(daala_comment *_dc){
  memset(_dc,0,sizeof(*_dc));
}

void daala_comment_clear(daala_comment *_dc){
  if(_dc!=NULL){
    int ci;
    for(ci=0;ci<_dc->comments;ci++)_ogg_free(_dc->user_comments[ci]);
    _ogg_free(_dc->user_comments);
    _ogg_free(_dc->comment_lengths);
    _ogg_free(_dc->vendor);
    memset(_dc,0,sizeof(*_dc));
  }
}
