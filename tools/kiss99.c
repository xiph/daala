/*Daala video codec.
  Copyright (C) 2002-2012 Daala project contributors.

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
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.*/
#include "kiss99.h"



void kiss99_srand(kiss99_ctx *_this,const unsigned char *_data,int _ndata){
  int i;
  _this->z=362436069;
  _this->w=521288629;
  _this->jsr=123456789;
  _this->jcong=380116160;
  for(i=3;i<_ndata;i+=4){
    _this->z^=_data[i-3];
    _this->w^=_data[i-2];
    _this->jsr^=_data[i-1];
    _this->jcong^=_data[i];
    kiss99_rand(_this);
  }
  if(i-3<_ndata)_this->z^=_data[i-3];
  if(i-2<_ndata)_this->w^=_data[i-2];
  if(i-1<_ndata)_this->jsr^=_data[i-1];
}

uint32_t kiss99_rand(kiss99_ctx *_this){
  uint32_t znew;
  uint32_t wnew;
  uint32_t mwc;
  uint32_t shr3;
  uint32_t cong;
  znew=36969*(_this->z&0xFFFF)+(_this->z>>16);
  wnew=18000*(_this->w&0xFFFF)+(_this->w>>16);
  mwc=(znew<<16)+wnew;
  shr3=_this->jsr^(_this->jsr<<17);
  shr3^=shr3>>13;
  shr3^=shr3<<5;
  cong=69069*_this->jcong+1234567;
  _this->z=znew;
  _this->w=wnew;
  _this->jsr=shr3;
  _this->jcong=cong;
  return (mwc^cong)+shr3;
}
