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
#if !defined(_kiss99_H)
# define _kiss99_H (1)
# include <stdint.h>



typedef struct kiss99_ctx kiss99_ctx;



struct kiss99_ctx{
  uint32_t z;
  uint32_t w;
  uint32_t jsr;
  uint32_t jcong;
};


void kiss99_srand(kiss99_ctx *_this,const unsigned char *_data,int _ndata);
uint32_t kiss99_rand(kiss99_ctx *_this);

#endif
