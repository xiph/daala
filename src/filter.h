/*
    Daala video codec
    Copyright (C) 2003-2010 Daala project contributors

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


#if !defined(_filter_H)
# define _filter_H (1)
# include "internal.h"

/*These are the pre/post filtering functions used by Daala.
  The idea is to pre/post filter in the spatial domain (the time domain in
   signal processing terms) to improve the energy compaction as well as reduce
   or eliminate blocking artifacts.*/

void od_pre_filter_4(ogg_int16_t _y[4],const ogg_int16_t _x[4]);
void od_post_filter_4(ogg_int16_t _x[4],const ogg_int16_t _y[4]);
void od_pre_filter_8(ogg_int16_t _y[8],const ogg_int16_t _x[8]);
void od_post_filter_8(ogg_int16_t _x[8],const ogg_int16_t _y[8]);
void od_pre_filter_16(ogg_int16_t _y[16],const ogg_int16_t _x[16]);
void od_post_filter_16(ogg_int16_t _x[16],const ogg_int16_t _y[16]);

#endif
