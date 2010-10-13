/*
    Daala video codec
    Copyright (C) 2010 Timothy B. Terriberry and Daala project contributors

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


#if !defined(_encint_H)
# define _encint_H (1)
# include "../include/daala/daaladec.h"
# include "state.h"



typedef struct daala_dec_ctx od_dec_ctx;



/*Constants for the packet state machine specific to the decoder.*/

/*Next packet to read: Data packet.*/
#define OD_PACKET_DATA        (0)



struct daala_dec_ctx{
  od_state        state;
  oggbyte_buffer  obb;
  int             packet_state;
};

#endif
