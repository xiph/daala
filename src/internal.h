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


#if !defined(_internal_H)
# define _internal_H (1)
# include <limits.h>
# include "../include/daala/codec.h"
# include "odintrin.h"


# define OD_VERSION_MAJOR (0)
# define OD_VERSION_MINOR (0)
# define OD_VERSION_SUB   (0)

# define OD_VENDOR_STRING "Xiph's experimental encoder library " __DATE__

/*Constants for the packet state machine common between encoder and decoder.*/

/*Next packet to emit/read: Codec info header.*/
#define OD_PACKET_INFO_HDR    (-3)
/*Next packet to emit/read: Comment header.*/
#define OD_PACKET_COMMENT_HDR (-2)
/*Next packet to emit/read: Codec setup header.*/
#define OD_PACKET_SETUP_HDR   (-1)
/*Next more packets to emit/read.*/
#define OD_PACKET_DONE        (INT_MAX)


# if defined(OD_ENABLE_ASSERTIONS)
#  include <stdio.h>
#  include <stdlib.h>
#  if __GNUC_PREREQ(2,5)
__attribute__((noreturn))
#  endif
void od_fatal_impl(const char *_str,const char *_file,int _line);

#  define OD_FATAL(_str) (od_fatal_impl(_str,__FILE__,__LINE__))

#  define OD_ASSERT(_cond) \
  do{ \
    if(!(_cond)){ \
      OD_FATAL("assertion failed: " #_cond); \
    } \
  } \
  while(0)

#  define OD_ASSERT2(_cond,_message) \
  do{ \
    if(!(_cond)){ \
      OD_FATAL("assertion failed: " #_cond "\n" _message); \
    } \
  } \
  while(0)

# else
#  define OD_ASSERT(_cond)
#  define OD_ASSERT2(_cond,_message)
# endif



#if 1
/*Currently this structure is only in Tremor, and is read-only.*/
typedef struct oggbyte_buffer oggbyte_buffer;

/*Simple libogg1-style buffer.*/
struct oggbyte_buffer{
  unsigned char *buf;
  unsigned char *ptr;
  long           storage;
};

/*Encoding functions.*/
void oggbyte_writeinit(oggbyte_buffer *_b);
void oggbyte_writetrunc(oggbyte_buffer *_b,long _bytes);
void oggbyte_write1(oggbyte_buffer *_b,unsigned _value);
void oggbyte_write4(oggbyte_buffer *_b,ogg_uint32_t _value);
void oggbyte_writecopy(oggbyte_buffer *_b,void *_source,long _bytes);
void oggbyte_writeclear(oggbyte_buffer *_b);
/*Decoding functions.*/
void oggbyte_readinit(oggbyte_buffer *_b,unsigned char *_buf,long _bytes);
int oggbyte_look1(oggbyte_buffer *_b);
int oggbyte_look4(oggbyte_buffer *_b,ogg_uint32_t *_val);
void oggbyte_adv1(oggbyte_buffer *_b);
void oggbyte_adv4(oggbyte_buffer *_b);
int oggbyte_read1(oggbyte_buffer *_b);
int oggbyte_read4(oggbyte_buffer *_b,ogg_uint32_t *_val);
/*Shared functions.*/
void oggbyte_reset(oggbyte_buffer *_b);
long oggbyte_bytes(oggbyte_buffer *_b);
unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b);

#endif

int od_ilog(ogg_uint32_t _v);
void **od_malloc_2d(size_t _height,size_t _width,size_t _sz);
void **od_calloc_2d(size_t _height,size_t _width,size_t _sz);
void od_free_2d(void *_ptr);

#endif
