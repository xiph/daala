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


#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "internal.h"

int od_ilog(ogg_uint32_t _v){
#if defined(OD_CLZ)
  return OD_CLZ0-OD_CLZ(_v)&-!!_v;
#else
  /*On a Pentium M, this branchless version tested as the fastest on
     1,000,000,000 random 32-bit integers, edging out a similar version with
     branches, and a 256-entry LUT version.*/
  int ret;
  int m;
  ret=!!_v;
  m=!!(_v&0xFFFF0000)<<4;
  _v>>=m;
  ret|=m;
  m=!!(_v&0xFF00)<<3;
  _v>>=m;
  ret|=m;
  m=!!(_v&0xF0)<<2;
  _v>>=m;
  ret|=m;
  m=!!(_v&0xC)<<1;
  _v>>=m;
  ret|=m;
  ret+=!!(_v&0x2);
  return ret;
#endif
}

void **od_malloc_2d(size_t _height,size_t _width,size_t _sz){
  size_t  rowsz;
  size_t  colsz;
  size_t  datsz;
  char   *ret;
  colsz=_height*sizeof(void *);
  rowsz=_sz*_width;
  datsz=rowsz*_height;
  /*Alloc array and row pointers.*/
  ret=(char *)_ogg_malloc(datsz+colsz);
  /*Initialize the array.*/
  if(ret!=NULL){
    size_t   i;
    void   **p;
    char    *datptr;
    p=(void **)ret;
    i=_height;
    for(datptr=ret+colsz;i-->0;p++,datptr+=rowsz)*p=(void *)datptr;
  }
  return (void **)ret;
}

void **od_calloc_2d(size_t _height,size_t _width,size_t _sz){
  size_t  colsz;
  size_t  rowsz;
  size_t  datsz;
  char   *ret;
  colsz=_height*sizeof(void *);
  rowsz=_sz*_width;
  datsz=rowsz*_height;
  /*Alloc array and row pointers.*/
  ret=(char *)_ogg_calloc(datsz+colsz,1);
  /*Initialize the array.*/
  if(ret!=NULL){
    size_t   i;
    void   **p;
    char    *datptr;
    p=(void **)ret;
    i=_height;
    for(datptr=ret+colsz;i-->0;p++,datptr+=rowsz)*p=(void *)datptr;
  }
  return (void **)ret;
}

void od_free_2d(void *_ptr){
  _ogg_free(_ptr);
}



#define BUFFER_INCREMENT (256)

void oggbyte_writeinit(oggbyte_buffer *_b){
  memset(_b,0,sizeof(*_b));
  _b->ptr=_b->buf=_ogg_malloc(BUFFER_INCREMENT);
  _b->storage=BUFFER_INCREMENT;
}

void oggbyte_writetrunc(oggbyte_buffer *_b,long _bytes){
  _b->ptr=_b->buf+_bytes;
}

void oggbyte_write1(oggbyte_buffer *_b,unsigned _value){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte>=_b->storage){
    _b->buf=_ogg_realloc(_b->buf,_b->storage+BUFFER_INCREMENT);
    _b->storage+=BUFFER_INCREMENT;
    _b->ptr=_b->buf+endbyte;
  }
  *(_b->ptr++)=(unsigned char)_value;
}

void oggbyte_write4(oggbyte_buffer *_b,ogg_uint32_t _value){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte+4>_b->storage){
    _b->buf=_ogg_realloc(_b->buf,_b->storage+BUFFER_INCREMENT);
    _b->storage+=BUFFER_INCREMENT;
    _b->ptr=_b->buf+endbyte;
  }
  *(_b->ptr++)=(unsigned char)_value;
  _value>>=8;
  *(_b->ptr++)=(unsigned char)_value;
  _value>>=8;
  *(_b->ptr++)=(unsigned char)_value;
  _value>>=8;
  *(_b->ptr++)=(unsigned char)_value;
}

void oggbyte_writecopy(oggbyte_buffer *_b,void *_source,long _bytes){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte+_bytes>_b->storage){
    _b->storage=endbyte+_bytes+BUFFER_INCREMENT;
    _b->buf=_ogg_realloc(_b->buf,_b->storage);
    _b->ptr=_b->buf+endbyte;
  }
  memmove(_b->ptr,_source,_bytes);
  _b->ptr+=_bytes;
}

void oggbyte_reset(oggbyte_buffer *_b){
  _b->ptr=_b->buf;
}

void oggbyte_writeclear(oggbyte_buffer *_b){
  _ogg_free(_b->buf);
  memset(_b,0,sizeof(*_b));
}

void oggbyte_readinit(oggbyte_buffer *_b,unsigned char *_buf,long _bytes){
  memset(_b,0,sizeof(*_b));
  _b->buf=_b->ptr=_buf;
  _b->storage=_bytes;
}

int oggbyte_look1(oggbyte_buffer *_b){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte>=_b->storage)return -1;
  else return _b->ptr[0];
}

int oggbyte_look4(oggbyte_buffer *_b,ogg_uint32_t *_val){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte+4>_b->storage){
    if(endbyte<_b->storage){
      *_val=_b->ptr[0];
      endbyte++;
      if(endbyte<_b->storage){
        *_val|=(ogg_uint32_t)_b->ptr[1]<<8;
        endbyte++;
        if(endbyte<_b->storage)*_val|=(ogg_uint32_t)_b->ptr[2]<<16;
      }
    }
    return -1;
  }
  else{
    *_val=_b->ptr[0];
    *_val|=(ogg_uint32_t)_b->ptr[1]<<8;
    *_val|=(ogg_uint32_t)_b->ptr[2]<<16;
    *_val|=(ogg_uint32_t)_b->ptr[3]<<24;
  }
  return 0;
}

void oggbyte_adv1(oggbyte_buffer *_b){
  _b->ptr++;
}

void oggbyte_adv4(oggbyte_buffer *_b){
  _b->ptr+=4;
}

int oggbyte_read1(oggbyte_buffer *_b){
  ptrdiff_t endbyte;
  endbyte=_b->ptr-_b->buf;
  if(endbyte>=_b->storage)return -1;
  else return *(_b->ptr++);
}

int oggbyte_read4(oggbyte_buffer *_b,ogg_uint32_t *_val){
  unsigned char *end;
  end=_b->buf+_b->storage;
  if(_b->ptr+4>end){
    if(_b->ptr<end){
      *_val=*(_b->ptr++);
      if(_b->ptr<end){
        *_val|=(ogg_uint32_t)*(_b->ptr++)<<8;
        if(_b->ptr<end)*_val|=(ogg_uint32_t)*(_b->ptr++)<<16;
      }
    }
    return -1;
  }
  else{
    *_val=(*_b->ptr++);
    *_val|=(ogg_uint32_t)*(_b->ptr++)<<8;
    *_val|=(ogg_uint32_t)*(_b->ptr++)<<16;
    *_val|=(ogg_uint32_t)*(_b->ptr++)<<24;
  }
  return 0;
}

long oggbyte_bytes(oggbyte_buffer *_b){
  return _b->ptr-_b->buf;
}

unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b){
  return _b->buf;
}



const char *daala_version_string(void){
  return OD_VENDOR_STRING;
}

ogg_uint32_t daala_version_number(void){
  return OD_VERSION_MAJOR<<16|OD_VERSION_MINOR<<8|OD_VERSION_SUB;
}

int daala_packet_isheader(ogg_packet *_op){
  return _op->bytes>0?_op->packet[0]>>7:0;
}

int daala_packet_iskeyframe(ogg_packet *_op){
  return _op->bytes<=0?0:_op->packet[0]&0x80?-1:!(_op->packet[0]&0x40);
}
