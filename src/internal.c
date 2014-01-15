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

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "internal.h"

/*Constants for use with OD_DIVU_SMALL().
  See \cite{Rob05} for details on computing these constants.
  @INPROCEEDINGS{Rob05,
    author="Arch D. Robison",
    title="{N}-bit Unsigned Division via {N}-bit Multiply-Add",
    booktitle="Proc. of the 17th IEEE Symposium on Computer Arithmetic
     (ARITH'05)",
    pages="131--139",
    address="Cape Cod, MA",
    month=Jun,
    year=2005
  }*/
ogg_uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2] = {
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xAAAAAAAB,         0 },
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xCCCCCCCD,         0 },
  { 0xAAAAAAAB,         0 },
  { 0x92492492, 0x92492492 },
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xE38E38E4,         0 },
  { 0xCCCCCCCD,         0 },
  { 0xBA2E8BA3,         0 },
  { 0xAAAAAAAB,         0 },
  { 0x9D89D89E,         0 },
  { 0x92492492, 0x92492492 },
  { 0x88888889,         0 },
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xF0F0F0F1,         0 },
  { 0xE38E38E4,         0 },
  { 0xD79435E5, 0xD79435E5 },
  { 0xCCCCCCCD,         0 },
  { 0xC30C30C3, 0xC30C30C3 },
  { 0xBA2E8BA3,         0 },
  { 0xB21642C9,         0 },
  { 0xAAAAAAAB,         0 },
  { 0xA3D70A3E,         0 },
  { 0x9D89D89E,         0 },
  { 0x97B425ED, 0x97B425ED },
  { 0x92492492, 0x92492492 },
  { 0x8D3DCB09,         0 },
  { 0x88888889,         0 },
  { 0x84210842, 0x84210842 },
  { 0xFFFFFFFF, 0xFFFFFFFF },
  { 0xF83E0F84,         0 },
  { 0xF0F0F0F1,         0 },
  { 0xEA0EA0EA, 0xEA0EA0EA },
  { 0xE38E38E4,         0 },
  { 0xDD67C8A6, 0xDD67C8A6 },
  { 0xD79435E5, 0xD79435E5 },
  { 0xD20D20D2, 0xD20D20D2 },
  { 0xCCCCCCCD,         0 },
  { 0xC7CE0C7D,         0 },
  { 0xC30C30C3, 0xC30C30C3 },
  { 0xBE82FA0C,         0 },
  { 0xBA2E8BA3,         0 },
  { 0xB60B60B6, 0xB60B60B6 },
  { 0xB21642C9,         0 },
  { 0xAE4C415D,         0 },
  { 0xAAAAAAAB,         0 },
  { 0xA72F053A,         0 },
  { 0xA3D70A3E,         0 },
  { 0xA0A0A0A1,         0 },
  { 0x9D89D89E,         0 },
  { 0x9A90E7D9, 0x9A90E7D9 },
  { 0x97B425ED, 0x97B425ED },
  { 0x94F2094F, 0x94F2094F },
  { 0x92492492, 0x92492492 },
  { 0x8FB823EE, 0x8FB823EE },
  { 0x8D3DCB09,         0 },
  { 0x8AD8F2FC,         0 },
  { 0x88888889,         0 },
  { 0x864B8A7E,         0 },
  { 0x84210842, 0x84210842 },
  { 0x82082082, 0x82082082 },
  { 0xFFFFFFFF, 0xFFFFFFFF }
};

#if defined(OD_ENABLE_ASSERTIONS)
void od_fatal_impl(const char *_str, const char *_file, int _line) {
  fprintf(stderr, "Fatal (internal) error in %s, line %d: %s\n",
   _file, _line, _str);
  abort();
}
#endif

int od_ilog(ogg_uint32_t _v) {
#if defined(OD_CLZ)
  return (OD_CLZ0-OD_CLZ(_v))&-!!_v;
#else
  /*On a Pentium M, this branchless version tested as the fastest on
     1,000,000,000 random 32-bit integers, edging out a similar version with
     branches, and a 256-entry LUT version.*/
  int ret;
  int m;
  ret = !!_v;
  m = !!(_v&0xFFFF0000)<<4;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xFF00)<<3;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xF0)<<2;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xC)<<1;
  _v >>= m;
  ret |= m;
  ret += !!(_v&0x2);
  return ret;
#endif
}

void **od_malloc_2d(size_t _height, size_t _width, size_t _sz) {
  size_t  rowsz;
  size_t  colsz;
  size_t  datsz;
  char *ret;
  colsz = _height*sizeof(void *);
  rowsz = _sz*_width;
  datsz = rowsz*_height;
  /*Alloc array and row pointers.*/
  ret = (char *)_ogg_malloc(datsz+colsz);
  /*Initialize the array.*/
  if (ret != NULL) {
    size_t   i;
    void **p;
    char *datptr;
    p = (void **)ret;
    i = _height;
    for (datptr = ret + colsz; i-- > 0; p++, datptr += rowsz)
      *p = (void *)datptr;
  }
  return (void **)ret;
}

void **od_calloc_2d(size_t _height, size_t _width, size_t _sz) {
  size_t  colsz;
  size_t  rowsz;
  size_t  datsz;
  char *ret;
  colsz = _height*sizeof(void *);
  rowsz = _sz*_width;
  datsz = rowsz*_height;
  /*Alloc array and row pointers.*/
  ret = (char *)_ogg_calloc(datsz + colsz, 1);
  /*Initialize the array.*/
  if (ret != NULL) {
    size_t   i;
    void **p;
    char *datptr;
    p = (void **)ret;
    i = _height;
    for (datptr = ret + colsz; i-- > 0; p++, datptr += rowsz)
      *p = (void *)datptr;
  }
  return (void **)ret;
}

void od_free_2d(void *_ptr) {
  _ogg_free(_ptr);
}



#define BUFFER_INCREMENT (256)

void oggbyte_writeinit(oggbyte_buffer *_b) {
  OD_CLEAR(_b, 1);
  _b->ptr = _b->buf = _ogg_malloc(BUFFER_INCREMENT);
  _b->storage = BUFFER_INCREMENT;
}

void oggbyte_writetrunc(oggbyte_buffer *_b, ptrdiff_t _bytes) {
  OD_ASSERT(_bytes >= 0);
  _b->ptr = _b->buf + _bytes;
}

void oggbyte_write1(oggbyte_buffer *_b, unsigned _value) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte >= _b->storage) {
    _b->buf = _ogg_realloc(_b->buf, _b->storage + BUFFER_INCREMENT);
    _b->storage += BUFFER_INCREMENT;
    _b->ptr = _b->buf+endbyte;
  }
  *(_b->ptr++) = (unsigned char)_value;
}

void oggbyte_write4(oggbyte_buffer *_b, ogg_uint32_t _value) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  if (endbyte+4 > _b->storage) {
    _b->buf = _ogg_realloc(_b->buf, _b->storage + BUFFER_INCREMENT);
    _b->storage += BUFFER_INCREMENT;
    _b->ptr = _b->buf + endbyte;
  }
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
}

void oggbyte_writecopy(oggbyte_buffer *_b, const void *_source,
 ptrdiff_t _bytes) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte+_bytes > _b->storage) {
    _b->storage = endbyte+_bytes+BUFFER_INCREMENT;
    _b->buf = _ogg_realloc(_b->buf, _b->storage);
    _b->ptr = _b->buf+endbyte;
  }
  memmove(_b->ptr, _source, _bytes);
  _b->ptr += _bytes;
}

void oggbyte_reset(oggbyte_buffer *_b) {
  _b->ptr = _b->buf;
}

void oggbyte_writeclear(oggbyte_buffer *_b) {
  _ogg_free(_b->buf);
  OD_CLEAR(_b, 1);
}

void oggbyte_readinit(oggbyte_buffer *_b, unsigned char *_buf,
 ptrdiff_t _bytes) {
  OD_ASSERT(_bytes >= 0);
  OD_CLEAR(_b, 1);
  _b->buf = _b->ptr = _buf;
  _b->storage = _bytes;
}

int oggbyte_look1(oggbyte_buffer *_b) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  if (endbyte >= _b->storage) return -1;
  else return _b->ptr[0];
}

int oggbyte_look4(oggbyte_buffer *_b, ogg_uint32_t *_val) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte > _b->storage-4) {
    if (endbyte < _b->storage) {
      *_val = _b->ptr[0];
      endbyte++;
      if (endbyte < _b->storage) {
        *_val |= (ogg_uint32_t)_b->ptr[1]<<8;
        endbyte++;
        if (endbyte < _b->storage) *_val |= (ogg_uint32_t)_b->ptr[2]<<16;
      }
    }
    return -1;
  }
  else {
    *_val = _b->ptr[0];
    *_val |= (ogg_uint32_t)_b->ptr[1]<<8;
    *_val |= (ogg_uint32_t)_b->ptr[2]<<16;
    *_val |= (ogg_uint32_t)_b->ptr[3]<<24;
  }
  return 0;
}

void oggbyte_adv1(oggbyte_buffer *_b) {
  _b->ptr++;
}

void oggbyte_adv4(oggbyte_buffer *_b) {
  _b->ptr += 4;
}

int oggbyte_read1(oggbyte_buffer *_b) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte >= _b->storage) return -1;
  else return *(_b->ptr++);
}

int oggbyte_read4(oggbyte_buffer *_b, ogg_uint32_t *_val) {
  unsigned char *end;
  end = _b->buf+_b->storage;
  if (_b->ptr+4 > end) {
    if (_b->ptr < end) {
      *_val = *(_b->ptr++);
      if (_b->ptr < end) {
        *_val |= (ogg_uint32_t)*(_b->ptr++)<<8;
        if (_b->ptr < end) *_val |= (ogg_uint32_t)*(_b->ptr++)<<16;
      }
    }
    return -1;
  }
  else {
    *_val = (*_b->ptr++);
    *_val |= (ogg_uint32_t)*(_b->ptr++)<<8;
    *_val |= (ogg_uint32_t)*(_b->ptr++)<<16;
    *_val |= (ogg_uint32_t)*(_b->ptr++)<<24;
  }
  return 0;
}

int oggbyte_readcopy(oggbyte_buffer *_b, void *_dest, ogg_uint32_t _bytes) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  OD_ASSERT(endbyte >= 0);
  OD_ASSERT(endbyte <= _b->storage);
  if ((size_t)(_b->storage - endbyte) < _bytes) return -1;
  memcpy(_dest, _b->ptr, _bytes);
  _b->ptr += _bytes;
  return 0;
}

ptrdiff_t oggbyte_bytes(oggbyte_buffer *_b) {
  return _b->ptr-_b->buf;
}

ptrdiff_t oggbyte_bytes_left(oggbyte_buffer *_b) {
  return _b->storage - oggbyte_bytes(_b);
}

unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b) {
  return _b->buf;
}



const char *daala_version_string(void) {
  return OD_VENDOR_STRING;
}

ogg_uint32_t daala_version_number(void) {
  return OD_VERSION_MAJOR<<16|OD_VERSION_MINOR<<8|OD_VERSION_SUB;
}

int daala_packet_isheader(ogg_packet *_op) {
  return _op->bytes > 0 ? _op->packet[0]>>7 : 0;
}

int daala_packet_iskeyframe(ogg_packet *_op) {
  return _op->bytes <= 0 ? 0 : _op->packet[0]&0x80 ?
   -1 : !(_op->packet[0]&0x40);
}
