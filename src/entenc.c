/* Copyright (c) 2001-2011 Timothy B. Terriberry
   Copyright (c) 2008-2009 Xiph.Org Foundation */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <stdlib.h>
#include <string.h>
#include "entenc.h"

/*A range encoder.
  See entdec.c and the references for implementation details \cite{Mar79,MNW98}.

  @INPROCEEDINGS{Mar79,
   author="Martin, G.N.N.",
   title="Range encoding: an algorithm for removing redundancy from a digitised
    message",
   booktitle="Video \& Data Recording Conference",
   year=1979,
   address="Southampton",
   month=Jul
  }
  @ARTICLE{MNW98,
   author="Alistair Moffat and Radford Neal and Ian H. Witten",
   title="Arithmetic Coding Revisited",
   journal="{ACM} Transactions on Information Systems",
   year=1998,
   volume=16,
   number=3,
   pages="256--294",
   month=Jul,
   URL="http://www.stanford.edu/class/ee398/handouts/papers/Moffat98ArithmCoding.pdf"
  }*/



static void od_ec_enc_normalize(od_ec_enc *_this,
 od_ec_window _low,unsigned _rng){
  int nbits_total;
  int d;
  int c;
  int s;
  nbits_total=_this->base.nbits_total;
  c=_this->base.cnt;
  d=16-OD_ILOG_NZ(_rng);
  s=c+d;
  /*TODO: Right now we flush every time we have at least one byte available.
    Instead we should use an ec_window and flush right before we're about to
     shift bits off the end of the window.
    For a 32-bit window this is about the same amount of work, but for a 64-bit
     window it should be a fair win.*/
  if(s>=0){
    ogg_uint16_t *buf;
    ogg_uint32_t  storage;
    ogg_uint32_t  offs;
    unsigned      m;
    buf=_this->precarry_buf;
    storage=_this->precarry_storage;
    offs=_this->base.offs;
    if(offs>=storage){
      storage=2*storage+2;
      buf=_ogg_realloc(buf,storage*sizeof(*buf));
      if(buf==NULL){
        _this->base.error=-1;
        _this->base.offs=0;
        return;
      }
      _this->precarry_buf=buf;
      _this->precarry_storage=storage;
    }
    c+=16;
    m=(1<<c)-1;
    if(s>=8){
      buf[offs++]=(ogg_uint16_t)(_low>>c);
      _low&=m;
      c-=8;
      m>>=8;
    }
    buf[offs++]=(ogg_uint16_t)(_low>>c);
    s=c+d-24;
    _low&=m;
    _this->base.offs=offs;
  }
  _this->base.nbits_total=nbits_total+d;
  _this->base.val=_low<<d;
  _this->base.rng=_rng<<d;
  _this->base.cnt=s;
}


void od_ec_enc_init(od_ec_enc *_this,ogg_uint32_t _size){
  od_ec_enc_reset(_this);
  _this->base.buf=(unsigned char *)_ogg_malloc(_size*sizeof(*_this->base.buf));
  _this->base.storage=_size;
  if(_size>0&&_this->base.buf==NULL){
    _this->base.storage=0;
    _this->base.error=-1;
  }
  _this->precarry_buf=
   (ogg_uint16_t *)_ogg_malloc(_size*sizeof(*_this->precarry_buf));
  _this->precarry_storage=_size;
  if(_size>0&&_this->precarry_buf==NULL){
    _this->precarry_storage=0;
    _this->base.error=-1;
  }
}

void od_ec_enc_clear(od_ec_enc *_this){
  _ogg_free(_this->precarry_buf);
  _ogg_free(_this->base.buf);
}


void od_ec_enc_reset(od_ec_enc *_this){
  _this->base.end_offs=0;
  _this->base.end_window=0;
  _this->base.nend_bits=0;
  /*We reserve one bit for termination.*/
  _this->base.nbits_total=1;
  _this->base.offs=0;
  _this->base.rng=0x8000;
  _this->base.cnt=-9;
  _this->base.val=0;
  _this->base.error=0;
}

void od_ec_encode_normalized(od_ec_enc *_this,
 unsigned _fl,unsigned _fh,unsigned _ft){
  od_ec_window l;
  unsigned     r;
  int          s;
  unsigned     d;
  unsigned     u;
  unsigned     v;
  l=_this->base.val;
  r=_this->base.rng;
  s=r>=_ft<<1;
  _ft<<=s;
  _fl<<=s;
  _fh<<=s;
  d=r-_ft;
  u=_fl+OD_MINI(_fl,d);
  v=_fh+OD_MINI(_fh,d);
  r=v-u;
  l+=u;
  od_ec_enc_normalize(_this,l,r);
}

void od_ec_encode(od_ec_enc *_this,unsigned _fl,unsigned _fh,unsigned _ft){
  int s;
  OD_ASSERT(_ft>=2);
  s=15-OD_ILOG_NZ(_ft-1);
  od_ec_encode_normalized(_this,_fl<<s,_fh<<s,_ft<<s);
}

void od_ec_encode_bin_normalized(od_ec_enc *_this,unsigned _fl,unsigned _fh){
  od_ec_window l;
  unsigned     r;
  unsigned     d;
  unsigned     u;
  unsigned     v;
  l=_this->base.val;
  r=_this->base.rng;
  d=r-32768U;
  u=_fl+OD_MINI(_fl,d);
  v=_fh+OD_MINI(_fh,d);
  r=v-u;
  l+=u;
  od_ec_enc_normalize(_this,l,r);
}

void od_ec_encode_bin(od_ec_enc *_this,
 unsigned _fl,unsigned _fh,unsigned _ftb){
  od_ec_encode_bin_normalized(_this,_fl<<15-_ftb,_fh<<15-_ftb);
}

/*The probability of having a "one" is 1/(1<<_logp).*/
void od_ec_enc_bit_logp(od_ec_enc *_this,int _val,unsigned _logp){
  od_ec_window l;
  unsigned     r;
  unsigned     v;
  l=_this->base.val;
  r=_this->base.rng;
  v=32768U-(1<<15-_logp);
  v+=OD_MINI(v,r-32768U);
  if(_val)l+=v;
  r=_val?r-v:v;
  od_ec_enc_normalize(_this,l,r);
}

void od_ec_enc_icdf(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ftb){
  if(_s>0){
    od_ec_encode_bin_normalized(_this,32768U-(_icdf[_s-1]<<15-_ftb),
     32768U-(_icdf[_s]<<15-_ftb));
  }
  else{
    od_ec_window l;
    unsigned     r;
    unsigned     d;
    l=_this->base.val;
    r=_this->base.rng;
    d=32768U-(_icdf[0]<<15-_ftb);
    r=d+OD_MINI(d,r-32768U);
    od_ec_enc_normalize(_this,l,r);
  }
}

void od_ec_enc_icdf_ft(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ft){
  od_ec_encode(_this,_s>0?_ft-_icdf[_s-1]:0,_ft-_icdf[_s],_ft);
}

void od_ec_enc_icdf16(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_icdf,unsigned _ftb){
  if(_s>0){
    od_ec_encode_bin_normalized(_this,32768U-(_icdf[_s-1]<<15-_ftb),
     32768U-(_icdf[_s]<<15-_ftb));
  }
  else{
    od_ec_window l;
    unsigned     r;
    unsigned     d;
    l=_this->base.val;
    r=_this->base.rng;
    d=32768U-(_icdf[0]<<15-_ftb);
    r=d+OD_MINI(d,r-32768U);
    od_ec_enc_normalize(_this,l,r);
  }
}

void od_ec_enc_icdf16_ft(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_icdf,unsigned _ft){
  od_ec_encode(_this,_s>0?_ft-_icdf[_s-1]:0,_ft-_icdf[_s],_ft);
}

void od_ec_enc_uint(od_ec_enc *_this,ogg_uint32_t _fl,ogg_uint32_t _ft){
  if(_ft>1U<<OD_EC_UINT_BITS){
    unsigned ft;
    unsigned fl;
    int      ftb;
    _ft--;
    ftb=OD_ILOG_NZ(_ft)-OD_EC_UINT_BITS;
    ft=(_ft>>ftb)+1<<15-OD_EC_UINT_BITS;
    fl=(unsigned)(_fl>>ftb)<<15-OD_EC_UINT_BITS;
    od_ec_encode_normalized(_this,fl,fl+(1<<15-OD_EC_UINT_BITS),ft);
    od_ec_enc_bits(_this,_fl&((ogg_uint32_t)1<<ftb)-1U,ftb);
  }
  else od_ec_encode(_this,_fl,_fl+1,_ft);
}

void od_ec_enc_bits(od_ec_enc *_this,ogg_uint32_t _fl,unsigned _bits){
  od_ec_window end_window;
  int          nend_bits;
  end_window=_this->base.end_window;
  nend_bits=_this->base.nend_bits;
  OD_ASSERT(_bits>0);
  if(nend_bits+_bits>OD_EC_WINDOW_SIZE){
    unsigned char *buf;
    ogg_uint32_t   storage;
    ogg_uint32_t   end_offs;
    buf=_this->base.buf;
    storage=_this->base.storage;
    end_offs=_this->base.end_offs;
    if(end_offs+(OD_EC_WINDOW_SIZE>>3)>=storage){
      unsigned char *new_buf;
      ogg_uint32_t   new_storage;
      new_storage=2*storage+(OD_EC_WINDOW_SIZE>>3);
      new_buf=(unsigned char *)_ogg_malloc(new_storage*sizeof(*new_buf));
      if(new_buf==NULL){
        _this->base.error=-1;
        _this->base.end_offs=0;
        return;
      }
      memcpy(new_buf+new_storage-end_offs,buf+storage-end_offs,
       end_offs*sizeof(*new_buf));
      _ogg_free(buf);
      _this->base.buf=new_buf;
      _this->base.storage=new_storage;
    }
    do{
      buf[storage-++(end_offs)]=(unsigned char)end_window;
      end_window>>=8;
      nend_bits-=8;
    }
    while(nend_bits>=8);
    _this->base.end_offs=end_offs;
  }
  end_window|=(od_ec_window)_fl<<nend_bits;
  nend_bits+=_bits;
  _this->base.end_window=end_window;
  _this->base.nend_bits=nend_bits;
  _this->base.nbits_total+=_bits;
}

void od_ec_enc_patch_initial_bits(od_ec_enc *_this,
 unsigned _val,unsigned _nbits){
  int      shift;
  unsigned mask;
  OD_ASSERT(_nbits<=8);
  shift=8-_nbits;
  mask=(1<<_nbits)-1<<shift;
  if(_this->base.offs>0){
    /*The first byte has been finalized.*/
    _this->precarry_buf[0]=
     (ogg_uint16_t)(_this->precarry_buf[0]&~mask|_val<<shift);
  }
  else if(9+_this->base.cnt+(_this->base.rng==0x8000)>_nbits){
    /*The first byte has yet to be output.*/
    _this->base.val=(_this->base.val&~((od_ec_window)mask<<16+_this->base.cnt))|
     (od_ec_window)_val<<16+_this->base.cnt+shift;
  }
  /*The encoder hasn't even encoded _nbits of data yet.*/
  else _this->base.error=-1;
}

unsigned char *od_ec_enc_done(od_ec_enc *_this,ogg_uint32_t *_nbytes){
  unsigned char *out;
  ogg_uint32_t   storage;
  ogg_uint16_t  *buf;
  ogg_uint32_t   offs;
  ogg_uint32_t   end_offs;
  int            nend_bits;
  od_ec_window   m;
  od_ec_window   e;
  od_ec_window   l;
  unsigned       r;
  int            c;
  int            s;
  if(_this->base.error)return NULL;
  /*We output the minimum number of bits that ensures that the symbols encoded
     thus far will be decoded correctly regardless of the bits that follow.*/
  l=_this->base.val;
  r=_this->base.rng;
  c=_this->base.cnt;
  s=9;
  m=0x7FFF;
  e=l+m&~m;
  while((e|m)>=l+r){
    s++;
    m>>=1;
    e=l+m&~m;
  }
  s+=c;
  offs=_this->base.offs;
  buf=_this->precarry_buf;
  if(s>0){
    unsigned m;
    storage=_this->precarry_storage;
    if(offs>=storage){
      storage=storage*2+2;
      buf=(ogg_uint16_t *)_ogg_realloc(buf,storage*sizeof(*buf));
      if(buf==NULL){
        _this->base.error=-1;
        return NULL;
      }
      _this->precarry_buf=buf;
      _this->precarry_storage=storage;
    }
    m=(1<<c+16)-1;
    do{
      buf[offs++]=(ogg_uint16_t)(e>>c+16);
      e&=m;
      s-=8;
      c-=8;
      m>>=8;
    }
    while(s>0);
  }
  /*Make sure there's enough room for the range-coder bits and the raw bits.*/
  out=_this->base.buf;
  storage=_this->base.storage;
  end_offs=_this->base.end_offs;
  e=_this->base.end_window;
  nend_bits=_this->base.nend_bits;
  s=-s;
  c=OD_MAXI(nend_bits-s>>3,0);
  if(offs+end_offs+c>storage){
    storage=offs+end_offs+c;
    out=(unsigned char *)_ogg_realloc(out,storage*sizeof(*out));
    if(out==NULL){
      _this->base.error=-1;
      return NULL;
    }
    _this->base.buf=out;
    _this->base.storage=storage;
  }
  /*If we have buffered raw bits, flush them as well.*/
  while(nend_bits>s){
    out[storage-++(end_offs)]=(unsigned char)e;
    e>>=8;
    nend_bits-=8;
  }
  *_nbytes=offs+end_offs;
  /*Perform carry propagation.*/
  out=out+storage-(offs+end_offs);
  c=0;
  end_offs=offs;
  while(offs-->0){
    c=buf[offs]+c;
    out[offs]=(unsigned char)c;
    c>>=8;
  }
  /*Add any remaining raw bits to the last byte.
    There is guaranteed to be enough room, because nend_bits<=s.*/
  if(nend_bits>0)out[end_offs-1]|=(unsigned char)e;
  return out;
}
