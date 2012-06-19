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

static int od_ec_write_byte(od_ec_enc *_this,unsigned _value){
  if(_this->offs+_this->end_offs>=_this->storage)return -1;
  _this->buf[_this->offs++]=(unsigned char)_value;
  return 0;
}

static int od_ec_write_byte_at_end(od_ec_enc *_this,unsigned _value){
  if(_this->offs+_this->end_offs>=_this->storage)return -1;
  _this->buf[_this->storage-++(_this->end_offs)]=(unsigned char)_value;
  return 0;
}

/*Outputs a symbol, with a carry bit.
  If there is a potential to propagate a carry over several symbols, they are
   buffered until it can be determined whether or not an actual carry will
   occur.
  If the counter for the buffered symbols overflows, then the stream becomes
   undecodable.
  This gives a theoretical limit of a few billion symbols in a single packet on
   32-bit systems.
  The alternative is to truncate the range in order to force a carry, but
   requires similar carry tracking in the decoder, needlessly slowing it down.*/
static void od_ec_enc_carry_out(od_ec_enc *_this,int _c){
  if(_c!=0xFF){
    /*No further carry propagation possible, flush buffer.*/
    int carry;
    carry=_c>>8;
    /*Don't output a byte on the first write.
      This compare should be taken care of by branch-prediction thereafter.*/
    if(_this->rem>=0)_this->error|=od_ec_write_byte(_this,_this->rem+carry);
    if(_this->ext>0){
      unsigned sym;
      sym=0xFF+carry&0xFF;
      do _this->error|=od_ec_write_byte(_this,sym);
      while(--(_this->ext)>0);
    }
    _this->rem=_c&0xFF;
  }
  else _this->ext++;
}

static void od_ec_enc_normalize(od_ec_enc *_this){
  /*If the range is too small, output some bits and rescale it.*/
  while(_this->rng<=0x800000){
    od_ec_enc_carry_out(_this,(int)(_this->val>>23));
    /*Move the next-to-high-order symbol into the high-order position.*/
    _this->val=_this->val<<8&0x7FFFFFFF;
    _this->rng<<=8;
    _this->nbits_total+=8;
  }
}

void od_ec_enc_init(od_ec_enc *_this,unsigned char *_buf,ogg_uint32_t _size){
  _this->buf=_buf;
  _this->end_offs=0;
  _this->end_window=0;
  _this->nend_bits=0;
  /*This is the offset from which od_ec_tell() will subtract partial bits.*/
  _this->nbits_total=33;
  _this->offs=0;
  _this->rng=0x80000000;
  _this->rem=-1;
  _this->val=0;
  _this->ext=0;
  _this->storage=_size;
  _this->error=0;
}

void od_ec_encode(od_ec_enc *_this,unsigned _fl,unsigned _fh,unsigned _ft){
  ogg_uint32_t r;
  r=_this->rng/_ft;
  if(_fl>0){
    _this->val+=_this->rng-r*(_ft-_fl);
    _this->rng=r*(_fh-_fl);
  }
  else _this->rng-=r*(_ft-_fh);
  od_ec_enc_normalize(_this);
}

void od_ec_encode_bin(od_ec_enc *_this,
 unsigned _fl,unsigned _fh,unsigned _bits){
  ogg_uint32_t r;
  r=_this->rng>>_bits;
  if(_fl>0){
    _this->val+=_this->rng-r*((1U<<_bits)-_fl);
    _this->rng=r*(_fh-_fl);
  }
  else _this->rng-=r*((1U<<_bits)-_fh);
  od_ec_enc_normalize(_this);
}

/*The probability of having a "one" is 1/(1<<_logp).*/
void od_ec_enc_bit_logp(od_ec_enc *_this,int _val,unsigned _logp){
  ogg_uint32_t r;
  ogg_uint32_t s;
  ogg_uint32_t l;
  r=_this->rng;
  l=_this->val;
  s=r>>_logp;
  r-=s;
  if(_val)_this->val=l+r;
  _this->rng=_val?s:r;
  od_ec_enc_normalize(_this);
}

void od_ec_enc_icdf(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ftb){
  ogg_uint32_t r;
  r=_this->rng>>_ftb;
  if(_s>0){
    _this->val+=_this->rng-r*_icdf[_s-1];
    _this->rng=r*(_icdf[_s-1]-_icdf[_s]);
  }
  else _this->rng-=r*_icdf[_s];
  od_ec_enc_normalize(_this);
}

void od_ec_enc_icdf_ft(od_ec_enc *_this,int _s,
 const unsigned char *_icdf,unsigned _ft){
  ogg_uint32_t r;
  r=_this->rng/_ft;
  if(_s>0){
    _this->val+=_this->rng-r*_icdf[_s-1];
    _this->rng=r*(_icdf[_s-1]-_icdf[_s]);
  }
  else _this->rng-=r*_icdf[_s];
  od_ec_enc_normalize(_this);
}

void od_ec_enc_icdf16_ft(od_ec_enc *_this,int _s,
 const unsigned short *_icdf,unsigned _ft){
  ogg_uint32_t r;
  r=_this->rng/_ft;
  if(_s>0){
    _this->val+=_this->rng-r*_icdf[_s-1];
    _this->rng=r*(_icdf[_s-1]-_icdf[_s]);
  }
  else _this->rng-=r*_icdf[_s];
  od_ec_enc_normalize(_this);
}

void od_ec_enc_icdf16(od_ec_enc *_this,int _s,
 const unsigned short *_icdf,unsigned _ftb){
  ogg_uint32_t r;
  r=_this->rng>>_ftb;
  if(_s>0){
    _this->val+=_this->rng-r*_icdf[_s-1];
    _this->rng=r*(_icdf[_s-1]-_icdf[_s]);
  }
  else _this->rng-=r*_icdf[_s];
  od_ec_enc_normalize(_this);
}

void od_ec_enc_uint(od_ec_enc *_this,ogg_uint32_t _fl,ogg_uint32_t _ft){
  unsigned  ft;
  unsigned  fl;
  int       ftb;
  if(_ft>1U<<OD_EC_UINT_BITS){
    _ft--;
    ftb=OD_ILOG_NZ(_ft)-OD_EC_UINT_BITS;
    ft=(_ft>>ftb)+1;
    fl=(unsigned)(_fl>>ftb);
    od_ec_encode(_this,fl,fl+1,ft);
    od_ec_enc_bits(_this,_fl&((ogg_uint32_t)1<<ftb)-1U,ftb);
  }
  else od_ec_encode(_this,_fl,_fl+1,_ft);
}

void od_ec_enc_bits(od_ec_enc *_this,ogg_uint32_t _fl,unsigned _bits){
  od_ec_window window;
  int          used;
  window=_this->end_window;
  used=_this->nend_bits;
  OD_ASSERT(_bits>0);
  if(used+_bits>OD_EC_WINDOW_SIZE){
    do{
      _this->error|=od_ec_write_byte_at_end(_this,(unsigned)window&0xFF);
      window>>=8;
      used-=8;
    }
    while(used>=8);
  }
  window|=(od_ec_window)_fl<<used;
  used+=_bits;
  _this->end_window=window;
  _this->nend_bits=used;
  _this->nbits_total+=_bits;
}

void od_ec_enc_patch_initial_bits(od_ec_enc *_this,
 unsigned _val,unsigned _nbits){
  int      shift;
  unsigned mask;
  OD_ASSERT(_nbits<=8);
  shift=8-_nbits;
  mask=(1<<_nbits)-1<<shift;
  if(_this->offs>0){
    /*The first byte has been finalized.*/
    _this->buf[0]=(unsigned char)(_this->buf[0]&~mask|_val<<shift);
  }
  else if(_this->rem>=0){
    /*The first byte is still awaiting carry propagation.*/
    _this->rem=_this->rem&~mask|_val<<shift;
  }
  else if(_this->rng<=0x80000000>>_nbits){
    /*The renormalization loop has never been run.*/
    _this->val=(_this->val&~((ogg_uint32_t)mask<<23))|
     (ogg_uint32_t)_val<<23+shift;
  }
  /*The encoder hasn't even encoded _nbits of data yet.*/
  else _this->error=-1;
}

void od_ec_enc_shrink(od_ec_enc *_this,ogg_uint32_t _size){
  OD_ASSERT(_this->offs+_this->end_offs<=_size);
  memmove(_this->buf+_size-_this->end_offs,
   _this->buf+_this->storage-_this->end_offs,_this->end_offs);
  _this->storage=_size;
}

void od_ec_enc_done(od_ec_enc *_this){
  od_ec_window window;
  int          used;
  ogg_uint32_t msk;
  ogg_uint32_t end;
  int          l;
  /*We output the minimum number of bits that ensures that the symbols encoded
     thus far will be decoded correctly regardless of the bits that follow.*/
  l=32-OD_ILOG_NZ(_this->rng);
  msk=0x7FFFFFFF>>l;
  end=_this->val+msk&~msk;
  if((end|msk)>=_this->val+_this->rng){
    l++;
    msk>>=1;
    end=_this->val+msk&~msk;
  }
  while(l>0){
    od_ec_enc_carry_out(_this,(int)(end>>23));
    end=end<<8&0x7FFFFFFF;
    l-=8;
  }
  /*If we have a buffered byte flush it into the output buffer.*/
  if(_this->rem>=0||_this->ext>0)od_ec_enc_carry_out(_this,0);
  /*If we have buffered extra bits, flush them as well.*/
  window=_this->end_window;
  used=_this->nend_bits;
  while(used>=8){
    _this->error|=od_ec_write_byte_at_end(_this,(unsigned)window&0xFF);
    window>>=8;
    used-=8;
  }
  /*Clear any excess space and add any remaining extra bits to the last byte.*/
  if(!_this->error){
    memset(_this->buf+_this->offs,0,
     _this->storage-_this->offs-_this->end_offs);
    if(used>0){
      /*If there's no range coder data at all, give up.*/
      if(_this->end_offs>=_this->storage)_this->error=-1;
      else{
        l=-l;
        /*If we've busted, don't add too many extra bits to the last byte; it
           would corrupt the range coder data, and that's more important.*/
        if(_this->offs+_this->end_offs>=_this->storage&&l<used){
          window&=(1<<l)-1;
          _this->error=-1;
        }
        _this->buf[_this->storage-_this->end_offs-1]|=(unsigned char)window;
      }
    }
  }
}
