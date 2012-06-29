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



/*Takes updated low and range values, renormalizes them so that
   32768<=rng<65536 (flushing bytes from low to the pre-carry buffer if
   necessary), and stores them back in the encoder context.
  _low: The new value of low.
  _rng: The new value of the range.*/
static void od_ec_enc_normalize(od_ec_enc *_this,
 od_ec_window _low,unsigned _rng){
  int d;
  int c;
  int s;
  c=_this->cnt;
  OD_ASSERT(_rng<=65535U);
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
    offs=_this->offs;
    if(offs+2>storage){
      storage=2*storage+2;
      buf=_ogg_realloc(buf,storage*sizeof(*buf));
      if(buf==NULL){
        _this->error=-1;
        _this->offs=0;
        return;
      }
      _this->precarry_buf=buf;
      _this->precarry_storage=storage;
    }
    c+=16;
    m=(1<<c)-1;
    if(s>=8){
      OD_ASSERT(offs<storage);
      buf[offs++]=(ogg_uint16_t)(_low>>c);
      _low&=m;
      c-=8;
      m>>=8;
    }
    OD_ASSERT(offs<storage);
    buf[offs++]=(ogg_uint16_t)(_low>>c);
    s=c+d-24;
    _low&=m;
    _this->offs=offs;
  }
  _this->low=_low<<d;
  _this->rng=_rng<<d;
  _this->cnt=s;
}


/*Initializes the encoder.
  _size: The initial size of the buffer, in bytes.*/
void od_ec_enc_init(od_ec_enc *_this,ogg_uint32_t _size){
  od_ec_enc_reset(_this);
  _this->buf=(unsigned char *)_ogg_malloc(_size*sizeof(*_this->buf));
  _this->storage=_size;
  if(_size>0&&_this->buf==NULL){
    _this->storage=0;
    _this->error=-1;
  }
  _this->precarry_buf=
   (ogg_uint16_t *)_ogg_malloc(_size*sizeof(*_this->precarry_buf));
  _this->precarry_storage=_size;
  if(_size>0&&_this->precarry_buf==NULL){
    _this->precarry_storage=0;
    _this->error=-1;
  }
}

/*Reinitializes the encoder.*/
void od_ec_enc_reset(od_ec_enc *_this){
  _this->end_offs=0;
  _this->end_window=0;
  _this->nend_bits=0;
  _this->offs=0;
  _this->low=0;
  _this->rng=0x8000;
  /*This is initialized to -9 so that it crosses zero after we've accumulated
     one byte + one carry bit.*/
  _this->cnt=-9;
  _this->error=0;
}

/*Frees the buffers used by the encoder.*/
void od_ec_enc_clear(od_ec_enc *_this){
  _ogg_free(_this->precarry_buf);
  _ogg_free(_this->buf);
}


/*Encodes a symbol given its scaled frequency information.
  The frequency information must be discernable by the decoder, assuming it
   has read only the previous symbols from the stream.
  You can change the frequency information, or even the entire source alphabet,
   so long as the decoder can tell from the context of the previously encoded
   information that it is supposed to do so as well.
  _fl: The cumulative frequency of all symbols that come before the one to be
        encoded.
  _fh: The cumulative frequency of all symbols up to and including the one to
        be encoded.
       Together with _fl, this defines the range [_fl,_fh) in which the
        decoded value will fall.
  _ft: The sum of the frequencies of all the symbols.
       This must be at least 16384, and no more than 32768.*/
static void od_ec_encode(od_ec_enc *_this,unsigned _fl,unsigned _fh,unsigned _ft){
  od_ec_window l;
  unsigned     r;
  int          s;
  unsigned     d;
  unsigned     u;
  unsigned     v;
  OD_ASSERT(_fl<_fh);
  OD_ASSERT(_fh<=_ft);
  OD_ASSERT(16384<=_ft);
  OD_ASSERT(_ft<=32768U);
  l=_this->low;
  r=_this->rng;
  OD_ASSERT(_ft<=r);
  s=r-_ft>=_ft;
  _ft<<=s;
  _fl<<=s;
  _fh<<=s;
  d=r-_ft;
  OD_ASSERT(d<_ft);
  u=_fl+OD_MINI(_fl,d);
  v=_fh+OD_MINI(_fh,d);
  r=v-u;
  l+=u;
  od_ec_enc_normalize(_this,l,r);
}

/*Equivalent to od_ec_encode() with _ft==32768.
  _fl: The cumulative frequency of all symbols that come before the one to be
        encoded.
  _fh: The cumulative frequency of all symbols up to and including the one to
        be encoded.*/
static void od_ec_encode_q15(od_ec_enc *_this,unsigned _fl,unsigned _fh){
  od_ec_window l;
  unsigned     r;
  unsigned     d;
  unsigned     u;
  unsigned     v;
  OD_ASSERT(_fl<_fh);
  OD_ASSERT(_fh<=32768U);
  l=_this->low;
  r=_this->rng;
  OD_ASSERT(32768U<=r);
  d=r-32768U;
  OD_ASSERT(d<32768U);
  u=_fl+OD_MINI(_fl,d);
  v=_fh+OD_MINI(_fh,d);
  r=v-u;
  l+=u;
  od_ec_enc_normalize(_this,l,r);
}

/*Encodes a symbol given its frequency information with an arbitrary scale.
  This operates just like od_ec_encode(), but does not require that _ft be at
   least 16384.
  _fl: The cumulative frequency of all symbols that come before the one to be
        encoded.
  _fh: The cumulative frequency of all symbols up to and including the one to
        be encoded.
  _ft: The sum of the frequencies of all the symbols.
       This must be at least 2 and no more than 32768.*/
static void od_ec_encode_unscaled(od_ec_enc *_this,
 unsigned _fl,unsigned _fh,unsigned _ft){
  int s;
  OD_ASSERT(_fl<_fh);
  OD_ASSERT(_fh<=_ft);
  OD_ASSERT(2<=_ft);
  OD_ASSERT(_ft<=32768U);
  s=15-OD_ILOG_NZ(_ft-1);
  od_ec_encode(_this,_fl<<s,_fh<<s,_ft<<s);
}

/*Encode a bit that has an _fz/_ft probability of being a zero.
  _val:  The value to encode (0 or 1).
  _fz:   The probability that _val is zero, scaled by _ft.
  _ft:   The total probability.
         This must be at least 16384 and no more than 32768.*/
void od_ec_encode_bool(od_ec_enc *_this,int _val,unsigned _fz,unsigned _ft){
  od_ec_window l;
  unsigned     r;
  int          s;
  unsigned     v;
  OD_ASSERT(0<_fz);
  OD_ASSERT(_fz<_ft);
  OD_ASSERT(16384<=_ft);
  OD_ASSERT(_ft<=32768U);
  l=_this->low;
  r=_this->rng;
  OD_ASSERT(_ft<=r);
  s=r-_ft>=_ft;
  _ft<<=s;
  _fz<<=s;
  OD_ASSERT(r-_ft<_ft);
  v=_fz+OD_MINI(_fz,r-_ft);
  if(_val)l+=v;
  r=_val?r-v:v;
  od_ec_enc_normalize(_this,l,r);
}

/*Equivalent to od_ec_encode_bool() with _ft==32768.
  _val:  The value to encode (0 or 1).
  _fz:   The probability that _val is zero, scaled by 32768.*/
void od_ec_encode_bool_q15(od_ec_enc *_this,int _val,unsigned _fz){
  od_ec_window l;
  unsigned     r;
  unsigned     v;
  OD_ASSERT(0<_fz);
  OD_ASSERT(_fz<32768U);
  l=_this->low;
  r=_this->rng;
  OD_ASSERT(32768U<=r);
  v=_fz+OD_MINI(_fz,r-32768U);
  if(_val)l+=v;
  r=_val?r-v:v;
  od_ec_enc_normalize(_this,l,r);
}

/*Encodes a symbol given a cumulative distribution function (CDF) table.
  _s:     The index of the symbol to encode.
  _cdf:   The CDF, such that symbol _s falls in the range
           [_s>0?_cdf[_s-1]:0,_cdf[_s]).
          The values must be monotonically non-decreasing, and the last value
           must be at least 16384, and no more than 32768.
  _nsyms: The number of symbols in the alphabet.
          This should be at most 16.*/
void od_ec_encode_cdf(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_cdf,int _nsyms){
  OD_ASSERT(_s>=0);
  OD_ASSERT(_s<_nsyms);
  od_ec_encode(_this,_s>0?_cdf[_s-1]:0,_cdf[_s],_cdf[_nsyms-1]);
}

/*Equivalent to od_ec_encode_cdf() with _cdf[_nsyms-1]==32768U.
  _s:     The index of the symbol to encode.
  _cdf:   The CDF, such that symbol _s falls in the range
           [_s>0?_cdf[_s-1]:0,_cdf[_s]).
          The values must be monotonically non-decreasing, and the last value
           must be exactly 32768.
  _nsyms: The number of symbols in the alphabet.
          This should be at most 16.*/
void od_ec_encode_cdf_q15(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_cdf,int _nsyms){
  OD_ASSERT(_s>=0);
  OD_ASSERT(_s<_nsyms);
  OD_ASSERT(_cdf[_nsyms-1]==32768U);
  od_ec_encode_q15(_this,_s>0?_cdf[_s-1]:0,_cdf[_s]);
}

/*Encodes a symbol given a cumulative distribution function (CDF) table.
  _s:     The index of the symbol to encode.
  _cdf:   The CDF, such that symbol _s falls in the range
           [_s>0?_cdf[_s-1]:0,_cdf[_s]).
          The values must be monotonically non-decreasing, and the last value
           must be at least 2, and no more than 32768.
  _nsyms: The number of symbols in the alphabet.
          This should be at most 16.*/
void od_ec_encode_cdf_unscaled(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_cdf,int _nsyms){
  OD_ASSERT(_s>=0);
  OD_ASSERT(_s<_nsyms);
  od_ec_encode_unscaled(_this,_s>0?_cdf[_s-1]:0,_cdf[_s],_cdf[_nsyms-1]);
}

/*Equivalent to od_ec_encode_cdf() with _cdf[_nsyms-1]==1<<_ftb.
  _s:     The index of the symbol to encode.
  _cdf:   The CDF, such that symbol _s falls in the range
           [_s>0?_cdf[_s-1]:0,_cdf[_s]).
          The values must be monotonically non-decreasing, and the last value
           must be exactly 1<<_ftb.
  _nsyms: The number of symbols in the alphabet.
          This should be at most 16.
  _ftb:   The number of bits of precision in the cumulative distribution.
          This must be no more than 15.*/
void od_ec_encode_cdf_unscaled_dyadic(od_ec_enc *_this,int _s,
 const ogg_uint16_t *_cdf,int _nsyms,unsigned _ftb){
  OD_ASSERT(_s>=0);
  OD_ASSERT(_s<_nsyms);
  OD_ASSERT(_ftb<=15);
  OD_ASSERT(_cdf[_nsyms-1]==1U<<_ftb);
  od_ec_encode_q15(_this,_s>0?_cdf[_s-1]<<15-_ftb:0,_cdf[_s]<<15-_ftb);
}

/*Encodes a raw unsigned integer in the stream.
  _fl: The integer to encode.
  _ft: The number of integers that can be encoded (one more than the max).
       This must be at least 2, and no more than 2**32-1.*/
void od_ec_enc_uint(od_ec_enc *_this,ogg_uint32_t _fl,ogg_uint32_t _ft){
  OD_ASSERT(_ft>=2);
  OD_ASSERT(_fl<_ft);
  if(_ft>1U<<OD_EC_UINT_BITS){
    int ft;
    int ftb;
    _ft--;
    ftb=OD_ILOG_NZ(_ft)-OD_EC_UINT_BITS;
    ft=(int)(_ft>>ftb)+1;
    od_ec_encode_cdf_q15(_this,(int)(_fl>>ftb),OD_UNIFORM_CDF_Q15(ft),ft);
    od_ec_enc_bits(_this,_fl&((ogg_uint32_t)1<<ftb)-1U,ftb);
  }
  else od_ec_encode_cdf_q15(_this,(int)_fl,OD_UNIFORM_CDF_Q15(_ft),(int)_ft);
}

/*Encodes a sequence of raw bits in the stream.
  _fl:  The bits to encode.
  _ftb: The number of bits to encode.
        This must be between 0 and 25, inclusive.*/
void od_ec_enc_bits(od_ec_enc *_this,ogg_uint32_t _fl,unsigned _ftb){
  od_ec_window end_window;
  int          nend_bits;
  OD_ASSERT(_ftb<=25);
  OD_ASSERT(_fl<(ogg_uint32_t)1<<_ftb);
  end_window=_this->end_window;
  nend_bits=_this->nend_bits;
  if(nend_bits+_ftb>OD_EC_WINDOW_SIZE){
    unsigned char *buf;
    ogg_uint32_t   storage;
    ogg_uint32_t   end_offs;
    buf=_this->buf;
    storage=_this->storage;
    end_offs=_this->end_offs;
    if(end_offs+(OD_EC_WINDOW_SIZE>>3)>=storage){
      unsigned char *new_buf;
      ogg_uint32_t   new_storage;
      new_storage=2*storage+(OD_EC_WINDOW_SIZE>>3);
      new_buf=(unsigned char *)_ogg_malloc(new_storage*sizeof(*new_buf));
      if(new_buf==NULL){
        _this->error=-1;
        _this->end_offs=0;
        return;
      }
      memcpy(new_buf+new_storage-end_offs,buf+storage-end_offs,
       end_offs*sizeof(*new_buf));
      storage=new_storage;
      _ogg_free(buf);
      _this->buf=buf=new_buf;
      _this->storage=storage;
    }
    do{
      OD_ASSERT(end_offs<storage);
      buf[storage-++(end_offs)]=(unsigned char)end_window;
      end_window>>=8;
      nend_bits-=8;
    }
    while(nend_bits>=8);
    _this->end_offs=end_offs;
  }
  OD_ASSERT(nend_bits+_ftb<=OD_EC_WINDOW_SIZE);
  end_window|=(od_ec_window)_fl<<nend_bits;
  nend_bits+=_ftb;
  _this->end_window=end_window;
  _this->nend_bits=nend_bits;
}

/*Overwrites a few bits at the very start of an existing stream, after they
   have already been encoded.
  This makes it possible to have a few flags up front, where it is easy for
   decoders to access them without parsing the whole stream, even if their
   values are not determined until late in the encoding process, without having
   to buffer all the intermediate symbols in the encoder.
  In order for this to work, at least _nbits bits must have already been
   encoded using probabilities that are an exact power of two.
  The encoder can verify the number of encoded bits is sufficient, but cannot
   check this latter condition.
  _val:   The bits to encode (in the least _nbits significant bits).
          They will be decoded in order from most-significant to least.
  _nbits: The number of bits to overwrite.
          This must be no more than 8.*/
void od_ec_enc_patch_initial_bits(od_ec_enc *_this,unsigned _val,int _nbits){
  int      shift;
  unsigned mask;
  OD_ASSERT(_nbits>=0);
  OD_ASSERT(_nbits<=8);
  OD_ASSERT(_val<1U<<_nbits);
  shift=8-_nbits;
  mask=(1U<<_nbits)-1<<shift;
  if(_this->offs>0){
    /*The first byte has been finalized.*/
    _this->precarry_buf[0]=
     (ogg_uint16_t)(_this->precarry_buf[0]&~mask|_val<<shift);
  }
  else if(9+_this->cnt+(_this->rng==0x8000)>_nbits){
    /*The first byte has yet to be output.*/
    _this->low=(_this->low&~((od_ec_window)mask<<16+_this->cnt))|
     (od_ec_window)_val<<16+_this->cnt+shift;
  }
  /*The encoder hasn't even encoded _nbits of data yet.*/
  else _this->error=-1;
}

/*Indicates that there are no more symbols to encode.
  All remaining output bytes are flushed to the output buffer.
  od_ec_enc_reset() should be called before using the encoder again.
  _bytes: Returns the size of the encoded data in the returned buffer.
  Return: A pointer to the start of the final buffer, or NULL if there was an
           encoding error.*/
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
  if(_this->error)return NULL;
  /*We output the minimum number of bits that ensures that the symbols encoded
     thus far will be decoded correctly regardless of the bits that follow.*/
  l=_this->low;
  r=_this->rng;
  c=_this->cnt;
  s=9;
  m=0x7FFF;
  e=l+m&~m;
  while((e|m)>=l+r){
    s++;
    m>>=1;
    e=l+m&~m;
  }
  s+=c;
  offs=_this->offs;
  buf=_this->precarry_buf;
  if(s>0){
    unsigned m;
    storage=_this->precarry_storage;
    if(offs+(s+7>>3)>storage){
      storage=storage*2+(s+7>>3);
      buf=(ogg_uint16_t *)_ogg_realloc(buf,storage*sizeof(*buf));
      if(buf==NULL){
        _this->error=-1;
        return NULL;
      }
      _this->precarry_buf=buf;
      _this->precarry_storage=storage;
    }
    m=(1<<c+16)-1;
    do{
      OD_ASSERT(offs<storage);
      buf[offs++]=(ogg_uint16_t)(e>>c+16);
      e&=m;
      s-=8;
      c-=8;
      m>>=8;
    }
    while(s>0);
  }
  /*Make sure there's enough room for the entropy-coded bits and the raw bits.*/
  out=_this->buf;
  storage=_this->storage;
  end_offs=_this->end_offs;
  e=_this->end_window;
  nend_bits=_this->nend_bits;
  s=-s;
  c=OD_MAXI(nend_bits-s+7>>3,0);
  if(offs+end_offs+c>storage){
    storage=offs+end_offs+c;
    out=(unsigned char *)_ogg_realloc(out,storage*sizeof(*out));
    if(out==NULL){
      _this->error=-1;
      return NULL;
    }
    memmove(out+storage-end_offs,out+_this->storage-end_offs,
     end_offs*sizeof(*out));
    _this->buf=out;
    _this->storage=storage;
  }
  /*If we have buffered raw bits, flush them as well.*/
  while(nend_bits>s){
    OD_ASSERT(end_offs<storage);
    out[storage-++(end_offs)]=(unsigned char)e;
    e>>=8;
    nend_bits-=8;
  }
  *_nbytes=offs+end_offs;
  /*Perform carry propagation.*/
  OD_ASSERT(offs+end_offs<=storage);
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
  OD_ASSERT(nend_bits==0||end_offs>0);
  if(nend_bits>0)out[end_offs-1]|=(unsigned char)e;
  /*Note: Unless there's an allocation error, if you keep encoding into the
     current buffer and call this function again later, everything will work
     just fine (you won't get a new packet out, but you will get a single
     buffer with the new data appended to the old).
    However, this function is O(N) where N is the amount of data coded so far,
     so calling it more than once for a given packet is a bad idea.*/
  return out;
}

/*Returns the number of bits "used" by the encoded symbols so far.
  This same number can be computed in either the encoder or the decoder, and is
   suitable for making coding decisions.
  Warning: The value returned by this function can decrease compared to an
   earlier call, even after encoding more data, if there is an encoding error
   (i.e., a failure to allocate enough space for the output buffer).
  Return: The number of bits.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
int od_ec_enc_tell(od_ec_enc *_this){
  /*The 10 here counteracts the offset of -9 baked into cnt, and adds 1 extra
     bit, which we reserve for terminating the stream.*/
  return (_this->offs+_this->end_offs)*8+_this->cnt+_this->nend_bits+10;
}

/*Returns the number of bits "used" by the encoded symbols so far.
  This same number can be computed in either the encoder or the decoder, and is
   suitable for making coding decisions.
  Warning: The value returned by this function can decrease compared to an
   earlier call, even after encoding more data, if there is an encoding error
   (i.e., a failure to allocate enough space for the output buffer).
  Return: The number of bits scaled by 2**OD_BITRES.
          This will always be slightly larger than the exact value (e.g., all
           rounding error is in the positive direction).*/
ogg_uint32_t od_ec_enc_tell_frac(od_ec_enc *_this){
  return od_ec_tell_frac(od_ec_enc_tell(_this),_this->rng);
}
