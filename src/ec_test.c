#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <string.h>

#define EC_MINI(_a,_b) ((_a)<(_b)?(_a):(_b))
#define EC_MAXI(_a,_b) ((_a)>(_b)?(_a):(_b))

# ifdef __GNUC_PREREQ
#  if __GNUC_PREREQ(3,4)
/*Note the casts to (int) below: this prevents EC_CLZ{32|64}_OFFS from
   "upgrading" the type of an entire expression to an (unsigned) size_t.*/
#   if INT_MAX>=2147483647
#    define EC_CLZ32_OFFS ((int)sizeof(unsigned)*CHAR_BIT)
#    define EC_CLZ32(_x) (__builtin_clz(_x))
#   elif LONG_MAX>=2147483647L
#    define EC_CLZ32_OFFS ((int)sizeof(unsigned long)*CHAR_BIT)
#    define EC_CLZ32(_x) (__builtin_clzl(_x))
#   endif
#  endif
# endif

# if defined(EC_CLZ32)
#  define EC_ILOGNZ_32(_v) (EC_CLZ32_OFFS-EC_CLZ32(_v))
#  define EC_ILOG_32(_v)   (EC_ILOGNZ_32(_v)&-!!(_v))
# else
#  error "Need __builtin_clz or equivalent."
# endif

typedef signed char    ec_bitree[32];
typedef unsigned char  ec_bitree_probs[16];
typedef unsigned short ec_probs[16];

#define EC_TREE_FTB (15)
#define EC_TREE_FT  (1<<EC_TREE_FTB)

typedef struct ec_tree       ec_tree;
typedef struct ec_code_entry ec_code_entry;



struct ec_tree{
  ec_bitree       t;
  ec_bitree_probs p;
  ec_probs        f;
  int             s;
};



struct ec_code_entry{
  int           tree;
  unsigned char v;
  unsigned char n;
  unsigned char l;
};



static ec_tree       *trees;
static int            ntrees;
static int            ctrees;
static ec_code_entry *entries;
static int            nentries;
static int            centries;



static int ec_tree_count_leaves(const signed char *_t,const unsigned char *_p,
 int _i){
  return _i<0?1:ec_tree_count_leaves(_t,_p,_t[2*_i+0])+
   ec_tree_count_leaves(_t,_p,_t[2*_i+1]);
}

static void tree_fill_freqs(unsigned short *_f,
 const signed char *_t,const unsigned char *_p,int _i,int _p0){
  int i0;
  int i1;
  int l;
  int p0;
  int p1;
  i0=_t[2*_i+0];
  p0=_p0*_p[_i]+128>>8;
  l=ec_tree_count_leaves(_t,_p,i0);
  p0=p0>l?p0:l;
  i1=_t[2*_i+1];
  l=ec_tree_count_leaves(_t,_p,i1);
  p0=p0<_p0-l?p0:_p0-l;
  if(i0<0)_f[-i0-1]=p0;
  else tree_fill_freqs(_f,_t,_p,i0,p0);
  p1=_p0-p0;
  if(i1<0)_f[-i1-1]=p1;
  else tree_fill_freqs(_f,_t,_p,i1,p1);
}

static void ec_tree_init(ec_tree *_tree,
 const signed char *_t,const unsigned char *_p,int _s){
  int l;
  memcpy(_tree->t,_t,2*_s*sizeof(*_t));
  memcpy(_tree->p,_p,_s*sizeof(*_p));
  tree_fill_freqs(_tree->f,_t,_p,0,EC_TREE_FT);
  for(l=1;l<=_s;l++)_tree->f[l]+=_tree->f[l-1];
  for(;l<16;l++)_tree->f[l]=_tree->f[l-1];
  _tree->s=_s;
}


static int tree_read(signed char *_t,unsigned char *_p){
  int s;
  int i;
  int j;
  i=j=0;
  s=1;
  do{
    int t0;
    int t1;
    int p0;
    if(scanf("%i%i%i",&t0,&t1,&p0)<3)return 0;
    if(t0<0)t0=-(++j);
    s=t0+1<s?s:t0+1;
    if(t1<0)t1=-(++j);
    s=t1+1<s?s:t1+1;
    _t[2*i+0]=(signed char)t0;
    _t[2*i+1]=(signed char)t1;
    _p[i]=p0;
  }
  while(++i<s);
  return s;
}

static int entry_read(void){
  signed char    t[64];
  unsigned char  p[32];
  int            s;
  int            tree;
  int            v;
  int            n;
  int            i;
  s=tree_read(t,p);
  if(s<=0)return s;
  for(tree=0;tree<ntrees;tree++){
    if(trees[tree].s==s
     &&memcmp(trees[tree].t,t,2*s*sizeof(*t))==0
     &&memcmp(trees[tree].p,p,s*sizeof(*p))==0){
      break;
    }
  }
  if(tree>=ntrees){
    if(ctrees<=ntrees){
      ctrees=ctrees<<1|1;
      trees=(ec_tree *)realloc(trees,ctrees*sizeof(*trees));
    }
    ec_tree_init(trees+tree,t,p,s);
    ntrees++;
  }
  if(centries<=nentries){
    centries=centries<<1|1;
    entries=(ec_code_entry *)realloc(entries,centries*sizeof(*entries));
  }
  entries[nentries].tree=tree;
  if(scanf("%i%i",&v,&n)<2)return 0;
  entries[nentries].v=(unsigned char)v;
  entries[nentries].n=(unsigned char)n;
  for(i=0;n-->0;)i=t[2*i+(v>>n&1)];
  entries[nentries].l=-i-1;
  nentries++;
  return 1;
}



typedef size_t ec_window;

typedef struct ec_enc ec_enc;
typedef struct ec_dec ec_dec;

#define EC_WINDOW_SZ (sizeof(ec_window)*CHAR_BIT)
/*This is meant to be a large, positive constant that can still be efficiently
   loaded as an immediate (on platforms like ARM, for example).
  Even relatively modest values like 100 would work fine.*/
#define EC_LOTS_OF_BITS (0x4000)



struct ec_enc{
  unsigned short *buf;
  size_t          cbuf;
  size_t          nbuf;
  unsigned        low;
  short           cnt;
  unsigned short  rng;
};



struct ec_dec{
  const unsigned char *buf;
  const unsigned char *end;
  ec_window            dif;
  short                cnt;
  unsigned short       rng;
};



#if defined(EC_MULTISYM)
static void ec_enc_init(ec_enc *_this){
  _this->buf=NULL;
  _this->cbuf=0;
  _this->nbuf=0;
  _this->low=0;
  _this->cnt=-9;
  _this->rng=0x8000;
}

static void ec_enc_clear(ec_enc *_this){
  free(_this->buf);
}


static int ec_encode(ec_enc *_this,unsigned _fl,unsigned _fh,unsigned _ftb){
  unsigned l;
  unsigned r;
  int      c;
  int      s;
  unsigned d;
  unsigned u;
  unsigned v;
#if defined(EC_MULT_FREE)
  unsigned ft;
#endif
  /*printf("0x%08X %2i 0x%04X  [0x%04X,0x%04X) {0x%04X}\n",
   _this->low,_this->cnt,_this->rng,_fl,_fh,1<<_ftb);*/
  l=_this->low;
  c=_this->cnt;
  r=_this->rng;
#if defined(EC_MULT_FREE)
  ft=1<<_ftb;
  s=r>=ft<<1;
  if(s)ft<<=1;
  _fl<<=s;
  _fh<<=s;
  d=r-ft;
  u=_fl+EC_MINI(_fl,d);
  v=_fh+EC_MINI(_fh,d);
#else
  u=r*_fl>>_ftb;
  v=r*_fh>>_ftb;
#endif
  r=v-u;
  l+=u;
  d=16-EC_ILOGNZ_32(r);
  s=c+d;
  /*TODO: Right now we flush every time we have at least one byte available.
    Instead we should use an ec_window and flush right before we're about to
     shift bits off the end of the window.
    For a 32-bit window this is about the same amount of work, but for a 64-bit
     window it should be a fair win.*/
  if(s>=0){
    unsigned short *buf;
    size_t          cbuf;
    size_t          nbuf;
    unsigned        m;
    buf=_this->buf;
    cbuf=_this->cbuf;
    nbuf=_this->nbuf;
    if(nbuf>=cbuf){
      cbuf=(cbuf<<1)+2;
      buf=(unsigned short *)realloc(buf,cbuf*sizeof(*buf));
      if(buf==NULL)return -ENOMEM;
    }
    c+=16;
    m=(1<<c)-1;
    if(s>=8){
      buf[nbuf++]=(unsigned short)(l>>c);
      l&=m;
      c-=8;
      m>>=8;
    }
    buf[nbuf++]=(unsigned short)(l>>c);
    s=c+d-24;
    l&=m;
    _this->buf=buf;
    _this->cbuf=cbuf;
    _this->nbuf=nbuf;
  }
  _this->low=l<<d;
  _this->cnt=s;
  _this->rng=r<<d;
  return 0;
}

static int ec_enc_done(ec_enc *_this,unsigned char **_out){
  unsigned short *buf;
  size_t          cbuf;
  size_t          nbuf;
  unsigned char  *out;
  unsigned        mask;
  unsigned        end;
  unsigned        l;
  unsigned        r;
  int             c;
  int             d;
  int             s;
  size_t          i;
  buf=_this->buf;
  cbuf=_this->cbuf;
  nbuf=_this->nbuf;
  l=_this->low;
  c=_this->cnt;
  r=_this->rng;
  /*Figure out the minimum number of bits that ensures that the symbols encoded
     thus far will be decoded correctly regardless of the bits that follow.*/
  d=0;
  mask=0x7FFF;
  end=l+mask&~mask;
  if((end|mask)>=l+r){
    d++;
    mask>>=1;
    end=l+mask&~mask;
  }
  /*Flush all remaining bits into the output buffer.*/
  s=c+d+8;
  if(s>=0){
    unsigned m;
    if(nbuf>=cbuf){
      cbuf=(cbuf<<1)+2;
      buf=(unsigned short *)realloc(buf,cbuf*sizeof(*buf));
      if(buf==NULL)return -ENOMEM;
    }
    m=(1<<c+16)-1;
    do{
      buf[nbuf++]=(unsigned short)(end>>c+16);
      end&=m;
      s-=8;
      c-=8;
      m>>=8;
    }
    while(s>=0);
    _this->buf=buf;
    _this->cbuf=cbuf;
  }
  /*Perform carry propagation.*/
  out=(unsigned char *)malloc(nbuf*sizeof(*out));
  if(out==NULL)return -ENOMEM;
  l=0;
  for(i=nbuf;i-->0;){
    l=buf[i]+l;
    out[i]=(unsigned char)l;
    l>>=8;
  }
  *_out=out;
  return nbuf;
}



static void ec_dec_init(ec_dec *_this,const unsigned char *_buf,size_t _sz){
  const unsigned char *end;
  ec_window            dif;
  int                  c;
  int                  s;
  end=_buf+_sz;
  dif=0;
  c=-15;
  for(s=EC_WINDOW_SZ-9;s>=0;){
    if(_buf>=end){
      c=EC_LOTS_OF_BITS;
      break;
    }
    c+=8;
    dif|=(ec_window)*_buf++<<s;
    s-=8;
  }
  _this->buf=_buf;
  _this->end=end;
  _this->dif=dif;
  _this->cnt=c;
  _this->rng=0x8000;
}

static int ec_decode(ec_dec *_this,const unsigned short _f[16],unsigned _ftb){
  ec_window dif;
  unsigned  r;
  unsigned  d;
  int       c;
  int       s;
  unsigned  u;
  unsigned  v;
  unsigned  q;
  int       l;
#if defined(EC_MULT_FREE)
  unsigned  fl;
  unsigned  fh;
  unsigned  ft;
#endif
  dif=_this->dif;
  c=_this->cnt;
  r=_this->rng;
#if defined(EC_MULT_FREE)
  ft=1<<_ftb;
  s=r>=ft<<1;
  if(s)ft<<=1;
  d=r-ft;
  q=EC_MAXI((int)(dif>>EC_WINDOW_SZ-16+1),(int)((dif>>EC_WINDOW_SZ-16)-d))>>s;
  fl=0;
  for(l=0;_f[l]<=q;l++)fl=_f[l];
  fh=_f[l];
  /*printf("0x%08X %2i 0x%04X  [0x%04X,0x%04X) {0x%04X}\n",
   (int)(_this->low>>EC_WINDOW_SZ-16),_this->cnt,_this->rng,fl,fh,1<<_ftb);*/
  fl<<=s;
  fh<<=s;
  u=fl+EC_MINI(fl,d);
  v=fh+EC_MINI(fh,d);
#else
  q=(unsigned)(dif>>EC_WINDOW_SZ-16);
  u=0;
  v=_f[0]*r>>_ftb;
  for(l=0;v<=q;){
    u=v;
    v=_f[++l]*r>>_ftb;
  }
#endif
  r=v-u;
  dif-=(ec_window)u<<EC_WINDOW_SZ-16;
  d=16-EC_ILOGNZ_32(r);
  s=c-d;
  if(s<0){
    const unsigned char *buf;
    const unsigned char *end;
    buf=_this->buf;
    end=_this->end;
    for(s=EC_WINDOW_SZ-9-(c+15);s>=0;){
      if(buf>=end){
        c=EC_LOTS_OF_BITS;
        break;
      }
      c+=8;
      dif|=(ec_window)*buf++<<s;
      s-=8;
    }
    s=c-d;
    _this->buf=buf;
  }
  _this->dif=dif<<d;
  _this->cnt=s;
  _this->rng=r<<d;
  return l;
}
#endif



#if defined(EC_BINARY)
static void ec_enc_bool_init(ec_enc *_this){
  _this->buf=NULL;
  _this->cbuf=0;
  _this->nbuf=0;
  _this->low=0;
  _this->cnt=-24;
  _this->rng=255;
}

static void ec_enc_bool_clear(ec_enc *_this){
  free(_this->buf);
}


static int ec_encode_bool(ec_enc *_this,int _b,unsigned _p){
  unsigned l;
  unsigned r;
  int      c;
  unsigned s;
  l=_this->low;
  c=_this->cnt;
  r=_this->rng;
  s=1+((r-1)*_p>>8);
  /*printf("0x%08X %3i 0x%02X  [0x%02X) %i {0x%02X}\n",
   _this->low,_this->cnt,_this->rng,s,_b,_p);*/
  r=_b?r-s:s;
  if(_b)l+=s;
  s=8-EC_ILOGNZ_32(r);
  r<<=s;
  c+=s;
  if(c>=0){
    unsigned short *buf;
    size_t          cbuf;
    size_t          nbuf;
    int             o;
    buf=_this->buf;
    cbuf=_this->cbuf;
    nbuf=_this->nbuf;
    if(nbuf>=cbuf){
      cbuf=(cbuf<<1)+1;
      buf=(unsigned short *)realloc(buf,cbuf*sizeof(*buf));
      if(buf==NULL)return -ENOMEM;
    }
    o=s-c;
    buf[nbuf++]=(unsigned short)(l>>24-o);
    l=l<<o&0xFFFFFF;
    s=c;
    c-=8;
    _this->buf=buf;
    _this->cbuf=cbuf;
    _this->nbuf=nbuf;
  }
  _this->low=l<<s;
  _this->cnt=c;
  _this->rng=r;
  return 0;
}

static int ec_encode_bool_tree(ec_enc *_this,
 const signed char *_t,const unsigned char *_p,int _v,int _n){
  int b;
  int i;
  int ret;
  for(i=0;_n-->0;i=_t[2*i+b]){
    b=_v>>_n&1;
    ret=ec_encode_bool(_this,b,_p[i]);
    if(ret<0)break;
  }
  return ret;
}

static int ec_enc_bool_done(ec_enc *_this,unsigned char **_out){
  unsigned short *buf;
  size_t          nbuf;
  unsigned char  *out;
  unsigned        l;
  size_t          i;
  /*Flush all remaining bits into the output buffer.
    TODO: Flush minimum amount only; don't update state.
    For right now we are copying the VP8 approach.*/
  for(i=0;i<32;i++)ec_encode_bool(_this,0,128);
  buf=_this->buf;
  nbuf=_this->nbuf;
  /*Perform carry propagation.*/
  out=(unsigned char *)malloc(nbuf*sizeof(*out));
  if(out==NULL)return -ENOMEM;
  l=0;
  for(i=nbuf;i-->0;){
    l=buf[i]+l;
    out[i]=(unsigned char)l;
    l>>=8;
  }
  *_out=out;
  return nbuf;
}



static void ec_dec_bool_init(ec_dec *_this,
 const unsigned char *_buf,size_t _sz){
  const unsigned char *end;
  ec_window            dif;
  int                  c;
  int                  s;
  end=_buf+_sz;
  dif=0;
  c=-8;
  for(s=EC_WINDOW_SZ-8;s>=0;){
    if(_buf>=end){
      c=EC_LOTS_OF_BITS;
      break;
    }
    c+=8;
    dif|=(ec_window)*_buf++<<s;
    s-=8;
  }
  _this->buf=_buf;
  _this->end=end;
  _this->dif=dif;
  _this->cnt=c;
  _this->rng=255;
}

static int ec_decode_bool(ec_dec *_this,unsigned _p){
  ec_window d;
  ec_window q;
  unsigned  r;
  int       c;
  int       s;
  int       v;
  d=_this->dif;
  c=_this->cnt;
  r=_this->rng;
  s=1+((r-1)*_p>>8);
  q=(ec_window)s<<EC_WINDOW_SZ-8;
  v=d>=q;
  /*printf("0x%08X %3i 0x%02X  [0x%02X) %i {0x%02X}\n",
   (int)(_this->dif>>EC_WINDOW_SZ-8),_this->cnt,_this->rng,s,v,_p);*/
  r=v?r-s:s;
  if(v)d-=q;
  s=8-EC_ILOGNZ_32(r);
  r<<=s;
  d<<=s;
  c-=s;
  if(c<0){
    const unsigned char *buf;
    const unsigned char *end;
    buf=_this->buf;
    end=_this->end;
    for(s=EC_WINDOW_SZ-8-(c+8);s>=0;){
      if(buf>=end){
        c=EC_LOTS_OF_BITS;
        break;
      }
      c+=8;
      d|=(ec_window)*buf++<<s;
      s-=8;
    }
    _this->buf=buf;
  }
  _this->dif=d;
  _this->cnt=c;
  _this->rng=r;
  return v;
}

static int ec_decode_bool_tree(ec_dec *_this,
 const signed char *_t,const unsigned char *_p){
  int i;
  i=0;
  do i=_t[2*i+ec_decode_bool(_this,_p[i])];
  while(i>0);
  return -i-1;
}
#endif



int main(void){
  ec_enc         enc;
  ec_dec         dec;
  unsigned char *out;
  int            nout;
  int            ei;
  int            ret;
  int            loop;
  while(entry_read());
  for(loop=0;loop<1024;loop++){
#if defined(EC_MULTISYM)
    ec_enc_init(&enc);
    for(ei=0;ei<nentries;ei++){
      const unsigned short *f;
      int                   l;
      f=trees[entries[ei].tree].f;
      l=entries[ei].l;
      ec_encode(&enc,l>0?f[l-1]:0,f[l],EC_TREE_FTB);
    }
    nout=ec_enc_done(&enc,&out);
    if(nout<0)return EXIT_FAILURE;
    /*printf("Encoded to %u bytes.\n",nout);*/
    ec_dec_init(&dec,out,nout);
    ec_enc_clear(&enc);
    ret=EXIT_SUCCESS;
    for(ei=0;ei<nentries;ei++){
      const unsigned short *f;
      int                   l;
      f=trees[entries[ei].tree].f;
      l=ec_decode(&dec,f,EC_TREE_FTB);
      if(l!=entries[ei].l){
        printf("%i should be %i\n",l,entries[ei].l);
        ret=EXIT_FAILURE;
      }
    }
    free(out);
#elif defined(EC_BINARY)
    ec_enc_bool_init(&enc);
    for(ei=0;ei<nentries;ei++){
      const signed char   *t;
      const unsigned char *p;
      int                  v;
      int                  n;
      t=trees[entries[ei].tree].t;
      p=trees[entries[ei].tree].p;
      v=entries[ei].v;
      n=entries[ei].n;
      ec_encode_bool_tree(&enc,t,p,v,n);
    }
    nout=ec_enc_bool_done(&enc,&out);
    if(nout<0)return EXIT_FAILURE;
    /*printf("Encoded to %u bytes.\n",nout);*/
    ec_dec_bool_init(&dec,out,nout);
    ec_enc_bool_clear(&enc);
    ret=EXIT_SUCCESS;
    for(ei=0;ei<nentries;ei++){
      const signed char   *t;
      const unsigned char *p;
      int                  l;
      t=trees[entries[ei].tree].t;
      p=trees[entries[ei].tree].p;
      l=ec_decode_bool_tree(&dec,t,p);
      if(l!=entries[ei].l){
        printf("%i should be %i\n",l,entries[ei].l);
        ret=EXIT_FAILURE;
      }
    }
    free(out);
#else
# error "Please define EC_MULTISYM or EC_BINARY"
#endif
  }
  return ret;
}
