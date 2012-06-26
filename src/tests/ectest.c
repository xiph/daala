#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../internal.c"
#include "../entcode.c"
#include "../entdec.c"
#include "../entenc.c"

#if !defined(M_LOG2E)
# define M_LOG2E (1.4426950408889634074)
#endif
#define DATA_SIZE  (10000000)
#define DATA_SIZE2 (10000)

int main(int _argc,char **_argv){
  od_ec_enc      enc;
  od_ec_dec      dec;
  long           nbits;
  long           nbits2;
  double         entropy;
  int            ft;
  int            ftb;
  int            sz;
  int            i;
  int            ret;
  unsigned int   sym;
  unsigned int   seed;
  unsigned char *ptr;
  ogg_uint32_t   ptr_sz;
  const char    *env_seed;
  ret=EXIT_SUCCESS;
  entropy=0;
  if(_argc>2){
    fprintf(stderr,"Usage: %s [<seed>]\n",_argv[0]);
    return EXIT_FAILURE;
  }
  env_seed=getenv("SEED");
  if(_argc>1)seed=atoi(_argv[1]);
  else if(env_seed)seed=atoi(env_seed);
  else seed=time(NULL);
  /*Testing encoding of raw bit values.*/
  od_ec_enc_init(&enc,DATA_SIZE);
  for(ft=2;ft<1024;ft++){
    for(i=0;i<ft;i++){
      entropy+=log(ft)*M_LOG2E;
      od_ec_enc_uint(&enc,i,ft);
    }
  }
  /*Testing encoding of raw bit values.*/
  for(ftb=1;ftb<16;ftb++){
    for(i=0;i<(1<<ftb);i++){
      entropy+=ftb;
      nbits=od_ec_tell(&enc);
      od_ec_enc_bits(&enc,i,ftb);
      nbits2=od_ec_tell(&enc);
      if(nbits2-nbits!=ftb){
        fprintf(stderr,"Used %li bits to encode %i bits directly.\n",
         nbits2-nbits,ftb);
        ret=EXIT_FAILURE;
      }
    }
  }
  nbits=od_ec_tell_frac(&enc.base);
  ptr=od_ec_enc_done(&enc,&ptr_sz);
  fprintf(stderr,
   "Encoded %0.2f bits of entropy to %0.2f bits (%0.3f%% wasted).\n",
   entropy,ldexp(nbits,-3),100*(nbits-ldexp(entropy,3))/nbits);
  fprintf(stderr,"Packed to %li bytes.\n",(long)ptr_sz);
  od_ec_dec_init(&dec,ptr,ptr_sz);
  for(ft=2;ft<1024;ft++){
    for(i=0;i<ft;i++){
      sym=od_ec_dec_uint(&dec,ft);
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ft of %i.\n",sym,i,ft);
        ret=EXIT_FAILURE;
      }
    }
  }
  for(ftb=1;ftb<16;ftb++){
    for(i=0;i<(1<<ftb);i++){
      sym=od_ec_dec_bits(&dec,ftb);
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ftb of %i.\n",sym,i,ftb);
        ret=EXIT_FAILURE;
      }
    }
  }
  nbits2=od_ec_tell_frac(&dec);
  if(nbits!=nbits2){
    fprintf(stderr,
     "Reported number of bits used was %0.2f, should be %0.2f.\n",
     ldexp(nbits2,-3),ldexp(nbits,-3));
    ret=EXIT_FAILURE;
  }
  srand(seed);
  fprintf(stderr,"Testing random streams... Random seed: %u (%.4X).\n",
   seed,rand()&65535);
  for(i=0;i<409600;i++){
    unsigned *data;
    unsigned *tell;
    unsigned  tell_bits;
    int       j;
    int zeros;
    ft=rand()/((RAND_MAX>>(rand()%11U))+1U)+10;
    sz=rand()/((RAND_MAX>>(rand()%9U))+1U);
    data=(unsigned *)malloc(sz*sizeof(*data));
    tell=(unsigned *)malloc((sz+1)*sizeof(*tell));
    od_ec_enc_reset(&enc);
    zeros=rand()%13==0;
    tell[0]=od_ec_tell_frac(&enc.base);
    for(j=0;j<sz;j++){
      if(zeros)data[j]=0;
      else data[j]=rand()%ft;
      od_ec_enc_uint(&enc,data[j],ft);
      tell[j+1]=od_ec_tell_frac(&enc.base);
    }
    if(!(rand()&1)){
      while(od_ec_tell(&enc)&7)od_ec_enc_uint(&enc,rand()&1,2);
    }
    tell_bits=od_ec_tell(&enc);
    ptr=od_ec_enc_done(&enc,&ptr_sz);
    if(tell_bits!=(unsigned)od_ec_tell(&enc)){
      fprintf(stderr,"od_ec_tell() changed after od_ec_enc_done(): "
       "%u instead of %u (Random seed: %u).\n",
       (unsigned)od_ec_tell(&enc),tell_bits,seed);
      ret=EXIT_FAILURE;
    }
    if(tell_bits+7>>3<ptr_sz){
      fprintf(stderr,"od_ec_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       ptr_sz,tell_bits+7>>3,seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,ptr_sz);
    if(od_ec_tell_frac(&dec)!=tell[0]){
      fprintf(stderr,"od_ec_tell() mismatch between encoder and decoder "
       "at symbol %i: %u instead of %u (Random seed: %u).\n",
       0,(unsigned)od_ec_tell_frac(&dec),tell[0],seed);
    }
    for(j=0;j<sz;j++){
      sym=od_ec_dec_uint(&dec,ft);
      if(sym!=data[j]){
        fprintf(stderr,"Decoded %i instead of %i with ft of %i "
         "at position %i of %i (Random seed: %u).\n",
         sym,data[j],ft,j,sz,seed);
        ret=EXIT_FAILURE;
      }
      if(od_ec_tell_frac(&dec)!=tell[j+1]){
        fprintf(stderr,"od_ec_tell() mismatch between encoder and decoder "
         "at symbol %i: %u instead of %u (Random seed: %u).\n",
         j+1,(unsigned)od_ec_tell_frac(&dec),tell[j+1],seed);
      }
    }
    free(tell);
    free(data);
  }
  /*Test compatibility between multiple different encode/decode routines.*/
  for(i=0;i<409600;i++){
    unsigned *fz;
    unsigned *ftb;
    unsigned *data;
    unsigned *tell;
    unsigned *enc_method;
    int       j;
    sz=rand()/((RAND_MAX>>(rand()%9U))+1U);
    fz=(unsigned *)malloc(sz*sizeof(*fz));
    ftb=(unsigned *)malloc(sz*sizeof(*ftb));
    data=(unsigned *)malloc(sz*sizeof(*data));
    tell=(unsigned *)malloc((sz+1)*sizeof(*tell));
    enc_method=(unsigned *)malloc(sz*sizeof(*enc_method));
    od_ec_enc_reset(&enc);
    tell[0]=od_ec_tell_frac(&enc.base);
    for(j=0;j<sz;j++){
      data[j]=rand()/((RAND_MAX>>1)+1);
      ftb[j]=(rand()%15)+1;
      fz[j]=rand()%32766>>15-ftb[j];
      fz[j]=OD_MAXI(fz[j],1);
      enc_method[j]=rand()/((RAND_MAX>>2)+1);
      switch(enc_method[j]){
        case 0:{
          od_ec_encode(&enc,data[j]?fz[j]:0,
           data[j]?1U<<ftb[j]:fz[j],1U<<ftb[j]);
        }break;
        case 1:{
          od_ec_encode_bin(&enc,data[j]?fz[j]:0,
           data[j]?1U<<ftb[j]:fz[j],ftb[j]);
        }break;
        case 2:{
          od_ec_enc_bool(&enc,data[j],fz[j]<<15-ftb[j]);
        }break;
        case 3:{
          ogg_uint16_t icdf[2];
          icdf[0]=(1U<<ftb[j])-fz[j];
          icdf[1]=0;
          od_ec_enc_icdf16(&enc,data[j],icdf,ftb[j]);
        }break;
      }
      tell[j+1]=od_ec_tell_frac(&enc.base);
    }
    ptr=od_ec_enc_done(&enc,&ptr_sz);
    if(od_ec_tell(&enc)+7U>>3<ptr_sz){
      fprintf(stderr,"od_ec_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       ptr_sz,od_ec_tell(&enc)+7>>3,seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,ptr_sz);
    if(od_ec_tell_frac(&dec)!=tell[0]){
      fprintf(stderr,"od_ec_tell() mismatch between encoder and decoder "
       "at symbol %i: %u instead of %u (Random seed: %u).\n",
       0,(unsigned)od_ec_tell_frac(&dec),tell[0],seed);
    }
    for(j=0;j<sz;j++){
      unsigned fs;
      int      dec_method;
      dec_method=rand()/((RAND_MAX>>2)+1);
      switch(dec_method){
        case 0:{
          fs=od_ec_decode(&dec,1U<<ftb[j]);
          sym=fs>=fz[j];
          od_ec_dec_update(&dec,sym?fz[j]:0,
           sym?1U<<ftb[j]:fz[j],1U<<ftb[j]);
        }break;
        case 1:{
          fs=od_ec_decode_bin(&dec,ftb[j]);
          sym=fs>=fz[j];
          od_ec_dec_update(&dec,sym?fz[j]:0,
           sym?1U<<ftb[j]:fz[j],1U<<ftb[j]);
        }break;
        case 2:{
          sym=od_ec_dec_bool(&dec,fz[j]<<15-ftb[j]);
        }break;
        case 3:{
          ogg_uint16_t icdf[2];
          icdf[0]=(1U<<ftb[j])-fz[j];
          icdf[1]=0;
          sym=od_ec_dec_icdf16(&dec,icdf,ftb[j]);
        }break;
      }
      if(sym!=data[j]){
        fprintf(stderr,"Decoded %i instead of %i with fz=%i and ftb=%i "
         "at position %i of %i (Random seed: %u).\n",
         sym,data[j],fz[j],ftb[j],j,sz,seed);
        fprintf(stderr,"Encoding method: %i, decoding method: %i\n",
         enc_method[j],dec_method);
        ret=EXIT_FAILURE;
      }
      if(od_ec_tell_frac(&dec)!=tell[j+1]){
        fprintf(stderr,"od_ec_tell() mismatch between encoder and decoder "
         "at symbol %i: %u instead of %u (Random seed: %u).\n",
         j+1,(unsigned)od_ec_tell_frac(&dec),tell[j+1],seed);
      }
    }
    free(enc_method);
    free(tell);
    free(data);
    free(ftb);
    free(fz);
  }
  od_ec_enc_reset(&enc);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,0,24576);
  od_ec_enc_patch_initial_bits(&enc,3,2);
  if(enc.base.error){
    fprintf(stderr,"od_ec_enc_patch_initial_bits() failed.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_patch_initial_bits(&enc,0,5);
  if(!enc.base.error){
    fprintf(stderr,
     "od_ec_enc_patch_initial_bits() didn't fail when it should have.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_reset(&enc);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,0,16384);
  od_ec_enc_bool(&enc,1,32256);
  od_ec_enc_bool(&enc,0,24576);
  od_ec_enc_patch_initial_bits(&enc,0,2);
  if(enc.base.error){
    fprintf(stderr,"od_ec_enc_patch_initial_bits() failed.\n");
    ret=EXIT_FAILURE;
  }
  ptr=od_ec_enc_done(&enc,&ptr_sz);
  if(ptr_sz!=2||ptr[0]!=63){
    fprintf(stderr,
     "Got %i when expecting 63 for od_ec_enc_patch_initial_bits().\n",ptr[0]);
    ret=EXIT_FAILURE;
  }
  od_ec_enc_clear(&enc);
  return ret;
}
