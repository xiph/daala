#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

int main(int _argc,char **_argv){
  od_ec_enc      enc;
  od_ec_enc      enc_bak;
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
  uint32_t   ptr_sz;
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
  /*Trigger resize during termination.*/
  for(ft=2;ft<1024;ft++){
    for(i=0;i<ft;i++){
      od_ec_enc_init(&enc, (ft + i)&1);
      od_ec_enc_uint(&enc,i,ft);
      nbits=od_ec_enc_tell_frac(&enc);
      ptr=od_ec_enc_done(&enc,&ptr_sz);
      od_ec_dec_init(&dec,ptr,ptr_sz);
      sym=od_ec_dec_uint(&dec,ft, "test");
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ft of %i.\n",sym,i,ft);
        ret=EXIT_FAILURE;
      }
      nbits2=od_ec_dec_tell_frac(&dec);
      if(nbits!=nbits2){
        fprintf(stderr,"enc_tell_frac == %li, dec_tell_frac == %li\n",
         nbits,nbits2);
        ret=EXIT_FAILURE;
      }
      if(dec.error){
        fprintf(stderr,"uint error decoding %i with ft of %i.\n",i,ft);
        ret=EXIT_FAILURE;
      }
      od_ec_enc_clear(&enc);
    }
  }
  /*Raw bits only w/ resize*/
  for(ftb=1;ftb<17;ftb++){
    for(i=0;i<(1<<ftb);i++){
      od_ec_enc_init(&enc, (ftb + i)&1);
      od_ec_enc_checkpoint(&enc_bak,&enc);
      od_ec_enc_bits(&enc,i,ftb);
      od_ec_enc_rollback(&enc,&enc_bak);
      od_ec_enc_bits(&enc,i,ftb);
      ptr=od_ec_enc_done(&enc,&ptr_sz);
      if (ptr_sz != ((unsigned)ftb + 7) >> 3) {
        fprintf(stderr,"Used %li bytes to encode %i bits directly.\n",
         (long)ptr_sz,ftb);
        ret=EXIT_FAILURE;
      }
      od_ec_dec_init(&dec,ptr,ptr_sz);
      sym=od_ec_dec_bits(&dec,ftb, "test");
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ftb of %i.\n",sym,i,ftb);
        ret=EXIT_FAILURE;
      }
      od_ec_enc_clear(&enc);
    }
  }
  /*Testing unsigned integer corruption*/
  od_ec_enc_init(&enc,2);
  od_ec_enc_uint(&enc,128,129);
  od_ec_enc_checkpoint(&enc_bak,&enc);
  od_ec_enc_uint(&enc,128,129);
  od_ec_enc_uint(&enc,128,129);
  od_ec_enc_uint(&enc,128,129);
  od_ec_enc_rollback(&enc,&enc_bak);
  ptr=od_ec_enc_done(&enc,&ptr_sz);
  if(ptr_sz!=1){
    fprintf(stderr,"Incorrect output size %li.\n",(long)ptr_sz);
    ret=EXIT_FAILURE;
  }
  for(i=0;i<256;i++){
    ptr[ptr_sz-1]=i;
    od_ec_dec_init(&dec,ptr,ptr_sz);
    sym=od_ec_dec_uint(&dec,129, "test");
    if(i>=228 && i!=240 && !dec.error){
      fprintf(stderr,"Failed to detect uint error with %i.\n",i);
      ret=EXIT_FAILURE;
    }
    if(sym>=255){
      fprintf(stderr,"Corrupt uint out of range %i>=255 for %d.\n",sym,i);
      ret=EXIT_FAILURE;
    }
  }
  od_ec_enc_clear(&enc);
  /*Testing encoding of unsigned integers.*/
  od_ec_enc_init(&enc,1);
  for(ft=2;ft<1024;ft++){
    for(i=0;i<ft;i++){
      entropy+=log(ft)*M_LOG2E;
      od_ec_enc_checkpoint(&enc_bak,&enc);
      od_ec_enc_uint(&enc,0,ft);
      od_ec_enc_rollback(&enc,&enc_bak);
      od_ec_enc_uint(&enc,i,ft);
      od_ec_enc_checkpoint(&enc_bak,&enc);
      od_ec_enc_uint(&enc,1,ft);
      od_ec_enc_rollback(&enc,&enc_bak);
    }
    if(ft==512)ptr=od_ec_enc_done(&enc,&ptr_sz);
  }
  /*Testing encoding of raw bit values.*/
  for(ftb=1;ftb<16;ftb++){
    for(i=0;i<(1<<ftb);i++){
      entropy+=ftb;
      nbits=od_ec_enc_tell(&enc);
      od_ec_enc_bits(&enc,i,ftb);
      nbits2=od_ec_enc_tell(&enc);
      if(nbits2-nbits!=ftb){
        fprintf(stderr,"Used %li bits to encode %i bits directly.\n",
         nbits2-nbits,ftb);
        ret=EXIT_FAILURE;
      }
    }
  }
  nbits=od_ec_enc_tell_frac(&enc);
  ptr=od_ec_enc_done(&enc,&ptr_sz);
  fprintf(stderr,
   "Encoded %0.2f bits of entropy to %0.2f bits (%0.3f%% wasted).\n",
   entropy,ldexp(nbits,-3),100*(nbits-ldexp(entropy,3))/nbits);
  fprintf(stderr,"Packed to %li bytes.\n",(long)ptr_sz);
  od_ec_dec_init(&dec,ptr,ptr_sz);
  for(ft=2;ft<1024;ft++){
    for(i=0;i<ft;i++){
      sym=od_ec_dec_uint(&dec,ft, "test");
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ft of %i.\n",sym,i,ft);
        ret=EXIT_FAILURE;
      }
    }
  }
  for(ftb=1;ftb<16;ftb++){
    for(i=0;i<(1<<ftb);i++){
      sym=od_ec_dec_bits(&dec,ftb, "test");
      if(sym!=(unsigned)i){
        fprintf(stderr,"Decoded %i instead of %i with ftb of %i.\n",sym,i,ftb);
        ret=EXIT_FAILURE;
      }
    }
  }
  nbits2=od_ec_dec_tell_frac(&dec);
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
    tell[0]=od_ec_enc_tell_frac(&enc);
    for(j=0;j<sz;j++){
      if(zeros)data[j]=0;
      else data[j]=rand()%ft;
      od_ec_enc_uint(&enc,data[j],ft);
      tell[j+1]=od_ec_enc_tell_frac(&enc);
      if ((rand() & 7) == 0) {
        od_ec_enc_checkpoint(&enc_bak,&enc);
        od_ec_enc_uint(&enc,rand()&1?0:ft-1,ft);
        od_ec_enc_rollback(&enc,&enc_bak);
      }
    }
    if(!(rand()&1)){
      while(od_ec_enc_tell(&enc)&7)od_ec_enc_uint(&enc,rand()&1,2);
    }
    tell_bits=od_ec_enc_tell(&enc);
    ptr=od_ec_enc_done(&enc,&ptr_sz);
    if(tell_bits!=(unsigned)od_ec_enc_tell(&enc)){
      fprintf(stderr,"od_ec_enc_tell() changed after od_ec_enc_done(): "
       "%u instead of %u (Random seed: %u).\n",
       (unsigned)od_ec_enc_tell(&enc),tell_bits,seed);
      ret=EXIT_FAILURE;
    }
    if (((tell_bits + 7) >> 3) < ptr_sz) {
      fprintf(stderr,"od_ec_enc_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       ptr_sz, (tell_bits + 7) >> 3, seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,ptr_sz);
    if(od_ec_dec_tell_frac(&dec)!=tell[0]){
      fprintf(stderr,"od_ec_dec_tell() mismatch between encoder and decoder "
       "at symbol %i: %u instead of %u (Random seed: %u).\n",
       0,(unsigned)od_ec_dec_tell_frac(&dec),tell[0],seed);
      ret=EXIT_FAILURE;
    }
    for(j=0;j<sz;j++){
      sym=od_ec_dec_uint(&dec,ft, "test");
      if(sym!=data[j]){
        fprintf(stderr,"Decoded %i instead of %i with ft of %i "
         "at position %i of %i (Random seed: %u).\n",
         sym,data[j],ft,j,sz,seed);
        ret=EXIT_FAILURE;
      }
      if(od_ec_dec_tell_frac(&dec)!=tell[j+1]){
        fprintf(stderr,"od_ec_dec_tell() mismatch between encoder and decoder "
         "at symbol %i: %u instead of %u (Random seed: %u).\n",
         j+1,(unsigned)od_ec_dec_tell_frac(&dec),tell[j+1],seed);
        ret=EXIT_FAILURE;
      }
    }
    free(tell);
    free(data);
  }
  /*Test compatibility between multiple different encode/decode routines.*/
  for(i=0;i<409600;i++){
    unsigned *fz;
    unsigned *ftbs;
    unsigned *data;
    unsigned *tell;
    unsigned *enc_method;
    int       j;
    sz=rand()/((RAND_MAX>>(rand()%9U))+1U);
    fz=(unsigned *)malloc(sz*sizeof(*fz));
    ftbs=(unsigned *)malloc(sz*sizeof(*ftbs));
    data=(unsigned *)malloc(sz*sizeof(*data));
    tell=(unsigned *)malloc((sz+1)*sizeof(*tell));
    enc_method=(unsigned *)malloc(sz*sizeof(*enc_method));
    od_ec_enc_reset(&enc);
    tell[0]=od_ec_enc_tell_frac(&enc);
    for(j=0;j<sz;j++){
      data[j]=rand()/((RAND_MAX>>1)+1);
      ftbs[j]=(rand()%15)+1;
      fz[j] = (rand() % 32766) >> (15 - ftbs[j]);
      fz[j]=OD_MAXI(fz[j],1);
      enc_method[j]=rand()&1;
      switch(enc_method[j]){
        case 0:{
          if (rand() & 1) {
            od_ec_encode_bool_q15(&enc, data[j], fz[j] << (15 - ftbs[j]));
          }
          else {
            od_ec_encode_bool(&enc, data[j], fz[j] << (15 - ftbs[j]), 32768);
          }
        }break;
        case 1:{
          uint16_t cdf[2];
          cdf[0]=fz[j];
          cdf[1]=1U<<ftbs[j];
          od_ec_encode_cdf_unscaled_dyadic(&enc,data[j],cdf,2,ftbs[j]);
        }break;
      }
      tell[j+1]=od_ec_enc_tell_frac(&enc);
    }
    ptr=od_ec_enc_done(&enc,&ptr_sz);
    if (((od_ec_enc_tell(&enc) + 7U) >> 3) < ptr_sz) {
      fprintf(stderr,"od_ec_enc_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       ptr_sz, (od_ec_enc_tell(&enc) + 7) >> 3, seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,ptr_sz);
    if(od_ec_dec_tell_frac(&dec)!=tell[0]){
      fprintf(stderr,"od_ec_dec_tell() mismatch between encoder and decoder "
       "at symbol %i: %u instead of %u (Random seed: %u).\n",
       0,(unsigned)od_ec_dec_tell_frac(&dec),tell[0],seed);
      ret=EXIT_FAILURE;
    }
    for(j=0;j<sz;j++){
      int      dec_method;
      dec_method=rand()&1;
      switch(dec_method){
        case 0:{
          if (rand() & 1) {
            sym = od_ec_decode_bool_q15(&dec, fz[j] << (15 - ftbs[j]), "test");
          }
          else {
            sym = od_ec_decode_bool(&dec, fz[j]<< (15 - ftbs[j]), 32768,
             "test");
          }
        }break;
        case 1:{
          uint16_t cdf[2];
          cdf[0]=fz[j];
          cdf[1]=1U<<ftbs[j];
          sym=od_ec_decode_cdf_unscaled_dyadic(&dec,cdf,2,ftbs[j], "test");
        }break;
      }
      if(sym!=data[j]){
        fprintf(stderr,"Decoded %i instead of %i with fz=%i and ftb=%i "
         "at position %i of %i (Random seed: %u).\n",
         sym,data[j],fz[j],ftbs[j],j,sz,seed);
        fprintf(stderr,"Encoding method: %i, decoding method: %i\n",
         enc_method[j],dec_method);
        ret=EXIT_FAILURE;
      }
      if(od_ec_dec_tell_frac(&dec)!=tell[j+1]){
        fprintf(stderr,"od_ec_dec_tell() mismatch between encoder and decoder "
         "at symbol %i: %u instead of %u (Random seed: %u).\n",
         j+1,(unsigned)od_ec_dec_tell_frac(&dec),tell[j+1],seed);
        ret=EXIT_FAILURE;
      }
    }
    free(enc_method);
    free(tell);
    free(data);
    free(ftbs);
    free(fz);
  }
  od_ec_enc_reset(&enc);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,0,24576);
  od_ec_enc_patch_initial_bits(&enc,3,2);
  if(enc.error){
    fprintf(stderr,"od_ec_enc_patch_initial_bits() failed.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_patch_initial_bits(&enc,0,5);
  if(!enc.error){
    fprintf(stderr,
     "od_ec_enc_patch_initial_bits() didn't fail when it should have.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_reset(&enc);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,0,16384);
  od_ec_encode_bool_q15(&enc,1,32256);
  od_ec_encode_bool_q15(&enc,0,24576);
  od_ec_enc_patch_initial_bits(&enc,0,2);
  if(enc.error){
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
