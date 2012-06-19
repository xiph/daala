#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
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
  ptr=(unsigned char *)malloc(DATA_SIZE);
  od_ec_enc_init(&enc,ptr,DATA_SIZE);
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
  nbits=od_ec_tell_frac(&enc);
  od_ec_enc_done(&enc);
  fprintf(stderr,
   "Encoded %0.2f bits of entropy to %0.2f bits (%0.3f%% wasted).\n",
   entropy,ldexp(nbits,-3),100*(nbits-ldexp(entropy,3))/nbits);
  fprintf(stderr,"Packed to %li bytes.\n",(long)od_ec_range_bytes(&enc));
  od_ec_dec_init(&dec,ptr,DATA_SIZE);
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
  /*Testing an encoder bust prefers range coder data over raw bits.
    This isn't a general guarantee, will only work for data that is buffered in
     the encoder state and not yet stored in the user buffer, and should never
     get used in practice.
    It's mostly here for code coverage completeness.*/
  /*Start with a 16-bit buffer.*/
  od_ec_enc_init(&enc,ptr,2);
  /*Write 7 raw bits.*/
  od_ec_enc_bits(&enc,0x55,7);
  /*Write 12.3 bits of range coder data.*/
  od_ec_enc_uint(&enc,1,2);
  od_ec_enc_uint(&enc,1,3);
  od_ec_enc_uint(&enc,1,4);
  od_ec_enc_uint(&enc,1,5);
  od_ec_enc_uint(&enc,2,6);
  od_ec_enc_uint(&enc,6,7);
  od_ec_enc_done(&enc);
  od_ec_dec_init(&dec,ptr,2);
  if(!enc.error
   /*The raw bits should have been overwritten by the range coder data.*/
   ||od_ec_dec_bits(&dec,7)!=0x05
   /*And all the range coder data should have been encoded correctly.*/
   ||od_ec_dec_uint(&dec,2)!=1
   ||od_ec_dec_uint(&dec,3)!=1
   ||od_ec_dec_uint(&dec,4)!=1
   ||od_ec_dec_uint(&dec,5)!=1
   ||od_ec_dec_uint(&dec,6)!=2
   ||od_ec_dec_uint(&dec,7)!=6){
    fprintf(stderr,"Encoder bust overwrote range coder data with raw bits.\n");
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
    od_ec_enc_init(&enc,ptr,DATA_SIZE2);
    zeros=rand()%13==0;
    tell[0]=od_ec_tell_frac(&enc);
    for(j=0;j<sz;j++){
      if(zeros)data[j]=0;
      else data[j]=rand()%ft;
      od_ec_enc_uint(&enc,data[j],ft);
      tell[j+1]=od_ec_tell_frac(&enc);
    }
    if(!(rand()&1)){
      while(od_ec_tell(&enc)&7)od_ec_enc_uint(&enc,rand()&1,2);
    }
    tell_bits=od_ec_tell(&enc);
    od_ec_enc_done(&enc);
    if(tell_bits!=(unsigned)od_ec_tell(&enc)){
      fprintf(stderr,"od_ec_tell() changed after od_ec_enc_done(): "
       "%u instead of %u (Random seed: %u).\n",
       (unsigned)od_ec_tell(&enc),tell_bits,seed);
      ret=EXIT_FAILURE;
    }
    if(tell_bits+7>>3<od_ec_range_bytes(&enc)){
      fprintf(stderr,"od_ec_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       od_ec_range_bytes(&enc),tell_bits+7>>3,seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,DATA_SIZE2);
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
    unsigned *logp1;
    unsigned *data;
    unsigned *tell;
    unsigned *enc_method;
    int       j;
    sz=rand()/((RAND_MAX>>(rand()%9U))+1U);
    logp1=(unsigned *)malloc(sz*sizeof(*logp1));
    data=(unsigned *)malloc(sz*sizeof(*data));
    tell=(unsigned *)malloc((sz+1)*sizeof(*tell));
    enc_method=(unsigned *)malloc(sz*sizeof(*enc_method));
    od_ec_enc_init(&enc,ptr,DATA_SIZE2);
    tell[0]=od_ec_tell_frac(&enc);
    for(j=0;j<sz;j++){
      data[j]=rand()/((RAND_MAX>>1)+1);
      logp1[j]=(rand()%15)+1;
      enc_method[j]=rand()/((RAND_MAX>>2)+1);
      switch(enc_method[j]){
        case 0:{
          od_ec_encode(&enc,data[j]?(1<<logp1[j])-1:0,
           (1<<logp1[j])-(data[j]?0:1),1<<logp1[j]);
        }break;
        case 1:{
          od_ec_encode_bin(&enc,data[j]?(1<<logp1[j])-1:0,
           (1<<logp1[j])-(data[j]?0:1),logp1[j]);
        }break;
        case 2:{
          od_ec_enc_bit_logp(&enc,data[j],logp1[j]);
        }break;
        case 3:{
          unsigned char icdf[2];
          icdf[0]=1;
          icdf[1]=0;
          od_ec_enc_icdf(&enc,data[j],icdf,logp1[j]);
        }break;
      }
      tell[j+1]=od_ec_tell_frac(&enc);
    }
    od_ec_enc_done(&enc);
    if(od_ec_tell(&enc)+7U>>3<od_ec_range_bytes(&enc)){
      fprintf(stderr,"od_ec_tell() lied: "
       "there's %i bytes instead of %i (Random seed: %u).\n",
       od_ec_range_bytes(&enc),od_ec_tell(&enc)+7>>3,seed);
      ret=EXIT_FAILURE;
    }
    od_ec_dec_init(&dec,ptr,DATA_SIZE2);
    if(od_ec_tell_frac(&dec)!=tell[0]){
      fprintf(stderr,"od_ec_tell() mismatch between encoder and decoder "
       "at symbol %i: %u instead of %u (Random seed: %u).\n",
       0,(unsigned)od_ec_tell_frac(&dec),tell[0],seed);
    }
    for(j=0;j<sz;j++){
      int fs;
      int dec_method;
      dec_method=rand()/((RAND_MAX>>2)+1);
      switch(dec_method){
        case 0:{
          fs=od_ec_decode(&dec,1<<logp1[j]);
          sym=fs>=(1<<logp1[j])-1;
          od_ec_dec_update(&dec,sym?(1<<logp1[j])-1:0,
           (1<<logp1[j])-(sym?0:1),1<<logp1[j]);
        }break;
        case 1:{
          fs=od_ec_decode_bin(&dec,logp1[j]);
          sym=fs>=(1<<logp1[j])-1;
          od_ec_dec_update(&dec,sym?(1<<logp1[j])-1:0,
           (1<<logp1[j])-(sym?0:1),1<<logp1[j]);
        }break;
        case 2:{
          sym=od_ec_dec_bit_logp(&dec,logp1[j]);
        }break;
        case 3:{
          unsigned char icdf[2];
          icdf[0]=1;
          icdf[1]=0;
          sym=od_ec_dec_icdf(&dec,icdf,logp1[j]);
        }break;
      }
      if(sym!=data[j]){
        fprintf(stderr,"Decoded %i instead of %i with logp1 of %i "
         "at position %i of %i (Random seed: %u).\n",
         sym,data[j],logp1[j],j,sz,seed);
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
    free(logp1);
  }
  od_ec_enc_init(&enc,ptr,DATA_SIZE2);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,0,2);
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
  od_ec_enc_done(&enc);
  if(od_ec_range_bytes(&enc)!=1||ptr[0]!=192){
    fprintf(stderr,
     "Got %i when expecting 192 for od_ec_enc_patch_initial_bits().\n",ptr[0]);
    ret=EXIT_FAILURE;
  }
  od_ec_enc_init(&enc,ptr,DATA_SIZE2);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,0,1);
  od_ec_enc_bit_logp(&enc,1,6);
  od_ec_enc_bit_logp(&enc,0,2);
  od_ec_enc_patch_initial_bits(&enc,0,2);
  if(enc.error){
    fprintf(stderr,"od_ec_enc_patch_initial_bits() failed.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_done(&enc);
  if(od_ec_range_bytes(&enc)!=2||ptr[0]!=63){
    fprintf(stderr,
     "Got %i when expecting 63 for od_ec_enc_patch_initial_bits().\n",ptr[0]);
    ret=EXIT_FAILURE;
  }
  od_ec_enc_init(&enc,ptr,2);
  od_ec_enc_bit_logp(&enc,0,2);
  for(i=0;i<48;i++){
    od_ec_enc_bits(&enc,0,1);
  }
  od_ec_enc_done(&enc);
  if(!enc.error){
    fprintf(stderr,"Raw bits overfill didn't fail when it should have.\n");
    ret=EXIT_FAILURE;
  }
  od_ec_enc_init(&enc,ptr,2);
  for(i=0;i<17;i++){
    od_ec_enc_bits(&enc,0,1);
  }
  od_ec_enc_done(&enc);
  if(!enc.error){
    fprintf(stderr,"17 raw bits encoded in two bytes.\n");
    ret=EXIT_FAILURE;
  }
  free(ptr);
  return ret;
}
