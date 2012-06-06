#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ogg/os_types.h>
#include "image.h"

#define HUE_MAX (0xFFFF*6)

static void rgba16_from_hue(od_rgba16_pixel _color,int _hue){
  int            h;
  int            i;
  int            f;
  int            y;
  h=_hue%HUE_MAX;
  if(h<0)h+=HUE_MAX;
  i=h/0xFFFF;
  f=h-i*0xFFFF;
  y=(0xFFFFU*f+0x7FFFU)/0xFFFFU;
  switch(i){
    case 0:{
      _color[0]=(unsigned short)0xFFFFU;
      _color[1]=(unsigned short)y;
      _color[2]=0;
      _color[3]=(unsigned short)0xFFFFU;
    }break;
    case 1:{
      _color[0]=(unsigned short)(0xFFFFU-y);
      _color[1]=(unsigned short)0xFFFFU;
      _color[2]=0;
      _color[3]=(unsigned short)0xFFFFU;
    }break;
    case 2:{
      _color[0]=0;
      _color[1]=(unsigned short)0xFFFFU;
      _color[2]=(unsigned short)y;
      _color[3]=(unsigned short)0xFFFFU;
    }break;
    case 3:{
      _color[0]=0;
      _color[1]=(unsigned short)(0xFFFFU-y);
      _color[2]=(unsigned short)0xFFFFU;
      _color[3]=(unsigned short)0xFFFFU;
    }break;
    case 4:{
      _color[0]=(unsigned short)y;
      _color[1]=0;
      _color[2]=(unsigned short)0xFFFFU;
      _color[3]=(unsigned short)0xFFFFU;
    }break;
    default:{
      _color[0]=(unsigned short)0xFFFFU;
      _color[1]=0;
      _color[2]=(unsigned short)(0xFFFFU-y);
      _color[3]=(unsigned short)0xFFFFU;
    }break;
  }
}

static int get_dimensions_from_filename(int *_nxblocks,int *_nyblocks,
 const char *_fin_name){
  const char *p;
  p=strstr(_fin_name,".intra.map");
  if(p!=NULL){
    const char *q;
    do{
      q=p;
      p=strstr(p+1,".intra.map");
    }
    while(p!=NULL);
    while(q>_fin_name&&*(q-1)>='0'&&*(q-1)<='9')q--;
    if(q>_fin_name&&*--q=='x'){
      while(q>_fin_name&&*(q-1)>='0'&&*(q-1)<='9')q--;
      if(sscanf(q,"%dx%d",_nxblocks,_nyblocks)==2){
        return 0;
      }
    }
  }
  return -1;
}

int main(int _argc,const char **_argv){
  int ai;
  for(ai=1;ai<_argc;ai++){
    char             fout_name[8192];
    unsigned char   *map;
    FILE            *fin;
    FILE            *fout;
    od_rgba16_pixel *colors;
    od_rgba16_image  image;
    int              nxblocks;
    int              nyblocks;
    int              mapi_max;
    int              mapi;
    int              bi;
    int              bj;
    /*Try to extract the map size from the file name.*/
    if(get_dimensions_from_filename(&nxblocks,&nyblocks,_argv[ai])<0){
      fprintf(stderr,"Could not determine dimensions of '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    map=(unsigned char *)_ogg_malloc(nxblocks*(size_t)nyblocks*sizeof(*map));
    fin=fopen(_argv[ai],"rb");
    if(fin==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    if(fread(map,nxblocks*(size_t)nyblocks,1,fin)<1){
      fprintf(stderr,"Error reading from input file '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    fclose(fin);
    /*Count how many different intra modes we tried.*/
    mapi_max=0;
    for(bj=0;bj<nyblocks;bj++){
      for(bi=0;bi<nxblocks;bi++){
        if(map[bj*nxblocks+bi]>mapi_max)mapi_max=map[bj*nxblocks+bi];
      }
    }
    mapi_max++;
    colors=(od_rgba16_pixel *)_ogg_malloc(mapi_max*sizeof(*colors));
    /*The first mode is always DC mode; use gray.*/
    colors[0][0]=(unsigned short)0x8080U;
    colors[0][1]=(unsigned short)0x8080U;
    colors[0][2]=(unsigned short)0x8080U;
    colors[0][3]=(unsigned short)0xFFFFU;
    /*Pull out fully saturated colors from the color wheel for all the
       directional modes.*/
    if(mapi_max>1){
      int dhue;
      dhue=HUE_MAX/(mapi_max-1);
      for(mapi=0;mapi<mapi_max;mapi++){
        rgba16_from_hue(colors[mapi+1],dhue*mapi);
      }
    }
    od_rgba16_image_init(&image,nxblocks,nyblocks);
    for(bj=0;bj<nyblocks;bj++){
      for(bi=0;bi<nxblocks;bi++){
        od_rgba16_image_draw_point(&image,bi,bj,colors[map[bj*nxblocks+bi]]);
      }
    }
    sprintf(fout_name,"%s.png",_argv[ai]);
    fout=fopen(fout_name,"wb");
    if(fout==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",fout_name);
      return EXIT_FAILURE;
    }
    od_rgba16_image_write_png(&image,fout);
    fclose(fout);
    od_rgba16_image_clear(&image);
    _ogg_free(colors);
    _ogg_free(map);
  }
  return EXIT_SUCCESS;
}
