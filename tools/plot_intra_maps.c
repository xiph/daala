/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
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
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ogg/os_types.h>
#include "image.h"

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
    intra_map_colors(colors,mapi_max);
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
