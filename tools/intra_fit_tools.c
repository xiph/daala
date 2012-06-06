#include <stdlib.h>
#include <string.h>
#include "intra_fit_tools.h"

/*Computes the starting offset and number of blocks which can be intra
   predicted with full context (i.e., all of their neighbors) available.*/
void get_intra_dims(const th_info *_ti,int _pli,
 int *_x0,int *_y0,int *_nxblocks,int *_nyblocks){
  int xshift;
  int yshift;
  xshift=_pli!=0&&!(_ti->pixel_fmt&1);
  yshift=_pli!=0&&!(_ti->pixel_fmt&2);
  /*An offset of 1 would be fine to provide enough context for VP8-style intra
     prediction, but for frequency-domain prediction, we'll want a full block,
     plus overlap.*/
  *_x0=(_ti->pic_x>>xshift)+2*B_SZ;
  *_y0=(_ti->pic_y>>yshift)+2*B_SZ;
  /*We take an extra block off the end to give enough context for above-right
     intra prediction.*/
  *_nxblocks=_ti->pic_width-(4*B_SZ<<xshift)>>B_SZ_LOG+xshift;
  *_nyblocks=_ti->pic_height-(4*B_SZ<<yshift)>>B_SZ_LOG+yshift;
}

char *get_map_filename(const char *_name,int _pli,int _nxblocks,int _nyblocks){
  char  *ret;
  char   fname[8192];
  size_t fname_len;
  sprintf(fname,"%s-pli%i-%i-%ix%i.intra.map",
   _name,_pli,B_SZ,_nxblocks,_nyblocks);
  fname_len=strlen(fname);
  ret=(char *)malloc(fname_len+1);
  memcpy(ret,fname,fname_len+1);
  return ret;
}


int apply_to_blocks(void *_ctx,plane_start_func _start,block_func _block,
 plane_finish_func _finish,int _argc,const char **_argv){
  int ai;
  for(ai=1;ai<_argc;ai++){
    video_input      vid;
    th_info          ti;
    th_ycbcr_buffer  ycbcr;
    FILE            *fin;
    int              pli;
    fin=fopen(_argv[ai],"rb");
    if(fin==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    if(video_input_open(&vid,fin)<0){
      fprintf(stderr,"Error reading video info from '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    video_input_get_info(&vid,&ti);
    if(video_input_fetch_frame(&vid,ycbcr,NULL)<0){
      fprintf(stderr,"Error reading first frame from '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    for(pli=0;pli<3;pli++){
      const unsigned char *data;
      int                  x0;
      int                  y0;
      int                  nxblocks;
      int                  nyblocks;
      int                  ret;
      int                  stride;
      int                  bi;
      int                  bj;
      get_intra_dims(&ti,pli,&x0,&y0,&nxblocks,&nyblocks);
      ret=(*_start)(_ctx,_argv[ai],&ti,pli,nxblocks,nyblocks);
      if(ret)return ret;
      data=ycbcr[pli].data;
      stride=ycbcr[pli].stride;
      for(bj=0;bj<nyblocks;bj++){
        for(bi=0;bi<nxblocks;bi++){
          (*_block)(_ctx,data+(y0+bj*B_SZ)*stride+x0+bi*B_SZ,stride,bi,bj);
        }
      }
      ret=(*_finish)(_ctx);
      if(ret)return ret;
    }
    video_input_close(&vid);
  }
  return EXIT_SUCCESS;
}
