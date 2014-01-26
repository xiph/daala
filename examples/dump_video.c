/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

/*This code was cargo-culted from libtheora's dump_video.c.*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if !defined(_REENTRANT)
#define _REENTRANT
#endif
#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#if !defined(_LARGEFILE_SOURCE)
#define _LARGEFILE_SOURCE
#endif
#if !defined(_LARGEFILE64_SOURCE)
#define _LARGEFILE64_SOURCE
#endif
#if !defined(_FILE_OFFSET_BITS)
#define _FILE_OFFSET_BITS 64
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*Yes, yes, we're going to hell.*/
#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#endif
#include <limits.h>
#include <math.h>
#include <signal.h>
#include "getopt.h"
#include "../src/logging.h"
#include "../include/daala/daaladec.h"

#define OD_SUPERBLOCK_SIZE (32)

const char *optstring = "o:rf";
struct option options [] = {
  {"output",required_argument,NULL,'o'},
  {"raw",no_argument, NULL,'r'}, /*Disable YUV4MPEG2 headers:*/
  {NULL,0,NULL,0}
};

/* Helper; just grab some more compressed bitstream and sync it for
   page extraction */
int buffer_data(FILE *in,ogg_sync_state *oy){
  char *buffer=ogg_sync_buffer(oy,4096);
  int bytes=fread(buffer,1,4096,in);
  ogg_sync_wrote(oy,bytes);
  return(bytes);
}

ogg_sync_state      oy;
ogg_page            og;
ogg_stream_state    vo;
ogg_stream_state    to;
daala_info          di;
daala_comment       dc;
daala_setup_info   *ds;
daala_dec_ctx      *dd;
od_img              img;

int              daala_p=0;
int              daala_processing_headers;
int              stateflag=0;

/* single frame video buffering */
int          videobuf_ready=0;
int          raw=0;

FILE* outfile = NULL;

int got_sigint=0;
static void sigint_handler (int signal) {
  (void)signal;
  got_sigint = 1;
}

/*Write out the planar YUV frame, uncropped.*/
static void video_write(void){
  int pli;
  int i;
  if(outfile){
    if(!raw)fprintf(outfile, "FRAME\n");
    for(pli=0;pli<img.nplanes;pli++){
      int plane_width;
      int plane_height;
      int xdec;
      int ydec;
      xdec=img.planes[pli].xdec;
      ydec=img.planes[pli].ydec;
      plane_width=(img.width+(1<<xdec)-1)>>xdec;
      plane_height=(img.height+(1<<ydec)-1)>>ydec;
      for(i=0;i<plane_height;i++){
        if(fwrite(img.planes[pli].data+img.planes[pli].ystride*i, 1,
         plane_width, outfile) < (size_t)plane_width){
          fprintf(stderr, "Error writing yuv frame");
          return;
        }
      }
    }
  }
}

/* dump the daala comment header */
static int dump_comments(daala_comment *_dc){
  int   i;
  int   len;
  FILE *out;
  out=stderr;
  fprintf(out,"Encoded by %s\n",_dc->vendor);
  if(_dc->comments){
    fprintf(out,"daala comment header:\n");
    for(i=0;i<_dc->comments;i++){
      if(_dc->user_comments[i]){
        len=_dc->comment_lengths[i]<INT_MAX?_dc->comment_lengths[i]:INT_MAX;
        fprintf(out,"\t%.*s\n",len,_dc->user_comments[i]);
      }
    }
  }
  return 0;
}

/* helper: push a page into the appropriate steam */
/* this can be done blindly; a stream won't accept a page
                that doesn't belong to it */
static int queue_page(ogg_page *page){
  if(daala_p)ogg_stream_pagein(&to,page);
  return 0;
}

static void usage(void){
  fprintf(stderr,
   "Usage: dumpvid [options] [<file.ogv>] [-o outfile.y4m]\n"
   "If no input file is given, stdin is used.\n"
   "Options:\n\n"
   "  -o --output <outfile.y4m> File name for decoded output. If\n"
   "                            this option is not given, the\n"
   "                            decompressed data is sent to stdout.\n"
   "  -r --raw                  Output raw YUV with no framing instead\n"
   "                            of YUV4MPEG2 (the default).\n");
  exit(EXIT_FAILURE);
}

int main(int argc,char *argv[]){

  ogg_packet op;

  int long_option_index;
  int c;

  int frames = 0;
  int pix_fmt = 1;
  ogg_int32_t pic_width = 0;
  ogg_int32_t pic_height = 0;
  ogg_int32_t fps_num = 0;
  ogg_int32_t fps_denom = 0;
  FILE *infile = stdin;
  outfile = stdout;
  daala_log_init();
#ifdef _WIN32 /* We need to set stdin/stdout to binary mode on windows. */
  /* Beware the evil ifdef. We avoid these where we can, but this one we
     cannot. Don't add any more, you'll probably go to hell if you do. */
  _setmode( _fileno( stdin ), _O_BINARY );
  _setmode( _fileno( stdout ), _O_BINARY );
#endif

  /* Process option arguments. */
  while((c=getopt_long(argc,argv,optstring,options,&long_option_index))!=EOF){
    switch(c){
    case 'o':
      if(strcmp(optarg,"-")!=0){
        outfile=fopen(optarg,"wb");
        if(outfile==NULL){
          fprintf(stderr,"Unable to open output file '%s'\n", optarg);
          exit(1);
        }
      }else{
        outfile=stdout;
      }
      break;

    case 'r':
      raw=1;
      break;

    default:
      usage();
    }
  }
  if(optind<argc){
    infile=fopen(argv[optind],"rb");
    if(infile==NULL){
      fprintf(stderr,"Unable to open '%s' for extraction.\n", argv[optind]);
      exit(1);
    }
    if(++optind<argc){
      usage();
      exit(1);
    }
  }
  /*Ok, Ogg parsing.
    The idea here is we have a bitstream that is made up of Ogg pages.
    The libogg sync layer will find them for us.
    There may be pages from several logical streams interleaved; we find the
     first daala stream and ignore any others.
    Then we pass the pages for our stream to the libogg stream layer which
     assembles our original set of packets out of them.
   start up Ogg stream synchronization layer */
  ogg_sync_init(&oy);

  /* init supporting Theora structures needed in header parsing */
  daala_comment_init(&dc);
  daala_info_init(&di);

  /*Ogg file open; parse the headers.
    Theora (like Vorbis) depends on some initial header packets for decoder
     setup and initialization.
    We retrieve these first before entering the main decode loop.*/

  /* Only interested in Daala streams */
  while(!stateflag){
    int ret=buffer_data(infile,&oy);
    if(ret==0)break;
    while(ogg_sync_pageout(&oy,&og)>0){
      int got_packet;
      ogg_stream_state test;

      /* is this a mandated initial header? If not, stop parsing */
      if(!ogg_page_bos(&og)){
        /* don't leak the page; get it into the appropriate stream */
        queue_page(&og);
        stateflag=1;
        break;
      }

      ogg_stream_init(&test,ogg_page_serialno(&og));
      ogg_stream_pagein(&test,&og);
      got_packet = ogg_stream_packetpeek(&test,&op);

      /* identify the codec: try daala */
      if((got_packet==1) && !daala_p && (daala_processing_headers=
       daala_decode_header_in(&di,&dc,&ds,&op))>=0){
        /* it is daala -- save this stream state */
        memcpy(&to,&test,sizeof(test));
        daala_p=1;
        /*Advance past the successfully processed header.*/
        if(daala_processing_headers)ogg_stream_packetout(&to,NULL);
      }else{
        /* whatever it is, we don't care about it */
        ogg_stream_clear(&test);
      }
    }
    /* fall through to non-bos page parsing */
  }

  /* we're expecting more header packets. */
  while(daala_p && daala_processing_headers){
    int ret;

    /* look for further daala headers */
    while(daala_processing_headers&&(ret=ogg_stream_packetpeek(&to,&op))){
      if(ret<0)continue;
      daala_processing_headers=daala_decode_header_in(&di,&dc,&ds,&op);
      if(daala_processing_headers<0){
        fprintf(stderr,"Error parsing Daala stream headers; "
         "corrupt stream?\n");
        exit(1);
      }
      else if(daala_processing_headers>0){
        /*Advance past the successfully processed header.*/
        ogg_stream_packetout(&to,NULL);
      }
      daala_p++;
    }

    /*Stop now so we don't fail if there aren't enough pages in a short
       stream.*/
    if(!(daala_p && daala_processing_headers))break;

    /* The header pages/packets will arrive before anything else we
       care about, or the stream is not obeying spec */

    if(ogg_sync_pageout(&oy,&og)>0){
      queue_page(&og); /* demux into the appropriate stream */
    }else{
      int ret=buffer_data(infile,&oy); /* someone needs more data */
      if(ret==0){
        fprintf(stderr,"End of file while searching for codec headers.\n");
        exit(1);
      }
    }
  }

  /* and now we have it all.  initialize decoders */
  if(daala_p){
    dump_comments(&dc);
    dd=daala_decode_alloc(&di,ds);
    fprintf(stderr,"Ogg logical stream %lx is Daala %dx%d %.02f fps video\n",
     to.serialno,di.pic_width,di.pic_height,
     di.timebase_numerator/(double)di.timebase_denominator*di.frame_duration);
  }else{
    /* tear down the partial daala setup */
    daala_info_clear(&di);
    daala_comment_clear(&dc);
  }

  /*Either way, we're done with the codec setup data.*/
  daala_setup_free(ds);


  if(!raw && outfile){
    static const char *CHROMA_TYPES[4]={"420jpeg",NULL,"422","444"};
    pic_width=di.pic_width;
    pic_height=di.pic_height;
    fps_num=di.timebase_numerator;
    fps_denom=di.timebase_denominator*di.frame_duration;
    /*calculate pixel_fmt based on the xdec & ydec values from one of the
      chroma planes.*/
    if(di.plane_info[1].xdec==1 && di.plane_info[1].ydec==1) {
     pix_fmt = 0;
    }else if(di.plane_info[1].xdec==1 && di.plane_info[1].ydec==0) {
      pix_fmt = 2;
    }else if(di.plane_info[1].xdec==0 && di.plane_info[1].ydec==0) {
      pix_fmt = 3;
    }
    if(pix_fmt>=4||pix_fmt==1){
      fprintf(stderr,"Unknown pixel format: %i\n",pix_fmt);
      exit(1);
    }
    /*Store header information*/
    fprintf(outfile,"YUV4MPEG2 C%s W%d H%d F%d:%d A%d:%d\n",
     CHROMA_TYPES[pix_fmt],pic_width,pic_height,fps_num,fps_denom,
     di.pixel_aspect_numerator,di.pixel_aspect_denominator);
  }

  /* install signal handler */
  signal (SIGINT, sigint_handler);

  /*Finally the main decode loop.

    It's one Daala packet per frame, so this is pretty straightforward if
     we're not trying to maintain sync with other multiplexed streams.

    The videobuf_ready flag is used to maintain the input buffer in the libogg
     stream state.
    If there's no output frame available at the end of the decode step, we must
     need more input data.
    We could simplify this by just using the return code on
     ogg_page_packetout(), but the flag system extends easily to the case where
     you care about more than one multiplexed stream (like with audio
     playback).
    In that case, just maintain a flag for each decoder you care about, and
     pull data when any one of them stalls.*/

  stateflag=0; /* playback has not begun */
  /* queue any remaining pages from data we buffered but that did not
      contain headers */
  while(ogg_sync_pageout(&oy,&og)>0){
    queue_page(&og);
  }

  while(!got_sigint){

    while(daala_p && !videobuf_ready){
      if(ogg_stream_packetout(&to,&op)>0){

        if(daala_decode_packet_in(dd,&img,&op)>=0){
          videobuf_ready=1;
          frames++;
        }

      }else
        break;
    }


    if(!videobuf_ready && feof(infile))break;

    if(!videobuf_ready){
      /* no data yet for somebody.  Grab another page */
      buffer_data(infile,&oy);
      while(ogg_sync_pageout(&oy,&og)>0){
        queue_page(&og);
      }
    }
    /* dumpvideo frame, and get new one */
    else if(outfile)video_write();

    videobuf_ready=0;
  }

  /* end of decoder loop -- close everything */

  if(daala_p){
    ogg_stream_clear(&to);
    daala_decode_free(dd);
    daala_comment_clear(&dc);
    /*TODO: Uncomment once implemetned
      daala_info_clear(&ti);*/
  }
  ogg_sync_clear(&oy);

  if(infile && infile!=stdin)fclose(infile);
  if(outfile && outfile!=stdout)fclose(outfile);

  fprintf(stderr, "\n\n%d frames\n", frames);
  fprintf(stderr, "\nDone.\n");

  return(0);

}
