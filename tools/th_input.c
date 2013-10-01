/*Daala video codec
Copyright (c) 2002-2013 Daala project contributors.  All rights reserved.

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

#include "vidinput.h"
#include <stdlib.h>
#include <string.h>
#include <ogg/ogg.h>
#include <theora/theoradec.h>

typedef struct th_input th_input;

struct th_input{
  ogg_sync_state    oy;
  int               theora_p;
  ogg_stream_state  to;
  th_info           ti;
  th_comment        tc;
  th_dec_ctx       *td;
};

/*Grab some more compressed bitstream and sync it for page extraction.*/
static int th_input_buffer_data(th_input *_th,FILE *_fin){
  char *buffer;
  int bytes;
  buffer=ogg_sync_buffer(&_th->oy,4096);
  bytes=fread(buffer,1,4096,_fin);
  ogg_sync_wrote(&_th->oy,bytes);
  return bytes;
}

/*Push a page into the appropriate steam.
  This can be done blindly; a stream won't accept a page that doesn't belong to
   it.*/
static void th_input_queue_page(th_input *_th,ogg_page *_og){
  if(_th->theora_p)ogg_stream_pagein(&_th->to,_og);
}

static int th_input_open_impl(th_input *_th,th_setup_info **_ts,FILE *_fin){
  ogg_packet op;
  ogg_page   og;
  int        nheaders_left;
  int        done_headers;
  ogg_sync_init(&_th->oy);
  th_info_init(&_th->ti);
  th_comment_init(&_th->tc);
  *_ts=NULL;
  _th->theora_p=0;
  nheaders_left=0;
  for(done_headers=0;!done_headers;){
    if(th_input_buffer_data(_th,_fin)==0)break;
    while(ogg_sync_pageout(&_th->oy,&og)>0){
      ogg_stream_state test;
      /*Is this a mandated initial header?
        If not, stop parsing.*/
      if(!ogg_page_bos(&og)){
        /*Don't leak the page; get it into the appropriate stream.*/
        th_input_queue_page(_th,&og);
        done_headers=1;
        break;
      }
      ogg_stream_init(&test,ogg_page_serialno(&og));
      ogg_stream_pagein(&test,&og);
      ogg_stream_packetpeek(&test,&op);
      /*Identify the codec: try Theora.*/
      if(!_th->theora_p){
        nheaders_left=th_decode_headerin(&_th->ti,&_th->tc,_ts,&op);
        if(nheaders_left>=0){
          /*It is Theora.*/
          memcpy(&_th->to,&test,sizeof(test));
          _th->theora_p=1;
          /*Advance past the successfully processed header.*/
          if(nheaders_left>0)ogg_stream_packetout(&_th->to,NULL);
          continue;
        }
      }
      /*Whatever it is, we don't care about it.*/
      ogg_stream_clear(&test);
    }
  }
  /*We're expecting more header packets.*/
  while(_th->theora_p&&nheaders_left>0){
    int ret;
    while(nheaders_left>0){
      ret=ogg_stream_packetpeek(&_th->to,&op);
      if(ret==0)break;
      if(ret<0)continue;
      nheaders_left=th_decode_headerin(&_th->ti,&_th->tc,_ts,&op);
      if(nheaders_left<0){
        fprintf(stderr,"Error parsing Theora stream headers; "
         "corrupt stream?\n");
        return -1;
      }
      /*Advance past the successfully processed header.*/
      else if(nheaders_left>0)ogg_stream_packetout(&_th->to,NULL);
      _th->theora_p++;
    }
    /*Stop now so we don't fail if there aren't enough pages in a short
       stream.*/
    if(!(_th->theora_p&&nheaders_left>0))break;
    /*The header pages/packets will arrive before anything else we care
       about, or the stream is not obeying spec.*/
    if(ogg_sync_pageout(&_th->oy,&og)>0)th_input_queue_page(_th,&og);
    /*We need more data.*/
    else if(th_input_buffer_data(_th,_fin)==0){
      fprintf(stderr,"End of file while searching for codec headers.\n");
      return -1;
    }
  }
  /*And now we have it all.
    Initialize the decoder.*/
  if(_th->theora_p){
    _th->td=th_decode_alloc(&_th->ti,*_ts);
    if(_th->td!=NULL){
      fprintf(stderr,"Ogg logical stream %lx is Theora %ix%i %.02f fps video.\n"
       "Encoded frame content is %ix%i with %ix%i offset.\n",
       _th->to.serialno,_th->ti.frame_width,_th->ti.frame_height,
       (double)_th->ti.fps_numerator/_th->ti.fps_denominator,
       _th->ti.pic_width,_th->ti.pic_height,_th->ti.pic_x,_th->ti.pic_y);
      return 1;
    }
  }
  return -1;
}

static void th_input_close(th_input *_th){
  if(_th->theora_p){
    ogg_stream_clear(&_th->to);
    th_decode_free(_th->td);
  }
  th_comment_clear(&_th->tc);
  th_info_clear(&_th->ti);
  ogg_sync_clear(&_th->oy);
}

static int th_input_open(th_input *_th,FILE *_fin){
  th_input       th;
  th_setup_info *ts;
  int            ret;
  ret=th_input_open_impl(&th,&ts,_fin);
  th_setup_free(ts);
  /*Clean up on failure.*/
  if(ret<0)th_input_close(&th);
  else memcpy(_th,&th,sizeof(th));
  return ret;
}

static void th_input_get_info(th_input *_th,th_info *_ti){
  memcpy(_ti,&_th->ti,sizeof(*_ti));
}

static int th_input_fetch_frame(th_input *_th,FILE *_fin,
 video_input_ycbcr _ycbcr,char _tag[5]){
  for(;;){
    ogg_page   og;
    ogg_packet op;
    int        ret;
    if(ogg_stream_packetout(&_th->to,&op)>0){
      ret=th_decode_packetin(_th->td,&op,NULL);
      if(ret>=0){
        th_ycbcr_buffer buffer;
        th_decode_ycbcr_out(_th->td,buffer);
        memcpy(_ycbcr,buffer,sizeof(*_ycbcr)*3);
        if(_tag!=NULL){
          _tag[0]=th_packet_iskeyframe(&op)?'K':'D';
          if(op.bytes>0)sprintf(_tag+1,"%02i ",op.packet[0]&0x3F);
          else strcpy(_tag+1,"-- ");
        }
        return ret+1;
      }
      else return -1;
    }
    while(ogg_sync_pageout(&_th->oy,&og)<=0){
      if(th_input_buffer_data(_th,_fin)==0)return feof(_fin)?0:-1;
    }
    th_input_queue_page(_th,&og);
  }
}

const video_input_vtbl TH_INPUT_VTBL={
  (video_input_open_func)th_input_open,
  (video_input_get_info_func)th_input_get_info,
  (video_input_fetch_frame_func)th_input_fetch_frame,
  (video_input_close_func)th_input_close
};
