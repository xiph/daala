#if !defined(_vidinput_H)
# define _vidinput_H (1)
# if !defined(_LARGEFILE_SOURCE)
#  define _LARGEFILE_SOURCE
# endif
# if !defined(_LARGEFILE64_SOURCE)
#  define _LARGEFILE64_SOURCE
# endif
# if !defined(_FILE_OFFSET_BITS)
#  define _FILE_OFFSET_BITS 64
# endif
# include <stdio.h>
# include "theora/theoradec.h"

# if defined(__cplusplus)
extern "C" {
# endif



typedef struct video_input      video_input;
typedef struct video_input_vtbl video_input_vtbl;

typedef void (*video_input_get_info_func)(void *_ctx,th_info *_ti);
typedef int (*video_input_fetch_frame_func)(void *_ctx,FILE *_fin,
 th_ycbcr_buffer _ycbcr,char _tag[5]);
typedef void (*video_input_close_func)(void *_ctx);



struct video_input_vtbl{
  video_input_get_info_func     get_info;
  video_input_fetch_frame_func  fetch_frame;
  video_input_close_func        close;
};

struct video_input{
  const video_input_vtbl *vtbl;
  void                   *ctx;
  FILE                   *fin;
};



int video_input_open(video_input *_vid,FILE *_fin);
void video_input_close(video_input *_vid);

void video_input_get_info(video_input *_vid,th_info *_ti);
int video_input_fetch_frame(video_input *_vid,
 th_ycbcr_buffer _ycbcr,char _tag[5]);

# if defined(__cplusplus)
}
# endif

#endif
