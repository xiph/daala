#if !defined(_image_H)
# define _image_H (1)
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

typedef unsigned short od_rgba16_pixel[4];

typedef struct od_rgba16_image od_rgba16_image;



struct od_rgba16_image{
  od_rgba16_pixel *data;
  int              stride;
  int              width;
  int              height;
};


int od_rgba16_image_read_png(od_rgba16_image *_this,FILE *_fin);
void od_rgba16_image_init(od_rgba16_image *_this,int _width,int _height);

void od_rgba16_image_draw_point(od_rgba16_image *_this,int _x,int _y,
 const od_rgba16_pixel _color);
void od_rgba16_image_draw_line(od_rgba16_image *_this,
 int _x0,int _y0,int _x1,int _y1,const od_rgba16_pixel _color);

int od_rgba16_image_write_png(const od_rgba16_image *_this,FILE *_fout);

void od_rgba16_image_clear(od_rgba16_image *_this);

#endif
