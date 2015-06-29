#include <stdio.h>
#include <stdlib.h>

#include "../src/internal.h"
#include "../src/zigzag.h"

#define HEADER \
 "<html>\n" \
 "  <head>\n" \
 "    <style type=\"text/css\">\n" \
 "      svg {\n" \
 "        width: %i;\n" \
 "        height: %i;\n" \
 "      }\n" \
 "      polygon {\n" \
 "        stroke: grey;\n" \
 "        stroke-width: 1;\n" \
 "      }\n" \
 "      polyline {\n" \
 "        stroke: grey;\n" \
 "        stroke-width: 1;\n" \
 "      }\n" \
 "    </style>\n" \
 "  </head>\n" \
 "  <body>\n" \
 "    <svg>\n"
#define FOOTER \
 "    </svg>\n" \
 "  </body>\n" \
 "</html>"

typedef unsigned short od_color[4];

#define HUE_MAX ((6*0x10000) - 1)

static void rgb_from_hue(od_color color, int hue) {
  int i;
  int y;
  hue %= HUE_MAX;
  if (hue < 0) hue += HUE_MAX;
  i = hue/0x10000;
  y = hue&0xFFFF;
  color[0] = color[1] = color[2] = color[3] = (unsigned short)0xFFFFU;
  switch (i) {
    case 0 : {
      color[1] = (unsigned short)y;
      color[2] = 0;
      break;
    }
    case 1 : {
      color[0] = (unsigned short)(0xFFFFU - y);
      color[2] = 0;
      break;
    }
    case 2 : {
      color[0] = 0;
      color[2] = (unsigned short)y;
      break;
    }
    case 3 : {
      color[0]=0;
      color[1] = (unsigned short)(0xFFFFU - y);
      break;
    }
    case 4 : {
      color[0] = (unsigned short)y;
      color[1] = 0;
      break;
    }
    default : {
      color[1] = 0;
      color[2] = (unsigned short)(0xFFFFU - y);
      break;
    }
  }
}

static void mode_colors(od_color *colors, int nmodes) {
  int i;
  for (i = 0; i < nmodes; i++) {
    rgb_from_hue(colors[i], HUE_MAX*i/nmodes);
  }
}

#define C_SZ (15)

typedef struct od_point od_point;

struct od_point {
  int x;
  int y;
};

static void draw_polygon(FILE *fp, od_point pts[], int n, int bs, od_color c) {
  int i;
  fprintf(fp,"      <polygon points=\"");
  for (i = 0; i < n; i++) {
    fprintf(fp, "%s%i,%i", i > 0 ? " " : "", (pts[i].x + (1 << bs))*C_SZ,
     (pts[i].y + (1 << bs))*C_SZ);
  }
  fprintf(fp,"\" fill=\"#%02x%02x%02x\"/>\n",c[0] >> 8, c[1] >> 8, c[2] >> 8);
}

static void draw_polyline(FILE *fp, od_point pts[], int n, int bs, int o) {
  int i;
  fprintf(fp,"      <polyline points=\"");
  for (i = 0; i < n; i++) {
    fprintf(fp, "%s%i,%i", i > 0 ? " " : "",
     (pts[i].x + (1 << bs))*C_SZ + (o*C_SZ/2),
     (pts[i].y + (1 << bs))*C_SZ + (o*C_SZ/2));
  }
  fprintf(fp,"\"/>\n");
}

#define draw_rect(fp, x0, y0, x1, y1, bs, c) \
  do { \
    int i; \
    od_point pts[4]; \
    pts[0].x = x0; \
    pts[0].y = y0; \
    pts[1].x = x1; \
    pts[1].y = y0; \
    pts[2].x = x1; \
    pts[2].y = y1; \
    pts[3].x = x0; \
    pts[3].y = y1; \
    draw_polygon(fp, pts, 4, bs, c); \
    for (i = x0 + 1; i < x1; i++) { \
      pts[0].x = i; \
      pts[0].y = y0; \
      pts[1].x = i; \
      pts[1].y = y1; \
      draw_polyline(fp, pts, 2, bs, 0); \
    } \
    for (i = y0 + 1; i < y1; i++) { \
      pts[0].x = x0; \
      pts[0].y = i; \
      pts[1].x = x1; \
      pts[1].y = i; \
      draw_polyline(fp, pts, 2, bs, 0); \
    } \
  } \
  while (0)

#define draw_notch(fp, x0, y0, x1, bs, c) \
  do { \
    int i; \
    od_point pts[6]; \
    pts[0].x = x0; \
    pts[0].y = x0; \
    pts[1].x = x0; \
    pts[1].y = y0; \
    pts[2].x = x1; \
    pts[2].y = y0; \
    pts[3].x = x1; \
    pts[3].y = x1; \
    pts[4].x = y0; \
    pts[4].y = x1; \
    pts[5].x = y0; \
    pts[5].y = x0; \
    draw_polygon(fp, pts, 6, bs, c); \
    for (i = y0 + 1; i < x1; i++) { \
      pts[0].x = i; \
      pts[0].y = i < x0 ? x0 : y0; \
      pts[1].x = i; \
      pts[1].y = x1; \
      draw_polyline(fp, pts, 2, bs, 0); \
      pts[0].x = i < x0 ? x0 : y0; \
      pts[0].y = i; \
      pts[1].x = x1; \
      pts[1].y = i; \
      draw_polyline(fp, pts, 2, bs, 0); \
    } \
  } \
  while (0)

#define draw_path(fp, zz, n, bs) \
  do { \
    int i; \
    od_point pts[2]; \
    for (i = 1; i < n; i++) { \
      pts[0].x = (zz)[i - 1][0]; \
      pts[0].y = (zz)[i - 1][1]; \
      pts[1].x = (zz)[i][0]; \
      pts[1].y = (zz)[i][1]; \
      draw_polyline(fp, pts, 2, bs, 1); \
    } \
  } \
  while (0)

typedef unsigned char index_pair[2];

const index_pair *const OD_ZIGZAGS[] = {
  OD_ZIGZAG4,
  OD_ZIGZAG8,
  OD_ZIGZAG16,
  OD_ZIGZAG32,
  OD_ZIGZAG64
};

int main(void) {
  FILE *fp;
  int width;
  int bs;
  od_color colors[3*OD_NBSIZES - 1];
  int ncolors;
  fp = fopen("zigzags.html", "w");
  width = C_SZ*(OD_BSIZE_MAX + (1 << OD_NBSIZES - 1));
  fprintf(fp, HEADER, width, width);
  ncolors = 3*OD_NBSIZES - 1;
  mode_colors(colors, ncolors);
  for (bs = 0; bs < OD_NBSIZES; bs++) {
    int n;
    n = 1 << (bs + OD_LOG_BSIZE0);
    fprintf(fp,"      <!-- %ix%i -->\n",n,n);
    /* Special case 4x4 block */
    if (bs == 0) {
      draw_rect(fp, 0, 0, 1, 1, bs, colors[0]);
      draw_notch(fp, 1, 0, 4, bs, colors[1*3]);
      draw_path(fp, OD_ZIGZAGS[0], 15, bs);
    }
    else {
      draw_rect(fp, n/2, 0, n, n/4, bs, colors[(3*bs - 1)*3%ncolors]);
      draw_path(fp, OD_ZIGZAGS[bs], n*n/8, bs);
      draw_rect(fp, 0, n/2, n/4, n, bs, colors[(3*bs + 1)*3%ncolors]);
      draw_path(fp, OD_ZIGZAGS[bs] + n*n/8, n*n/8, bs);
      draw_notch(fp, n/2, n/4, n, bs, colors[3*bs*3%ncolors]);
      draw_path(fp, OD_ZIGZAGS[bs] + n*n/4, n*n/2, bs);
    }
  }
  fprintf(fp, FOOTER);
  fclose(fp);
  return EXIT_SUCCESS;
}
