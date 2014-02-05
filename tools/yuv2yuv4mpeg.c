#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*This is a small utility used to convert many of the standard test sequences
   into a common format, the YUV4MPEG format used by mjpegtools.
  This is the format that the example encoder takes as input.

  The set of file name patterns/extentions it supports is intentionally very
   limited to keep the code simple.
  However, it is not too difficult to rename existing file names to fit these
   patterns.
  E.g., to add leading 0's to 1- and 2-digit frame numbers with /bin/sh,
    for i in tennis?.yuv ; do mv "$i" "tennis00${i##tennis}" ; done
    for i in tennis??.yuv ; do mv "$i" "tennis0${i##tennis}" ; done*/
int main(int _argc,char *_argv[]){
  FILE          *in_y;
  FILE          *in_u;
  FILE          *in_v;
  FILE          *in_a;
  FILE          *in_f;
  FILE          *out_y4m;
  int            y_w;
  int            y_h;
  int            c_w;
  int            c_h;
  char           ip;
  int            fps_n;
  int            fps_d;
  int            par_n;
  int            par_d;
  char          *chroma;
  int            have_chroma;
  int            have_alpha;
  int            frame_digits;
  int            separate_planes;
  unsigned char *buf;
  char           fname[8192];
  int            frame_num;
  int            argi;
  if(_argc<2){
    fprintf(stderr,
     "Usage: yuv2yuv4mpeg <sequence_qcif> [-p | -d<digits>] [<options>]\n"
     "Converts one or more raw data files in a sequence into a single\n"
     " YUV4MPEG file with appropriate headers to describe the video format.\n"
     "\n"
     "The type of input read is controlled by one or more of:\n"
     " -p              Input is split by component.\n"
     "                 Opens up three files, <sequence>.y, <sequence>.u,\n"
     "                  and <sequence>.v containing raw, unpadded image\n"
     "                  data for each component, respectively.\n"
     "                 If only a Y' component is present, <sequence>.y is\n"
     "                  used.\n"
     "                 If an alpha component is present, <sequence>.a is\n"
     "                  used.\n"
     "                 Which components are present is determined from the\n"
     "                  -c flag described below.\n"
     " -d<digits>      Input is split by frame.\n"
     "                 Opens up <sequence><ddd>.yuv for frame <ddd>\n"
     "                  containing raw, unpadded image data for all three\n"
     "                  components.\n"
     "                 <digits> specifies the number of decimal digits in\n"
     "                  the frame number <ddd>.\n"
     " If both of these options are be specified, the three files\n"
     "  <sequence><ddd>.y, <sequence><ddd>.u, and <sequence><ddd>.v are\n"
     "  opened for each frame.\n"
     " If neither of these options are specified, a single <sequence>.yuv\n"
     "  file is opened, containing, the raw, unpadded image data for all\n"
     "  components of all frames.\n"
     " This must be organized in a planar fashion, with each complete frame\n"
     "  in one block, in the order Y, U (Cb), V (Cr), alpha.\n"
     " Which components are present is determined from the -c flag described\n"
     "  below.\n"
     "\n"
     "The default input video format is 29.97 fps progressive QCIF.\n"
     "Options:\n"
     " -w<width>\n"
     " -h<height>\n"
     " -fn<frame rate numerator>\n"
     " -fd<frame rate denominator>\n"
     " -i<p|t|b>       Progressive or Top-field first or Bottom-field first.\n"
     " -c<chroma>      Chroma subsampling format. One of:\n"
     "                   420jpeg\n"
     "                   420mpeg2\n"
     "                   420paldv\n"
     "                   411\n"
     "                   422\n"
     "                   444\n"
     "                   444alpha\n"
     "                   mono\n"
     " -an<pixel aspect numerator>\n"
     " -ad<pixel aspect denominator>\n"
     " Note: Most common frame sizes do not use square pixels:\n"
     "  Name  Width  Height  -an   -ad\n"
     "        720    576     128   117\n"
     "  CIF   352    288     128   117\n"
     "  QCIF  176    144     128   117\n"
     "        720    480     4320  4739\n"
     "  SIF   352    240     4320  4739\n"
     " The HD TV sizes (1920x1080, 1280x720) really do use square pixels.\n"
     " -o<offset>      The first frame number to use in -d mode.\n"
     "                 If unspecified, 0 is used.\n");
    return -1;
  }
  frame_digits=frame_num=separate_planes=0;
  in_y=in_u=in_v=in_a=in_f=NULL;
  y_w=176;
  y_h=144;
  fps_n=30000;
  fps_d=1001;
  ip='p';
  par_n=128;
  par_d=117;
  chroma="420jpeg";
  for(argi=2;argi<_argc;argi++){
    if(_argv[argi][0]!='-'){
      fprintf(stderr,"Error parsing arguments: must start with '-'.\n");
      return -1;
    }
    switch(_argv[argi][1]){
      case 'p':separate_planes=1;break;
      case 'd':frame_digits=atoi(_argv[argi]+2);break;
      case 'w':y_w=atoi(_argv[argi]+2);break;
      case 'h':y_h=atoi(_argv[argi]+2);break;
      case 'i':ip=_argv[argi][2];break;
      case 'c':chroma=_argv[argi]+2;break;
      case 'f':{
        if(_argv[argi][2]=='n'){
          fps_n=atoi(_argv[argi]+3);
          break;
        }
        else if(_argv[argi][2]=='d'){
          fps_d=atoi(_argv[argi]+3);
          break;
        }
      }
      case 'a':{
        if(_argv[argi][2]=='n'){
          par_n=atoi(_argv[argi]+3);
          break;
        }
        else if(_argv[argi][2]=='d'){
          par_d=atoi(_argv[argi]+3);
          break;
        }
      }
      case 'o':frame_num=atoi(_argv[argi]+2);break;
      default:{
        fprintf(stderr,"Error parsing arguments: unknown switch '%c'.\n",
         _argv[argi][1]);
        return -1;
      }
    }
  }
  if(strncmp(chroma,"420",3)==0){
    c_w=y_w+1>>1;
    c_h=y_h+1>>1;
    have_chroma=1;
  }
  else if(strncmp(chroma,"411",3)==0){
    c_w=y_w+3>>2;
    c_h=y_h;
    have_chroma=1;
  }
  else if(strcmp(chroma,"422")==0){
    c_w=y_w+1>>1;
    c_h=y_h;
    have_chroma=1;
  }
  else if(strncmp(chroma,"444",3)==0){
    c_w=y_w;
    c_h=y_h;
    have_chroma=1;
  }
  else if(strncmp(chroma,"mono",4)==0){
    c_w=c_h=0;
    have_chroma=0;
  }
  else{
    fprintf(stderr,"Error parsing arguments: Unsupported chroma mode: %s.\n",
     chroma);
    return -1;
  }
  have_alpha=strstr(chroma,"alpha")!=NULL;
  if(frame_digits<0)frame_digits=0;
  if(strlen(_argv[1])>8178-frame_digits){
    fprintf(stderr,"File name too long.\n");
    return -1;
  }
  /*Open input and output files.*/
  if(separate_planes){
    if(frame_digits>0){
      sprintf(fname,"%s%0*d.y",_argv[1],frame_digits,frame_num);
      in_y=fopen(fname,"rb");
      if(have_chroma){
         sprintf(fname,"%s%0*d.u",_argv[1],frame_digits,frame_num);
         in_u=fopen(fname,"rb");
         sprintf(fname,"%s%0*d.v",_argv[1],frame_digits,frame_num);
         in_v=fopen(fname,"rb");
      }
      if(have_alpha){
         sprintf(fname,"%s%0*d.a",_argv[1],frame_digits,frame_num);
         in_a=fopen(fname,"rb");
      }
    }
    else{
      sprintf(fname,"%s.y",_argv[1]);
      in_y=fopen(fname,"rb");
      if(have_chroma){
        sprintf(fname,"%s.u",_argv[1]);
        in_u=fopen(fname,"rb");
        sprintf(fname,"%s.v",_argv[1]);
        in_v=fopen(fname,"rb");
      }
      if(have_alpha){
        sprintf(fname,"%s.a",_argv[1]);
        in_a=fopen(fname,"rb");
      }
    }
  }
  else{
    if(frame_digits>0){
      sprintf(fname,"%s%0*d.yuv",_argv[1],frame_digits,frame_num);
      in_f=fopen(fname,"rb");
    }
    else{
      sprintf(fname,"%s.yuv",_argv[1]);
      in_f=fopen(fname,"rb");
    }
    in_y=in_f;
    if(have_chroma)in_u=in_v=in_f;
    if(have_alpha)in_a=in_f;
  }
  if(in_y==NULL||have_chroma&&(in_u==NULL||in_v==NULL)||
   have_alpha&&in_a==NULL){
    fprintf(stderr,"Error opening input file(s).\n");
    return -1;
  }
  sprintf(fname,"%s.y4m",_argv[1]);
  out_y4m=fopen(fname,"wb");
  if(out_y4m==NULL){
    fprintf(stderr,"Error opening output file.\n");
    return -1;
  }
  /*Start output.*/
  fprintf(out_y4m,"YUV4MPEG2 W%i H%i F%i:%i I%c A%i:%i",
   y_w,y_h,fps_n,fps_d,ip,par_n,par_d);
  if(strcmp(chroma,"420jpeg")!=0)fprintf(out_y4m," C%s",chroma);
  fprintf(out_y4m,"\n");
  buf=(unsigned char *)malloc((size_t)(y_w*y_h));
  for(;;){
    if(fread(buf,y_w,y_h,in_y)<y_h)break;
    fprintf(out_y4m,"FRAME\n");
    fwrite(buf,y_w,y_h,out_y4m);
    if(have_chroma){
      fread(buf,c_w,c_h,in_u);
      fwrite(buf,c_w,c_h,out_y4m);
      fread(buf,c_w,c_h,in_v);
      fwrite(buf,c_w,c_h,out_y4m);
    }
    if(have_alpha){
      fread(buf,y_w,y_h,in_a);
      fwrite(buf,y_w,y_h,out_y4m);
    }
    if(frame_digits>0){
      frame_num++;
      if(separate_planes){
        sprintf(fname,"%s%0*d.y",_argv[1],frame_digits,frame_num);
        fclose(in_y);
        in_y=fopen(fname,"rb");
        if(have_chroma){
          sprintf(fname,"%s%0*d.u",_argv[1],frame_digits,frame_num);
          fclose(in_u);
          in_u=fopen(fname,"rb");
          sprintf(fname,"%s%0*d.v",_argv[1],frame_digits,frame_num);
          fclose(in_v);
          in_v=fopen(fname,"rb");
        }
        if(have_alpha){
          sprintf(fname,"%s%0*d.a",_argv[1],frame_digits,frame_num);
          fclose(in_a);
          in_a=fopen(fname,"rb");
        }
      }
      else{
        sprintf(fname,"%s%0*d.yuv",_argv[1],frame_digits,frame_num);
        fclose(in_f);
        in_f=fopen(fname,"rb");
        in_y=in_f;
        if(have_chroma)in_u=in_v=in_f;
        if(have_alpha)in_a=in_f;
      }
      if(in_y==NULL||have_chroma&&(in_u==NULL||in_v==NULL)||
       have_alpha&&in_a==NULL){
        break;
      }
    }
  }
  return 0;
}
