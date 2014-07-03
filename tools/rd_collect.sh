#!/bin/bash
set -e

CODECS="<daala|vp8|vp9|x264|x265|libjpeg|mozjpeg|theora>"

if [ $# == 0 ]; then
  echo "usage: DAALA_ROOT=<build_dir> $0 $CODECS *.y4m"
  exit 1
fi

export CODEC=$1
shift

if [ $# == 0 ] || [[ ! $CODECS =~ [\<\|]$CODEC[\|\>] ]]; then
  echo "usage: DAALA_ROOT=<build_dir> $0 $CODECS *.y4m"
  exit 1
fi

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ ! -d $DAALA_ROOT ]; then
  echo "Please set DAALA_ROOT to the location of your libvpx git clone"
  exit 1
fi

case $CODEC in
  daala)
    if [ ! -f $DAALA_ROOT/config.h ]; then
      echo "File not found $DAALA_ROOT/config.h"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if ! grep -Fxq "#define OD_LOGGING_ENABLED 1" "$DAALA_ROOT/config.h"; then
      echo "Logging not enabled, re-run configure with --enable-logging"
      exit 1
    fi

    if ! grep -Fxq "#define OD_DUMP_IMAGES 1" "$DAALA_ROOT/config.h"; then
      echo "Image dumping not enabled, re-run configure without --disable-dump-images"
      exit 1
    fi

    if [ -z "$ENCODER_EXAMPLE" ]; then
      export ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
    fi

    if [ ! -x "$ENCODER_EXAMPLE" ]; then
      echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname "$0")/rd_collect_daala.sh
    ;;
  vp8 | vp9)
    if [ -z $LIBVPX_ROOT ] || [ ! -d $LIBVPX_ROOT ]; then
      echo "Please set LIBVPX_ROOT to the location of your libvpx git clone"
      exit 1
    fi

    if [ -z "$VPXENC" ]; then
      export VPXENC=$LIBVPX_ROOT/vpxenc
    fi

    if [ -z "$VPXDEC" ]; then
      export VPXDEC=$LIBVPX_ROOT/vpxdec
    fi

    if [ ! -x "$VPXENC" ]; then
      echo "Executable not found VPXENC=$VPXENC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    if [ ! -x "$VPXDEC" ]; then
      echo "Executable not found VPXDEC=$VPXDEC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_libvpx.sh
    ;;
  x264)
    if [ -z $X264_ROOT ] || [ ! -d $X264_ROOT ]; then
      echo "Please set X264_ROOT to the location of your x264 git clone"
      exit 1
    fi

    if [ -z "$X264" ]; then
      export X264=$X264_ROOT/x264
    fi

    if [ ! -x "$X264" ]; then
      echo "Executable not found X264=$X264"
      echo "Do you have the right X264_ROOT=$X264_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_x264.sh
    ;;
  x265)
    if [ -z $X265_ROOT ] || [ ! -d $X265_ROOT ]; then
      echo "Please set X265_ROOT to the location of your x265 hg checkout"
      exit 1
    fi

    if [ -z "$X265" ]; then
      export X265=$X265_ROOT/build/linux/x265
    fi

    if [ ! -x "$X265" ]; then
      echo "Executable not found X265=$X265"
      echo "Do you have the right X265_ROOT=$X265_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_x265.sh
    ;;
  libjpeg)
    if [ -z "$YUVJPEG" ]; then
      export YUVJPEG=$DAALA_ROOT/tools/yuvjpeg
    fi

    if [ -z "$JPEGYUV" ]; then
      export JPEGYUV=$DAALA_ROOT/tools/jpegyuv
    fi

    if [ ! -x "$YUVJPEG" ]; then
      echo "Executable not found YUVJPEG=$YUVJPEG"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if [ ! -x "$JPEGYUV" ]; then
      echo "Executable not found JPEGYUV=$JPEGYUV"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_jpeg.sh
    ;;
  mozjpeg)
    if [ -z $MOZJPEG_ROOT ] || [ ! -d $MOZJPEG_ROOT ]; then
      echo "Please set MOZJPEG_ROOT to the location of your mozjpeg git clone"
      exit 1
    fi

    if [ -z "$YUVJPEG" ]; then
      export YUVJPEG=$MOZJPEG_ROOT/yuvjpeg
    fi

    if [ -z "$JPEGYUV" ]; then
      export JPEGYUV=$MOZJPEG_ROOT/jpegyuv
    fi

    if [ ! -x "$YUVJPEG" ]; then
      echo "Executable not found YUVJPEG=$YUVJPEG"
      echo "Do you have the right MOZJPEG_ROOT=$MOZJPEG_ROOT"
      exit 1
    fi

    if [ ! -x "$JPEGYUV" ]; then
      echo "Executable not found JPEGYUV=$JPEGYUV"
      echo "Do you have the right MOZJPEG_ROOT=$MOZJPEG_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_jpeg.sh
    ;;
  theora)
    if [ -z $THEORA_ROOT ] || [ ! -d $THEORA_ROOT ]; then
      echo "Please set THEORA_ROOT to the location of your theora svn checkout"
      exit 1
    fi

    if [ -z "$ENCODER_EXAMPLE" ]; then
      export ENCODER_EXAMPLE=$THEORA_ROOT/examples/encoder_example
    fi

    if [ -z "$DUMP_VIDEO" ]; then
      export DUMP_VIDEO=$THEORA_ROOT/examples/dump_video
    fi

    if [ ! -x "$ENCODER_EXAMPLE" ]; then
      echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
      echo "Do you have the right THEORA_ROOT=$THEORA_ROOT"
      exit 1
    fi

    if [ ! -x "$DUMP_VIDEO" ]; then
      echo "Executable not found DUMP_VIDEO=$DUMP_VIDEO"
      echo "Do you have the right THEORA_ROOT=$THEORA_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_theora.sh
    ;;
  *)
    echo "Unknown codec: $CODEC"
    exit 1
esac

if [ -z "$PLANE" ]; then
  export PLANE=0
fi

if [ $PLANE != 0 ] && [ $PLANE != 1 ] && [ $PLANE != 2 ] &&
  [ $PLANE != -1 ]; then
  echo "Invalid plane $PLANE. Must be 0, 1, 2, or -1 (all planes)."
  exit 1
fi

# TODO refactor these out of the daala project into a metrics project

if [ -z "$YUV2YUV4MPEG" ]; then
  export YUV2YUV4MPEG=$DAALA_ROOT/tools/yuv2yuv4mpeg
fi

if [ -z "$DUMP_PSNR" ]; then
  export DUMP_PSNR=$DAALA_ROOT/tools/dump_psnr
fi

if [ -z "$DUMP_PSNRHVS" ]; then
  export DUMP_PSNRHVS=$DAALA_ROOT/tools/dump_psnrhvs
fi

if [ -z "$DUMP_SSIM" ]; then
  export DUMP_SSIM=$DAALA_ROOT/tools/dump_ssim
fi

if [ -z "$DUMP_FASTSSIM" ]; then
  export DUMP_FASTSSIM=$DAALA_ROOT/tools/dump_fastssim
fi

if [ ! -x "$YUV2YUV4MPEG" ]; then
  echo "Executable not found YUV2YUV4MPEG=$YUV2YUV4MPEG"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNR" ]; then
  echo "Executable not found DUMP_PSNR=$DUMP_PSNR"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNRHVS" ]; then
  echo "Executable not found DUMP_PSNRHVS=$DUMP_PSNRHVS"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_SSIM" ]; then
  echo "Executable not found DUMP_FASTSSIM=$DUMP_SSIM"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_FASTSSIM" ]; then
  echo "Executable not found DUMP_FASTSSIM=$DUMP_FASTSSIM"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ -z "$CORES" ]; then
  CORES=`grep -i processor /proc/cpuinfo | wc -l`
  #echo "CORES not set, using $CORES"
fi

find $@ -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB
