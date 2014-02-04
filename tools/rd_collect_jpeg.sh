#!/bin/bash
set -e

if [ $# == 0 ]; then
  echo "usage: BUILD_ROOT=<build_dir> $0 *.y4m"
  exit 1
fi

if [ -z $BUILD_ROOT ]; then
  BUILD_ROOT=.
fi

if [ -z "$PLANE" ]; then
  PLANE=0
fi

if [ $PLANE != 0 ] && [ $PLANE != 1 ] && [ $PLANE != 2 ]; then
  echo "Invalid plane $PLANE. Must be 0, 1 or 2."
  exit 1
fi

if [ -z "$YUVJPEG" ]; then
  YUVJPEG=$BUILD_ROOT/tools/yuvjpeg
fi

if [ -z "$JPEGYUV" ]; then
  JPEGYUV=$BUILD_ROOT/tools/jpegyuv
fi

if [ -z "$YUV2YUV4MPEG" ]; then
  YUV2YUV4MPEG=$BUILD_ROOT/tools/yuv2yuv4mpeg
fi

if [ -z "$DUMP_PSNRHVS" ]; then
  DUMP_PSNR=$BUILD_ROOT/tools/dump_psnr
fi

if [ -z "$DUMP_PSNRHVS" ]; then
  DUMP_PSNRHVS=$BUILD_ROOT/tools/dump_psnrhvs
fi

if [ -z "$DUMP_SSIM" ]; then
  DUMP_SSIM=$BUILD_ROOT/tools/dump_ssim
fi

if [ -z "$DUMP_FASTSSIM" ]; then
  DUMP_FASTSSIM=$BUILD_ROOT/tools/dump_fastssim
fi

if [ ! -x "$YUVJPEG" ]; then
  echo "Executable not found YUVJPEG=$YUVJPEG"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$JPEGYUV" ]; then
  echo "Executable not found JPEGYUV=$JPEGYUV"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$YUV2YUV4MPEG" ]; then
  echo "Executable not found YUV2YUV4MPEG=$YUV2YUV4MPEG"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNR" ]; then
  echo "Executable not found DUMP_PSNR=$DUMP_PSNR"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNRHVS" ]; then
  echo "Executable not found DUMP_PSNRHVS=$DUMP_PSNRHVS"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_SSIM" ]; then
  echo "Executable not found DUMP_SSIM=$DUMP_SSIM"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_FASTSSIM" ]; then
  echo "Executable not found DUMP_FASTSSIM=$DUMP_FASTSSIM"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

RD_COLLECT_SUB=$(dirname "$0")/rd_collect_sub_jpeg.sh

if [ -z "$CORES" ]; then
  CORES=`grep -i processor /proc/cpuinfo | wc -l`
  #echo "CORES not set, using $CORES"
fi

find $@ -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB $PLANE $YUVJPEG $JPEGYUV $YUV2YUV4MPEG $DUMP_PSNR $DUMP_PSNRHVS $DUMP_SSIM $DUMP_FASTSSIM
