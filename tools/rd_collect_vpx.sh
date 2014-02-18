#!/bin/bash
set -e

if [ $# -lt 2 ]; then
  echo "usage: BUILD_ROOT=<build_dir> $0 <vp8|vp9> *.y4m"
  exit 1
fi

CODEC=$1
shift

if [ -z $BUILD_ROOT ]; then
  BUILD_ROOT=.
fi

if [ -z $LIBVPX_ROOT ]; then
  LIBVPX_ROOT=$BUILD_ROOT/../libvpx
fi

if [ -z "$PLANE" ]; then
  PLANE=0
fi

if [ $PLANE != 0 ] && [ $PLANE != 1 ] && [ $PLANE != 2 ]; then
  echo "Invalid plane $PLANE. Must be 0, 1 or 2."
  exit 1
fi

if [ -z "$VPXENC" ]; then
  VPXENC=$LIBVPX_ROOT/vpxenc
fi

if [ -z "$VPXDEC" ]; then
  VPXDEC=$LIBVPX_ROOT/vpxdec
fi

if [ -z "$DUMP_PSNR" ]; then
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

if [ ! -x "$VPXENC" ]; then
  echo "Executable not found VPXENC=$VPXENC"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -x "$VPXDEC" ]; then
  echo "Executable not found VPXDEC=$VPXDEC"
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

RD_COLLECT_SUB=$(dirname "$0")/rd_collect_sub_vpx.sh

if [ -z "$CORES" ]; then
  CORES=`grep -i processor /proc/cpuinfo | wc -l`
  #echo "CORES not set, using $CORES"
fi

find $@ -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB $PLANE $CODEC $VPXENC $VPXDEC $DUMP_PSNR $DUMP_PSNRHVS $DUMP_SSIM $DUMP_FASTSSIM
