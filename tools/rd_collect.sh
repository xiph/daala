#!/bin/bash
set -e

if [ $# == 0 ]; then
  echo "usage: BUILD_ROOT=<build_dir> $0 *.y4m"
  exit 1
fi

if [ -z $BUILD_ROOT ]; then
  BUILD_ROOT=.
fi

if ! grep -Fxq "#define OD_LOGGING_ENABLED 1" "$BUILD_ROOT/config.h"; then
  echo "Logging not enabled, re-run configure with --enable-logging"
  exit 1
fi

if ! grep -Fxq "#define OD_DUMP_IMAGES 1" "$BUILD_ROOT/config.h"; then
  echo "Image dumping not enabled, re-run configure without --disable-dump-images"
  exit 1
fi

if [ -z "$PLANE" ]; then
  PLANE=0
fi

if [ $PLANE != 0 ] && [ $PLANE != 1 ] && [ $PLANE != 2 ]; then
  echo "Invalid plane $PLANE. Must be 0, 1 or 2."
  exit 1
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$BUILD_ROOT/examples/encoder_example
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

if [ ! -x "$ENCODER_EXAMPLE" ]; then
  echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
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

RD_COLLECT_SUB=$(dirname "$0")/rd_collect_sub.sh

if [ -z "$CORES" ]; then
  CORES=`grep -i processor /proc/cpuinfo | wc -l`
  #echo "CORES not set, using $CORES"
fi

find $@ -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB $PLANE $ENCODER_EXAMPLE $DUMP_PSNRHVS $DUMP_SSIM $DUMP_FASTSSIM
