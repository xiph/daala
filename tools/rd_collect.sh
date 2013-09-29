#!/bin/bash

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

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$BUILD_ROOT/examples/encoder_example
fi

if [ -z "$DUMP_VIDEO" ]; then
  DUMP_VIDEO=$BUILD_ROOT/examples/dump_video
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

if [ ! -f "$ENCODER_EXAMPLE" ]; then
  echo "File not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -f "$DUMP_VIDEO" ]; then
  echo "File not found DUMP_VIDEO=$DUMP_VIDEO"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -f "$DUMP_PSNRHVS" ]; then
  echo "File not found DUMP_PSNRHVS=$DUMP_PSNRHVS"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -f "$DUMP_SSIM" ]; then
  echo "File not found DUMP_SSIM=$DUMP_SSIM"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ ! -f "$DUMP_FASTSSIM" ]; then
  echo "File not found DUMP_FASTSSIM=$DUMP_FASTSSIM"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

if [ -z "$TMP_DIR" ]; then
  TMP_DIR=.images.tmp-$BASHPID
fi

RD_COLLECT_SUB=$(dirname "$0")/rd_collect_sub.sh

if [ -z "$CORES" ]; then
  CORES=`grep -i processor /proc/cpuinfo | wc -l`
  #echo "CORES not set, using $CORES"
fi

mkdir -p $TMP_DIR

find $@ -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB $ENCODER_EXAMPLE $DUMP_VIDEO $DUMP_PSNRHVS $DUMP_SSIM $DUMP_FASTSSIM $TMP_DIR

rm -r $TMP_DIR 2> /dev/null
