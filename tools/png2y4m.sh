#!/bin/bash

# Converts PNG files to Y4M.
#
# This script is different from png2y4m.c in that it:
#   - Uses JPEG color levels instead of video levels
#   - Uses the ITU-RBT.601 matrix instead of ITU-RBT.709 matrix
#   - Does not dither the input image

if [ $# == 0 ]; then
  echo "usage: BUILD_ROOT=<build_dir> $0 *.png"
  exit 1
fi

if [ -z $BUILD_ROOT ]; then
  BUILD_ROOT=.
fi

if [ -z "$YUV2YUV4MPEG" ]; then
  YUV2YUV4MPEG=$BUILD_ROOT/tools/yuv2yuv4mpeg
fi

if [ ! -f "$YUV2YUV4MPEG" ]; then
  echo "File not found YUV2YUV4MPEG=$YUV2YUV4MPEG"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

FILES=$(find $@ -type f -name "*.png")

for f in $FILES; do
  BASENAME=${f%.*}
  convert png:$BASENAME.png -sampling-factor 4:2:0 -depth 8 $BASENAME.yuv
  SIZE=$(identify $BASENAME.png | sed -re 's/^.*PNG\ +([0-9]+x[0-9]+).*$/\1/')
  WIDTH=$(echo $SIZE | cut -dx -f1)
  HEIGHT=$(echo $SIZE | cut -dx -f2)
  $YUV2YUV4MPEG $BASENAME -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  rm $BASENAME.yuv
done
