#!/bin/bash
set -e

PLANE=$1
YUVJPEG=$2
JPEGYUV=$3
YUV2YUV4MPEG=$4
DUMP_PSNR=$5
DUMP_PSNRHVS=$6
DUMP_SSIM=$7
DUMP_FASTSSIM=$8
FILE=$9

BASENAME=$(basename $FILE)
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME
tail -n+3 $FILE > $BASENAME-in.yuv
WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

for x in {0..100}; do
  $YUVJPEG $x "$WIDTH"x$HEIGHT $BASENAME-in.yuv $BASENAME.jpeg
  $JPEGYUV $BASENAME.jpeg $BASENAME.yuv
  $YUV2YUV4MPEG $BASENAME -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  PIXELS=$(($WIDTH*$HEIGHT))
  SIZE=$(wc -c $BASENAME.jpeg | cut -d\  -f 1)
  PSNR=$($DUMP_PSNR $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  rm $BASENAME.jpeg $BASENAME.yuv $BASENAME.y4m
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done

rm $BASENAME-in.yuv
