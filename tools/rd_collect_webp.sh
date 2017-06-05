#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh webp *.y4m"
  exit 1
fi

FILE=$1

BASENAME=$(basename $FILE)
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME
tail -n+3 $FILE > $BASENAME-in.yuv

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
PIXELS=$(($WIDTH*$HEIGHT))

for x in $(seq 100 -1 0); do
  $CWEBP -q $x -s $WIDTH $HEIGHT $BASENAME-in.yuv -o $BASENAME.webp 2> /dev/null
  $DWEBP $BASENAME.webp -yuv -o $BASENAME.yuv 2> /dev/null
  $YUV2YUV4MPEG $BASENAME -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  SIZE=$(stat -c %s $BASENAME.webp)
  PSNR=$($DUMP_PSNR $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  CIEDE=$($DUMP_CIEDE $FILE $BASENAME.y4m 2> /dev/null | grep Total | cut -d' ' -f2-)
  rm $BASENAME.webp $BASENAME.yuv $BASENAME.y4m
  echo -$x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
  #tail -1 $BASENAME.out
done

rm $BASENAME-in.yuv
