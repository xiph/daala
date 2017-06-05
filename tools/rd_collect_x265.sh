#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh x265 *.y4m"
  exit 1
fi

FILE=$1

BASENAME=$(basename $FILE)-$CODEC
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

RANGE=$(seq 1 51)
QSTR="--preset slow --frame-threads 1 --min-keyint 256 --keyint 256 --no-scenecut --crf=\$x"

for x in $RANGE; do
  $X265 -r $BASENAME.y4m $(echo $QSTR | sed 's/\$x/'$x'/g') -o $BASENAME.x265 $FILE 2> $BASENAME-$x-enc.out > /dev/null
  SIZE=$(wc -c $BASENAME.x265 | awk '{ print $1 }')
  $DUMP_PSNR $FILE $BASENAME.y4m > $BASENAME-$x-psnr.out 2> /dev/null
  FRAMES=$(cat $BASENAME-$x-psnr.out | grep ^0 | wc -l)
  PIXELS=$(($WIDTH*$HEIGHT*$FRAMES))
  PSNR=$(cat $BASENAME-$x-psnr.out | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  CIEDE=$($DUMP_CIEDE $FILE $BASENAME.y4m 2> /dev/null | grep Total | cut -d' ' -f2-)
  rm $BASENAME.x265 $BASENAME.y4m $BASENAME-$x-enc.out $BASENAME-$x-psnr.out
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
  #tail -1 $BASENAME.out
done
