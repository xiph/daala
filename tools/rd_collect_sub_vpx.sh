#!/bin/bash
set -e

PLANE=$1
CODEC=$2
VPXENC=$3
VPXDEC=$4
DUMP_PSNR=$5
DUMP_PSNRHVS=$6
DUMP_SSIM=$7
DUMP_FASTSSIM=$8
FILE=$9

BASENAME=$(basename $FILE)
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME
WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

if [ $CODEC == vp8 ]; then
  RANGE=$(seq 50 100 5000)
  QUALITY=--target-bitrate
else
  RANGE=$(seq 3 63)
  QUALITY="--end-usage=q --cq-level"
fi

for x in $RANGE; do
  $VPXENC --codec=$CODEC --good --cpu-used=0 $QUALITY=$x -o $BASENAME.vpx $FILE 2> $BASENAME-$x-enc.out
  $VPXDEC --codec=$CODEC -o $BASENAME.y4m $BASENAME.vpx
  PIXELS=$(($WIDTH*$HEIGHT))
  SIZE=$(wc -c $BASENAME.vpx | awk '{ print $1 }')
  PSNR=$($DUMP_PSNR $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  rm $BASENAME.vpx $BASENAME.y4m $BASENAME-$x-enc.out
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done
