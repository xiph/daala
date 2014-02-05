#!/bin/bash
set -e

PLANE=$1
ENCODER_EXAMPLE=$2
DUMP_PSNRHVS=$3
DUMP_SSIM=$4
DUMP_FASTSSIM=$5
FILE=$6

BASENAME=$(basename $FILE)
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

for x in {2..40}; do
  OD_LOG_MODULES='encoder:10' OD_DUMP_IMAGES_SUFFIX=$BASENAME $ENCODER_EXAMPLE -v $x $FILE -o /dev/null 2> $BASENAME-$x-enc.out
  PIXELS=$(cat $BASENAME-$x-enc.out | grep "Plane $PLANE" | sed -re 's/^.*Pixels:\ +([0-9]+).*$/\1/')
  PSNR=$(cat $BASENAME-$x-enc.out | grep "Plane $PLANE" | sed -re 's/^.*PSNR:\ +([0-9\.]+).*$/\1/')
  SIZE=$(cat $BASENAME-$x-enc.out | grep 'Output' | sed -re 's/^.*Bytes:\ +([0-9]+).*$/\1/')
  PSNRHVS=$($DUMP_PSNRHVS $FILE 00000000$BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  if [ -f "$DUMP_SSIM" ]; then
    SSIM=$($DUMP_SSIM $FILE 00000000$BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  fi
  if [ -f "$DUMP_FASTSSIM" ]; then
    FASTSSIM=$($DUMP_FASTSSIM -c $FILE 00000000$BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  fi
  rm 00000000$BASENAME.y4m $BASENAME-$x-enc.out
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done
