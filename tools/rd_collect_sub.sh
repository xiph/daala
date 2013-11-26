#!/bin/bash

PLANE=$1
ENCODER_EXAMPLE=$2
DUMP_VIDEO=$3
DUMP_PSNRHVS=$4
DUMP_SSIM=$5
DUMP_FASTSSIM=$6
TMP_DIR=$7
FILE=$8

BASENAME=$(basename $FILE)
rm $BASENAME.out 2> /dev/null
echo $BASENAME

for x in {2..40}; do
  OD_LOG_MODULES='encoder:10' $ENCODER_EXAMPLE -v $x $FILE -o $TMP_DIR/$BASENAME-$x.ogv 2> $TMP_DIR/$BASENAME-$x-enc.out
  PIXELS=$(cat $TMP_DIR/$BASENAME-$x-enc.out | grep "Plane $PLANE" | sed -re 's/^.*Pixels:\ +([0-9]+).*$/\1/')
  PSNR=$(cat $TMP_DIR/$BASENAME-$x-enc.out | grep "Plane $PLANE" | sed -re 's/^.*PSNR:\ +([0-9\.]+).*$/\1/')
  SIZE=$(cat $TMP_DIR/$BASENAME-$x-enc.out | grep 'Output' | sed -re 's/^.*Bytes:\ +([0-9]+).*$/\1/')
  $DUMP_VIDEO $TMP_DIR/$BASENAME-$x.ogv -o $TMP_DIR/$BASENAME-$x.y4m 2> /dev/null
  PSNRHVS=$($DUMP_PSNRHVS $FILE $TMP_DIR/$BASENAME-$x.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  if [ -f "$DUMP_SSIM" ]; then
    SSIM=$($DUMP_SSIM $FILE $TMP_DIR/$BASENAME-$x.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  fi
  if [ -f "$DUMP_FASTSSIM" ]; then
    FASTSSIM=$($DUMP_FASTSSIM -c $FILE $TMP_DIR/$BASENAME-$x.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  fi
  rm $TMP_DIR/$BASENAME-$x.ogv $TMP_DIR/$BASENAME-$x.y4m $TMP_DIR/$BASENAME-$x-enc.out
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done
