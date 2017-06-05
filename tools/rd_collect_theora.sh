#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh <theora> *.y4m"
  exit 1
fi

FILE=$1

BASENAME=$(basename $FILE)-$CODEC
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

for x in $(seq 63 -1 0); do
    V=$(echo $x | awk '{print $1*0.15873015873015}');
    $ENCODER_EXAMPLE -o $BASENAME.ogv -v $V -z 0 $FILE 2> $BASENAME-$x-enc.err >/dev/null
    $DUMP_VIDEO $BASENAME.ogv -c -o $BASENAME.y4m 2>/dev/null
    SIZE=$(wc -c $BASENAME.ogv | awk '{ print $1 }')
    $DUMP_PSNR $FILE $BASENAME.y4m > $BASENAME-$x-psnr.out 2> /dev/null
    FRAMES=$(cat $BASENAME-$x-psnr.out | grep ^0 | wc -l)
    PIXELS=$(($WIDTH*$HEIGHT*$FRAMES))
    PSNR=$(cat $BASENAME-$x-psnr.out | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
    PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
    SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
    FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
    CIEDE=$($DUMP_CIEDE $FILE $BASENAME.y4m 2> /dev/null | grep Total | cut -d' ' -f2-)
    rm $BASENAME.ogv $BASENAME.y4m  $BASENAME-$x-enc.err $BASENAME-$x-psnr.out
    echo -$x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
    #tail -1 $BASENAME.out
done
