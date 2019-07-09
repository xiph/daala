#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh bpg *.y4m"
  exit 1
fi

FILE=$1
QP=$2

BASENAME=$(basename $FILE)-$CODEC-$QP
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME
$Y4M2PNG $FILE -o $BASENAME-in.png 2> /dev/null

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
PIXELS=$(($WIDTH*$HEIGHT))

$BPGENC -q $QP $BASENAME-in.png -o $BASENAME.bpg
$BPGDEC $BASENAME.bpg -o $BASENAME-out.png -b 16
$PNG2Y4M --no-dither $BASENAME-out.png -o $BASENAME.y4m 2> /dev/null
SIZE=$(stat -c %s $BASENAME.bpg)
PSNR=$($DUMP_PSNR $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
CIEDE=$($DUMP_CIEDE $FILE $BASENAME.y4m 2> /dev/null | grep Total | cut -d' ' -f2-)
rm $BASENAME.bpg $BASENAME-out.png $BASENAME.y4m
echo $QP $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
#tail -1 $BASENAME.out

rm $BASENAME-in.png
