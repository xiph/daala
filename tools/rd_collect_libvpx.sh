#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh <vp8|vp9> *.y4m"
  exit 1
fi

FILE=$1

BASENAME=$(basename $FILE)-$CODEC
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

RANGE=$(seq 1 63)
case $CODEC in
vp8)
  QSTR="--end-usage=cq --cq-level=\$x"
  ;;
vp9)
  QSTR="--end-usage=q --cq-level=\$x"
  ;;
esac

for x in $RANGE; do
  $VPXENC --codec=$CODEC --good --cpu-used=0 --ivf $(echo $QSTR | sed 's/\$x/'$x'/g') -o $BASENAME.ivf $FILE 2> $BASENAME-$x-enc.out
  $VPXDEC --codec=$CODEC -o $BASENAME.y4m $BASENAME.ivf
  SIZE=$(wc -c $BASENAME.ivf | awk '{ print $1 }')
  $DUMP_PSNR $FILE $BASENAME.y4m > $BASENAME-$x-psnr.out 2> /dev/null
  FRAMES=$(cat $BASENAME-$x-psnr.out | grep ^0 | wc -l)
  PIXELS=$(($WIDTH*$HEIGHT*$FRAMES))
  PSNR=$(cat $BASENAME-$x-psnr.out | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  rm $BASENAME.ivf $BASENAME.y4m $BASENAME-$x-enc.out $BASENAME-$x-psnr.out
  echo $x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done
