#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh libaom *.y4m"
  exit 1
fi

FILE=$1
QP=$2

BASENAME=$(basename $FILE)-$CODEC-$QP
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')

case $CODEC in
libaom)
  QSTR="--ivf --frame-parallel=0 --auto-alt-ref=2 --passes=2 --threads=1 --kf-min-dist=1000 --kf-max-dist=1000 --lag-in-frames=25 --end-usage=q --cq-level=\$x"
  ;;
libaom-rt)
  QSTR="--ivf --frame-parallel=0 --passes=1 --threads=1 --kf-min-dist=1000 --kf-max-dist=1000 --lag-in-frames=0 --end-usage=q --cq-level=\$x"
  ;;
esac

ENCTIME=$BASENAME-enctime.out
TIMER='time -v --output='"$ENCTIME"
$TIMER $AOMENC $EXTRA_OPTS $(echo $QSTR | sed 's/\$x/'$QP'/g') -o $BASENAME.ivf $FILE 2> $BASENAME-enc.out
$AOMDEC -o $BASENAME.y4m $BASENAME.ivf
SIZE=$(wc -c $BASENAME.ivf | awk '{ print $1 }')
$DUMP_PSNR $FILE $BASENAME.y4m > $BASENAME-psnr.out 2> /dev/null
FRAMES=$(cat $BASENAME-psnr.out | grep ^0 | wc -l)
PIXELS=$(($WIDTH*$HEIGHT*$FRAMES))
TIME=$(cat $ENCTIME | grep User | cut -d\  -f4)
PSNR=$(cat $BASENAME-psnr.out | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
CIEDE=$($DUMP_CIEDE $FILE $BASENAME.y4m 2> /dev/null | grep Total | cut -d' ' -f2-)
rm $BASENAME.ivf $BASENAME.y4m $BASENAME-enc.out $BASENAME-psnr.out $ENCTIME
echo $QP $PIXELS $SIZE $TIME $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
#tail -1 $BASENAME.out
