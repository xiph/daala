#!/bin/bash
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh vtm *.y4m"
  exit 1
fi

FILE=$1
QP=$2

BASENAME=$(basename $FILE)-$CODEC-$QP
rm $BASENAME.out 2> /dev/null || true
echo $BASENAME

WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
FRAMES=$(grep FRAME $FILE | wc -l)
FPS=$(ffmpeg -i $FILE 2>&1 | sed -n "s/.*, \(.*\) fp.*/\1/p")

QSTR="-wdt $WIDTH -hgt $HEIGHT -f $FRAMES -fr $FPS -c $VTM_ROOT/cfg/encoder_randomaccess_vtm.cfg -q $QP --OutputBitDepth=8"

ENCTIME=$BASENAME-enctime.out
TIMER='time -v --output='"$ENCTIME"
$Y4M2YUV -o $BASENAME.yuv $FILE
$TIMER $VTM $EXTRA_OPTS $QSTR -o $BASENAME-out.yuv -b $BASENAME.vvc -i $BASENAME.yuv > $BASENAME-enc.out
$YUV2YUV4MPEG $BASENAME-out -w$WIDTH -h$HEIGHT
SIZE=$(wc -c $BASENAME.vvc | awk '{ print $1 }')
$DUMP_PSNR $FILE $BASENAME-out.y4m > $BASENAME-psnr.out 2> /dev/null
FRAMES=$(cat $BASENAME-psnr.out | grep ^0 | wc -l)
PIXELS=$(($WIDTH*$HEIGHT*$FRAMES))
TIME=$(cat $ENCTIME | grep User | cut -d\  -f4)
PSNR=$(cat $BASENAME-psnr.out | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME-out.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
SSIM=$($DUMP_SSIM $FILE $BASENAME-out.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME-out.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
CIEDE=$($DUMP_CIEDE $FILE $BASENAME-out.yuv 2> /dev/null | grep Total | cut -d' ' -f2-)
rm $BASENAME-out.yuv $BASENAME.vvc $BASENAME.yuv $BASENAME-enc.out $BASENAME-enctime.out $BASENAME-psnr.out $BASENAME-out.y4m
echo $QP $PIXELS $SIZE $TIME $PSNR $PSNRHVS $SSIM $FASTSSIM $CIEDE >> $BASENAME.out
#tail -1 $BASENAME.out
