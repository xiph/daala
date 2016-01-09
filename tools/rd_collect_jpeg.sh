#!/bin/sh
set -e

if [ -z $RD_COLLECT_SUB ]; then
  echo "Please use: $(dirname $0)/rd_collect.sh <libjpeg|mozjpeg> *.y4m"
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
  $YUVJPEG $x "$WIDTH"x$HEIGHT $BASENAME-in.yuv $BASENAME.jpeg
  $JPEGYUV $BASENAME.jpeg $BASENAME.yuv
  $YUV2YUV4MPEG $BASENAME -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  SIZE=$(wc -c $BASENAME.jpeg | awk '{ print $1 }')
  PSNR=$($DUMP_PSNR $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  PSNRHVS=$($DUMP_PSNRHVS $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  SSIM=$($DUMP_SSIM $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  FASTSSIM=$($DUMP_FASTSSIM -c $FILE $BASENAME.y4m 2> /dev/null | grep Total | tr -s ' ' | cut -d\  -f $((4+$PLANE*2)))
  rm $BASENAME.jpeg $BASENAME.yuv $BASENAME.y4m
  echo -$x $PIXELS $SIZE $PSNR $PSNRHVS $SSIM $FASTSSIM >> $BASENAME.out
  #tail -1 $BASENAME.out
done

rm $BASENAME-in.yuv
