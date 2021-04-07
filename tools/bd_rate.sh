#!/bin/sh
set -e

if [ $# == 0 ]; then
  echo "usage: BUILD_ROOT=<build_dir> $0 <RD-1.out> <RD-2.out>"
  exit 1
fi

if [ -z "$BUILD_ROOT" ]; then
  BUILD_ROOT=.
fi

if [ -z "$BJONTEGAARD" ]; then
  BJONTEGAARD=$BUILD_ROOT/tools/bjontegaard
fi

if [ ! -f "$BJONTEGAARD" ]; then
  echo "File not found BJONTEGAARD=$BJONTEGAARD"
  echo "Do you have the right BUILD_ROOT=$BUILD_ROOT"
  exit 1
fi

IFS=' ' read -r VALUES <<< "$(head -n 1 $1)"
VALUES=($VALUES)
NUM_VALUES_1=${#VALUES[@]}
IFS=' ' read -r VALUES <<< "$(head -n 1 $2)"
VALUES=($VALUES)
NUM_VALUES_2=${#VALUES[@]}

CIEDE_ROW=9

if [ $NUM_VALUES_1 != $NUM_VALUES_2 ]; then
  echo "Number of metrics in out files don't match"
  exit 1
fi

N1=$(cat $1 | wc -l)
N2=$(cat $2 | wc -l)

AREA1=$(cut -d\  -f 2 $1 | xargs | sed 's/ /,/g')
AREA2=$(cut -d\  -f 2 $2 | xargs | sed 's/ /,/g')
SIZE1=$(cut -d\  -f 3 $1 | xargs | sed 's/ /,/g')
SIZE2=$(cut -d\  -f 3 $2 | xargs | sed 's/ /,/g')
PSNR1=$(cut -d\  -f 5 $1 | xargs | sed 's/ /,/g')
PSNR2=$(cut -d\  -f 5 $2 | xargs | sed 's/ /,/g')
PSNRHVS1=$(cut -d\  -f 6 $1 | xargs | sed 's/ /,/g')
PSNRHVS2=$(cut -d\  -f 6 $2 | xargs | sed 's/ /,/g')
SSIM1=$(cut -d\  -f 7 $1 | xargs | sed 's/ /,/g')
SSIM2=$(cut -d\  -f 7 $2 | xargs | sed 's/ /,/g')
FASTSSIM1=$(cut -d\  -f 8 $1 | xargs | sed 's/ /,/g')
FASTSSIM2=$(cut -d\  -f 8 $2 | xargs | sed 's/ /,/g')

PSNR_RATE=$($BJONTEGAARD 0 $N1 $AREA1 $SIZE1 $PSNR1 $N2 $AREA2 $SIZE2 $PSNR2)
PSNR_DSNR=$($BJONTEGAARD 1 $N1 $AREA1 $SIZE1 $PSNR1 $N2 $AREA2 $SIZE2 $PSNR2)
PSNRHVS_RATE=$($BJONTEGAARD 0 $N1 $AREA1 $SIZE1 $PSNRHVS1 $N2 $AREA2 $SIZE2 $PSNRHVS2)
PSNRHVS_DSNR=$($BJONTEGAARD 1 $N1 $AREA1 $SIZE1 $PSNRHVS1 $N2 $AREA2 $SIZE2 $PSNRHVS2)
SSIM_RATE=$($BJONTEGAARD 0 $N1 $AREA1 $SIZE1 $SSIM1 $N2 $AREA2 $SIZE2 $SSIM2)
SSIM_DSNR=$($BJONTEGAARD 1 $N1 $AREA1 $SIZE1 $SSIM1 $N2 $AREA2 $SIZE2 $SSIM2)
FASTSSIM_RATE=$($BJONTEGAARD 0 $N1 $AREA1 $SIZE1 $FASTSSIM1 $N2 $AREA2 $SIZE2 $FASTSSIM2)
FASTSSIM_DSNR=$($BJONTEGAARD 1 $N1 $AREA1 $SIZE1 $FASTSSIM1 $N2 $AREA2 $SIZE2 $FASTSSIM2)

echo "           RATE (%) DSNR (dB)"
echo "    PSNR" $(echo $PSNR_RATE     | cut -d\  -f 3) $(echo $PSNR_DSNR     | cut -d\  -f 3)
echo " PSNRHVS" $(echo $PSNRHVS_RATE  | cut -d\  -f 3) $(echo $PSNRHVS_DSNR  | cut -d\  -f 3)
echo "    SSIM" $(echo $SSIM_RATE     | cut -d\  -f 3) $(echo $SSIM_DSNR     | cut -d\  -f 3)
echo "FASTSSIM" $(echo $FASTSSIM_RATE | cut -d\  -f 3) $(echo $FASTSSIM_DSNR | cut -d\  -f 3)

if [ $NUM_VALUES_1 -ge $CIEDE_ROW ]; then
  CIEDE1=$(cut -d\  -f 9 $1 | xargs | sed 's/ /,/g')
  CIEDE2=$(cut -d\  -f 9 $2 | xargs | sed 's/ /,/g')
  CIEDE_RATE=$($BJONTEGAARD 0 $N1 $AREA1 $SIZE1 $CIEDE1 $N2 $AREA2 $SIZE2 $CIEDE2)
  CIEDE_DSNR=$($BJONTEGAARD 1 $N1 $AREA1 $SIZE1 $CIEDE1 $N2 $AREA2 $SIZE2 $CIEDE2)
  echo "   CIEDE" $(echo $CIEDE_RATE    | cut -d\  -f 3) $(echo $CIEDE_DSNR    | cut -d\  -f 3)
fi
