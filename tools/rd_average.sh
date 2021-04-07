#!/bin/sh

if [ $# == 0 ]; then
  echo "usage: OUTPUT=<label> $0 *.out"
  exit 1
fi

TOTAL=total.out

if [ -n "$OUTPUT" ]; then
  TOTAL="$OUTPUT.out"
fi

if [ -e "$TOTAL" ]; then
  echo "ERROR: $TOTAL already exists and will be included in average, please remove it first"
  exit 1
fi

IFS=' ' read -r VALUES <<< "$(head -n 1 $1)"
VALUES=($VALUES)
NUM_VALUES=${#VALUES[@]}
CIEDE_ROW=9

AWK_SUM='size[$1]+=$2;bytes[$1]+=$3;time[$1]+=$4;psnr[$1]+=$2*$5;psnrhvs[$1]+=$2*$6;ssim[$1]+=$2*$7;fastssim[$1]+=$2*$8;'
AWK_DIV='psnr[i]/size[i],psnrhvs[i]/size[i],ssim[i]/size[i],fastssim[i]/size[i]'

if [ $NUM_VALUES -ge $CIEDE_ROW ]; then
  AWK_SUM+='ciede[$1]+=$2*$9;'
  AWK_DIV+=',ciede[i]/size[i]'
fi

AWK_CMD="{$AWK_SUM}END{for(i in size)print i,size[i],bytes[i],time[i],$AWK_DIV;}"

awk "$AWK_CMD" $@ | sort -n > $TOTAL
