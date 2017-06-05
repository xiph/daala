#!/bin/bash

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

awk '{size[$1]+=$2;bytes[$1]+=$3;psnr[$1]+=$2*$4;psnrhvs[$1]+=$2*$5;ssim[$1]+=$2*$6;fastssim[$1]+=$2*$7;ciede[$1]+=$2*$8;}END{for(i in size)print i,size[i],bytes[i],psnr[i]/size[i],psnrhvs[i]/size[i],ssim[i]/size[i],fastssim[i]/size[i],ciede[i]/size[i];}' $@ | sort -n > $TOTAL
