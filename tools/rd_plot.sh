#!/bin/bash

# Use this to average data from multiple runs
#awk '{size[FNR]+=$2;bytes[FNR]+=$3;psnr[FNR]+=$2*$4;psnrhvs[FNR]+=$2*$5;ssim[FNR]+=$2*$6;fastssim[FNR]+=$2*$7;}END{for(i=1;i<=FNR;i++)print i+1,size[i],bytes[i],psnr[i]/size[i],psnrhvs[i]/size[i],ssim[i]/size[i],fastssim[i]/size[i];}' *.out > total.out

if [ -z "$IMAGE" ]; then
  IMAGE=psnrhvs.png
fi

if [ $# == 0 ]; then
  echo "no parameters"
  exit 1
fi

CMDS="$CMDS set term png;"
CMDS="$CMDS set output \"$IMAGE\";"
CMDS="$CMDS set xlabel 'Bits/Pixel';"
CMDS="$CMDS set ylabel 'dB';"
CMDS="$CMDS set key bot right;"

CMDS="$CMDS plot"
for FILE in "$@"; do
  BASENAME=$(basename $FILE)
  CMDS="$CMDS $PREFIX '$FILE' using (\$3*8/\$2):4 with lines title '${BASENAME%.*} (PSNR)'"
  PREFIX=","
  CMDS="$CMDS $PREFIX '$FILE' using (\$3*8/\$2):5 with lines title '${BASENAME%.*} (PSNR-HSV)'"
  PREFIX=","
  #CMDS="$CMDS $PREFIX '$FILE' using (\$3*8/\$2):6 with lines title '${BASENAME%.*} (SSIM)'"
  #PREFIX=","
  #CMDS="$CMDS $PREFIX '$FILE' using (\$3*8/\$2):7 with lines title '${BASENAME%.*} (FAST SSIM)'"
  #PREFIX=","
done
echo $CMDS
CMDS="$CMDS;"

gnuplot -e "$CMDS" 2> /dev/null
