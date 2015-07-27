#!/bin/bash
set -e

while getopts 's:c:k:r:E:D:d:Y:y:n:' OPTIONS; do
  case $OPTIONS in
    s) SIZE="$OPTARG";;
    c) CODEC="$OPTARG";;
    k) KEYINT="$OPTARG";;
    r) LIBVPX_ROOT="$OPTARG";;
    E) VPXENC="$OPTARG";;
    D) VPXDEC="$OPTARG";;
    d) DAALA_ROOT="$OPTARG";;
    Y) Y4M2PNG="$OPTARG";;
    y) YUV2YUV4MPEG="$OPTARG";;
    n) FRAMES="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [ -z $CODEC ]; then
  CODEC=vp8
fi

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $LIBVPX_ROOT ] || [ ! -d $LIBVPX_ROOT ]; then
  echo "Please set LIBVPX_ROOT to the location of your libvpx git clone"
  exit 1
fi

if [ -z "$Y4M2PNG" ]; then
  Y4M2PNG=$DAALA_ROOT/tools/y4m2png
fi

if [ -z "$VPXENC" ]; then
  export VPXENC=$LIBVPX_ROOT/vpxenc
fi

if [ -z "$VPXDEC" ]; then
  export VPXDEC=$LIBVPX_ROOT/vpxdec
fi

if [ ! -x "$Y4M2PNG" ]; then
  echo "Y4M2PNG not found at '$Y4M2PNG'."
  exit 1
fi

if [ ! -x "$VPXENC" ]; then
  echo "vpx encoder not found at '$VPXENC'."
  exit 1
fi

if [ ! -x "$VPXDEC" ]; then
  echo "vpx decoder not found at '$VPXDEC'."
  exit 1
fi

if [ -z "$FRAMES" ]; then
  FRAMES=1
fi

if [ -z "$KEYINT" ]; then
  KEYINT=256;
fi

FILE=$1
echo $FILE
BASENAME=$(basename $FILE)
# With libvpx, the lowest quantizer number yields the highest quality and
# vice versa. Here, MAX_QUALITY produces the best looking image, so it's the
# lowest number.
MAX_QUALITY=3
MIN_QUALITY=63
while (( $MIN_QUALITY - $MAX_QUALITY > 1 )); do
  QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
  VPX_FILE=$BASENAME-$QUALITY.$CODEC.tmp
  $VPXENC --codec=$CODEC --good --cpu-used=0 -y --min-q=$QUALITY --max-q=$QUALITY --kf-max-dist=$KEYINT -o $VPX_FILE $FILE 2> /dev/null
  VPX_SIZE=$(wc -c $VPX_FILE | awk '{ print $1 }')
  if (($VPX_SIZE > $SIZE)); then
    MAX_QUALITY=$QUALITY
    MAX_QUALITY_SIZE=$VPX_SIZE
  else
    MIN_QUALITY=$QUALITY
    MIN_QUALITY_SIZE=$VPX_SIZE
  fi
done

if [ $MIN_QUALITY -eq 63 ]; then
  $VPXENC --codec=$CODEC --good --cpu-used=0 -y --kf-max-dist=$KEYINT --min-q=$MIN_QUALITY --max-q=$MIN_QUALITY -o $BASENAME-$MIN_QUALITY.$CODEC.tmp $FILE 2> /dev/null
  MIN_QUALITY_SIZE=$(stat -c %s $BASENAME-$MIN_QUALITY.$CODEC.tmp)
fi

if [ $MAX_QUALITY -eq 3 ]; then
  $VPXENC --codec=$CODEC --good --cpu-used=0 -y --kf-max-dist=$KEYINT --min-q=$MAX_QUALITY --max-q=$MAX_QUALITY -o $BASENAME-$MAX_QUALITY.$CODEC.tmp $FILE 2> /dev/null
  MAX_QUALITY_SIZE=$(stat -c %s $BASENAME-$MAX_QUALITY.$CODEC.tmp)
fi

if (( $MAX_QUALITY_SIZE - $SIZE < $SIZE - $MIN_QUALITY_SIZE )); then
  VPX_SIZE=$MAX_QUALITY_SIZE
  BEST_QUALITY=$MAX_QUALITY
else
  BEST_QUALITY=$MIN_QUALITY
  VPX_SIZE=$MIN_QUALITY_SIZE
fi

BEST_FILE=$BASENAME-$BEST_QUALITY.$CODEC
$VPXDEC --codec=$CODEC -o $BEST_FILE.y4m $BEST_FILE.tmp
mv $BEST_FILE.tmp $BEST_FILE

if [ $FRAMES -eq 1 ]; then
  $Y4M2PNG -o $BEST_FILE.png $BEST_FILE.y4m
  rm $BEST_FILE.y4m
fi

rm $BASENAME-*.$CODEC.tmp
