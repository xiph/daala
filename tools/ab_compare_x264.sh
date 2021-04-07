#!/bin/sh
set -e

while getopts 's:k:r:d:Y:y:n:' OPTIONS; do
  case $OPTIONS in
    s) SIZE="$OPTARG";;
    k) KEYINT="$OPTARG";;
    r) X264_ROOT="$OPTARG";;
    d) DAALA_ROOT="$OPTARG";;
    Y) Y4M2PNG="$OPTARG";;
    y) YUV2YUV4MPEG="$OPTARG";;
    n) FRAMES="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $X264_ROOT ] || [ ! -d $X264_ROOT ]; then
  echo "Please set X264_ROOT to the location of your x264 git clone"
  exit 1
fi

if [ -z "$X264" ]; then
  export X264=$X264_ROOT/x264
fi

if [ -z "$YUV2YUV4MPEG" ]; then
  YUV2YUV4MPEG=$DAALA_ROOT/tools/yuv2yuv4mpeg
fi

if [ -z "$Y4M2PNG" ]; then
  Y4M2PNG=$DAALA_ROOT/tools/y4m2png
fi

if [ ! -x "$X264" ]; then
  echo "x264 encoder not found at '$X264'."
  exit 1
fi

if [ ! -x "$YUV2YUV4MPEG" ]; then
  echo "YUV2YUV4MPEG not found at '$YUV2YUV4MPEG'."
  exit 1
fi

if [ ! -x "$Y4M2PNG" ]; then
  echo "Y4M2PNG not found at '$Y4M2PNG'."
  exit 1
fi

if [ -z "$FRAMES" ]; then
  FRAMES=1
fi

if [ -z "$KEYINT" ]; then
  KEYINT=256
fi

FILE=$1
echo $FILE
BASENAME=$(basename $FILE)
WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
QSTR="--preset placebo --crf=\$x"
# With x264, the lowest quantizer number yields the highest quality and vice
# versa. Here, MAX_QUALITY produces the best looking image, so it's the
# lowest number.
MIN_QUALITY=51
MAX_QUALITY=1
while (( $MIN_QUALITY - $MAX_QUALITY > 1 )); do
  QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
  X264_FILE=$BASENAME-$QUALITY.x264.tmp
  $X264 --dump-yuv $X264_FILE.yuv $(echo $QSTR | sed 's/\$x/'$QUALITY'/g') --keyint $KEYINT -o $X264_FILE $FILE 2> /dev/null > /dev/null
  X264_SIZE=$(wc -c $X264_FILE | awk '{ print $1 }')
  if (($X264_SIZE > $SIZE)); then
    MAX_QUALITY=$QUALITY
    MAX_QUALITY_SIZE=$X264_SIZE
  else
    MIN_QUALITY=$QUALITY
    MIN_QUALITY_SIZE=$X264_SIZE
  fi
done

if [ $MIN_QUALITY -eq 51 ]; then
  X264_FILE="$BASENAME-$MIN_QUALITY.x264.tmp"
  $X264 --dump-yuv "$X264_FILE.yuv" $(echo $QSTR | sed 's/\$x/'$MIN_QUALITY'/g') --keyint $KEYINT -o "$X264_FILE" $FILE 2> /dev/null > /dev/null
  MIN_QUALITY_SIZE=$(stat -c %s "$X264_FILE")
fi

if [ $MAX_QUALITY -eq 1 ]; then
  X264_FILE="$BASENAME-$MAX_QUALITY.x264.tmp"
  $X264 --dump-yuv "$X264_FILE.yuv" $(echo $QSTR | sed 's/\$x/'$MAX_QUALITY'/g') --keyint $KEYINT -o "$X264_FILE" $FILE 2> /dev/null > /dev/null
  MAX_QUALITY_SIZE=$(stat -c %s "$X264_FILE")
fi

if (( $MAX_QUALITY_SIZE - $SIZE < $SIZE - $MIN_QUALITY_SIZE )); then
  BEST_QUALITY=$MAX_QUALITY
else
  BEST_QUALITY=$MIN_QUALITY
fi

BEST_FILE=$BASENAME-$BEST_QUALITY.x264
mv $BEST_FILE.tmp $BEST_FILE
$YUV2YUV4MPEG $BEST_FILE.tmp -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2

if [ $FRAMES -eq 1 ]; then
  $Y4M2PNG -o $BEST_FILE.png $BEST_FILE.tmp.y4m
  rm $BEST_FILE.tmp.y4m
else
  mv $BEST_FILE.tmp.y4m $BEST_FILE.y4m
fi

rm $BASENAME-*.x264.tmp $BASENAME-*.x264.tmp.yuv
