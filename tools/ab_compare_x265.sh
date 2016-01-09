#!/bin/sh
set -e

while getopts 's:k:r:d:Y:n:' OPTIONS; do
  case $OPTIONS in
    s) SIZE="$OPTARG";;
    k) KEYINT="$OPTARG";;
    r) X265_ROOT="$OPTARG";;
    d) DAALA_ROOT="$OPTARG";;
    Y) Y4M2PNG="$OPTARG";;
    n) FRAMES="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $X265_ROOT ] || [ ! -d $X265_ROOT ]; then
  echo "Please set X265_ROOT to the location of your x265 hg checkout"
  exit 1
fi

if [ -z "$X265" ]; then
  export X265=$X265_ROOT/build/linux/x265
fi

if [ -z "$Y4M2PNG" ]; then
  Y4M2PNG=$DAALA_ROOT/tools/y4m2png
fi

if [ ! -x "$X265" ]; then
  echo "x265 encoder not found at '$X265'."
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
QSTR="--crf=\$x"
# With x265, the lowest quantizer number yields the highest quality and vice
# versa. Here, MAX_QUALITY produces the best looking image, so it's the
# lowest number.
MAX_QUALITY=1
MIN_QUALITY=51
while (( $MIN_QUALITY - $MAX_QUALITY > 1 )); do
  QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
  X265_FILE=$BASENAME-$QUALITY.x265.tmp
  $X265 -r $X265_FILE.y4m $(echo $QSTR | sed 's/\$x/'$QUALITY'/g') --keyint $KEYINT -o $X265_FILE $FILE 2> /dev/null > /dev/null
  X265_SIZE=$(wc -c $X265_FILE | awk '{ print $1 }')
  if (($X265_SIZE > $SIZE)); then
    MAX_QUALITY=$QUALITY
    MAX_QUALITY_SIZE=$X265_SIZE
  else
    MIN_QUALITY=$QUALITY
    MIN_QUALITY_SIZE=$X265_SIZE
  fi
done

if [ $MIN_QUALITY -eq 51 ]; then
  X265_FILE="$BASENAME-$MIN_QUALITY.x265.tmp"
  $X265 -r "$X265_FILE.y4m" $(echo $QSTR | sed 's/\$x/'$MIN_QUALITY'/g') --keyint $KEYINT -o "$X265_FILE" $FILE 2> /dev/null > /dev/null
  MIN_QUALITY_SIZE=$(stat -c %s "$X265_FILE")
fi

if [ $MAX_QUALITY -eq 1 ]; then
  X265_FILE="$BASENAME-$MAX_QUALITY.x265.tmp"
  $X265 -r "$X265_FILE.y4m" $(echo $QSTR | sed 's/\$x/'$MAX_QUALITY'/g') --keyint $KEYINT -o "$X265_FILE" $FILE 2> /dev/null > /dev/null
  MAX_QUALITY_SIZE=$(stat -c %s "$X265_FILE")
fi

if (( $MAX_QUALITY_SIZE - $SIZE < $SIZE - $MIN_QUALITY_SIZE )); then
  BEST_QUALITY=$MAX_QUALITY
else
  BEST_QUALITY=$MIN_QUALITY
fi

mv $BASENAME-$BEST_QUALITY.x265.tmp $BASENAME-$BEST_QUALITY.x265

if [ $FRAMES -eq 1 ]; then
  $Y4M2PNG -o $BASENAME-$BEST_QUALITY.x265.png $BASENAME-$BEST_QUALITY.x265.tmp.y4m
else
  mv $BASENAME-$BEST_QUALITY.x265.tmp.y4m $BASENAME-$BEST_QUALITY.x265.y4m
fi

rm $BASENAME-*.x265.tmp $BASENAME-*.x265.tmp.y4m
