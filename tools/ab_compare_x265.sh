#!/bin/bash
set -e

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $X265_ROOT ] || [ ! -d $X265_ROOT ]; then
  echo "Please set X265_ROOT to the location of your x265 hg checkout"
  exit 1
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
fi

if [ -z "$DUMP_VIDEO" ]; then
  DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
fi

if [ -z "$X265" ]; then
  export X265=$X265_ROOT/build/linux/x265
fi

if [ -z "$Y4M2PNG" ]; then
  Y4M2PNG=$DAALA_ROOT/tools/y4m2png
fi

if [ ! -x "$ENCODER_EXAMPLE" ]; then
  echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_VIDEO" ]; then
  echo "Executable not found DUMP_VIDEO=$DUMP_VIDEO"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$X265" ]; then
  echo "Executable not found X265=$X265"
  echo "Do you have the right X265_ROOT=$X265_ROOT"
  exit 1
fi

if [ ! -x "$Y4M2PNG" ]; then
  echo "Executable not found Y4M2PNG=$Y4M2PNG"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ -z $V ]; then
  V=29
fi

for FILE in $@; do
  echo $FILE
  BASENAME=$(basename $FILE)
  $ENCODER_EXAMPLE -v $V $FILE -o $BASENAME-$V.ogv 2> /dev/null
  $DUMP_VIDEO -o $BASENAME-$V.ogv.y4m $BASENAME-$V.ogv 2> /dev/null
  $Y4M2PNG -o $BASENAME-$V.ogv.png $BASENAME-$V.ogv.y4m
  OGV_SIZE=$(stat -c %s $BASENAME-$V.ogv)
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
    $X265 -r $X265_FILE.y4m $(echo $QSTR | sed 's/\$x/'$QUALITY'/g') -o $X265_FILE $FILE 2> /dev/null > /dev/null
    X265_SIZE=$(stat -c %s $X265_FILE)
    if (($X265_SIZE > $OGV_SIZE)); then
      MAX_QUALITY=$QUALITY
      MAX_QUALITY_SIZE=$X265_SIZE
    else
      MIN_QUALITY=$QUALITY
      MIN_QUALITY_SIZE=$X265_SIZE
    fi
  done

  if [ $MIN_QUALITY -eq 51 ]; then
    $X265 -r $BASENAME-$MIN_QUALITY.x265.tmp.y4m $(echo $QSTR | sed 's/\$x/'$MIN_QUALITY'/g') -o $BASENAME-$MIN_QUALITY.x265.tmp $FILE 2> /dev/null > /dev/null
    MIN_QUALITY_SIZE=$(stat -c %s $BASENAME-$MIN_QUALITY.tmp)
  fi

  if [ $MAX_QUALITY -eq 1 ]; then
    $X265 -r $BASENAME-$MAX_QUALITY.x265.tmp.y4m $(echo $QSTR | sed 's/\$x/'$MAX_QUALITY'/g') -o $X265_FILE $FILE 2> /dev/null > /dev/null
    MAX_QUALITY_SIZE=$(stat -c %s $BASENAME-$MAX_QUALITY.tmp)
  fi

  if (( $MAX_QUALITY_SIZE - $OGV_SIZE < $OGV_SIZE - $MIN_QUALITY_SIZE )); then
    BEST_QUALITY=$MAX_QUALITY
  else
    BEST_QUALITY=$MIN_QUALITY
  fi

  mv $BASENAME-$BEST_QUALITY.x265.tmp $BASENAME-$BEST_QUALITY.x265
  $Y4M2PNG -o $BASENAME-$BEST_QUALITY.x265.png $BASENAME-$BEST_QUALITY.x265.tmp.y4m
  rm $BASENAME-$BEST_QUALITY.x265.tmp.y4m
  rm $BASENAME-*.x265.tmp $BASENAME-*.x265.tmp.y4m $BASENAME-$V.ogv.y4m
done
