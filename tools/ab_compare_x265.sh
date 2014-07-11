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
  for i in {51..1}; do
    X265_FILE=$BASENAME-$i.x265.tmp
    $X265 -r $X265_FILE.y4m $(echo $QSTR | sed 's/\$x/'$i'/g') -o $X265_FILE $FILE 2> /dev/null > /dev/null
    X265_SIZE=$(stat -c %s $X265_FILE)
    if (($X265_SIZE > $OGV_SIZE)); then
      if [[ -z $X265_LAST_SIZE ]]; then
        mv $X265_FILE $BASENAME-$i.x265
        mv $X265_FILE.y4m $BASENAME-$i.x265.y4m
      else
        if (($X265_SIZE - $OGV_SIZE < $OGV_SIZE - $X265_LAST_SIZE)); then
          mv $X265_FILE $BASENAME-$i.x265
          mv $X265_FILE.y4m $BASENAME-$i.x265.y4m
        else
          mv $X265_LAST_FILE $BASENAME-$i.x265
          mv $X265_LAST_FILE.y4m $BASENAME-$i.x265.y4m
        fi
      fi
      $Y4M2PNG -o $BASENAME-$i.x265.png $BASENAME-$i.x265.y4m
      rm $BASENAME-$i.x265.y4m
      break
    fi
    X265_LAST_SIZE=$X265_SIZE
    X265_LAST_FILE=$X265_FILE
  done
  rm $BASENAME-*.x265.tmp $BASENAME-*.x265.tmp.y4m $BASENAME-$V.ogv.y4m
done
