#!/bin/bash
set -e

# TODO - add a test for this
export CODEC=$1
shift

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $LIBVPX_ROOT ] || [ ! -d $LIBVPX_ROOT ]; then
  echo "Please set LIBVPX_ROOT to the location of your libvpx git clone"
  exit 1
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
fi

if [ -z "$DUMP_VIDEO" ]; then
  DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
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

if [ ! -x "$Y4M2PNG" ]; then
  echo "Executable not found Y4M2PNG=$Y4M2PNG"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$VPXENC" ]; then
  echo "Executable not found VPXENC=$VPXENC"
  echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
  exit 1
fi

if [ ! -x "$VPXDEC" ]; then
  echo "Executable not found VPXDEC=$VPXDEC"
  echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
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
  # With libvpx, the lowest quantizer number yields the highest quality and
  # vice versa. Here, MAX_QUALITY produces the best looking image, so it's the
  # lowest number.
  MAX_QUALITY=3
  MIN_QUALITY=63
  while (( $MIN_QUALITY - $MAX_QUALITY > 1 )); do
    QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
    VPX_FILE=$BASENAME-$QUALITY.$CODEC.tmp
    $VPXENC --codec=$CODEC --good --cpu-used=0 -y --min-q=$QUALITY --max-q=$QUALITY -o $VPX_FILE $FILE 2> /dev/null
    VPX_SIZE=$(stat -c %s $VPX_FILE)
    if (($VPX_SIZE > $OGV_SIZE)); then
      MAX_QUALITY=$QUALITY
      MAX_QUALITY_SIZE=$VPX_SIZE
    else
      MIN_QUALITY=$QUALITY
      MIN_QUALITY_SIZE=$VPX_SIZE
    fi
  done

  if [ $MIN_QUALITY -eq 63 ]; then
    $VPXENC --codec=$CODEC --good --cpu-used=0 -y --min-q=$MIN_QUALITY --max-q=$MIN_QUALITY -o $BASENAME-$MIN_QUALITY.$CODEC.tmp $FILE 2> /dev/null
    MIN_QUALITY_SIZE=$(stat -c %s $BASENAME-$MIN_QUALITY.$CODEC.tmp)
  fi

  if [ $MAX_QUALITY -eq 3 ]; then
    $VPXENC --codec=$CODEC --good --cpu-used=0 -y --min-q=$MAX_QUALITY --max-q=$MAX_QUALITY -o $BASENAME-$MAX_QUALITY.$CODEC.tmp $FILE 2> /dev/null
    MAX_QUALITY_SIZE=$(stat -c %s $BASENAME-$MAX_QUALITY.$CODEC.tmp)
  fi

  if (( $MAX_QUALITY_SIZE - $OGV_SIZE < $OGV_SIZE - $MIN_QUALITY_SIZE )); then
    VPX_SIZE=$MAX_QUALITY_SIZE
    BEST_QUALITY=$MAX_QUALITY
  else
    BEST_QUALITY=$MIN_QUALITY
    VPX_SIZE=$MIN_QUALITY_SIZE
  fi

  BEST_FILE=$BASENAME-$BEST_QUALITY.$CODEC
  $VPXDEC --codec=$CODEC -o $BEST_FILE.y4m $BEST_FILE.tmp
  $Y4M2PNG -o $BEST_FILE.png $BEST_FILE.y4m
  mv $BEST_FILE.tmp $BEST_FILE
  rm $BEST_FILE.y4m $BASENAME-*.$CODEC.tmp $BASENAME-$V.ogv.y4m
done
