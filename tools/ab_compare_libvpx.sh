#!/bin/sh

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
  for i in {63..3}; do
    VPX_FILE=$BASENAME-$i.$CODEC.tmp
    $VPXENC --codec=$CODEC --good --cpu-used=0 -y --min-q=$i --max-q=$i -o $VPX_FILE $FILE 2> /dev/null
    VPX_SIZE=$(stat -c %s $VPX_FILE)
    if (($VPX_SIZE > $OGV_SIZE)); then
      if [[ -z $VPX_LAST_SIZE ]]; then
        mv $VPX_FILE $BASENAME-$i.$CODEC
      else
        if (($VPX_SIZE - $OGV_SIZE < $OGV_SIZE - $VPX_LAST_SIZE)); then
          mv $VPX_FILE $BASENAME-$i.$CODEC
        else
          mv $VPX_LAST_FILE $BASENAME-$i.$CODEC
        fi
      fi
      $VPXDEC --codec=$CODEC -o $BASENAME-$i.$CODEC.y4m $BASENAME-$i.$CODEC
      $Y4M2PNG -o $BASENAME-$i.$CODEC.png $BASENAME-$i.$CODEC.y4m
      rm $BASENAME-$i.$CODEC.y4m
      break
    fi
    VPX_LAST_SIZE=$VPX_SIZE
    VPX_LAST_FILE=$VPX_FILE
  done
  rm $BASENAME-*.$CODEC.tmp $BASENAME-$V.ogv.y4m
done
