#!/bin/bash
set -e

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
fi

if [ -z "$DUMP_VIDEO" ]; then
  DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
fi

if [ -z "$YUVJPEG" ]; then
  YUVJPEG=$DAALA_ROOT/tools/yuvjpeg
fi

if [ -z "$JPEGYUV" ]; then
  JPEGYUV=$DAALA_ROOT/tools/jpegyuv
fi

if [ -z "$YUV2YUV4MPEG" ]; then
  YUV2YUV4MPEG=$DAALA_ROOT/tools/yuv2yuv4mpeg
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

if [ ! -x "$YUVJPEG" ]; then
  echo "Executable not found YUVJPEG=$YUVJPEG"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$JPEGYUV" ]; then
  echo "Executable not found JPEGYUV=$JPEGYUV"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ ! -x "$YUV2YUV4MPEG" ]; then
  echo "Executable not found YUV2YUV4MPEG=$YUV2YUV4MPEG"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
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
  tail -n+3 $FILE > $BASENAME-in.yuv
  WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
  HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
  MAX_QUALITY=100
  MIN_QUALITY=0
  while (( $MAX_QUALITY - $MIN_QUALITY > 1 )); do
    QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
    JPEG_FILE=$BASENAME-$QUALITY.jpeg.tmp
    $YUVJPEG $QUALITY "$WIDTH"x$HEIGHT $BASENAME-in.yuv $JPEG_FILE
    JPEG_SIZE=$(stat -c %s $JPEG_FILE)
    if (($JPEG_SIZE > $OGV_SIZE)); then
      MAX_QUALITY=$QUALITY
      MAX_QUALITY_SIZE=$JPEG_SIZE
    else
      MIN_QUALITY=$QUALITY
      MIN_QUALITY_SIZE=$JPEG_SIZE
    fi
  done

  if [ $MIN_QUALITY -eq 0 ]; then
    $YUVJPEG $MIN_QUALITY "$WIDTH"x$HEIGHT $BASENAME-in.yuv $BASENAME-$MIN_QUALITY.jpeg.tmp
    MIN_QUALITY_SIZE=$(stat -c %s $BASENAME-$MIN_QUALITY.jpeg.tmp)
  fi

  if [ $MAX_QUALITY -eq 100 ]; then
    $YUVJPEG $MAX_QUALITY "$WIDTH"x$HEIGHT $BASENAME-in.yuv $BASENAME-$MAX_QUALITY.jpeg.tmp
    MAX_QUALITY_SIZE=$(stat -c %s $BASENAME-$MAX_QUALITY.jpeg.tmp)
  fi

  if (( $MAX_QUALITY_SIZE - $OGV_SIZE < $OGV_SIZE - $MIN_QUALITY_SIZE )); then
    BEST_QUALITY=$MAX_QUALITY
  else
    BEST_QUALITY=$MIN_QUALITY
  fi

  BEST_FILE=$BASENAME-$BEST_QUALITY.jpeg
  mv $BEST_FILE.tmp $BEST_FILE
  $JPEGYUV $BEST_FILE $BEST_FILE.yuv
  $YUV2YUV4MPEG $BEST_FILE -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  $Y4M2PNG -o $BEST_FILE.png $BEST_FILE.y4m
  rm $BEST_FILE.yuv $BASENAME-*.jpeg.tmp
  rm $BEST_FILE.y4m $BASENAME-in.yuv $BASENAME-$V.ogv.y4m
done
