#!/bin/sh

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
  for i in $(seq 0 100); do
    JPEG_FILE=$BASENAME-$i.jpeg.tmp
    $YUVJPEG $i "$WIDTH"x$HEIGHT $BASENAME-in.yuv $JPEG_FILE
    JPEG_SIZE=$(stat -c %s $BASENAME-$i.jpeg.tmp)
    if (($JPEG_SIZE > $OGV_SIZE)); then
      if [[ -z $JPEG_LAST_SIZE ]]; then
        mv $JPEG_FILE $BASENAME-$i.jpeg
      else
        if (($JPEG_SIZE - $OGV_SIZE < $OGV_SIZE - $JPEG_LAST_SIZE)); then
          mv $JPEG_FILE $BASENAME-$i.jpeg
        else
          mv $JPEG_LAST_FILE $BASENAME-$i.jpeg
        fi
      fi
      $JPEGYUV $BASENAME-$i.jpeg $BASENAME-$i.jpeg.yuv
      $YUV2YUV4MPEG $BASENAME-$i.jpeg -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
      $Y4M2PNG -o $BASENAME-$i.jpeg.png $BASENAME-$i.jpeg.y4m
      rm $BASENAME-$i.jpeg.yuv $BASENAME-$i.jpeg.y4m
      break
    fi
    JPEG_LAST_SIZE=$JPEG_SIZE
    JPEG_LAST_FILE=$JPEG_FILE
  done
  rm $BASENAME-*.jpeg.tmp $BASENAME-in.yuv $BASENAME-$V.ogv.y4m
done
