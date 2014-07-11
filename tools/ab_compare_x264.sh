#!/bin/bash
set -e

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $X264_ROOT ] || [ ! -d $X264_ROOT ]; then
  echo "Please set X264_ROOT to the location of your x264 git clone"
  exit 1
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
fi

if [ -z "$DUMP_VIDEO" ]; then
  DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
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

if [ ! -x "$X264" ]; then
  echo "Executable not found X264=$X264"
  echo "Do you have the right X264_ROOT=$X264_ROOT"
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
  WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
  HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
  QSTR="--preset placebo --crf=\$x"
  for i in {51..1}; do
    X264_FILE=$BASENAME-$i.x264.tmp
    $X264 --dump-yuv $X264_FILE.yuv $(echo $QSTR | sed 's/\$x/'$i'/g') -o $X264_FILE $FILE 2> /dev/null > /dev/null
    X264_SIZE=$(stat -c %s $X264_FILE)
    if (($X264_SIZE > $OGV_SIZE)); then
      if [[ -z $X264_LAST_SIZE ]]; then
        mv $X264_FILE $BASENAME-$i.x264
        mv $X264_FILE.yuv $BASENAME-$i.x264.yuv
      else
        if (($X264_SIZE - $OGV_SIZE < $OGV_SIZE - $X264_LAST_SIZE)); then
          mv $X264_FILE $BASENAME-$i.x264
          mv $X264_FILE.yuv $BASENAME-$i.x264.yuv
        else
          mv $X264_LAST_FILE $BASENAME-$i.x264
          mv $X264_LAST_FILE.yuv $BASENAME-$i.x264.yuv
        fi
      fi
      $YUV2YUV4MPEG $BASENAME-$i.x264 -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
      $Y4M2PNG -o $BASENAME-$i.x264.png $BASENAME-$i.x264.y4m
      rm $BASENAME-$i.x264.yuv $BASENAME-$i.x264.y4m
      break
    fi
    X264_LAST_SIZE=$X264_SIZE
    X264_LAST_FILE=$X264_FILE
  done
  rm $BASENAME-*.x264.tmp $BASENAME-*.x264.tmp.yuv $BASENAME-$V.ogv.y4m
done
