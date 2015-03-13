#!/bin/bash
set -e

while getopts 's:d:J:j:Y:y:' OPTIONS; do
  case $OPTIONS in
    s) SIZE="$OPTARG";;
    d) DAALA_ROOT="$OPTARG";;
    J) YUVJPEG="$OPTARG";;
    j) JPEGYUV="$OPTARG";;
    Y) Y4M2PNG="$OPTARG";;
    y) YUV2YUV4MPEG="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [ -z "$SIZE" ]; then
  echo "No file size given."
  exit 1
fi

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
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

if [ ! -x "$YUVJPEG" ]; then
  echo "YUVJPEG not found at '$YUVJPEG'."
#  exit 1
fi

if [ ! -x "$JPEGYUV" ]; then
  echo "JPEGYUV not found at '$JPEGYUV'."
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

for FILE in $@; do
  echo $FILE
  BASENAME=$(basename $FILE)
  tail -n+3 $FILE > $BASENAME-in.yuv
  WIDTH=$(head -1 $FILE | cut -d\  -f 2 | tr -d 'W')
  HEIGHT=$(head -1 $FILE | cut -d\  -f 3 | tr -d 'H')
  # With JPEG, the highest quality number creates the best-looking image.
  MAX_QUALITY=100
  MIN_QUALITY=0
  while (( $MAX_QUALITY - $MIN_QUALITY > 1 )); do
    QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
    JPEG_FILE=$BASENAME-$QUALITY.jpeg.tmp
    $YUVJPEG $QUALITY "$WIDTH"x$HEIGHT $BASENAME-in.yuv $JPEG_FILE
    JPEG_SIZE=$(stat -c %s $JPEG_FILE)
    if (($JPEG_SIZE > $SIZE)); then
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

  if (( $MAX_QUALITY_SIZE - $SIZE < $SIZE - $MIN_QUALITY_SIZE )); then
    BEST_QUALITY=$MAX_QUALITY
  else
    BEST_QUALITY=$MIN_QUALITY
  fi

  BEST_FILE=$BASENAME-$BEST_QUALITY.jpeg
  mv $BEST_FILE.tmp $BEST_FILE
  $JPEGYUV $BEST_FILE $BEST_FILE.yuv
  $YUV2YUV4MPEG $BEST_FILE -w$WIDTH -h$HEIGHT -an0 -ad0 -c420mpeg2
  $Y4M2PNG -o $BEST_FILE.png $BEST_FILE.y4m
  rm $BASENAME-in.yuv  $BASENAME-*.jpeg.tmp
  rm $BEST_FILE.yuv $BEST_FILE.y4m
done
