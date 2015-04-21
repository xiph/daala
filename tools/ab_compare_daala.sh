#!/bin/bash
set -e

while getopts 's:v:bk:d:E:D:y:n:' OPTIONS; do
  case $OPTIONS in
    s) SIZE="$OPTARG";;
    v) V="$OPTARG";;
    b) B_MODE='true';;
    k) KEYINT="$OPTARG";;
    # -b is whether or not we are running this script to produce images that
    # are for the secondary Daala in an a-b comparison, as opposed to the
    # first.
    d) DAALA_ROOT="$OPTARG";;
    E) ENCODER_EXAMPLE="$OPTARG";;
    D) DUMP_VIDEO="$OPTARG";;
    Y) Y4M2PNG="$OPTARG";;
    n) FRAMES="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [[ -n "$V" && -n "$SIZE" ]]; then
  echo "Both quality setting and size specified."
  exit 1
fi

if [[ -z "$V" && -z "$SIZE" ]]; then
  echo "Neither quality setting and size specified."
  exit 1
fi

if [ -z "$DAALA_ROOT" ]; then
  DAALA_ROOT=.
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

if [ ! -x "$ENCODER_EXAMPLE" ]; then
  echo "Example encoder not found at '$ENCODER_EXAMPLE.'"
  exit 1
fi

if [ ! -x "$DUMP_VIDEO" ]; then
  echo "Video dumper not found at '$DUMP_VIDEO'."
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
BASENAME=$(basename "$FILE")
if [ -n "$V" ]; then
  echo "$FILE"
  $ENCODER_EXAMPLE -v $V "$FILE" -o "$BASENAME-$V.ogv" 2> /dev/null
  $DUMP_VIDEO -o "$BASENAME-$V.ogv.y4m" "$BASENAME-$V.ogv" 2> /dev/null

  if [ $FRAMES -eq 1 ]; then
    $Y4M2PNG -o "$BASENAME-$V.ogv.png" "$BASENAME-$V.ogv.y4m"
    rm "$BASENAME-$V.ogv.y4m"
  fi

else
  # With Daala, the lowest quantizer number yields the highest quality and
  # vice versa. Here, MAX_QUALITY produces the best looking image, so it's the
  # lowest number.
  MAX_QUALITY=0
  MIN_QUALITY=511
  while (( $MIN_QUALITY - $MAX_QUALITY > 1 )); do
    QUALITY=$(( ($MIN_QUALITY + $MAX_QUALITY) / 2 ))
    if [ "$B_MODE" == 'true' ]; then
      OGV_FILE="$BASENAME-b-$QUALITY.ogv.tmp"
    else
      OGV_FILE="$BASENAME-$QUALITY.ogv.tmp"
    fi
    $ENCODER_EXAMPLE -v $QUALITY "$FILE" -k $KEYINT -o "$OGV_FILE" 2> /dev/null
    OGV_SIZE=$(stat -c %s "$OGV_FILE")
    if (($OGV_SIZE > $SIZE)); then
      MAX_QUALITY=$QUALITY
      MAX_QUALITY_SIZE=$OGV_SIZE
    else
      MIN_QUALITY=$QUALITY
      MIN_QUALITY_SIZE=$OGV_SIZE
    fi
  done

  if [ $MIN_QUALITY -eq 511 ]; then
    if [ "$B_MODE" == 'true' ]; then
      FILENAME="$BASENAME-b-$MIN_QUALITY.ogv.tmp"
    else
      FILENAME="$BASENAME-$MIN_QUALITY.ogv.tmp"
    fi
    $ENCODER_EXAMPLE -v $MIN_QUALITY "$FILE" -k $KEYINT -o "$FILENAME" 2> /dev/null
    MIN_QUALITY_SIZE=$(stat -c %s "$FILENAME")
  fi

  if [ $MAX_QUALITY -eq 0 ]; then
    if [ "$B_MODE" == 'true' ]; then
      FILENAME="$BASENAME-b-$MAX_QUALITY.ogv.tmp"
    else
      FILENAME="$BASENAME-$MAX_QUALITY.ogv.tmp"
    fi
    $ENCODER_EXAMPLE -v $QUALITY "$FILE" -k $KEYINT -o "$FILENAME" 2> /dev/null
    MAX_QUALITY_SIZE=$(stat -c %s "$FILENAME")
  fi

  if (( $MAX_QUALITY_SIZE - $SIZE < $SIZE - $MIN_QUALITY_SIZE )); then
    VPX_SIZE=$MAX_QUALITY_SIZE
    BEST_QUALITY=$MAX_QUALITY
  else
    BEST_QUALITY=$MIN_QUALITY
    VPX_SIZE=$MIN_QUALITY_SIZE
  fi

  if [ "$B_MODE" == 'true' ]; then
   BEST_FILE="$BASENAME-b-$BEST_QUALITY.ogv"
  else
   BEST_FILE="$BASENAME-$BEST_QUALITY.ogv"
  fi

  mv "$BEST_FILE.tmp" "$BEST_FILE"
  $DUMP_VIDEO -o "$BEST_FILE.y4m" "$BEST_FILE" 2>/dev/null
  if [ $FRAMES -eq 1 ]; then
    $Y4M2PNG -o "$BEST_FILE.png" "$BEST_FILE.y4m"
    rm "$BASENAME"-*.ogv.tmp "$BEST_FILE.y4m"
  else
    rm "$BASENAME"-*.ogv.tmp
  fi
fi
