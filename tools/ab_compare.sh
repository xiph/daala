#!/bin/sh
set -e

USAGE=\
'This script creates images or videos that are compressed by different codecs
but have nearly the same size, facilitating comparison among codecs.

Usage:
ab_compare.sh [OPTIONS] file1.y4m file2.y4m ... fileN.y4m

Options:
-c codec1,codec2,codec3 ...
Specify which codecs to compare. Availaible codecs are daala, daala2 (a
second build of Daala, useful for comparing two versions of Daala against
each other), vp8, vp9, x264, x265, and jpeg.

-v quality
Match the images compressed by other codecs to the size of the Daala file
compressed with that quality.

-b bpp
Compress the images at this many bits per pixel.

-s size
Compress the images to this size in bytes.

-k interval
Set the interval between keyframes. If not given, defaults to 256.

-d daala_root
The root directory of the main build of Daala. If not given, assumed to be ./ .

-D daala2_root
The root directory of a secondary build of Daala. Only needed if
using a second version of daala.

-l libvpx_root
A directory containing a build of libvpx. Only needed if using vp8 or vp9.

-x x264_root
A directory containing a build of x264. Only needed if using x264.

-X x265_root
A directory containing a build of x265. Only needed if using x265.

-a ab_root
Specify the directory containing the other ab_compare scripts. If not given,
assumed to be ./tools/ .

-h
Print this help text.'

while getopts 'c:s:v:b:k:ha:d:D:l:x:X:' OPTIONS; do
  case $OPTIONS in
    c) CODECS="$OPTARG";;
    s) SIZE="$OPTARG";;
    v) V="$OPTARG";;
    b) BPP="$OPTARG";;
    k) KEYINT="$OPTARG";;
    h) echo "$USAGE"; exit 0;;
    a) AB_ROOT="$OPTARG";;
    d) DAALA_ROOT="$OPTARG";;
    D) DAALA2_ROOT="$OPTARG";;
    l) LIBVPX_ROOT="$OPTARG";;
    x) X264_ROOT="$OPTARG";;
    X) X265_ROOT="$OPTARG";;
  esac
done
shift $(($OPTIND - 1))

if [ -z "$1" ]; then
  echo "No input files given."
  exit 0
fi

if [ -z "$AB_ROOT" ]; then
  AB_ROOT="./tools/"
fi

if [ -z "$DAALA_ROOT" ]; then
  DAALA_ROOT='.'
fi

if [ ! -x "$AB_ROOT" ]; then
  echo "Comparison scripts not found at '$AB_ROOT'. Please give their location"
  echo "with -a."
  exit 1
fi

if [ ! -x "$LIBVPX_ROOT" ] && echo "$CODECS" | grep 'vp[89]' > /dev/null; then
  echo "Libvpx build not found at '$LIBVPX_ROOT'. Please give its location"
  echo "with -l."
  exit 1
fi

if [ ! -x "$X264_ROOT" ] && echo "$CODECS" | grep 'x264' > /dev/null; then
  echo "x264 build not found at $X264_ROOT. Please give its location"
  echo "with -x."
  exit 1
fi

if [ ! -x "$X265_ROOT" ] && echo "$CODECS" | grep 'x265' > /dev/null; then
  echo "x265 build not found at $X265_ROOT. Please give its location"
  echo "with -X."
  exit 1
fi

if [ ! -x "$DAALA_ROOT" ] && echo "$CODECS" | grep 'daala$\|daala,' > /dev/null; then
  echo "Daala build not found at $DAALA_ROOT. Please give its location"
  echo "with -d."
  exit 1
fi

if [ ! -x "$DAALA2_ROOT" ] && echo "$CODECS" | grep 'daala2' > /dev/null; then
  echo "Second Daala build not found at $DAALA2_ROOT. Please give its location"
  echo "with -D."
  exit 1
fi

if [ -z "$KEYINT" ]; then
  KEYINT=256;
fi

if [ -z "$CODECS" ]; then
  CODECS="daala vp8 x264 x265"
else
  CODECS=$( echo $CODECS | tr '[,]' '[ ]')
fi

if echo "$CODECS" | grep 'daala$\|daala ' > /dev/null; then
  USE_DAALA="true"
else
  USE_DAALA="false"
fi

if [[ -n "$V" && "$USE_DAALA" == "false" ]]; then
  echo "If you give a Daala quality setting, Daala must be one of the
codecs used in the comparison"
  exit 1
fi

if [[ -z "$V" && -z "$SIZE" && -z "$BPP" ]]; then
  echo "No method of determining file size specified."
  exit 1
fi

if [[ -n "$V" && (-n "$SIZE" || -n "$BPP") || -n "$SIZE" && -n "$BPP" ]]; then
  echo "More than one way of determining file size given."
  exit 1
fi

for FILE in $@; do
  HEADER=$(head -1 "$FILE")
  WIDTH=$(echo "$HEADER" | cut -d\  -f 2 | tr -d 'W')
  HEIGHT=$(echo "$HEADER" | cut -d\  -f 3 | tr -d 'H')

  # CHROMA_FACTOR is the ratio of the size of all the image data to the size of
  # the luma plane's data. It's 1.5 for 4:2:0, 2 for 4:2:2 and 3 for 4:4:4. We
  # need to know this number in order to calculate the number of frames from
  # the file size.
  # If the header says nothing about subsampling, it's 4:2:0.
  if ! echo "$HEADER" | grep ' C' > /dev/null; then
    CHROMA_FACTOR=1.5
  else
    if echo "$HEADER" | grep ' C420' > /dev/null; then
      CHROMA_FACTOR=1.5
    elif echo "$HEADER" | grep ' C422' > /dev/null; then
      CHROMA_FACTOR=2.0
    elif echo "$HEADER" | grep ' C444' > /dev/null; then
      CHROMA_FACTOR=3.0
    else
      echo "Chroma subsampling unknown, assuming 4:2:0."
      CHROMA_FACTOR=1.5
    fi
  fi

  INSIZE=$(ls -l "$FILE" | awk '{print $5}')
  # This ignores the size of the headers, but should be good enough for most
  # cases.
  FRAMES=$( echo "scale=5; $INSIZE / ( $WIDTH * $HEIGHT * $CHROMA_FACTOR )" | bc )
  FRAMES=$(echo "scale=0;$FRAMES/1" | bc)

  if [ -n "$BPP" ]; then
  SIZE=$(echo "scale=5;( $BPP * $HEIGHT * $WIDTH * $FRAMES) / 8 + 0.5" | bc)
  SIZE=$(echo "scale=0;$SIZE/1" | bc)
  fi

  if [ $USE_DAALA == "true" ]; then
    if [ -n "$V" ]; then
      $AB_ROOT/ab_compare_daala.sh -d "$DAALA_ROOT" -n "$FRAMES" -v "$V" "$FILE"
      SIZE=$(wc -c $(basename "$FILE")-$V.ogv | awk '{ print $1 }')
    else
      $AB_ROOT/ab_compare_daala.sh -d "$DAALA_ROOT" -n "$FRAMES" -s "$SIZE" "$FILE"
    fi
  fi

  for CODEC in $CODECS; do
    case $CODEC in
      daala) : ;;
      jpeg) $AB_ROOT/ab_compare_jpeg.sh -d "$DAALA_ROOT" -n "$FRAMES" -s "$SIZE" "$FILE";;
      vp8) $AB_ROOT/ab_compare_libvpx.sh -d "$DAALA_ROOT" -r "$LIBVPX_ROOT" -n "$FRAMES" -k "$KEYINT" -s "$SIZE" -c vp8 "$FILE";;
      vp9) $AB_ROOT/ab_compare_libvpx.sh -d "$DAALA_ROOT" -r "$LIBVPX_ROOT" -n "$FRAMES" -k "$KEYINT" -s "$SIZE" -c vp9 "$FILE";;
      x264) $AB_ROOT/ab_compare_x264.sh -d "$DAALA_ROOT" -r "$X264_ROOT" -n "$FRAMES" -k "$KEYINT" -s "$SIZE" "$FILE";;
      x265) $AB_ROOT/ab_compare_x265.sh -d "$DAALA_ROOT" -r "$X265_ROOT" -n "$FRAMES" -k "$KEYINT" -s "$SIZE" "$FILE";;
      daala2) $AB_ROOT/ab_compare_daala.sh -d "$DAALA_ROOT" -E "$DAALA2_ROOT/examples/encoder_example" -D "$DAALA2_ROOT/examples/dump_video" -n "$FRAMES" -s "$SIZE" -b "$FILE";;
      *) echo "Unknown codec: $CODEC"
    esac
  done
done
