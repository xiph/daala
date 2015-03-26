#!/bin/bash
set -e

if [ $# == 0 ]; then
  echo "usage: DAALA_ROOT=<build_dir> $0 *.y4m"
  exit 1
fi

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $CHECKSUM]; then
  CHECKSUM=md5sum
fi

if [ -z "$ENCODER_EXAMPLE" ]; then
  ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
fi

if [ ! -x "$ENCODER_EXAMPLE" ]; then
  echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

if [ -z "$DUMP_VIDEO" ]; then
  export DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
fi

if [ ! -x "$DUMP_VIDEO" ]; then
  echo "Executable not found DUMP_VIDEO=$DUMP_VIDEO"
  echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
  exit 1
fi

for FILE in "$@"; do
  BASENAME=$(basename $FILE)
  rm $BASENAME.out 2> /dev/null || true
  echo $BASENAME
  $ENCODER_EXAMPLE -v 0 $FILE -o $BASENAME.ogv
  $DUMP_VIDEO $BASENAME.ogv -o $BASENAME_out.y4m
  OUTSUM=`tail -n+3 $BASENAME_out.y4m | $CHECKSUM`
  INSUM=`tail -n+3 $FILE | $CHECKSUM`
  if [ "$INSUM" == "$OUTSUM" ]
  then
    echo "Lossless matches input for $FILE"
  else
    echo "Lossless MISMATCHES input for $FILE"
    # Skip deletion of generated files for debugging
    exit 1
  fi
  rm $BASENAME_out.y4m $BASENAME.ogv > /dev/null || true
done
