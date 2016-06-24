#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

if [ -z "$VIDEOS"]; then
  VIDEOS=/usr/local/share/videos
fi

./autogen.sh
CFLAGS='-O2 -g' ./configure --enable-assertions --enable-check-asm --enable-logging --enable-accounting
make clean
make distcheck
make docs
make
./examples/encoder_example -k 4 ${VIDEOS}/claire_qcif-2frames.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
./examples/encoder_example -k 4 ${VIDEOS}/tos444.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
