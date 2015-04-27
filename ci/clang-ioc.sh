#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

./autogen.sh
CC=clang CFLAGS='-Wall -O3 -g -fsanitize=undefined -fno-sanitize-recover' LDFLAGS='-fsanitize=undefined -fno-sanitize-recover' ./configure
make clean all
./examples/encoder_example -k 4 ${VIDEOS}/claire_qcif-2frames.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
