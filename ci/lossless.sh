#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

VIDEOS=/usr/local/share/videos

./autogen.sh
CFLAGS='-O2 -g' ./configure --enable-assertions --enable-check-asm --enable-logging --enable-accounting
make clean all
DAALA_ROOT=. tools/check_lossless.sh ${VIDEOS}/claire_qcif-2frames.y4m
DAALA_ROOT=. tools/check_lossless.sh ${VIDEOS}/tos444.y4m
