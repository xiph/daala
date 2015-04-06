#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

./autogen.sh
CC=clang ./configure
make clean
make
