#!/bin/bash
# continuous integration test script
# run this from the top-level source directory

LCOV_PATH=/opt/lcov/bin
LCOV=${LCOV_PATH}/lcov
GENHTML=${LCOV_PATH}/genhtml

cd unix
rm -Rf objs
make clean
CFLAGS="-g3 -fprofile-arcs -ftest-coverage -UOD_ENABLE_ASSERTIONS" make
CFLAGS="-g3 -fprofile-arcs -ftest-coverage -UOD_ENABLE_ASSERTIONS" make check
${LCOV} -c -b `pwd` -d `pwd` -o makecheck.info
${GENHTML} -s -o coverage/ makecheck.info
