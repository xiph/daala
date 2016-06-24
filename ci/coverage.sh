#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

# In Jenkins scripts like this should generally be sourced so that they
#  run in the jenkins controlled shell environment so that errors get noticed.

# run the lcov code-coverage analysis tool

LCOV_PATH=/opt/lcov/bin
LCOV=${LCOV_PATH}/lcov
GENHTML=${LCOV_PATH}/genhtml
if [ -z "$VIDEOS" ]; then
  VIDEOS=/usr/local/share/videos
fi

cd unix && rm -Rf objs
make clean
CFLAGS="-g3 -fprofile-arcs -ftest-coverage -UOD_ENABLE_ASSERTIONS" make
${LCOV} -z -d `pwd` -b `pwd`
${LCOV} -c -i -b `pwd` -d `pwd` -t baseline -o baseline.info
CFLAGS="-g3 -fprofile-arcs -ftest-coverage -UOD_ENABLE_ASSERTIONS" make check
${LCOV} -c -b `pwd` -d `pwd` -t check -o makecheck.info
${LCOV} -z -d `pwd` -b `pwd`
rm -f out.ogv
./encoder_example -k 4 ${VIDEOS}/claire_qcif-2frames.y4m -o out.ogv
./dump_video out.ogv -o /dev/null
${LCOV} -c -b `pwd` -d `pwd` -t encoderdecoder -o encdec.info
${LCOV} -a baseline.info -a makecheck.info -a encdec.info -o daala_coverage.info
${GENHTML} -s -o coverage/ daala_coverage.info
