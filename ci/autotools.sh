#!/bin/bash
# continuous integration test script
# run this from the top-level source directory

OGG_PATH=/srv/jenkins/jobs/libogg/workspace
VIDEOS=/usr/local/share/videos

./autogen.sh
CFLAGS='-O2 -DOD_CHECKASM -g' ./configure --enable-assertions --enable-logging PKG_CONFIG_PATH=${OGG_PATH}
make distcheck PKG_CONFIG_PATH=${OGG_PATH}
make docs
./examples/encoder_example -k 4 ${VIDEOS}/claire_qcif-2frames.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
