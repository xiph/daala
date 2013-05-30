#!/bin/bash
# continuous integration test script
# run this from the top-level source directory

OGG_PATH=/srv/jenkins/jobs/libogg/workspace

./autogen.sh
./configure PKG_CONFIG_PATH=${OGG_PATH}
make distcheck PKG_CONFIG_PATH=${OGG_PATH}
make docs
