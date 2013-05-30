#!/bin/bash
# continuous integration test script
# run this from the top-level source directory

# Run clang's scan-build static analysis tool

cd unix
rm -Rf objs
make clean
scan-build -o ./scan-build make
cd scan-build
echo == Updating $JENKINS_URL/job/$JOB_NAME/ws/scan-build/current/
rm -f current && ln -s `ls -td ????-* | head -1` current
cd ..
# Validate the generated code as proxy for testing under clang.
make check
