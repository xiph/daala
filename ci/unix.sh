#!/bin/bash
# continuous integration test script
# run this from the top-level source directory

make -C unix check
make -C unix clean
