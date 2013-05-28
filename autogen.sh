#!/bin/sh
# Run this to set up the build system: configure, makefiles, etc.
set -e

srcdir=`dirname $0`
test -n "$srcdir" && cd "$srcdir"

echo "Updating build configuration files for daala, please wait...."

ACLOCAL_FLAGS="-I m4"
autoreconf -if

if test "$NO_CONFIGURE" = "" ; then
  ./configure $@
fi

