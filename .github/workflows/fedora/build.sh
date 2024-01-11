#!/bin/sh

set -e

dnf install -y gcc make autoconf automake libtool
dnf install -y atlas-devel fftw-devel plplot-devel curl-devel

# test dependencies
dnf install -y python3-numpy.x86_64 python3-astropy.x86_64

USER=root
export USER

cd "$WORK_DIR"

./autogen.sh
./configure
make

# run the tests
make check
