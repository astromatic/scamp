#!/bin/sh

set -e

dnf install -y gcc make autoconf automake libtool
dnf install -y atlas-devel fftw-devel plplot-devel curl-devel

USER=root
export USER

cd "$WORK_DIR"

./autogen.sh
./configure
make
