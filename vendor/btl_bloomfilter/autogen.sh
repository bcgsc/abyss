#!/bin/sh
set -ex
aclocal
autoconf
autoheader
automake -a
