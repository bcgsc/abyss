#!/bin/sh

MIN_AUCONF_VERSION="2.62"
MIN_AUMAKE_VERSION="1.9.6"

AUCONF_VERSION=`autoconf --version 2>/dev/null | awk '/^autoconf/ { print $4 }'`
AUMAKE_VERSION=`automake --version 2>/dev/null | awk '/^automake/ { print $4 }'`

##
# Compare version strings.
#
# Echo '0' if the two version strings are equal, '1' if the $1 represents a
# later version than $2, or '-1' if $1 represents an earlier version of than
# $2.
#
# E.g.:
#
#       compare_version 1.2     1.2
#       # => 0
#       compare_version 1.10    1.2
#       # => 1
#       compare_version 1.2.10  1.3.0
#       # => -1
#
function compare_version() {
    echo $1 $2 | awk '
    { split($1, a, ".");
      split($2, b, ".");
      x=0;
      for (i = 1; i <= 4; i++)
          if (a[i] < b[i]) {
              x =-1;
              break;
          } else if (a[i] > b[i]) {
              x = 1;
              break;
          }
      print x;
   }'
}

if [ -z "${AUCONF_VERSION}" ]; then
    echo "autoconf not found" 1>&2
elif [ -z "${AUMAKE_VERSION}" ]; then
    echo "automake not found" 1>&2
elif [ `compare_version $AUCONF_VERSION $MIN_AUCONF_VERSION` = "-1" ]; then
    echo "autoconf version $MIN_AUCONF_VERSION or later required" 2>&1
elif [ `compare_version $AUMAKE_VERSION $MIN_AUMAKE_VERSION` = "-1" ]; then
    echo "automake version $MIN_AUMAKE_VERSION or later required" 2>&1
else
    autotools_ok=1
fi

if [ "${autotools_ok}" != 1 ]; then
    echo "
        Building ABySS from source requires a recent version of autotools.

        To compile, either install the latest versions of automake and
        autoconf and re-run autogen.sh, or download the pre-configured source
        from https://github.com/bcgsc/abyss/releases

    " 2>&1
    exit 1
fi

set -ex
aclocal
autoconf
autoheader
automake -a
