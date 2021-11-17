#! /bin/sh
# Generate test data for double- and quadruple-precision tests.
#
# Usage:
#	./gendat.sh xmin xmax function-name
#
# [31-Jul-2000]
LIBDIR=`dirname $0`
gawk -f  $LIBDIR/gendat.awk -v xmin=$1 -v xmax=$2 -v fun=$3
