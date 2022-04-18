#! /bin/sh
# Generate test data for single-precision tests.
#
# Usage:
#	./gendats.sh xmin xmax function-name
#
# [31-Jul-2000]
LIBDIR=`dirname $0`
gawk -f  $LIBDIR/gendats.awk -v xmin=$1 -v xmax=$2 -v fun=$3
