#! /bin/sh
########################################################################
### Convert one or more Fortran files to man pages in the subdirectory
### tmpman.  PostScript versions are generated too.
###
### Usage:
###	fortoman.sh Fortran-file(s)
###
### [04-Nov-2003] -- add check for multiple directories
### [21-Dec-2000] -- add SEEALSO support
### [18-Dec-2000]
########################################################################

AWK=${AWK:-mawk}
MAN2PS=${MAN2PS:-man2ps}
test -d tmpman || mkdir tmpman
for f in "$@"
do
	g=`basename $f .f`
	echo $g
	srcdir=common
	for d in common gampsi jacobi laguerre
	do
		case $d in
		gampsi) AUTHORFILE=AUTHOR.gam ;;
		*)	AUTHORFILE=AUTHOR.std ;;
		esac
		if test -f $d/$f
		then
			$AWK -f fortoman.awk \
				-v AUTHORFILE=$AUTHORFILE \
				-v INCLUDEFILES="gampsi.h gjl.h" \
				-v SEEALSO=SEEALSO.txt \
				$d/$f >tmpman/$g.man
			$MAN2PS < tmpman/$g.man > tmpman/$g.ps 2>/dev/null
		fi
	done
done
