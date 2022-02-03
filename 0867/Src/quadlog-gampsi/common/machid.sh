#! /bin/sh
#=======================================================================
# Print on stdout a string of the form VENDOR-HARDWARE-OS-LEVEL,
# where each of the four parts are guaranteed to not contain hyphens,
# and all characters are suitable for use in a UNIX filename within a
# single directory.
#
# Usage:
#	./machid.sh
#
# [01-Jun-2000]
#=======================================================================

simplify()
{
    echo $1 | sed -e 's/[^A-Za-z0-9.]/_/g'
}

os=`uname -s || echo UNKNOWN`
os=`simplify "$os"`

case $os in
	AIX)			# IBM RS/6000
		vendor=IBM
		hardware=RS6000
		level=`oslevel`
		;;

	HP_UX)			# Hewlett-Packard 9000
		vendor=HP
		hardware=`uname -m`
		level=`uname -r`
		;;

	IRIX64)			# Silicon Graphics
		vendor=SGI
		hardware=`uname -m`
		level=`uname -r`
		;;

	Linux)			# GNU/Linux
		hardware=`uname -m`
		level=`uname -r`
		case $hardware in
			alpha)
				vendor=Compaq_DEC
				;;
			i*)
				vendor=Intel
				;;
			ppc)
				vendor=Apple
				;;
			sparc)
				vendor=Sun
				;;
			*)
				vendor=UNKNOWN
				;;
		esac
		;;

	Mach)			# NeXT
		vendor=NeXT
		hardware=`uname -m`
		level=`uname -r`
		;;

	OSF1)			# Compaq/DEC OSF/1 == Tru64
		vendor=Compaq_DEC
		hardware=`uname -m`
		level=`uname -r`
		;;

	Rhapsody)		# Apple MacOS 10
		vendor=Apple
		hardware=`uname -m`
		level=`uname -r`
		;;

	SunOS)			# Sun SunOS or Sun Solaris
		vendor=Sun
		hardware=`uname -m`
		level=`uname -r`
		;;

	UNKNOWN|*)
		vendor=UNKNOWN
		hardware=UNKNOWN
		level=UNKNOWN
		;;
esac

vendor=`simplify "$vendor"`
hardware=`simplify "$hardware"`
level=`simplify "$level"`

echo $vendor-$hardware-$os-$level
