#/c!/bin/bash
#
# Shell script that checks both the basic dispmodule and the non-default-kinds add-on
# modules using the compilers g95, gfortran, f95 from Nag, ifort, af95 from Absoft and
# pathf95 from Pathscale (or those of these that are installed). The shell script can
# be run from Cygwin on MS Windows if the compilers are set up to work with Cygwin.
# There are two checks which are made:
#
#     A. The modules are compiled with swiches that provide evidence of conformance
#        with the Fortran 95 standard and/or switches that warn when Fortran 95 extensions
#        are used.
#     B. Test programs are run via makefile targets ("check", "check-all-kinds",
#        "check-naninf", and, for Nag f95, "check-naninf-ieee". For the "check" target
#        all available checking compile flags are applied (or at least it has been
#        attempted to include them all). Note that some new versions of g95 are
#        ieee-compliant and check-naninf-ieee works with them.
#
# To use the script one can issue the command: ./test-all.sh > result.log 2>&1

platform=windows
if test `uname` = SunOS
    then make=gmake; platform=unix
    else make=make
fi
if test `uname` = Linux
    then platform=unix
fi

echo "OPERATING SYSTEM:"
echo "*****************"
echo   "    uname -a:        "`uname -a`
if test -e /etc/issue; then 
  echo "    cat /etc/issue:  "`cat /etc/issue | grep -v '^$' | head -1`; 
fi
echo

# Standard conformance:
echo "                                     EVIDENCE OF STANDARD CONFORMANCE"
echo "                                     ================================"
echo
set -x
if test -x "`which g95 2>/dev/null`"; then
    g95 -c -O0 -std=f95 -pedantic dispmodule.f90 disp_i?mod.f90 disp_l?mod.f90
fi
if test -x "`which gfortran 2>/dev/null`"; then
    gfortran -c -O0 -std=f95 -pedantic dispmodule.f90 disp_i?mod.f90 disp_l?mod.f90
fi
if test -x "`which ifort 2>/dev/null`"; then
    if test $platform = windows
	then ifort -c -nologo -Od -stand:f95 dispmodule.f90 disp_*mod.f90
	else ifort -c -nologo -O0 -stand:f95 dispmodule.f90 disp_*mod.f90
    fi
fi
set +x

# Check:
#   main module with full checks,
#   non-default kind modules with default checks
#   display of NaN-s and infinities

echo
echo "                                          EVIDENCE OF PORTABILITY"
echo "                                          ======================="

if test -x "`which g95 2>/dev/null`"; then
    echo -e \\nCOMPILER: `g95 --version | head -1`
    echo -----------------------------------------------------------------
    set -x
    $make clean
    $make check compiler=g95 fflags="-c -Wall -Wextra -fbounds-check -fimplicit-none"
    $make check-all-kinds compiler=g95
    $make check-naninf compiler=g95
    set +x
fi

if test -x "`which gfortran 2>/dev/null`"; then
    echo -e \\nCOMPILER: `gfortran --version | head -1`
    echo -----------------------------------------------------------------
    set -x
    $make clean
    $make check compiler=gfortran fflags="-c -Wall -W -fbounds-check -fimplicit-none"
    $make check-all-kinds compiler=gfortran
    $make check-naninf compiler=gfortran
    set +x
fi

if test -x "`which af95 2>/dev/null`"; then
    echo -e \\nCOMPILER: `af95 -V | head -2`
    echo -----------------------------------------------------------------
    set -x
    $make clean
    $make check compiler=absoft
    $make check-all-absoft-kinds compiler=absoft
    $make check-naninf compiler=absoft
    set +x
fi

if test -x "`which pathf95 2>/dev/null`"; then
    echo -e \\nCOMPILER: `pathf95 -v | head -1`
    echo -----------------------------------------------------------------
    set -x
    $make clean
    $make check compiler=pathf95
    $make check-all-kinds compiler=pathf95
    $make check-naninf compiler=pathf95
    $make check-naninf-ieee compiler=pathf95
    set +x
fi

if test -x "`which f95 2>/dev/null`"
    then if f95 -V 2>&1 | grep NAG > /dev/null; then
	echo -e \\nCOMPILER: `f95 -V 2>&1 | head -1`
	echo -----------------------------------------------------------------
	set -x
	$make clean
	$make check compiler=nag fflags="-c -C=all -C=undefined"
        # (must get rid of -C=undefined compiled module):
	$make clean
	$make check-all-kinds compiler=nag
	$make check-naninf compiler=nag
	$make check-naninf-ieee compiler=nag
	set +x
    fi
fi

if test -x "`which ifort 2>/dev/null`"; then
    if test $platform = windows
        then echo -e \\nCOMPILER: `ifort /logo 2>&1 | head -1`
        else echo -e \\nCOMPILER: `ifort --version 2>&1 | head -1`
    fi
    echo -----------------------------------------------------------------
    set -x
    $make clean platform=$platform compiler=ifort
    if test $platform = windows
        then fflags="/c /Od /nologo /warn:all /check:bounds,format,pointers,uninit"
	else fflags="-c -O0 -nologo -warn all -check bounds,format,pointers,uninit"
    fi
    $make check platform=$platform compiler=ifort fflags="$fflags"
    $make check-all-kinds platform=$platform compiler=ifort
    $make check-quadprec platform=$platform compiler=ifort
    $make check-naninf platform=$platform compiler=ifort
    set +x
fi
echo -e \\n\\n\\n