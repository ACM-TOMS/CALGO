# Process this file with autoconf to produce a configure script.

# Set various version strings
m4_define([genial_major_version],    [2])
m4_define([genial_minor_version],    [1])
m4_define([genial_subminor_version], [0])
m4_define([genial_version],          [genial_major_version.genial_minor_version.genial_subminor_version])

# Before making a release, the LT_VERSION string should be modified.
# The string is of the form C:R:A.
# - If interfaces have been changed or added, but binary compatibility has
#   been preserved, change to C+1:0:A+1
# - If binary compatibility has been broken (eg removed or changed interfaces)
#   change to C+1:0:0
# - If the interface is the same as the previous version, change to C:R+1:A
m4_define([genial_lt_version],       [0:0:0])


# init autotools
AC_PREREQ(2.59)
AC_INIT([Genial],[genial_version],[laurent@ient.rwth-aachen.de])
AC_CONFIG_SRCDIR([genial/utility.h])
AC_CONFIG_AUX_DIR([autotools])
AC_CONFIG_MACRO_DIR([autotools/aclocal])

# put CVS revision of this configure.in into the resulting configure script
AC_REVISION($Revision: 1.9 $)

# Detect the canonical host and target build environment
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

# setup automake - we need >= 1.6 because of aclocal support
AM_CONFIG_HEADER([genial/autoconfig.h])
AM_INIT_AUTOMAKE([1.6 foreign dist-zip])

# substitute versionnumbers
AC_SUBST(MMAJOR,     [genial_major_version])
AC_SUBST(MMINOR,     [genial_minor_version])
AC_SUBST(MSUBMINOR,  [genial_subminor_version])
AC_SUBST(VERSION,    [genial_version])
AC_SUBST(LT_VERSION, [genial_lt_version])

# checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

# check if the user wants to build with stlport
# the configure.ac from OpenOffice has a similar option, 
# but it works quite different
#AC_ARG_WITH([stlport],
#            [AS_HELP_STRING([--with-stlport],
#                            [specify if we are using STLport])], 
#            [with_stlport=$withval], 
#            [with_stlport="yes"])
#AC_MSG_CHECKING([if we are using stlport])
#AC_PREPROC_IFELSE(
#   [AC_LANG_PROGRAM([[#include<vector>]],[[
#ifndef _STLPORT_VERSION
#  #error _STLPORT_VERSION is not defined
#  // For some compilers (eg. SGI's CC), #error is not enough ...
#  please, do fail!
#endif
#   ]])],
#   [
#      ac_have_stlport="yes"
#      AC_MSG_RESULT([yes])
#   ],
#   [
#      ac_have_stlport="no"
#      AC_MSG_RESULT([no])
#   ])
#if test x"$with_stlport" = x"no" ; then
#   if test x"$ac_have_stlport" = x"yes" ; then
#      AC_MSG_WARN([[configuring --with-stlport=no but stlport is in the include path]])
#      with_stlport="stlport"
#   fi
#else
#   if test x"$ac_have_stlport" = x"no" ; then
#      AC_MSG_ERROR([[Could not find or compile header <vector> from stlport.
#
#There are two options now:
#
#a) ./configure CPPFLAGS=-I/absolutepathtostlportheaders --with-stlport
#b) ./configure --without-stlport
#
#The latter requires that you have a decent vendor-supplied STL implementation.
#
#]])
#   fi
#   if test x"$with_stlport" = x"yes" ; then
#      with_stlport="stlport"
#   fi
#fi
#if test x"$ac_have_stlport" = x"yes" ; then
#   AC_CHECK_LIB([$with_stlport],
#                [__get_temporary_buffer],
#                [],
#                [AC_MSG_ERROR([[Could not link to library <$with_stlport>.
#
#This might mean, your STLport shared library has a different name.
#If your STLport library is called "libstlport4.so", for example,
#you can supply that name with
#
#   ./configure --with-stlport=stlport4
#
#Also note, that the given library should reside somewhere in the
#library search path for your compiler. You can extend the search path
#via the LDFLAGS parameter
#
#   ./configure --with-stlport=stlport4 LDFLAGS="-L/usr/exotic/path/lib"
#
#]])])
#fi


# cure annoying default arguments for CFLAGS and CXXFLAGS
if test x"$CXXFLAGS" = x"-g" ; then
    CXXFLAGS=""
fi
if test x"$CXXFLAGS" = x"-g -O2" ; then
    CXXFLAGS=""
fi


# See if the user wants aggressive optimizations of the code
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],[[build debug version <default=no>]])],
              [debug=$enableval], 
              [debug="no"])
if test x"$debug" != x"no"; then
    if test x"$debug" = x"yes"; then
        if test x"$GXX" = x"yes"; then
            CXXFLAGS="-ggdb -O0 $CXXFLAGS"
        else
            CXXFLAGS="-g $CXXFLAGS"
        fi
    else
        CXXFLAGS="$debug $CXXFLAGS"
    fi
    CXXFLAGS="-DDEBUG -D_DEBUG $CXXFLAGS"
else
    case $target in
    i686-*-*)
        if test x"$GXX" = x"yes" ; then
            CXXFLAGS="-g -march=i686 -ffast-math -fomit-frame-pointer $CXXFLAGS"
        fi
        ;;
    *-*-*)
        if test x"GXX" = x"yes" ; then
            CXXFLAGS="-fomit-frame-pointer $CXXFLAGS"
        fi
        ;;
    esac
    CXXFLAGS="-O3 -DNDEBUG $CXXFLAGS"
fi



AC_CONFIG_FILES([
  Makefile
  genial/Makefile
  src/Makefile
  samples/Makefile
])
AC_OUTPUT
