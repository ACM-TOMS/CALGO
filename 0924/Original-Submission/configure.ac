AC_INIT([libTIDES], [1.3], [tides@unizar.es],[tides],[http://gme.unizar.es/software/tides])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR([configure.ac])

AC_PROG_RANLIB
AC_PROG_CC

AX_CC_MAXOPT

AC_ARG_ENABLE([fortran-tests],
  [  --enable-fortran-tests     enable fortran tests],
  [fortran_tests=${enableval}],
  [fortran_tests=no])

if test "x${fortran_tests}" = xyes; then
  AC_PROG_FC
  AC_FC_SRCEXT([f])
fi
AM_CONDITIONAL([HAVE_FORTRAN_TESTS],[test "x${fortran_tests}" = xyes])

my_save_cflags="$CFLAGS"
CFLAGS=-m64
AC_MSG_CHECKING([whether CC supports -m64])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_CFLAGS=-m64],
    [AC_MSG_RESULT([no])]
)
CFLAGS="$my_save_cflags"
AC_SUBST([AM_CFLAGS])

AC_ARG_ENABLE([multiple-precision],
  [  --disable-multiple-precision     disable multiple precision],
  [multiple_precision=${enableval}],
  [multiple_precision=yes])

if test "x${multiple_precision}" = xyes; then
  TIDES_WITH_GMP
  TIDES_LIB_GMP
  TIDES_WITH_MPFR
  TIDES_LIB_MPFR
fi
AM_CONDITIONAL([HAVE_MP],[test "x${multiple_precision}" = xyes])

AC_C_INLINE

if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE(HAVE_INLINE,1,"adding inline flag")
  AC_SUBST(HAVE_INLINE)
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile libTIDES/Makefile test/Makefile])
AC_OUTPUT
echo \
"-------------------------------------------------

 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'.
 C Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
 F Compiler: '${FC} ${FCFLAGS}'

 Package features:
   Multiple precision: ${multiple_precision}

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything

--------------------------------------------------"
