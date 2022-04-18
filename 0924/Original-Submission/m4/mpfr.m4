AC_DEFUN([TIDES_LIB_MPFR], [
    AC_MSG_CHECKING([whether mpfr library is installed])
    AH_TEMPLATE([HAVE_MPFR],[Define if you have the mpfr library])
    LIBS="$LIBS -lmpfr -lgmp -lm"
AC_TRY_LINK(
    [
#include <gmp.h>
#include <mpfr.h>
    ],
    [mpfr_t x; mpfr_init (x);],
    [CONFIG_MPFR="#define MML_SUPPORTS_MPFR"],
    AC_MSG_RESULT([no]);AC_MSG_ERROR([missing library: mpfr]))
    AC_MSG_RESULT([yes])
])

AC_DEFUN([TIDES_WITH_MPFR], [
    AC_ARG_WITH(mpfr,
	[  --with-mpfr=DIR         Mpfr installation directory],
	with_mpfr=$withval)

    if test "$with_mpfr" != ""; then
        CPPFLAGS="-I$with_mpfr/include $CPPFLAGS"
        LIBS="-L$with_mpfr/lib $LIBS"
    fi
])
