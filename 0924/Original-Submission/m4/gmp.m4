AC_DEFUN([TIDES_LIB_GMP], [
    AC_MSG_CHECKING([whether gmp library is installed])
    AH_TEMPLATE([HAVE_GMP],[Define if you have the gmp library])
    LIBS="$LIBS -lgmp -lm"
    AC_TRY_LINK(
        [#include <gmp.h>],
        [mpz_t x; mpz_init(x)],
        [AC_DEFINE([HAVE_GMP],[])],
        AC_MSG_RESULT([no]);AC_MSG_ERROR([missing library: gmp]))
    AC_MSG_RESULT([yes])
])

AC_DEFUN([TIDES_WITH_GMP], [
    AC_ARG_WITH(gmp,
	[  --with-gmp=DIR          Gmp installation directory],
	with_gmp=$withval)

    if test "$with_gmp" != ""; then
        CPPFLAGS="-I$with_gmp/include $CPPFLAGS"
        LIBS="-L$with_gmp/lib $LIBS"
    fi
])
