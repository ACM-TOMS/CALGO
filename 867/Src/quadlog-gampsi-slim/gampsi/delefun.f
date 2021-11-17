      PROGRAM delefun
************************************************************************
*     (Test of special arguments to elementary functions)
*     Test Fortran elementary functions with special arguments that
*     generate exceptional return values.
*     [10-Jun-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           dacos,       dasin,       datan,       dcos
      INTRINSIC           dcosh,       dexp,        dlog,        dlog10
      INTRINSIC           dsin,        dsinh,       dsqrt,       dtan
      INTRINSIC           dtanh
*
*     Built-in functions
*
      DOUBLE PRECISION    dacos,       dasin,       datan,       dcos
      DOUBLE PRECISION    dcosh,       dexp,        dlog,        dlog10
      DOUBLE PRECISION    dsin,        dsinh,       dsqrt,       dtan
      DOUBLE PRECISION    dtanh
*
*     External functions
*
      EXTERNAL            f
*
      DOUBLE PRECISION    f
*
*     Parameter variables
*
      DOUBLE PRECISION    ZERO
      PARAMETER           (ZERO = 0.0d+00)
      DOUBLE PRECISION    ONE
      PARAMETER           (ONE = 1.0d+00)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      DOUBLE PRECISION    inf,         nan
*
      inf = ZERO
      inf = ONE/inf
*
      nan = inf - inf
*
      WRITE (stdout, 10000) 'dacos( ',  inf,   ' ) = ', f(dacos,inf)
      WRITE (stdout, 10000) 'dasin( ',  inf,   ' ) = ', f(dasin,inf)
      WRITE (stdout, 10000) 'datan( ',  inf,   ' ) = ', f(datan,inf)
      WRITE (stdout, 10000) 'dcos( ',   inf,   ' ) = ', f(dcos,inf)
      WRITE (stdout, 10000) 'dcosh( ',  inf,   ' ) = ', f(dcosh,inf)
      WRITE (stdout, 10000) 'dexp( ',   inf,   ' ) = ', f(dexp,inf)
      WRITE (stdout, 10000) 'dlog( ',   inf,   ' ) = ', f(dlog,inf)
      WRITE (stdout, 10000) 'dlog10( ', inf,   ' ) = ', f(dlog10,inf)
      WRITE (stdout, 10000) 'dsin( ',   inf,   ' ) = ', f(dsin,inf)
      WRITE (stdout, 10000) 'dsinh( ',  inf,   ' ) = ', f(dsinh,inf)
      WRITE (stdout, 10000) 'dsqrt( ',  inf,   ' ) = ', f(dsqrt,inf)
      WRITE (stdout, 10000) 'dtan( ',   inf,   ' ) = ', f(dtan,inf)
      WRITE (stdout, 10000) 'dtanh( ',  inf,   ' ) = ', f(dtanh,inf)
*
      WRITE (stdout, 10000) 'dacos( ',  nan,   ' ) = ', f(dacos,nan)
      WRITE (stdout, 10000) 'dasin( ',  nan,   ' ) = ', f(dasin,nan)
      WRITE (stdout, 10000) 'datan( ',  nan,   ' ) = ', f(datan,nan)
      WRITE (stdout, 10000) 'dcos( ',   nan,   ' ) = ', f(dcos,nan)
      WRITE (stdout, 10000) 'dcosh( ',  nan,   ' ) = ', f(dcosh,nan)
      WRITE (stdout, 10000) 'dexp( ',   nan,   ' ) = ', f(dexp,nan)
      WRITE (stdout, 10000) 'dlog( ',   nan,   ' ) = ', f(dlog,nan)
      WRITE (stdout, 10000) 'dlog10( ', nan,   ' ) = ', f(dlog10,nan)
      WRITE (stdout, 10000) 'dsin( ',   nan,   ' ) = ', f(dsin,nan)
      WRITE (stdout, 10000) 'dsinh( ',  nan,   ' ) = ', f(dsinh,nan)
      WRITE (stdout, 10000) 'dsqrt( ',  nan,   ' ) = ', f(dsqrt,nan)
      WRITE (stdout, 10000) 'dtan( ',   nan,   ' ) = ', f(dtan,nan)
      WRITE (stdout, 10000) 'dtanh( ',  nan,   ' ) = ', f(dtanh,nan)

*
10000 FORMAT (A, 1X, 1P, E21.08E4, 1X, A, E41.15E4)
*
      END


      DOUBLE PRECISION FUNCTION f(fun,x)
************************************************************************
*     (Function evaluation)
*     Evaluate and return fun(x). This separate function is necessary to
*     foil clever optimizers (e.g., SGI f90) which diagnose special
*     arguments to intrinsics at compile time as fatal errors.
*     [10-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            fun
*
      DOUBLE PRECISION             fun
*
*     Argument variables
*
      DOUBLE PRECISION             x
*
      f = fun(x)
*
      END
