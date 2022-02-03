      PROGRAM elefun
************************************************************************
*     (Test of special arguments to elementary functions)
*     Test Fortran elementary functions with special arguments that
*     generate exceptional return values.
*     [10-Jun-2000]
************************************************************************
*
*     Intrinsic functions
*
      INTRINSIC           acos,        alog,        alog10,      asin
      INTRINSIC           atan,        cos,         cosh,        exp
      INTRINSIC           sin,         sinh,        sqrt,        tan
      INTRINSIC           tanh
*
*     Built-in functions
*
      REAL                acos,        alog,        alog10,      asin
      REAL                atan,        cos,         cosh,        exp
      REAL                sin,         sinh,        sqrt,        tan
      REAL                tanh
*
*     External functions
*
      EXTERNAL            f
*
      REAL                f
*
*     Parameter variables
*
      REAL                zero
      PARAMETER           (zero = 0.0)
      REAL                one
      PARAMETER           (one = 1.0)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      REAL                inf,         nan
*
      inf = zero
      inf = one/inf
*
      nan = inf - inf
*
      WRITE (stdout, 10000) 'acos( ',  inf,   ' ) = ', f(acos,inf)
      WRITE (stdout, 10000) 'asin( ',  inf,   ' ) = ', f(asin,inf)
      WRITE (stdout, 10000) 'atan( ',  inf,   ' ) = ', f(atan,inf)
      WRITE (stdout, 10000) 'cos(  ',  inf,   ' ) = ', f(cos,inf)
      WRITE (stdout, 10000) 'cosh( ',  inf,   ' ) = ', f(cosh,inf)
      WRITE (stdout, 10000) 'exp(  ',  inf,   ' ) = ', f(exp,inf)
      WRITE (stdout, 10000) 'alog( ',  inf,   ' ) = ', f(alog,inf)
      WRITE (stdout, 10000) 'alog10( ',inf,   ' ) = ', f(alog10,inf)
      WRITE (stdout, 10000) 'sin( ',   inf,   ' ) = ', f(sin,inf)
      WRITE (stdout, 10000) 'sinh( ',  inf,   ' ) = ', f(sinh,inf)
      WRITE (stdout, 10000) 'sqrt( ',  inf,   ' ) = ', f(sqrt,inf)
      WRITE (stdout, 10000) 'tan( ',   inf,   ' ) = ', f(tan,inf)
      WRITE (stdout, 10000) 'tanh( ',  inf,   ' ) = ', f(tanh,inf)
*
      WRITE (stdout, 10000) 'acos( ',  nan,   ' ) = ', f(acos,nan)
      WRITE (stdout, 10000) 'asin( ',  nan,   ' ) = ', f(asin,nan)
      WRITE (stdout, 10000) 'atan( ',  nan,   ' ) = ', f(atan,nan)
      WRITE (stdout, 10000) 'cos( ',   nan,   ' ) = ', f(cos,nan)
      WRITE (stdout, 10000) 'cosh( ',  nan,   ' ) = ', f(cosh,nan)
      WRITE (stdout, 10000) 'exp( ',   nan,   ' ) = ', f(exp,nan)
      WRITE (stdout, 10000) 'alog( ',  nan,   ' ) = ', f(alog,nan)
      WRITE (stdout, 10000) 'alog10( ', nan,   ' ) = ', f(alog10,nan)
      WRITE (stdout, 10000) 'sin( ',   nan,   ' ) = ', f(sin,nan)
      WRITE (stdout, 10000) 'sinh( ',  nan,   ' ) = ', f(sinh,nan)
      WRITE (stdout, 10000) 'sqrt( ',  nan,   ' ) = ', f(sqrt,nan)
      WRITE (stdout, 10000) 'tan( ',   nan,   ' ) = ', f(tan,nan)
      WRITE (stdout, 10000) 'tanh( ',  nan,   ' ) = ', f(tanh,nan)

*
10000 FORMAT (A, 1X, 1P, E21.08E4, 1X, A, E41.15E4)
*
      END


      REAL FUNCTION f(fun,x)
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
      REAL                fun
*
*     Argument variables
*
      REAL             x
*
      f = fun(x)
*
      END
