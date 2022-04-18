      PROGRAM qelefun
************************************************************************
*     (Test of special arguments to elementary functions)
*     Test Fortran elementary functions with special arguments that
*     generate exceptional return values.
*
*     Because of bugs in the Cray, SGI, and Sun Fortran 90 and 95
*     compilers that prevent certain quadruple-precision intrinsics
*     from being passed as arguments, we have to provide wrappers for
*     all of them, sigh...
*
*     [10-Jun-2000]
************************************************************************
*
*     External functions
*
      EXTERNAL            f,           qqacos,      qqasin,      qqatan
      EXTERNAL            qqcos,       qqcosh,      qqexp,       qqlog
      EXTERNAL            qqlog10,     qqsin,       qqsinh,      qqsqrt
      EXTERNAL            qqtan,       qqtanh
*
      REAL*16             f,           qqacos,      qqasin,      qqatan
      REAL*16             qqcos,       qqcosh,      qqexp,       qqlog
      REAL*16             qqlog10,     qqsin,       qqsinh,      qqsqrt
      REAL*16             qqtan,       qqtanh
*
*     Parameter variables
*
      REAL*16             ZERO
      PARAMETER           (ZERO = 0.0q+00)
      REAL*16             ONE
      PARAMETER           (ONE = 1.0q+00)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      REAL*16             inf,         nan
*
      inf = ZERO
      inf = ONE/inf
*
      nan = inf - inf
*
      WRITE (stdout, 10000) 'qacos( ',  inf,   ' ) = ', f(qqacos,inf)
      WRITE (stdout, 10000) 'qasin( ',  inf,   ' ) = ', f(qqasin,inf)
      WRITE (stdout, 10000) 'qatan( ',  inf,   ' ) = ', f(qqatan,inf)
      WRITE (stdout, 10000) 'qcos( ',   inf,   ' ) = ', f(qqcos,inf)
      WRITE (stdout, 10000) 'qcosh( ',  inf,   ' ) = ', f(qqcosh,inf)
      WRITE (stdout, 10000) 'qexp( ',   inf,   ' ) = ', f(qqexp,inf)
      WRITE (stdout, 10000) 'qlog( ',   inf,   ' ) = ', f(qqlog,inf)
      WRITE (stdout, 10000) 'qlog10( ', inf,   ' ) = ', f(qqlog10,inf)
      WRITE (stdout, 10000) 'qsin( ',   inf,   ' ) = ', f(qqsin,inf)
      WRITE (stdout, 10000) 'qsinh( ',  inf,   ' ) = ', f(qqsinh,inf)
      WRITE (stdout, 10000) 'qsqrt( ',  inf,   ' ) = ', f(qqsqrt,inf)
      WRITE (stdout, 10000) 'qtan( ',   inf,   ' ) = ', f(qqtan,inf)
      WRITE (stdout, 10000) 'qtanh( ',  inf,   ' ) = ', f(qqtanh,inf)
*
      WRITE (stdout, 10000) 'qacos( ',  nan,   ' ) = ', f(qqacos,nan)
      WRITE (stdout, 10000) 'qasin( ',  nan,   ' ) = ', f(qqasin,nan)
      WRITE (stdout, 10000) 'qatan( ',  nan,   ' ) = ', f(qqatan,nan)
      WRITE (stdout, 10000) 'qcos( ',   nan,   ' ) = ', f(qqcos,nan)
      WRITE (stdout, 10000) 'qcosh( ',  nan,   ' ) = ', f(qqcosh,nan)
      WRITE (stdout, 10000) 'qexp( ',   nan,   ' ) = ', f(qqexp,nan)
      WRITE (stdout, 10000) 'qlog( ',   nan,   ' ) = ', f(qqlog,nan)
      WRITE (stdout, 10000) 'qlog10( ', nan,   ' ) = ', f(qqlog10,nan)
      WRITE (stdout, 10000) 'qsin( ',   nan,   ' ) = ', f(qqsin,nan)
      WRITE (stdout, 10000) 'qsinh( ',  nan,   ' ) = ', f(qqsinh,nan)
      WRITE (stdout, 10000) 'qsqrt( ',  nan,   ' ) = ', f(qqsqrt,nan)
      WRITE (stdout, 10000) 'qtan( ',   nan,   ' ) = ', f(qqtan,nan)
      WRITE (stdout, 10000) 'qtanh( ',  nan,   ' ) = ', f(qqtanh,nan)
*
10000 FORMAT (A, 1X, 1P, E21.08E4, 1X, A, E41.15E4)
*
      END


      REAL*16 FUNCTION f(fun,x)
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
      REAL*16             fun
*
*     Argument variables
*
      REAL*16             x
*
      f = fun(x)
*
      END


      REAL*16 FUNCTION qqacos(x)
*
*     External functions
*
      REAL*16             qacos
*
*     Argument variables
*
      REAL*16             x
*
      qqacos = qacos(x)
*
      END


      REAL*16 FUNCTION qqasin(x)
*
*     External functions
*
      REAL*16             qasin
*
*     Argument variables
*
      REAL*16             x
*
      qqasin = qasin(x)
*
      END


      REAL*16 FUNCTION qqatan(x)
*
      REAL*16             qatan
*
*     Argument variables
*
      REAL*16             x
*
      qqatan = qatan(x)
*
      END


      REAL*16 FUNCTION qqcos(x)
*
      REAL*16             qcos
*
*     Argument variables
*
      REAL*16             x
*
      qqcos = qcos(x)
*
      END


      REAL*16 FUNCTION qqcosh(x)
*
      REAL*16             qcosh
*
*     Argument variables
*
      REAL*16             x
*
      qqcosh = qcosh(x)
*
      END


      REAL*16 FUNCTION qqexp(x)
*
      REAL*16             qexp
*
*     Argument variables
*
      REAL*16             x
*
      qqexp = qexp(x)
*
      END


      REAL*16 FUNCTION qqlog(x)
*
      REAL*16             qlog
*
*     Argument variables
*
      REAL*16             x
*
      qqlog = qlog(x)
*
      END


      REAL*16 FUNCTION qqlog10(x)
*
      REAL*16             qlog10
*
*     Argument variables
*
      REAL*16             x
*
      qqlog10 = qlog10(x)
*
      END


      REAL*16 FUNCTION qqsin(x)
*
      REAL*16             qsin
*
*     Argument variables
*
      REAL*16             x
*
      qqsin = qsin(x)
*
      END


      REAL*16 FUNCTION qqsinh(x)
*
      REAL*16             qsinh
*
*     Argument variables
*
      REAL*16             x
*
      qqsinh = qsinh(x)
*
      END


      REAL*16 FUNCTION qqsqrt(x)
*
      REAL*16             qsqrt
*
*     Argument variables
*
      REAL*16             x
*
      qqsqrt = qsqrt(x)
*
      END


      REAL*16 FUNCTION qqtan(x)
*
      REAL*16             qtan
*
*     Argument variables
*
      REAL*16             x
*
      qqtan = qtan(x)
*
      END


      REAL*16 FUNCTION qqtanh(x)
*
      REAL*16             qtanh
*
*     Argument variables
*
      REAL*16             x
*
      qqtanh = qtanh(x)
*
      END
