c Although there is a rigorous theory justifying a Wolfe line search,
c the performance of the Approximate Wolfe line search if often much
c better. Nonetheless, the user can completely turn off the Approximate
c Wolfe line search by setting AWolfe to false and AWolfeFac to 0.
c The output is the following:

c Termination status:  4
c Line search fails, too many secant steps
c Possible causes of this error message:
c   - your tolerance (grad_tol = 0.1000D-07) may be too strict
c absolute largest component of gradient: 0.2923D-06
c function value:   -653.07867273306
c cg iterations:  28
c function evaluations:  146
c gradient evaluations:  132

c Hence, due to numerical errors, it was not possible to achieve the specified
c 1.e-8 error tolerance. By decreasing the error tolerance to 1.e-6 (the first
c argument of cg_descent), the problem is solved:

c Termination status:  0
c Convergence tolerance for gradient satisfied
c absolute largest component of gradient: 0.7722D-06
c function value:   -653.07867273306
c cg iterations:  26
c function evaluations:  49
c gradient evaluations:  33

c On the other hand, if we turn on the Approximate Wolfe line search by
c resetting AWolfeFac to its default value 1.e-3, we
c obtain slightly faster convergence:

c Termination status:  0
c Convergence tolerance for gradient satisfied
c absolute largest component of gradient: 0.9577D-06
c function value:   -653.07867273306
c cg iterations:  25
c function evaluations:  48
c gradient evaluations:  31

      integer maxsize

      parameter (maxsize = 100000)

      double precision x (maxsize), d (maxsize), g (maxsize),
     &                 xtemp (maxsize), gtemp (maxsize), gnorm, f

      integer i, n, status, iter, nfunc, ngrad

      external myvalue, mygrad

      n = 100

c starting guess

      do 10 i = 1, n
          x (i) = 1.d0
10    continue

      call cg_descent (1.d-8, x, n, myvalue, mygrad, status, gnorm, f,
     &                 iter, nfunc, ngrad, d, g, xtemp, gtemp)
      
      end

      subroutine myvalue (f, x, n)
      integer i, n
      double precision x (*), f, t
      f = 0.d0
      do 10 i = 1, n
          t = i
          t = dsqrt (t)
          f = f +  dexp (x (i)) - t*x (i)
10    continue
      return
      end

      subroutine mygrad (g, x, n)
      integer i, n
      double precision g (*), x (*), t
      do 10 i = 1, n
          t = i
          t = dsqrt (t)
          g (i) = dexp (x (i)) -  t
10    continue
      return
      end

c cg.parm (approximate Wolfe line search turned off):

c .1d0      delta        (Wolfe line search parameter)
c .9d0      sigma        (Wolfe line search parameter)
c 1.d-6     eps          (perturbation parameter for computing fpert)
c .66d0     gamma        (required decay factor in interval)
c 5.0d0     rho          (interval growth factor used to get bracketing interval)
c .01d0     eta          (lower bound for cg's beta_k)
c .01d0     psi0         (factor used in starting guess for iteration 1)
c .1d0      psi1         (factor previous step multiplied by in QuadStep)
c 2.d0      psi2         (factor previous step is multipled by for startup)
c 1.d-12    QuadCutOff   (QuadStep if relative change in f > QuadCutOff)
c 0.d-12    StopFact     (factor multiplying starting |grad|_infty in StopRule)
c 0.d-3     AWolfeFac    (AWolfe = F => set AWolfe = T if |f-f0| < Awolfe_fac*Ck)
c 1.0d0     restart_fac  (restart cg in restart_fac*n iterations)
c 500.d0    maxit_fac    (terminate in maxit_fac*n iterations)
c 0.d0      feps         (stop when value change <= feps*|f|)
c .7d0      Qdecay       (used in Qk update: Qk = Qdecay*Qk + 1)
c 50        nexpand      (number of grow/shrink allowed in bracket)
c 50        nsecant      (number of secant steps allowed in line search)
c .true.    PertRule     (F => eps, T => eps*Ck)
c .true.    QuadStep     (use initial quad interpolation in line search)
c .false.   PrintLevel   F (no print) T (intermediate results)
c .true.    PrintFinal   F (no print) T (print error messages, final error)
c .true.    StopRule     T (|grad|_infty <= max(tol,|grad|_0*StopFact) F (... <= tol*(1+|f|))
c .false.   AWolfe       F (Wolfe -- see AWolfeFac above) T (approx Wolfe)
c .false.   Step         F (no initial line search guess) T (guess in gnorm)
c .false.   debug        F (no debugging) T (check for no increase in f)
