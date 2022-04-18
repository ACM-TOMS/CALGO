c The operation of the code is mostly controlled by the parameter file
c cg.parm. In the following example, the parameter QuadStep is set to
c false. When QuadStep is true, the trial step in each iteration is
c computed as the minimizer of a quadratic interpolant along the search
c direction. In performing the quadstep, we hope to find a suitable
c line search point right away, completely by-passing the secant iteration.
c However, as the iterates approach a minimizer, the numerical accuracy of
c the minimizer of the quadratic interpolant becomes worse. When the relative
c change in the function values for two consecutive iterations reach
c QuadCutOff, then the code completely turns off the quadstep. The user
c can turn off the quadstep by making QuadStep false in cg.parm. By leaving
c QuadStep true, but increasing QuadCutOff (default 1.d-12), the code turns
c off the QuadStep sooner. Note that the code is usually faster when the
c quadstep is employed. For the test problem, turning off QuadStep
c (setting its value to false in cg.parm) yields the following output.
c Notice that the number of gradient evaluations increases while the
c number of function evaluations decreases (compared to the previous run).
c This behavior is typical.

c Termination status:  0
c Convergence tolerance for gradient satisfied
c absolute largest component of gradient: 0.6273D-08
c function value:   -653.07867273306
c cg iterations:  35
c function evaluations:  39
c gradient evaluations:  74

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

c cg.parm (QuadStep false):

c  .1d0      delta        (Wolfe line search parameter)
c  .9d0      sigma        (Wolfe line search parameter)
c  1.d-6     eps          (perturbation parameter for computing fpert)
c  .66d0     gamma        (required decay factor in interval)
c  5.0d0     rho          (interval growth factor used to get bracketing interval)
c  .01d0     eta          (lower bound for cg's beta_k)
c  .01d0     psi0         (factor used in starting guess for iteration 1)
c  .1d0      psi1         (factor previous step multiplied by in QuadStep)
c  2.d0      psi2         (factor previous step is multipled by for startup)
c  1.d-12    QuadCutOff   (QuadStep if relative change in f > QuadCutOff)
c  0.d-12    StopFact     (factor multiplying starting |grad|_infty in StopRule)
c  1.d-3     AWolfeFac    (AWolfe = F => set AWolfe = T if |f-f0| < Awolfe_fac*Ck)
c  1.0d0     restart_fac  (restart cg in restart_fac*n iterations)
c  500.d0    maxit_fac    (terminate in maxit_fac*n iterations)
c  0.d0      feps         (stop when value change <= feps*|f|)
c  .7d0      Qdecay       (used in Qk update: Qk = Qdecay*Qk + 1)
c  50        nexpand      (number of grow/shrink allowed in bracket)
c  50        nsecant      (number of secant steps allowed in line search)
c  .true.    PertRule     (F => eps, T => eps*Ck)
c  .false.   QuadStep     (use initial quad interpolation in line search)
c  .false.   PrintLevel   F (no print) T (intermediate results)
c  .true.    PrintFinal   F (no print) T (print error messages, final error)
c  .true.    StopRule     T (|grad|_infty <= max(tol,|grad|_0*StopFact) F (... <= tol*(1+|f|))
c  .false.   AWolfe       F (Wolfe -- see AWolfeFac above) T (approx Wolfe)
c  .false.   Step         F (no initial line search guess) T (guess in gnorm)
c  .false.   debug        F (no debugging) T (check for no increase in f)
