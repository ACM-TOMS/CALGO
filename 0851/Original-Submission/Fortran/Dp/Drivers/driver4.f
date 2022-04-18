c After the first iteration, the line search might take too large a
c step, when it tries to bracket a local minimizer along the search
c direction. In this expansion phase, the expansion factor is rho.
c The default value of rho is 5. Thus, in each step of the expansion
c phase, the bracketing interval grows by a factor of 5. If the cost
c function contains a log or an exponential function, this
c expansion might result in either an overflow or an attempt to compute
c the log of a negative number. A less aggressive expansion may avoid
c this problem. In the following example, we choose a small initial
c stepsize (Step is true, and initial stepsize is 1.e-5), QuadStep is
c turned off, and rho is 1.5. The code has to do a number of expansions
c to reach a suitable interval bracketing the minimizer in the initial
c search direction.

c Termination status:  0
c Convergence tolerance for gradient satisfied
c absolute largest component of gradient: 0.6748D-08
c function value:   -653.07867273306
c cg iterations:  34
c function evaluations:  58
c gradient evaluations:  92

c We then set rho = 5. In this case, fewer expansions are needed to bracket
c the minimizer in the initial search direction:

c Termination status:  0
c Convergence tolerance for gradient satisfied
c absolute largest component of gradient: 0.6028D-08
c function value:   -653.07867273306
c cg iterations:  35
c function evaluations:  42
c gradient evaluations:  77

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

      gnorm = 1.e-5
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

c cg.parm (QuadStep turned off, initial stepsize is specified, expansion factor
c rho is 1.5):

c .1d0      delta        (Wolfe line search parameter)
c .9d0      sigma        (Wolfe line search parameter)
c 1.d-6     eps          (perturbation parameter for computing fpert)
c .66d0     gamma        (required decay factor in interval)
c 1.5d0     rho          (interval growth factor used to get bracketing interval)
c .01d0     eta          (lower bound for cg's beta_k)
c .01d0     psi0         (factor used in starting guess for iteration 1)
c .1d0      psi1         (factor previous step multiplied by in QuadStep)
c 2.d0      psi2         (factor previous step is multipled by for startup)
c 1.d-12    QuadCutOff   (QuadStep if relative change in f > QuadCutOff)
c 0.d-12    StopFact     (factor multiplying starting |grad|_infty in StopRule)
c 1.d-3     AWolfeFac    (AWolfe = F => set AWolfe = T if |f-f0| < Awolfe_fac*Ck)
c 1.0d0     restart_fac  (restart cg in restart_fac*n iterations)
c 500.d0    maxit_fac    (terminate in maxit_fac*n iterations)
c 0.d0      feps         (stop when value change <= feps*|f|)
c .7d0      Qdecay       (used in Qk update: Qk = Qdecay*Qk + 1)
c 50        nexpand      (number of grow/shrink allowed in bracket)
c 50        nsecant      (number of secant steps allowed in line search)
c .true.    PertRule     (F => eps, T => eps*Ck)
c .false.   QuadStep     (use initial quad interpolation in line search)
c .false.   PrintLevel   F (no print) T (intermediate results)
c .true.    PrintFinal   F (no print) T (print error messages, final error)
c .true.    StopRule     T (|grad|_infty <= max(tol,|grad|_0*StopFact) F (... <= tol*(1+|f|))
c .false.   AWolfe       F (Wolfe -- see AWolfeFac above) T (approx Wolfe)
c .true.    Step         F (no initial line search guess) T (guess in gnorm)
c .false.   debug        F (no debugging) T (check for no increase in f)
