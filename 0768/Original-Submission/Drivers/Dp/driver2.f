      
      program example2
     
c     TENSOLVE finds roots of systems of n nonlinear equations in
c     n unknowns, or minimizers of the sum of squares of m > n
c     nonlinear functions in n unknowns, using tensor methods.
c     
c     This example illustrates the use of TENSOLVE to solve a
c     nonlinear least-squares problem defined by subroutines
c     fwood and jwood (included below).
     
      integer           maxm, maxn, maxp, m, n, itnlim, jacflg, method,
     +                  global, ipr, lunc, lnem, lnen, lin, msg, termcd
      double precision  gradtl, steptl, ftol, stepmx, dlt
      parameter         (maxm = 100, maxn = 30, maxp = 6)
      parameter         (lin = 3, lunc = 14, lnem = 51, lnen = 19)
      integer           iwrkn(maxn,lin)
      double precision  x0(maxn), xp(maxn), fp(maxm), gp(maxn),
     +                  typx(maxn), typf(maxm), 
     +                  wrknen(maxn,lnen), wrkunc(maxp,lunc),
     +                  wrknem(maxm,lnem)
      external          fwood, jwood
      
c     Set dimensions of the problem.
      
      m      = 6
      n      = 4

c     Set values for the initial point.
      
      x0(1)  = -300.0d0
      x0(2)  = -100.0d0
      x0(3)  = -300.0d0
      x0(4)  = -100.0d0

c     Set default values for the TENSOLVE parameters.

      call tsdflt(m     , n     , itnlim, jacflg, gradtl, steptl,
     +            ftol  , method, global, stepmx, dlt   ,
     +            typx  , typf  , ipr   , msg    )

c     Alter some of the parameters.

      gradtl = 1.0d-5
      ftol   = 1.0d-9
      steptl = 1.0d-9
      global = 1

c     Call TENSOLVE.

      call tsneci(maxm  , maxn  , maxp  , x0    , m     , n     ,
     +            typx  , typf  , itnlim, jacflg, gradtl, steptl, 
     +            ftol  , method, global, stepmx, dlt   , ipr   , 
     +            wrkunc, lunc  , wrknem, lnem  , wrknen, lnen  ,
     +            iwrkn , lin   , fwood , jwood , msg   , 
     +            xp    , fp    , gp    , termcd )
      
c     end of example2 main program
      end


      subroutine fwood ( x, f, m, n )
      integer            m, n
      double precision   x(n), f(m)

c     fwood defines function values for the Wood function.
      
      f(1) = 10.0d0 * (x(2) - x(1)**2)
      f(2) = 1.0d0 - x(1)
      f(3) = sqrt(90.0d0) * (x(4) - x(3)**2)
      f(4) = 1.0d0 - x(3)
      f(5) = (x(2) + x(4) - 2.0d0) * sqrt(10.0d0)
      f(6) = (x(2) - x(4)) / sqrt(10.0d0)

c     end of fwood.
      end
      

      subroutine jwood ( x, jac, maxm, m, n )
      integer            maxm, m, n
      double precision   x(n), jac(maxm,n)

c     jwood defines Jacobian values for the Wood function.

      integer            i, j
      double precision   tval

      do 20 j = 1, n
         do 10 i = 1, m
            jac(i,j) = 0.0d0
 10      continue
 20   continue
      
      jac(1,1) = -20.0d0 * x(1)
      jac(1,2) =  10.0d0
      
      jac(2,1) = -1.0d0
      
      tval     =  sqrt(90.0d0)
      jac(3,3) = -2.0d0 * tval * x(3)
      jac(3,4) =  tval
      
      jac(4,3) = -1.0d0
      
      tval     =  sqrt(10.0d0)
      jac(5,2) =  tval
      jac(5,4) =  tval
      
      tval     =  1.0d0/tval
      jac(6,2) =  tval
      jac(6,4) = -tval

c     end of jwood.
      end


