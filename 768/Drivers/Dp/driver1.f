
      program example1
     
c     TENSOLVE finds roots of systems of n nonlinear equations in
c     n unknowns, or minimizers of the sum of squares of m > n
c     nonlinear functions in n unknowns, using tensor methods.
c     
c     This example illustrates the use of TENSOLVE to solve a
c     nonlinear equation problem defined by the subroutine
c     frosen (included below).
     
      integer           maxm, maxn, maxp, m, n, lunc, lnem, lnen, 
     +                  lin, msg, termcd
      parameter         (maxm = 100, maxn = 30, maxp = 6)
      parameter         (lin = 3, lunc = 14, lnem = 51, lnen = 19)
      integer           iwrkn(maxn,lin)
      double precision  x0(maxn), xp(maxn), fp(maxm), gp(maxn),
     +                  wrknen(maxn,lnen), wrkunc(maxp,lunc),
     +                  wrknem(maxm,lnem)
      external          frosen
      
c     Set dimensions of the problem.

      m      = 2
      n      = 2

c     Set values for the initial point.      

      x0(1)  = -1.20d0
      x0(2)  = 1.0d0

c     Call TENSOLVE.

      call tsnesi(maxm  , maxn, maxp  , x0  , m     , n   ,
     +            wrkunc, lunc, wrknem, lnem, wrknen, lnen,
     +            iwrkn , lin , frosen, msg ,
     +            xp    , fp  , gp    , termcd )
 
       
c     end of example1 main program.
      end


      subroutine frosen ( x, f, m, n )
      integer             m, n
      double precision    x(n), f(m)

c     frosen defines function values for the Rosenbrock function.
      
      f(1) = 10.0d0 * (x(2) - x(1)**2)
      f(2) = 1.0d0 - x(1)
      
c     end of frosen. 
      end
      

