      program spgma1

C     SPG Test driver. This master file uses SPG to solve
C     an easy (very easy) bound constrained problem.
C
C     This version 02 FEB 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.
C     Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

C     PARAMETERS
      integer nmax
      parameter (nmax=100)

C     COMMON ARRAYS
      double precision l(nmax),u(nmax)

C     LOCAL SCALARS
      double precision pginfn,pgtwon,eps,eps2,f
      integer fcnt,flag,gcnt,i,iter,m,maxfc,maxit,n
      logical output

C     LOCAL ARRAYS
      double precision x(nmax)

C     EXTERNAL SUBROUTINES
      external spg

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

C     COMMON BLOCKS
      common /bounds/l,u

C     DEFINE DIMENSION OF THE PROBLEM

      n = 10

C     DEFINE BOUNDS

      do i = 1,n
          l(i)= -100.0d0
          u(i)=   50.0d0
      end do

C     DEFINE INITIAL POINT

      do i = 1,n
          x(i) = 60.0d0
      end do

C     SET UP OPTIMIZER PARAMETERS 

      output= .true.
      maxit = 1000
      maxfc = 2000
      eps = 0.0d0
      eps2 = 1.0d-6
      m = 10

C     CALL THE OPTIMIZER

      call spg(n,x,m,eps,eps2,maxit,maxfc,output,f,pginfn,pgtwon,iter,
     +         fcnt,gcnt,flag)

C     WRITE STATISTICS

      write (*,fmt=9010) f,pginfn,sqrt(pgtwon),flag
      write (*,fmt=9020) iter,fcnt,gcnt

      stop

 9010 format (/' F= ',8X,1P,D17.10,/' PGINFNORM= ',1X,1P,D16.10,
     +       /' PGTWONORM= ',1X,1P,D16.10,/' FLAG= ',6X,I1)
 9020 format (/' ITER= ',1X,I10,/' FCNT= ',1X,I10,/' GCNT= ',1X,I10)

      end


      subroutine evalf(n,x,f,inform)

C     This subroutine computes the objective function.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point at which the function will be evaluated.
C
C     On Return
C
C     f     double precision,
C           function value at x,
C
C     inform integer,
C           termination parameter:
C           0= the function was successfully evaluated,
C           1= some error occurs in the function evaluation.

C     SCALAR ARGUMENTS
      double precision f
      integer n,inform

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i

      inform = 0

      f = 0.0d0
      do i = 1,n
          f = f + x(i)**2
      end do

      end


      subroutine evalg(n,x,g,inform)

C     This subroutine computes the gradient of the 
C     objective function.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point at which the gradient will be evaluated.
C
C     On Return
C
C     g     double precision g(n),
C           gradient vector at x,
C
C     inform integer,
C           termination parameter:
C           0= the gradient was successfully evaluated,
C           1= some error occurs in the gradient evaluation.

C     SCALAR ARGUMENTS
      integer n,inform

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

      inform = 0

      do i = 1,n
          g(i) = 2.0d0 * x(i)
      end do

      end


      subroutine proj(n,x,inform)

C     This subroutine computes the projection of an arbitrary 
C     point onto the feasible set. Since the feasible set can be 
C     described in many ways, its information is in a common 
C     block.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point that will be projected.
C
C     On Return
C
C     x     double precision x(n),
C           projected point,
C
C     inform integer,
C           termination parameter:
C           0= the projection was successfully done,
C           1= some error occurs in the projection.

C     PARAMETERS
      integer nmax
      parameter (nmax=100)

C     COMMON ARRAYS
      double precision l(nmax),u(nmax)

C     SCALAR ARGUMENTS
      integer n,inform

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      integer i

C     COMMON BLOCKS
      common /bounds/l,u

      inform = 0

      do i = 1,n
          x(i) = max(l(i),min(x(i),u(i)))
      end do

      end
