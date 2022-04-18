      subroutine fsub(ncomp, x, u, f)
****************************************************************
*   USER PROVIDED FUNCTION
****************************************************************
*   FIRST ORDER ODE SYSTEM f' = F(t,u)
*
*   Solve the BVP [Example 4b]
*      U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2   x in ]0,L[
*      U(0,s) = s/(s^2+9)^2                                       <==> U(0,s) - s/(s^2+9)^2 = 0
*      U(L,s) = [(s*L+1)*sin(3*L) + 3*L*cos(3*L)] /6/(s^2+9) +
*             + [s*cos(3*L) - 3*sin(3*L)] / (s^2+9)^2             <==> U(L,s) - ... = 0
*
*   U1 = U(x,s)
*   Complex ODE system
*       U1' = U2             ==>   U2' = U"
*       U2' = s^2*U1 - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
* i.e.
*       U1' = U2
*       U2' = s^2*U1 - s*x*sin(3*x)/6 - sin(3*x)/6 - x*cos(3*x)/2
*
*
*  Complex boundary values
*       U1(0) = s/(s^2+9)^2
*       U1(L) = [(s*L+1)*sin(3*L)+3*L*cos(3*L)]/6/(s^2+9) + [s*cos(3*L)-3*sin(3*L)]/(s^2+9)^2
*
*
*   Real ODE system u' = f(x,u)   [U1=u1+i*u2, U2=u3+i*u4]
*       u1' = u3
*       u2' = u4
*       u3' = Re[s^2*U1 - s*x*sin(3*x)/6 - sin(3*x)/6 - x*cos(3*x)/2]
*       u4' = Im[s^2*U1 - s*x*sin(3*x)/6 - sin(3*x)/6 - x*cos(3*x)/2]
* i.e.
*       u1' = u3
*       u2' = u4
*       u3' = Re(s^2)*u1 - Im(s^2)*u2 - (Re(s)*x+1)*sin(3*x)/6 - x*cos(3*x)/2
*       u4' = Im(s^2)*u1 + Re(s^2)*u2 - Im(s)*x*sin(3*x)/6
*
*  Real boundary values
*       u1(0) = Re(s/(s^2+9)^2)   <==>   u1(0) - Re(s/(s^2+9)^2) = 0
*       u2(0) = Im(s/(s^2+9)^2)   <==>   u2(0) - Im(s/(s^2+9)^2) = 0
*          U1(L) = [(s*L+1)*sin(3*L)+3*L*cos(3*L)]/6/(s^2+9) + [s*cos(3*L)-3*sin(3*L)]/(s^2+9)^2
*       u1(L) = Re(U1(L))         <==>   u1(L) - Re(U1(L)) = 0
*       u2(L) = Im(U1(L))         <==>   u2(L) - Im(U1(L)) = 0
*
****************************************************************
!     subroutine fsub(ncomp, x, u, f)
      implicit double precision (a-h,o-z)
      dimension u(*), f(*)

      intrinsic sin, cos

!     maximum number of OpenMP threads and thread id
      integer maxTHRDS, myid, omp_get_thread_num
      parameter (maxTHRDS=20)

!     problem's and algorithm's parameters
      double precision L
      double complex st, s, ss
      dimension st(maxTHRDS)
      common/probs/ st, L

      common/counts/nfcall
      nfcall = nfcall + 1

!     OpenMP
#ifdef _OPENMP
        myid = omp_get_thread_num()
#else
        myid = 0
#endif

!     each thread has its parameter
      s = st(myid+1) ! FORTRAN index starts from 1

      f(1) = u(3) ! real part
      f(2) = u(4) ! imag part
	  ss = s*s
!     U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
      f(3) = RealPart(ss)*u(1) - ImagPart(ss)*u(2) -        ! real part of s^2*U
     +       (RealPart(s)*x+1)*sin(3*x)/6 - x*cos(3*x)/2    ! real part
      f(4) = ImagPart(ss)*u(1) + RealPart(ss)*u(2) -        ! imag part of s^2*U
     +	     ImagPart(s)*x*sin(3*x)/6                       ! imag part

      return
      end



      subroutine dfsub(ncomp, x, u, df)
****************************************************************
*   USER PROVIDED FUNCTION
****************************************************************
*  dfsub CALCULATES THE JACOBIAN MATRIX OF F AS DEFINED BY THE SUBROUTINE FSUB
*  part of sample test function for code twpbvp, for two-point
*  boundary value problems.
*
*   Real ODE system u' = f(x,u)   [U1=u1+i*u2, U2=u3+i*u4]
*       u1' = u3
*       u2' = u4
*       u3' = Re(s^2)*u1 - Im(s^2)*u2 - (Re(s)*x+1)*sin(3*x)/6 - x*cos(3*x)/2
*       u4' = Im(s^2)*u1 + Re(s^2)*u2 - Im(s)*x*sin(3*x)/6
*
*      |   0         0        1        0 |
*      |   0         0        0        1 |
*  J = | Re(s^2)  -Im(s^2)    0        0 |
*      | Im(s^2)   Re(s^2)    0        0 |
*
****************************************************************
!     subroutine dfsub(ncomp, x, u, df)
      implicit double precision (a-h,o-z)
      dimension u(*), df(ncomp,*)

!     maximum number of OpenMP threads and thread id
      integer maxTHRDS, myid, omp_get_thread_num
      parameter (maxTHRDS=20)

!     problem's and algorithm's parameters
      double precision L
      double complex st, s, ss
      dimension st(maxTHRDS)
      common/probs/ st, L

!     OpenMP
#ifdef _OPENMP
        myid = omp_get_thread_num()
#else
        myid = 0
#endif

!     each thread has its parameter
      s = st(myid+1) ! FORTRAN index starts from 1
      ss = s*s

      df(1,1) = 0.0d+0 ! 1st row
      df(1,2) = 0.0d+0
      df(1,3) = 1.0d+0
      df(1,4) = 0.0d+0

      df(2,1) = 0.0d+0 ! 2nd row
      df(2,2) = 0.0d+0
      df(2,3) = 0.0d+0
      df(2,4) = 1.0d+0

      df(3,1) =  RealPart(ss) ! 3rd row
      df(3,2) = -ImagPart(ss)
      df(3,3) = 0.0d+0
      df(3,4) = 0.0d+0

      df(4,1) = ImagPart(ss) ! 4th row
      df(4,2) = RealPart(ss)
      df(4,3) = 0.0d+0
      df(4,4) = 0.0d+0

      return
      end



      subroutine gsub(i, ncomp, u, g)
****************************************************************
*   USER PROVIDED FUNCTION
****************************************************************
*   i-TH BOUNDARY CONDITION
*  part of sample test function for code twpbvp, for two-point
*  boundary value problems.
*
*  Real boundary values
*  i=1  u1(0) = Re[s/(s^2+9)^2]   <==>   u1(0) - Re[s/(s^2+9)^2] = 0
*    2  u2(0) = Im[s/(s^2+9)^2]   <==>   u2(0) - Im[s/(s^2+9)^2] = 0
*    3  u1(L) = Re(U1(L))         <==>   u1(L) - Re(U1(L)) = 0
*    4  u2(L) = Im(U1(L))         <==>   u2(L) - Im(U1(L)) = 0
*  where
*  U1(L) = [(s*L+1)*sin(3*L)+3*L*cos(3*L)]/6/(s^2+9) + [s*cos(3*L)-3*sin(3*L)]/(s^2+9)^2
****************************************************************
      implicit double precision (a-h,o-z)
      dimension u(*)

!     maximum number of OpenMP threads and thread id
      integer maxTHRDS, myid, omp_get_thread_num
      parameter (maxTHRDS=20)

!     problem's and algorithm's parameters
      double precision L
      double complex st, s
      dimension st(maxTHRDS)
      common/probs/ st, L

!     OpenMP
#ifdef _OPENMP
        myid = omp_get_thread_num()
#else
        myid = 0
#endif

!     each thread has its parameter
      s = st(myid+1) ! FORTRAN index starts from 1

      SELECT CASE (i)
      CASE (1)
          g = u(1) - RealPart( s/(s*s+9.0d0)**2 )
      CASE (2)
          g = u(2) - ImagPart( s/(s*s+9.0d0)**2 )
      CASE (3)
          g = u(1) - RealPart( (3*L*cos(3*L) + (s*L+1)*sin(3*L))/6
     +             /(s*s+9) + (s*cos(3*L) - 3*sin(3*L))/(s*s+9)**2 )
      CASE (4)
          g = u(2) - ImagPart( (3*L*cos(3*L) + (s*L+1)*sin(3*L))/6
     +             /(s*s+9) + (s*cos(3*L) - 3*sin(3*L))/(s*s+9)**2 )
      CASE DEFAULT
          WRITE(*,*)  "This value is not valid"
      END SELECT

      return
      end



      subroutine dgsub(i, ncomp, u, dg)
****************************************************************
*   USER PROVIDED FUNCTION
****************************************************************
*   dgsub CALCULATES THE i-TH JACOBIAN MATRIX CORRESPONDING TO THE FUNCTION gsub
*  part of sample test function for code twpbvp, for two-point
*  boundary value problems.
*
****************************************************************
      implicit double precision (a-h,o-z)
      dimension u(*), dg(*)

      select case (i)
      case (1,3)
          dg(1) = 1.0d+0
          dg(2) = 0.0d+0
      case (2,4)
          dg(1) = 0.0d+0
          dg(2) = 1.0d+0
      case default
          write(*,*)  "This value is not valid"
      end select
      dg(3) = 0.0d+0
      dg(4) = 0.0d+0

      return
      end
