      subroutine fsub(ncomp, x, u, f)
****************************************************************
*   USER PROVIDED FUNCTION
****************************************************************
*   FIRST ORDER ODE SYSTEM f' = F(t,u)
*
*   Solve the BVP [Example 3b]
*       U" = s*U - u(x,0+)   x in ]0,L[
*                  u(x,0+) = x*(x-1)
*       U(0,s)=2/s^2              <==> U(0,1)-2/s^2 = 0
*       U(L,s)=2/s^2 + L*(L-1)/s  <==> U(L,s)-2/s^2-L*(L-1)/s = 0
*
*   U1 = U(x,s)
*   Complex ODE system
*       U1' = U2             ==>   U2' = U"
*       U2' = s*U1 - x*(x-1)
*
*  Complex boundary values
*       U1(0) = 2/s^2
*       U1(L) = 2/s^2 + L*(L-1)/s
*
*   Real ODE system U' = f(x,U)   [U1=u1+i*u2, U2=u3+i*u4]
*       u1' = u3
*       u2' = u4
*       u3' = Re(s)*u1 - Im(s)*u2 - x*(x-1)
*       u4' = Im(s)*u1 + Re(s)*u2
*
*  Real boundary values
*       u1(0) = Re(2/s^2)               <==> u1(0) - Re(2/s^2) = 0
*       u2(0) = Im(2/s^2)               <==> u2(0) - Im(2/s^2) = 0
*       u3(L) = Re( 2/s^2 + L*(L-1)/s ) <==> u3(L) - Re( 2/s^2 + L*(L-1)/s ) = 0 
*       u4(L) = Im( 2/s^2 + L*(L-1)/s ) <==> u4(L) - Im( 2/s^2 + L*(L-1)/s ) = 0
*
****************************************************************
!     subroutine fsub(ncomp, x, u, f)
      implicit double precision (a-h,o-z)
      dimension u(*), f(*)

!     problem's and algorithm's parameters
      double precision L
      double complex s
      common/probs/ s, L

      common/counts/nfcall
      nfcall = nfcall + 1

      f(1) = u(3) ! real part
      f(2) = u(4) ! imag part
      f(3) = RealPart(s)*u(1) - ImagPart(s)*u(2) - x*(x-1) ! real part of s*U1
      f(4) = ImagPart(s)*u(1) + RealPart(s)*u(2)           ! imag part of s*U1
 
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
*   Real ODE system U' = f(x,U)   [U1=u1+i*u2, U2=u3+i*u4]
*       u1' = u3
*       u2' = u4
*       u3' = Re(s)*u1 - Im(s)*u2 - x*(x-1)
*       u4' = Im(s)*u1 + Re(s)*u2
*
*      |   0      0      1      0 |
*      |   0      0      0      1 |
*  J = | Re(s) -Im(s)    0      0 |
*      | Im(s)  Re(s)    0      0 |
*
****************************************************************
!     subroutine dfsub(ncomp, x, u, df)
      implicit double precision (a-h,o-z)
      dimension u(*), df(ncomp,*)

!     problem's and algorithm's parameters
      double precision L
      double complex s
      common/probs/ s, L

      df(1,1) = 0.0d+0
      df(1,2) = 0.0d+0
      df(1,3) = 1.0d+0
      df(1,4) = 0.0d+0

      df(2,1) = 0.0d+0
      df(2,2) = 0.0d+0
      df(2,3) = 0.0d+0
      df(2,4) = 1.0d+0

      df(3,1) =  RealPart(s)
      df(3,2) = -ImagPart(s)
      df(3,3) = 0.0d+0
      df(3,4) = 0.0d+0

      df(4,1) = ImagPart(s)
      df(4,2) = RealPart(s)
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
*  Real boundary values
*  i=1  u1(0) = Re(2/s^2)               <==> u1(0) - Re(2/s^2) = 0
*    2  u2(0) = Im(2/s^2)               <==> u2(0) - Im(2/s^2) = 0
*    3  u3(L) = Re( 2/s^2 + L*(L-1)/s ) <==> u3(L) - Re( 2/s^2 + L*(L-1)/s ) = 0
*    4  u4(L) = Im( 2/s^2 + L*(L-1)/s ) <==> u4(L) - Im( 2/s^2 + L*(L-1)/s ) = 0
****************************************************************
      implicit double precision (a-h,o-z)
      dimension u(*)
!     problem's and algorithm's parameters
      double precision L
      double complex s, ss
      common/probs/ s, L

      ss = s**2
      SELECT CASE (i)
      CASE (1)
          g = u(1) - RealPart(2.0d+0/ss)
      CASE (2)
          g = u(2) - ImagPart(2.0d+0/ss)
      CASE (3)
          g = u(1) - RealPart(2.0d+0/ss + L*(L-1)/s)
      CASE (4)
          g = u(2) - ImagPart(2.0d+0/ss + L*(L-1)/s)
      CASE DEFAULT
          WRITE(*,*)  "Hmmmm, I don't know"
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
*  Real boundary values
*  i=1  u1(0) = Re(2/s^2)               <==> u1(0) - Re(2/s^2) = 0
*    2  u2(0) = Im(2/s^2)               <==> u2(0) - Im(2/s^2) = 0
*    3  u1(L) = Re( 2/s^2 + L*(L-1)/s ) <==> u1(L) - Re( 2/s^2 + L*(L-1)/s ) = 0 
*    4  u2(L) = Im( 2/s^2 + L*(L-1)/s ) <==> u2(L) - Im( 2/s^2 + L*(L-1)/s ) = 0
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
          write(*,*)  "Hmmmm, I don't know"
      end select
      dg(3) = 0.0d+0
      dg(4) = 0.0d+0

      return
      end
