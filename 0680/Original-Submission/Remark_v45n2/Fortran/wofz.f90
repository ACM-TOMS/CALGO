  MODULE wofzmod

    CONTAINS
!      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 16, NO. 1, PP. 47.
! Modified March 2015:
!  + Made elemental
!
!  + Parameterized precision
!
!  + Generic intrinsics
!
!  + Parameterized constants
!
!  + Some machine dependent parameters generated in terms of numerical
!    enquiry functions.
!
! Routine gives 12+ significent figures on a wide range of values.
!
      ELEMENTAL SUBROUTINE wofz(z,w,flag)
!  Given a complex number z = (xi,yi), this subroutine computes
!  the value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
!  where erfc is the complex complementary error-function and i
!  means sqrt(-1).
!  The accuracy of the algorithm for z in the 1st and 2nd quadrant
!  is 14 significant digits; in the 3rd and 4th it is 13 significant
!  digits outside a circular region with radius 0.126 around a zero
!  of the function.
!  All real variables in the program are double precision.
!  The code contains a few compiler-dependent parameters :
!     rmaxreal = the maximum value of rmaxreal equals the root of
!                rmax = the largest number which can still be
!                implemented on the computer in double precision
!                Floating-point arithmetic
!     rmaxexp  = ln(rmax) - ln(2)
!     rmaxgoni = the largest possible argument of a double precision
!                goniometric function (cos, sin, ...)
!  The reason why these parameters are needed as they are defined will
!  be explained in the code by means of comments
!
!  PARAMETER LIST
!  --------------
!      Z (COMPLEX (wp))  -- input value
!      W (COMPLEX (wp))  -- value of w(z) = exp(-z**2)*erfc(-i*z)
!   FLAG (LOGICAL)       -- an error flag indicating whether overflow will
!                           occur or not; type logical;
!                           the values of this variable have the following
!                           meaning :
!                           FLAG=.false. : no error condition
!                           FLAG=.true.  : overflow will occur, the routiNE
!                                          becomes inactive
!  Furthermore the parameter factor equals 2/sqrt(pi)
!  The routine is not underflow-protected but any variable can be
!  put to 0 upon underflow;
!  Reference - GPM Poppe, CMJ Wijers; More efficient computation of
!  the complex error-function, ACM Trans. Math. Software.
! .. Use Statements ..
        USE set_rk, ONLY : wp=>rk
! ..
! .. Parameters ..
        REAL (wp), PARAMETER :: zero = 0.0E0_wp, two = 2.0E0_wp
        REAL (wp), PARAMETER :: half = 0.5E0_wp, one = 1.0E0_wp, &
                                factor = one/SQRT(ATAN(one))
        REAL (wp), PARAMETER :: rmaxexp = LOG(HUGE(one))-LOG(two)
        REAL (wp), PARAMETER :: rmaxgoni = 3.53711887601422D+15
        REAL (wp), PARAMETER :: rmaxreal = SQRT(HUGE(one))
        REAL (wp), PARAMETER :: xScale = 6.3E0_wp, yScale = 4.4E0_wp
        REAL (wp), PARAMETER :: bound = 0.085264E0_wp
        REAL (wp), PARAMETER :: ps(3) = [0.85E0_wp, 6.0E0_wp, 72.0E0_wp]
        REAL (wp), PARAMETER :: te(9) = [3.0E0_wp, 1442.0E0_wp, 26.0E0_wp, &
                 77.0E0_wp, 1.88E0_wp, 7.0E0_wp, 34.0E0_wp, 16.0E0_wp, 26.0E0_wp]
! ..
! .. Scalar Arguments ..
        COMPLEX (wp), INTENT (OUT) :: w
        COMPLEX (wp), INTENT (IN) :: z
        LOGICAL, INTENT (OUT) :: flag
! ..
! .. Local Scalars ..
        REAL (wp) :: c, daux, h, h2, qlambda, qrho, rx, ry, sx, sy, tx, ty, u, &
          u1, u2, v, v1, v2, w1, x, xabs, xabsq, xaux, xi, xquad, xsum, y, &
          yabs, yi, yquad, ysum
        INTEGER :: i, j, kapn, n, np1, nu
        LOGICAL :: a, b
! ..
! .. Intrinsic Functions ..
        INTRINSIC aimag, cmplx, abs, cos, exp, sin, sqrt, int, real
! ..
        flag = .FALSE.

        xi = real(z,kind=wp)
        yi = aimag(z)
        xabs = abs(xi)
        yabs = abs(yi)
        x = xabs/xScale
        y = yabs/yScale


!     The following if-statement protects
!     qrho = (x**2 + y**2) against overflow

        IF ((xabs>rmaxreal) .OR. (yabs>rmaxreal)) GO TO 100

        qrho = x**2 + y**2

        xabsq = xabs**2
        xquad = xabsq - yabs**2
        yquad = two*xabs*yabs

        a = qrho < bound

        IF (a) THEN

!  If (qrho.lt.0.085264d0) then the Faddeeva-function is evaluated
!  using a power-series (Abramowitz/Stegun, equation (7.1.5), p.297)
!  n is the minimum number of terms needed to obtain the required
!  accuracy

          qrho = (one-ps(1)*y)*sqrt(qrho)
          n = nint(ps(2)+ps(3)*qrho)
          j = 2*n + 1
          xsum = one/REAL(j, KIND=wp)
          ysum = zero
          DO i = n, 1, -1
            j = j - 2
            xaux = (xsum*xquad-ysum*yquad)/REAL(i, KIND=wp)
            ysum = (xsum*yquad+ysum*xquad)/REAL(i, KIND=wp)
            xsum = xaux + one/REAL(j, KIND=wp)
          END DO
          u1 = -factor*(xsum*yabs+ysum*xabs) + one
          v1 = factor*(xsum*xabs-ysum*yabs)
          daux = exp(-xquad)
          u2 = daux*cos(yquad)
          v2 = -daux*sin(yquad)

          u = u1*u2 - v1*v2
          v = u1*v2 + v1*u2

        ELSE

!  If (qrho.gt.1.o) then w(z) is evaluated using the Laplace
!  continued fraction
!  nu is the minimum number of terms needed to obtain the required
!  accuracy

!  If ((qrho.gt.0.085264d0).and.(qrho.lt.1.0)) then w(z) is evaluated
!  by a truncated Taylor expansion, where the Laplace continued fraction
!  is used to calculate the derivatives of w(z)
!  kapn is the minimum number of terms in the Taylor expansion needed
!  to obtain the required accuracy
!  nu is the minimum number of terms of the continued fraction needed
!  to calculate the derivatives with the required accuracy


          IF (qrho>one) THEN
            h = zero
            kapn = 0
            qrho = sqrt(qrho)
            nu = int(te(1)+(te(2)/(te(3)*qrho+te(4))))
          ELSE
            qrho = (one-y)*sqrt(one-qrho)
            h = te(5)*qrho
            h2 = two*h
            kapn = nint(te(6)+te(7)*qrho)
            nu = nint(te(8)+te(9)*qrho)
          END IF

          b = (h>zero)

          IF (b) qlambda = h2**kapn

          rx = zero
          ry = zero
          sx = zero
          sy = zero

          DO n = nu, 0, -1
            np1 = n + 1
            tx = yabs + h + np1*rx
            ty = xabs - np1*ry
            c = half/(tx**2+ty**2)
            rx = c*tx
            ry = c*ty
            IF ((b) .AND. (n<=kapn)) THEN
              tx = qlambda + sx
              sx = rx*tx - ry*sy
              sy = ry*tx + rx*sy
              qlambda = qlambda/h2
            END IF
          END DO

          IF (h==zero) THEN
            u = factor*rx
            v = factor*ry
          ELSE
            u = factor*sx
            v = factor*sy
          END IF

          IF (yabs==zero) u = exp(-xabs**2)

        END IF

!  Evaluation of w(z) in the other quadrants

        IF (yi<zero) THEN

          IF (a) THEN
            u2 = two*u2
            v2 = two*v2
          ELSE
            xquad = -xquad


!  The following if-statement protects 2*exp(-z**2) against overflow

            IF ((yquad>rmaxgoni) .OR. (xquad>rmaxexp)) GO TO 100

            w1 = two*exp(xquad)
            u2 = w1*cos(yquad)
            v2 = -w1*sin(yquad)
          END IF

          u = u2 - u
          v = v2 - v
          IF (xi>zero) v = -v
        ELSE
          IF (xi<zero) v = -v
        END IF

        w = cmplx(u,v,kind=wp)
        RETURN

100     flag = .TRUE.
        RETURN

      END SUBROUTINE wofz
    END MODULE wofzmod

