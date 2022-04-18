MODULE wofz_v2_mod
  USE set_rk, ONLY : rk
  PUBLIC :: daw715


  REAL (rk), PARAMETER :: zero = 0.0E0_rk, one = 1.0E0_rk, two = 2.0E0_rk
  REAL (rk), PARAMETER :: half = 0.5E0_rk, one_third= 0.333333333333333_rk, &
    two_thirds=0.6666666666666667_rk, four_tenths=0.4_rk, two_seventh=0.285714285714285_rk, &
    factor = 1.128379167095513e+000_rk   !2/sqrt(pi)   !one/SQRT(ATAN(one))
  REAL (rk), PARAMETER :: two_sqrt_pi = 1.128379167095513e+000_rk   !2/sqrt(pi)
  REAL (rk), PARAMETER :: rmaxexp = LOG(HUGE(one))-LOG(two)
  REAL (rk), PARAMETER :: rmaxgoni = 3.53711887601422e+15_rk
  REAL(rk), PARAMETER ::  sixtenth_Rmax=1.844674407370955e+019_rk  !!sqrt(sqrt(sqrt(Huge(1.0_rk))))
  REAL(rk), PARAMETER ::  Root4_Rmax=sixtenth_Rmax*sixtenth_Rmax*sixtenth_Rmax*sixtenth_Rmax
  REAL (rk), PARAMETER :: rmaxreal =Root4_Rmax*Root4_Rmax!SQRT(HUGE(one))
  REAL (rk), PARAMETER :: xScale = 6.3E0_rk, yScale = 4.4E0_rk
  REAL (rk), PARAMETER :: bound = 0.085264E0_rk
  REAL (rk), PARAMETER :: ps(3) = (/0.85E0_rk, 6.0E0_rk, 72.0E0_rk/)
  REAL (rk), PARAMETER :: te(9) = (/3.0E0_rk, 1442.0E0_rk, 26.0E0_rk, &
    77.0E0_rk, 1.88E0_rk, 7.0E0_rk, 34.0E0_rk, 16.0E0_rk, 26.0E0_rk/)
  ! constants used in the present correction
  REAL (rk), PARAMETER :: yabs_corr=0.031623_rk, xabs_corr_min=1.8396_rk, &
    xabs_corr_max=20.0_rk

  COMPLEX(rk), PARAMETER :: j1= (0.0e0_rk,1.0e0_rk)

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
  ELEMENTAL SUBROUTINE wofz_v2(z,w,flag)
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
    !      Z (COMPLEX (rk))  -- input value
    !      W (COMPLEX (rk))  -- value of w(z) = exp(-z**2)*erfc(-i*z)
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

    ! ..
    ! .. Parameters ..


    ! ..
    ! .. Scalar Arguments ..
    COMPLEX (rk), INTENT (OUT) :: w
    COMPLEX (rk), INTENT (IN) :: z
    LOGICAL, INTENT (OUT) :: flag
    ! ..
    ! .. Local Scalars ..
    REAL (rk) :: c, daux, h, h2, qlambda, qrho, rx, ry, sx, sy, tx, ty, u, &
      u1, u2, v, v1, v2, w1, x, xabs, xabsq, xaux, xi, xquad, xsum, y, &
      yabs, yi, yquad, ysum
    INTEGER :: i, j, kapn, n, np1, nu
    LOGICAL :: a, b
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC AIMAG, CMPLX, ABS, COS, EXP, SIN, SQRT, INT, REAL
    ! ..
    flag = .FALSE.

    xi = REAL(z,KIND=rk)
    yi = AIMAG(z)
    xabs = ABS(xi)
    yabs = ABS(yi)

    x = xabs/xScale
    y = yabs/yScale


    !     The following if-statement protects
    !     qrho = (x**2 + y**2) against overflow

    !    IF ((xabs>rmaxreal) .OR. (yabs>rmaxreal)) GO TO 100 !original statement
    !*** modified
    IF ((xabs>rmaxreal) .OR. (yabs>rmaxreal)) THEN
      flag = .TRUE.
      RETURN
    END IF
    !***


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

      qrho = (one-ps(1)*y)*SQRT(qrho)
      n = NINT(ps(2)+ps(3)*qrho)
      j = 2*n + 1
      xsum = one/REAL(j, KIND=rk)
      ysum = zero
      DO i = n, 1, -1
        j = j - 2
        xaux = (xsum*xquad-ysum*yquad)/REAL(i, KIND=rk)
        ysum = (xsum*yquad+ysum*xquad)/REAL(i, KIND=rk)
        xsum = xaux + one/REAL(j, KIND=rk)
      END DO
      u1 = -factor*(xsum*yabs+ysum*xabs) + one
      v1 = factor*(xsum*xabs-ysum*yabs)
      daux = EXP(-xquad)
      u2 = daux*COS(yquad)
      v2 = -daux*SIN(yquad)

      u = u1*u2 - v1*v2
      v = u1*v2 + v1*u2

      !%%%%%%%%%%%%% present modification %%%%%%%%%
      !  If (1.8396<=xabs<=20.0) and yabs<0.0316227766 the Faddeeva-function is evaluated
      !  using "upward" truncated Taylor series (Taylor expansion about z0=x)
      !
    ELSEIF (yabs<=yabs_corr .and. xabs>=xabs_corr_min .and. xabs<=xabs_corr_max) THEN
      w=use_taylor_daw(xi,yabs,7)
      IF (yi<zero)THEN
        w=two*EXP(yi*yi-xi*xi+j1*two*yabs*xi)-(REAL(w,KIND=rk)-j1*AIMAG(w))
      ENDIF
      RETURN
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        qrho = SQRT(qrho)
        nu = INT(te(1)+(te(2)/(te(3)*qrho+te(4))))
      ELSE
        qrho = (one-y)*SQRT(one-qrho)
        h = te(5)*qrho
        h2 = two*h
        kapn = NINT(te(6)+te(7)*qrho)
        nu = NINT(te(8)+te(9)*qrho)
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

      IF (yabs==zero) u = EXP(-xabs**2)

    END IF

    !  Evaluation of w(z) in the other quadrants

    IF (yi<zero) THEN

      IF (a) THEN
        u2 = two*u2
        v2 = two*v2
      ELSE
        xquad = -xquad


        !  The following if-statement protects 2*exp(-z**2) against overflow

        !       IF ((yquad>rmaxgoni) .OR. (xquad>rmaxexp)) GO TO 100  !original statement
        !*** modified
        IF ((yquad>rmaxgoni) .OR. (xquad>rmaxexp))THEN
          flag = .TRUE.
          RETURN
        END IF
        !***


        w1 = two*EXP(xquad)
        u2 = w1*COS(yquad)
        v2 = -w1*SIN(yquad)
      END IF

      u = u2 - u
      v = v2 - v
      IF (xi>zero) v = -v
    ELSE
      IF (xi<zero) v = -v
    END IF

    w = CMPLX(u,v,KIND=rk)
    RETURN

    100     flag = .TRUE.
    RETURN

  END SUBROUTINE wofz_v2





  !!:--
  ELEMENTAL FUNCTION use_taylor_daw ( xx, yy, nterms ) RESULT ( w )
    REAL(rk), INTENT(IN) :: xx, yy
    INTEGER, INTENT(IN) :: nterms
    COMPLEX(rk) ::  w
    REAL(rk) ::  x_sqr,y_sqr, ysqr_minus_xsqr,  two_xy, &
      dw,dw_1,dw_2,dw_3,dw_4,dw_5,dw_6,dw_7, y_4

    x_sqr=xx*xx
    y_sqr=yy*yy
    ysqr_minus_xsqr=y_sqr-x_sqr

    dw=daw715(xx)
    dw_1=(one-two*xx*dw)
    two_xy=two*xx*yy

    SELECT CASE (nterms)
      CASE (1); !    nterms==1
        w=EXP(ysqr_minus_xsqr)*(one-j1*two_xy)
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy)
      CASE (2); !    nterms==2
        dw_2=-(xx*dw_1+dw)
        w=EXP(ysqr_minus_xsqr)*(one-j1*two_xy)
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr)

      CASE (3); !    nterms==3
        dw_2=-(xx*dw_1+dw)
        dw_3=-two_thirds*(xx*dw_2+dw_1)
        w=EXP(ysqr_minus_xsqr)*(COS(two_xy)-j1*SIN(two_xy))
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy)
      CASE (4); !    nterms==4
        dw_2=-(xx*dw_1+dw)
        dw_3=-two_thirds*(xx*dw_2+dw_1)
        dw_4=-half*(xx*dw_3+dw_2)
        w=EXP(ysqr_minus_xsqr)*(COS(two_xy)-j1*SIN(two_xy))
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy+&
          dw_4*y_sqr*y_sqr);
      CASE (5); !    nterms==5
        dw_2=-(xx*dw_1+dw)
        dw_3=-two_thirds*(xx*dw_2+dw_1)
        dw_4=-half*(xx*dw_3+dw_2)
        dw_5=-four_tenths*(xx*dw_4+dw_3)
        w=EXP(ysqr_minus_xsqr)*(COS(two_xy)-j1*SIN(two_xy))
        y_4=y_sqr*y_sqr;
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy+&
          dw_4*y_4+j1*dw_5*yy*y_4)
      CASE (6); !    nterms==6
        dw_2=-(xx*dw_1+dw)
        dw_3=-two_thirds*(xx*dw_2+dw_1)
        dw_4=-half*(xx*dw_3+dw_2)
        dw_5=-four_tenths*(xx*dw_4+dw_3)
        dw_6=-one_third*(xx*dw_5+dw_4)
        w=EXP(ysqr_minus_xsqr)*(COS(two_xy)-j1*SIN(two_xy))
        y_4=y_sqr*y_sqr;
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy+&
          dw_4*y_4+j1*dw_5*yy*y_4-dw_6*y_sqr*y_4)
      CASE (7); !    nterms==7
        dw_2=-(xx*dw_1+dw)
        dw_3=-two_thirds*(xx*dw_2+dw_1)
        dw_4=-half*(xx*dw_3+dw_2)
        dw_5=-four_tenths*(xx*dw_4+dw_3)
        dw_6=-one_third*(xx*dw_5+dw_4)
        dw_7=-two_seventh*(xx*dw_6+dw_5)
        w=EXP(ysqr_minus_xsqr)*(COS(two_xy)-j1*SIN(two_xy))
        y_4=y_sqr*y_sqr;
        w=w+(two_sqrt_pi*j1)*(dw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy+&
          dw_4*y_4+j1*dw_5*yy*y_4-dw_6*y_sqr*y_4-j1*dw_7*y_4*y_sqr*yy)
    ENDSELECT
  END FUNCTION use_taylor_daw

  ELEMENTAL FUNCTION daw715 ( xx ) RESULT(daw)

    !*****************************************************************************80
    !
    !! DAW EVALUATES DAWSON'S INTEGRAL
    !
    !  DISCUSSION:
    !
    !    THIS FUNCTION EVALUATES DAWSON'S INTEGRAL,
    !      F(X) = EXP ( -X^2 ) * INTEGRAL ( 0 <= T <= X ) EXP ( T^2 ) DT
    !    FOR A REAL ARGUMENT X.
    !
    !    THE CALLING SEQUENCE FOR THIS FUNCTION IS
    !      Y=DAW(X)
    !    THE MAIN COMPUTATION USES RATIONAL CHEBYSHEV APPROXIMATIONS
    !    PUBLISHED IN MATH. COMP. 24, 171-178 (1970) BY CODY, PACIOREK
    !    AND THACHER.  THIS TRANSPORTABLE PROGRAM IS PATTERNED AFTER THE
    !    MACHINE-DEPENDENT FUNPACK PROGRAM DDAW(X), BUT CANNOT MATCH THAT
    !    VERSION FOR EFFICIENCY OR ACCURACY.  THIS VERSION USES RATIONAL
    !    APPROXIMATIONS THAT ARE THEORETICALLY ACCURATE TO ABOUT 19
    !    SIGNIFICANT DECIMAL DIGITS.  THE ACCURACY ACHIEVED DEPENDS ON THE
    !    ARITHMETIC SYSTEM, THE COMPILER, THE INTRINSIC FUNCTIONS, AND
    !    PROPER SELECTION OF THE MACHINE-DEPENDENT CONSTANTS.
    !
    !  MODIFIED:
    !
    !    10 JANUARY 2016
    !
    !  AUTHOR:
    !
    !    ORIGINAL FORTRAN77 VERSION BY WILLIAM CODY.
    !    FORTRAN90 VERSION BY JOHN BURKARDT.
    !
    ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS.  LET
    !
    !   XINF   = LARGEST POSITIVE MACHINE NUMBER
    !   XMIN   = THE SMALLEST POSITIVE MACHINE NUMBER.
    !   EPS    = SMALLEST POSITIVE NUMBER SUCH THAT 1+EPS > 1.
    !    APPROXIMATELY  BETA**(-P), WHERE BETA IS THE MACHINE
    !    RADIX AND P IS THE NUMBER OF SIGNIFICANT BASE-BETA
    !    DIGITS IN A FLOATING-POINT NUMBER.
    !
    ! THEN THE FOLLOWING MACHINE-DEPENDENT CONSTANTS MUST BE DECLARED
    !   IN DATA STATEMENTS.  IEEE VALUES ARE PROVIDED AS A DEFAULT.
    !
    !   XMAX   = ABSOLUTE ARGUMENT BEYOND WHICH DAW(X) UNDERFLOWS.
    !    XMAX = MIN(0.5/XMIN, XINF).
    !   XSMALL = ABSOLUTE ARGUMENT BELOW DAW(X)  MAY BE REPRESENTED
    !    BY X.  WE RECOMMEND XSMALL = SQRT(EPS).
    !   XLARGE = ARGUMENT BEYOND WHICH DAW(X) MAY BE REPRESENTED BY
    !    1/(2X).  WE RECOMMEND XLARGE = 1/SQRT(EPS).
    !
    ! ERROR RETURNS
    !
    !  THE PROGRAM RETURNS 0.0 FOR |X| > XMAX.
    !

    REAL(rk), INTENT(IN) :: xx
    REAL ( rk):: daw
    INTEGER:: i
    REAL ( rk ) &
      frac,sump,sumq,w2,x,y
    !  DIMENSION P1(10),P2(10),P3(10),P4(10),Q1(10),Q2(9),Q3(9),Q4(9)
    !
    !  MATHEMATICAL CONSTANTS.
    !
    REAL (rk),PARAMETER ::six25=6.25E0_rk
    REAL (rk),PARAMETER ::one225=12.25E0_rk
    REAL (rk),PARAMETER ::two5=25.0E0_rk
    !
    !  MACHINE-DEPENDENT CONSTANTS
    !
    REAL (rk),PARAMETER ::xsmall=1.05E-08_rk
    REAL (rk),PARAMETER ::xlarge=9.49E+07_rk
    REAL (rk),PARAMETER ::xmax=HUGE(1.0_rk)
    !
    !  COEFFICIENTS FOR R(9,9) APPROXIMATION FOR  |X| < 2.5
    !
    REAL (rk), PARAMETER :: p1(10) =(/-2.69020398788704782410E-12_rk, 4.18572065374337710778E-10_rk, &
      -1.34848304455939419963E-08_rk, 9.28264872583444852976E-07_rk, &
      -1.23877783329049120592E-05_rk, 4.07205792429155826266E-04_rk, &
      -2.84388121441008500446E-03_rk, 4.70139022887204722217E-02_rk, &
      -1.38868086253931995101E-01_rk, 1.00000000000000000004E+00_rk/)

    REAL (rk), PARAMETER :: q1(10)=(/1.71257170854690554214E-10_rk, 1.19266846372297253797E-08_rk, &
      4.32287827678631772231E-07_rk, 1.03867633767414421898E-05_rk, &
      1.78910965284246249340E-04_rk, 2.26061077235076703171E-03_rk, &
      2.07422774641447644725E-02_rk, 1.32212955897210128811E-01_rk, &
      5.27798580412734677256E-01_rk, 1.00000000000000000000E+00_rk/)
    !
    !  COEFFICIENTS FOR R(9,9) APPROXIMATION IN J-FRACTION FORM
    !     FOR  X IN [2.5, 3.5)
    !
    REAL (rk), PARAMETER :: p2(10)=(/-1.70953804700855494930E+00_rk,-3.79258977271042880786E+01_rk, &
      2.61935631268825992835E+01_rk, 1.25808703738951251885E+01_rk, &
      -2.27571829525075891337E+01_rk, 4.56604250725163310122E+00_rk, &
      -7.33080089896402870750E+00_rk, 4.65842087940015295573E+01_rk, &
      -1.73717177843672791149E+01_rk, 5.00260183622027967838E-01_rk/)

    REAL (rk), PARAMETER :: q2(9)=(/1.82180093313514478378E+00_rk, 1.10067081034515532891E+03_rk, &
      -7.08465686676573000364E+00_rk, 4.53642111102577727153E+02_rk, &
      4.06209742218935689922E+01_rk, 3.02890110610122663923E+02_rk, &
      1.70641269745236227356E+02_rk, 9.51190923960381458747E+02_rk, &
      2.06522691539642105009E-01_rk/)
    !
    !  COEFFICIENTS FOR R(9,9) APPROXIMATION IN J-FRACTION FORM
    !     FOR  X IN [3.5, 5.0]
    !
    REAL (rk), PARAMETER :: p3(10)=(/-4.55169503255094815112E+00_rk,-1.86647123338493852582E+01_rk, &
      -7.36315669126830526754E+00_rk,-6.68407240337696756838E+01_rk, &
      4.84507265081491452130E+01_rk, 2.69790586735467649969E+01_rk, &
      -3.35044149820592449072E+01_rk, 7.50964459838919612289E+00_rk, &
      -1.48432341823343965307E+00_rk, 4.99999810924858824981E-01_rk/)

    REAL (rk), PARAMETER :: q3(9)=(/ 4.47820908025971749852E+01_rk, 9.98607198039452081913E+01_rk, &
      1.40238373126149385228E+01_rk, 3.48817758822286353588E+03_rk, &
      -9.18871385293215873406E+00_rk, 1.24018500009917163023E+03_rk, &
      -6.88024952504512254535E+01_rk,-2.31251575385145143070E+00_rk, &
      2.50041492369922381761E-01_rk/)
    !
    !  COEFFICIENTS FOR R(9,9) APPROXIMATION IN J-FRACTION FORM
    !  FOR  |X| > 5.0
    !
    REAL (rk), PARAMETER :: p4(10)=(/-8.11753647558432685797E+00_rk,-3.84043882477454453430E+01_rk, &
      -2.23787669028751886675E+01_rk,-2.88301992467056105854E+01_rk, &
      -5.99085540418222002197E+00_rk,-1.13867365736066102577E+01_rk, &
      -6.52828727526980741590E+00_rk,-4.50002293000355585708E+00_rk, &
      -2.50000000088955834952E+00_rk, 5.00000000000000488400E-01_rk/)

    REAL (rk), PARAMETER :: q4(9)=(/ 2.69382300417238816428E+02_rk, 5.04198958742465752861E+01_rk, &
      6.11539671480115846173E+01_rk, 2.08210246935564547889E+02_rk, &
      1.97325365692316183531E+01_rk,-1.22097010558934838708E+01_rk, &
      -6.99732735041547247161E+00_rk,-2.49999970104184464568E+00_rk, &
      7.49999999999027092188E-01_rk/)

    x = ABS(xx)
    y=xx*xx
    IF (x > xlarge) THEN

      IF (x <= xmax) THEN
        daw = half / xx
      ELSE
        daw = zero
      END IF

    ELSE IF (x < xsmall) THEN

      daw = xx

    ELSEIF (5.0_rk<=x .and. xlarge>=x)THEN
      !
      !  5.0 <= ABS(X) <= XLARGE
      !
      w2 = one / y
      frac = zero
      DO i = 1, 9
        frac = q4(i) / (p4(i) + y + frac)
      END DO
      frac = p4(10) + frac
      daw = (half + half * w2 * frac) / xx

    ELSE

      IF (y < six25) THEN
        !
        !  ABS(X) < 2.5
        !
        sump = p1(1)
        sumq = q1(1)
        DO i = 2, 10
          sump = sump * y + p1(i)
          sumq = sumq * y + q1(i)
        END DO
        daw = xx * sump / sumq

      ELSE IF (y < one225) THEN
        !
        !  2.5 <= ABS(X) < 3.5
        !
        frac = zero
        DO i = 1, 9
          frac = q2(i) / (p2(i) + y + frac)
        END DO
        daw = (p2(10) + frac) / xx

      ELSE IF (y < two5) THEN
        !
        !  3.5 <= ABS(X) < 5.0
        !
        frac = zero
        DO i = 1, 9
          frac = q3(i) / (p3(i) + y + frac)
        END DO
        daw = (p3(10) + frac) / xx
      END IF

    END IF

    RETURN
  END FUNCTION daw715

END MODULE wofz_v2_mod
