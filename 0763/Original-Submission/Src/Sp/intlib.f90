!*** basicops.f
      SUBROUTINE ADD (A,B,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!***********************************************************************
!
! Called by --
!
!  any routine requiring interval addition.
!
!***********************************************************************
!
! Function --
!
!  This routine adds the interval A and the interval B.  It
!  simulates directed roundings with the routine RNDOUT;  the interval
!  result should contain the interval which would have been obtained
!  with exact interval arithmetic.  However, in general it will not
!  be the smallest possible machine-representable such containing
!  interval.  See the documentation in subroutine RNDOUT for more
!  detailed information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  A        is the first operand to the addition.
!           (INPUT)
!
!  B        is the second operand to the addition.
!           (INPUT)
!
!  RESULT   is the interval-arithmetic sum of A and B.
!           (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      LOGICAL RNDDWN, RNDUP
!
!***********************************************************************
!
! Internal variable descriptions --
!
! RNDDWN    is set to .TRUE. if RNDOUT has to round down, and is set
!           to .FALSE. otherwise.
!
! RNDUP     is set to .TRUE. if RNDOUT has to round up, and is set
!           to .FALSE. otherwise.
!
!***********************************************************************
!
! Common block declarations -- none
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines -- none
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
! RNDOUT
!
!***********************************************************************
!
! User-supplied functions and subroutines -- none
!
!***********************************************************************
!
! I/O functions --  none
!
!***********************************************************************
!
! Internal constant declarations -- none
!
!***********************************************************************
!
! Beginning of executable statements --
!
      RNDDWN = (A(1).NE.0D0).AND.(B(1).NE.0D0)
      RNDUP = (A(2).NE.0D0).AND.(B(2).NE.0D0)
!
      RESULT(1) = A(1) + B(1)
      RESULT(2) = A(2) + B(2)
!
      CALL RNDOUT(RESULT,RNDDWN,RNDUP)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE CANCEL (A,B,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   March 13, 1990
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!***********************************************************************
!
! Called by --
!
!  Any routine requiring an interval addend to be canceled (removed)
!  from a previously accumulated interval sum.
!
!***********************************************************************
!
! Function --
!
!  Given an interval B, and a previously accumulated interval sum A
!  for which B was an addend, this routine returns an interval which
!  contains, and is hopefully close to, the sum of the other addends
!  for A.  Directed roundings are simulated with the routine RNDOUT;
!  the interval result should contain the interval which would have
!  been obtained with exact interval arithmetic.  However, in general
!  it will not be the smallest possible machine-representable such
!  containing interval.  See the documentation in subroutine RNDOUT for
!  more detailed information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  A        is the first operand to the subtraction.
!           (INPUT)
!
!  B        is the second operand to the subtraction
!           (INPUT)
!
!  RESULT   is the interval-arithmetic value of A - B.
!           (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      LOGICAL RNDDWN, RNDUP
!
!***********************************************************************
!
! Internal variable descriptions --
!
!
! RNDDWN    is set to .TRUE. if RNDOUT has to round down, and is set
!           to .FALSE. otherwise.
!
! RNDUP     is set to .TRUE. if RNDOUT has to round up, and is set
!           to .FALSE. otherwise.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
! RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      RNDDWN = (B(1).NE.0D0)
      RNDUP = (B(2).NE.0D0)
!
      RESULT(1) = A(1) - B(1)
      RESULT(2) = A(2) - B(2)
!
      CALL RNDOUT(RESULT,RNDDWN,RNDUP)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IDIV (AA,BB,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   May 21, 1992
!
! Part of the interval arithmetic elementary function package.
!
!***********************************************************************
!
! Function --
!
!  This routine performs interval division.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION AA(2), BB(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --

      DOUBLE PRECISION C(2)
!
!***********************************************************************
!
!  Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ERRTST, MULT, RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 3
      IERR = 0
!
!  Do usual interval division if zero is not in the denominator --
!
      IF (BB(1).GT.ZERO(2)) THEN
         C(1) = ONE(1)/BB(2)
         C(2) = ONE(2)/BB(1)
         CALL RNDOUT(C,.TRUE.,.TRUE.)
         CALL MULT(AA,C,RESULT)
      ELSE IF (BB(2).LT.ZERO(1)) THEN
         C(1) = ONE(2)/BB(2)
         C(2) = ONE(1)/BB(1)
         CALL RNDOUT(C,.TRUE.,.TRUE.)
         CALL MULT(AA,C,RESULT)
      ELSE
         IERR = 6
         CALL ERRTST(BB)
         RESULT(1) = NEGINF
         RESULT(2) = POSINF
      END IF
!
      RETURN
!
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE MULT (A, B, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   November 17, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  MULT multiplies the interval A and the interval B.  It
!  puts the result into output parameter RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP, A1, B1
      DOUBLE PRECISION AA(2), BB(2)
!
!***********************************************************************
!
! Internal Constant Declarations --
!
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Pictures for cases.
!
      IF ((ZERO .LE. A(1)) .AND. (ZERO .LE. B(1))) THEN
!  Case 1   ---------------- 0 -----------------
!        A:                    [==========]
!        B:                         [===========]
          RESULT(1) = A(1) * B(1)
          RESULT(2) = A(2) * B(2)
      ELSE IF ((A(1) .LT. ZERO) .AND. (ZERO .LT. A(2))                  &
     &                          .AND. (ZERO .LE. B(1))) THEN
!  Case 2   ---------------- 0 -----------------
!        A:            [==================]
!        B:                         [===========]
          RESULT(1) = A(1) * B(2)
          RESULT(2) = A(2) * B(2)
      ELSE IF ((A(2) .LE. ZERO) .AND. (ZERO .LE. B(1))) THEN
!  Case 3   ---------------- 0 -----------------
!        A:    [==========]
!        B:                         [===========]
          B1 = B(1)
          RESULT(1) = A(1) * B(2)
          RESULT(2) = A(2) * B1
      ELSE IF ((ZERO .LE. A(1)) .AND. (B(1) .LT. ZERO)                  &
     &                          .AND. (ZERO .LT. B(2))) THEN
!  Case 4   ---------------- 0 -----------------
!        A:                     [==========]
!        B:               [===========]
          RESULT(1) = A(2) * B(1)
          RESULT(2) = A(2) * B(2)
      ELSE IF ((A(2) .LE. ZERO) .AND. (B(1) .LT. ZERO)                  &
     &                          .AND. (ZERO .LT. B(2))) THEN
!  Case 5   ---------------- 0 -----------------
!        A:    [==========]
!        B           [===========]
          A1 = A(1)
          B1 = B(1)
          RESULT(1) = A(1) * B(2)
          RESULT(2) = A1 * B1
      ELSE IF ((ZERO .LE. A(1)) .AND. (B(2) .LE. ZERO)) THEN
!  Case 6   ---------------- 0 -----------------
!        A:                    [==========]
!        B:  [===========]
          A1 = A(1)
          RESULT(1) = A(2) * B(1)
          RESULT(2) = A1 * B(2)
      ELSE IF ((A(1) .LT. ZERO) .AND. (ZERO .LT. A(2))                  &
     &                          .AND. (B(2) .LE. ZERO)) THEN
!  Case 7   ---------------- 0 -----------------
!        A:              [==========]
!        B:  [===========]
          A1 = A(1)
          B1 = B(1)
          RESULT(1) = A(2) * B(1)
          RESULT(2) = A1 * B1
      ELSE IF ((A(2) .LE. ZERO) .AND. (B(2) .LE. ZERO)) THEN
!  Case 8   ---------------- 0 -----------------
!        A:    [==========]
!        B:   [===========]
          A1 = A(1)
          B1 = B(1)
          RESULT(1) = A(2) * B(2)
          RESULT(2) = A1 * B1
      ELSE
!  Case 9   ---------------- 0 -----------------
!        A:              [==========]
!        B:              [===========]
!  Must check two cases.
          AA(1) = A(1)
          AA(2) = A(2)
          BB(1) = B(1)
          BB(2) = B(2)
          RESULT(1) = AA(1) * BB(2)
          TEMP = AA(2) * BB(1)
          IF (TEMP .LT. RESULT(1)) THEN
              RESULT(1) = TEMP
          ELSE
          END IF
              RESULT(2) = AA(1) * BB(1)
          TEMP = AA(2) * BB(2)
          IF (TEMP .GT. RESULT(2)) THEN
              RESULT(2) = TEMP
          ELSE
          END IF
      END IF
!
      CALL RNDOUT(RESULT,.TRUE.,.TRUE.)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE RNDOUT (X,RNDDWN,RNDUP)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
!   Modified August, 1991 and April, 1992 by Kaisheng Du and
!   R. Baker Kearfott to be more accurate near the underflow
!   threshold.
!
! Part of the generalized bisection package
! (interval arithmetic subpackage)
!
!***********************************************************************
!
! Function --
!
!  This routine is intended to simulate directed roundings in a
!  reasonably transportable way.  It is called for each elementary
!  operation involving intervals.  The endpoints of the result interval
!  are first computed with the machine's usual floating point
!  arithmetic.
!
!  If RNDDWN = .TRUE., then this routine decreases the left
!  endpoint of that approximate result by the absolute value of
!  that endpoint times a rigorous estimate for the maximum relative
!  error in an elementary operation.
!
!  If RNDUP = .TRUE., then this routine increases the right
!  endpoint of that approximate result by the absolute value of
!  that endpoint times a rigorous estimate for the maximum relative
!  error in an elementary operation.
!
!  For this routine to work properly, a machine-dependent parameter
!  must be installed in the routine SIMINI.  See the documentation in
!  that routine for details.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2)
      LOGICAL RNDDWN, RNDUP
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
! X         is the interval to be adjusted.
!           (I/O)
!
! RNDDWN    is set to .TRUE. if the left endpoint is to be adjusted, and
!           is set to .FALSE. otherwise.
!           (INPUT)
!
! RNDUP     is set to .TRUE. if the right endpoint is to be adjusted,
!           and is set to .FALSE. otherwise.
!           (INPUT)
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
!  This common block holds machine parameters which are set in
!  SIMINI and used here.
!
!  Variable descriptions
!
!  MXULP       (machine epsilon)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TTINY2      2 * (smallest representable positive machine number)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TOL0        TTINY2 / MXULP
!
      DOUBLE PRECISION TINY, TEST
      COMMON /MACH2/ TINY, TEST
!
!  See SIMINI for an explanation of the common block MACH2.
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IF (RNDDWN) THEN
         IF (X(1).GE.TEST) THEN
            X(1) = (1.D0 - MXULP ) * X(1)
         ELSE IF (X(1).LE.-TEST) THEN
            X(1) = (1D0 + MXULP ) * X(1)
         ELSE IF (X(1).LE.0.D0) THEN
            X(1) = -TEST
         ELSE
            X(1) = 0.D0
         END IF
      END IF
!
      IF (RNDUP) THEN
         IF (X(2).GE.TEST) THEN
            X(2) = (1.D0 + MXULP )* X(2)
         ELSE IF (X(2).LE.-TEST) THEN
            X(2) = (1.D0 - MXULP ) * X(2)
         ELSE IF(X(2).GE.0D0) THEN
            X(2) = TEST
         ELSE
            X(2) = 0.D0
         ENDIF
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE SCLADD (R,A,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!***********************************************************************
!
! Called by --
!
!  Any routine requiring addition of a point value to an interval.
!
!***********************************************************************
!
! FUNCTION --
!
!  This routine adds the interval A to the point R.  It simulates
!  directed roundings with the routine RNDOUT;  the interval
!  result should contain the interval which would have been obtained
!  with exact interval arithmetic.  However, in general it will not
!  be the smallest possible machine-representable such containing
!  interval.  See the documentation in subroutine RNDOUT for more
!  detailed information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION R, A(2), RESULT(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  R        is the point to be added to the interval.
!           (INPUT)
!
!  A        is the interval to be added to the point.
!           (INPUT)
!
!  RESULT   is the interval-arithmetic sum of R and B.
!           (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      LOGICAL RNDDWN, RNDUP
!
!***********************************************************************
!
! Internal variable descriptions --
!
! RNDDWN    is set to .TRUE. if RNDOUT has to round down, and is set
!           to .FALSE. otherwise.
!
! RNDUP     is set to .TRUE. if RNDOUT has to round up, and is set
!           to .FALSE. otherwise.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      RNDDWN = (R.NE.0D0).AND.(A(1).NE.0D0)
      RNDUP = (R.NE.0D0).AND.(A(2).NE.0D0)
!
      RESULT(1) = R + A(1)
      RESULT(2) = R + A(2)
!
      CALL RNDOUT(RESULT,RNDDWN,RNDUP)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE SCLMLT (R,A,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
! Modified by Manuel Novoa III on March 13, 1990 to clean the code
!   slightly, and to remove the need for MAX and MIN.
!
!***********************************************************************
!
! Called by --
!
!  Any routine requiring multiplication of an interval and a point
!  value.
!
!***********************************************************************
!
! Function --
!
!  This routine multiplies the interval A and the point R.  It
!  simulates directed roundings with the routine RNDOUT;  the interval
!  result should contain the interval which would have been obtained
!  with exact interval arithmetic.  However, in general it will not
!  be the smallest possible machine-representable such containing
!  interval.  See the documentation in subroutine RNDOUT for more
!  detailed information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION R, A(2), RESULT(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  R        is the point to be multiplied to the interval.
!           (INPUT)
!
!  A        is the interval to be multiplied to the point.
!           (INPUT)
!
!  RESULT   is the interval-arithmetic product R * B.
!           (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      LOGICAL RNDDWN, RNDUP
      DOUBLE PRECISION T1, T2
!
!***********************************************************************
!
! Internal variable descriptions --
!
! RNDDWN    is set to .TRUE. if RNDOUT is to round down, and is set
!           to .FALSE. otherwise.
!
! RNDUP     is set to .TRUE. if RNDOUT is to round up, and is set
!           to .FALSE. otherwise.
!
! T1 and T2 are temporary variables.
!
!***********************************************************************
!
! Common block declarations -- none
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --  none
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! User-supplied functions and subroutines -- none
!
!***********************************************************************
!
! I/O functions --  none
!
!***********************************************************************
!
! Internal constant declarations -- none
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IF ((R.EQ.0D0).OR.((A(1).EQ.0D0).AND.(A(2).EQ.0D0))) THEN
         RESULT(1) = 0D0
         RESULT(2) = 0D0
         RETURN
      END IF
!
      T1 = A(1)
      T2 = A(2)
      RNDDWN = .TRUE.
      RNDUP = .TRUE.
!
      IF (T1.EQ.0D0) THEN
         IF (R.LT.0D0) THEN
            RESULT(1) = R * T2
            RESULT(2) = 0D0
            RNDUP = .FALSE.
         ELSE
            RESULT(1) = 0D0
            RESULT(2) = R * T2
            RNDDWN = .FALSE.
         END IF
      ELSE IF (T2.EQ.0D0) THEN
         IF (R.LT.0D0) THEN
            RESULT(1) = 0D0
            RESULT(2) = R * T1
            RNDDWN = .FALSE.
         ELSE
            RESULT(1) = R * T1
            RESULT(2) = 0D0
            RNDUP = .FALSE.
         END IF
      ELSE IF (R.GT.0D0) THEN
         RESULT(1) = R * T1
         RESULT(2) = R * T2
      ELSE
         RESULT(1) = R * T2
         RESULT(2) = R * T1
      END IF
!
      CALL RNDOUT(RESULT,RNDDWN,RNDUP)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE SUB (A,B,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!***********************************************************************
!
! Called by --
!
!  Any routine requiring interval subtraction.
!
!***********************************************************************
!
! Function --
!
!  This routine subtracts the interval B from the interval A.  It
!  simulates directed roundings with the routine RNDOUT;  the interval
!  result should contain the interval which would have been obtained
!  with exact interval arithmetic.  However, in general it will not
!  be the smallest possible machine-representable such containing
!  interval.  See the documentation in subroutine RNDOUT for more
!  detailed information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  A        is the first operand to the subtraction.
!           (INPUT)
!
!  B        is the second operand to the subtraction
!           (INPUT)
!
!  RESULT   is the interval-arithmetic value of A - B.
!           (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TA1, TA2, TB1, TB2
      LOGICAL RNDDWN, RNDUP
!
!***********************************************************************
!
! Internal variable descriptions --
!
! TA1, TA2, TB1, and TB2 are temporaries.
!
! RNDDWN    is set to .TRUE. if RNDOUT has to round down, and is set
!           to .FALSE. otherwise.
!
! RNDUP     is set to .TRUE. if RNDOUT has to round up, and is set
!           to .FALSE. otherwise.
!
!***********************************************************************
!
! Common block declarations -- none
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines -- none
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
! RNDOUT
!
!***********************************************************************
!
! User-supplied functions and subroutines -- none
!
!***********************************************************************
!
! I/O functions --  none
!
!***********************************************************************
!
! Internal constant declarations -- none
!
!***********************************************************************
!
! Beginning of executable statements --
!
      TA1 = A(1)
      TA2 = A(2)
      TB1 = B(1)
      TB2 = B(2)
!
      RNDDWN = (TB2.NE.0D0)
      RNDUP = (TB1.NE.0D0)
!
      RESULT(1) = TA1 - TB2
      RESULT(2) = TA2 - TB1
!
      CALL RNDOUT(RESULT,RNDDWN,RNDUP)
!
      RETURN
      END

!*** utilfuns.f
      SUBROUTINE ICAP (A, B, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  DICAP computes the intersection of the two intervals A and B, and
!  places  the result in RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
!
!***********************************************************************
!
! Common block declarations --
!
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!   MIN, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!   ERRTST
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 13
      IERR = 0
!
      T(1) = MAX (A(1), B(1))
      T(2) = MIN (A(2), B(2))
      RESULT(1) = T(1)
      RESULT(2) = T(2)
      IF (RESULT(1).GT.RESULT(2)) THEN
         IERR=13
         CALL ERRTST(RESULT)
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      LOGICAL FUNCTION IDISJ (A, B)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  This function returns .TRUE. if the intervals A and B are disjoint,
!  and .FALSE. otherwise.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IDISJ = (A(2) .LT. B(1)) .OR. (B(2) .LT. A(1))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IHULL (A, B, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IHULL returns the convex-hull of the interval A and the interval B
!  in RESULT.  If one of the intervals is empty (lower bound is greater
!  than upper bound), then just the upper interval is returned.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IF ( A(1).GT.A(2) ) THEN
         IF (B(1).GT.B(2)) THEN
            T(1) = MAX(A(2),B(2))
            T(2) = MIN(A(1),B(2))
         ELSE
            T(1) = B(1)
            T(2) = B(2)
         END IF
      ELSE IF ( B(1).GT.B(2) ) THEN
         T(1) = A(1)
         T(2) = A(2)
      ELSE
         T(1) = MIN (A(1), B(1))
         T(2) = MAX (A(2), B(2))
         RESULT(1) = T(1)
         RESULT(2) = T(2)
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      LOGICAL FUNCTION IILEI (A, B)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IILEI sets its value to .TRUE. if interval A is in the closure of
!  interval B.  The value of IILEI is .FALSE. otherwise.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IILEI = (A(1) .GE. B(1)) .AND. (A(2) .LE. B(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      LOGICAL FUNCTION IILTI (A, B)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IILE sets its value to .TRUE. if interval A is in the interior of
!  interval B.  The value of IILE is .FALSE. otherwise.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), B(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IILTI = (A(1) .GT. B(1)) .AND. (A(2) .LT. B(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION IINF (A)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IINF returns the lower endpoint of the interval A.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IINF = A(1)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION IMID (X)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IMID returns a floating point approximation to the midpoint
!  of the interval A, using available floating point arithmetic.  The
!  value returned by this routine can be considered to DEFINE the
!  midpoint.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2)
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IMID = X(1) + (X(2) - X(1)) / TWO(1)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION IMIG (A)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IMIG returns the mignitude of the interval A. Since ABS is not
!  assumed to give an exact result, the result is rounded down.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2)
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  ABS, MIN
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IF ( ((A(1).GT.0D0) .AND. (A(2).GT.0D0))                          &
     & .OR.((A(1).LT.0D0) .AND. (A(2).LT.0D0)) ) THEN
         TMP(1) = ABS(A(2))
         CALL RNDOUT(TMP,.TRUE.,.FALSE.)
         TMP(2)=TMP(1)
         TMP(1) = ABS(A(1))
         CALL RNDOUT(TMP,.TRUE.,.FALSE.)
         IMIG = MIN (TMP(1), TMP(2))
      ELSE
         IMIG = 0D0
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE INEG (A, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  INEG performs interval unary negation.  The result is rounded out in
!  case the negatives of the endpoints are not representable.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), RESULT(2)
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      T(1) = A(1)
      T(2) = A(2)
      RESULT(1) = -T(2)
      RESULT(2) = -T(1)
      CALL RNDOUT (RESULT, .TRUE.,.TRUE.)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION INTABS (A)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  INTABS returns the maximum absolute value of the interval A.
!  Since ABS is not assumed to give an exact result, the result is
!  rounded up.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2)
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  ABS, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      TMP(2) = ABS(A(1))
      CALL RNDOUT(TMP,.FALSE.,.TRUE.)
      TMP(1)=TMP(2)
      TMP(2) = ABS(A(2))
      CALL RNDOUT(TMP,.FALSE.,.TRUE.)
!
      INTABS = MAX (TMP(1), TMP(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      LOGICAL FUNCTION IRLEI (A, B)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IRLEI sets its value to .TRUE. if double precision value A is in the
!  closure of the interval B. The value of IRLEI is set to .FALSE.
!  otherwise.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A, B(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IRLEI = (A .GE. B(1)) .AND. (A .LE. B(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      LOGICAL FUNCTION IRLTI (A, B)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IRLEI sets its value to .TRUE. if double precision value A is in the
!  interior of the interval B. The value of IRLEI is set to .FALSE.
!  otherwise.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A, B(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IRLTI = (A .GT. B(1)) .AND. (A .LT. B(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION ISUP (A)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IINF returns the upper endpoint of the interval A.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      ISUP = A(2)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IVL1 (R, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IVL1 constructs an interval RESULT from the double precision variable
!  R.  This is done by using simulated directed roundings with RNDOUT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION R, RESULT(2)
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T
!
!***********************************************************************
!
! Beginning of executable statements --
!
      T = R
      RESULT(1) = T
      RESULT(2) = T
      IF (T.NE.0D0) CALL RNDOUT (RESULT, .TRUE.,.TRUE.)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IVL2 (R1, R2, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!
!   Hong Jiang
!   Dept. of Math., Stat. and Computer Science
!   Marquette University
!   Milwaukee, WI 53233
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IVL2 constructs an interval RESULT, roughly equal to [R1,R2], from
!  the double precision variables R1 and R2.  This is done by using
!  simulated directed roundings with RNDOUT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION R1, R2, RESULT(2)
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T(2)
      LOGICAL RNDUP, RNDDWN
!
!***********************************************************************
!
! Common block declarations --
!
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 15
      IERR = 0
!
      T(1) = R1
      T(2) = R2
      RESULT(1) = T(1)
      RESULT(2) = T(2)
!
      IF (RESULT(2).LT.RESULT(1)) THEN
         IERR = 1
         CALL ERRTST(RESULT)
      END IF
!
      RNDDWN = T(1).NE.0D0
      RNDUP = T(2).NE.0D0
      CALL RNDOUT (RESULT, RNDDWN, RNDUP)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IVLABS (A, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IVLABS returns a rigorous bound on the range of the absolute value
!  of the interval A.  The intrinsic function ABS is assumed to be
!  accurate to the same accuracy as the elementary operations.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2), TA(2)
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  ABS
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      TA(1) = A(1)
      TA(2) = A(2)

      TMP(1) = ABS(TA(1))
      TMP(2) = ABS(TA(2))

      IF (TMP(1).LE.TMP(2)) THEN
         RESULT(1) = TMP(1)
         RESULT(2) = TMP(2)
      ELSE
         RESULT(1) = TMP(2)
         RESULT(2) = TMP(1)
      END IF

      CALL RNDOUT(RESULT,.TRUE.,.TRUE.)

      IF ( ( TA(1).LE.0 .AND. TA(2).GE.0 ) .OR. RESULT(1).LT. 0)        &
     &   RESULT(1) = 0D0

      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IVLI (A, RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   January 8, 1993
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IVLI places the contents of interval A in RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2), RESULT(2)
!
!***********************************************************************
!
! Beginning of executable statements --
!
      RESULT(1) = A(1)
      RESULT(2) = A(2)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      DOUBLE PRECISION FUNCTION IWID (A)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
! (Utility function)
!
!***********************************************************************
!
! Function --
!
!  IWID returns the width of the interval A, rounded up.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION A(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2)
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      TMP(2) = A(2)-A(1)
      CALL RNDOUT(TMP,.FALSE.,.TRUE.)
!
      IWID = TMP(2)
!
      RETURN
      END

!*** elemfuns.f
      SUBROUTINE IACOS (X,IRCCOS)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   July 12, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine computes an interval enclosure for the arccos over
!  the interval X and returns the result in the interval IRCCOS.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), IRCCOS(2)
!
!***********************************************************************
!
! Common blocks --
!
        DOUBLE PRECISION MXULP, TTINY2, TOL0
        COMMON /MACH1/MXULP, TTINY2, TOL0
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ERRTST, IASIN, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 8
      IERR = 0
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
      IF (X(2).GT.ONE(1)) THEN
         IERR = 9
         CALL ERRTST(X)
         IRCCOS(1) = NEGINF
         IRCCOS(2) = POSINF
         RETURN
      ELSE IF (X(1).LT.-ONE(1)+MXULP) THEN
         IERR = 10
         CALL ERRTST(X)
         IRCCOS(1) = NEGINF
         IRCCOS(2) = POSINF
         RETURN
      END IF
!
      CALL IASIN(X, IRCCOS)
!
      CALL SUB(PI2,IRCCOS,IRCCOS)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IACOT(X,IRCCOT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   June 16, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!   This program finds bounds on the value of ARCCOT(X) as X
!   ranges of the interval XX, and places the resulting interval
!   in IRCCOT.
!
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), IRCCOT(2)
!
!***********************************************************************
!
! Common blocks--
!
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  IATAN, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!

      CALL IATAN(X, IRCCOT)
      CALL SUB(PI2,IRCCOT,IRCCOT)

      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IASIN(XX,IRCSIN)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   July 7, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine computes an interval enclosure for the arcsin over
!  the interval XX and returns the result in the interval IRCSIN.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XX(2), IRCSIN(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION LB, UB
      DOUBLE PRECISION X(2), RSLT(2), TMP(2)
      LOGICAL OVER, NEGATV, FLIP
!
!***********************************************************************
!
! Common blocks --
!
        DOUBLE PRECISION MXULP, TTINY2, TOL0
        COMMON /MACH1/MXULP, TTINY2, TOL0
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --  none
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ASNSER, ERRTST, ISQRT, MULT, RNDOUT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 9
      IERR = 0
!
      IF (XX(2).LT.XX(1)) THEN
         IERR = 1
         CALL ERRTST(XX)
      END IF
!
      IF (XX(2).GT.ONE(1)) THEN
         IERR = 9
         CALL ERRTST(XX)
         RETURN
      ELSE IF (XX(1).LT.-ONE(1)+MXULP) THEN
         IERR = 10
         CALL ERRTST(XX)
         RETURN
      END IF
!
      LB = XX(1)
      UB = XX(2)
!
! For the lower end point ---
!----------------------------------------------------------------------
!
!  If LB < 0, convert it to positive --

      IF (LB .LT. ZERO(1)) THEN
         OVER = .TRUE.
         NEGATV = .TRUE.
         TMP(2) = -LB
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         LB = TMP(2)
      ELSE
         OVER = .FALSE.
         NEGATV =.FALSE.
      ENDIF
!
!  Transform LB to a small interval X such that 0 <= X <= 0.5 --
!
      X(1) = LB
      X(2) = LB
      IF (NEGATV) CALL RNDOUT(X,.TRUE.,.TRUE.)
!
      IF (LB .LE. OD2F(2) ) THEN
         FLIP = .FALSE.
      ELSE
         FLIP = .TRUE.
!
!   X <-- SQRT( (1 - X)/2 )
!
         CALL SUB(ONE,X,X)
         CALL MULT(OD2F,X,X)
         IF(X(1).LT.0D0) THEN
            X(1)=0D0
         END IF
         CALL ISQRT(X,X)
      ENDIF
!
!  Find a bound on the arcsin at the lower end point --
!
        CALL ASNSER (X, OVER, RSLT)
!
!  Undo the argument transformation, if it was applied --
!
      IF (FLIP) THEN
         CALL MULT(TWO,RSLT,RSLT)
         CALL SUB(PI2,RSLT,RSLT)
      END IF

      IF (NEGATV) THEN
         IRCSIN(1) = -RSLT(2)
         CALL RNDOUT(IRCSIN,.TRUE.,.FALSE.)
      ELSE
         IRCSIN(1) = RSLT(1)
      ENDIF
!
! For the upper end point ---
!----------------------------------------------------------------------
!
      IF (UB .LT. 0D0) THEN
         OVER =.FALSE.
         NEGATV = .TRUE.
         TMP(2) = -UB
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         UB = TMP(2)
      ELSE
         OVER = .TRUE.
         NEGATV = .FALSE.
      ENDIF
!
      X(1) = UB
      X(2) = UB
      IF (NEGATV) CALL RNDOUT(X,.TRUE.,.TRUE.)
!
      IF (UB .LE. OD2F(2) ) THEN
         FLIP = .FALSE.
      ELSE
         FLIP = .TRUE.
!
!   X <-- SQRT( (1 - X)/2 )
!
         CALL SUB(ONE,X,X)
         CALL MULT(OD2F,X,X)
         IF(X(1).LT.0D0) THEN
            X(1)=0D0
         END IF
         CALL ISQRT(X,X)
      ENDIF
!
      CALL ASNSER(X, OVER, RSLT)
!
      IF (FLIP) THEN
         CALL MULT(TWO,RSLT,RSLT)
         CALL SUB(PI2,RSLT,RSLT)
      END IF

      IF (NEGATV) THEN
         IRCSIN(2) = -RSLT(1)
         CALL RNDOUT(IRCSIN,.FALSE.,.TRUE.)
      ELSE
         IRCSIN(2) = RSLT(2)
      ENDIF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
!
      SUBROUTINE ASNSER (X, OVER, RSLT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   July 7, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine computes an interval enclosure for the arcsin over
!  the interval X and returns the result in the interval RSLT.  It is
!  assumed that X is nearly a point, between 0 and .5.  This routine is
!  normally called from IASIN.  The argument OVER indicates whether
!  the upper end point or lower end point is important.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2)
      LOGICAL OVER
      DOUBLE PRECISION RSLT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP1(2), TEMP2(2), TEMP3(2), SUM(2)
      DOUBLE PRECISION T
      INTEGER I
      DOUBLE PRECISION D2IM1(2), D2I(2), D2IP1(2)
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/MXULP, TTINY2, TOL0
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  ABS, DBLE, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ADD, ERRTST, IDIV, MULT, RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 10
      IERR = 0
!
      SUM(1) = X(1)
      SUM(2) = X(2)
!
      TEMP1(1) = X(1) * X(1)
      TEMP1(2) = X(2) * X(2)
      CALL RNDOUT(TEMP1, .TRUE., .TRUE.)
!
      TEMP2(1) = X(1)
      TEMP2(2) = X(2)
!
      DO 10 I = 1,100
!
         D2IM1(1) = DBLE(2*I-1)
         D2IM1(2) = D2IM1(1)
         CALL RNDOUT(D2IM1,.TRUE.,.TRUE.)
         D2I(1) = DBLE(2*I)
         D2I(2) = D2I(1)
         CALL RNDOUT(D2I,.TRUE.,.TRUE.)
         D2IP1(1) = DBLE(2*I+1)
         D2IP1(2) = D2IP1(1)
         CALL RNDOUT(D2IP1,.TRUE.,.TRUE.)
!
!  TEMP2 <-- (2i-1)TEMP1 TEMP2 / (2i)
!
         CALL MULT(D2IM1,TEMP2,TEMP2)
         CALL MULT(TEMP2,TEMP1,TEMP2)
         CALL IDIV(TEMP2,D2I,TEMP2)
!
         CALL IDIV(TEMP2,D2IP1,TEMP3)
!
         SUM(1) = SUM(1) + TEMP3(1)
         SUM(2) = SUM(2) + TEMP3(2)
         CALL RNDOUT(SUM,.TRUE.,.TRUE.)
         T = MAX( ABS(SUM(1)), TOL0 )
         IF (TEMP3(2)/T .LT. MXULP) GOTO 20
   10 END DO
         IERR=11
         CALL ERRTST(TEMP3)
         RETURN
   20 CONTINUE
!
      IF (OVER) THEN
         CALL MULT(TWO,TEMP3,TEMP3)
         CALL ADD(SUM,TEMP3,RSLT)
      ELSE
         RSLT(1) = SUM(1)
         RSLT(2) = SUM(2)
      ENDIF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IATAN(XX,IRCTAN)

!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   June 16, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!   This program finds bounds on the value of ARCTAN(X) as X
!   ranges over the interval XX, and places the resulting interval
!   in IRCTAN.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XX(2), IRCTAN(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION LB, UB, X(2), RSLT(2), TMP(2)
      INTEGER CODE
      LOGICAL EVEN, NEGATV, RLB, RUB
!
!***********************************************************************
!
! Common blocks--
!
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/MXULP, TTINY2, TOL0
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ADD, ATNRED, ATNSER, RNDOUT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 11
      IERR = 0
!
      IF (XX(2).LT.XX(1)) THEN
         IERR = 1
         CALL ERRTST(XX)
      END IF
!
      LB = XX(1)
      UB = XX(2)
!
!  Find the logical value EVEN, used in subroutine ATNSER to determine
!  the number of terms to be computed in the series --
!
      IF (LB .GT. 1.0D0) THEN
         EVEN = .TRUE.
      ELSE IF (LB .GE. 0D0) THEN
         EVEN = .FALSE.
      ELSE IF (LB .GT. -1D0) THEN
         EVEN = .FALSE.
      ELSE
         EVEN = .TRUE.
      ENDIF
!
!  If LB < 0, convert it to positive and update even --

      IF (LB .LT. 0D0) THEN
         EVEN = .NOT. EVEN
         NEGATV = .TRUE.
         TMP(2) = -LB
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         LB = TMP(2)
      ELSE
         NEGATV =.FALSE.
      ENDIF
!
!  Transform LB to an interval X such that 0 <= X < 1/sqrt(3).  Here, X
!  is a small interval containing the lower bound of the the interval
!  [LB,UB] after transformation --
!
      CALL ATNRED(LB, CODE, X)
!
!  Find the value of the arctangent for the interval X --
!
      CALL ATNSER(X, EVEN, RSLT)
!
!  Undo the transformation and store the result --
!
      IF (CODE .EQ. 11) THEN
      ELSE IF (CODE .EQ. 12 ) THEN
         CALL ADD(PI6,RSLT,RSLT)
      ELSE IF (CODE .EQ. 21 ) THEN
         CALL SUB(PI2,RSLT,RSLT)
      ELSE IF (CODE .EQ. 22 ) THEN
         CALL SUB(PI3,RSLT,RSLT)
      END IF
!
      IF (NEGATV) THEN
         IRCTAN(1) = -RSLT(2)
         RLB = .TRUE.
      ELSE
         IRCTAN(1) = RSLT(1)
         RLB = .FALSE.
      ENDIF
!
!  Similarly obtain an enclosure for the arctangent at the upper
!  end point --
!
      IF (UB .GT. 1.0D0) THEN
         EVEN = .FALSE.
      ELSE IF (UB .GE. 0D0) THEN
         EVEN = .TRUE.
      ELSE IF (UB .GT. -1D0) THEN
         EVEN = .TRUE.
      ELSE
         EVEN = .FALSE.
      ENDIF
!
      IF (UB .LT. 0D0) THEN
         EVEN = .NOT. EVEN
         NEGATV = .TRUE.
         TMP(2) = -UB
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         UB = TMP(2)
      ELSE
         NEGATV = .FALSE.
      ENDIF
!
      CALL ATNRED(UB, CODE, X)
      CALL ATNSER(X, EVEN, RSLT)
!
      IF (CODE .EQ. 11) THEN
      ELSE IF (CODE .EQ. 12 ) THEN
         CALL ADD(PI6,RSLT,RSLT)
      ELSE IF (CODE .EQ. 21 ) THEN
         CALL SUB(PI2,RSLT,RSLT)
      ELSE IF (CODE .EQ. 22 ) THEN
         CALL SUB(PI3,RSLT,RSLT)
      END IF
!
      IF (NEGATV) THEN
         IRCTAN(2) = -RSLT(1)
         RUB = .TRUE.
      ELSE
         IRCTAN(2) = RSLT(2)
         RUB = .FALSE.
      ENDIF
!
      CALL RNDOUT(IRCTAN,RLB,RUB)
!
      RETURN
      END
!**********************************************************************
!**********************************************************************
      SUBROUTINE ATNRED(PP,CODE,X)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   June 16, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This is an auxiliary routine for IATAN.  It performs an argument
!  reduction in preparation for computing a series approximation to the
!  arctangent function.
!
!  The subroutine transforms a given nonnegative number P to a
!  nonnegative number less than the reciprocal of the square root of 3,
!  then converts it to an interval.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION PP, X(2)
      INTEGER CODE
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
! PP     A floating point representation of an endpoint, to be reduced.
!        (INPUT)
!
! X      A small interval enclosure for the reduced quantity
!        corresponding to PP.
!        (OUTPUT)
!
! CODE   A flag indicating how the argument was reduced, as follows:
!
!          CODE = 11      no transformation
!          CODE = 12      1 / sqrt(3) < P < = 1
!          CODE = 21      P > 1 and 0 < 1/P < 1 / sqrt(3)
!          CODE = 22      P > 1 and 1/P >= 1 / sqrt(3)
!
!        (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION P(2), TMP(2)
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ADD, IDIV, MULT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
      P(1) = PP
      P(2) = PP
!
      IF (P(1) .LE. ODSQT3(1)) THEN
         CODE = 11
      ELSE IF (P(1) .LE. ONE(2)) THEN
!
!  P <-- (sqrt(3)P - 1)/(sqrt(3)+P):
!
         CALL MULT(P,SQT3,TMP)
         CALL SUB(TMP,ONE,TMP)
         CALL ADD(SQT3,P,P)
         CALL IDIV(TMP,P,P)
         CODE = 12
      ELSE
!
!  P <-- 1/P:
!
         CALL IDIV(ONE,P,P)
         IF (P(1) .LE. ODSQT3(1)) THEN
            CODE = 21
         ELSE
            CALL MULT(P,SQT3,TMP)
            CALL SUB(TMP,ONE,TMP)
            CALL ADD(SQT3,P,P)
            CALL IDIV(TMP,P,P)
            CODE = 22
         ENDIF
      ENDIF
!
      X(1) = P(1)
      X(2) = P(2)
!
      RETURN
      END
!**********************************************************************
!**********************************************************************
      SUBROUTINE ATNSER (X,EVEN,RSLT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   By Chenyi Hu
!       and
!   Abdulhamid Awad
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!       and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   June 16, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This subroutine computes an enclosure for arctan X and puts the
!  result in RSLT.  It is assumed that X is an interval of small width
!  between 0 and 1 / sqrt(3).  The argument EVEN is used to determine
!  whether an even or odd number of terms should be taken.  This routine
!  is normally called from IATAN.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RSLT(2)
      LOGICAL EVEN
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP1(2), TEMP2(2), TEMP3(2), SUM(2), XS(2)
      DOUBLE PRECISION DI(2), T
      INTEGER I, SIGN
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  DBLE
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ADD, IDIV, MULT, POWER, RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 12
      IERR = 0
!
      IF (X(1).LT.0D0) THEN
         IERR=12
         CALL ERRTST(X)
         RETURN
      END IF
!
      CALL POWER(X,2,XS)
!
      IF ((X(1) .EQ. 0D0) .AND. (X(2) .EQ. 0D0)) THEN
         RSLT(1) = ZERO(1)
         RSLT(2) = ZERO(2)
         RETURN
      ENDIF
!
      SUM(1) = X(1)
      SUM(2) = X(2)
!
      SIGN = -1
      CALL POWER(X,3,TEMP1)
!
      DO 10 I = 1,10000
         SIGN = -SIGN
         CALL MULT(TWO,TEMP1,TEMP2)
         DI(1) = DBLE (4*I*I - 1)
         DI(2) = DI(1)
         CALL RNDOUT(DI,.TRUE.,.TRUE.)
         CALL IDIV(TEMP2,DI,TEMP2)
!
         IF (SIGN .GT. 0) THEN
            SUM(1) = SUM(1) + TEMP2(1)
            SUM(2) = SUM(2) + TEMP2(2)
         ELSE
            SUM(1) = SUM(1) - TEMP2(2)
            SUM(2) = SUM(2) - TEMP2(1)
         ENDIF
         CALL RNDOUT(SUM,.TRUE.,.TRUE.)
!
         TEMP1(1) = TEMP1(1) * XS(1)
         TEMP1(2) = TEMP1(2) * XS(2)
         CALL RNDOUT(TEMP1, .TRUE., .TRUE.)
!
         T = MAX( ABS(SUM(1)),2D0*TOL0 )
         IF (TEMP1(2)/T .GT. MXULP) GOTO 10
         IF (EVEN .AND. (SIGN .LT. 0)) GOTO 10
         IF (EVEN .OR. (SIGN .LE. 0)) GOTO 20
   10 END DO
         IERR=11
         CALL ERRTST(TEMP2)
         RETURN
   20 CONTINUE
!
      CALL ADD(XS,ONE,TEMP3)
      CALL IDIV(SUM,TEMP3,RSLT)
!
        RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ICOS(X,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine returns the interval value of the cosine function
!  (evaluated over the interval X) in the interval variable
!  RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION V1(2), V2(2)
      DOUBLE PRECISION VAL1(2),VAL2(2)
      DOUBLE PRECISION TMP(2), T1, T2
!
      INTEGER N1, N2
!
      LOGICAL EVEN
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!     ABS, ERRTST, IDINT, MAX, MIN, MOD
!
!***********************************************************************
!
! Package-supplied functions and subroutines -- none
!
!     EXTERNAL ERRTST, MULT, RCOS, RNDOUT
!
!***********************************************************************
!
! User-supplied functions and subroutines -- none
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 2
      IERR = 0
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
!  Compute the value at the left endpoint using interval arithmetic --
!
      V1(1) = X(1)
      V1(2) = X(1)
!
!  If X1/PI is out of range, return the interval [-1,1] --
!
      CALL MULT(V1,A,TMP)
      T1 =ABS(TMP(1))
      IF(T1.GE.MAXX) THEN
         IERR = 5
         CALL ERRTST(X)
         RESULT(1) = -1D0
         RESULT(2) =  1D0
         RETURN
      ENDIF
!
!  Compute the value at the right endpoint using interval arithmetic --
!
      V2(1) = X(2)
      V2(2) = X(2)
      CALL RNDOUT(V2,.TRUE.,.TRUE.)
!
!  If X2/pi is out of range, return the interval [-1,1] --
!
      CALL MULT(V2,A,TMP)
      T2 = ABS(TMP(2))
      IF(T2.GE.MAXX) THEN
         IERR = 5
         CALL ERRTST(X)
         RESULT(1) = -1D0
         RESULT(2) =  1D0
         RETURN
      ENDIF
!
!  If both X1/pi and X2/pi are within range, calculate the function
!  values at the left and right endpoints using interval arithmetic --
!
      CALL RCOS(V1,VAL1)
!
      CALL RCOS(V2,VAL2)
!
!  Calculate the number of half-periods from zero for the left
!  and right endpoints in order to normalize to [-pi/2,pi/2] --
!
      N1 = IDINT(T1)
      N2 = IDINT(T2)
      IF(V1(1).LT.0D0) N1 = -N1-1
      IF(V2(2).LT.0D0) N2 = -N2-1
!
!   In even half periods, the function is decreasing --
!
      EVEN = MOD(N1,2).EQ.0
!
!   If X1 and X2 are in the same half period, then the lower bound
!   and upper bound on the range occur at the endpoints --
!
      IF(N1.EQ.N2) THEN
         RESULT(1) = MIN(VAL1(1),VAL2(1))
         RESULT(2) = MAX(VAL1(2),VAL2(2))
         RETURN
      ENDIF
!
!   Consider the case X2 is in the half period adjacent to X1 --
!
      IF(N2.EQ.N1+1) THEN
!
!    If X1 is in the increasing half period, then the upper bound should
!    equal 1;  otherwise, the lower bound should be -1 --
!
         IF(EVEN) THEN
            RESULT(1) = -1.D0
            RESULT(2) = MAX(VAL1(2),VAL2(2))
         ELSE
           RESULT(1) = MIN(VAL1(1),VAL2(1))
           RESULT(2) = 1.D0
         ENDIF
!
!   If X1 and X2 are not in adjecent half periods, then the lower
!   bound should be -1 and the upper bound should be 1 --
!
      ELSE
         RESULT(1) = -1.D0
         RESULT(2) =  1.D0
      ENDIF
!
      CALL RNDOUT(RESULT,.TRUE.,.TRUE.)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IEXP(X,RESULT)
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This routine returns the interval exp(X) in RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP(2), TVAL(2), T
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions --
!
!     ABS
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
      EXTERNAL ERRTST, REXP
      DOUBLE PRECISION D1MACH
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 1
      IERR = 0
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
      IF (X(1).EQ.ZERO(1)) THEN
         RESULT(1) =ONE(1)
      ELSE
         T = ABS(X(1))
         IF (T.GT.MXLGM1) THEN
            IF(X(1).GT.ZERO(1))THEN
               IERR = 3
            ELSE
               IERR = 2
            END IF
            CALL ERRTST(X)
         END IF
         IF (IERR.EQ.2) THEN
            RESULT(1) = ZERO(1)
         ELSE IF (IERR.EQ.3) THEN
            RESULT(1) = POSINF
            RESULT(2) = POSINF
            RETURN
         ELSE
            TEMP(1) = X(1)
            TEMP(2) = X(1)
            CALL REXP(TEMP,TVAL)
            RESULT(1) = TVAL(1)
         END IF
      END IF
!
      IF (X(2).EQ.ZERO(2)) THEN
         RESULT(2) =ONE(2)
      ELSE
         IF (X(2).GT.MXLGM1) THEN
            IERR = 4
            CALL ERRTST(X)
         END IF
         IF (IERR.EQ.4) THEN
            RESULT(2) = D1MACH(2)
            RETURN
         ELSE
            TEMP(1) = X(2)
            TEMP(2) = X(2)
            CALL REXP(TEMP,TVAL)
            RESULT(2) = TVAL(2)
         END IF
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE IIPOWR(XX, YY, RESULT)
!
!  Kaisheng Du and R. B. Kearfott, Summer, 1991
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!**********************************************************************
!
! Function --
!
!  This routine returns the interval X**Y, provided that X is a
!  positive interval.
!  RESULT.
!
!**********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XX(2), YY(2), RESULT(2)
!
!**********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION X(2), Y(2)
      LOGICAL LBZERO

!**********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!   ABS, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!   ERRTST, IEXP, ILOG, MULT
!
!**********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 14
      IERR = 0
!
      X(1) = XX(1)
      X(2) = XX(2)
      Y(1) = YY(1)
      Y(2) = YY(2)
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
      LBZERO = .FALSE.
      IF (X(1).EQ.0D0) THEN
         IF (Y(1).LT.0D0) THEN
            IF(Y(2).GE.0D0) THEN
               IERR = 14
               RESULT(1) = NEGINF
               RESULT(2) = POSINF
               RETURN
            END IF
         END IF
         LBZERO = .TRUE.
         X(1) = X(2)
      END IF

      IF (X(1).LE.ZERO(2)) THEN
         IERR = 7
         CALL ERRTST(X)
         RETURN
      END IF
!
      CALL ILOG(X,RESULT)
!
      CALL MULT(RESULT,Y,RESULT)
!
      IF (ABS(RESULT(1)).GT.MXLGM1) THEN
         IF(RESULT(1).GT.ZERO(1)) THEN
            IERR = 3
         ELSE
            IERR = 2
         END IF
         CALL ERRTST(RESULT)
      END IF
      IF (RESULT(2).GT.MXLGM1) THEN
         IERR = 4
         CALL ERRTST(RESULT)
      END IF
!
      CALL IEXP(RESULT,RESULT)
      IF(LBZERO) THEN
         IF (Y(1).LT.0D0) THEN
            RESULT(2) = POSINF
            IERR = 15
         ELSE
            RESULT(1) = 0D0
         END IF
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ILOG(X,RESULT)
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This routine returns the interval log(X) in RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP(2), TVAL(2)
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!     ERRTST, RLOG
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 4
      IERR = 0
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
      IF (X(1).LE.ZERO(2)) THEN
         IERR = 7
         CALL ERRTST(X)
         RETURN
      END IF
!
      TEMP(1) = X(1)
      TEMP(2) = X(1)
      CALL RLOG(TEMP,TVAL)
      RESULT(1) = TVAL(1)
!
      TEMP(1) = X(2)
      TEMP(2) = X(2)
      CALL RLOG(TEMP,TVAL)
      RESULT(2) = TVAL(2)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ISIN(X,RESULT)
!
!  Kaisheng Du and R. B. Kearfott, Summer, 1991
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!**********************************************************************
!
! Function --
!
!  This routine returns the interval value of the sine function
!  (evaluated over the interval X) in the interval variable
!  RESULT.
!
!**********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RESULT(2)
!
!**********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Fortran supplied functions and subroutines --
!
!      ABS, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!      ICOS, MULT, POWER, SUB
!
!**********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP1(2), TMP2(2), TMP3(2)
      LOGICAL SMALL1, SMALL2
!
!**********************************************************************
!
!
! Beginning of executable statements --
!
!  Assure good relative accuracy near zero by the case separately;
!  otherwise, shift the argument and use the interval cosine function --
!
      SMALL1=.FALSE.
      SMALL2=.FALSE.
      IF ( MAX( ABS(X(1)), ABS(X(2)) ) .LT. PI2(1) ) THEN
         IF(ABS(X(1)).LT.CBTEP) THEN
            SMALL1 = .TRUE.
            TMP1(1) = X(1)
            TMP1(2) = X(1)
            IF (X(1).GT.ZERO(2)) THEN
               CALL POWER(TMP1,3,TMP2)
               CALL MULT(TMP2,OD3F,TMP2)
               CALL SUB(TMP1,TMP2,TMP1)
            END IF
         END IF
         IF(ABS(X(2)).LT.CBTEP) THEN
            SMALL2 = .TRUE.
            TMP3(1) = X(2)
            TMP3(2) = X(2)
            IF (X(2).LT.ZERO(1)) THEN
               CALL POWER(TMP3,3,TMP2)
               CALL MULT(TMP2,OD3F,TMP2)
               CALL SUB(TMP3,TMP2,TMP3)
            END IF
         END IF
      END IF
!
      IF (SMALL1 .AND. SMALL2) THEN
         RESULT(1) = TMP1(1)
         RESULT(2) = TMP3(2)
         RETURN
      ELSE
         CALL SUB(X,PI2,RESULT)
         CALL ICOS(RESULT,RESULT)
      END IF
!
      IF (SMALL1) RESULT(1) = MAX(TMP1(1),RESULT(1))
      IF (SMALL2) RESULT(2) = MIN(TMP3(2),RESULT(2))
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ISINH(XX,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   Chenyi Hu
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   October 20, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This routine returns an interval inclusion for the hyperbolic sine
!  over the interval XX in RESULT.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XX(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION AA, BB, X(2), RSLT(2)
      DOUBLE PRECISION TMP(2), TMP2(2)
      INTEGER CODE
      LOGICAL OVER, NEGATV
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/MXULP, TTINY2, TOL0
!
!***********************************************************************
!
! Package-supplied functions and subroutines --

!  ERRTST, IDIV, IEXP, ISHSER, ISNRED, ISNVAL, MULT, RNDOUT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 6
      IERR = 0
!
      IF (XX(2).LT.XX(1)) THEN
         IERR = 1
         CALL ERRTST(XX)
      END IF
!
      AA = XX(1)
      BB = XX(2)
!
!-------------------------------------------------------------------
!  Compute the value of the hyperbolic sine at the lower endpoint A.
!-------------------------------------------------------------------
!
      IF (AA .LT. ZERO(1)) THEN
         OVER = .TRUE.
         NEGATV = .TRUE.
         TMP(2) = -AA
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         AA = TMP(2)
      ELSE
         OVER = .FALSE.
         NEGATV =.FALSE.
      ENDIF
!
!  Transform AA to an interval X such that 0 <= X <= 1 --
!
      CALL ISNRED(AA, CODE, X)
!
      IF (CODE .EQ. 20) THEN
!
!  Use the exponential function for arguments larger than 27 --
!
         CALL IEXP(X, TMP)
         CALL IDIV(ONE,TMP,TMP2)
         CALL SUB(TMP,TMP2,TMP)
         CALL MULT(OD2F,TMP,RSLT)
      ELSE
!
!  Evaluate the series X + X^3/3! + X^5/5! + ... for 0 <= X <= 1 --
!
         CALL ISHSER(X, OVER, RSLT)
         CALL ISNVAL(CODE, RSLT)
       END IF
!
!  Undo the possible argument reduction from ISNRED --
!
      IF (NEGATV) THEN
         RESULT(1) = - RSLT(2)
         CALL RNDOUT(RESULT,.TRUE.,.FALSE.)
      ELSE
         RESULT(1) = RSLT(1)
      ENDIF
!
!-------------------------------------------------------------------
!  Compute the value of the hyperbolic sine at the upper endpoint B
!  (completely analogous to computation for the lower endpoint A).
!-------------------------------------------------------------------
!
      IF (BB .LT. ZERO(1)) THEN
         OVER =.FALSE.
         NEGATV = .TRUE.
         TMP(2) = -BB
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         BB = TMP(2)
      ELSE
         OVER = .TRUE.
         NEGATV = .FALSE.
      ENDIF
!
      CALL ISNRED(BB, CODE, X)
!
      IF (CODE .EQ. 20) THEN
         CALL IEXP(X, TMP)
         CALL IDIV(ONE,TMP,TMP2)
         CALL SUB(TMP,TMP2,TMP)
         CALL MULT(OD2F,TMP,RSLT)
      ELSE
         CALL ISHSER(X, OVER, RSLT)
         CALL ISNVAL(CODE, RSLT)
      END IF
!
      IF (NEGATV) THEN
         RESULT(2) = - RSLT(1)
         CALL RNDOUT(RESULT,.FALSE.,.TRUE.)
      ELSE
         RESULT(2) = RSLT(2)
      ENDIF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ISNRED(P,CODE,X)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   Chenyi Hu
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   October 20, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This subroutine transforms the nonnegative number P to a nonnegative
!  number which less than 1, and then converts the result to an
!  interval.  The following code numbers are set to record how the
!  transformation was done:
!
!     CODE = 11      0 < P < =  1
!     CODE = 12      1 < P < =  3
!     CODE = 13      3 < P < =  9
!     CODE = 14      9 < P < = 27
!     CODE = 20          P >   27
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION P
      INTEGER CODE
      DOUBLE PRECISION X(2)
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Beginning of executable statements --
!
      X(1) = P
      X(2) = P
!
      IF (P .LE. ONE(1) ) THEN
         CODE = 11
         GO TO 10
      ELSE IF (P .LE. THREE(1)) THEN
         CALL MULT(X,THIRD,X)
         CODE = 12
         GO TO 10
      ELSE IF (P .LE. NINE(1)) THEN
         CALL MULT(X,NINTH,X)
         CODE = 13
         GO TO 10
      ELSE IF (P .LE. TWOT7(1)) THEN
         CALL MULT(X,TT7TH,X)
         CODE = 14
         GO TO 10
      ELSE
         CODE = 20
         GO TO 10
      ENDIF
!
   10 CALL RNDOUT(X, .TRUE., .TRUE.)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
        SUBROUTINE ISHSER (X, OVER, RSLT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   Chenyi Hu
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!         and
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   October 20, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!   This subroutine uses a series to bound the value of the hyperbolic
!   sine function over the nonnegative interval X of small width.  The
!   result is returned in RSLT.  The logical variable OVER indicates
!   whether negation was applied during the argument reduction.  This
!   routine is not meant to be called by the user.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2)
      LOGICAL OVER
      DOUBLE PRECISION RSLT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP1(2), TEMP2(2), SUM(2), ERR(2), T
      INTEGER I
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!  ABS, DBLE, FLOAT, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  MULT, POWER, RNDOUT
!
!***********************************************************************
!
! Beginning of executable statements --
!
      SUM(1) = X(1)
      SUM(2) = X(2)

      I = 1
!
!  Use temp1 to hold X^2 --
!
      CALL POWER(X,2,TEMP1)

      TEMP2(1) = X(1)
      TEMP2(2) = X(2)
!
   10 CONTINUE
      I = I + 2
         TEMP2(1) = TEMP1(1) * TEMP2(1)
         TEMP2(2) = TEMP1(2) * TEMP2(2)
         CALL RNDOUT(TEMP2, .TRUE., .TRUE.)
!
         TEMP2(1) = TEMP2(1)/DBLE(FLOAT(I))
         TEMP2(2) = TEMP2(2)/DBLE(FLOAT(I))
         CALL RNDOUT(TEMP2, .TRUE., .TRUE.)
!
         TEMP2(1) = TEMP2(1)/DBLE(FLOAT(I-1))
         TEMP2(2) = TEMP2(2)/DBLE(FLOAT(I-1))
         CALL RNDOUT(TEMP2, .TRUE., .TRUE.)
!
         SUM(1) = SUM(1) + TEMP2(1)
         SUM(2) = SUM(2) + TEMP2(2)
         CALL RNDOUT(SUM,.TRUE.,.TRUE.)
!
         CALL MULT(TWO,TEMP2,ERR)
         T = MAX(ABS(SUM(1)),ABS(SUM(2)),2D0*TOL0)
         IF ((ERR(1)/T) .GT. MXULP) GOTO 10
!
      IF (OVER) THEN
         RSLT(1) = SUM(1) + ERR(1)
         RSLT(2) = SUM(2) + ERR(2)
         CALL RNDOUT(RSLT, .TRUE., .TRUE.)
      ELSE
         RSLT(1) = SUM(1)
         RSLT(2) = SUM(2)
      ENDIF
!
      RETURN
      END

! **********************************************************
      SUBROUTINE ISNVAL (CODE,RSLT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Chenyi Hu
!
!   Computer and Mathematical Sciences Department
!   University of Houston-Downtown
!   Houston, TX 77002
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   October 20, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!   This routine transforms the result from the reduced argument
!   computation for the arcsine function back to the value corresponding
!   to the original argument.  The type of argument reduction which was
!   applied is given in the integer variable CODE.  Both the reduced
!   value input to this routine and the final value output from this
!   routine are stored in RSLT.
!
!***********************************************************************
!
! Argument declarations --
!
      INTEGER CODE
      DOUBLE PRECISION RSLT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TEMP(2)
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  ADD, MULT, RNDOUT
!
!
!***********************************************************************
!
! Beginning of executable statements --
!
      IF (CODE .EQ. 11) RETURN
!
   20 CONTINUE
      TEMP(1) = RSLT(1) * RSLT(1)
      TEMP(2) = RSLT(2) * RSLT(2)
      CALL RNDOUT(TEMP, .TRUE., .TRUE.)
!
!  rslt <-- rslt * (3 + 4 rslt^2)
!
      CALL MULT(FOUR,TEMP,TEMP)
      CALL ADD(THREE,TEMP,TEMP)
      CALL MULT(RSLT,TEMP,RSLT)
!
      IF (CODE .EQ. 12) THEN
         RETURN
      ELSE IF (CODE .EQ. 13) THEN
         CODE = 12
         GO TO 20
      ELSE IF (CODE .EQ. 14) THEN
         CODE = 13
         GO TO 20
      ENDIF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE ISQRT(X,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This subroutine returns the interval value of the square root of X
!  in RESULT.  It uses the routine RSQRT to get rigorous values at the
!  endpoints.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION V(2), RL(2), RR(2)
!
!***********************************************************************
!
!  Common block declarations --
!
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Package-supplied functions and subroutines --

!     ERRTST, RSQRT
!
!***********************************************************************
!
! Beginning of executable statements --
!
!  Identifying code for this routine --
!
      IROUT = 5
      IERR = 0
!
      IF (X(2).LT.X(1)) THEN
         IERR = 1
         CALL ERRTST(X)
      END IF
!
      IF (X(1).LT.0D0) THEN
         IERR = 7
         CALL ERRTST(X)
         RETURN
      ELSE IF (X(1).EQ.0D0) THEN
         RL(1) = 0D0
         RL(2) = 0D0
      ELSE
         V(1) = X(1)
         V(2) = X(1)
         CALL RSQRT(V,RL)
      END IF
!
      IF (X(2).EQ.0D0) THEN
         RR(1) = 0D0
         RR(2) = 0D0
      ELSE
         V(1) = X(2)
         V(2) = X(2)
         CALL RSQRT(V,RR)
      END IF
!
      RESULT(1) = RL(1)
      RESULT(2) = RR(2)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE POWER(AA,NDUM,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Manuel Novoa III
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   September 29, 1987
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!   Modified by Manuel Novoa III on March 13, 1990 to fix a bug which
!   caused overestimation for N odd and zero in the interval's interior.
!
!   Modified August, 1991 and April, 1992 by Kaisheng Du and
!   R. Baker Kearfott to allow computation of negative integer powers.
!
!***********************************************************************
!
! Called by --
!
!  Any routine requiring computation of a positive integer power of
!  an interval.
!
!***********************************************************************
!
! Function --
!
!  This routine computes the NDUM-th power, of the interval A, where
!  NDUM is an integer.  It simulates directed roundings with the routine
!  RNDOUT; the interval result should contain the interval which would
!  have been obtained with exact interval arithmetic.  However, in
!  general it will not be the smallest possible machine-representable
!  such containing interval.  See the documentation in subroutine RNDOUT
!  for more detailed information.
!
!  This routine can clearly be made more efficient and more readable on
!  Fortran systems for which the intrinsic represented by A**N is
!  optimally accurate.  See the documentation in the subsidiary routine
!  RRPOWR for more information.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION AA(2), RESULT(2)
      INTEGER NDUM
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION B(2), TMP(2)
      LOGICAL EVEN ,L
      DOUBLE PRECISION TEMP
      INTEGER N
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!    MAX, MOD
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!    ERRTST, MULT, RNDOUT, RRPOWR
!
!***********************************************************************
!
! Beginning of executable statements --
!
!
!  Identifying code for this routine --
!
      IROUT = 7
      IERR = 0
!
      IF (AA(2).LT.AA(1)) THEN
         IERR = 1
         CALL ERRTST(AA)
      END IF
!
!  If N is equal to 0 then return the point interval [1,1] --
!
      N = NDUM
!
      IF (N.EQ.0) THEN
         RESULT(1) = ONE(1)
         RESULT(2) = ONE(2)
         RETURN
      END IF
!
!     If N = 1 then return AA.
!
      IF (N.EQ.1) THEN
         RESULT(1) = AA(1)
         RESULT(2) = AA(2)
         RETURN
      END IF
!
!  If N < 0, let N = -N.
!
      L = (N.LT.0)
      IF(L) THEN
         IF (AA(1).LT.ZERO(1) .AND. AA(2).GT.ZERO(2)) THEN
            IERR=8
            CALL ERRTST(AA)
            RESULT(1) = NEGINF
            RESULT(2) = POSINF
            RETURN
         END IF
         N = -N
      END IF
!
!  If N > 1, check cases for 0 in or not in the interval --
!
      EVEN = (MOD(N,2).EQ.0)
!
      IF (AA(1).GT.ZERO(2)) THEN
         CALL RRPOWR(AA(1),N,TMP)
         B(1) = TMP(1)
         CALL RRPOWR(AA(2),N,TMP)
         B(2) = TMP(2)
      ELSE IF (AA(2).LT.ZERO(1)) THEN
         IF (EVEN) THEN
            CALL RRPOWR(AA(2),N,TMP)
            B(1) = TMP(1)
            CALL RRPOWR(AA(1),N,TMP)
            B(2) = TMP(2)
         ELSE
            CALL RRPOWR(AA(1),N,TMP)
            B(1) = TMP(1)
            CALL RRPOWR(AA(2),N,TMP)
            B(2) = TMP(2)
         END IF
      ELSE IF (EVEN) THEN
         B(1) = ZERO(1)
         TMP(2) = -AA(1)
         CALL RNDOUT(TMP,.FALSE.,.TRUE.)
         TEMP = MAX(TMP(2),AA(2))
         CALL RRPOWR(TEMP,N,TMP)
         B(2) = TMP(2)
      ELSE
         CALL RRPOWR(AA(1),N,TMP)
         B(1) = TMP(1)
         CALL RRPOWR(AA(2),N,TMP)
         B(2) = TMP(2)
      END IF
!
      IF (L) THEN
         IF(B(1).GT.0D0)  THEN
            TEMP = B(1)
            B(1) = ONE(1)/B(2)
            B(2) = ONE(2)/TEMP
            CALL RNDOUT(B,.TRUE.,.TRUE.)
         ELSE
            TEMP = B(1)
            B(1) = ONE(2)/B(2)
            B(2) = ONE(1)/TEMP
            CALL RNDOUT(B,.TRUE.,.TRUE.)
         END IF
      END IF
!
      RESULT(1) = B(1)
      RESULT(2) = B(2)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE RRPOWR(AA,NDUM,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   November 14, 1992
!
! Part of the generalized bisection package
! (interval arithmetic subpackage).
!
!***********************************************************************
!
! Called by --
!
!  POWER
!
!***********************************************************************
!
! Function --
!
!  This routine computes the NDUM-th power, of the interval point AA,
!  using simulated directed roundings to rigorously bound the roundoff
!  error.  The results is placed in RESULT.  The reason for this routine
!  is to allow rigor on machines on which exponentiation by an integer
!  is not accurate to within one digit in the last place.  The POWER
!  routine can certainly be made more efficient without this routine
!  if the Fortran intrinsic A**NDUM is optimally accurate.  This routine
!  can also be made more efficient by rewriting it in Fortran 90, and
!  making good use of recursion.  (See ``Programmer's Guide to Fortran
!  90, W. S. Brainerd, C. H. Goldberg, J. C. Adams, McGraw Hill, 1990
!  p. 222 ff.)
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION AA, RESULT(2)
      INTEGER NDUM
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION X(2)
      INTEGER I
!
!***********************************************************************
!
! Beginning of executable statements --
!
      X(1) = AA
      X(2) = AA
!
      RESULT(1) = X(1)
      RESULT(2) = X(2)
      IF (NDUM.EQ.1) RETURN
!
      DO 10 I = 2,NDUM,1
         CALL MULT(RESULT,X,RESULT)
   10 END DO
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE RCOS(XA,VAL)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This routine is used to compute one of the endpoints of the
!  cosine of the interval X.  It is assumed that X is an interval of
!  very small width.  In particular, it is assumed that the length of
!  the interval is at most pi/16.  It is also assumed that the
!  absolute value of the argument is no greater than MAXX / pi,
!  where MAXX is the largest representable integer in the machine.
!  Since this routine doesn't check these conditions, it is assumed
!  that the calling routine has done so.
!
!  This routine is usually called by ISIN, which returns interval
!  values of the sin for general interval arguments.
!
!  This routine translates the interval to [0,pi/8], and possibly
!  uses a double angle formula twice.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XA(2), VAL(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  XA    is the interval argument
!       (INPUT)
!
!  VAL   is the resulting interval cosine value
!        (OUTPUT)
!
!***********************************************************************
!
! Internal variable declarations --
!
      LOGICAL L2
      INTEGER N
      DOUBLE PRECISION  DN, ERTRM(2), TEM(2), TMP, V(2), VALT(2),       &
     &                  VSQ(2), X(2)
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MXULP, TOL0, TTINY2
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
!  The above common block holds machine parameters which are set in
!  SIMINI and used here.
!
!  Variable descriptions
!
!  MXULP       (machine epsilon)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TTINY2      2 * (smallest representable positive machine number)
!                * (maximum error in ULP's of the floating pt. op's)
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines
!
!     DBLE, MOD, NINT
!
!***********************************************************************
!
! Package-supplied functions and subroutines
!
!     ADD, MULT, POWER, RNDOUT, SCLMLT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
      X(1) = XA(1)
      X(2) = XA(2)
      IF (X(2).LT.0D0) THEN
         TMP = X(2)
         X(2) = -X(1)
         X(1) = -TMP
      END IF
!
!  Reduce the angle to between 0 and pi using translation by a
!  multiple of pi/2, making sure that the resulting interval is either
!  entirely positive or entirely negative.
!
      CALL MULT (X,A,TEM)
      N = NINT((TEM(1)+TEM(2))/2D0)
      DN = DBLE(N)
!
      CALL SCLMLT (DN,PI,V)
      CALL SUB (X,V,V)
!
      IF (V(2).LT.0D0) THEN
         TMP = V(2)
         V(2) = -V(1)
         V(1) = -TMP
         CALL RNDOUT(V,.TRUE.,.TRUE.)
      END IF
!
!  Further reduce the angle to between 0 and pi/8, using angle halving.
!
      L2 = V(2) .GE. PI8(1)
      IF (L2) THEN
         V(1) = V(1)/4D0
         V(2) = V(2)/4D0
         CALL RNDOUT(V,.TRUE.,.TRUE.)
      END IF
!
!  Compute value of the cosine of the reduced angle.
!
      CALL POWER(V,2,VSQ)
!
!  val  = one + vsq * (-od2f + vsq * (od4f + vsq * (-od6f
!        + vsq * (od8f + vsq * (-od10f + vsq * od12f)))))
!
      CALL MULT (VSQ, OD12F, VAL)
      CALL SUB (VAL, OD10F, VAL)
      CALL MULT (VSQ, VAL, VAL)
      CALL ADD (OD8F, VAL, VAL)
      CALL MULT (VSQ, VAL, VAL)
      CALL SUB (VAL, OD6F, VAL)
      CALL MULT (VSQ,VAL, VAL)
      CALL ADD (OD4F, VAL, VAL)
      CALL MULT (VSQ, VAL, VAL)
      CALL SUB (VAL, OD2F, VAL)
      CALL MULT (VSQ, VAL, VAL)
      CALL ADD (ONE, VAL, VAL)
!
!  valt = - v **14 * od14f
!
      CALL POWER (V, 14, VALT)
      CALL MULT (VALT, OD14F, VALT)
      TMP = VALT(2)
      VALT(2) = -VALT(1)
      VALT(1) = -TMP
      CALL RNDOUT(VALT,.TRUE.,.FALSE.)

      ERTRM(1) = VALT(1)
      ERTRM(2) = 0D0
!
      CALL ADD (VAL,ERTRM,VAL)
!
!  Use double angle formulas to get back to original angle in
!  [0, pi/2]. ( val <-- 8val^4 - 8val^2 + 1 )
!
      IF (L2) THEN
         CALL POWER(VAL, 2, TEM)
         CALL SUB  (TEM, ONE, VAL)
         CALL MULT (VAL, TEM, VAL)
         CALL MULT (VAL, EIGHT, VAL)
         CALL ADD  (ONE, VAL, VAL)
      END IF
!
!  Change the sign of the result if the shift was by an odd multiple of
!  pi --
!
      IF (MOD(N,2).NE.0) THEN
         TMP = VAL(1)
         VAL(1) = -VAL(2)
         VAL(2) = -TMP
      END IF
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE REXP(X,VAL)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!***********************************************************************
!
! Function --
!
!  This routine is used to compute the interval value exp(X).
!  It is assumed that the width of the interval X is small.
!
!  This routine is usually called by the routine IEXP,
!  to compute interval values of endpoints.
!
!  This routine translates the midpoint to [0,1/16], by first making
!  it positive,  then subtracting the integer part,  and finally
!  possibly dividing by 2 a number of times (To translate
!  back to the original interval, it then squares the result a number of
!  times, multiplies by e**N for some N, and possibly takes a
!  reciprocal.)
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), VAL(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2)
      INTEGER N
      DOUBLE PRECISION ERTRM(2), G(2), VALT(2)
!
!***********************************************************************
!
! Common block declarations --
!
!
!  This common block holds machine parameters which are set in
!  SIMINI and used here.
!
!  Variable descriptions
!
!  MXULP       (machine epsilon)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TTINY2      2 * (smallest representable positive machine number)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TOL0        TTINY2 / MXULP
!
      DOUBLE PRECISION MXULP, TOL0, TTINY2
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!     DBLE, DMAX1, DMIN1, NINT
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!     ADD, MULT, POWER, RNDOUT, SCLMLT
!
!***********************************************************************
!
! Beginning of executable statements --
!
!   Reduce the arguments with the formula X = G + N /4, where G is
!   between -.25 and .25 --
!
      CALL MULT (FOUR,X,TMP(1))
      N=NINT(TMP(1))
      TMP(1) = DBLE(N)
      TMP(2) = DBLE(N)
      CALL RNDOUT(TMP,.TRUE.,.TRUE.)
      CALL IDIV(TMP,FOUR,TMP)
      CALL SUB(X, TMP, G)
!
!  Compute Taylor series approximation --
!
!     val = one + g*(one + g*(od2f + g*(od3f + g*(od4f + g*(od5f +
!    *     g*(od6f + g*(od7f + g*(od8f +g*(od9f +g*(od10f +g*od11f
!    *     )))))))))) --
!
      CALL MULT (G,OD11F,TMP)
      CALL ADD (OD10F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD9F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD8F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD7F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD6F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD5F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD4F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD3F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (OD2F,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (ONE,TMP,TMP)
      CALL MULT (G,TMP,TMP)
      CALL ADD (ONE,TMP,VAL)
!
!  Compute the error term and add to the approximation --
!
!     VALT = G**12 * OD12F * A3
!
      CALL POWER (G,12,TMP)
      CALL MULT (TMP,OD12F,TMP)
      CALL MULT (TMP,A3,VALT)
!
      ERTRM(1) = DMIN1(0D0,VALT(1))
      ERTRM(2) = DMAX1(0D0,VALT(2))
!
      CALL ADD (VAL,ERTRM,VAL)
!
!  Translate back to original interval -- val=val*e14**n --
!
      CALL POWER (E14,N,TMP)
      CALL MULT (VAL,TMP,VAL)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE RLOG(XX,RESULT)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Rebecca Yun
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This subroutine returns the interval value of the logarithm over the
!  interval XX in RESULT.  It is assumed that XX is relatively small,
!  and that it is possible.  The usual use of this routine is from the
!  routine ILOG, which checks the argument values and calls this routine
!  to get rigorous bounds on values at the end points.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION XX(2), RESULT(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      INTEGER K
      DOUBLE PRECISION ERRTRM(2), VAL(2), X(2), XI(2), Y(2)
!
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!     DBLE, INT, LOG
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!   ADD, IDIV, MULT, POWER, RNDOUT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
      X(1) = XX(1)
      X(2) = XX(2)
!
!  Argument reduction: rewrite X <-- R^K * X
!               where  1/R <= X <= R = EXP(1/16)  ---
!
!   It is O.K. to use the Fortran LOG function here, since it it not
!   critical that X be exactly between 1/R and R ---
!
      IF (X(1).GE.ONE(1)) THEN
         K=INT(LOG((X(1)+X(2))/TWO(1))*SXTEEN(1))
      ELSE
         K = INT(LOG((X(1)+X(2))/TWO(1))*SXTEEN(1))-1
      END IF
!
      IF (K.LT.JTINY2) K=JTINY2
!
!    X <-- X/(R^K) ---
!
      CALL POWER(ESXTNT,K,VAL)
      CALL IDIV(X,VAL,X)
!
!  Initialize the values ---
!
      CALL SUB (X,ONE,Y)
!
      RESULT(1) = ZERO(1)
      RESULT(2) = ZERO(2)
!
!  Use Horner's scheme to evaluate the thirteenth degree Taylor
!  polynomial centered at one --
!
      CALL MULT(Y,THRTTH,RESULT)
      CALL SUB(RESULT,TWLVTH,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(ELEVTH,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL SUB(RESULT,TENTH,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(NINTH,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL SUB(RESULT,EIGHTH,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(SEVNTH,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL SUB(RESULT,SIXTH,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(FIFTH,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL SUB(RESULT,FOURTH,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(THIRD,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL SUB(RESULT,OD2F,RESULT)
      CALL MULT(Y,RESULT,RESULT)
      CALL ADD(ONE,RESULT,RESULT)
      CALL MULT(Y,RESULT,RESULT)
!
!  Compute the error term ---
!
      IF (X(1).LT.ONE(1)) THEN
         XI(1) = X(1)
      ELSE
         XI(1) = ONE(1)
      END IF
!
!
      IF(X(2).GT.ONE(2)) THEN
         XI(2) = X(2)
      ELSE
         XI(2) = ONE(2)
      END IF
!
!    ERRTRM <-- Y^14 / [14*XI^14]
!
      CALL IDIV(Y,XI,ERRTRM)
      CALL POWER(ERRTRM,14,ERRTRM)
      CALL MULT(ERRTRM,FORTTH,ERRTRM)
!
!    RESULT <-- RESULT - ERRTRM
!
      CALL SUB (RESULT,ERRTRM,RESULT)
!
!  Undo the argument reduction: RESULT <-- RESULT+ K/16 --
!
      Y(1) = DBLE(K)
      Y(2) = Y(1)
      CALL RNDOUT(Y,.TRUE.,.TRUE.)
      CALL MULT(Y,SXTNTH,Y)
      CALL ADD (RESULT,Y,RESULT)
!
      RETURN
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE RSQRT(X,VAL)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine returns the square root of the interval X in VAL.  It is
!  assumed that the interval X has relatively small width, and that its
!  lower bound is positive.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2), VAL(2)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION AM1(2), TMP1(2), TMP2(2), ERRTRM(2)
      DOUBLE PRECISION XTMP(2), XMID(2), XOLD(2)
      DOUBLE PRECISION CMIN, CMAX, T
      LOGICAL L, LSMALL
      INTEGER ITER, N
!
!***********************************************************************
!
! Common blocks --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines --
!
!   LOG10, MAX, MIN, NINT
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!   ADD, IDIV, MULT, POWER, RNDOUT, SUB
!
!***********************************************************************
!
! Beginning of executable statements --
!
      XTMP(1) = X(1)
      XTMP(2) = X(2)
!
!  Reduce to between 1/5 and 5 by dividing by the appropriate power
!  of 10.  Use the Fortran logarithm function to get the power.  This is
!  O.K., provided the argument to this function is within range,
!  since it is not crucial that the argument be exactly in [1/5,5] --
!
      N = NINT( LOG10( (XTMP(1)+XTMP(2))/2D0 ) )
!
!  Treat values between the smallest number which can be reciprocated
!  and the smallest representable number (if these two values are
!  different) as a special case --
!
      IF (N.LT.ITINY2) THEN
         N=ITINY2
         XTMP(1) = TINY2
         XTMP(2) = TINY2
         LSMALL = .TRUE.
      ELSE
         LSMALL = .FALSE.
      END IF
      CALL POWER(TEN,-N,TMP1)
      CALL MULT(XTMP,TMP1,XTMP)
!
      L = .FALSE.
!
      IF (XTMP(2).LT.ONE(1)) L = .TRUE.
!
!  Take the reciprocal if initial value is greater than one, in
!  order to have stable behavior in the interval Newton method --
!
      IF (L .AND. .NOT.LSMALL) THEN
         T = ONE(1)/XTMP(2)
         XTMP(2) = ONE(2)/XTMP(1)
         XTMP(1) = T
         CALL RNDOUT(XTMP,.TRUE.,.TRUE.)
      END IF
!
!  The reduced argument can be less than zero, depending on the
!  rounding properties and size of the original argument --
!
      IF (XTMP(1).LE.ZERO(2)) THEN
         WRITE(6,*) 'IN RSQRT, ARGUMENT WAS DETERMINED TO BE TOO'
         WRITE(6,*) 'LARGE OR TOO SMALL.'
         STOP
      END IF
!
!  Make an initial estimate for the result using a degree two
!  Taylor series with remainder term, provided the interval is
!  sufficiently bounded away from zero
!  (i.e. V <-- 1 + 1/2 (c-1) - 1/8 (c-1)^2 + 1/32 I (c-1)^3 )
!
         VAL(1) = ONE(1)
         VAL(2) = ONE(2)
         CALL SUB(XTMP,ONE,AM1)
         CALL MULT(AM1,OD2F,TMP1)
         CALL ADD(VAL,TMP1,VAL)
         CALL POWER(AM1,2,TMP1)
         CALL MULT(TMP1,EIGHTH,TMP1)
         CALL SUB(VAL,TMP1,VAL)
!
         CMIN = MIN(XTMP(1),ONE(1))
         CMAX = MAX(XTMP(2),ONE(2))
         TMP1(1) = ONE(1)/CMAX
         TMP1(2) = ONE(2)/CMIN
         CALL RNDOUT(TMP1,.TRUE.,.TRUE.)
         CALL POWER(TMP1,2,TMP2)
         CALL POWER(TMP1,3,TMP1)
         ERRTRM(1) = MIN(TMP1(1),TMP2(1))
         ERRTRM(2) = MAX(TMP1(2),TMP2(2))
         CALL MULT(ERRTRM,SXTNTH,ERRTRM)
         CALL POWER(AM1,3,TMP1)
         CALL MULT(ERRTRM,TMP1,ERRTRM)
!
         CALL ADD(ERRTRM,VAL,VAL)
         VAL(1) = MAX(VAL(1),TINY2)
!
!  Do an interval Newton iteration until it becomes stationary --
!
      DO 10    ITER = 1,100
!
         XOLD(1) = VAL(1)
         XOLD(2) = VAL(2)
         XMID(1) = (VAL(1)+VAL(2))/2D0
         XMID(2) = XMID(1)
!
!       VAL <-- MID(VAL) - (MID(VAL)^2-A) / (TWO*VAL) --
!
         CALL POWER (XMID,2,TMP1)
         CALL SUB (TMP1,XTMP,TMP1)
         CALL MULT (TWO,VAL,TMP2)
         CALL IDIV(TMP1,TMP2,TMP2)
         CALL SUB (XMID,TMP2,VAL)
!
         VAL(1) = MAX(VAL(1),XOLD(1))
         VAL(2) = MIN(VAL(2),XOLD(2))
!
         IF (VAL(1).EQ.XOLD(1) .AND. VAL(2).EQ.XOLD(2)) GO TO 20
!
   10 END DO
   20 CONTINUE
!
!  Take the reciprocal of the result, if the argument was
!  reciprocated --
!
      IF (L.AND. .NOT. LSMALL) THEN
         T = VAL(1)
         VAL(1) = ONE(1)/VAL(2)
         VAL(2) = ONE(2)/T
         CALL RNDOUT(VAL,.TRUE.,.TRUE.)
      END IF
!
!  Transform back to the original position by multiplying by the
!  appropriate power of the square root of 10 --
!
      CALL POWER(SQT10,N,TMP1)
      CALL MULT(TMP1,VAL,VAL)
!
      IF (LSMALL) VAL(1)=ZERO(1)
!
      RETURN
      END

!*** miscmach.f
      SUBROUTINE ERRTST(X)
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Function --
!
!  This routine prints error conditions which have been signalled
!  in other routines.  It halts execution, changes X, or returns,
!  depending on the particular error and on the way the error handling
!  control flags are set.
!
!***********************************************************************
!
! Argument declarations --
!
      DOUBLE PRECISION X(2)
!
!***********************************************************************
!
! Argument descriptions -- (INPUT  = set on entry and not alterable)
!                          (OUTPUT = to be set by the routine)
!                          (I/O    = set on entry but alterable)
!
!  X  is an interval which depends on the error set (I/O)
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION T
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Beginning of executable statements --
!
      ISIG = 0
!
      IF (IERR.EQ.1) THEN
         ISIG = 2
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 2 ERROR:  LOWER BOUND ON INTERVAL IS GREATER',      &
     &    ' THAN UPPER BOUND'
            WRITE(IERPUN,100) 'X:',X
         END IF
         T = X(1)
         X(1) = X(2)
         X(2) = T
      ELSE IF (IERR.EQ.2) THEN
         ISIG = 1
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 1 ERROR: LOWER BOUND ON X FOR IEXP WOULD UNDERFLOW'
         END IF
      ELSE IF (IERR.EQ.3) THEN
         ISIG = 2
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 2 ERROR: LOWER BOUND ON X FOR IEXP WOULD OVERFLOW'
         END IF
      ELSE IF (IERR.EQ.4) THEN
         ISIG = 2
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 2 ERROR: UPPER BOUND ON X FOR IEXP WOULD OVERFLOW'
         END IF
      ELSE IF (IERR.EQ.5) THEN
         ISIG = 0
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'WARNING: Loss of accuracy in trig function due to the',      &
     &    ' argument range.'
            WRITE(IERPUN,*) 'X:',X
         END IF
      ELSE IF (IERR.EQ.6) THEN
         ISIG = 2
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 2 ERROR: ZERO IN DENOMINATOR IN ORDINARY INTERVAL', &
     &    ' DIVISION.'
            WRITE(IERPUN,*) 'DENOMINATOR:',X
         END IF
      ELSE IF (IERR.EQ.7) THEN
         ISIG = 3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 3 ERROR: ARGUMENT CONTAINS ZERO IN ELEMENTARY',     &
     &    ' FUNCTION.'
            WRITE(IERPUN,*) 'ARGUMENT:',X
         END IF
      ELSE IF (IERR.EQ.8) THEN
         ISIG = 2
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 2 ERROR: NEGATIVE POWER OF A ZERO-CONTAINING'
            WRITE(IERPUN,*)                                             &
     &    'INTERVAL IS UNDEFINED IN ORDINARY INTERVAL ARITHMETIC.'
            WRITE(IERPUN,*) 'ARGUMENT:',X
         END IF
      ELSE IF (IERR.EQ.9) THEN
         ISIG = 3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 3 ERROR: ARGUMENT TO INVERSE TRIG FUNCTION'
            WRITE(IERPUN,*)                                             &
     &    'MAY CONTAIN NUMBERS WHICH ARE GREATER THAN 1.'
            WRITE(IERPUN,*) 'ARGUMENT:',X
         END IF
      ELSE IF (IERR.EQ.10) THEN
         ISIG = 3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 3 ERROR: ARGUMENT TO INVERSE TRIG FUNCTION'
            WRITE(IERPUN,*)                                             &
     &    'MAY CONTAIN NUMBERS WHICH ARE LESS THAN -1.'
            WRITE(IERPUN,*) 'ARGUMENT:',X
         END IF
      ELSE IF (IERR.EQ.11) THEN
         ISIG = 3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &    'SEVERITY 3 ERROR: APPROXIMATING SERIES FOR REDUCED ARGUMENT',&
     &    ' DID NOT CONVERGE.'
            WRITE(IERPUN,*) 'CURRENT TERM IN SERIES:',X
         END IF
      ELSE IF (IERR.EQ.12) THEN
         ISIG = 3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*) 'INTERNAL ERROR;  NEGATIVE ARGUMENT:'
            WRITE(IERPUN,*) X
            WRITE(IERPUN,*) 'POSSIBLY, AN ARRAY WAS DIMENSIONED',       &
     &                   ' INCORRECTLY OR MEMORY WAS NOT PROPERLY'
            WRITE(IERPUN,*) 'ALLOCATED.  IF THE ERROR PERSISTS, THEN:'
            WRITE(IERPUN,*) 'Contact R. Baker Kearfott, Dept. Math.,'
            WRITE(IERPUN,*) 'USL Box 4-1010, Lafayette, LA 70504-1010.'
            WRITE(IERPUN,*) 'email: rbk@usl.edu'
         END IF
      ELSE IF (IERR.EQ.13) THEN
         ISIG = 0
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &        'Warning:  Disjoint intervals in intersection.'
            X(1) = POSINF
            X(2) = NEGINF
         END IF
      ELSE IF (IERR.EQ.14) THEN
         ISIG = 1
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &        'Power function with operands containing 0^0 occured.'
            WRITE(IERPUN,*)                                             &
     &        'Result is sent to [NEGINF,POSINF].'
         END IF
      ELSE IF (IERR.EQ.15) THEN
         ISIG = 1
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*)                                             &
     &        'Power function with operands of the form'
            WRITE(IERPUN,*)                                             &
     &        '[0,pos]^[pos,pos] occurred. Result is set to'
            WRITE(IERPUN,*)                                             &
     &        '[pos,POSINF].'
         END IF
      ELSE
         ISIG=3
         IF (IPRTCL.LE.ISIG) THEN
            WRITE(IERPUN,*) 'INTERNAL ERROR;  UNKNOWN ERROR TYPE,',IERR,&
     &                      ' IN ERRTST.'
            WRITE(IERPUN,*) 'POSSIBLY, AN ARRAY WAS DIMENSIONED',       &
     &                   ' INCORRECTLY OR MEMORY WAS NOT PROPERLY'
            WRITE(IERPUN,*) 'ALLOCATED.  IF THE ERROR PERSISTS, THEN:'
            WRITE(IERPUN,*) 'Contact R. Baker Kearfott, Dept. Math.,'
            WRITE(IERPUN,*) 'USL Box 4-1010, Lafayette, LA 70504-1010.'
            WRITE(IERPUN,*) 'email: rbk@usl.edu'
         END IF
      END IF

!
      IF (IPRTCL.LE.ISIG) THEN
         IF (IROUT.EQ.1) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IEXP.'
         ELSE IF (IROUT.EQ.2) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ICOS.'
         ELSE IF (IROUT.EQ.3) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IDIV.'
         ELSE IF (IROUT.EQ.4) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ILOG.'
         ELSE IF (IROUT.EQ.5) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ISQRT.'
         ELSE IF (IROUT.EQ.6) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ISINH.'
         ELSE IF (IROUT.EQ.7) THEN
            WRITE(IERPUN,*) 'Error occurred in routine POWER.'
         ELSE IF (IROUT.EQ.8) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IACOS.'
         ELSE IF (IROUT.EQ.9) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IASIN.'
         ELSE IF (IROUT.EQ.10) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ASNSER, ',       &
     &                      'an auxiliary routine called by IASIN.'
         ELSE IF (IROUT.EQ.11) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IATAN.'
         ELSE IF (IROUT.EQ.12) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ATNSER, ',       &
     &                      'an auxiliary routine called by IATAN.'
         ELSE IF (IROUT.EQ.13) THEN
            WRITE(IERPUN,*) 'Error occurred in routine ICAP, ',         &
     &                      'a utility function routine.'
         ELSE IF (IROUT.EQ.14) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IIPOWR.'
         ELSE IF (IROUT.EQ.15) THEN
            WRITE(IERPUN,*) 'Error occurred in routine IVL2, ',         &
     &                      'a utility function routine.'
         ELSE
            ISIG = 3
            WRITE(IERPUN,*) 'INTERNAL ERROR;  UNKNOWN ROUTINE NUMBER,', &
     &                      IROUT, ' IN ERRTST.'
            WRITE(IERPUN,*) 'POSSIBLY, AN ARRAY WAS DIMENSIONED',       &
     &                   ' INCORRECTLY OR MEMORY WAS NOT PROPERLY'
            WRITE(IERPUN,*) 'ALLOCATED.  IF THE ERROR PERSISTS, THEN:'
            WRITE(IERPUN,*) 'Contact R. Baker Kearfott, Dept. Math.,'
            WRITE(IERPUN,*) 'USL Box 4-1010, Lafayette, LA 70504-1010.'
            WRITE(IERPUN,*) 'email: rbk@usl.edu'
         END IF
      END IF
!
      IF (ISIG.GE.ISEVER) STOP
!
      RETURN
!
  100 FORMAT(1X,A,3X,D23.15,2X,D23.15)
      END
!***********************************************************************
!***********************************************************************
      SUBROUTINE SIMINI
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
! Written by:
!
!   R. Baker Kearfott
!
!         and
!
!   Kaisheng Du
!
!   Department of Mathematics
!   U.S.L. Box 4-1010
!   Lafayette, LA 70504
!
!   April 11, 1992
!
! Part of the interval elementary function library
!
!
!***********************************************************************
!
! Called by -- GENBIS
!
!***********************************************************************
!
! Function --
!
!  This routine sets certain machine parameters used to simulate
!  directed roundings in a reasonably transportable way.  In
!  particular, it sets the amount by which to decrease the left endpoint
!  and increase the right endpoint of an interval computed using usual
!  floating point arithmetic to guarantee that the resulting interval
!  will contain the result which would have been obtained with true
!  interval arithmetic.
!
!
!  In addition to setting certain parameters, certain constants, such as
!  pi and the natural logarithm base e, are set here.  The data
!  for the representation for these constants which appears here
!  should be accurate at least to the number of significant figures
!  present in the machine.  If it is possible to input a binary form
!  for the double precision representation of this data, then that
!  form should be given to be correct to all representable digits.
!
!  This routine assumes that the four elementary floating point
!  operations, and unary negation, will give results with a maximum
!  error of one ULP (unit in the last place).  If this is not so, change
!  the value of MAXERR in the data statement below to the maximum number
!  of ULP's by which a floating point result can differ from the true
!  result (for '+', '-', '*', '/' and conversion from integer
!  to double precision).
!
!  An additional assumption is that the standard routines MIN and MAX
!  return exact values corresponding to one of their arguments.
!
!  When determining the maximum error of the result A op B, where
!  A and B are floating point numbers, we assume that A and B are
!  represented exactly.  For example, if A and B are almost equal,
!  then it is not unreasonable to assume that A - B, where the
!  subtraction is a floating point subtraction, is within a few
!  units of the last place of the true result.
!
!  Throughout the elementary function routines, it is assumed that
!  storing the double precision expressions 0D0 leads to an exact
!  floating point representation of 0. It is also assumed that a
!  floating point assignment statement (such as A=B) causes exactly the
!  same value to be in A as in B.
!
!  Note:  On some machines, an underflow naturally occurs when this
!         routine is executed.  (It would happen in the computation of
!         TINY2.)  This is not a problem.
!
!  If SIMINI is installed correctly, then the conclusions this
!  package prints out will have mathematical rigor.
!
!***********************************************************************
!
! Common block declarations --
!
      DOUBLE PRECISION MXULP, TTINY2, TOL0
      COMMON /MACH1/ MXULP, TTINY2, TOL0
!
!  This common block holds machine parameters which are set here
!  and used in RNDOUT.
!
!  Variable descriptions
!
!  MXULP       (machine epsilon)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TTINY2      2 * (smallest representable positive machine number)
!                * (maximum error in ULP's of the floating pt. op's)
!
!  TOL0        TTINY2 / MXULP
!
      DOUBLE PRECISION TINY, TEST
      COMMON /MACH2/ TINY, TEST
!
!  Common block MACH2 stores machine constants used in RNDOUT.  TINY
!  is the smallest representable machine number, while TEST is the
!  smallest number which can be safely rounded to something other than
!  zero.
!
      DOUBLE PRECISION MAXX, MXLGM1, NEGINF, POSINF, A(2), PI(2), E(2), &
     &   PI2(2), PI3(2), PI4(2), PI6(2), PI8(2), E14(2), A3(2), ONE(2), &
     &   OD2F(2), OD3F(2), OD4F(2), OD5F(2), OD6F(2), OD7F(2), OD8F(2), &
     &   OD9F(2), OD10F(2), OD11F(2), OD12F(2), OD14F(2), EIGHT(2),     &
     &   ZERO(2), TWO(2), THREE(2), FOUR(2), SXTNTH(2), NINE(2), TEN(2),&
     &   TWOT7(2), SQT10(2), ESXTNT(2), SXTEEN(2), THIRD(2), FOURTH(2), &
     &   FIFTH(2), SIXTH(2), SEVNTH(2), EIGHTH(2), NINTH(2), TENTH(2),  &
     &   ELEVTH(2), TWLVTH(2), THRTTH(2), FORTTH(2), TT7TH(2), TINY2,   &
     &   CBTEP, SQT3(2), ODSQT3(2)
!
      COMMON /MTHCNS/  MAXX, MXLGM1, NEGINF, POSINF, A, PI, E, PI2, PI3,&
     &   PI4, PI6, PI8, ONE, OD2F, OD3F, E14, A3, OD4F, OD5F, OD6F,     &
     &   OD7F, OD8F, OD9F, OD10F, OD11F, OD12F, OD14F, EIGHT, ZERO, TWO,&
     &   THREE, FOUR, SXTNTH, NINE, TEN, TWOT7, SQT10, ESXTNT, SXTEEN,  &
     &   THIRD, FOURTH, FIFTH, SIXTH, SEVNTH, EIGHTH, NINTH, TENTH,     &
     &   ELEVTH, TWLVTH, THRTTH, FORTTH, TT7TH, TINY2, CBTEP, SQT3,     &
     &   ODSQT3
!
      INTEGER ITINY2, JTINY2
!
      COMMON /IMATH/   ITINY2, JTINY2
!
!
!  The above common blocks hold mathematical constants which are used in
!  the elementary function routines.
!
!  Variable descriptions
!
!  MAXX        a double precision representation of the largest
!              representable integer
!
!  MXLGM1      an approximation to the logarithm of .25 times the
!              largest representable machine number, minus 1. This
!              should be a rigorous lower bound on the logarithm of the
!              largest representable machine number.  Its default
!              computation in SIMINI is to use the Fortran LOG function.
!              This should be changed if the Fortran LOG function is not
!              sufficiently accurate.
!
!
!  A           an interval enclosure for 1/pi
!
!  PI          an interval enclosure for pi
!
!  PI2         an interval enclosure for pi/2
!
!  PI3         an interval enclosure for pi/3
!
!  PI4         an interval enclosure for pi/4
!
!  PI6         an interval enclosure for pi/6
!
!  PI8         an interval enclosure for pi/8
!
!  E           an interval enclosure for E
!
!  E14         an interval enclosure for e^{1/4}
!
!  ESXTNT      an interval enclosure for e^{1/16}
!
!  SQT10       an interval enclosure for SQRT(10)
!
!  TINY2       the maximum of the smallest machine number and the
!              reciprocal of the largest machine number.  This quantity
!              is checked to avoid overflow in certain places.
!
!  CBTEP       an approximation to the cube root of six times the
!              largest distance between numbers.
!
!  SQT3        an interval enclosure for the square root of 3.
!
!  ODSQT3      the reciprocal of the square root of 3 times
!              1+ 100*MXULP, used in argument reduction in the
!              arctangent routine.
!
!  ITINY2      the logarithm base 10 of TINY2 truncated to an integer.
!
!  JTINY2      sixteen times the logarithm base e of TINY2 truncated to
!              an integer.
!
!  See the statements in the routine SIMINI for the definitions of
!  the other constants.
!
      INTEGER ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
      COMMON /ERFLGS/ ISIG, IERR, IROUT, IPRTCL, ISEVER, IERPUN
!
! The above common block stores signalling information for error
! conditions.  In particular,
!
! ISIG   is set to 0 at the beginning of SIMINI, and is reset
!        to the severity level of the error (1, 2, or 3) when an error
!        condition occurs.  The user may reset the flag to zero after an
!        error condition, depending on the error, if ISEVER is set to
!        allow execution after an error of its severity.
!
! IERR   is the number of the error condition, if an error occurred in
!        the last routine in INTLIB with error checking.  If no errors
!        occurred in the last such routine called, then IERR is zero.
!        The specific error conditions associated with particular error
!        numbers is defined in the routine ERRTST.
!
! IROUT  is the code number of the package routine in which the error
!        occurred
!
! IPRTCL controls the level of error which is printed.  IPRCTL=0 prints
!        all levels, while IPRCTL=1 prints only errors of level 1 or
!        greater.  If IPRTCL=4, then no error information is printed.
!
! ISEVER gives the error level which will stop execution.  ISEVER=0
!        causes any error to stop execution, while ISEVER>=3 causes only
!        errors of severity 3 to stop execution.  Generally, severity 3
!        errors correspond to when it is impossible to assign a result
!        with any meaningful interpretation.
!
! IERPUN is the Fortran unit number to which errors should be printed.
!
!***********************************************************************
!
! Fortran-supplied functions and subroutines -- DBLE, INT, MAX
!
!***********************************************************************
!
! Package-supplied functions and subroutines --
!
!  D1MACH, I1MACH (the SLATEC routines for machine constants)
!
!  RNDOUT
!
      DOUBLE PRECISION D1MACH
      INTEGER I1MACH
!
!***********************************************************************
!
! User-supplied functions and subroutines -- none
!
!***********************************************************************
!
! I/O functions -- none
!
!***********************************************************************
!
! Internal variable declarations --
!
      DOUBLE PRECISION TMP(2)
!
!***********************************************************************
!
! Internal constant declarations --
!
      INTEGER MAXERR
      DATA MAXERR/1/
!
!***********************************************************************
!
! Internal constant descriptions --
!
!  MAXERR      is the maximum number of ULP's (units in the last
!              place) by which a result of one of the floating point
!              operations (+, -, *, /, ** N) can differ from the
!              true result.  (See explanation above.)
!
!********** WARNING: The value of MAXERR is machine dependent and
!                    must be manually set.
!
!***********************************************************************
!
! Beginning of executable statements --
!
      TINY = D1MACH(1)
      MXULP = DBLE(MAXERR) * D1MACH(4)
      TEST = TINY/(1D0-3D0*MXULP)
      TTINY2 = 2D0 * DBLE(MAXERR) * TINY
      TOL0 = TTINY2 / MXULP
!
      TMP(1) = 1D0/D1MACH(2)
      TMP(2) = TMP(1)
      CALL RNDOUT (TMP,.TRUE.,.TRUE.)
      TINY2 = MAX(TINY,TMP(2))
      CBTEP = (6D0*D1MACH(4))**(1D0/3D0)
      ITINY2 = INT(LOG10(TINY2))
      JTINY2 = 16*INT(LOG(TINY2))
!
      MAXX = DBLE(I1MACH(9))
!
      MXLGM1 = LOG(D1MACH(2)) -LOG(4D0) - 1D0
      NEGINF = -D1MACH(2)
      POSINF =  D1MACH(2)
!
      A(1) = 0.31830988618379067153776752674502864D+00
      A(2) = 0.31830988618379067153776752674502864D+00
      CALL RNDOUT(A,.TRUE.,.TRUE.)
!
      PI(1) = 0.31415926535897932384626433832795028D+01
      PI(2) = 0.31415926535897932384626433832795028D+01
      CALL RNDOUT(PI,.TRUE.,.TRUE.)

      E(1) =  0.27182818284590452353602874713526625D+01
      E(2) =  0.27182818284590452353602874713526625D+01
      CALL RNDOUT(E,.TRUE.,.TRUE.)
!
      PI2(1) = 0.15707963267948966192313216916397499D+01
      PI2(2) = 0.15707963267948966192313216916397499D+01
      CALL RNDOUT(PI2,.TRUE.,.TRUE.)
!
      PI3(1) = 1.0471975511965977461542144610931676D0
      PI3(2) = 1.0471975511965977461542144610931676D0
      CALL RNDOUT(PI3,.TRUE.,.TRUE.)
!
      PI4(1) = 0.7853981633974483096156608458198757D0
      PI4(2) = 0.7853981633974483096156608458198757D0
      CALL RNDOUT(PI4,.TRUE.,.TRUE.)
!
      PI6(1) = 0.5235987755982988730771072305465838D0
      PI6(2) = 0.5235987755982988730771072305465838D0
      CALL RNDOUT(PI6,.TRUE.,.TRUE.)
!
      PI8(1) = 0.39269908169872415480783042290993785D+00
      PI8(2) = 0.39269908169872415480783042290993785D+00
      CALL RNDOUT(PI8,.TRUE.,.TRUE.)
!
      E14(1) = 0.12840254166877414840734205680624368D+01
      E14(2) = 0.12840254166877414840734205680624368D+01
      CALL RNDOUT(E14,.TRUE.,.TRUE.)
!
      A3(1) = 0.11331484530668263168290072278117947D+01
      A3(2) = 0.11331484530668263168290072278117947D+01
      CALL RNDOUT(A3,.TRUE.,.TRUE.)
!
      ONE  (1) = 1D0
      OD2F (1) = 1D0 /           2D0
      OD3F (1) = 1D0 /           6D0
      OD4F (1) = 1D0 /          24D0
      OD5F (1) = 1D0 /         120D0
      OD6F (1) = 1D0 /         720D0
      OD7F (1) = 1D0 /        5040D0
      OD8F (1) = 1D0 /       40320D0
      OD9F (1) = 1D0 /      362880D0
      OD10F(1) = 1D0 /     3628800D0
      OD11F(1) = 1D0 /    39916800D0
      OD12F(1) = 1D0 /   479001600D0
      OD14F(1) = 1D0 / 87178291200D0
!
      ONE  (2) = ONE  (1)
      OD2F (2) = OD2F (1)
      OD3F (2) = OD3F (1)
      OD4F (2) = OD4F (1)
      OD5F (2) = OD5F (1)
      OD6F (2) = OD6F (1)
      OD7F (2) = OD7F (1)
      OD8F (2) = OD8F (1)
      OD9F (2) = OD9F (1)
      OD10F(2) = OD10F(1)
      OD11F(2) = OD11F(1)
      OD12F(2) = OD12F(1)
      OD14F(2) = OD14F(1)
!
      CALL RNDOUT(  ONE,.TRUE.,.TRUE.)
      CALL RNDOUT( OD2F,.TRUE.,.TRUE.)
      CALL RNDOUT( OD4F,.TRUE.,.TRUE.)
      CALL RNDOUT( OD6F,.TRUE.,.TRUE.)
      CALL RNDOUT( OD8F,.TRUE.,.TRUE.)
      CALL RNDOUT(OD10F,.TRUE.,.TRUE.)
      CALL RNDOUT(OD12F,.TRUE.,.TRUE.)
      CALL RNDOUT(OD14F,.TRUE.,.TRUE.)
!
      EIGHT(1) = 8D0
      EIGHT(2) = 8D0
      CALL RNDOUT(EIGHT,.TRUE.,.TRUE.)
!
      ZERO(1) = 0D0
      ZERO(2) = 0D0
!
      TWO(1) = 2D0
      TWO(2) = 2D0
      CALL RNDOUT(TWO,.TRUE.,.TRUE.)
!
      THREE(1) = 3D0
      THREE(2) = 3D0
      CALL RNDOUT(THREE,.TRUE.,.TRUE.)
!
      FOUR(1) = 4D0
      FOUR(2) = 4D0
      CALL RNDOUT(FOUR,.TRUE.,.TRUE.)
!
      NINE(1) = 9D0
      NINE(2) = 9D0
      CALL RNDOUT(NINE,.TRUE.,.TRUE.)
!
      SXTNTH(1) = 1D0/16D0
      SXTNTH(2) = SXTNTH(1)
      CALL RNDOUT(SXTNTH,.TRUE.,.TRUE.)
!
      TEN(1) = 10D0
      TEN(2) = 10D0
      CALL RNDOUT(TEN,.TRUE.,.TRUE.)
!
      TWOT7(1) = 27D0
      TWOT7(2) = 27D0
      CALL RNDOUT(TWOT7,.TRUE.,.TRUE.)
!
      SQT10(1) = 3.1622776601683793319988935444327185D0
      SQT10(2) = 3.1622776601683793319988935444327185D0
      CALL RNDOUT(SQT10,.TRUE.,.TRUE.)
!
      ESXTNT(1) = 1.0644944589178594295633905946428909D0
      ESXTNT(2) = 1.0644944589178594295633905946428909D0
      CALL RNDOUT(ESXTNT,.TRUE.,.TRUE.)
!
      SXTEEN(1) = 16D0
      SXTEEN(2) = 16D0
      CALL RNDOUT(SXTEEN, .TRUE.,.TRUE.)
!
      THIRD(1) = 1D0/3D0
      THIRD(2) = 1D0/3D0
      CALL RNDOUT(THIRD,.TRUE.,.TRUE.)
!
      FOURTH(1) = .25D0
      FOURTH(2) = .25D0
      CALL RNDOUT(FOURTH,.TRUE.,.TRUE.)
!
      FIFTH(1) = .2D0
      FIFTH(2) = .2D0
      CALL RNDOUT(FIFTH,.TRUE.,.TRUE.)
!
      SIXTH(1) = 1D0/6D0
      SIXTH(2) = 1D0/6D0
      CALL RNDOUT(SIXTH,.TRUE.,.TRUE.)
!
      SEVNTH(1) = 1D0/7D0
      SEVNTH(2) = 1D0/7D0
      CALL RNDOUT(SEVNTH,.TRUE.,.TRUE.)
!
      EIGHTH(1) = .125D0
      EIGHTH(2) = .125D0
      CALL RNDOUT(EIGHTH,.TRUE.,.TRUE.)
!
      NINTH(1) = 1D0/9D0
      NINTH(2) = 1D0/9D0
      CALL RNDOUT(NINTH,.TRUE.,.TRUE.)
!
      TENTH(1) = .1D0
      TENTH(2) = .1D0
      CALL RNDOUT(TENTH,.TRUE.,.TRUE.)
!
      ELEVTH(1) = 1D0/11D0
      ELEVTH(2) = 1D0/11D0
      CALL RNDOUT(ELEVTH,.TRUE.,.TRUE.)
!
      TWLVTH(1) = 1D0/12D0
      TWLVTH(2) = 1D0/12D0
      CALL RNDOUT(TWLVTH,.TRUE.,.TRUE.)
!
      THRTTH(1) = 1D0/13D0
      THRTTH(2) = 1D0/13D0
      CALL RNDOUT(THRTTH,.TRUE.,.TRUE.)
!
      FORTTH(1) = 1D0/14D0
      FORTTH(2) = 1D0/14D0
      CALL RNDOUT(FORTTH,.TRUE.,.TRUE.)
!
      TT7TH(1) = 1D0/27D0
      TT7TH(2) = 1D0/27D0
      CALL RNDOUT(TT7TH,.TRUE.,.TRUE.)
!
      SQT3(1) = 1.7320508075688772935274463415058724D0
      SQT3(2) = 1.7320508075688772935274463415058724D0
      CALL RNDOUT(SQT3,.TRUE.,.TRUE.)
!
      ODSQT3(1) = 0.57735026918962576450914878050195746D0
      ODSQT3(2) = 0.57735026918962576450914878050195746D0
      CALL RNDOUT(ODSQT3,.TRUE.,.TRUE.)
      CALL SCLMLT(1D0+100D0*MXULP,ODSQT3,ODSQT3)
!
!  Set default values for the error checking routine --
!
      ISIG = 0
      IERR = 0
      IPRTCL = 0
      ISEVER = 3
      IERPUN = 6
!
      RETURN
      END
