      REAL FUNCTION GAMINC ( A, X1, X2, GAM )
C
C  COMPUTE THE DIFFERENCE BETWEEN TWO MODIFIED INCOMPLETE
C  GAMMA FUNCTIONS FOR (A,X1) AND (A,X2), THEN MULTIPLY BY
C  EXP(X1).  THAT IS, COMPUTE THE INTEGRAL OF ABS(X)**(A-1.)
C  * EXP(X1-X) FROM X1 TO X2.  IF X1 .GT. X2, THEN X1 - X2 MUST
C  BE .LE. EXPLIM.
C  EXPLIM CAN BE A MACHINE DEPENDENT CONSTANT WHICH PREVENTS
C  EXPONENTIATION OVER- AND UNDERFLOWS.  IT IS USED HERE TO
C  SUPPRESS THE CALCULATION OF MIGAM(A,X2) WHEN THE VALUE OF
C  MIGAM(A,X2) IS INSIGNIFICANT.  THIS USAGE REQUIRES X2 +
C  EXPLIM .GE. X1.  (MIGAM IS AN ABBREVIATION FOR MODIFIED IN-
C  COMPLETE GAMMA FUNCTION.)
C  GAM IS THE COMPLETE GAMMA FUNCTION OF A SUPPLIED BY THE
C  CALLING PROGRAM.
C
C  FOR X .GT. 5., GAM - MIGAM(A,X) IS COMPUTED WITH A CONTINUED
C  FRACTION APPROXIMATION.  FOR ABS(X) .LE. 1.0, THE INTEGRAL
C  IS TRANSFORMED AND EXP(Q) IS APPROXIMATED WITH A CHEBYSHEV
C  SERIES SO THAT THE NEW INTEGRAL MAY BE DONE ANALYTICALLY.
C  FOR X .GT. -12, AND X .LT. 5. (ABS (X) .GT. 1.0), A CONTIN-
C  UED FRACTION APPROXIMATION IS USED.  FINALLY FOR X .LE.
C  -12., THE ASYMPTOTIC EXPANSION IS USED.
C
C  SGN IS A SWITCH WHICH, IF NONZERO, INDICATES WHETHER GAM
C  SHOULD BE ADDED OR SUBTRACTED FROM AN INTERMEDIATE RESULT.
C
C
C     .. Scalar Arguments ..
      REAL A,GAM,X1,X2
C     ..
C     .. Local Scalars ..
      REAL AZ,EXPDIF,EXPLIM,GAM1,SGN,TIM,Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP,SIGN
      DATA EXPLIM / 20. /
      Z = X1
      SGN = 0.
      TIM = -1.
      EXPDIF = 1.0
5     IF ( Z .NE. 0. ) GO TO 10
      GAM1 = 0.
      SGN = SGN + TIM
      GO TO 40
10    IF ( Z .LE. 5. ) GO TO 20
C  USE EQUATION 10.
      GAM1 = - EXPDIF * Z**A / ( Z + (1.-A)/(1.+1./(Z+(2.-A)/(1.+2.
     &  /(Z+(3.-A)/(1.+3./(Z+1.7) ))))))
      GO TO 40
20    AZ = ABS ( Z )
      IF ( Z .LE. -12. ) GO TO 30
      SGN = SGN + TIM
C  USE EQUATION 17.
      IF ( AZ .LE. 1. ) GAM1 = EXPDIF * Z / A *AZ**(A-1.)
     &  *(1.       +Z/(A+1.) *(.9999999+Z/(A+2.)
     &  *(.9999999 +Z/(A+3.) *(1.000008+Z/(A+4.)
     &  *(1.000005 +Z/(A+5.) *(.9994316+Z/(A+6.)
     &  *(.9995587 +Z/(A+7.) *(1.031684+Z/(A+8.)
     &  *1.028125))))))))
C  USE EQUATIONS 11 AND 12.  EVALUATION MUST BE DONE
C  IN DOUBLE PRECISION IF COMPUTER HAS 32 OR FEWER BITS
C  PER WORD.  ON AN IBM 360, D. P. EVALUATION IS FORCED
C  BY THE D. P. CONSTANTS IN CONTINUATION CARD 9.
      IF ( AZ .GT. 1. ) GAM1 = EXPDIF * Z / A * AZ**(A-1.)
     &  /(1.- A    *Z/( A     *(A+ 1.+   Z/((A+ 2.)
     &  *(1.-(A+1.)*Z/((A+ 2.)*(A+ 3.+2.*Z/((A+ 4.)
     &  *(1.-(A+2.)*Z/((A+ 4.)*(A+ 5.+3.*Z/((A+ 6.)
     &  *(1.-(A+3.)*Z/((A+ 6.)*(A+ 7.+4.*Z/((A+ 8.)
     &  *(1.-(A+4.)*Z/((A+ 8.)*(A+ 9.+5.*Z/((A+10.)
     &  *(1.-(A+5.)*Z/((A+10.)*(A+11.+6.*Z/((A+12.)
     &  *(1.-(A+6.)*Z/((A+12.)*(A+13.+7.*Z/((A+14.)
     &  *(1.00150-A*8.95D-5 + Z*(-.0337062+A*.0004182
     &  +Z*(000999294-A*.000104103))) )))) )))) ))))
     &  )))) )))) )))) ))))
      GO TO 40
C  USE EQUATION 10 AND SHANK'S E1 PROCESS ONCE.
30    GAM1 = - EXPDIF * AZ**(A-1.)*(1.+(A-1.)*(1.+(A-2.)*
     &  (1.+(A-3.)*(1.*(A-4.)*(1.*(A-5.)/(Z-A+6.))
     &  /Z)/Z)/Z)/Z)
40    IF ( TIM. GT. 0. ) GO TO 55
      GAMINC = GAM1
      IF ( ABS ( X1 - X2 ) .GT. EXPLIM ) GO TO 50
C  IF TRUE, CONTRIBUTION AT X2 IS .LT. 1.E-7 * (CONTR AT X1),
C  PROVIDED X2 .GT. X1.
      Z = X2
      EXPDIF = EXP ( X1 - X2 )
      TIM = 1.
      GO TO 5
50    GAM1 = 0.
55    GAMINC = GAM1 - GAMINC
      IF ( SGN .NE. 0. ) GAMINC = GAMINC - SIGN ( GAM * EXP (X1), SGN )
      RETURN
      END
