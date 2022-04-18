C     PROGRAM TVARU
C
C     TEST PROGRAM FOR THE SUBROUTINE PVARU
C
C
C      CALL TYTIM1(ITIME)
C
C     LOOP FOR 20 TEST PROBLEMS
C
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NEXT,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,FMIN,GMAX
      INTEGER I,IERR,ITERM,NA,NF
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RA(8000),RPAR(7),X(50)
      INTEGER IPAR(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL PVARU,TIUD19
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      DO 30 NEXT = 1,20
C
C     CHOICE OF INTEGER AND REAL PARAMETERS
C
          DO 10 I = 1,7
              IPAR(I) = 0
   10     CONTINUE
          DO 20 I = 1,7
              RPAR(I) = 0.0D0
   20     CONTINUE
          IPAR(1) = 1
          IPAR(7) = 1
C
C     PROBLEM DIMENSION
C
          NF = 50
          NA = 0
C
C     INITIATION OF X AND CHOICE OF RPAR(7)
C
          CALL TIUD19(NF,X,FMIN,RPAR(7),NEXT,IERR)
          IF (IERR.NE.0) GO TO 30
C
C     THE USER SUPPLIED VALUES
C
          IF (IPAR(1).EQ.0) THEN
              IF (NEXT.EQ.25) IPAR(4) = 5
              RPAR(7) = 1.0D3
              IF (NEXT.EQ.1) RPAR(7) = 1.0D0
              IF (NEXT.EQ.2) RPAR(7) = 1.0D0
              IF (NEXT.EQ.3) RPAR(7) = 1.0D0
              IF (NEXT.EQ.8) RPAR(7) = 5.0D0
              IF (NEXT.EQ.9) RPAR(7) = 1.0D0
              IF (NEXT.EQ.10) RPAR(7) = 1.0D0
              IF (NEXT.EQ.11) RPAR(7) = 1.0D0
              IF (NEXT.EQ.13) RPAR(7) = 1.0D-1
              IF (NEXT.EQ.14) RPAR(7) = 2.0D-1
              IF (NEXT.EQ.15) RPAR(7) = 1.0D0
              IF (NEXT.EQ.16) RPAR(7) = 1.0D0
              IF (NEXT.EQ.17) RPAR(7) = 1.0D1
              IF (NEXT.EQ.18) RPAR(7) = 1.0D0
              IF (NEXT.EQ.19) RPAR(7) = 1.0D1
              IF (NEXT.EQ.23) RPAR(7) = 2.0D0
              IF (NEXT.EQ.24) RPAR(7) = 1.0D1
              IF (NEXT.EQ.25) RPAR(7) = 5.0D0

          ELSE
              IF (NEXT.EQ.21) IPAR(4) = 3
              IF (NEXT.EQ.22) IPAR(4) = 4
              IF (NEXT.EQ.25) IPAR(4) = 5
              RPAR(5) = 1.0D-9
              IF (NEXT.EQ.1) RPAR(5) = 1.0D0
              IF (NEXT.EQ.2) RPAR(5) = 2.0D0
              IF (NEXT.EQ.3) RPAR(5) = 2.0D0
              IF (NEXT.EQ.5) RPAR(5) = 1.0D0
              IF (NEXT.EQ.7) RPAR(5) = 2.0D0
              IF (NEXT.EQ.8) RPAR(5) = 1.0D-2
              IF (NEXT.EQ.10) RPAR(5) = 1.0D0
              IF (NEXT.EQ.13) RPAR(5) = 2.5D-1
              IF (NEXT.EQ.14) RPAR(5) = 1.0D-3
              IF (NEXT.EQ.15) RPAR(5) = 1.0D0
              IF (NEXT.EQ.16) RPAR(5) = 1.0D-3
              IF (NEXT.EQ.17) RPAR(5) = 2.5D-1
              IF (NEXT.EQ.18) RPAR(5) = 2.0D0
              IF (NEXT.EQ.19) RPAR(5) = 1.0D-1
              IF (NEXT.EQ.21) RPAR(5) = 1.0D-1
              IF (NEXT.EQ.23) RPAR(5) = 1.0D-5
              IF (NEXT.EQ.24) RPAR(5) = 1.0D-1
              IF (NEXT.EQ.25) RPAR(5) = 1.0D-1
              RPAR(7) = 1.0D3
              IF (NEXT.EQ.1) RPAR(7) = 1.0D0
              IF (NEXT.EQ.3) RPAR(7) = 1.0D0
              IF (NEXT.EQ.6) RPAR(7) = 1.0D0
              IF (NEXT.EQ.7) RPAR(7) = 1.0D0
              IF (NEXT.EQ.8) RPAR(7) = 2.0D-1
              IF (NEXT.EQ.9) RPAR(7) = 1.0D0
              IF (NEXT.EQ.10) RPAR(7) = 1.0D0
              IF (NEXT.EQ.11) RPAR(7) = 1.0D0
              IF (NEXT.EQ.12) RPAR(7) = 1.0D0
              IF (NEXT.EQ.13) RPAR(7) = 5.0D-1
              IF (NEXT.EQ.14) RPAR(7) = 2.0D-1
              IF (NEXT.EQ.15) RPAR(7) = 1.0D0
              IF (NEXT.EQ.16) RPAR(7) = 2.0D1
              IF (NEXT.EQ.17) RPAR(7) = 1.0D1
              IF (NEXT.EQ.18) RPAR(7) = 1.0D0
              IF (NEXT.EQ.19) RPAR(7) = 1.0D1
              IF (NEXT.EQ.23) RPAR(7) = 1.0D0
              IF (NEXT.EQ.24) RPAR(7) = 5.0D0
              IF (NEXT.EQ.25) RPAR(7) = 1.0D1
          END IF
C
C     SOLUTION
C
          CALL PVARU(NF,NA,X,RA,IPAR,RPAR,F,GMAX,ITERM)
   30 CONTINUE
C      CALL TYTIM2(ITIME)
      STOP

      END
C
C     USER SUPPLIED SUBROUTINE (CALCULATION OF F AND G)
C
      SUBROUTINE FUNDER(NF,X,F,G)
C
C     FUNCTION EVALUATION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NEXT
C     ..
C     .. External Subroutines ..
      EXTERNAL TFFU19,TFGU19
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
C     ..
      CALL TFFU19(NF,X,F,NEXT)
C
C     GRADIENT EVALUATION
C
      CALL TFGU19(NF,X,G,NEXT)
      RETURN

      END
*
*     EMPTY SUBROUTINES
*
      SUBROUTINE FUN(NF,KA,X,FA)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FA
      INTEGER KA,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(*)
C     ..
      RETURN

      END
      SUBROUTINE DER(NF,KA,X,GA)
C     .. Scalar Arguments ..
      INTEGER KA,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GA(*),X(*)
C     ..
      RETURN

      END
      SUBROUTINE HES(NF,X,H)
C     .. Scalar Arguments ..
      INTEGER NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(*),X(*)
C     ..
      RETURN

      END
