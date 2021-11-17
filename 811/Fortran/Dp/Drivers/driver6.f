*     PROGRAM TNEWU
*
*     TEST PROGRAM FOR THE SUBROUTINE PNEWU
*
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NEXT,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,FMIN,GMAX
      INTEGER I,IERR,IHES,ITERM,NA,NF
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RA(75000),RPAR(9),X(50)
      INTEGER IA(200),IPAR(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL PNEWU,TIUD19
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
C      CALL TYTIM1(ITIME)
*
*     LOOP FOR 20 TEST PROBLEMS
*
      DO 30 NEXT = 1,20
*
*     CHOICE OF INTEGER AND REAL PARAMETERS
*
          DO 10 I = 1,7
              IPAR(I) = 0
   10     CONTINUE
          DO 20 I = 1,9
              RPAR(I) = 0.0D0
   20     CONTINUE
          IPAR(1) = 2
          IPAR(2) = 4
          IF (NEXT.EQ.25) IPAR(3) = 50
          IF (NEXT.EQ.10) IPAR(4) = 4
          IF (NEXT.EQ.15) IPAR(4) = 4
          IPAR(7) = 1
          RPAR(8) = 1D-10
          IF (IPAR(1).EQ.1) THEN
              IF (NEXT.EQ.1) RPAR(8) = 0.5D0
              IF (NEXT.EQ.2) RPAR(8) = 1.0D-4
              IF (NEXT.EQ.3) RPAR(8) = 0.25D0
              IF (NEXT.EQ.4) RPAR(8) = 1.0D-2
              IF (NEXT.EQ.5) RPAR(8) = 1.0D-1
              IF (NEXT.EQ.8) RPAR(8) = 1.0D-1
              IF (NEXT.EQ.10) RPAR(8) = 1.0D-4
              IF (NEXT.EQ.13) RPAR(8) = 4.0D-2
              IF (NEXT.EQ.16) RPAR(8) = 1.0D-4
              IF (NEXT.EQ.17) RPAR(8) = 0.25D0
              IF (NEXT.EQ.18) RPAR(8) = 1.0D-2
              IF (NEXT.EQ.25) RPAR(8) = 5.0D-2

          ELSE
              IF (NEXT.EQ.1) RPAR(8) = 1.3D0
              IF (NEXT.EQ.2) RPAR(8) = 1.0D-3
              IF (NEXT.EQ.3) RPAR(8) = 1.0D-1
              IF (NEXT.EQ.5) RPAR(8) = 0.25D0
              IF (NEXT.EQ.8) RPAR(8) = 1.0D-4
              IF (NEXT.EQ.10) RPAR(8) = 0.25D0
              IF (NEXT.EQ.13) RPAR(8) = 0.25D0
              IF (NEXT.EQ.14) RPAR(8) = 1.0D-3
              IF (NEXT.EQ.15) RPAR(8) = 1.0D-3
              IF (NEXT.EQ.16) RPAR(8) = 1.0D-2
              IF (NEXT.EQ.17) RPAR(8) = 0.25D0
              IF (NEXT.EQ.18) RPAR(8) = 1.0D-2
              IF (NEXT.EQ.21) RPAR(8) = 0.5D0
              IF (NEXT.EQ.25) RPAR(8) = 5.0D-2
          END IF
*
*     PROBLEM DIMENSION
*
          NF = 50
          NA = 0
*
*     INITIATION OF X AND CHOICE OF RPAR(9)
*
          CALL TIUD19(NF,X,FMIN,RPAR(9),NEXT,IERR)
          IF (IERR.NE.0) GO TO 30
          IF (NEXT.EQ.3) RPAR(9) = 1.0D0
          IF (NEXT.EQ.4) RPAR(9) = 1.0D4
          IF (NEXT.EQ.5) RPAR(9) = 1.0D0
          IF (NEXT.EQ.8) RPAR(9) = 1.0D4
          IF (NEXT.EQ.10) RPAR(9) = 1.0D1
          IF (NEXT.EQ.17) RPAR(9) = 1.0D0
          IF (NEXT.EQ.21) RPAR(9) = 5.0D1
          IF (NEXT.EQ.22) RPAR(9) = 2.0D1
          IF (NEXT.EQ.24) RPAR(9) = 1.0D4
          IHES = 1
*
*     SOLUTION
*
          CALL PNEWU(NF,NA,X,IA,RA,IPAR,RPAR,F,GMAX,IHES,ITERM)
   30 CONTINUE
C      CALL TYTIM2(ITIME)
      STOP

      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF F AND G)
*
      SUBROUTINE FUNDER(NF,X,F,G)
*
*     FUNCTION EVALUATION
*
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
*
*     GRADIENT EVALUATION
*
      CALL TFGU19(NF,X,G,NEXT)
      RETURN

      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF H)
*
      SUBROUTINE HES(NF,X,H)
*
*     HESSIAN EVALUATION
*
C     .. Scalar Arguments ..
      INTEGER NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NEXT
C     ..
C     .. External Subroutines ..
      EXTERNAL TFHD19
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
C     ..
      CALL TFHD19(NF,X,H,NEXT)
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
