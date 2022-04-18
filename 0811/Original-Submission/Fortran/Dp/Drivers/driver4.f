*     PROGRAM TMINU
*
*     TEST PROGRAM FOR THE SUBROUTINE PMINU
*
C     .. Scalars in Common ..
      INTEGER NADD,NDECF,NEXT,NFG,NFH,NFV,NIT,NRED,NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,FMIN,GMAX
      INTEGER I,IERR,IEXT,ITERM,NA,NF
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AF(200),RA(4000),RPAR(7),X(40)
      INTEGER IA(200),IPAR(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL PMINU,TIUD06
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
C      CALL TYTIM1(ITIME)
*
*     LOOP FOR 25 TEST PROBLEMS
*
      DO 30 NEXT = 1,25
*
*     CHOICE OF INTEGER AND REAL PARAMETERS
*
          DO 10 I = 1,6
              IPAR(I) = 0
   10     CONTINUE
          DO 20 I = 1,7
              RPAR(I) = 0.0D0
   20     CONTINUE
          IPAR(2) = 2
          IPAR(7) = 1
*
*     PROBLEM DIMENSION
*
          NF = 20
          NA = 65
*
*     INITIATION OF X AND CHOICE OF RPAR(7)
*
          CALL TIUD06(NF,NA,X,FMIN,RPAR(7),NEXT,IEXT,IERR)
          IF (IERR.NE.0) GO TO 30
*
*     SOLUTION
*
          CALL PMINU(NF,NA,X,AF,IA,RA,IPAR,RPAR,F,GMAX,IEXT,ITERM)
   30 CONTINUE
C      CALL TYTIM2(ITIME)
      STOP

      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF FA)
*
      SUBROUTINE FUN(NF,KA,X,FA)
*
*     FUNCTION EVALUATION
*
C     .. Scalar Arguments ..
      DOUBLE PRECISION FA
      INTEGER KA,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NEXT
C     ..
C     .. External Subroutines ..
      EXTERNAL TAFU06
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
C     ..
      CALL TAFU06(NF,KA,X,FA,NEXT)
      RETURN

      END
*
*     USER SUPPLIED SUBROUTINE (CALCULATION OF GA)
*
      SUBROUTINE DER(NF,KA,X,GA)
*
*     GRADIENT EVALUATION
*
C     .. Scalar Arguments ..
      INTEGER KA,NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GA(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER NEXT
C     ..
C     .. External Subroutines ..
      EXTERNAL TAGU06
C     ..
C     .. Common blocks ..
      COMMON /PROB/NEXT
C     ..
      CALL TAGU06(NF,KA,X,GA,NEXT)
      RETURN

      END
*
*     EMPTY SUBROUTINES
*
      SUBROUTINE FUNDER(NF,X,F,G)
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(*),X(*)
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
