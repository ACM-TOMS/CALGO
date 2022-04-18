C     PROGRAM TNEWL
C
C     TEST PROGRAM FOR THE SUBROUTINE PNEWL
C
C      CALL TYTIM1(ITIME)
C
C     LOOP FOR 10 TEST PROBLEMS
C
C     .. Scalars in Common ..
      INTEGER IEXT,KAP,LAP,NAA,NADD,NDECF,NEXT,NFG,NFH,NFV,NIT,NRED,
     +        NREM,NRES
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,FMIN,GMAX
      INTEGER I,IERR,IHES,ITERM,NA,NB,NC,NF
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION CF(20),CG(200),CL(20),CU(20),RA(19000),RPAR(9),
     +                 X(40),XL(40),XU(40)
      INTEGER IA(200),IC(20),IPAR(7),IX(40)
C     ..
C     .. External Subroutines ..
      EXTERNAL PNEWL,TILD22
C     ..
C     .. Common blocks ..
      COMMON /PROB/IEXT,NEXT,NAA,KAP,LAP
      COMMON /STAT/NDECF,NRES,NRED,NREM,NADD,NIT,NFV,NFG,NFH
C     ..
      DO 30 NEXT = 1,10
C
C     CHOICE OF INTEGER AND REAL PARAMETERS
C
          DO 10 I = 1,7
              IPAR(I) = 0
   10     CONTINUE
          DO 20 I = 1,9
              RPAR(I) = 0.0D0
   20     CONTINUE
          IPAR(7) = 1
          RPAR(8) = 1.0D-10
          IF (NEXT.EQ.3) RPAR(8) = 0.5D0
          IF (NEXT.EQ.5) RPAR(8) = 1.0D-2
          IF (NEXT.EQ.7) RPAR(8) = 1.0D0
          IF (NEXT.EQ.8) RPAR(8) = 0.5D0
          IF (NEXT.EQ.9) RPAR(8) = 1.0D0
          IF (NEXT.EQ.11) RPAR(8) = 0.5D0
          IF (NEXT.EQ.12) RPAR(8) = 1.0D-4
          IF (NEXT.EQ.13) RPAR(8) = 1.0D-2
          IF (NEXT.EQ.14) RPAR(8) = 1.0D-3
          IF (NEXT.EQ.15) RPAR(8) = 5.0D-4
C
C     PROBLEM DIMENSION
C
          NF = 20
          NA = 0
          NB = 20
          NC = 15
          NAA = 165
C
C     INITIATION OF X AND CHOICE OF RPAR(9)
C
          CALL TILD22(NF,NAA,NB,NC,X,IX,XL,XU,IC,CL,CU,CG,FMIN,RPAR(9),
     +                NEXT,IEXT,IERR)
          IF (IERR.NE.0) GO TO 30
          IF (NEXT.EQ.1) RPAR(9) = 1.0D2
          IF (NEXT.EQ.2) RPAR(9) = 1.0D0
          IF (NEXT.EQ.3) RPAR(9) = 1.0D0
          IF (NEXT.EQ.5) RPAR(9) = 5.0D-2
          IF (NEXT.EQ.7) RPAR(9) = 5.0D-2
          IF (NEXT.EQ.8) RPAR(9) = 1.0D1
          IF (NEXT.EQ.9) RPAR(9) = 1.0D2
          IF (NEXT.EQ.10) RPAR(9) = 1.0D2
          IF (NEXT.EQ.11) RPAR(9) = 1.0D1
          IF (NEXT.EQ.12) RPAR(9) = 1.0D2
          IF (NEXT.EQ.13) RPAR(9) = 1.0D2
C
C     SOLUTION
C
          IHES = 1
          CALL PNEWL(NF,NA,NB,NC,X,IX,XL,XU,CF,IC,CL,CU,CG,IA,RA,IPAR,
     +               RPAR,F,GMAX,IHES,ITERM)
   30 CONTINUE
C      CALL TYTIM2(ITIME)
      STOP

      END
C
C     USER SUPPLIED SUBROUTINE (CALCULATION OF FA)
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
      INTEGER IEXT,KAP,LAP,NA,NEXT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FTEMP,FVAL
      INTEGER K,KA
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVNEG,TAFU22,TAGU22
C     ..
C     .. Common blocks ..
      COMMON /PROB/IEXT,NEXT,NA,KAP,LAP
C     ..
      DO 10 KA = 1,NA
          CALL TAFU22(NF,KA,X,F,NEXT)
          IF (IEXT.EQ.0 .AND. F.GE.0.0D0 .OR. IEXT.LT.0) THEN
              FTEMP = F
              K = 1

          ELSE
              FTEMP = -F
              K = -1
          END IF

          IF (KA.EQ.1 .OR. FVAL.LT.FTEMP) THEN
              FVAL = FTEMP
              KAP = KA
              LAP = K
          END IF

   10 CONTINUE
      F = FVAL
C
C     GRADIENT EVALUATION
C
      CALL TAGU22(NF,KAP,X,G,NEXT)
      IF (LAP.GE.0) THEN

      ELSE
          CALL MXVNEG(NF,G,G)
      END IF

      RETURN

      END
      SUBROUTINE HES(NF,X,H)
C
C     HESSIAN EVALUATION
C
C     .. Scalar Arguments ..
      INTEGER NF
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(*),X(*)
C     ..
C     .. Scalars in Common ..
      INTEGER IEXT,KAP,LAP,NA,NEXT
C     ..
C     .. External Subroutines ..
      EXTERNAL MXVNEG,TAHD22
C     ..
C     .. Common blocks ..
      COMMON /PROB/IEXT,NEXT,NA,KAP,LAP
C     ..
      CALL TAHD22(NF,KAP,X,H,NEXT)
      IF (LAP.GE.0) THEN

      ELSE
          CALL MXVNEG(NF* (NF+1)/2,H,H)
      END IF

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
