!     TEST PROGRAM FOR THE SUBROUTINE PSEDS
!
      INTEGER NF,NA,MA,IX(1000),IAG(2001),JAG(8000),IPAR(7),ISPAS,IPRNT,
     &ITERM
      DOUBLE PRECISION X(1000),XL(1000),XU(1000),AF(2000),RPAR(9),F,
     &GMAX
      INTEGER NEXT,IERR,I,ITIME
      INTEGER NITER,NFVAL,NSUCC
      COMMON /PROB/ NEXT
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      NITER=0
      NFVAL=0
      NSUCC=0
      CALL TYTIM1 (ITIME)
!
!     LOOP FOR 22 TEST PROBLEMS
!
      DO 30 NEXT=1,22
!
!     CHOICE OF INTEGER AND REAL PARAMETERS
!
        DO 10 I=1,7
          IPAR(I)=0
   10   CONTINUE
        DO 20 I=1,9
          RPAR(I)=0.0D0
   20   CONTINUE
        ISPAS=2
        IPRNT=1
!
!     PROBLEM DIMENSION
!
        NF=1000
        NA=2000
!
!     INITIATION OF X AND CHOICE OF RPAR(1) AND RPAR(6)
!
        CALL TIUB14 (NF, NA, MA, X, IAG, JAG, RPAR(6), RPAR(1), NEXT,
     &   IERR)
        IF (IERR.NE.0) GO TO 30
        CALL MXVINS (NF, 3, IX)
        CALL MXVSET (NF, -1.0D0, XL)
        CALL MXVSET (NF, 1.0D0, XU)
        IF (RPAR(6).EQ.0.0D0) IPAR(4)=1
        RPAR(1)=0.0D0
        IF (NEXT.EQ.5) RPAR(1)=1.0D1
        IF (NEXT.EQ.10) RPAR(1)=1.0D0
        IF (NEXT.EQ.12) RPAR(1)=1.0D0
!
!     SOLUTION
!
        CALL PSEDS (NF, NA, MA, X, IX, XL, XU, AF, IAG, JAG, IPAR, RPAR,
     &    F, GMAX, ISPAS, IPRNT, ITERM)
        NITER=NITER+NIT
        NFVAL=NFVAL+NFV
        IF (ITERM.GT.0.AND.ITERM.LT.9) NSUCC=NSUCC+1
   30 CONTINUE
      WRITE (6,40) NITER,NFVAL,NFVAL,NSUCC
   40 FORMAT (' NITER =',I5,3X,' NFVAL =',I5,3X,' NGVAL =',I5,3X,' NSUCC
     & =',I5)
      CALL TYTIM2 (ITIME)
      STOP
      END
!     USER SUPPLIED SUBROUTINE (CALCULATION OF FA)
!
      SUBROUTINE FUN (NF, KA, X, FA)
      INTEGER NF,KA
      DOUBLE PRECISION X(*),FA
      INTEGER NEXT
      COMMON /PROB/ NEXT
!
!     FUNCTION EVALUATION
!
      CALL TAFU14 (NF, KA, X, FA, NEXT)
      RETURN
      END
!     USER SUPPLIED SUBROUTINE (CALCULATION OF GA)
!
      SUBROUTINE DFUN (NF, KA, X, GA)
      INTEGER NF,KA
      DOUBLE PRECISION X(*),GA(*)
      INTEGER NEXT
      COMMON /PROB/ NEXT
!
!     GRADIENT EVALUATION
!
      CALL TAGU14 (NF, KA, X, GA, NEXT)
      RETURN
      END
!     EMPTY SUBROUTINES
!
      SUBROUTINE OBJ (NF, X, FF)
      INTEGER NF
      DOUBLE PRECISION X(*),FF
      NF=1
      FF=X(1)
      END
      SUBROUTINE DOBJ (NF, X, GF)
      INTEGER NF
      DOUBLE PRECISION X(*),GF(*)
      NF=1
      GF(1)=X(1)
      RETURN
      END
