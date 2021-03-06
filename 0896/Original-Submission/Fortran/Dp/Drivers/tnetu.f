!     TEST PROGRAM FOR THE SUBROUTINE PNETU
!
      INTEGER NF,IPAR(7),IHES,IPRNT,ITERM
      DOUBLE PRECISION X(1000),RPAR(9),F,GMAX
      INTEGER NEXT,IERR,I,ITIME
      INTEGER NITER,NFVAL,NGVAL,NSUCC
      COMMON /PROB/ NEXT
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      NITER=0
      NFVAL=0
      NGVAL=0
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
        IHES=0
        IPRNT=1
!
!     PROBLEM DIMENSION
!
        NF=1000
!
!     INITIATION OF X AND CHOICE OF RPAR(1) AND RPAR(6)
!
        CALL TIUD14 (NF, X, RPAR(6), RPAR(1), NEXT, IERR)
        IF (IERR.NE.0) GO TO 30
        IF (RPAR(6).EQ.0.0D0) IPAR(4)=1
        IF (NEXT.EQ.2) IPAR(4)=0
        RPAR(1)=0.0D0
        IF (NEXT.EQ. 2) RPAR(1)=1.0D1
        IF (NEXT.EQ.10) RPAR(1)=3.5D1
        IF (NEXT.EQ.18) RPAR(1)=1.0D1
!
!     SOLUTION
!
        CALL PNETU (NF, X, IPAR, RPAR, F, GMAX, IHES, IPRNT, ITERM)
        NITER=NITER+NIT
        NFVAL=NFVAL+NFV
        NGVAL=NGVAL+NFG
        IF (ITERM.GT.0.AND.ITERM.LT.9) NSUCC=NSUCC+1
   30 CONTINUE
      WRITE (6,40) NITER,NFVAL,NGVAL,NSUCC
   40 FORMAT (' NITER =',I5,3X,' NFVAL =',I5,3X,' NGVAL =',I5,3X,' NSUCC
     & =',I5)
      CALL TYTIM2 (ITIME)
      STOP
      END
!     USER SUPPLIED SUBROUTINE (CALCULATION OF FF)
!
      SUBROUTINE OBJ (NF, X, FF)
      INTEGER NF
      DOUBLE PRECISION X(*),FF
      INTEGER NEXT
      COMMON /PROB/ NEXT
!
!     FUNCTION EVALUATION
!
      CALL TFFU14 (NF, X, FF, NEXT)
      RETURN
      END
!     USER SUPPLIED SUBROUTINE (CALCULATION OF GF)
!
      SUBROUTINE DOBJ (NF, X, GF)
      INTEGER NF
      DOUBLE PRECISION X(*),GF(*)
      INTEGER NEXT
      COMMON /PROB/ NEXT
!
!     GRADIENT EVALUATION
!
      CALL TFGU14 (NF, X, GF, NEXT)
      RETURN
      END
!     EMPTY SUBROUTINES
!
      SUBROUTINE HVEC (NF, X, D, HD)
      INTEGER NF
      DOUBLE PRECISION X(*),D(*),HD(*)
      NF=1
      D(1)=X(1)
      HD(1)=X(1)
      RETURN
      END
      SUBROUTINE FUN (NF, KA, X, FA)
      INTEGER NF,KA
      DOUBLE PRECISION X(*),FA
      KA=NF
      FA=X(1)
      RETURN
      END
      SUBROUTINE DFUN (NF, KA, X, GA)
      INTEGER NF,KA
      DOUBLE PRECISION X(*),GA(*)
      KA=NF
      GA(1)=X(1)
      RETURN
      END
