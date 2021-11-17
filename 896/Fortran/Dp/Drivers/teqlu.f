!     TEST PROGRAM FOR THE SUBROUTINE PEQLU
!
      DOUBLE PRECISION F,FMIN,GMAX
      INTEGER I,IERR,IDER,ISPAS,IPRNT,ITERM,ITIME,M,MM,N
      INTEGER NITER,NFVAL,NSUCC,NITCG
      DOUBLE PRECISION AF(5000),RPAR(9),X(5000)
      INTEGER IAG(5001),IPAR(7),JAG(100000)
      INTEGER NEXT
      COMMON /PROB/ NEXT
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      NITER=0
      NFVAL=0
      NSUCC=0
      NITCG=0
      CALL TYTIM1 (ITIME)
!
!     LOOP FOR 30 TEST PROBLEMS
!
      DO 30 NEXT=1,30
!
!     CHOICE OF INTEGER AND REAL PARAMETERS
!
        DO 10 I=1,7
          IPAR(I)=0
   10   CONTINUE
        DO 20 I=1,9
          RPAR(I)=0.0D0
   20   CONTINUE
        IDER=0
        ISPAS=2
        IPRNT=1
!
!     PROBLEM DIMENSION
!
        N=3000
        MM=10000
!
!     INITIATION OF X, DETERMINATION IAG AND JAG AND CHOICE OF RPAR(1)
!
        CALL TIUB18 (N, N, MM, X, IAG, JAG, FMIN, RPAR(1), NEXT, IERR)
        IF (NEXT.EQ.5) RPAR(8)=1.0D-3
        IF (NEXT.EQ.29) RPAR(8)=1.0D-2
        IF (NEXT.EQ.30) RPAR(8)=1.0D-2
        RPAR(1)=0.0D0
        IF (NEXT.EQ.5) RPAR(1)=1.0D4
        IF (NEXT.EQ.7) RPAR(1)=1.0D2
        IF (NEXT.EQ.9) RPAR(1)=1.0D1
        IF (NEXT.EQ.18) RPAR(1)=1.0D0
        IF (NEXT.EQ.22) RPAR(1)=1.0D4
        IF (IERR.NE.0) GO TO 30
!
!     SOLUTION
!
        CALL PEQLU (N, M, X, AF, IAG, JAG, IPAR, RPAR, F, GMAX, IDER,
     &   ISPAS, IPRNT, ITERM)
        NITER=NITER+NIT
        NFVAL=NFVAL+NFV
        NITCG=NITCG+NIN
        IF (ITERM.EQ.3) NSUCC=NSUCC+1
   30 CONTINUE
      WRITE (6,40) NITER,NFVAL,NITCG,NSUCC
   40 FORMAT (' NITER =',I5,3X,' NFVAL =',I5,3X,' NITCG =',I5,3X,' NSUCC
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
      CALL TAFU18 (NF, KA, X, FA, NEXT)
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
      CALL TAGU18 (NF, KA, X, GA, NEXT)
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
