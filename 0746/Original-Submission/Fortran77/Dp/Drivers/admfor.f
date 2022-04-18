C
C   GENERATION OF FORTRAN CODES FOR AUTOMATIC DIFFERENTIATION
C
      INTEGER LRSYM,LISYM
      PARAMETER (LRSYM=15000,LISYM=15000)
C
      DOUBLE PRECISION RSYM(LRSYM)
      INTEGER ISYM(LISYM)
      INTEGER LMRSYM,LMISYM,LARSYM,LAISYM,IERR,DEBFIL
      INTEGER LROW,LEVDEB,NVAR,NFUNC
C
      LROW=0
      DEBFIL=0
      LEVDEB=0
C
C   OPEN FILES
C
      OPEN(2,FILE='ADM.FUN',STATUS='UNKNOWN')
      OPEN(3,FILE='ADM.SYM',STATUS='UNKNOWN')
      OPEN(4,FILE='adm_gen.f',STATUS='UNKNOWN')
C
C   PARSE FUNCTION INPUT
C
      LMRSYM=LRSYM
      LMISYM=LISYM
      CALL SYMINP(2,3,RSYM,LRSYM,ISYM,LISYM,LARSYM,LAISYM,IERR,LROW,1,
     /            NVAR,NFUNC,DEBFIL,LEVDEB)
C      CALL SYMPRP(3,RSYM,LRSYM,ISYM,LISYM,LARSYM,LAISYM,IERR,1,NV,NF)
      IF (IERR.GT.0) GOTO 830
C
C   GENERATE FORTRAN CODE FOR FUNCTIONS AND GRADIENTS
C
      CALL SYMFOR(4,RSYM,LMRSYM,ISYM,LMISYM,IERR)
      IF (IERR.GT.0) GOTO 830
      GOTO 9999
C
C   ERROR
C
  830 CALL SYMERR(LROW,IERR)
C
C   CLOSE FILES
C
 9999 CONTINUE
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
C
C   END OF MAIN PROGRAMM
C
      STOP
      END
C