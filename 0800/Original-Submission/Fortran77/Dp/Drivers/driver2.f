*     DHABL EXAMPLE PROGRAM TEXT.
*
C     .. Parameters ..
*        NIN  is the standard input.
*        NOUT is the standard output, e.g., a terminal
      INTEGER NIN,NOUT
      PARAMETER (NIN=5,NOUT=6)
      INTEGER NMAX
      PARAMETER (NMAX=20)
      INTEGER LDA,LDQG
      PARAMETER (LDA=NMAX,LDQG=NMAX)
      INTEGER LRWORK
      PARAMETER (LRWORK=NMAX)
C     ..
C     .. Local Scalars ..
      INTEGER I,IERR,J,N
      CHARACTER JOBSCL
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX),D(NMAX),QG(LDQG,NMAX+1),RWORK(LRWORK)
C     ..
C     .. External Subroutines ..
      EXTERNAL DHABL
C     ..
C     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*
      WRITE (NOUT,FMT=9000)
*     Skip the heading in the data file and read the data.
      READ (NIN,FMT='()')
      READ (NIN,FMT=*) N,JOBSCL
      IF (N.LE.0 .OR. N.GT.NMAX) THEN
          WRITE (NOUT,FMT=9010) N

      ELSE
          READ (NIN,FMT=*) ((A(I,J),J=1,N),I=1,N)
          READ (NIN,FMT=*) ((QG(J,I+1),I=J,N),J=1,N)
          READ (NIN,FMT=*) ((QG(I,J),I=J,N),J=1,N)
      END IF
*     .. scale the Hamiltonian matrix ..
      CALL DHABL(JOBSCL,N,A,LDA,QG,LDQG,D,RWORK,IERR)
      IF (IERR.NE.0) THEN
          WRITE (NOUT,FMT=9020) IERR

      ELSE
*     .. show the scaled Hamiltonian matrix..
          WRITE (NOUT,FMT=9030)
          DO 10 I = 1,N
              WRITE (NOUT,FMT=9060) (A(I,J),J=1,N),
     +          (QG(J,I+1),J=1,I-1), (QG(I,J+1),J=I,N)
   10     CONTINUE
          DO 20 I = 1,N
              WRITE (NOUT,FMT=9060) (QG(I,J),J=1,I-1), (QG(J,I),J=I,N),
     +          (-A(J,I),J=1,N)
   20     CONTINUE
*     .. show the scaling factors..
          IF (LSAME(JOBSCL,'A')) THEN
              WRITE (NOUT,FMT=9040)
              WRITE (NOUT,FMT=9060) (D(I),I=1,N)

          ELSE IF (LSAME(JOBSCL,'B')) THEN
              WRITE (NOUT,FMT=9050)
              WRITE (NOUT,FMT=9060) D(1)
          END IF

      END IF

      STOP
*
 9000 FORMAT (' DHABL EXAMPLE PROGRAM RESULTS',/,1X)
 9010 FORMAT (/,' N is out of range.',/,' N = ',I5)
 9020 FORMAT (' IERR on exit from DHABL = ',I2)
 9030 FORMAT (/,' The scaled Hamiltonian is ')
 9040 FORMAT (/,' The scaling factors are ')
 9050 FORMAT (/,' The scaling factor tau is ')
 9060 FORMAT (1X,8 (F10.4))
      END
