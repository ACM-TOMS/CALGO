*     DHAEVS EXAMPLE PROGRAM TEXT.
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
      PARAMETER (LRWORK=NMAX* (NMAX+1))
C     ..
C     .. Local Scalars ..
      INTEGER I,IERR,J,N
      CHARACTER JOBSCL
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX),QG(LDQG,NMAX+1),RWORK(LRWORK),
     +                 WI(NMAX),WR(NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL DHAEVS
C     ..
*
      WRITE (NOUT,FMT=9000)
*     Skip the heading in the data file and read the data.
*     NOTE: input must define a square-reduced Hamiltonian matrix.
      READ (NIN,FMT='()')
      READ (NIN,FMT=*) N,JOBSCL
      IF (N.LE.0 .OR. N.GT.NMAX) THEN
          WRITE (NOUT,FMT=9010) N

      ELSE
          READ (NIN,FMT=*) ((A(I,J),J=1,N),I=1,N)
          READ (NIN,FMT=*) ((QG(J,I+1),I=J,N),J=1,N)
          READ (NIN,FMT=*) ((QG(I,J),I=J,N),J=1,N)
*       Compute the eigenvalues.
          CALL DHAEVS(JOBSCL,N,A,LDA,QG,LDQG,WR,WI,RWORK,IERR)
*
          IF (IERR.NE.0) THEN
              WRITE (NOUT,FMT=9020) IERR

          ELSE
*       Show the computed eigenvalues.
              WRITE (NOUT,FMT=9030)
              DO 10 I = 1,N
                  WRITE (NOUT,FMT=9040) WR(I),' + (',WI(I),')i'
   10         CONTINUE
              DO 20 I = N,1,-1
                  WRITE (NOUT,FMT=9040) - WR(I),' + (',-WI(I),')i'
   20         CONTINUE
          END IF

      END IF

      STOP
*
 9000 FORMAT (' DHAEVS EXAMPLE PROGRAM RESULTS',/,1X)
 9010 FORMAT (/,' N is out of range.',/,' N = ',I5)
 9020 FORMAT (' IERR on exit from DHAEVS = ',I2)
 9030 FORMAT (/,' The eigenvalues are ')
 9040 FORMAT (1X,F8.4,A,F8.4,A)
      END
