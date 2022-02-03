*     DHASRD EXAMPLE PROGRAM TEXT.
*
C     .. Parameters ..
*        NIN  is the standard input.
*        NOUT is the standard output, e.g., a terminal
      INTEGER NIN,NOUT
      PARAMETER (NIN=5,NOUT=6)
      INTEGER NMAX
      PARAMETER (NMAX=20)
      INTEGER LDA,LDQG,LDU
      PARAMETER (LDA=NMAX,LDQG=NMAX,LDU=NMAX)
      INTEGER LRWORK
      PARAMETER (LRWORK= (NMAX+NMAX)* (NMAX+NMAX+1))
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Local Scalars ..
      INTEGER I,IERR,IJ,J,JI,N,POS,WPOS
      CHARACTER COMPU
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX),QG(LDQG,NMAX+1),RWORK(LRWORK),
     +                 U(LDU,NMAX)
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DGEMM,DHASRD,DSYMV
C     ..
*
      WRITE (NOUT,FMT=9000)
*     Skip the heading in the data file and read the data.
      READ (NIN,FMT='()')
      READ (NIN,FMT=*) N,COMPU
      IF (N.LE.0 .OR. N.GT.NMAX) THEN
          WRITE (NOUT,FMT=9010) N

      ELSE
          READ (NIN,FMT=*) ((A(I,J),J=1,N),I=1,N)
          READ (NIN,FMT=*) ((QG(J,I+1),I=J,N),J=1,N)
          READ (NIN,FMT=*) ((QG(I,J),I=J,N),J=1,N)
      END IF
*     .. square-reduce by symplectic orthogonal similarity ..
      CALL DHASRD(COMPU,N,A,LDA,QG,LDQG,U,LDU,RWORK,IERR)
      IF (IERR.NE.0) THEN
          WRITE (NOUT,FMT=9020) IERR

      ELSE
*     .. show the square-reduced Hamiltonian ..
          WRITE (NOUT,FMT=9030)
          DO 10 I = 1,N
              WRITE (NOUT,FMT=9050) (A(I,J),J=1,N),
     +          (QG(J,I+1),J=1,I-1), (QG(I,J+1),J=I,N)
   10     CONTINUE
          DO 20 I = 1,N
              WRITE (NOUT,FMT=9050) (QG(I,J),J=1,I-1), (QG(J,I),J=I,N),
     +          (-A(J,I),J=1,N)
   20     CONTINUE
*     .. show the square of H ..
          WRITE (NOUT,FMT=9040)
          WPOS = (NMAX+NMAX)* (NMAX+NMAX)
*                                                 T
*     .. compute N11 = A*A + G*Q and set N22 = N11  ..
          CALL DGEMM('N','N',N,N,N,ONE,A,LDA,A,LDA,ZERO,RWORK,N+N)
          DO 40 I = 1,N
              CALL DCOPY(N-I+1,QG(I,I),1,RWORK(WPOS+I),1)
              DO 30 J = 1,I - 1
                  RWORK(WPOS+J) = QG(I,J)
   30         CONTINUE
              CALL DSYMV('U',N,ONE,QG(1,2),LDQG,RWORK(WPOS+1),1,ONE,
     +                   RWORK((I-1)* (N+N)+1),1)
              POS = N* (N+N) + N + I
              CALL DCOPY(N,RWORK((I-1)* (N+N)+1),1,RWORK(POS),N+N)
   40     CONTINUE
          DO 50 I = 1,N
              CALL DSYMV('U',N,-ONE,QG(1,2),LDQG,A(I,1),LDA,ZERO,
     +                   RWORK((N+I-1)* (N+N)+1),1)
              CALL DSYMV('L',N,ONE,QG,LDQG,A(1,I),1,ZERO,
     +                   RWORK((I-1)* (N+N)+N+1),1)
   50     CONTINUE
          DO 70 J = 1,N
              DO 60 I = J,N
                  IJ = (N+J-1)* (N+N) + I
                  JI = (N+I-1)* (N+N) + J
                  RWORK(IJ) = RWORK(IJ) - RWORK(JI)
                  RWORK(JI) = -RWORK(IJ)
                  IJ = N + I + (J-1)* (N+N)
                  JI = N + J + (I-1)* (N+N)
                  RWORK(IJ) = RWORK(IJ) - RWORK(JI)
                  RWORK(JI) = -RWORK(IJ)
   60         CONTINUE
   70     CONTINUE
          DO 80 I = 1,N + N
              WRITE (NOUT,FMT=9050) (RWORK(I+ (J-1)* (N+N)),J=1,N+N)
   80     CONTINUE
      END IF

      STOP
*
 9000 FORMAT (' DHASRD EXAMPLE PROGRAM RESULTS',/,1X)
 9010 FORMAT (/,' N is out of range.',/,' N = ',I5)
 9020 FORMAT (' IERR on exit from DHASRD = ',I2)
 9030 FORMAT (/,' The square-reduced Hamiltonian is ')
 9040 FORMAT (/,' The square of the square-reduced Hamiltonian is ')
 9050 FORMAT (1X,8 (F10.4))
      END
