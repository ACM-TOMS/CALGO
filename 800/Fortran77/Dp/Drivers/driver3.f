*       test driver for DHABL.f
*     *****************************************************************
*     This program scales several ad-hoc Hamiltonian matrices, then
*     verifies that the scaled matrices are similar to the originals
*     as they should be.
*     *****************************************************************
C     .. Parameters ..
*        NOUT is the standard output, e.g., a terminal
      INTEGER NOUT
      PARAMETER (NOUT=6)
      INTEGER NMAX
      PARAMETER (NMAX=15)
      INTEGER LDA,LDQG
      PARAMETER (LDA=NMAX,LDQG=NMAX)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CHARACTER*(*) EQULS
      PARAMETER (EQULS='=============================================')
C     ..
C     .. Local Scalars ..
*
      DOUBLE PRECISION ANRM,BASE,EPS,ERRA,ERRG,ERRQ,GNRM,GPNRM,QNRM,
     +                 QPNRM,R,SIMERR,TAU,Y
      INTEGER I,IERR,IFUDGE,INFO,IREP,J,N,NFAILS
      LOGICAL BBAL,FAILED,NORM,SYMP
      CHARACTER JOBSCL
C     ..
C     .. Local Arrays ..
*
      DOUBLE PRECISION A(LDA,LDA),AP(LDA,LDA),D(NMAX),QG(LDQG,LDQG+1),
     +                 QPGP(LDQG,LDQG+1),RWORK(NMAX)
      INTEGER ISEED(4)
C     ..
C     .. External Functions ..
*
      DOUBLE PRECISION DLAMCH,DLANGE,DLANSY,DLAPY2
      EXTERNAL DLAMCH,DLANGE,DLANSY,DLAPY2
C     ..
C     .. External Subroutines ..
*
      EXTERNAL DHABL,DLACPY,DLARNV,DLASCL,DRSCL,DSCAL
C     ..
C     .. Intrinsic Functions ..
*
      INTRINSIC DBLE,MAX
C     ..
C     .. Data statements ..
*     . (seed for random number generation) .
*
*     ******************************************************************
      DATA ISEED/86,1967,2001,1995/
C     ..
      EPS = DLAMCH('P')
      BASE = DLAMCH('B')
*
      WRITE (NOUT,FMT=*) EQULS
      WRITE (NOUT,FMT=*) EQULS
      WRITE (NOUT,FMT=*) 'TESTING SYMPLECTIC SCALING'

      NFAILS = 0

      DO 100 IREP = 1,2
          DO 90 N = -1,NMAX
              WRITE (NOUT,FMT=*)
              WRITE (NOUT,FMT=*) EQULS
              WRITE (NOUT,FMT=*) 'DHAEV TESTS: ORDER N = ',N
              IF (N.LT.0) THEN
                  WRITE (NOUT,FMT=*)
     +              'OF COURSE, THIS IS NOT A VALID TEST CASE.'
                  WRITE (NOUT,FMT=*)
     +              'IT CHECKS THAT DHAEV RETURNS WITH AN ERROR.'

              END IF
*
*             .. set up Hamiltonian ..
              DO 10 I = 1,N
                  CALL DLARNV(3,ISEED,N,A(1,I))
                  CALL DLARNV(3,ISEED,N,QG(1,I))
   10         CONTINUE
              DO 20 I = 1,N
                  CALL DRSCL(N,TWO** (-I),A(I,1),LDA)
                  CALL DSCAL(N,TWO** (-I),A(1,I),1)
   20         CONTINUE
              IF (N.GT.0) THEN
                  CALL DLARNV(3,ISEED,N,QG(1,N+1))
                  CALL DLASCL('U',0,0,DBLE(10*N),ONE,N,N,QG(1,2),LDQG,
     +                        INFO)
                  CALL DLASCL('L',0,0,ONE,DBLE(10*N),N,N,QG,LDQG,INFO)

              END IF

              CALL DLACPY('A',N,N,A,LDA,AP,LDA)
              CALL DLACPY('A',N,N+1,QG,LDQG,QPGP,LDQG)
*
*             .. symplectic scaling ..
              IF (IREP.EQ.1) THEN
                  WRITE (NOUT,FMT=*) 'SYMPLECTIC SCALING'
                  JOBSCL = 'A'
                  SYMP = .TRUE.
                  NORM = .FALSE.

              ELSE
                  WRITE (NOUT,FMT=*) 'NORM SCALING'
                  JOBSCL = 'B'
                  SYMP = .FALSE.
                  NORM = .TRUE.

              END IF

              CALL DHABL(JOBSCL,N,AP,LDA,QPGP,LDQG,D,RWORK,IERR)
*
*             **********************************************************
*             .. Has DHABL something to tell? =====
              IF (IERR.EQ.0) THEN
                  FAILED = .FALSE.
                  WRITE (NOUT,FMT=*) 'NO ERRORS DETECTED.'
                  IF (N.GT.0) TAU = D(1)

              ELSE
                  WRITE (NOUT,FMT=*) 'ERROR: DHABL RETURNS IERR = ',IERR
                  FAILED = .TRUE.

              END IF

              IF (.NOT.FAILED) THEN
*                 *****************************************************
*                 .. Nonsingularity check ..
                  DO 30 I = 1,N
                      IF (D(I).EQ.ZERO) THEN
                          FAILED = .TRUE.
                          WRITE (NOUT,FMT=*) 'ERROR IN DHABL: D(',I,
     +                      ') = 0'

                      END IF

                      IF (NORM) GO TO 40
   30             CONTINUE
   40             CONTINUE

              END IF

              IF (.NOT.FAILED .AND. N.GT.0) THEN
*                 *****************************************************
*                 .. block scaling check ..
*
                  IF (NORM) THEN
                      ANRM = DLANGE('1',N,N,A,LDA,RWORK)
                      GNRM = DLANSY('1','U',N,QG(1,2),LDQG,RWORK)
                      QNRM = DLANSY('1','L',N,QG,LDQG,RWORK)
                      Y = MAX(ONE,ANRM,GNRM,QNRM)
                      IF (Y.LE.TAU*BASE .AND. Y.GE.TAU/BASE) THEN
                          WRITE (NOUT,FMT=*) 'SCALE OK:'

                      ELSE
                          WRITE (NOUT,FMT=*)
     +                      'ERROR IN DHABL: OUT OF SCALE'
                          FAILED = .TRUE.

                      END IF

                      WRITE (NOUT,FMT=*) '             TAU   = ',TAU
                      WRITE (NOUT,FMT=*) '             ANRM  = ',ANRM
                      WRITE (NOUT,FMT=*) '             GNRM  = ',GNRM
                      WRITE (NOUT,FMT=*) '             QNRM  = ',QNRM

                  ELSE
                      GPNRM = DLANSY('1','U',N,QPGP(1,2),LDQG,RWORK)
                      QPNRM = DLANSY('1','L',N,QPGP,LDQG,RWORK)
                      BBAL = (GPNRM.LT.QPNRM*BASE) .AND.
     +                       (QPNRM/BASE.LT.GPNRM)
                      BBAL = BBAL .OR. (QPNRM.EQ.ZERO) .OR.
     +                       (GPNRM.EQ.ZERO)
                      IF (BBAL) THEN
                          WRITE (NOUT,FMT=*) 'SCALE OK:    GPNRM = ',
     +                      GPNRM
                          WRITE (NOUT,FMT=*) '             QPNRM = ',
     +                      QPNRM

                      ELSE
                          WRITE (NOUT,FMT=*)
     +                      'ERROR IN DHABL: OUT OF SCALE'
                          WRITE (NOUT,FMT=*) 'GP NORM = ',GPNRM
                          WRITE (NOUT,FMT=*) 'QP NORM = ',QPNRM
                          FAILED = .TRUE.

                      END IF

                  END IF

              END IF

              IF (.NOT.FAILED) THEN
*                 *****************************************************
*                 .. Similarity check (modulo rounding) ..
*
                  ERRA = ZERO
                  DO 60 J = 1,N
                      DO 50 I = 1,N
                          IF (NORM) R = A(I,J) - AP(I,J)*TAU
                          IF (SYMP) R = A(I,J) - AP(I,J)*D(I)/D(J)
                          IF (A(I,J).NE.ZERO) THEN
                              R = R/A(I,J)

                          ELSE
                              R = R/ANRM

                          END IF

                          ERRA = DLAPY2(ERRA,R)
   50                 CONTINUE
   60             CONTINUE

                  ERRQ = ZERO
                  ERRG = ZERO
                  DO 80 J = 1,N
                      DO 70 I = J,N
                          IF (NORM) R = QG(I,J) - QPGP(I,J)
                          IF (SYMP) R = QG(I,J) - QPGP(I,J)/D(I)/D(J)
                          IF (QG(I,J).NE.ZERO) THEN
                              R = R/QG(I,J)

                          ELSE
                              R = R/QNRM

                          END IF

                          ERRQ = DLAPY2(ERRQ,R)

                          IF (NORM) R = QG(J,I+1) - QPGP(J,I+1)*TAU*TAU
                          IF (SYMP) R = QG(J,I+1) -
     +                                  QPGP(J,I+1)*D(I)*D(J)
                          IF (QG(J,I+1).NE.ZERO) THEN
                              R = R/QG(J,I+1)

                          ELSE
                              R = R/GNRM

                          END IF

                          ERRG = DLAPY2(ERRG,R)
   70                 CONTINUE
   80             CONTINUE
                  SIMERR = DLAPY2(DLAPY2(ERRA,ERRG),DLAPY2(ERRQ,ERRA))
                  IFUDGE = 2*N
*
                  WRITE (NOUT,FMT=*)
                  WRITE (NOUT,FMT=*) 'SIMILARITY CHECK:'
                  WRITE (NOUT,FMT=*) 'RELATIVE ERROR  = ',SIMERR
*
                  IF (SIMERR.LE.DBLE(IFUDGE)*EPS) THEN

                      WRITE (NOUT,FMT=*)
     +                  'OK: THE OUTPUT APPEARS TO BE SIMILAR TO THE'
                      WRITE (NOUT,FMT=*)
     +                  '    INPUT ''MODULO ROUNDING ERRORS'' AS IT '
                      WRITE (NOUT,FMT=*) '    SHOULD BE.'

                  ELSE
                      FAILED = .TRUE.
                      WRITE (NOUT,FMT=*)
     +                  'SUSPICIOUS: THE RELATIVE ERROR OF SIMILARITY'
                      WRITE (NOUT,FMT=*)
     +                  '            IS LARGER THAN EXPECTED. '

                  END IF

              END IF

              IF (FAILED .AND. N.GE.0) NFAILS = NFAILS + 1
              WRITE (NOUT,FMT=*)

   90     CONTINUE
          WRITE (NOUT,FMT=*)
          WRITE (NOUT,FMT=*) EQULS
          WRITE (NOUT,FMT=*) EQULS
          WRITE (NOUT,FMT=*)
  100 CONTINUE

*     ******************************************************************
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*) EQULS
      IF (NFAILS.EQ.0) THEN
          WRITE (NOUT,FMT=*) 'NO SUSPICIOUS TEST RESULTS.'

      ELSE
          WRITE (NOUT,FMT=*) 'SIGH..., ',NFAILS,
     +      ' SUSPICIOUS TEST RESULTS.'

      END IF

      STOP

      END
