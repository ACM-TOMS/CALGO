*       test driver for DHAEVD.f
*     *****************************************************************
*
*     This program sets up several ad-hoc Hamiltonian matrices and then
*     calls DHAEVD to find their eigenvalues.  The program verifies
*     that the nonzero smallest singular values sigma_min(H - lambda I)
*     are can be accounted for by rounding errors.
*
*     ******************************************************************
C     .. Parameters ..
*        NOUT is the standard output, e.g., a terminal

      INTEGER NOUT
      PARAMETER (NOUT=6)
      INTEGER NMAX,NMIN
      PARAMETER (NMAX=15,NMIN=-1)
      INTEGER LD
      PARAMETER (LD=NMAX)
      INTEGER LDA,LDQG,LDH
      PARAMETER (LDA=LD,LDQG=LD,LDH=2*LD)
      INTEGER LWORK
      PARAMETER (LWORK=5*LDH)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      CHARACTER*(*) EQULS
      PARAMETER (EQULS=' =============================================')
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ANRM,GNRM,HNRM,QNRM,RATMAX,RTEPS
      INTEGER I,IERR,INFO,J,K,MULT,N,NFAILS,NUNSRT
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX HP(LDH,LDH),U(1),VT(1),ZWORK(LWORK)
      DOUBLE PRECISION A(LDA,LDA),H(LDH,LDH),QG(LDQG,LDQG+1),
     +                 RWORK(10*LD),S(LDH),WI(2*LD),WORK(LD,LD),WR(2*LD)
      INTEGER ISEED(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,DLANGE,DLANSY,DLAPY2
      EXTERNAL DLAMCH,DLANGE,DLANSY,DLAPY2
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,DHAEVD,DLACPY,DLARNV,DLASCL,ZGESVD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX,MAX,SQRT
C     ..
C     .. Data statements ..
*     . (seed for random number generation) .
*
      DATA ISEED/86,1967,2001,1995/
C     ..
*     ******************************************************************

      RTEPS = SQRT(DLAMCH('P'))
*
      NFAILS = 0
      NUNSRT = 0
      MULT = 0
*         ==== 2 test runs: with and without Hessenberg scaling =====
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*) EQULS
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*) EQULS

      DO 90 N = NMIN,NMAX
          WRITE (NOUT,FMT=*)
          WRITE (NOUT,FMT=*) EQULS
          WRITE (NOUT,FMT=*) 'DHAEVD TEST: ORDER N = ',N
          IF (N.LT.0) THEN
              WRITE (NOUT,FMT=*)
     +          'OF COURSE, THIS IS NOT A VALID TEST CASE.'
              WRITE (NOUT,FMT=*)
     +          'IT CHECKS THAT DHAEVD RETURNS WITH ERROR IERR = -1.'

          END IF
*
*             ==== set up Hamiltonian ====
          DO 10 I = 1,N
              CALL DLARNV(3,ISEED,N,A(1,I))
              CALL DLARNV(3,ISEED,N,QG(1,I))
   10     CONTINUE
          IF (N.GE.0) THEN
              CALL DLARNV(3,ISEED,N,QG(1,N+1))

              CALL DLACPY('A',N,N,A,LDA,H,LDH)
              CALL DLACPY('L',N,N,QG,LDQG,H(N+1,1),LDH)
              CALL DLACPY('U',N,N,QG(1,2),LDQG,H(1,N+1),LDH)
          END IF

          DO 20 J = 1,N
              CALL DCOPY(N,A(1,J),1,H(N+J,N+1),LDH)
              CALL DCOPY(J,QG(J,1),LDQG,H(N+1,J),1)
              CALL DCOPY(J,QG(1,J+1),1,H(J,N+1),LDH)
   20     CONTINUE
          IF (N.GT.0) CALL DLASCL('G',0,0,ONE,-ONE,N,N,H(N+1,N+1),LDH,
     +                            INFO)
*
*             ==== Frobenius norm of Hamiltonian Matrix ====
          ANRM = DLANGE('F',N,N,A,LD,RWORK)
          GNRM = DLANSY('F','L',N,QG,LDQG,RWORK)
          QNRM = DLANSY('F','U',N,QG(1,2),LDQG,RWORK)
          HNRM = DLAPY2(DLAPY2(ANRM,GNRM),DLAPY2(ANRM,QNRM))
*
*             ==== find eigenvalues ====
          CALL DHAEVD(N,A,LDA,QG,LDQG,WR,WI,WORK,IERR)
*
*             **********************************************************
*             ==== Test 0: Error exit from DHAEVD? =====

          IF (IERR.NE.0) WRITE (NOUT,FMT=*) 'DHAEVD RETURNS IERR = ',
     +        IERR
          IF (IERR.GT.0) WRITE (NOUT,FMT=*
     +        ) 'DHSEQR FAILED TO CONVERGE WHILE COMPUTING THE',J,
     +        'TH EIGENVALUE.  (THIS IS VERY RARE.)'
          IF (IERR.EQ.-1) WRITE (NOUT,FMT=*) 'N = ',N,' IS NEGATIVE.'
          IF (IERR.EQ.-3) WRITE (NOUT,FMT=*) 'LDA = ',LDA,
     +        ' IS LESS THAN N = ',N,'.'
          IF (IERR.EQ.-5) WRITE (NOUT,FMT=*) 'LDQG = ',LDQG,
     +        ' IS LESS THAN N = ',N,'.'

          IF (IERR.EQ.0) THEN
              WRITE (NOUT,FMT=*) 'NO ERRORS WERE DETECTED.'
*                 ******************************************************
*                 ==== Test 1: Sort check =====
              DO 30 I = 1,N - 1
                  IF ((WR(I).LT.WR(I+1)) .OR.
     +                ((WR(I).EQ.ZERO).AND. (WR(I+1).EQ.ZERO).AND.
     +                (WI(I).LT.WI(I+1)))) THEN
                      WRITE (NOUT,FMT=*) 'EIGENVALUE OUT OF ORDER!'
                      WRITE (NOUT,FMT=*) 'WR(',I,') = ',WR(I)
                      WRITE (NOUT,FMT=*) 'WR(',I + 1,') = ',WR(I+1)
                      NUNSRT = NUNSRT + 1
                      GO TO 50

                  END IF

   30         CONTINUE
              DO 40 I = 1,N - 1
                  IF (WR(I).EQ.WR(I+1) .AND. WR(I).NE.ZERO) THEN
                      IF (WI(I).GT.ZERO .AND. WI(I).NE.-WI(I+1)) THEN
                          WRITE (NOUT,FMT=*)
     +                      'COMPLEX CONJUGATES NOT ADJACENT OR ',
     +                      'OUT OF ORDER.'
                          WRITE (NOUT,FMT=*) 'WR(',I,') = ',WR(I),
     +                      '  WI(',I,') = ',WI(I)
                          WRITE (NOUT,FMT=*) 'WR(',I + 1,') = ',WR(I+1),
     +                      '  WI(',I + 1,') = ',WI(I+1)
                          NUNSRT = NUNSRT + 1
                          GO TO 50

                      END IF

                      IF (WR(I).EQ.WR(I+1) .AND. WI(I).EQ.WI(I+1)) THEN
                          WRITE (NOUT,FMT=*)
     +                      'PECULIAR! THERE ARE MULTIPLE EIGENVALUES.'
                          WRITE (NOUT,FMT=*) 'WR(',I,') = ',WR(I),
     +                      '  WI(',I,') = ',WI(I)
                          WRITE (NOUT,FMT=*) 'WR(',I + 1,') = ',WR(I+1),
     +                      '  WI(',I + 1,') = ',WI(I+1)

                          MULT = MULT + 1

                      END IF

                  END IF

   40         CONTINUE

   50         CONTINUE
*                 ==== Test 2: eigenvalue check ====

              RATMAX = ZERO
              DO 80 K = 1,N
                  DO 70 J = 1,2*N
                      DO 60 I = 1,2*N
                          HP(I,J) = H(I,J)
   60                 CONTINUE
                      HP(J,J) = HP(J,J) - DCMPLX(WR(K),WI(K))
   70             CONTINUE
                  CALL ZGESVD('N','N',2*N,2*N,HP,LDH,S,U,1,VT,1,ZWORK,
     +                        LWORK,RWORK,INFO)
                  RATMAX = MAX(RATMAX,S(2*N)/HNRM)
                  IF (S(2*N).GT.RTEPS*HNRM) THEN
                      WRITE (NOUT,FMT=*)
                      WRITE (NOUT,FMT=*) 'SUSPICIOUS RATIO'
                      WRITE (NOUT,FMT=*) 'COMPUTED EIGENVALUE = ',WR(K),
     +                  ' + ',WI(K),'I'
                      WRITE (NOUT,FMT=*) 'SIGMA-MIN(H-LAMBDA*I) = ',
     +                  S(2*N)
                      WRITE (NOUT,FMT=*)

                  END IF

   80         CONTINUE

              WRITE (NOUT,FMT=*) 'MAXIMUM RATIO =  ',RATMAX
              IF (RATMAX.GT.RTEPS) NFAILS = NFAILS + 1

          END IF

   90 CONTINUE

*     ******************************************************************
      WRITE (NOUT,FMT=*)
      WRITE (NOUT,FMT=*) EQULS
      WRITE (NOUT,FMT=*)

      IF (NFAILS.EQ.0 .AND. NUNSRT.EQ.0) WRITE (NOUT,
     +    FMT=*) 'NO SUSPICIOUS TEST RESULTS.'
      IF (MULT.NE.0) WRITE (NOUT,FMT=*) 'PECULIAR! THERE WERE ',MULT,
     +    ' MULTIPLE EIGENVALUES.'
      IF (NUNSRT.NE.0) WRITE (NOUT,FMT=*) 'SIGH..., THERE WERE ',NUNSRT,
     +    ' EIGENVALUE SORTING FAILURES.'
      IF (NFAILS.NE.0) WRITE (NOUT,FMT=*) 'SIGH..., THERE WERE ',NFAILS,
     +    ' SUSPICIOUS TEST RESULTS.'

      STOP

      END
