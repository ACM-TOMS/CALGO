      PROGRAM MAIN
C  ****************************************************************************
C                                                                             *
C  IMPORTANT NOTES                                                            *
C  ***************                                                            *
C                                                                             *
C  1. IWORK doesn't have to be as big as it is set here unless you are        *
C     running at seriously tight tolerances.                                  *
C                                                                             *
C  2. CPU times are measured using the routine DTIME, which is available      *
C     on SUN and Silicon Graphics machines (although it doesn't always        *
C     seem to do what the documentation says it does). For other machines,    *
C     DTIME will have to be replaced.                                         *
C                                                                             *
C  3. SLEUTH itself does not need the array WORKS to be as big as it is.      *
C     It needs an array of length 1136. The big size of WORKS is due to       *
C     the computation of eigenfunctions later on, in the call to SL4EFN.      *
C     SL4EFN requires 27 DOUBLE PRECISION storage allocations for each        *
C     meshpoint used, plus a further 56 DOUBLE PRECISION storage allocations. *
C     We could have written SL4EFN to require very few storage allocations    *
C     but then the eigenfunction computation process would either have been   *
C     more expensive or else less flexible in the way that the user can       *
C     specify the evaluation points.                                          *
C                                                                             *
C                                                                             *
C  END IMPORTANT NOTES                                                        *
C  ****************************************************************************
C     .. Parameters ..
C
      INTEGER IWORK
      PARAMETER (IWORK=40960)
      INTEGER IWORKS
      PARAMETER (IWORKS=27*IWORK+56)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ALFA,BETA,M,RPARM
      INTEGER IPROB
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,ELAM,ERR,TOL,X
      REAL ELAPSE
      INTEGER IFAIL,IK,IMATCH,K,MULTI,NMESH,NXTRAP
      LOGICAL FSTPAS,HOT,SYM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AEND(10),BEND(10),TARRAY(2),WORK(0:IWORK,1:6),
     &                 WORKS(IWORKS),Y1(1:4),Y2(1:4),YSTOR1(4,0:IWORK),
     &                 YSTOR2(4,0:IWORK)
      INTEGER KINDEX(10,10),NUMEIG(10)
      LOGICAL SYMARR(10)
C     ..
C     .. External Functions ..
      REAL DTIME
      EXTERNAL DTIME
C     ..
C     .. External Subroutines ..
      EXTERNAL SL4BCS,SL4COF,SL4EFN,SLEUTH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
C     .. Common blocks ..
      COMMON /COMBLK/IPROB
      COMMON /PARM/ALFA,RPARM,M,BETA
C     ..

C This code runs the test-problems described in the paper
C `The Code SLEUTH for solving Fourth Order Sturm-Liouville
C  Problems', by Greenberg and Marletta. All the input is
C hardwired (see below). The output for each test is as
C follows:
C
C ** The number of the problem attempted
C ** The index of the eigenvalue sought
C ** The error flag
C ** The eigenvalue approximation and error assessment
C ** The CPU time required
C ** The number of Richardson extrapolations used
C ** The number of mesh intervals used
C ** The value of the eigenfunction at the centre of the
C    interval [A,B] on which the problem is posed.
C

C Arrays of information about the problems (endpoints,
C number of eigenvalues to be computed and their
C indices):

C For problem 1:
      AEND(1) = -2.D0*ATAN(1.D0)
      BEND(1) = -AEND(1)
      SYMARR(1) = .TRUE.
      NUMEIG(1) = 3
      KINDEX(1,1) = 2
      KINDEX(1,2) = 3
      KINDEX(1,3) = 4
C For problem 2:
      AEND(2) = 0.D0
      BEND(2) = 4.D0*ATAN(1.D0)
      SYMARR(2) = .FALSE.
      NUMEIG(2) = 3
      KINDEX(2,1) = 0
      KINDEX(2,2) = 9
      KINDEX(2,3) = 99
C For problem 3:
      AEND(3) = -100.D0
      BEND(3) = 100.D0
      SYMARR(3) = .TRUE.
      NUMEIG(3) = 3
      KINDEX(3,1) = 0
      KINDEX(3,2) = 9
      KINDEX(3,3) = 99
C For problem 4:
      AEND(4) = 0.D0
      BEND(4) = 100.D0
      SYMARR(4) = .FALSE.
      NUMEIG(4) = 3
      KINDEX(4,1) = 0
      KINDEX(4,2) = 9
      KINDEX(4,3) = 19
C For problem 5:
      AEND(5) = -0.999999D0
      BEND(5) = 0.999999D0
      SYMARR(5) = .TRUE.
      NUMEIG(5) = 5
      KINDEX(5,1) = 0
      KINDEX(5,2) = 1
      KINDEX(5,3) = 2
      KINDEX(5,4) = 3
      KINDEX(5,5) = 9
C For problem 6:
      AEND(6) = 1.0D-7
      BEND(6) = 1.D0
      SYMARR(6) = .FALSE.
      NUMEIG(6) = 5
      KINDEX(6,1) = 0
      KINDEX(6,2) = 1
      KINDEX(6,3) = 2
      KINDEX(6,4) = 3
      KINDEX(6,5) = 9
C For problem 7:
      AEND(7) = -0.999999D0
      BEND(7) = 0.999999D0
      SYMARR(7) = .TRUE.
      NUMEIG(7) = 5
      KINDEX(7,1) = 0
      KINDEX(7,2) = 1
      KINDEX(7,3) = 2
      KINDEX(7,4) = 3
      KINDEX(7,5) = 9
C For problem 8:
      AEND(8) = 1.0D-6
      BEND(8) = 30.D0
      SYMARR(8) = .FALSE.
      NUMEIG(8) = 5
      KINDEX(8,1) = 0
      KINDEX(8,2) = 1
      KINDEX(8,3) = 2
      KINDEX(8,4) = 3
      KINDEX(8,5) = 9
C For problem 9:
      AEND(9) = 1.0D-6
      BEND(9) = 0.999999D0
      SYMARR(9) = .FALSE.
      NUMEIG(9) = 5
      KINDEX(9,1) = 0
      KINDEX(9,2) = 1
      KINDEX(9,3) = 2
      KINDEX(9,4) = 3
      KINDEX(9,5) = 9
C For problem 10:
      AEND(10) = -0.999999D0
      BEND(10) = 0.999999D0
      SYMARR(10) = .TRUE.
      NUMEIG(10) = 5
      KINDEX(10,1) = 0
      KINDEX(10,2) = 1
      KINDEX(10,3) = 2
      KINDEX(10,4) = 3
      KINDEX(10,5) = 9

C Tolerance:
      TOL = 1.0D-6

      FSTPAS = .TRUE.

      DO 30 IPROB = 1,10
          WRITE (6,FMT=*)
          WRITE (7,FMT=*)
          WRITE (6,FMT=*) 'Problem Number:',IPROB
          WRITE (6,FMT=*) '======================================'
          WRITE (7,FMT=*) 'Problem Number:',IPROB
          WRITE (7,FMT=*) '======================================'
   10     IF (IPROB.EQ.1) THEN
C Problem 1 is done with two different values for the parameter
C BETA:
              IF (FSTPAS) THEN
                  BETA = 10.D0
              ELSE
                  BETA = 30.D0
              END IF
          END IF

          IF (IPROB.EQ.7) THEN
              ALFA = 1.D0
          END IF

          IF (IPROB.EQ.8) THEN
              RPARM = 1.D0
          END IF

          IF (IPROB.EQ.9) THEN
              ALFA = 1.D0
              M = 1.D0
          END IF

C For each problem, compute the eigenvalues with the indices
C K specified:
          DO 20 IK = 1,NUMEIG(IPROB)
              K = KINDEX(IPROB,IK)
              WRITE (6,FMT=*) 'Eigenvalue index:',K
              WRITE (6,FMT=*) '-------------------------------------'
              WRITE (7,FMT=*) 'Eigenvalue index:',K
              WRITE (7,FMT=*) '-------------------------------------'
              A = AEND(IPROB)
              B = BEND(IPROB)
              SYM = SYMARR(IPROB)

              ELAPSE = DTIME(TARRAY)
C
C Pretty random-looking initial guess for the eigenvalue;
C did not use ELAM = 0 since some of the problems have
C lambda = 0 as an eigenvalue, which would make life too
C easy for the code.
C
              ELAM = 0.123456D0
              ERR = 1.D-1
              IFAIL = 0

              CALL SLEUTH(A,B,ELAM,ERR,K,SYM,SL4COF,SL4BCS,TOL,NMESH,
     &                    IMATCH,NXTRAP,WORK,IWORK,WORKS,IFAIL)

              ELAPSE = DTIME(TARRAY) - ELAPSE

              IF (IFAIL.EQ.0.D0) THEN
                  WRITE (6,FMT=*) 'Successful exit from SLEUTH'
                  WRITE (6,FMT=*) 'Eigenvalue approximation:',ELAM
                  WRITE (6,FMT=*) 'Error estimate:',ERR
                  WRITE (6,FMT=*) 'Number of extrapolations:',NXTRAP
                  WRITE (6,FMT=*) 'CPU Time:',ELAPSE
                  WRITE (6,FMT=*) 'Number of meshpoints:',NMESH
                  WRITE (6,FMT=*) 'Matchpoint:',IMATCH

                  WRITE (7,FMT=*) 'Successful exit from SLEUTH'
                  WRITE (7,FMT=*) 'Eigenvalue approximation:',ELAM
                  WRITE (7,FMT=*) 'Error estimate:',ERR
                  WRITE (7,FMT=*) 'Number of extrapolations:',NXTRAP
                  WRITE (7,FMT=*) 'CPU Time:',ELAPSE
                  WRITE (7,FMT=*) 'Number of meshpoints:',NMESH
                  WRITE (7,FMT=*) 'Matchpoint:',IMATCH

              ELSE
                  WRITE (6,FMT=*) 'Failure exit from SLEUTH, IFAIL=',
     &              IFAIL
                  WRITE (6,FMT=*) 'Eigenvalue approximation:',ELAM
                  WRITE (6,FMT=*) 'Error estimate:',ERR
                  WRITE (6,FMT=*) 'Number of extrapolations:',NXTRAP
                  WRITE (6,FMT=*) 'CPU Time:',ELAPSE

                  WRITE (7,FMT=*) 'Failure exit from SLEUTH, IFAIL=',
     &              IFAIL
                  WRITE (7,FMT=*) 'Eigenvalue approximation:',ELAM
                  WRITE (7,FMT=*) 'Error estimate:',ERR
                  WRITE (7,FMT=*) 'Number of extrapolations:',NXTRAP
                  WRITE (7,FMT=*) 'CPU Time:',ELAPSE
              END IF

              IF (IFAIL.EQ.0) THEN
C Compute the eigenfunction at the midpoint:
                  HOT = .FALSE.
                  X = (WORK(0,1)+WORK(NMESH,1))*0.5D0
                  CALL SL4EFN(X,Y1,Y2,HOT,WORK,IWORK,NMESH,IMATCH,ELAM,
     &                        MULTI,WORKS,YSTOR1,YSTOR2,SL4BCS,IFAIL)

                  IF (IFAIL.NE.0) THEN
                      WRITE (6,FMT=*)
     &                  'Failure exit from SL4EFN, IFAIL =',IFAIL
                      WRITE (7,FMT=*)
     &                  'Failure exit from SL4EFN, IFAIL =',IFAIL
                  END IF
                  WRITE (6,FMT=*) 'Multiplicity was',MULTI
                  WRITE (7,FMT=*) 'Multiplicity was',MULTI
                  WRITE (6,FMT=*) X,Y1(1)
                  WRITE (7,FMT=*) X,Y1(1)
              END IF
   20     CONTINUE

          IF (IPROB.EQ.1 .AND. FSTPAS) THEN
              FSTPAS = .FALSE.
              GO TO 10
          END IF

   30 CONTINUE

      STOP

      END

C -------------------------------------------------------------------------

      SUBROUTINE SL4COF(X,P,S,Q,W)

C     .. Scalars in Common ..
      DOUBLE PRECISION ALFA,BETA,M,RPARM
      INTEGER IPROB
C     ..
C     .. Common blocks ..

      COMMON /COMBLK/IPROB
      COMMON /PARM/ALFA,RPARM,M,BETA
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION P,Q,S,W,X
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,SIN
C     ..
      GO TO (10,20,30,40,50,60,70,80,90,
     &       100,110,120,130) IPROB

C Coffey Evans squared
   10 W = 1.D0
      P = 1.D0
      S = 2.D0* ((BETA*SIN(2.D0*X))**2) - 4.D0*BETA*COS(2.D0*X)
      Q = ((BETA**2)* (SIN(2.D0*X)**2)-2.D0*BETA*COS(2.D0*X))**2 -
     &    8.D0* (BETA**2)* ((COS(2.D0*X)**2)- (SIN(2.D0*X)**2)) -
     &    8.D0*BETA*COS(2.D0*X)
      RETURN

C Made-up-problem (oscillatory q(x)) squared
   20 P = 1.D0
      W = 1.D0
      S = 2.D0*COS(X) + 4.D0*COS(2.D0*X) + 6.D0*COS(3.D0*X)
      Q = (COS(X)+2.D0*COS(2.D0*X)+3.D0*COS(3.D0*X))**2 + COS(X) +
     &    8.D0*COS(2.D0*X) + 27.D0*COS(3.D0*X)
      RETURN

C Modified harmonic oscillator squared
   30 P = 1.D0
      W = 1.D0
      S = 2.D0* ((X**2)+ (X**4))
      Q = ((X**2)+ (X**4))**2 - 2.D0 - 12.D0* (X**2)
      RETURN

C Harmonic oscillator squared
   40 P = 1.D0
      W = 1.D0
      S = 2.D0* (X**2)
      Q = X**4 - 2
      RETURN

C Legendre equation squared
   50 P = (1-X**2.D0)**2.D0
      W = 1.D0
      Q = 0.D0
      S = 2.D0* (1-X**2.D0)
      RETURN

C Bessel equation squared
   60 P = X
      S = 9.D0/X
      Q = 0.D0
      W = X
      RETURN

C Littlejohn Problem 1
   70 W = 1.D0
      P = 1.D0 - X*X
      P = P*P
      S = 4.D0*ALFA* (1.D0-X*X) + 8.D0
      Q = 0.D0
      RETURN

C Littlejohn Problem 2
   80 W = EXP(-X)
      P = X*X*W
      S = ((2.D0*RPARM+2.D0)*X+2.D0)*W
      Q = 0.D0
      RETURN

C Littlejohn Problem 3
   90 P = X*X* ((1.D0-X)** (2.D0+ALFA))
      S = 2.D0* ((1.D0-X)** (1.D0+ALFA))* ((ALFA+M+1.D0)*X+1.D0)
      Q = 0.D0
      W = (1.D0-X)**ALFA
      RETURN

C Littlejohn Problem 1 with ALFA = 0
  100 W = 1.D0
      P = 1.D0 - X*X
      P = P*P
      S = 8.D0
      Q = 0.D0
      RETURN

C Additional Problem 1
  110 P = 1.D0
      S = X + X
      Q = X*X
      W = 1.D0
      RETURN

C Additional Problem 2
  120 P = 1.D0
      S = X + X
      Q = X*X
      W = 1.D0
      RETURN

C Additional Problem 3
  130 P = 1.D0
      S = 0.D0
      Q = X*X
      W = 1.D0
      RETURN

      END


      SUBROUTINE SL4BCS(IEND,ISING,X,U,V,ELAM)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ELAM,X
      INTEGER IEND
      LOGICAL ISING
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION U(2,2),V(2,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ALFA,BETA,M,RPARM
      INTEGER IPROB
C     ..
C     .. Common blocks ..

      COMMON /COMBLK/IPROB
      COMMON /PARM/ALFA,RPARM,M,BETA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,SIN
C     ..
      GO TO (10,20,30,40,50,60,70,80,90,
     &       100,110,120,130) IPROB

   10 U(1,1) = 0.D0
      U(1,2) = 0.D0
      U(2,1) = 1.D0
      U(2,2) = 0.D0
      V(1,1) = 1.D0
      V(1,2) = 1.D0
      V(2,1) = 0.D0
      V(2,2) = 0.D0
      RETURN
   20 IF (IEND.EQ.0) THEN
          U(1,1) = 0.D0
          U(1,2) = 0.D0
          U(2,1) = 1.D0
          U(2,2) = 0.D0
          V(1,1) = 1.D0
          V(1,2) = 1.D0
          V(2,1) = 0.D0
          V(2,2) = 0.D0
      ELSE
          ALPHA = SIN(X) + 4.D0*SIN(2.D0*X) + 9.D0*SIN(3.D0*X)
          U(1,1) = 0.D0
          U(1,2) = 1.D0
          U(2,1) = 0.D0
          U(2,2) = 0.D0
          V(1,1) = 0.D0
          V(1,2) = ALPHA
          V(2,1) = 1.D0
          V(2,2) = 0.D0
      END IF
      RETURN
   30 U(1,1) = 0.D0
      U(1,2) = 0.D0
      U(2,1) = 1.D0
      U(2,2) = 0.D0
      V(1,1) = 1.D0
      V(1,2) = 1.D0
      V(2,1) = 0.D0
      V(2,2) = 0.D0
      RETURN

C40        U(1,1) = 0.D0
C          U(1,2) = 0.D0
C          U(2,1) = 1.D0
C          U(2,2) = 0.D0
C          V(1,1) = 1.D0
C          V(1,2) = 1.D0
C          V(2,1) = 0.D0
C          V(2,2) = 0.D0
C      RETURN

   40 U(1,1) = 1.D0
      U(1,2) = 0.D0
      U(2,1) = 0.D0
      U(2,2) = 1.D0
      V(1,1) = 0.D0
      V(1,2) = 0.D0
      V(2,1) = 0.D0
      V(2,2) = 0.D0
      RETURN

   50 U(1,1) = 1.D0
      U(1,2) = 0.D0
      U(2,1) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = 0.D0
      V(1,2) = 0.D0
      V(2,1) = 0.D0
      V(2,2) = 1.D0
      RETURN
   60 IF (IEND.EQ.0) THEN
          ALPHA = (1.D0/X)* (2.D0+ELAM*X*X)/ (1.D0+ (ELAM/4.D0)*X**2.D0)
C          WRITE (10,FMT=*) 'EP = ',X,' AND ALPHA = ',ALPHA
C          ALPHA = 2.D0/X
          U(1,1) = 1.D0
          U(1,2) = 0.D0
          U(2,1) = ALPHA
          U(2,2) = 0.D0
          V(1,1) = 0.D0
          V(1,2) = 1.D0
          V(2,1) = (1.D0/ALPHA)* ((8.D0/X)+8.D0*ALPHA- (ALPHA**2.D0))
          V(2,2) = -1.D0/ALPHA
      ELSE
          U(1,1) = 0.D0
          U(1,2) = 0.D0
          U(2,1) = 1.D0
          U(2,2) = 0.D0
          V(1,1) = 0.D0
          V(1,2) = 1.D0
          V(2,1) = -1.D0
          V(2,2) = 0.D0
      END IF
      RETURN

C Lambda-dependent BC for Littlejohn Problem 1
   70 U(1,1) = 1.D0
      U(2,1) = ELAM/ (8.D0*ALFA*X)
      U(1,2) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = ELAM/ (ALFA*X)
      V(2,1) = 0.D0
      V(1,2) = -ELAM/ (8.D0*ALFA*X)
      V(2,2) = 1.D0
      RETURN

C Lambda-dependent BC for Littlejohn Problem 2
   80 IF (IEND.EQ.0) THEN
C          U(1,1) = 1.D0
C          U(2,1) = -0.5D0*ELAM/RPARM
C          U(1,2) = 0.D0
C          U(2,2) = 0.D0
C          V(1,1) = -ELAM/RPARM
C          V(2,1) = 0.D0
C          V(1,2) = 0.5D0*ELAM/RPARM
C          V(2,2) = 1.D0
          U(1,1) = 1.D0
          U(2,1) = 0.D0
          U(1,2) = -3.D0 + 3.D0*X*X/10.D0
          U(2,2) = 3.D0*X/10.D0
          V(1,1) = 0.D0
          V(2,1) = 0.D0
          V(1,2) = 3.D0*X*X*EXP(-X)
          V(2,2) = 3.D0*X*X*EXP(-X)

      ELSE
          U(1,1) = 1.D0
          U(1,2) = 0.D0
          U(2,1) = 0.D0
          U(2,2) = 1.D0
          V(1,1) = 0.D0
          V(2,1) = 0.D0
          V(1,2) = 0.D0
          V(2,2) = 0.D0
      END IF
      RETURN

C Lambda-dependent BC for Littlejohn Problem 3
   90 IF (IEND.EQ.0) THEN
          U(1,1) = 1.D0
          U(2,1) = -0.5D0*ELAM/M
          U(1,2) = 0.D0
          U(2,2) = 0.D0

          V(1,1) = -ELAM/M
          V(2,1) = 0.D0
          V(1,2) = 0.5D0*ELAM/M
          V(2,2) = 1.D0

      ELSE
          U(1,1) = 1.D0
          U(2,1) = 0.D0
          U(1,2) = 0.D0
          U(2,2) = 1.D0

          V(1,1) = 0.D0
          V(2,1) = 0.D0
          V(1,2) = 0.D0
          V(2,2) = 0.D0
      END IF

      RETURN

C Lambda-independent BC for Littlejohn Problem 1 with ALFA = 0.
  100 U(1,1) = 0.D0
      U(2,1) = 1.D0
      U(1,2) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = 0.D0
      V(2,1) = 0.D0
      V(1,2) = 1.D0
      V(2,2) = 0.D0

      RETURN

C Additional problem 1
  110 U(1,1) = 1.D0
      U(1,2) = 0.D0
      U(2,1) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = -1.D0
      V(2,1) = 0.D0
      V(1,2) = 0.D0
      V(2,2) = 1.D0
      RETURN

C Additional problem 2
  120 U(1,1) = 0.D0
      U(1,2) = 0.D0
      U(2,1) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = 1.D0
      V(2,1) = 0.D0
      V(1,2) = 0.D0
      V(2,2) = 1.D0
      RETURN

C Additional problem 3
  130 U(1,1) = 0.D0
      U(1,2) = 0.D0
      U(2,1) = 0.D0
      U(2,2) = 0.D0
      V(1,1) = 1.D0
      V(2,1) = 0.D0
      V(1,2) = 0.D0
      V(2,2) = 1.D0
      RETURN

      END
