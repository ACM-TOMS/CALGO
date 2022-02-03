c FILE driver_example1.f
c Example of driver for the subroutine frb_step
c InputData file: init_example1.dat
c Output data file: out_example1.dat


c     Initialization
C     .. Local Scalars ..
      DOUBLE PRECISION ANORMPI,H,HIJ,TFIN
      INTEGER I,J,K,L,N,NSTAGES,NSTEPS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AINERT(3),AN6(8),BN6(7),PIIN(3),PIOUT(3),
     +                 QIN(3,3),QOUT(3,3),U0(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL FRB_STEP,POTENTIAL_STEP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT,INT
C     ..
      OPEN (11,FILE='init_example1.dat')
      REWIND (11)


C     read coefficients splitting method
C     RK-Nystrom, SRKN6^a_14 6th order coeffs, 14  stages (Bl&Moan)
C     Ef =0.63
      READ (11,FMT=*) (AN6(I),I=1,8)
      READ (11,FMT=*) (BN6(I),I=1,7)
      NSTAGES = 7

c     read timestep, integration time
      READ (11,FMT=*) H
      READ (11,FMT=*) TFIN

c     read inertia matrix, initial angular moment, initial attitude
      READ (11,FMT=*) (PIIN(I),I=1,3)
c     read Qin columnwise
      READ (11,FMT=*) ((QIN(I,J),I=1,3),J=1,3)
      READ (11,FMT=*) (AINERT(I),I=1,3)
      READ (11,FMT=*) (U0(I),I=1,3)
      CLOSE (11)

c     re-adjust the timestep to fit into the time interval
      NSTEPS = INT(TFIN/H)
      H = TFIN/DBLE(NSTEPS)


c     open output file
      OPEN (11,FILE='out_example1.dat')
      REWIND (11)

      DO I = 1,3
          WRITE (11,FMT=*) PIIN(I)
      END DO
c     write Qin columnwise
      DO J = 1,3
          DO I = 1,3
              WRITE (11,FMT=*) QIN(I,J)
          END DO
      END DO
      WRITE (11,FMT=*) 0.0D0*H

c     display the data
      WRITE (6,FMT=*) 'Piin =', (PIIN(I),I=1,3)
c     display Qin _row_wise
      WRITE (6,FMT=*) 'Qin =', ((QIN(I,J),J=1,3),I=1,3)
      WRITE (6,FMT=*) 'h = ',H,'Tfin = ',TFIN
      WRITE (6,FMT=*) 'u0 =', (U0(I),I=1,3)

c     main computational routine
c     As we require output at the end of each step, we do 
c     not take into account techniques for reducing the computational
c     efforts of the splitting method.

      DO N = 1,NSTEPS

c     first internal stage
          ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
          DO K = 1,3
              PIIN(K) = PIIN(K)/ANORMPI
          END DO

          HIJ = AN6(1)*ANORMPI*H
          CALL FRB_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN,0)

          DO K = 1,3
              PIOUT(K) = ANORMPI*PIOUT(K)
              PIIN(K) = PIOUT(K)
              DO L = 1,3
                  QIN(K,L) = QOUT(K,L)
              END DO
          END DO

c     end first internal stage

c     other stages (forward)
          DO J = 1,NSTAGES
              HIJ = BN6(J)*H
              CALL POTENTIAL_STEP(PIOUT,QOUT,HIJ,U0,PIIN,QIN)
              DO K = 1,3
                  PIIN(K) = PIOUT(K)
                  DO L = 1,3
                      QIN(K,L) = QOUT(K,L)
                  END DO
              END DO

              ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
              DO K = 1,3
                  PIIN(K) = PIIN(K)/ANORMPI
              END DO
              HIJ = AN6(J+1)*ANORMPI*H
              CALL FRB_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN,0)

              DO K = 1,3
                  PIOUT(K) = ANORMPI*PIOUT(K)
                  PIIN(K) = PIOUT(K)
                  DO L = 1,3
                      QIN(K,L) = QOUT(K,L)
                  END DO
              END DO
          END DO

c     other stages (backward)
          DO J = NSTAGES,1,-1
              HIJ = BN6(J)*H
              CALL POTENTIAL_STEP(PIOUT,QOUT,HIJ,U0,PIIN,QIN)
              DO K = 1,3
                  PIIN(K) = PIOUT(K)
                  DO L = 1,3
                      QIN(K,L) = QOUT(K,L)
                  END DO
              END DO

              ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
              DO K = 1,3
                  PIIN(K) = PIIN(K)/ANORMPI
              END DO
              HIJ = AN6(J)*ANORMPI*H
              CALL FRB_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN,0)
              DO K = 1,3
                  PIOUT(K) = ANORMPI*PIOUT(K)
                  PIIN(K) = PIOUT(K)
                  DO L = 1,3
                      QIN(K,L) = QOUT(K,L)
                  END DO
              END DO

          END DO

c     write the data in the output file
          DO I = 1,3
              WRITE (11,FMT=*) PIOUT(I)
          END DO
c     write Qout columwise on output file
          DO J = 1,3
              DO I = 1,3
                  WRITE (11,FMT=*) QOUT(I,J)
              END DO
          END DO
          WRITE (11,FMT=*) N*H

      END DO

      CLOSE (11)

      END

c     ----- end driver program -------

      SUBROUTINE POTENTIAL_STEP(PIOUT,QOUT,H,U0,PIIN,QIN)


c     Compute u = transpose(Q)*u0
C     .. Scalar Arguments ..
      DOUBLE PRECISION H
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PIIN(3),PIOUT(3),QIN(3,3),QOUT(3,3),U0(3)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION U(3)
C     ..
      DO I = 1,3
          U(I) = 0.0D0
          DO J = 1,3
              U(I) = U(I) + QIN(J,I)*U0(J)
          END DO
      END DO

      PIOUT(1) = PIIN(1) + H*U(2)
      PIOUT(2) = PIIN(2) - H*U(1)
      PIOUT(3) = PIIN(3)

      DO I = 1,3
          DO J = 1,3
              QOUT(I,J) = QIN(I,J)
          END DO
      END DO

      END
