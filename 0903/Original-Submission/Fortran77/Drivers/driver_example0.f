c FILE driver_example0.f
c Example of driver for the subroutine frb_step. Single step of length h
c for the free rigid body.
c InputData file: init_example0.dat
c Output data file: out_example0.dat


c     Initialization
C     .. Local Scalars ..
      DOUBLE PRECISION ANORMPI,H,HIJ
      INTEGER I,J,K
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AINERT(3),PIIN(3),PIOUT(3),QIN(3,3),QOUT(3,3)
C     ..
C     .. External Subroutines ..
      EXTERNAL FRB_STEP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      OPEN (11,FILE='init_example0.dat')
      REWIND (11)

c     read timestep
      READ (11,FMT=*) H

c     read inertia matrix, initial angular moment, initial attitude
      READ (11,FMT=*) (AINERT(I),I=1,3)
      READ (11,FMT=*) (PIIN(I),I=1,3)
c     OBS: read Qin columnwise
      READ (11,FMT=*) ((QIN(I,J),I=1,3),J=1,3)
      CLOSE (11)

c     open output file
      OPEN (11,FILE='out_example0.dat')
      REWIND (11)

      DO I = 1,3
          WRITE (11,FMT=*) PIIN(I)
      END DO
c     OBS: write Qout columnwise
      DO J = 1,3
          DO I = 1,3
              WRITE (11,FMT=*) QIN(I,J)
          END DO
      END DO
      WRITE (11,FMT=*) 0.0D0*H

c     display the data
C      WRITE (6,FMT=*) 'Initial data'
C      WRITE (6,FMT=*) 'Piin = ', (PIIN(I),I=1,3)
C      WRITE (6,FMT=*) 'Qin =', ((QIN(I,J),J=1,3),I=1,3)
C      WRITE (6,FMT=*) 'h = ',H

c     main computational routine

      ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
      DO K = 1,3
          PIIN(K) = PIIN(K)/ANORMPI
      END DO

      HIJ = ANORMPI*H
      CALL FRB_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN, 0)

      DO K = 1,3
          PIOUT(K) = ANORMPI*PIOUT(K)
      END DO

c     write the data in the output file
      DO I = 1,3
          WRITE (11,FMT=*) PIOUT(I)
      END DO
c     OBS: write Qout columnwise
      DO J = 1,3
          DO I = 1,3
              WRITE (11,FMT=*) QOUT(I,J)
          END DO
      END DO
      CLOSE (11)

c     display the data
C      WRITE (6,FMT=*) 'Computational results:'
C      WRITE (6,FMT=*) 'Piout =', (PIOUT(I),I=1,3)
C      WRITE (6,FMT=*) 'Qout =', ((QOUT(I,J),J=1,3),I=1,3)

      END
