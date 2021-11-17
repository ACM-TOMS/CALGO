c FILE driver_example2.f
c Example of driver for the subroutine quat_step. Single step
c of length h for the free rigid body.
c InputData file: init_example2.dat
c Output data file: out_example2.dat


c     Initialization
C     .. Local Scalars ..
      DOUBLE PRECISION ANORMPI,H,HIJ
      INTEGER I,K
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AINERT(3),PIIN(3),PIOUT(3),QIN(4),QOUT(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL QUAT_STEP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      OPEN (11,FILE='init_example2.dat')
      REWIND (11)

c     read timestep
      READ (11,FMT=*) H

c     read inertia matrix, initial angular moment, initial attitude
      READ (11,FMT=*) (AINERT(I),I=1,3)
      READ (11,FMT=*) (PIIN(I),I=1,3)
      READ (11,FMT=*) (QIN(I),I=1,4)
      CLOSE (11)

c     open output file
      OPEN (11,FILE='out_example2.dat')
      REWIND (11)

      DO I = 1,3
          WRITE (11,FMT=*) PIIN(I)
      END DO
      DO I = 1,4
          WRITE (11,FMT=*) QIN(I)
      END DO
      WRITE (11,FMT=*) 0.0D0*H

c     display the data
      WRITE (6,FMT=*) 'Initial data'
      WRITE (6,FMT=*) 'Piin = ', (PIIN(I),I=1,3)
      WRITE (6,FMT=*) 'Qin =', (QIN(I),I=1,4)
      WRITE (6,FMT=*) 'h = ',H

c     main computational routine

      ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
      DO K = 1,3
          PIIN(K) = PIIN(K)/ANORMPI
      END DO

      HIJ = ANORMPI*H
      CALL QUAT_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN,0)

      DO K = 1,3
          PIOUT(K) = ANORMPI*PIOUT(K)
      END DO

c     write the data in the output file
      DO I = 1,3
          WRITE (11,FMT=*) PIOUT(I)
      END DO
      DO I = 1,4
          WRITE (11,FMT=*) QOUT(I)
      END DO
      CLOSE (11)

c     display the data
      WRITE (6,FMT=*) 'Computational results:'
      WRITE (6,FMT=*) 'Piout =', (PIOUT(I),I=1,3)
      WRITE (6,FMT=*) 'Qout =', (QOUT(I),I=1,4)

      END
