c FILE driver_example3.f
c Example of driver for the subroutine frb_step
c InputData file: init_example3.dat
c Output data file: out_example3.dat


c     Initialization
C     .. Local Scalars ..
      DOUBLE PRECISION ANORMPI,H,HIJ,TFIN,TMPNP
      INTEGER I,J,K,L,N,NSTEPS,NP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AINERT(3),PIIN(3),PIOUT(3),QIN(3,3),QOUT(3,3)
C     ..
C     .. External Subroutines ..
      EXTERNAL FRB_STEP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT
C     ..
      OPEN (11,FILE='init_example3.dat')
      REWIND (11)

c     read timestep, integration time
      READ (11,FMT=*) H
      READ (11,FMT=*) TFIN

c     read inertia matrix, initial angular moment, initial attitude
      READ (11,FMT=*) (PIIN(I),I=1,3)
c     read Qin columnwise
      READ (11,FMT=*) ((QIN(I,J),I=1,3),J=1,3)
      READ (11,FMT=*) (AINERT(I),I=1,3)
c     read order of the method
      READ (11,FMT=*) TMPNP
      
      CLOSE (11)

      NP = INT(TMPNP)

c     re-adjust the timestep to fit into the time interval
      NSTEPS = INT(TFIN/H)
      H = TFIN/DBLE(NSTEPS)


c     open output file
      OPEN (11,FILE='out_example3.dat')
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
C      WRITE (6,FMT=*) 'Piin =', (PIIN(I),I=1,3)
c     display Qin _row_wise
C      WRITE (6,FMT=*) 'Qin =', ((QIN(I,J),J=1,3),I=1,3)
C      WRITE (6,FMT=*) 'h = ',H,'Tfin = ',TFIN

c     main computational routine
c     As we require output at the end of each step, we do not take 
c     into account techniques for reducing the computational efforts 
c     of the splitting method.

      DO N = 1,NSTEPS

          ANORMPI = SQRT(PIIN(1)**2+PIIN(2)**2+PIIN(3)**2)
          DO K = 1,3
              PIIN(K) = PIIN(K)/ANORMPI
          END DO

          HIJ = ANORMPI*H
          CALL FRB_STEP(PIOUT,QOUT,HIJ,AINERT,PIIN,QIN,NP)

          DO K = 1,3
              PIOUT(K) = ANORMPI*PIOUT(K)
              PIIN(K) = PIOUT(K)
              DO L = 1,3
                  QIN(K,L) = QOUT(K,L)
              END DO
          END DO

c     write the data in the output file
          DO I = 1,3
              WRITE (11,FMT=*) PIOUT(I)
          END DO
c     write Qout column_wise on output file
          DO J = 1,3
              DO I = 1,3
                  WRITE (11,FMT=*) QOUT(I,J)
              END DO
          END DO
          WRITE (11,FMT=*) N*H

      END DO

      CLOSE (11)

      END
