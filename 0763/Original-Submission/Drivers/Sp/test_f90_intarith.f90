PROGRAM TEST_INTERVAL_ARITHMETIC

!  This routine tests the Fortran 90 interface to the elementary
!  interval arithmetic portion of INTLIB.

USE INTERVAL_ARITHMETIC

      TYPE(INTERVAL) A, B, C
      DOUBLE PRECISION D
      INTEGER N

      CALL SIMINI

      A = INTERVAL(1,2)
      B = INTERVAL(3,4)
      N = 3

      OPEN(6,FILE='TEST_F90_INTARITH.OUT')

      C = A+B;          WRITE(6,'(2(1X,ES12.3E2))') C
      C = A**N;         WRITE(6,'(2(1X,ES12.3E2))') C
      C = A**B;         WRITE(6,'(2(1X,ES12.3E2))') C
      C = COS(A);       WRITE(6,'(2(1X,ES12.3E2))') C
      C = A.CH.B;       WRITE(6,'(2(1X,ES12.3E2))') C
      D = COS(A%LOWER); WRITE(6,'((1X,ES12.3E2))') D
      D = COS(A%UPPER); WRITE(6,'((1X,ES12.3E2))') D

      WRITE(6,'(2(1X,L1))') A.SB.C, A.SB.B

      STOP

END PROGRAM TEST_INTERVAL_ARITHMETIC
