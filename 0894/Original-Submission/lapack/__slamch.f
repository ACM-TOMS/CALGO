*     The minimal SLAMCH required from RECSYCT and xLACON.
*     This is written by the author, because we could not compile
*     the distributed version of SLAMCH. Please do not use this code
*     for the other purposes.
*
      FUNCTION SLAMCH(CMACH)
      CHARACTER CMACH
      REAL      SLAMCH, ONE
      ONE = 1
      SLAMCH = 0
      IF      (CMACH .EQ. 'P') then
        SLAMCH = EPSILON(ONE)
        RETURN
      ELSE IF (CMACH .EQ. 'S') then
        SLAMCH = TINY(ONE)
        RETURN
      ELSE
        WRITE(*,*) 'THIS IS NOT THE FULL VERSION OF SLAMCH OF LAPACK.'
        STOP 999
      END IF
      RETURN
      END FUNCTION SLAMCH
