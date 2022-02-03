*     The minimal DLAMCH required from RECSYCT and xLACON.
*     This is written by the author, because we could not compile
*     the distributed version of DLAMCH. Please do not use this code
*     for the other purposes.
*
      FUNCTION DLAMCH(CMACH)
      CHARACTER CMACH
      REAL*8    DLAMCH, ONE
      ONE = 1
      DLAMCH = 0
      IF      (CMACH .EQ. 'P') then
        DLAMCH = EPSILON(ONE)
        RETURN
      ELSE IF (CMACH .EQ. 'S') then
        DLAMCH = TINY(ONE)
        RETURN
      ELSE
        WRITE(*,*) 'THIS IS NOT THE FULL VERSION OF DLAMCH OF LAPACK.'
        STOP 999
      END IF
      RETURN
      END FUNCTION DLAMCH
