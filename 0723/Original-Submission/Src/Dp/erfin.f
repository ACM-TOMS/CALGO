      SUBROUTINE ERFIN
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1994-11-11 CLL Typing all variables.
C>> 1985-09-23 ERFIN  Lawson  Initial code.
C
      integer idelta, ialpha
      COMMON/M77ERR/IDELTA,IALPHA
      SAVE /M77ERR/
C
 1003 FORMAT(1X,72('$')/' ')
      PRINT 1003
      IF (IALPHA.GE.2) STOP
      RETURN
      END
