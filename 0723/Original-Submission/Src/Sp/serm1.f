      SUBROUTINE SERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1994-10-20 SERM1  Krogh  Changes to use M77CON
C>> 1994-04-20 SERM1  CLL Edited to make DP & SP files similar.
C>> 1985-08-02 SERM1  Lawson  Initial code.
c--S replaces "?": ?ERM1, ?ERV1
C
      CHARACTER*(*) SUBNAM,MSG,LABEL
      CHARACTER*1 FLAG
      integer INDIC, LEVEL
      REAL             VALUE
C
      CALL ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      CALL SERV1(LABEL,VALUE,FLAG)
C
      RETURN
      END
