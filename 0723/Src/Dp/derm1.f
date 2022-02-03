      SUBROUTINE DERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
C     .  Copyright (C) 1989, California Institute of Technology.
C     .  U. S. Government sponsorship under
C     .  NASA contract NAS7-918 is acknowledged.
C>> 1994-10-20 DERM1  Krogh  Changes to use M77CON
C>> 1994-04-20 DERM1  CLL Edited to make DP & SP files similar.
C>> 1985-08-02 DERM1  Lawson  Initial code.
c--D replaces "?": ?ERM1, ?ERV1
C
      CHARACTER*(*) SUBNAM,MSG,LABEL
      CHARACTER*1 FLAG
      integer INDIC, LEVEL
      DOUBLE PRECISION VALUE
C
      CALL ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      CALL DERV1(LABEL,VALUE,FLAG)
C
      RETURN
      END
