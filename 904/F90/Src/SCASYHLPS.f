CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     This file contains two helper functions that casts integers to
C     logicals and vice versa.
C
      FUNCTION ILG2NT( ILOG )
      IMPLICIT NONE
      LOGICAL ILOG
      INTEGER ILG2NT
      IF( ILOG ) THEN
         ILG2NT = 1
      ELSE
         ILG2NT = 0
      END IF
      RETURN
      END
C
      FUNCTION INT2LG( INT )
      IMPLICIT NONE
      INTEGER INT
      LOGICAL INT2LG
      IF( INT.GT.0 ) THEN
         INT2LG = .TRUE.
      ELSE
         INT2LG = .FALSE.
      END IF
      RETURN
      END
C     
C     End of SCASYHLPS
C     
C *** Last line of SCASYHLPS ***
