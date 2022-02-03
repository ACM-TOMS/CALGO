CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DUBEXMA( EXROW, EXCOL, ROWS, COLS, IA, JA, A, 
     $                    LLDA, NBC, MB, NB, EXA, EXTMAT, LLDEXT )
C
C  -- ScaLAPACK-style tool routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 21, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar arguments ..
      LOGICAL              EXROW, EXCOL
      INTEGER              ROWS, COLS, IA, JA, LLDA, NBC, MB, NB, LLDEXT
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION     A( * ), EXA( * ) , EXTMAT( * )
C     ..
C     
C     Purpose and description
C     =======================
C     To simplify solving with, updating with and communication of 
C     the extended matrices we do unbuild them in this subroutine.
C
C     .. Local Scalars ..
      INTEGER I, LBI, LBJ, POS0, POS, SKIP
C     ..
C     .. External Subroutines ..
      EXTERNAL           DLACPY
C     ..
C     .. External Functions ..
      INTEGER ICEIL
      EXTERNAL ICEIL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Executable Statements ..
C
C     Compute some local block indicies
C
      IF( EXROW.OR.EXCOL ) THEN
         LBI = ICEIL( IA, MB )
         LBJ = ICEIL( JA, NB )
         POS0 =  ((LBI-1)*NBC+(LBJ-1)) * (MB+NB+1) + 1
      END IF
C
C     Choose branch depending on how the matrix is extended.
C     Here we must take into account the possibility that the first
C     element of the extension row or column shall not be addressed,
C     because of the extensions of neighbouring blocks. We do this in
C     the following way:
C
C     if we have ( IA mod MB ) /= 1, and / or
C     if we have ( JA mod NB ) /= 1
C
C     we address the extension column and the extension row from
C     their second element, respectively. In fact the value of the
C     modulo-expressions above can just have two values, namely 1 or 2.
C 
      IF( EXROW.AND.EXCOL ) THEN
         CALL DLACPY('All', ROWS-1, COLS-1, EXTMAT, LLDEXT,
     $               A((JA-1)*LLDA + IA), LLDA )
         IF( MOD( JA-1, NB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS = POS0 + SKIP
         CALL DLACPY('All', 1, COLS-1, EXTMAT(ROWS), LLDEXT,
     $               EXA( POS ), 1 )
C
         POS = POS0 + NB
         EXA( POS ) = EXTMAT( (COLS-1)*LLDEXT + ROWS)
C
         IF( MOD( IA-1, MB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS =  POS0 + NB + 1 + SKIP
         CALL DLACPY('All', ROWS-1, 1, EXTMAT((COLS-1)*LLDEXT+1),
     $               LLDEXT, EXA( POS ), MB + NB + 1 )    
      ELSEIF( EXROW ) THEN
         CALL DLACPY('All', ROWS-1, COLS, EXTMAT, LLDEXT,
     $               A((JA-1)*LLDA + IA), LLDA )
         IF( MOD( JA-1, NB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS = POS0 + SKIP
         CALL DLACPY('All', 1, COLS, EXTMAT(ROWS), LLDEXT,
     $               EXA( POS ), 1 )
      ELSEIF( EXCOL ) THEN
         CALL DLACPY('All', ROWS, COLS-1, EXTMAT, LLDEXT,
     $               A((JA-1)*LLDA + IA), LLDA )
         IF( MOD( IA-1, MB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS =  POS0 + NB + 1 + SKIP
         CALL DLACPY('All', ROWS, 1, EXTMAT((COLS-1)*LLDEXT+1),
     $               LLDEXT, EXA( POS ), MB + NB + 1 )    
      ELSE
         CALL DLACPY('All', ROWS, COLS, EXTMAT, LLDEXT,
     $               A((JA-1)*LLDA + IA), LLDA ) 
      END IF
C
      END
C
C     End of DUBEXMA
C
C *** Last line of DUBEXMA ***
