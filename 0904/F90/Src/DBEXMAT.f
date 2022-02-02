CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DBEXMAT( EXROW, EXCOL, ROWS, COLS, LIA, LJA, A, 
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
      INTEGER              ROWS, COLS, LIA, LJA, LLDA, NBC, MB, NB, 
     $                     LLDEXT
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION     A( * ), EXA( * ) , EXTMAT( * )
C     ..
C     
C     Purpose and description
C     =======================
C     To simplify solving with and updating with the extended 
C     matrices we do build them in this subroutine.
C
C     .. Local Scalars ..
      INTEGER I, LBI, LBJ, POS0, POS, SKIP, LENG
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
C     Compute some local indicies
C
      IF( EXROW .OR. EXCOL ) THEN
         LBI = ICEIL( LIA, MB )
         LBJ = ICEIL( LJA, NB )
         POS0 = ((LBI-1)*NBC + (LBJ-1))*(MB+NB+1) + 1
      END IF
C
C     First copy over the part of the possible extended matrix that
C     is stored in A. After that go for the extension elements:
C     Then check if we shall add any elements from EXA to the matrix
C     If so, then we add them. Here we must take into account the 
C     possibility that the first element of the extension row or
C     column shall not be addressed, because of the extensions of
C     neighbouring blocks. We do this in the following way:
C
C     if we have ( LIA mod MB ) /= 1, and / or
C     if we have ( LJA mod NB ) /= 1
C
C     we address the extension column and the extension row from
C     their second element, respectively. In fact the value of the
C     modulo-expressions above can just have two values, namely 1 or 2.
C         
C
      IF( EXROW.AND.EXCOL ) THEN
         CALL DLACPY( 'All', ROWS-1, COLS-1, A( (LJA-1) * LLDA + LIA ),
     $               LLDA, EXTMAT, LLDEXT )
C
         IF( MOD( LJA-1, NB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS = POS0 + SKIP 
         CALL DLACPY('All', 1, COLS-1, EXA( POS ), 1, EXTMAT( ROWS ),
     $               LLDEXT)
C 
         POS = POS0 + NB 
         EXTMAT( (COLS-1)*LLDEXT + ROWS ) = EXA( POS )
C
         IF( MOD( LIA-1, MB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS =  POS0 + NB + 1 + SKIP
         CALL DLACPY('All', ROWS-1, 1, EXA( POS ), MB+NB+1,
     $               EXTMAT( (COLS-1)*LLDEXT + 1 ), LLDEXT)
      ELSEIF( EXROW ) THEN
         CALL DLACPY( 'All', ROWS-1, COLS, A( (LJA-1) * LLDA + LIA ),
     $                LLDA, EXTMAT, LLDEXT )
C
         IF( MOD( LJA-1, NB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS = POS0 + SKIP 
         CALL DLACPY('All', 1, COLS, EXA( POS ), 1, EXTMAT( ROWS ),
     $               LLDEXT)
      ELSEIF( EXCOL ) THEN
         CALL DLACPY( 'All', ROWS, COLS-1, A( (LJA-1) * LLDA + LIA ),
     $                LLDA, EXTMAT, LLDEXT )
C
         IF( MOD( LIA-1, MB ).NE.0 ) THEN
            SKIP = 1
         ELSE
            SKIP = 0
         END IF
         POS =  POS0 + NB + 1 + SKIP
         CALL DLACPY('All', ROWS, 1, EXA( POS ), MB+NB+1,
     $               EXTMAT( (COLS-1)*LLDEXT + 1 ), LLDEXT)
      ELSE
         CALL DLACPY( 'All', ROWS, COLS, A( (LJA-1) * LLDA + LIA ),
     $                LLDA, EXTMAT, LLDEXT )
      END IF
C     
      END
C
C     End of DBEXMAT
C
C *** Last line of DBEXMAT ***
