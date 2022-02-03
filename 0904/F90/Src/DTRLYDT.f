CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRLYDT( JOB, OP, I, N, NB, A, LDA, C, LDC, WORK, 
     $     LWORK, ROWS, XRST, SCALE, INFO )
C
C  -- LAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     October 24, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          JOB, OP
      INTEGER            I, NB, INFO, LDA, LDC, N, ROWS, XRST, LWORK
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), C( * ), WORK( * )
C     ..
C     
      INTEGER UPLOSIGN
      LOGICAL LSAME
      EXTERNAL LSAME
C
      IF( LSAME( OP, 'N' ) ) THEN
         UPLOSIGN = 0
      ELSE
         UPLOSIGN = 1
      END IF
C
      CALL RECLYDT77( UPLOSIGN, SCALE, N, A, LDA, C, LDC, INFO )
      CALL DLATCPY( 'Upper', N-1, N-1, C(1+LDC), LDC,
     $     C(2), LDC )
C
      ROWS = N
      XRST = 1
C     
      END
C     
C     End of DTRLYDT
C     
C *** Last line of DTRLYDT ***
