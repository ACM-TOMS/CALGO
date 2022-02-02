CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DLAPRNT( M, N, A, IA, JA, LDA, CMATNM, NOUT )
C
C  -- ScaLAPACK-style routine (preliminary version ) --
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     January 18, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER            IA, JA, M, N, NOUT, LDA
C     ..
C     .. Array Arguments ..
      CHARACTER*(*)      CMATNM
      DOUBLE PRECISION   A( * )
C     ..
C
C  Purpose
C  =======
C
C  DLAPRNT prints to the standard output a matrix sub( A )
C  denoting A(IA:IA+M-1,JA:JA+N-1). 
C
C  Arguments
C  =========
C
C  M       INTEGER
C          The number of rows to be operated on i.e the number of rows
C          of the submatrix sub( A ). M >= 0.
C
C  N       INTEGER
C          The number of columns to be operated on i.e the number of
C          columns of the submatrix sub( A ). N >= 0.
C
C  A       DOUBLE PRECISION pointer into an array of dimension (LDA,N) 
C          containing the matrix sub( A ).
C
C  IA      INTEGER
C          The row index indicating the first row of sub( A ).
C
C  JA      INTEGER
C          The column index indicating the first column of sub( A ).
C
C  LDA     INTEGER 
C          The leading dimension of the matrix sub( A ).
C
C  CMATNM  CHARACTER*(*)
C          Identifier of the matrix to be printed.
C
C  NOUT    INTEGER
C          The unit number for output file. NOUT = 6, ouput to screen,
C          NOUT = 0, output to stderr.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     .. Local Scalars ..
      INTEGER            I, J, INDEX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
C
      DO 10 I = IA, IA+M-1
         DO 20 J = JA, JA+N-1
            INDEX = (J-1)*LDA + I
            WRITE( NOUT, FMT = 9999 )
     $           CMATNM, I, J, A( INDEX )
 20      CONTINUE
 10   CONTINUE
C     
 9999 FORMAT(A,'(',I6,',',I6,')=',D30.18,';')
C     
      RETURN
C     
      END
C     
C     End of PTRGCSYD
C     
C *** Last line of PTRGCSYD ***
