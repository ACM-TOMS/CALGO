      PROGRAM makdat
************************************************************************
*     (Make Data)
*     Read a *.out file, discarding comments, and write a new output
*     file containing pairs (x,f(x)) suitable for sorting and plotting.
*     This step must be done in quadruple-precision, so awk is not
*     available for the job.
*     [06-Jun-2000]
************************************************************************
*
*     Parameter variables
*
      REAL*16             two
      PARAMETER           (two = 2.0q+00)
*
      INCLUDE 'stdio.inc'
*
*     Local variables
*
      CHARACTER*256       line
      CHARACTER*80        fofx
*
      INTEGER             ignore,      p
*
      REAL*16             f
*
  100 READ (stdin, '(A)', END=200) line
      IF (line(1:1) .NE. '#') THEN
           READ (line,*) ignore, f, p, fofx
           WRITE (stdout,'(1X, 1P, E25.15E4, 1X, A)')
     X          f * two**p, fofx
      END IF
      GO TO 100
*
  200 CONTINUE
      END
