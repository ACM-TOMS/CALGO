      DOUBLE PRECISION FUNCTION DSECND( )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*     01-16-04:  Replace with F95 routine CPU_TIME (eca)
*
*  Purpose
*  =======
*
*  DSECND returns the user time for a process in seconds.
*
* =====================================================================
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Executable Statements ..
*
      CALL CPU_TIME( T1 )
      DSECND = T1
      RETURN
*
*     End of DSECND
*
      END
