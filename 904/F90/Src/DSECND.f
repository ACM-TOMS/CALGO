      DOUBLE PRECISION FUNCTION DSECND( )
*
*  -- Modified LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*  Purpose
*  =======
*
*  DSECND returns the user time for a process in seconds.
*  This version gets the time from the system function CLOCK().
*  No, sorry, this routine does not return anything but zero!
*
* =====================================================================
*
*      USE DFPORT
*
*     ..
*     .. Executable Statements ..
*
      DSECND = 0.0d+00
      RETURN
*
*     End of DSECND
*
      END
