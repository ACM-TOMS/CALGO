*
*     Provide symbolic names for three standard (in UNIX, MS DOS, and
*     MS Windows) I/O units.  The values 0, 5, and 6 are common in
*     many, but NOT all, UNIX systems.  The extra value, stddat, is
*     available for use as an additional I/O unit.
*
      INTEGER             stddat
      PARAMETER           (stddat = 1)
      INTEGER             stderr
      PARAMETER           (stderr = 0)
      INTEGER             stdin
      PARAMETER           (stdin = 5)
      INTEGER             stdout
      PARAMETER           (stdout = 6)
