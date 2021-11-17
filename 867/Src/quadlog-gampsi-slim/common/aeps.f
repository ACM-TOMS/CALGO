      REAL FUNCTION aeps(x)
************************************************************************
*     (Single-precision machine epsilon)
*     Return the smallest positive number such that (x + aeps(x))
*     differs from x.
*     (27-Dec-1999)
************************************************************************
*
*     Argument variables
*
      REAL             x
*
*     Local variables
*
      REAL             y,           yhalf,       xyhalf
*
*     The subroutine store() is intended to trick the compiler into
*     storing y, thereby avoiding incorrect results on architectures
*     with floating-point registers that have higher precision than
*     memory words (e.g., Intel x86).
*
*     The reduction by factors of two is correct for all architectures
*     with floating-point base 2; this includes all commercially
*     significant architectures in the 1990s, except for the IBM
*     mainframe S/360 -- S/390 families, which have a floating-point
*     base of 16.
*
      y = x
  100 yhalf = 0.5 * y
      xyhalf = (x + yhalf)
      CALL astore(xyhalf)
      IF (xyhalf .EQ. x) GO TO 200
      y = yhalf
      GO TO 100
  200 aeps = y
      END
