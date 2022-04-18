
      INTEGER FUNCTION EXOR(IIN,JIN)
C
C     THIS FUNCTION CALCULATES THE EXCLUSIVE-OR OF ITS
C     TWO INPUT PARAMETERS
C
C     .. Scalar Arguments ..
      INTEGER IIN,JIN
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      I = IIN
      J = JIN
      K = 0
      L = 1
C
   10 IF (I.EQ.J) THEN
          EXOR = K
          RETURN

      END IF
C
C     CHECK THE CURRENT RIGHT-HAND BITS OF I AND J.
C     IF THEY DIFFER, SET THE APPROPRIATE BIT OF K.
C
      IF (MOD(I,2).NE.MOD(J,2)) K = K + L
      I = I/2
      J = J/2
      L = 2*L
      GO TO 10

      END
      DOUBLE PRECISION FUNCTION UNI()
*
*     Random number generator, adapted from F. James
*     "A Review of Random Number Generators"
*      Comp. Phys. Comm. 60(1990), pp. 329-344.
*
C     .. Parameters ..
      DOUBLE PRECISION TWOM24
      PARAMETER (TWOM24=1D0/16777216.0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CARRY
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION SEEDS(24)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Save statement ..
      SAVE I,J,CARRY,SEEDS
C     ..
C     .. Data statements ..
      DATA I,J,CARRY/24,10,0.0d0/
      DATA SEEDS/0.8804418d0,0.2694365d0,0.0367681d0,0.4068699d0,
     +0.4554052d0,
     +     0.2880635d0,0.1463408d0,0.2390333d0,0.6407298d0,
     +0.1755283d0,0.7132940d0,
     +     0.4913043d0,0.2979918d0,0.1396858d0,0.3589528d0,
     +0.5254809d0,0.9857749d0,
     +     0.4612127d0,0.2196441d0,0.7848351d0,0.4096100d0,
     +0.9807353d0,0.2689915d0,
     +     0.5140357d0/
C     ..


      UNI = SEEDS(I) - SEEDS(J) - CARRY
      IF (UNI.LT.0) THEN
          UNI = UNI + 1
          CARRY = TWOM24

      ELSE
          CARRY = 0
      END IF

      SEEDS(I) = UNI
      I = 24 - MOD(25-I,24)
      J = 24 - MOD(25-J,24)
      RETURN

      END
