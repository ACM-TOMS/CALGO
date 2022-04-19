!***********************************************************************
!***********************************************************************
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
!
!  I/O UNIT NUMBERS.
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS.
!
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!
!  INTEGERS.
!
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!
!    I1MACH( 7) = A, THE BASE.
!
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!
!    I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
!  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
!  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
!
!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
!  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
!  FOR IMACH(1) - IMACH(4).
!
      INTEGER IMACH(16),OUTPUT
!
      EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   32 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
       DATA IMACH(11) /   24 /
       DATA IMACH(12) / -125 /
       DATA IMACH(13) /  128 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!
!      DATA IMACH( 1) /    7 /
!      DATA IMACH( 2) /    2 /
!      DATA IMACH( 3) /    2 /
!      DATA IMACH( 4) /    2 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   33 /
!      DATA IMACH( 9) / Z1FFFFFFFF /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -256 /
!      DATA IMACH(13) /  255 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) / -256 /
!      DATA IMACH(16) /  255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -50 /
!      DATA IMACH(16) /  76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -32754 /
!      DATA IMACH(16) /  32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   60 /
!      DATA IMACH( 6) /   10 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   48 /
!      DATA IMACH( 9) / 00007777777777777777B /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   48 /
!      DATA IMACH(12) / -974 /
!      DATA IMACH(13) / 1070 /
!      DATA IMACH(14) /   96 /
!      DATA IMACH(15) / -927 /
!      DATA IMACH(16) / 1070 /
!
!     MACHINE CONSTANTS FOR CONVEX C-1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   53 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!
!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /   102 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    63 /
!      DATA IMACH( 9) /  777777777777777777777B /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    47 /
!      DATA IMACH(12) / -8189 /
!      DATA IMACH(13) /  8190 /
!      DATA IMACH(14) /    94 /
!      DATA IMACH(15) / -8099 /
!      DATA IMACH(16) /  8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
!
!      DATA IMACH( 1) /   11 /
!      DATA IMACH( 2) /   12 /
!      DATA IMACH( 3) /    8 /
!      DATA IMACH( 4) /   10 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) /32767 /
!      DATA IMACH(10) /   16 /
!      DATA IMACH(11) /    6 /
!      DATA IMACH(12) /  -64 /
!      DATA IMACH(13) /   63 /
!      DATA IMACH(14) /   14 /
!      DATA IMACH(15) /  -64 /
!      DATA IMACH(16) /   63 /
!
!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
!
!      DATA IMACH( 1) /       5 /
!      DATA IMACH( 2) /       6 /
!      DATA IMACH( 3) /       0 /
!      DATA IMACH( 4) /       6 /
!      DATA IMACH( 5) /      24 /
!      DATA IMACH( 6) /       3 /
!      DATA IMACH( 7) /       2 /
!      DATA IMACH( 8) /      23 /
!      DATA IMACH( 9) / 8388607 /
!      DATA IMACH(10) /       2 /
!      DATA IMACH(11) /      23 /
!      DATA IMACH(12) /    -127 /
!      DATA IMACH(13) /     127 /
!      DATA IMACH(14) /      38 /
!      DATA IMACH(15) /    -127 /
!      DATA IMACH(16) /     127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /   43 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   63 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
       DATA IMACH(16) /  1024 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z7FFFFFFF /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /
!
!     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
!
!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
!
!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   6 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z'7FFFFFFF' /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  62 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  62 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   54 /
!      DATA IMACH(15) / -101 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   62 /
!      DATA IMACH(15) / -128 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     16-BIT INTEGER ARITHMETIC.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) / 32767 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
!
!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!     SUPPLIED BY IGOR BRAY.
!
!      DATA IMACH( 1) /            1 /
!      DATA IMACH( 2) /            1 /
!      DATA IMACH( 3) /            2 /
!      DATA IMACH( 4) /            1 /
!      DATA IMACH( 5) /           32 /
!      DATA IMACH( 6) /            4 /
!      DATA IMACH( 7) /            2 /
!      DATA IMACH( 8) /           31 /
!      DATA IMACH( 9) / :17777777777 /
!      DATA IMACH(10) /            2 /
!      DATA IMACH(11) /           23 /
!      DATA IMACH(12) /         -127 /
!      DATA IMACH(13) /         +127 /
!      DATA IMACH(14) /           47 /
!      DATA IMACH(15) /       -32895 /
!      DATA IMACH(16) /       +32637 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!
!      DATA IMACH( 1) /     0 /
!      DATA IMACH( 2) /     0 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     0 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    6 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /
!
!     MACHINE CONSTANTS FOR VAX.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 999
      I1MACH=IMACH(I)
!/6S
!/7S
      IF(I.EQ.6) I1MACH=1
!/
      RETURN
  999 WRITE(OUTPUT,1999) I
 1999 FORMAT(' I1MACH - I OUT OF BOUNDS',I10)
      STOP
      END

