C.....................................................................
C     Double Precision
C
C     This function generates random numbers, RND, such that
C               -1 < RND < 1
C      which are added to given function measurements, in order
C      to construct data sets for testing subroutine L2WPMA.
C.....................................................................
      DOUBLE PRECISION FUNCTION RND(ISEED)
C
C     Function RND can generate 65536 "random" numbers before
C      repeating itself. It will work on most computers for
C      which 2**31 - 1 .LT. MAXINT.
C
C.... I N P U T ....
C     ISEED     Integer variable between 1 and 65535 .
C
C.... O U T P U T ....
C     RND       The return of the function evaluation,
C                between -1 and 1.
C
C.... M E T H O D (Linear congruential) ....
C      Each number in the random sequence, R(K) say,
C      is calculated from its predecessor, R(K - 1), using the formula:
C               R(K)=(MULTIPLIER * R(K - 1) + INCREMENT) MOD 2**16
C      where MULTIPLIER = 25173 and INCREMENT = 13849.
C
C..... Ref: P.GROGONO, Programming in PASCAL, Addison-Wesley Pub.Co.,
C          1980, P.118.
C.....................................................................
C     .. Scalar Arguments ..
      INTEGER ISEED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,MOD
C     ..
      RND = DBLE(ISEED)/65535.
      RND = 2*RND - 1
      ISEED = 25173*ISEED + 13849
      ISEED = MOD(ISEED,65536)
C
      RETURN
      END
