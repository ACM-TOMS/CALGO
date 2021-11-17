      SUBROUTINE POLYAN ( C, CM, N )
C  POLYAN OBTAINS INFORMATION ABOUT THE LOCATION
C  OF THE ROOTS OF A POLYNOMIAL BY USING
C  BOUND, RADIUS AND HRWTZR.
C  C IS AN N ELEMENT ARRAY CONTAINING THE COEFFICIENTS
C  NORMALIZED SO THAT THE LEADING COEFFICIENT (WHICH
C  IS NOT INCLUDED IN C) IS +1.0.
C  CM IS A WORKING ARRAY THE SAME SIZE AS C.
C  N = DEGREE OF POLYNOMIAL.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL C(N),CM(N)
C     ..
C     .. Local Scalars ..
      REAL RIN,RNL,RNU,ROUT,RPL,RPU,X
      INTEGER I,NI,NM1
C     ..
C     .. External Functions ..
      REAL BOUND,RADIUS
      LOGICAL HRWTZR
      EXTERNAL BOUND,RADIUS,HRWTZR
C  TEST FOR ZERO ROOT
      IF ( C(N) .EQ. 0.0 ) GO TO 50
C  COEFFICIENTS FOR RECIPROCAL POLYNOMIAL ARE PUT IN CM.
      CM(N) = 1. / C(N)
      NM1 = N - 1
      DO 5 I = 1, NM1
        NI = N - I
        CM(I) = CM(N) * C(NI)
5     CONTINUE
      ROUT = RADIUS ( C, N )
      RIN = 1. / RADIUS ( CM, N )
      WRITE(*, 201 ) RIN, ROUT
201   FORMAT ( " ROOTS ARE IN AN ANNULUS OF INNER RADIUS", E10.3, " AND 
     +OUTER RADIUS", E10.3 )
      RPU = BOUND ( C, N )
      IF ( RPU .NE. 0.0 ) GO TO 10
      WRITE(*, 202 )
202   FORMAT ( " THERE ARE NO REAL POSITIVE ROOTS" )
      GO TO 20
10    RPL = 1. / BOUND ( CM, N )
      WRITE(*, 203 ) RPL, RPU
C  COEFFICIENTS FOR NEGATIVE RECIPROCAL ARE PUT IN CM.
203   FORMAT  ( " THE POSITIVE ROOTS (IF ANY) ARE BETWEEN", E10.3, " AND
     +", E10.3 )
20    DO 25 I = 1, N, 2
        CM(I) = - CM(I)
25    CONTINUE
      RNU = BOUND ( CM, N )
      IF ( RNU .NE. 0.0 ) GO TO 30
      WRITE(*, 204 )
204   FORMAT  ( " THERE ARE NO NEGATIVE REAL ROOTS" )
      GO TO 40
C  COEFFICIENTS FOR NEGATIVE ROOTS ARE PUT IN CM.
30    X = -1.0
      DO 35 I = 1, N
        CM(I) = X * C(I)
        X = -X
35    CONTINUE
      RNU = -1. / RNU
      RNL = - BOUND ( CM, N )
      WRITE(*, 205 ) RNU, RNL
205   FORMAT   ( " THE REAL NEGATIVE ROOTS(IF ANY) ARE BETWEEN",  E10.3,
     + " AND", E10.3 )
40    IF ( HRWTZR ( C, N ) ) WRITE(*, 206 )
206   FORMAT   ( " THERE ARE NO ROOTS WITH POSITIVE REAL PARTS" )
      IF ( HRWTZR ( CM, N ) ) WRITE(*, 207 )
207   FORMAT   ( " THERE ARE NO ROOTS WITH NEGATIVE REAL PARTS" )
      RETURN
50    WRITE(*, 208 )
208   FORMAT  ( " POLYNOMIAL HAS A ZERO ROOT-REDUCE DEGREE" )
      RETURN
      END
      REAL FUNCTION RADIUS ( C, N )
C  RADIUS RETURNS AN UPPER LIMIT FOR THE MODULUS
C  OF THE ROOTS OF AN N DEGREE POLYNOMIAL.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL C(N)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
      RADIUS = ABS ( C(1) )
      DO 10 I = 2, N
        IF ( ABS ( C(I) ) .GT. RADIUS ) RADIUS = ABS ( C(I) )
10    CONTINUE
      RADIUS = 1. + RADIUS
      RETURN
      END
      REAL FUNCTION BOUND ( C, N )
C  BOUND RETURNS AN UPPER LIMIT FOR THE
C  POSITIVE REAL ROOTS OF AN N DEGREE POLYNOMIAL.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL C(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT
C     ..
      M = 0
      BOUND = 0.0
      DO 9 I  =  1, N
        IF ( M .GT. 0 ) GO TO 10
        IF ( C(I) .LT. 0.0 ) M = I
10      IF ( C(I) .LT. BOUND ) BOUND = C(I)
9     CONTINUE
      IF ( M .EQ. 0 ) RETURN
      BOUND = 1. + ( -BOUND )** ( 1. / FLOAT ( M ) )
      RETURN
      END
      LOGICAL FUNCTION HRWTZR ( C, N )
C  HRWTZR RETURNS .TRUE. IF ALL THE ROOTS HAVE
C  NEGATIVE REAL PARTS, OTHERWISE .FALSE. IS RETURNED.
C  IF A REAL PART IS ZERO, THEN .FALSE. IS RETURNED.
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL C(N)
C     ..
C     .. Local Scalars ..
      REAL C1
      INTEGER I,K,M
C     ..
      HRWTZR = .FALSE.
C%%
C%% Modifications as suggested by Driessen and Hunt,
C%% CACM Volume 16, Number 9, page 579, September 1973.
C%%
      IF ( C(1) .LE. 0.0 .OR. C(N) .LE. 0. ) RETURN
      C1 = C(1)
      M = N - 1
      DO 30 I = 2, M
        DO 20 K = I, M, 2
          C(K) = C(K) - C(K+1) / C1
20    CONTINUE
        C1 = C(I) / C1
        IF ( C1 .LE. 0.0 ) RETURN
30    CONTINUE
      HRWTZR = .TRUE.
      RETURN
      END

