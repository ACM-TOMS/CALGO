***********************************************************************
*                                                                     *
*     matran.f: Random matrix generator                               *
*                                                                     *
***********************************************************************
      SUBROUTINE MATRAN( N, NZR, ISEED, A, COLPTR, ROWIND )
*     ..
*     .. Scalar Arguments ..
      INTEGER            N, NZR 
*     ..
*     .. Array Arguments ..
      INTEGER            COLPTR( * ), ROWIND( * ), ISEED( 4 )
      DOUBLE PRECISION   A( * )                                  
*
*
*  Purpose
*  =======
*
*  Generates a sparse nonsymmetric matrix with nonzeros placed 
*  randomly according to user specified the number of nonzeros 
*  per column. The random numbers are from the uniform distribution
*  on (-1,1).  The generated matrix is stored in compressed
*  column format.                 
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          the order of matrix A 
*
*  NZR     (input) INTEGER
*          The number of nonzero entries per column.
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  A       (output) DOUBLE PRECISION array, dimension (N*NZR) 
*          the numerical values of the generated matrix. 
*
*  COLPTR  (output) INTEGER array, dimension (N+1)
*          the column start pointers
*  
*  ROWIND  (output) INTEGER array, dimension (N*NZR) 
*          The row indices 
*
*  ====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 ) 
*     .. 
*     .. Local Scalars ..
      INTEGER            IR, IRT, J, K, K1, KM   
      DOUBLE PRECISION   GN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARAN
*     ..
*     .. Intrinsic functions ..
      INTRINSIC          ABS
*
*     .. Executable statements ..
*
      IR = 0                                                            
*
      DO 10 J = 1,N                                                     
*
*        set up the column start pointer
*
         COLPTR(J) = (J-1)*NZR +1 
*
         DO 20 K = 1, NZR 
*
*           Select index of nonzero in column J   
*
   30       GN = TWO*DLARAN( ISEED ) - ONE  
*
            IRT = N*ABS(GN) + 1   
            IF( IRT.GT.N ) 
     $          GO TO 30    
*
            K1 = K-1  
            IF( K1.EQ.0 ) 
     $         GO TO 60 
*
            DO 40 KK = 1, K1  
               KM = IR - KK + 1 
               IF( IRT.EQ.ROWIND( KM ) ) 
     $            GO TO 30  
   40       CONTINUE    
*
   60       IR = IR + 1 
            ROWIND( IR ) = IRT  
            A( IR ) = TWO*DLARAN( ISEED ) - ONE
*
   20    CONTINUE
*
   10 CONTINUE                                                          
*
      COLPTR( N+1 ) = IR + 1 
*
      RETURN                                                            
*
*     End of MATRAN
      END                                                               

*
*
      DOUBLE PRECISION FUNCTION DLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      DLARAN = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $         ( DBLE( IT4 ) ) ) ) )
      RETURN
*
*     End of DLARAN
*
      END
