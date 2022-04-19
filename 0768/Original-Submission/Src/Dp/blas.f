
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(*),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)              
C                                                         
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.              
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.    
C     JACK DONGARRA, LINPACK, 3/11/78.                    
C                                                         
      DOUBLE PRECISION DX(*), DY(*), DA                   
      INTEGER I, INCX, INCY, IX, IY, M, MP1, N            
C                                                         
      IF (N .LE. 0) RETURN                                
      IF (DA .EQ. 0.0D0) RETURN                           
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20         
C                                                         
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS  
C          NOT EQUAL TO 1                                 
C                                                         
      IX = 1                                              
      IY = 1                                              
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1               
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1               
      DO 10 I = 1, N                                      
         DY(IY) = DY(IY) + DA*DX(IX)                      
         IX = IX + INCX                                   
         IY = IY + INCY                                   
   10 CONTINUE                                            
C                                                            
      RETURN                                              
C                                                         
C        CODE FOR BOTH INCREMENTS EQUAL TO 1              
C                                                         
                                                          
C        CLEAN-UP LOOP                                    
C                                                         
   20 CONTINUE                                            
      M = MOD(N,4)                                        
      IF (M .EQ. 0) GO TO 40                              
      DO 30 I = 1, M                                      
         DY(I) = DY(I) + DA*DX(I)                         
   30 CONTINUE                                            
      IF (N .LT. 4) RETURN                                
   40 CONTINUE                                            
      MP1 = M + 1                                         
      DO 50 I = MP1, N, 4                                 
         DY(I) = DY(I) + DA*DX(I)                         
         DY(I+1) = DY(I+1) + DA*DX(I+1)                   
         DY(I+2) = DY(I+2) + DA*DX(I+2)                   
         DY(I+3) = DY(I+3) + DA*DX(I+3)                   
   50 CONTINUE                                            
C                                                         
      RETURN                                              
C                                                         
      END                                                 
C                                                          
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)                  
C                                                          
C     COPIES A VECTOR, X, TO A VECTOR, Y.                  
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.     
C     JACK DONGARRA, LINPACK, 3/11/78.                     
C                                                          
      DOUBLE PRECISION DX(*), DY(*)                        
      INTEGER I, INCX, INCY, IX, IY, M, MP1, N             
C                                                          
      IF (N .LE. 0) RETURN                                 
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20          
C                                                          
C                                                          
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS   
C          NOT EQUAL TO 1                                  
C                                                          
      IX = 1                                               
      IY = 1                                               
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1                
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1                
      DO 10 I = 1, N                                       
         DY(IY) = DX(IX)                                   
         IX = IX + INCX                                    
         IY = IY + INCY                                    
   10 CONTINUE                                             
C                                                          
      RETURN                                               
C                                                          
C        CODE FOR BOTH INCREMENTS EQUAL TO 1               
C                                                          
C                                                          
C        CLEAN-UP LOOP                                     
C                                                          
   20 CONTINUE                                             
      M = MOD(N,7)                                         
      IF (M .EQ. 0) GO TO 40                               
      DO 30 I = 1, M                                       
         DY(I) = DX(I)                                     
   30 CONTINUE                                             
      IF (N .LT. 7) RETURN                                 
   40 CONTINUE                                             
      MP1 = M + 1                                          
      DO 50 I = MP1, N, 7                                  
         DY(I) = DX(I)                                     
         DY(I+1) = DX(I+1)                                 
         DY(I+2) = DX(I+2)                                 
         DY(I+3) = DX(I+3)                                 
         DY(I+4) = DX(I+4)                                 
         DY(I+5) = DX(I+5)                                 
         DY(I+6) = DX(I+6)                                    
   50 CONTINUE                                             
C                                                          
      RETURN                                               
C                                                          
      END                                                  

      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END


      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT,N,INCX,I,J,NN
      DOUBLE PRECISION   DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END


      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGEMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
*
*     Y := ALPHA*A*X + BETA*Y,   OR   Y := ALPHA*A'*X + BETA*Y,
*
*  WHERE ALPHA AND BETA ARE SCALARS, X AND Y ARE VECTORS AND A IS AN
*  M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  TRANS  - CHARACTER*1.
*           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
*           FOLLOWS:
*
*              TRANS = 'N' OR 'N'   Y := ALPHA*A*X + BETA*Y.
*
*              TRANS = 'T' OR 'T'   Y := ALPHA*A'*X + BETA*Y.
*
*              TRANS = 'C' OR 'C'   Y := ALPHA*A'*X + BETA*Y.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
*           M MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
*           CONTAIN THE MATRIX OF COEFFICIENTS.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCX ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCX ) ) OTHERWISE.
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE
*           VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY, BETA SPECIFIES THE SCALAR BETA. WHEN BETA IS
*           SUPPLIED AS ZERO THEN Y NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCY ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ) OTHERWISE.
*           BEFORE ENTRY WITH BETA NON-ZERO, THE INCREMENTED ARRAY Y
*           MUST CONTAIN THE VECTOR Y. ON EXIT, Y IS OVERWRITTEN BY THE
*           UPDATED VECTOR Y.
*
*  INCY   - INTEGER.
*           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           Y. INCY MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 2 BLAS ROUTINE.
*
*  -- WRITTEN ON 22-OCTOBER-1986.
*     JACK DONGARRA, ARGONNE NATIONAL LAB.
*     JEREMY DU CROZ, NAG CENTRAL OFFICE.
*     SVEN HAMMARLING, NAG CENTRAL OFFICE.
*     RICHARD HANSON, SANDIA NATIONAL LABS.
*
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     SET  LENX  AND  LENY, THE LENGTHS OF THE VECTORS X AND Y, AND SET
*     UP THE START POINTS IN  X  AND  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
*     FIRST FORM  Y := BETA*Y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        FORM  Y := ALPHA*A*X + Y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        FORM  Y := ALPHA*A'*X + Y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     END OF DGEMV .
*
      END


      SUBROUTINE DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED 3/93 TO RETURN IF INCX .LE. 0.
C
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, M, MP1, N, NINCX
C
      IF (N .LE. 0 .OR. INCX .LE. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
         DX(I) = DA*DX(I)
   10 CONTINUE

      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 CONTINUE
      M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 5
         DX(I) = DA*DX(I)
         DX(I+1) = DA*DX(I+1)
         DX(I+2) = DA*DX(I+2)
         DX(I+3) = DA*DX(I+3)
         DX(I+4) = DA*DX(I+4)
   50 CONTINUE

      RETURN

      END


      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(*), DY(*), DTEMP
      INTEGER I, INCX, INCY, IX, IY, M, MP1, N
C
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         DTEMP = DX(IX)
         DX(IX) = DY(IY)
         DY(IY) = DTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE

      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 CONTINUE
      M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         DTEMP = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 3
         DTEMP = DX(I)
         DX(I) = DY(I)
         DY(I) = DTEMP
         DTEMP = DX(I+1)
         DX(I+1) = DY(I+1)
         DY(I+1) = DTEMP
         DTEMP = DX(I+2)
         DX(I+2) = DY(I+2)
         DY(I+2) = DTEMP
   50 CONTINUE

      RETURN

      END


      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED 3/93 TO RETURN IF INCX .LE. 0.
C
      DOUBLE PRECISION DX(*), DMAX
      INTEGER I, INCX, IX, N
C
      IDAMAX = 0
      IF (N .LT. 1 .OR. INCX .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 30
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 20 I = 2, N
         IF (DABS(DX(IX)) .LE. DMAX) GO TO 10
         IDAMAX = I
         DMAX = DABS(DX(IX))
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE

      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   30 CONTINUE
      DMAX = DABS(DX(1))
      DO 40 I = 2, N
         IF (DABS(DX(I)) .LE. DMAX) GO TO 40
         IDAMAX = I
         DMAX = DABS(DX(I))
   40 CONTINUE

      RETURN

      END

      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          CA, CB
*     ..
*
*  PURPOSE
*  =======
*
*  LSAME RETURNS .TRUE. IF CA IS THE SAME LETTER AS CB REGARDLESS OF
*  CASE.
*
*  ARGUMENTS
*  =========
*
*  CA      (INPUT) CHARACTER*1
*  CB      (INPUT) CHARACTER*1
*          CA AND CB SPECIFY THE SINGLE CHARACTERS TO BE COMPARED.
*
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ICHAR
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST IF THE CHARACTERS ARE EQUAL
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     NOW TEST FOR EQUIVALENCE IF BOTH CHARACTERS ARE ALPHABETIC.
*
      ZCODE = ICHAR( 'Z' )
*
*     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
*     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
*     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
*     ICHAR('A') ON AN EBCDIC MACHINE.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
*        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     END OF LSAME
*
      END

      SUBROUTINE XERBLA ( SRNAME, INFO )
*     ..    SCALAR ARGUMENTS ..
      INTEGER            INFO
      CHARACTER*6        SRNAME
*     ..
*
*  PURPOSE
*  =======
*
*  XERBLA  IS AN ERROR HANDLER FOR THE LEVEL 2 BLAS ROUTINES.
*
*  IT IS CALLED BY THE LEVEL 2 BLAS ROUTINES IF AN INPUT PARAMETER IS
*  INVALID.
*
*  INSTALLERS SHOULD CONSIDER MODIFYING THE STOP STATEMENT IN ORDER TO
*  CALL SYSTEM-SPECIFIC EXCEPTION-HANDLING FACILITIES.
*
*  PARAMETERS
*  ==========
*
*  SRNAME - CHARACTER*6.
*           ON ENTRY, SRNAME SPECIFIES THE NAME OF THE ROUTINE WHICH
*           CALLED XERBLA.
*
*  INFO   - INTEGER.
*           ON ENTRY, INFO SPECIFIES THE POSITION OF THE INVALID
*           PARAMETER IN THE PARAMETER-LIST OF THE CALLING ROUTINE.
*
*
*  AUXILIARY ROUTINE FOR LEVEL 2 BLAS.
*
*  WRITTEN ON 20-JULY-1986.
*
*     .. EXECUTABLE STATEMENTS ..
*
      WRITE (*,99999) SRNAME, INFO
*
      STOP
*
99999 FORMAT ( ' ** ON ENTRY TO ', A6, ' PARAMETER NUMBER ', I2,
     $         ' HAD AN ILLEGAL VALUE' )
*
*     END OF XERBLA.
*
      END

