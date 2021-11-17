*--------1--------2--------3---------4---------5---------6---------7-*
*                  BHESS and Auxiliary Routines                      *
*                                                                    *
*     Fed a square matrix N by N matrix dimensioned A(LDA,N+1)       *
*     and with the N by N matrix having upper left corner A(1,1)     *
*     BHESS overwrites the upper N by N left corner with a           *
*     similar banded Hessenberg matrix and the multipliers           *
*     required to accomplish the elimination.                        *
*                                                                    *
*     The code is essentially Fortran 77 but contains minor          *
*     extensions that are accepted by all Fortran 90 compilers       *
*     (and indeed by all Fortran 77 compilers we have found).        * 
*--------------------------------------------------------------------*
*                                                                    *
*               CONTENTS                                             *
*                                                                    *
*     BHESS  -- Performs H = inv(Z)* A * Z                           *
*               A a general square matrix.                           *
*               H an upper Hessenberg matrix.                        *
*               Z is a product of elementary matrices.               *
*               H and the multipliers for Z overwrite A              * 
*                                                                    *
*     BHPIVT  - auxiliary routine used by BHESS to determine which   *
*               rows can be eliminated, and to choose permutations   *  
*                                                                    *
*     BHAP1   - BHESS gives H = inv(Z)*A*Z.                          *
*               BHAP1 performs x^T <-- x^T * inv(Z)                  *
*               where X is stored as the two vectors XRE and XIM     * 
*                                                                    *
*     BHAP2   - x^T <-- x^T * Z  where x is a vector                 *
*                                                                    *
*     BHAP3   - (XRE,XIM) <-- Z * (XRE,XIM)                          *
*                                                                    *
*     BHAPC   - Condition number estimator for Z                     *
*                                                                    *
*     BHRPCK  - Zeros out the multipliers so that only the upper     *
*               Hessenberg matrix remains                            *
*                                                                    *
*     BHUNPK  - Puts the upper Hessenberg matrix is diagonal         *
*               storage                                              *
*                                                                    *
*               For some other useful routines, see the anonymous    *
*               ftp site                                             *
*                                                                    *
*     ftp cs.fit.edu                                                 *
*     cd pub/howell                                                  * 
*                                                                    *
*     In particular, the br eigenvalue algorithm                     *
*     (SIMAX, summer 1999) is often successful in extracting         *
*     eigenvalues of the matrix returned by BHESS in time            *
*     proportional to the square of the matrix size.                 *
*                                                                    *
*--------------------------------------------------------------------*
*
*     Code authored by Gary Howell and Nadia Diaa
*        Authors offer it for public use on condition that
*        our authorship recognized.  We hope the code will be useful
*        but do not warrant its reliability. 
*
*        Note: much of the code is a revision of ATOTRI
*            by Al Geist (SIMAX 1991).  The BHESS version differs 
*            in that multipliers are bounded and the returned
*            Hessenberg matrix is not necessarily tridiagonal. 
*
*     Thanks to Charles Fulton, David Watkins, 
*        Michael Saunders, and Al Geist for suggestions.
*
*     First submitted to TOMS in 1994.  
*
*--------------------------------------------------------------------*


      SUBROUTINE BHESS(LDA,A,N,PIVOTS,NU,KRVECT,KROW,V,TOL,INFO)

*     ----------------------
*     
*     .. Scalar Arguments ..
*     ----------------------

      INTEGER LDA, N, INFO
      DOUBLE PRECISION TOL

*     --------------------------------------------
*      
*     .. Array Arguments ..
*     --------------------------------------------

      INTEGER PIVOTS(*), NU(*), KRVECT(*), KROW(*)
      DOUBLE PRECISION A(LDA,*), V(*)
 
*     --------------------------------------------
* 
*     .. External Functions ..
*     --------------------------------------------

      DOUBLE PRECISION DDOT

*     ---------------------------------------------------------
*
*    DDOT is from the DBLAS-1 library.
*    ----------------------------------------------------------
*     
*  Purpose:
*    This subroutine reduces an n-by-n real general matrix A to
*    small-band Hessenberg form using elementary Gaussian 
*    similarity transformations.
*
*    At each step k the permutation that minimizes the maximum entry
*    in the transformation matrix that reduces column k and then row k
*    is applied.
*
*               July 9, 2004
*       Added quick returns for N less than 2 or N greater than LDA 
*
*               September 15, 1997
*       Cache performance is improved by performing both the row
*       and column similarity transformations on a given column
*       in one access.  Execution time relative to LAPACK DGEHRD is
*       (on a RS-6000) 85% for size 200 and 80% for size 1200.
*
*               June 25, 1995
*       This version of BHESS does not use BLAS level 2.  Outer loops 
*       corresponding to the BLAS 2 operation DGER and DGMEMV are in
*       the column variable.  This version is twice as fast on various
*       IBM workstations.
*
*
*  Arguments:
*   LDA      -integer
*             leading dimension of A
*
*   A        -double precision array of dimension (LDA,N+1)
*             On entry, the first N columns of A contain the matrix 
*             being reduced.
*             On exit the first N columns of A are overwritten by 
*             its banded Hessenberg form and the multipliers used 
*             to produce the small band form.  
*             The extra (N+1st) column of A is used as a work 
*             vector.  
*
*   N        -integer
*             N specifies the order of the matrix A.
*             N must be nonnegative.
*             not modified.
*
*   PIVOTS   -integer vector of dimension (LDA)
*             On exit pivots contains the pivot sequence used during
*             the reduction (permutation vector).
*
*   NU       -integer vector of dimension N.
*             NU is initialized as all ones
*             on exit NU(KR) is the number of consecutive
*             presumed nonzero entries to right of 
*             diagonal in row KR. 
*
*   KRVECT -dummy integer vector of dimension N.
*             on entry: n/a as it is a dummy var passed only so we
*                       can use variable dimension N.
*             on exit: n/a
*             KRVECT is initialized as all ones.
*             Whenever row KR is zeroed, KRVECT(KR) becomes 0.
*             KRVECT is used internally to keep track of which rows
*             are eligible to be zeroed.
*
*   KROW     -integer vector of dimension N.  
*             On exit KROW(K) contains the row eliminated after 
*             column K was eliminated.  If the row elimination was 
*             never performed, then KROW(K) is zero.  
*
*   TOL        not modified
*            - when the L2 norms of rows and columns of multipliers are equal,
*              then TOL corresponds to the maximal standard deviation
*              of multipliers from zero. More generally TOL is the maximal
*              allowed root mean square of standard deviations
*              from zero of row and column multipliers. 
*
*   V        -double precision dummy vector
*
*
*   INFO     -integer
*             On exit, INFO is set to
*             0  normal return.
*             1  if an error is detected.
*             -1 if N less than 3
*             -2 if N greater than LDA 
*
*     .. Local Scalars ..
*     --------------------------------------------------------------------

      INTEGER FIRST1, I, J, K, KR, L, PIVP, ERR, CNT
      DOUBLE PRECISION TEMP, MINMLT, ONE, INPROD
      PARAMETER (ONE = 1.0)

*     --------------------------------------------------------------------
*    
*     FIRST1 is the first 1 in KRVECT, i.e., the index of the first
*     nonzero row.
*
      
*     .. External Subroutines ..
*     --------------------------------------------------------------------

      EXTERNAL BHPIVT, DAXPY, DSCAL

*     DBLAS DDOT,DAXPY,DSCAL
      
*     
*=============Executable Statements ===============================
*     --------------
*
*     ------------------------------------
*     Quick returns for wrong sized matrices
*     ------------------------------------
      
      IF (N. LE. 2) THEN 
         INFO = -1
         PRINT*,' Matrix size N must be greater than 2' 
         RETURN
      ELSEIF (N. GE. LDA) THEN
         INFO = -2
         PRINT*,' Matrix size N must be smaller than', LDA
         RETURN
      ENDIF 

*     --------------
*     Initialization 
*     --------------

      INFO = 0
      CNT = 0
      L = 1
      PIVOTS(1) = 1
      PIVOTS(N) = N
      K = 1
      KR = 1
      FIRST1 = 1
      DO  I = 1,N-1
         NU(I) = 1
         KRVECT(I) = 1
         KROW(I) = 0
      END DO
      NU(N) = 0
      KRVECT(N) = 1
      KROW(N) = 0

*     ------------------------------------
*     For each column of the matrix
*     ------------------------------------

      DO  1000 K=1, N-2

*        ------------------------------------
*        We will come to 10 when we use same column with a different row
*         i.e., do not increment column counter K.
*        ------------------------------------

 10      CONTINUE

*        ---------------------------------
*          Find a suitable pivot.  Integer variable ERR returned by 
*          BHPIVT contains information as to whether a row elimination
*          should be performed.         
*        ------------------------------------------------------

         CALL BHPIVT ( LDA,A,N,K,KR,TOL,PIVP,INPROD,MINMLT,ERR )
         PIVOTS(K+1) = PIVP

*        ------------------------------------------------------
*          Check for deflation , i.e. both row & column are zero.
*          Column elimination is not necessary.         
*          (instead of going directly to 1000, update a few 
*          variables first and then fall into 1000 continue statement)
*        ------------------------------------------------------

         IF( ERR .EQ. 1 ) THEN
            L = K+1
            GO TO 888
         ENDIF

*        -----------------------------------------------------------
*          Check for (inner product = zero) OR (inner product < tolerance)
*        -----------------------------------------------------------

         IF ((ERR .EQ. -1) .OR. (ERR .EQ. -2)) THEN

*           --------------------------------------------------------
*            Row KR is not eligible so try the next available row
*               (if any)
*            First, update NU(KR) to indicate failure to zero out row KR
*           --------------------------------------------------------

            IF (KR .EQ. K) THEN

*              -------------------------------------------------------
*               We have tried all eligible rows with col K, 
*               all resulted in ERR =-1 or -2. 
*               This means we cannot zero out any rows, so just 
*               go ahead and zero out col K anyway.
*               update subsequent entries of KR to reflect wider
*               bandwidth resulting from failure to zero out row KR.
*               Zero out column K if needed     
*              --------------------------------------------------------

               IF( ERR .NE. 3 ) THEN

*                ------------------------------------------------------
*                 *********** Zero out column K **********************
*                 but no row is zeroed in conjunction with column K.
*
*                  Interchange row and column
*                ------------------------------------------------------

                  IF( PIVP .NE. K+1 ) THEN
                     CALL DSWAP (N,A(PIVP,1),LDA,A(K+1,1),LDA)
                     CALL DSWAP (N,A(1,PIVP),1,A(1,K+1),1) 
                  ENDIF
                  TEMP = A(K+1,K)
                  CALL DSCAL((N-(K+2)+1),ONE/TEMP,A(K+2,K),1)

*                 ----------------------------------------------------
*                   The next block of code combines the rank 1 update
*                     and the matrix vector multiply needed to 
*                     accomplish an elimination of column K
*                     and preserve similarity.
*                     Each column K+1 to N is accessed only once for the 
*                     conbined rank 1 update and matrix vector multiply.                
*                 --------------------
*                    Rank 1 update
*                 --------------------
*                    Update is addition of -u*vtr where 
*                    u is K+2 to N entries of column K
*                    vtr is K+1 to N entries of row K+1
*                    The upper left corner of updated block is A(K+2,K).
*                 ---------------------
*                    Matrix vector product
*                 ---------------------
*                    For KR(I) not zero  
*                     column K+1 of A starting at row FIRST1 
*                     is incremented by A multiplied by the column 
*                     K of A (entries K+2 to N)
*                     where the entries of A used consist of the block 
*                     starting at column K+2 and the row K.  
*
*                 -----------------------------------------------------

                  CALL DAXPY(N-K-1,-A(K+1,K+1),A(K+2,K),1,A(K+2,K+1),1)
                  DO  J = K+2,N
                     CALL DAXPY(N-K-1,-A(K+1,J),A(K+2,K),1,A(K+2,J),1)
                     CALL DAXPY(N-K-1,A(J,K),A(K+2,J),1,A(K+2,K+1),1)
                     DO  I = FIRST1,K+1
                        IF ( KRVECT(I) .NE. 0 ) THEN
                           A(I,K+1) = A(I,K+1) + A(J,K)*A(I,J)
                        ENDIF
                     END DO
                  END DO
               ENDIF

*              --------------------------------------------------------
*               End of (ERR .NE. 3)
*               Eliminated column K so drop to bottom of main loop
*               and start work on column K+1 
*              --------------------------------------------------------

               GO TO 999
            ELSE

*           ----------------------------------------------
*             KR is smaller than K,
*             could not zero out row KR      
*             find next nonzero row to try with col K
*           ----------------------------------------------

               DO  I = KR+1,K
                  IF (KRVECT(I) .EQ. 1) GO TO 451 !  Row I has 
                                                  !  not been eliminated.
               END DO
 451           CONTINUE
               KR = I

*              ----------------------------------------------------
*               Only change KR, not FIRST1,
*               now go find pivot using same col K with new row KR,
*               so bypass loop increment of K by going to 10
*              ----------------------------------------------------

               GO TO 10   ! where we will see if Row I can be
                          ! given a paired elimination with Column K 
            ENDIF
         ENDIF

*        -------------------------------------------------------------
*          Check if row-column inner product large enough.
*          If so (ERR .NE. 2), then zero out row KR along with column K
*        ------------------------------------------------------------

         IF( ERR .NE. 2 ) THEN

*           ------------------------------------------------------
*            zero out row KR and set KROW(K) as KR
*           *****************************************************
*            Interchange row and column (using pivot chosen in BHPIVT)
*           -----------------------------------------------------

            IF( PIVP .NE. K+1 ) THEN
               CALL DSWAP (N,A(PIVP,1),LDA,A(K+1,1),LDA)
               CALL DSWAP (N,A(1,PIVP),1,A(1,K+1),1) 
            ENDIF
            KROW(K)=KR
            TEMP = A(KR,K+1)
            KRVECT(KR) = 0

*           -----------------------------------------------
*          The vector V(I) is the row multipliers.  These overwrite
*            the KR row of the matrix and are stored as V(I) to
*            avoid reaccessing columns of the matrix A.
*           -----------------------------------------------

            DO  I = K+2,N
               V(I) = A(KR,I)/TEMP
               A(KR,I) = V(I)
            END DO
            DO I = FIRST1,K
               IF( KRVECT(I) .NE. 0 ) A(I,N+1)=A(I,K+1)
            END DO

*           ------------------------------------------------
*           A backup copy of column K+1 of the matrix is used.
*             The copy stored in column N+1 .
*             The copy is updated by multiples of columns as part
*             of the matrix vector multiply associated with the row
*             elimination.  The other is used for the the rank one
*             update used to eliminate row KR.
*          -----------------------------------------------------

            CALL DCOPY(N-K,A(K+1,K+1),1,A(K+1,N+1),1) 

*           ------------------------------------------------
*           Row K+1 of columns K and K+1
*             is updated by the matrix vector
*             multiply associated with 
*             eliminating column K.
*           ------------------------------------------------    

            A(K+1,N+1) = A(K+1,N+1) + DDOT(N-K-1,A(K+2,N+1),1,V(K+2),1)
            A(K+1,K) = A(K+1,K) + DDOT(N-K-1,A(K+2,K),1,V(K+2),1)

*           ------------------------------------------------
*           The multipliers to eliminate column K are stored
*             in the entry they are used to eliminate
*           ------------------------------------------------

            TEMP = A(K+1,K)
            IF ( ABS(TEMP) .GT. 1.D-12 ) THEN  
               CALL DSCAL((N-(K+2)+1),ONE/TEMP,A(K+2,K),1)
            ENDIF 

*           -------------------------------------------------
*             The following loop on columns accesses each column
*               only once in the course of combined a combined
*               similarity transformation to eliminate both
*               column K and row KR
*           -------------------------------------------------

            DO J = N,K+1,-1

*              -----------------------------------------------
*               The rank one update for eliminating row KR.
*             -----------------------------------------------

               IF ( J .NE. K+1 ) THEN
                  DO I = FIRST1,K
                     IF ( KRVECT(I) .NE. 0 ) THEN
                        A(I,J) = A(I,J) - A(I,K+1)*V(J)
                     ENDIF
                  END DO
                  CALL DAXPY(N-K,-V(J),A(K+1,K+1),1,A(K+1,J),1)
               ENDIF

*              ----------------------------------------------
*              Increment row K+1 by the inner product of multipliers
*                and succeeding rows.
*              ----------------------------------------------

               IF ( J .NE. K+1 ) THEN
                  TEMP = DDOT(N-K-1,A(K+2,J),1,V(K+2),1)
                  A(K+1,J) = A(K+1,J) + TEMP
               ENDIF

*              ----------------------------------------------
*              Increment column K+1 by a multiple of column J
*                and perform rank one update corresponding to
*                elimination of column K
*              ----------------------------------------------

               IF ( J .NE. K+1 ) THEN
                  CALL DAXPY(N-K-1,A(J,K),A(K+2,J),1,A(K+2,N+1),1)
                  CALL DAXPY(N-K-1,-A(K+1,J),A(K+2,K),1,A(K+2,J),1)
                  DO I = FIRST1,K+1
                     IF ( KRVECT(I) .NE. 0 ) THEN
                        A(I,N+1) = A(I,N+1) + A(J,K)*A(I,J)
                     ENDIF
                  END DO
               ENDIF
            END DO

*           ----------------------------------------------------
*            Overwrite the old column K+1 with the new version which
*              has been incrementing in column N+1 and do the
*              last bit of the rank one update for column elimination.
*           -----------------------------------------------------

            DO I = FIRST1,K+1
               IF ( KRVECT(I) .NE. 0 ) THEN
                  A(I,K+1) = A(I,N+1)
               ENDIF
            END DO
            CALL DCOPY(N-K-1,A(K+2,N+1),1,A(K+2,K+1),1) 
            CALL DAXPY(N-K-1,-A(K+1,K+1),A(K+2,K),1,A(K+2,K+1),1)
         
         ENDIF
 888     CONTINUE

*        ---------------------------------------------------
*          We have zeroed row KR so reflect this information
*          by updating KRVECT(KR) and then FIRST1
*        ---------------------------------------------------

         NU(KR)=K+1-KR
         KRVECT(KR) = 0
         DO I = FIRST1, K + 1
            IF ( KRVECT(I) .EQ. 1 ) GO TO 895
         END DO
 895     CONTINUE
         FIRST1 = I      
 999     CONTINUE

*        -----------------------------------------------------
*          Going to get a new col K, so we will want to try 
*           it out with first nonzero row 
*        -----------------------------------------------------

         KR = FIRST1
 1000 CONTINUE
      DO  I = 1,N-1
         IF( KRVECT(I) .EQ. 1 ) NU(I) = N-I
      END DO
      RETURN
*
      END      ! End of BHESS

*     -------------------------------------
*
*
*     ============================================================
*
*     ------------------------------------------------------------

      SUBROUTINE BHPIVT (LDA,A,N,K,KR,TOL,PIVP,INPROD,MINMLT,ERR )

*     ------------------------------------------------------------
*
*     .. Scalar Arguments ..
*     ------------------------------------------------------------

      INTEGER LDA, N, K, KR, PIVP, ERR
      DOUBLE PRECISION MINMLT, TOL, INPROD

*     ------------------------------------
*
*     .. Array Arguments ..
*     ------------------------------------

      DOUBLE PRECISION  A(LDA,*)

*     ------------------------------------
*
*  Purpose:
*
*    This subroutine tests whether row KR is eligible to eliminated
*    at in conjunction with column K.  If so, it finds the pivot 
*    that minimizes the maximum entry
*    in the matrix N, where N=N_rN_c and N_c reduces column k and
*    inv(N_r) reduces row k.  In this event return err=0.
*
*    If the row is not eligible for elimination with the kth column, then the
*    pivot is chosen to correspond to the largest entry below the subdiagonal
*    in column k.  Then return err=-2.  
*
*    Routine also checks if the problem can be deflated at step k
*    returning err=1 if so. If inner product of the row and column is
*    zero then a pivot is selected based on maximum row/col entries
*    and err is set to 2. 
*
*  Arguments:
*    A   -double precision array of dimension (n,n)
*         On entry A contains the partially reduced matrix.
*         Not modified
*
*    N   -integer
*         N specifies the order of the matrix A.
*         N must be at least zero.
*         Not modified.
*
*    K   -integer
*         K specifies the column under consideration.
*         Not modified
*
*    KR  -integer
*         KR specifies the row under consideration.
*         Not modified
*
*   TOL      -double precision scalar is between 0 and +1
*            on entry contains allowable cosine of angle 
*               between col k and row kr
*            not modified
*
*  PIVP  -integer
*         On return PIVP contains the pivot index
*         necessary to minimize the maximum entry in N.
*
* MINMLT -double precision
*         On return MINMLT contains the minimized maximum entry.
*
*   ERR  -integer
*         On exit ERR is set to
*==>      -2  inner product is too small.
*         -1  inner product = 0, recovery step required.(obsolete)
*          0  normal return.
*          1  problem deflated, row and column are zero
*          2  problem deflated, row is zero
*          3  problem deflated, column is zero
*
*     .. Local Scalars ..
*     -----------------------------------------------

      INTEGER I, PIVC, PIVR
      DOUBLE PRECISION MAXCOL, MAXROW, NMXCOL, NMXROW
      DOUBLE PRECISION MAXNC, MAXNR, MXDIAG, V1
      DOUBLE PRECISION TEMC, TEMR, TEMP,SUMC,SUMR
      DOUBLE PRECISION BIG,TOLK

*     -----------------------------------------------
*
*     .. Parameters ..
*     -----------------------------------------------

      PARAMETER (BIG = 1.0D8)

*     -----------------------
*
*
*     .. External Functions ..
*     none
*
*     .. External Subroutines ..
*     none
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS
c     INTRINSIC MAX
*
*=================== Executable Statements =================================
*
*     Initialize returned values
*     -------------------------

      PIVP = K+1
      MINMLT = BIG
      ERR = 0
      TOLK=TOL*(N-K)
      TOLK=DMIN1(1.D0/TOLK,1.D0)

*     -----------------------------------------------------
*     Find maximum and next-to-max entries in row and col k
*     and compute inner product.
*     -----------------------------------------------------

      PIVC = K+1
      PIVR = K+1
      NMXCOL = 0.D0
      NMXROW = 0.D0
      MAXCOL = A(K+1,K)
      MAXROW = A(KR,K+1)
      SUMC   = MAXCOL*MAXCOL
      SUMR   = MAXROW*MAXROW
      INPROD = MAXCOL * MAXROW
      MAXCOL = DABS(MAXCOL)
      MAXROW = DABS(MAXROW)
      DO  I=K+2, N
         TEMC = A(I,K)
         TEMR = A(KR,I)
         SUMC   = SUMC + TEMC*TEMC
         SUMR   = SUMR + TEMR*TEMR
         INPROD = INPROD + TEMC * TEMR
         TEMC = DABS(TEMC)
         TEMR = DABS(TEMR)
         IF( TEMC .GT. NMXCOL ) THEN
            IF( TEMC .GT. MAXCOL ) THEN
               NMXCOL = MAXCOL
               MAXCOL = TEMC
               PIVC = I
            ELSE
               NMXCOL = TEMC
            ENDIF
         ENDIF
         IF( TEMR .GT. NMXROW ) THEN
            IF( TEMR .GT. MAXROW ) THEN
               NMXROW = MAXROW
               MAXROW = TEMR
               PIVR = I
            ELSE
               NMXROW = TEMR
            ENDIF
         ENDIF
      END DO 

*     --------------------------------------------------
*     Check for deflation
*     --------------------------------------------------

      IF( MAXCOL .EQ. 0.D0 .AND. MAXROW .EQ. 0.D0 ) THEN
         ERR = 1
         RETURN
      ENDIF
      IF( MAXROW .EQ. 0.D0 ) THEN
         ERR = 2
         PIVP = PIVC
         RETURN
      ENDIF
      IF( MAXCOL .EQ. 0.D0 ) THEN
         ERR = 3
         PIVP = PIVR
         RETURN
      ENDIF

*     ---------------------------------------------------
*     Check for inner product too small.
*     ---------------------------------------------------

      SUMC = DSQRT(SUMC)
      SUMR = DSQRT(SUMR)
      IF( DABS(INPROD/SUMC/SUMR) .LT. TOLK ) THEN

*        ---------------------------------------------------
*         Inner product of row and column is too small for elimination
*         of row so return pivot to correspond to largest column element
*        ---------------------------------------------------

         PIVP = PIVC
         ERR = -2
         RETURN
      ENDIF

*     ---------------------------------------------------
*     Calculate maximum multipliers over all permutations i and
*     determine the minimum triple (maxnc,maxnr,mxdiag)_i.
*     In this version, the multipliers are selected for the 
*     row elimination first, then column elimination.   
*     ---------------------------------------------------

      DO  I=K+1, N
         V1 = A(KR,I)
         IF( V1 .EQ. 0.D0 ) THEN
            MAXNR = BIG
c        ELSE IF( I .EQ. PIVC ) THEN
         ELSE IF( I .EQ. PIVr ) THEN
            MAXNR = DABS(NMXROW/V1)
         ELSE
            MAXNR = DABS(MAXROW/V1)
         ENDIF
c        IF( I .EQ. PIVR ) THEN
         IF( I .EQ. PIVc ) THEN
            MAXNC = DABS(NMXCOl*V1/INPROD)
         ELSE
            MAXNC = DABS(MAXCOl*V1/INPROD)
         ENDIF
         MXDIAG = ABS(V1*A(I,K)/INPROD)
         TEMP = MAXNR
         IF( TEMP .LT. MAXNC ) TEMP = MAXNC
         IF( TEMP .LT. MXDIAG ) TEMP = MXDIAG
         IF( TEMP .LT. MINMLT ) THEN
            MINMLT = TEMP
            PIVP = I
         ENDIF
      END DO 
* 
      RETURN
      END     ! End of BHPIVT

*     -------------------------------------
*
*====================================================

*=========================================================================

      SUBROUTINE BHAP1(N,Z,LDZ,XRE,XIM,IPVT,KROW,IFLAG)
*
      INTEGER N, LDZ, IPVT(*), KROW(*), IFLAG
      DOUBLE PRECISION Z(LDZ,*), XRE(*), XIM(*)

*---------------------------------------------------------
*
*
*  Purpose:
*     This routine applies the accumulated elementary
*     transformations inv(Z) to a row vector.  The matrix
*     Z is the product of the elementary transformations
*     and permutation matrices that near-tridiagonalize the
*     original matrix A.  The computation is of the form
*         inv(Z)*A*Z = H, where
*     inv(Z) = Lc(n-2) * inv(Lr(n-2)) *  P(n-2) *  
*              Lc(n-3)*inv(Lr(n-3) * P(n-3) * ... Lc(1)*inv(Lr(1))*P(1)
*     Each of the elementary transformations is stored as a row or
*     column of the supplied array in this routine.
*     Multiplying a row eigenvector using BHAP1 results in
*       in an eigenvector of A.  
*
*     July 9, 2004 -- removed unused dummy variable W. 
*
*  Arguments:
*
*              N specifies the order of the matrix Z.
*              N must be nonnegative.
*              not modified.
*
*     Z       -double precision array, dimension (LDZ,N)
*              Z contains information about the reduction to
*              tridiagonal form.
*
*     LDZ     -integer
*              The leading dimension of the array Z.
*              LDZ >= max(1,N).
*
*     XRE     -double precision array, dimension (N)
*              The real part of the vector used in applying the
*              transformations.
*
*     XIM     -double precision array, dimension (N)
*              The imaginary part of the vector used in applying
*              the transformations.
*
*     IPIV    -integer array, dimension (N)
*              Contains the pivot sequence used during the reduction
*              to tridiagonal form.
*
*     KROW    -integer array, dimension (N)
*              Which row (if any) was eliminated after 
*              column elimination k.
*
*     IFLAG   -integer
*              You can put an error flag if you want.
*              Not currently used.  
*
*     .. Parameters ..
*
*     ------------------------------------------

      DOUBLE PRECISION TWO,ZERO
      PARAMETER (TWO = 2.0, ZERO = 0.0)

*     -------------------------------------------
*     ..
*     .. Local Scalars ..
*     -------------------------------------------

      INTEGER K, L, KP1
      DOUBLE PRECISION MRE, MIM, DDOT
      EXTERNAL DDOT

*     -------------------------------------------
*     ..
*     .. Executable Statements ..
*
*     Apply the permutations.
*     -----------------------------
 
      DO  K = N-2,1,-1
         KP1 = K+1
*        -------------------------------------------------------
*           Multiply by Lc(k).
*        -------------------------------------------------------

         XRE(KP1) = XRE(KP1) - DDOT(N-K-1,XRE(K+2),1,Z(K+2,K),1)
         XIM(KP1) = XIM(KP1) - DDOT(N-K-1,XIM(K+2),1,Z(K+2,K),1)

*        --------------------------
*           Multiply by inv(Lr(k)).
*        --------------------------

         IF ( KROW(K).NE.0 ) THEN

*          -------------------------------------------------------
*           The stride of LDZ in the DAXPY is poor for cache 
*           efficiency, but is a consequence of storing the 
*           multipliers in the space they were used to eliminate.  
*           If multiplications by inv(Z)are to be performed many 
*           times, consider packed  triangular storage for
*           upper triangular multipliers.  The packed version would
*           store the row elimination multipliers sequentially.  
*          --------------------------------------------------------

            CALL DAXPY(N-K-1,XRE(KP1),Z(KROW(K),K+2),LDZ,XRE(K+2),1)
            CALL DAXPY(N-K-1,XIM(KP1),Z(KROW(K),K+2),LDZ,XIM(K+2),1)

         ENDIF
         
      END DO 
 
      DO 70  K = N-2,1,-1
         KP1 = K+1
         L = IPVT(KP1)
         MRE = XRE(L)
         MIM = XIM(L)
         IF ( L .EQ. KP1 ) GO TO 70
         XRE(L) = XRE(KP1)
         XIM(L) = XIM(KP1)
         XRE(KP1) = MRE
         XIM(KP1) = MIM
 70   CONTINUE
*
      RETURN  
      END      ! End of BHAP1

*=======================================================================


      SUBROUTINE BHAP2(N,Z,LDZ,X,IPVT,KROW,IFLAG)

*     ----------------------------------------------
*
*     Global variables
*     ---------------------------------------

      INTEGER N, LDZ, IPVT(*), KROW(*), IFLAG
      DOUBLE PRECISION Z(LDZ,*), X(*)

*     ---------------------------------------
*
*     July 9, 2004 -- removed unused dummy variable W 
*
*  Purpose:
*     This routine is used to multiply the row vector X  
*     by the accumulated elementary transformations stored in
*     the array Z.  See routine BHAP1 for further comments about Z.
*
*  Arguments:
*     N       -integer
*              The number of rows and columns in the matrix Z.
*              N >= 0.
*
*     Z       -double precision array, dimension (LDZ,N)
*              Z contains information about the reduction to
*              tridiagonal form.
*
*     LDZ     -integer
*              The leading dimension of the array Z.
*              LDZ >= max(1,N).
*
*     X       -double precision array, dimension (N)
*              The vector used in applying the transformations.
*
*     IPIV    -integer array, dimension (N)
*              Contains the pivot sequence used during the
*              reduction to tridiagonal form.
*
*     KROW    -integer array, dimension (N)
*              KR(K) was the row eliminated after the kth column 
*              eliminated.  If no row was eliminated KR(K) = 0.  
*
*     IFLAG   -integer
*              Not currently used.  Could use to signal exceptions.  
*
*     .. Parameters ..
*     --------------------------------

      DOUBLE PRECISION TWO, ZERO
      PARAMETER (TWO = 2.0, ZERO = 0.0)

*     ---------------------------------
*     ..
*     .. Local Scalars ..
*     ---------------------------------

      INTEGER  K, L, KP1
      DOUBLE PRECISION M, DDOT
      EXTERNAL DDOT

*     ---------------------------------
*     ..
*     .. Executable Statements ..
*
*     Apply the permutations.
*     ---------------------------------

      DO 30  K = 1, N-2
         KP1 = K+1
         L = IPVT(KP1)
         M = X(L)
         IF (L .EQ. KP1) GO TO 30
         X(L) = X(KP1)
         X(KP1) = M
 30   CONTINUE 
      DO  K=1, N-2
         KP1 = K+1

*        -------------------------------------------------------
*          Postmultiply the row vector x by Lr(k).
*        -------------------------------------------------------

         IF(KROW(K).NE.0)THEN
            CALL DAXPY(N-K-1,-X(KP1),Z(KROW(K),K+2),LDZ,X(K+2),1) 
         ENDIF

*        --------------------------------------------------------
*          Postmultiply the row vector x by inv(Lc(k)).
*        --------------------------------------------------------

         X(KP1) = X(KP1) + DDOT(N-K-1,Z(K+2,K),1,X(K+2),1)

      END DO
*
      RETURN
      END    ! End of BHAP2

*===================================================================

      SUBROUTINE BHAP3(N,Z,LDZ,XRE,XIM,IPVT,KROW,IFLAG)

*     ----------------------------------------------------
*
*     Global Variables
*     ----------------------------------------------------

      INTEGER N, LDZ, IPVT(*), KROW(*), IFLAG
      DOUBLE PRECISION Z(LDZ,*), XRE(*), XIM(*)

*     ----------------------------------------------------
*
*     July 9, 2004 -- removed unused dummy variable W. 
*
*  Purpose:
*     This routine used to multiply the accumulated elementary
*     transformations stored in the array Z by the column vector
*     (XRE,XIM).  See the routine BHAP1 for comments about the
*      array Z.
*
*  Arguments:
*     N       -integer
*              The number of rows and columns in the matrix Z.
*              N >= 0.
*
*     Z       -double precision array, dimension (LDZ,N)
*              Z contains information about the reduction to
*              tridiagonal form.
*
*     LDZ     -integer
*              The leading dimension of the array Z.
*              LDZ >= max(1,N).
*
*     XRE     -double precision array, dimension (N)
*              The real part of the vector used in applying the
*              transformations.
*
*     XIM     -double precision array, dimension (N)
*              The imaginary part of the vector used in applying the
*              transformations.
*
*     IPIV    -integer array, dimension (N)
*              Contains the pivot sequence used during the reduction
*              to tridiagonal form.
*
*     KROW    -integer array, dimension (N)
*              KROW(K) was the row eliminated after the kth column 
*              eliminated.  If no row was eliminated KROW(K) = 0.  
*
*     IFLAG   -integer
*              Not currently used.  Can be used to throw
*              an exception.   
*
*     .. Parameters ..
*     ------------------------------------------
      DOUBLE PRECISION TWO,ZERO
      PARAMETER (TWO = 2.0, ZERO = 0.0)

*     ------------------------------------------
*     ..
*     .. Local Scalars ..
*     ------------------------------------------

      INTEGER K, L, KP1
      DOUBLE PRECISION MRE, MIM, DDOT
      EXTERNAL DDOT

*     -------------------------------------------
*
*     .. Executable Statements ..
*     -------------------------------------------

      DO  K = N-2, 1, -1

*        ----------------------------------------
*          Multiply by Lr(k).
*       -----------------------------------------

         KP1 = K+1

*        ----------------------------------------
*          Multiply by inv(Lc(k)).
*        ----------------------------------------

         CALL DAXPY(N-K-1,XRE(KP1),Z(K+2,K),1,XRE(K+2),1)
         CALL DAXPY(N-K-1,XIM(KP1),Z(K+2,K),1,XIM(K+2),1)
   
         IF(KROW(K).NE.0)THEN

*           ---------------------------------------------
*           A Transposed version of Z would avoid the stride
*               of LDZ here and be more cache efficient
*               but would require additional storage
*               A user who needs many eigenvectors
*               might consider allocating the storage 
*               space.  
*           ---------------------------------------------

            XRE(KP1) = XRE(KP1) - 
     +                 DDOT(N-K-1,Z(KROW(K),K+2),LDZ,XRE(K+2),1)
            XIM(KP1) = XIM(KP1) - 
     +                 DDOT(N-K-1,Z(KROW(K),K+2),LDZ,XIM(K+2),1)
         ENDIF
      END DO 

*     ---------------------------
*     Apply the permutations.
*     ---------------------------

      DO 40 K = N-2, 1, -1
         KP1 = K+1
         L = IPVT(KP1)
         MRE = XRE(L)
         MIM = XIM(L)
         IF (L .EQ. KP1) GO TO 40
         XRE(L) = XRE(KP1)
         XIM(L) = XIM(KP1)
         XRE(KP1) = MRE
         XIM(KP1) = MIM
 40   CONTINUE 
*
      RETURN  
      END     ! END OF BHAP3

*     --------------------------
*
*===============================================================


      SUBROUTINE BHAPC(N,Z,LDZ,XRE,XIM,IPVT,KROW,CONDIT,IFLAG)

*     -----------------------------------------------------------
*
      INTEGER N, LDZ, IPVT(*), KROW(*), IFLAG
      DOUBLE PRECISION Z(LDZ,*), XRE(*), XIM(*)

*     -----------------------------------------------------------
*
*     July 9, 2004 removed unused dummy variable W. 
*
*  Purpose:
*     This routine provides an estimator for condition of the
*     accumulated similarity transformations.  It chooses entries
*     of absolute value 1 for u to maximize the maximal entry of
*     Z*u.  
*
*     This is the LINPACK style of condition estimator, usually
*     considered to be more robust than the LAPACK style, but
*     in the average case less accurate.   
*
*  Arguments:
*     N       -integer
*              The number of rows and columns in the matrix Z.
*              N >= 0.
*
*     Z       -double precision array, dimension (LDZ,N)
*              Z contains information about the reduction to
*              tridiagonal form.
*
*     LDZ     -integer
*              The leading dimension of the array Z.
*              LDZ >= max(1,N).
*
*     XRE     -double precision array, dimension (N)
*              These vectors are set to one and minus one
*              as the computation proceeds. 
*
*     XIM     -double precision array, dimension (N)
*              Set to one or minus one as the computation
*              proceeds.  
*
*     IPIV    -integer array, dimension (N)
*              Contains the pivot sequence used during the reduction
*              to tridiagonal form.
*
*     KROW    -integer array, dimension (N)
*              KROW(K) was the row eliminated after the kth column 
*              eliminated.  If no row was eliminated KROW(K) = 0.  
*
*     IFLAG   -integer
*              Not currently used.  Can be used to throw
*              an exception.  
*
*     .. Parameters ..
*     ---------------------------------

      DOUBLE PRECISION TWO,ZERO
      PARAMETER (TWO = 2.0, ZERO = 0.0)

*     ---------------------------------
*     ..
*     .. Local Scalars ..
*     ---------------------------------

      INTEGER I, K, KP1
      DOUBLE PRECISION CONDIT, DDOT, DNRM2
      EXTERNAL DNRM2, DDOT

*     --------------------------------------------------
*
*     .. Executable Statements ..
*     --------------------------------------------------
      DO I = 1,N
         XRE(N) = 1.
         XIM(N) = - 1.
      END DO
      DO  K = N-2, 1, -1

*        --------------------
*
*          Multiply by Lr(k).
*        --------------------

         KP1 = K+1
         XRE(KP1) = 1.
         XIM(KP1) = -1.

*        ---------------------
*
*          Multiply by inv(Lc(k)).
*        --------------------

         CALL DAXPY(N-K-1,XRE(KP1),Z(K+2,K),1,XRE(K+2),1)
         CALL DAXPY(N-K-1,XIM(KP1),Z(K+2,K),1,XIM(K+2),1)

         IF(KROW(K).NE.0)THEN
            XRE(KP1) = XRE(KP1) - 
     +                    DDOT(N-K-1,Z(KROW(K),K+2),LDZ,XRE(K+2),1)
            XIM(KP1) = XIM(KP1) - 
     +                    DDOT(N-K-1,Z(KROW(K),K+2),LDZ,XIM(K+2),1)
         ENDIF
         IF (DNRM2(N-K,XRE(KP1),1).GT.DNRM2(N-K,XIM(KP1),1)) THEN 
            CALL DCOPY(N-K,XRE(KP1),1,XIM(KP1),1) 
         ELSE
            CALL DCOPY(N-K,XIM(KP1),1,XRE(KP1),1) 
         ENDIF
      END DO 

*     -----------------------------------
*
*     Apply the permutations.
*     -----------------------------------

      CONDIT = 0.
      DO  I = 1,N
         CONDIT = MAX(CONDIT,ABS(XRE(I)) )
      END DO 
      CONDIT = CONDIT*CONDIT
*
      RETURN
      END     ! End of BHAPC 

*     ------------------------------------
*
*=================================================

     
      SUBROUTINE BHRPCK (LDX,N,B,NU)

*     --------------------------
*
*       Array arguments
*     --------------------------

       DOUBLE PRECISION B(LDX,*)
       INTEGER NU(*)

*      -------------------------
*
*      Scalar arguments         
*      -------------------------

       INTEGER N, LDX

*      -----------------------------------------------------------
* 
*  Purpose:     
*       this subroutine wipes out the mutipliers in B, leaving only 
*       the banded Hessenberg form.
*
* Arguments
*
* LDX   -integer
*       leading dimension of B. (Parameter)
*       not modified.
*
* N     -integer
*       B is an N by N matrix.
*       not modified.
*
* B     -double precision array of dimension (LDX,N)
*       LDX should be greater than or equal to N.
*       On input B has multipliers packed into it.
*       On exit the entries of B corresponding to multipliers
*       have been zeroed.
*       modified.
*
* NU    -integer array of dimension N.
*       NU(I) is the number of nonzero entries of B to 
*       the right of the diagonal in row I.
*       not modified.
*
*      ----------------------------------------------------------
         
       INTEGER I,J
       DO  J = NU(1)+2,N
          B(1,J) = 0.D0
       END DO 
       DO  I = 2,N
          DO  J = 1,I-2
             B(I,J) = 0.D0
          END DO 
          DO  J = I+1+NU(I),N
             B(I,J) = 0.D0
          END DO 
      END DO 
* 
      RETURN
      END    ! End of BHRPCK

*     ---------------------------------------------------------
*
*==============================================================
   
      SUBROUTINE BHUNPK(LDX,B,BH,N,NU,NUMON,NUMAX)

*     --------------------------------------------
*
*       .. Scalar Arguments ..
*       ------------------------------------------

        INTEGER LDX, N, NUMAX

*       ----------------------
*
*       ..  Array Arguments ..  
*       ----------------------

        INTEGER NU(*), NUMON(*)
        DOUBLE PRECISION B(LDX,*), BH(LDX,*)

*       ------------------------------------
*
*  Purpose:
*       This subroutine extracts a copy of a banded Hessenberg matrix
*       from B and returns it in column storage in BH.
*       A copy of B is stored with subdiagonal in the first column
*       of BH, diagonal in second column, etc.  Extra zeros
*       fill out enough columns to allow storage of a backup
*       copy of BH if so desired.
*       Also returned is NUMON, which is the smallest integer array 
*       greater than or equal to NU.  It is useful in limiting the
*       range of operations for an LU decomposition of the banded Hessenberg
*       matrix.
*
*  Arguments:
*
*  LDX  - integer
*       leading dimension of B and BH
*       not modified.
*
*  B    - double precision array of dimension (LDX,N)
*       On entry B is an array of a banded Hessenberg matrix
*       and Gaussian multipliers used in reducing a general 
*       matrix to banded Hessenberg form
*       not modified.
*
*  BH   - double precision array of dimension (LDX,NUB)
*       On exit the first column of BH corresponds to
*       the subdiagonal of B, the second column of B to the diagonal,
*       etc.  
*       Sufficient subsequent columns of zeros of BH are provided
*       to allow a back-up copy of B to be stored.
*       modified.
*
*  N    - integer.
*       N is the number of rows and column of B and the number of
*       rows of BH
*       not modified.
*
*  NU   - integer array of dimension N.
*       NU(I) is the number of bands beyond the diagonal in the 
*       ith row of B.
*       not modified.
*
* NUMON - integer array of dimension N.
*       on return NUMON is the smallest monotone increasing array larger than
*       NU.  It is useful for constructing an LU decomposition of 
*       BH.
*       modified.
*
* NUMAX - integer
*       on return NUMAX is the largest element of NU, i.e.,
*       NUMAX+2 is the maximal bandwidth of B
*       modified.
*                       
*  Internal variables
*     ----------------

      INTEGER NUB, NUT, I, J

*     ----------------------------------------------------------------
*
* NUB   -integer
*       NUB is used in the loop producing zeros at the end of each row.
*
* NUT   -integer
*       NUT is used in the loop producing zeros at the end of each row.
*
*     ----------------------------------------------------------------
              
      NUMAX = NU(1)
      NUMON(1) = 1
      DO  I = 2,N
         NUMAX   =  MAX(NUMAX,NU(I))
         NUMON(I) = MAX(NUMON(I-1),NU(I))
      END DO 
      DO  I = N-NUMAX,N
         NUMON(I) = N-I
      END DO 
      NUT = NUMAX       +2
      DO  J = 1,1+NU(1)
         BH(1,J+1) = B(1,J)
      END DO 
      BH(1,1) = 0.0
      DO  J = NU(1)+3,NUT*2
         BH(1,J) = 0.0
      END DO 
      DO  I = 2,N
         NUB = MIN(I+NU(I),N)
         DO  J = I-1,NUB
            BH(I,J-(I-2)) = B(I,J)
         END DO 
         DO  J = NUB+1,NUT*2
            BH(I,J) = 0.0
         END DO 
      END DO
* 
      RETURN
      END    ! End of BHUNPK

*     -------------------------------
*
*=================================
