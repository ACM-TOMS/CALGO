*{SIGMA/TESTING/TestZGESVDQ.f}
      PROGRAM Z_C_GESVDQ_test
*
* SIGMA library, TESTING section updated February 2016. 
* Developed and coded by Zlatko Drmac, Department of Mathematics
* University of Zagreb, Croatia, drmac@math.hr
* Submitted to ACM TOMS 
*
********************************************************************************
*  Systematic stress testing of the complex versions of XGESVDQ.               *
*  The procedures under testing are CGESVDQ and ZGESVDQ.                       *
*  For details on the method and on the implementation see                     *
*  [1] Z. Drmac: A QR-preconditioned QR SVD method for computing the SVD       *
*                with high accuracy, TR-2016-HRZZ-9345-1,                      *
*                Department of Mathematics, University of Zagreb.              *
*                Submitted to ACM TOMS.                                        *
*  [2] Z. Drmac: xGESVDQ: A software implementation of the QR-preconditioned   *
*                QR SVD method for computing the singular value decomposition  *
*                TR-2016-HRZZ-9345-2, Department of Mathematics,               *
*                University of Zagreb. Submitted to ACM TOMS.                  *
*  The purpose of the test is to explore to what extent does the implementation*
*  obey the theoretical results of error analysis and perturbation theory as   * 
*  explained in [1].                                                           *
********************************************************************************
      IMPLICIT NONE 
      INTEGER     LDA, LDA2, LCWORK, LDU, LDV, LDVT, LDWORK, LIWORK
      PARAMETER ( LDA = 1001, LDA2 = 1001, LDU = LDA, 
     $            LDV = LDA2, LDVT = LDA2 )
      PARAMETER ( LCWORK = 5*LDA*LDA, LDWORK = 5*LDA2*LDA2, 
     $            LIWORK = 5 * LDA )
      INTEGER     LDCHKP
      PARAMETER   ( LDCHKP = 65000 )
*     .. 0 and 1 
      REAL                 ZERO,          ONE
      PARAMETER          ( ZERO = 0.0E0,  ONE = 1.0E0)
      DOUBLE PRECISION    DZERO,         DONE
      PARAMETER         ( DZERO = 0.0D0, DONE = 1.0D0 )
      COMPLEX         CZERO,                  CONE
      PARAMETER     ( CZERO = (0.0E0, 0.0E0), CONE = (1.0E0, 0.0E0) )   
      COMPLEX*16      ZZERO,                  ZONE
      PARAMETER     ( ZZERO = (0.0D0, 0.0D0), ZONE = (1.0D0, 0.0D0) )
*
*     .. Scalar Arguments ..
      LOGICAL            LSVEC, RSVEC
      CHARACTER          JOBA, JOBRWPVT, JOBT, JOBU, JOBUQ, JOBV, JOBVQ 
      INTEGER            dsc, ACCMODE, i, ICASE, ICASEPIV, ICASER, isc, 
     $                   IDIST, INFO, it, j,  M,  M100, MODE, N, SMODE, 
     $                   NEXP, NUMRANK, iexp, N1, PASSED, FAILED
      DOUBLE PRECISION   ANORM, COND, SCOND       
      DOUBLE PRECISION   ZGESVJ_R, ZGESVDQ_R    
      DOUBLE PRECISION   SVERR,  SVERR1, CGESVDQ_R 
      DOUBLE PRECISION   DMINV,  DMAXV,  OFFMAXV
      DOUBLE PRECISION   DMINU,  DMAXU,  OFFMAXU
      DOUBLE PRECISION   DMINVZ, DMAXVZ, OFFMAXVZ, SVECERR_U, SVECERR_V,
     $                   SVECERR_UZ, SVECERR_VZ,  SVERRZ, SVERRZ1     
      DOUBLE PRECISION   DMINUZ, DMAXUZ, OFFMAXUZ
*      
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   DSJ(LDA2), DSQ(LDA)
      REAL               SQ(LDA)
      COMPLEX*16     A( LDA, LDA2 ), A1(LDA,LDA2)
      COMPLEX*16     VT(LDVT,LDVT), ZU1( LDU, LDU ), ZVT(LDVT,LDVT) 
*
      COMPLEX*16     AJ(LDA,LDA2),  UJ(LDU,LDU), VJ(LDA2,LDA2)
*     ..       
      COMPLEX*16     ZA(LDA,LDA2), ZU(LDU,LDU), ZV(LDV,LDV)
      COMPLEX        CA(LDA,LDA2), CU(LDU,LDU), CV(LDV,LDV) 
      COMPLEX        CVT(LDVT,LDVT)
*     .. workspaces      
      COMPLEX*16     ZWORK(LCWORK)
      COMPLEX            CWORK(LCWORK)
      DOUBLE PRECISION   DWORK(LDWORK)
      REAL               RWORK(LDWORK)
      INTEGER            IWORK(LIWORK)
*     .. check points      
      DOUBLE PRECISION   CHKP1(LDCHKP,14)
      DOUBLE PRECISION   ZCHKP1(LDCHKP,14), ZCHKP2(LDCHKP,14)
      INTEGER            SJEME(LDCHKP,10)
*     .. dummies 
      DOUBLE PRECISION  DUMMY
      DOUBLE PRECISION  IVIV(1)
*     .. External Subroutines (SIGMA library)
*     .. new SVD subrutines under examination
      EXTERNAL         CGESVDQ, ZGESVDQ 
*     .. auxiliary subroutines
      EXTERNAL         ZCHKUN, ZCMPSV, ZMGEN, XGEC2Z, XGEZ2C  
*     .. External Functions (SIGMA library) ..            
      EXTERNAL         ZRESVD
      DOUBLE PRECISION ZRESVD
*     
*     .. External Subroutines (BLAS, LAPACK)...................................
      EXTERNAL         DLASET, ZDSCAL, ZLACPY, ZGESVD, SGESVJ 
*     .. External Functions (BLAS, LAPACK)......................................   
      DOUBLE PRECISION DZNRM2
      LOGICAL                                   LSAME      
      REAL                     SLAMCH
      DOUBLE PRECISION                 DLAMCH
      EXTERNAL         DZNRM2, SLAMCH, DLAMCH , LSAME 
*     .. External MKL
*      EXTERNAL         mkl_free_buffers
*
*..............................................................................
*     
      WRITE(*,*)' Test program for {C,Z}GESVDQ - implementation test'
      WRITE(*,*)' For A (M x N), M >= N, {C,Z}GESVDQ computes the SVD '
      WRITE(*,*)' A = U * S * V^H .'
      WRITE(*,*)' Data type: COMPLEX  and DOUBLE COMPLEX '
      WRITE(*,*)' Reference routine: DBLE COMPLEX Jacobi SVD, ZGESVJX.'
      WRITE(*,*)
      WRITE(*,*) '>> Job description ( U and/or V ) : '
      WRITE(*,*)
      WRITE(*,*) ' 1. SIGMA only       '
      WRITE(*,*) ' 2. SIGMA, U(MxN)    '
      WRITE(*,*) ' 3. SIGMA, V(NxN)    '
      WRITE(*,*) ' 4. SIGMA, U(MxN), V(NxN) '
      WRITE(*,*) ' 5. SIGMA, U(MxM), V(NxN) '
      WRITE(*,*) ' 6. SIGMA, U(MxM)    '
      WRITE(*,*) ' 7. SIGMA, U(MxNR), V(NxNR), '
      WRITE(*,*) '    (here NR denotes the numerical rank)'
      WRITE(*,*) ' 8. SIGMA, V(NxNR)   '
      WRITE(*,*) ' 9. SIGMA, U(NxNR)   '
      WRITE(*,713,ADVANCE="NO")
 713  FORMAT(' Enter the case, 1, 2, 3, 4, 5, 6, 7, 8 or 9 >> ')
      READ(*,*) ICASE
      IF ( ICASE .EQ. 1 ) THEN
            JOBU  = 'N'
            JOBV  = 'N'
            JOBUQ = 'N'
            JOBVQ = 'N'
      ELSE IF ( ICASE .EQ. 2 ) THEN
            JOBU  = 'U'
            JOBV  = 'N'
            JOBUQ = 'S'           
            JOBVQ = 'N'
      ELSE IF ( ICASE .EQ. 3 ) THEN
            JOBU  = 'N'
            JOBV  = 'V'
            JOBUQ = 'N'
            JOBVQ = 'A'
      ELSE IF ( ICASE .EQ. 4 ) THEN
            JOBU  = 'U'           
            JOBV  = 'V'
            JOBUQ = 'S'
            JOBVQ = 'A'
       ELSE IF ( ICASE .EQ. 5 ) THEN
            JOBU  = 'U'            
            JOBV  = 'V'
            JOBUQ = 'A'
            JOBVQ = 'A'
       ELSE IF ( ICASE .EQ. 6 ) THEN
            JOBU  = 'U'           
            JOBV  = 'N'
            JOBUQ = 'A'
            JOBVQ = 'N'
       ELSE IF ( ICASE .EQ. 7 ) THEN
            JOBU  = 'U'
            JOBV  = 'V'
            JOBUQ = 'R'
            JOBVQ = 'R'
       ELSE IF ( ICASE .EQ. 8 ) THEN
            JOBU  = 'N'
            JOBV  = 'V'
            JOBUQ = 'N'
            JOBVQ = 'R'
       ELSE IF ( ICASE .EQ. 9 ) THEN
            JOBU  = 'U'
            JOBV  = 'N'
            JOBUQ = 'R'
            JOBVQ = 'N'
       ELSE           
            WRITE(*,*) ' :( BAD OPTION! '
            STOP
      END IF
      LSVEC  = LSAME( JOBUQ, 'A' ) .OR. LSAME( JOBUQ, 'S' ) .OR. 
     $         LSAME( JOBUQ, 'R')
      RSVEC  = LSAME( JOBVQ, 'A' ) .OR. LSAME( JOBVQ, 'R' )
*      
      WRITE(*,*)'.  . ..'
      WRITE(*,*)'Choose accuracy level:                                '
      WRITE(*,*)'(1) The requested accuracy corresponds to the backward'
      WRITE(*,*)'    error bounded by                                  '
      WRITE(*,*)'   || delta A ||_F <= f(m,n) * EPS * || A ||_F,       '
      WRITE(*,*)'    where EPS = {S,D}LAMCH(Epsilon). This authorises  '
      WRITE(*,*)'    {C,Z}GESVDQ to truncate the computed triangular   '
      WRITE(*,*)'    factor in the rank revealing QR factorization     '
      WRITE(*,*)'    whenever the truncated part is below the          '
      WRITE(*,*)'    threshold EPS * ||A||_F.                          '
      WRITE(*,*)'(2) Similarly as with (1), but the truncation is more '
      WRITE(*,*)'    gentle: it is allowed only when there is a drop on'
      WRITE(*,*)'    the diagonal of the triangular factor.            '
      WRITE(*,*)'    This is a medium level of accuracy.               '
      WRITE(*,*)'(3) High accuracy requested. No numerical rank        '
      WRITE(*,*)'    determination based on the rank revealing QR      '
      WRITE(*,*)'    factorization is attempted.                       '
      WRITE(*,*)'(4) Accuracy level same as in (3), and in addition an '
      WRITE(*,*)'    estimate of the scaled cond. number is computed.  '
      WRITE(*,717,ADVANCE="NO")
717   FORMAT(' Choose  (1), (2),  (3) or (4) >> ')
      READ(*,*) ACCMODE
      IF ( ACCMODE.EQ. 1 ) THEN
                JOBA = 'A'
      ELSE IF ( ACCMODE .EQ. 2 ) THEN
                JOBA = 'M'
      ELSE IF ( ACCMODE .EQ. 3 ) THEN
                JOBA = 'H'
      ELSE IF ( ACCMODE .EQ. 4 ) THEN
                JOBA = 'E'
      ELSE          
               WRITE(*,*) ' :( Bad option!  :)'
               STOP
      END IF
*      
      WRITE(*,*)'.  .  . ..'
      WRITE(*,*)'Select the pivoting in the initial QR factorization:  '
      WRITE(*,*)'(1) Column-pivoted QR factorization as preconditioner.'
      WRITE(*,*)'(2) Row sorting prior to column pivoted QRF. '
      WRITE(*,*)'(This may give more accurate results if the rows of A '
      WRITE(*,*)' vary in length due to both ery large and very small  '
      WRITE(*,*)' weighting factors, causing large s_cond(A).          '
      WRITE(*,718,ADVANCE="NO")
 718  FORMAT(' Choose  (1) or (2) >> ')
      READ(*,*) ICASEPIV
      IF ( ICASEPIV.EQ. 1 ) THEN
                JOBRWPVT = 'N'
      ELSE IF ( ICASEPIV .EQ. 2 ) THEN
                JOBRWPVT = 'P'
      ELSE
               WRITE(*,*) ' :( Bad option!  :)'
               STOP
      END IF
      WRITE(*,*)'.    .   .  .  .  . ...'
      WRITE(*,*) 'After the QR factorization, proceed with computing:  '
      WRITE(*,*) '(1) The SVD of the upper triangular factor R         '
      WRITE(*,*) '    (This is recommended when using xGESVDQ.)        '
      WRITE(*,*) '(2) The SVD of the conjugate-transposed R**H         '
      WRITE(*,*) '    (This involves more data movement and it is left '
      WRITE(*,*) '     as optional for experimenting in the R&D.)      '
      WRITE(*,718)
      READ(*,*) ICASER
      IF ( ICASER.EQ. 1 ) THEN
                JOBT = 'N'
      ELSE IF ( ICASER .EQ. 2 ) THEN
                JOBT = 'T'
      ELSE
               WRITE(*,*) ' :( Bad option!  :)'
               STOP
      END IF
      WRITE(*,*) ' Maximal allowed size of the matrix is ', 
     $                LDA, ' x ', LDA2
      WRITE(*,*) 'The matrix must be square or tall rectangular.'
      WRITE(*,600,ADVANCE="NO")
 600  FORMAT('  >>>> Number of rows of A;    M = ')
      READ(*,*) M
      WRITE(*,601,ADVANCE="NO")
 601  FORMAT('  >>>> Number of columns of A; M >= N = ')
      READ(*,*) N
      IF ( ( M .GT. LDA ) .OR. ( N .GT. LDA2 ) ) THEN
         WRITE(*,*) ' :( The matrix is too big. '
         WRITE(*,*) ' Read the instructions. STOP.'
         STOP
      END IF
      IF ( ( M .LT. N ) .OR. ( M .LE. 0 ) .OR. ( N.LE. 0 ) ) THEN
           WRITE(*,*) ' :( The dimensions are invalid. STOP. '
           STOP
      END IF
      WRITE(*,*)        
      WRITE(*,603,ADVANCE="NO")
 603  FORMAT(' Number of tests in a specific class (1 or 2) = ')
      READ(*,*) NEXP
      IF ( ( NEXP .LT. 1 ) .OR. ( NEXP .GT. 2 ) ) THEN 
          WRITE(*,*) ':( ..bad input; bailing out.', 
     $               ' Keep calm and try again.'
          STOP
      END IF      
*     ................................................................. 
      it = 0       
      SCOND = 1.0D0
*     Initialize the seed for the psuedo-random number generator.
      ISEED(1) = 13
      ISEED(2) = 3
      ISEED(3) = 2002
      ISEED(4) = 2015     
*     .. go   go  go  go go ...
      DO 1203 isc = 1, 8
*       loop over scaled condition number from 1.0D1 to 1.0D8
      SCOND = SCOND * 1.0D1
*     .. now the scaled condition number is SCOND = 10**isc  
      COND = 1.0D1      
      DO 2103 dsc = 1, 8
          COND = COND * 1.0D2
      DO 123  MODE  = -6, 6 
*     .. try different modes of distribution of the initial diagonal          
      DO 1224 SMODE = -6, 6 
*     .. try different modes for additional diagonal scaling       
      DO 1225 IDIST = 1, 3    
      DO 1234 iexp = 1, NEXP
*     .. for each triple (SCOND, MODE, SMODE) try few random examples     
          WRITE(*,*) 'Grid point: ', isc, dsc, MODE, SMODE, IDIST, iexp 
      it = it + 1 
      DO 51 j = 1, 4
          SJEME(it,j) = ISEED(j)
 51   CONTINUE 
      SJEME(it,5)  = isc
      SJEME(it,6)  = dsc
      SJEME(it,7)  = MODE
      SJEME(it,8)  = SMODE
      SJEME(it,9)  = IDIST
      SJEME(it,10) = iexp
      WRITE(*,*) 'Test case # ', it 
*     .. generate the matrix  as DOUBLE COMPLEX      
      CALL ZMGEN( M, N, A, LDA, SCOND, COND, SMODE,    
     $     MODE, IDIST, ISEED, DWORK, ZWORK, INFO )
      WRITE(*,*) 'Test matrix successfully generated.'
      WRITE(*,*)'Its genetic material is saved', 
     $         ' for post festum analysis.'
*    
*     .. compute the scaled condition number using ZGESVD. It should
*     be close to the desired values given on input to ZMGEN     
      CALL ZLACPY( 'A', M, N, A, LDA, A1,  LDA )
      DO 122 j = 1, N 
         ANORM = DZNRM2( M, A1(1,j), 1 )
          IF (ANORM .NE. DZERO) CALL ZDSCAL(M, DONE/ANORM, A1(1,j),1)
122   CONTINUE        
      CALL ZGESVD( 'N', 'N', M, N, A1, LDA, DSJ, ZU1, LDU, 
     $     VT, LDVT, ZWORK, LCWORK, DWORK, INFO )   
      N1 = 0 
      DO 221 j = 1, N 
         IF ( DSJ(j) .GT. DZERO ) THEN
             N1 = N1 + 1
         ELSE
             GO TO 222
         END IF
 221  CONTINUE 
 222  CONTINUE 
      IF ( N1 .GT. 0 ) THEN 
         CHKP1(it,1)  = DSJ(1) / DSJ(N1)
         ZCHKP1(it,1) = DSJ(1) / DSJ(N1)
         ZCHKP2(it,1) = DSJ(1) / DSJ(N1) 
      ELSE
         CHKP1(it,1)  = -DONE
         ZCHKP1(it,1) = -DONE
         ZCHKP2(it,1) = -DONE
      END IF      
* 
*.......................................................................
* RUNNING THE DOUBLE COMPLEX ONE SIDED JACOBI >> ZGESVJ <<
* THIS IS THE REFERENCE ROUTINE 
*-----------------------------------------------------------------------     
       CALL ZLACPY( 'A', M, N, A, LDA, AJ,  LDA )
       WRITE(*,*) '>> Computing the SVD (ZGESVJX, reference values),...'
*>>>>>>.. SVD routine used for reference values        
       CALL ZGESVJX( 'G', JOBU, JOBV, M, N, AJ, LDA, DSJ, N, VJ,  
     $                   LDV, ZWORK, LCWORK, DWORK, LDWORK, INFO ) 
*>>>>>>
       WRITE(*,*) '   ZGESVJX ... done.'
*
       CALL ZLACPY( 'A', M, N, AJ, LDA, UJ, LDU )
*
*      After the call to a reference routine, the reference values are
*      UJ  :: M x N matrix of the  left singular vectors, DOUBLE COMPLEX
*      VJ  :: N x N matrix of the right singular vectors, DOUBLE COMPLEX  
*      DSJ :: N x 1 vector of the singular values,        DOUBLE PRECISION
*
*      Reference values are also tested: orthogonality of the singular 
*      vectors and the residual are routinely checked.
*
       WRITE(*,*) '>> Test of the singular vector matrices V, U :'
*
       IF ( RSVEC ) THEN 
           CALL ZCHKUN( 'C', N, N, VJ, LDV, DMINVZ, DMAXVZ,  
     $          OFFMAXVZ, ZWORK, N, INFO )
           WRITE(*,*)
           WRITE(*,1464) '   -> Max_i<>j (V^* * V)_ij  = ', OFFMAXVZ
           WRITE(*,1464) '   -> Max_i || V(:,i) ||     = ', DMAXVZ
           WRITE(*,1464) '   -> Min_i || V(:,i) ||     = ', DMINVZ
           WRITE(*,*) '..'
       ELSE
           WRITE(*,*)'The right sing. vector matrix V was not computed.'
           OFFMAXVZ = -DONE
           DMAXVZ   = -DONE
           DMINVZ   = -DONE
       END IF      
*      
       IF ( LSVEC ) THEN 
           CALL ZCHKUN( 'C', M, N, UJ, LDA, DMINUZ, DMAXUZ,  
     $                     OFFMAXUZ,  ZWORK, N, INFO )
           WRITE(*,*)
           WRITE(*,1464) '   -> Max_i<>j (U^* * U)_ij  = ', OFFMAXUZ
           WRITE(*,1464) '   -> Max_i || U(:,i) ||     = ', DMAXUZ
           WRITE(*,1464) '   -> Min_i || U(:,i) ||     = ', DMINUZ
           WRITE(*,*) '..' 
       ELSE
           WRITE(*,*) 'The left sing. vector matrix U was not computed.'
           OFFMAXUZ = -DONE
           DMAXUZ   = -DONE
           DMINUZ   = -DONE
       END IF
*          
       CALL ZLACPY( 'A', M, N, A, LDA, A1, LDA )
       IF ( LSVEC .AND. RSVEC ) THEN 
           ZGESVJ_R = ZRESVD( M, N, N, A1, LDA, AJ, LDA, VJ, LDV,
     $                 'N', DSJ, 'F', DWORK )       
           WRITE(*,1464) '  -> computed rezidual      = ', ZGESVJ_R    
       ELSE
           ZGESVJ_R = -DONE
       END IF   
*       
       ZCHKP1(it,2)  = - -DONE
       ZCHKP1(it,3)  = ZGESVJ_R
       ZCHKP1(it,4)  = -DONE
       ZCHKP1(it,5)  = OFFMAXUZ
       ZCHKP1(it,6)  = DMINUZ
       ZCHKP1(it,7)  = DMAXUZ
       ZCHKP1(it,8)  = OFFMAXVZ
       ZCHKP1(it,9)  = DMINVZ
       ZCHKP1(it,10) = DMAXVZ
       ZCHKP1(it,11) = -DONE
* 
************************************************************************
************************************************************************
*.......................................................................
*.......................................................................
* RUNNING THE PRECONDITIONED ZGESVD >> ZGESVDQ <<
* THIS IS THE ROUTINE UNDER EXAMINATION
*----------------------------------------------------------------------- 
*
       CALL ZLACPY( 'A', M, N, A, LDA, ZA,  LDA )
       WRITE(*,*) '>> Computing the SVD (ZGESVDQ), ...'
*>>>>>>.. SVD routine under examination (DOUBLE COMPLEX)     
       CALL ZGESVDQ( JOBA, JOBRWPVT, JOBT, JOBUQ, JOBVQ, M, N, ZA, LDA,
     $      DSQ, ZU, LDU, ZVT, LDV, NUMRANK, IWORK, 
     $      ZWORK, LCWORK, DWORK, LDWORK, INFO ) 
*>>>>>>
       DO 30 j = 1, N
         DO 40 i = 1, N 
            ZV(i,j) = CONJG(ZVT(j,i))
 40      CONTINUE
 30    CONTINUE               
       WRITE(*,*) '   ZGESVDQ ... done.'
       IF ( INFO .NE. 0 ) THEN
           WRITE(*,*) ':$ bad input; bailing out.'
           STOP
       END IF
*
*      After the call to the routine under examination, the values to be 
*      tested  are
*      ZU  :: M x {N,M} matrix of the  left singular vectors, DOUBLE COMPLEX
*      ZV  :: N x N matrix of the right singular vectors,     DOUBLE COMPLEX  
*             (possibly N x NUMRANK; always only N x NUMRANK part tested)
*      DSQ :: N x 1 vector of the singular values,            DOUBLE PRECISION
*       
       IF ( ACCMODE .EQ. 4 ) THEN 
           ZCHKP2(it,2) = DWORK(1) 
       ELSE
           ZCHKP2(it,2) = -DONE 
       END IF
*       
       IF ( RSVEC ) THEN 
           WRITE(*,*) '>> Test of the right singular vector matrix V:'
           IF ( LSAME(JOBVQ, 'R') ) THEN
               N1 = NUMRANK
           ELSE
               N1 = N
           END IF
*          .. test orthogonality
           CALL ZCHKUN( 'C', N, N1, ZV, LDV, DMINVZ, DMAXVZ,
     $          OFFMAXVZ,  ZWORK, N, INFO )
*          .. compare with reference values and compute max error*gap           
           CALL ZCMPSV( 'M', N,NUMRANK, VJ,LDV, DSJ, ZV,LDV, SVECERR_VZ,
     $          DUMMY, IVIV, ZWORK, INFO ) 
           WRITE(*,*)
           WRITE(*,1464) '   -> Max_i<>j (V^* * V)_ij  = ', OFFMAXVZ
           WRITE(*,1464) '   -> Max_i || V(:,i) ||     = ', DMAXVZ
           WRITE(*,1464) '   -> Min_i || V(:,i) ||     = ', DMINVZ
           WRITE(*,*) '..'
       ELSE
*          WRITE(*,*) 'The right singular vector matrix V was not computed.'
           OFFMAXVZ   = -DONE
           DMAXVZ     = -DONE
           DMINVZ     = -DONE
           SVECERR_VZ = -DONE
       END IF      
* 
       IF ( LSVEC ) THEN 
           WRITE(*,*) '>> Test of the left singular vector matrix U :' 
           IF ( LSAME(JOBUQ, 'A') ) THEN
               N1 = M
           ELSE IF ( LSAME(JOBUQ, 'S') ) THEN
               N1 = N
           ELSE
               N1 = NUMRANK
           END IF
*          .. test orthogonality
           CALL ZCHKUN( 'C', M, N1, ZU, LDU, DMINUZ, DMAXUZ, 
     $          OFFMAXUZ,  ZWORK, N, INFO )
*          .. compare with reference values and compute max error*gap            
           CALL ZCMPSV( 'M', M,NUMRANK, UJ,LDU, DSJ, ZU,LDU, SVECERR_UZ,
     $          DUMMY, IVIV, ZWORK, INFO )
           WRITE(*,*)
           WRITE(*,1464) '   -> Max_i<>j (U^* * U)_ij  = ', OFFMAXUZ
           WRITE(*,1464) '   -> Max_i || U(:,i) ||     = ', DMAXUZ
           WRITE(*,1464) '   -> Min_i || U(:,i) ||     = ', DMINUZ
       ELSE
*          WRITE(*,*) 'The left singular vector matrix U was not computed.'
           OFFMAXUZ   = -DONE
           DMAXUZ     = -DONE
           DMINUZ     = -DONE
           SVECERR_UZ = -DONE
       END IF
*           
       SVERRZ  = DZERO
       SVERRZ1 = DZERO       
       DO 12 j = 1, NUMRANK 
           SVERRZ  = MAX(SVERRZ, ABS((DSJ(j) -  DSQ(j)) / DSJ(j)))
           SVERRZ1 = MAX(SVERRZ1,ABS((DSJ(j) -  DSQ(j)) / DSJ(1))) 
 12     CONTINUE  
       WRITE(*,*) 'Maximal relative error in the singular values was ', 
     $            SVERRZ
       WRITE(*,*) 'Maximal error in the singular values relative to ', 
     $            ' ||A||_2 was ', SVERRZ1
*       
        IF ( LSVEC .AND. RSVEC ) THEN 
*          .. compute the residual in double complex precision           
           CALL ZLACPY( 'A', M, N, A, LDA, A1, LDA ) 
           ZGESVDQ_R =  ZRESVD( M, N, NUMRANK, A1, LDA, ZU, LDU, ZVT, 
     $                 LDVT, 'C', DSQ, 'F', DWORK ) 
           WRITE(*,1464) '  -> computed rezidual      =', ZGESVDQ_R
       ELSE
           ZGESVDQ_R = -DONE
       END IF       
*          
       ZCHKP2(it,3)  = ZGESVDQ_R
       ZCHKP2(it,4)  = SVERRZ
       ZCHKP2(it,5)  = OFFMAXUZ
       ZCHKP2(it,6)  = DMINUZ
       ZCHKP2(it,7)  = DMAXUZ
       ZCHKP2(it,8)  = OFFMAXVZ
       ZCHKP2(it,9)  = DMINVZ
       ZCHKP2(it,10) = DMAXVZ
       ZCHKP2(it,11) = SVECERR_VZ
       ZCHKP2(it,12) = SVECERR_UZ
       ZCHKP2(it,13) = SVERRZ1 
       ZCHKP2(it,14) = NUMRANK    
       
************************************************************************
************************************************************************
*.......................................................................
*.......................................................................
* RUNNING THE PRECONDITIONED CGESVD >> CGESVDQ <<
* THIS IS THE ROUTINE UNDER EXAMINATION
*----------------------------------------------------------------------- 
*
       CALL XGEZ2C( M, N, A, LDA, CA, LDA )
       WRITE(*,*) '>> Computing the SVD (CGESVDQ), ...'
*>>>>>>.. SVD routine under examination (COMPLEX)  
       CALL CGESVDQ( JOBA, JOBRWPVT, JOBT, JOBUQ, JOBVQ, M, N, CA, LDA,
     $      SQ, CU, LDU, CVT, LDV, NUMRANK, IWORK, 
     $      CWORK, LCWORK, RWORK, LDWORK, INFO ) 
*>>>>>>
       DO 50 j = 1, N
         DO 60 i = 1, N 
            CV(i,j) = CONJG(CVT(j,i))
 60      CONTINUE
 50    CONTINUE
       WRITE(*,*) '   CGESVDQ ... done.'
       IF ( INFO .NE. 0 ) THEN
           WRITE(*,*) ':$ bad input; bailing out.'
           STOP
       END IF
*
*      After the call to the routine under examination, the values to be 
*      tested  are
*      CU  :: M x {N,M} matrix of the  left singular vectors, COMPLEX
*      CV  :: N x N matrix of the right singular vectors,     COMPLEX  
*             (possibly N x NUMRANK; always only N x NUMRANK part tested)
*      SQ :: N x 1 vector of the singular values,             REAL
*       
       IF ( ACCMODE .EQ. 4 ) THEN 
           CHKP1(it,2) = RWORK(1) 
       ELSE
           CHKP1(it,2) = -DONE 
       END IF
*       
       IF ( RSVEC ) THEN 
           WRITE(*,*) '>> Test of the right singular vector matrix V:'
           IF ( LSAME(JOBVQ, 'R') ) THEN
               N1 = NUMRANK
           ELSE
               N1 = N
           END IF
           CALL XGEC2Z( N, N1, CV, LDV, VT, LDVT )
*          .. test orthogonality
           CALL ZCHKUN( 'C', N, N1, VT, LDVT, DMINV, DMAXV,
     $          OFFMAXV, ZWORK, N, INFO )
*          .. compare with reference values and compute max error*gap            
           CALL ZCMPSV( 'M', N,NUMRANK, VJ,LDV, DSJ, VT,LDVT, SVECERR_V,
     $          DUMMY, IVIV, ZWORK, INFO )
           WRITE(*,*)
           WRITE(*,1463) '   -> Max_i<>j (V^* * V)_ij  = ', OFFMAXV
           WRITE(*,1463) '   -> Max_i || V(:,i) ||     = ', DMAXV
           WRITE(*,1463) '   -> Min_i || V(:,i) ||     = ', DMINV
           WRITE(*,*) '..'
       ELSE
*          WRITE(*,*) 'The right singular vector matrix V was not computed.'
           OFFMAXV = -DONE
           DMAXV   = -DONE
           DMINV   = -DONE
           SVECERR_V = -DONE
       END IF      
* 
       IF ( LSVEC ) THEN 
           WRITE(*,*) '>> Test of the left singular vector matrix U :' 
           IF ( LSAME(JOBUQ, 'A') ) THEN
               N1 = M
           ELSE IF ( LSAME(JOBUQ, 'S') ) THEN
               N1 = N
           ELSE
               N1 = NUMRANK
           END IF
           CALL XGEC2Z( M, N1,  CU, LDU, ZU1, LDU )
*          .. test orthogonality
           CALL ZCHKUN( 'C', M, N1, ZU1, LDU, DMINU, DMAXU, 
     $          OFFMAXU, ZWORK, N, INFO )
*          .. compare with reference values and compute max error*gap            
           CALL ZCMPSV( 'M', M,NUMRANK, UJ,LDU, DSJ, ZU1,LDU, SVECERR_U,
     $          DUMMY, IVIV, ZWORK, INFO )
           WRITE(*,*)
           WRITE(*,1463) '   -> Max_i<>j (U^* * U)_ij  = ', OFFMAXU
           WRITE(*,1463) '   -> Max_i || U(:,i) ||     = ', DMAXU
           WRITE(*,1463) '   -> Min_i || U(:,i) ||     = ', DMINU
       ELSE
*          WRITE(*,*) 'The left singular vector matrix U was not computed.'
           OFFMAXU = -ONE
           DMAXU   = -ONE
           DMINU   = -ONE
           SVECERR_U = -DONE
       END IF
*           
       SVERR  = DZERO
       SVERR1 = DZERO
       CALL DLASET('G', N , 1 , 0.0D0, 0.0D0, DWORK(1), N )
       DO 102 j = 1, NUMRANK 
           DWORK(j) = DBLE(SQ(j))
           SVERR  = MAX(SVERR, ABS(((DSJ(j) -  DBLE(SQ(j))) / DSJ(j))))
           SVERR1 = MAX(SVERR1,ABS(((DSJ(j) -  DBLE(SQ(j))) / DSJ(1)))) 
 102     CONTINUE  
       WRITE(*,*) 'Maximal relative error in the singular values was ', 
     $            SVERR
       WRITE(*,*) 'Maximal error in the singular values relative to ', 
     $            ' ||A||_2 was ', SVERR1
       
        IF ( LSVEC .AND. RSVEC ) THEN 
*          .. compute the residual in double complex precision           
           CALL XGEC2Z( M, N, CU, LDU, ZU1, LDU )
           CALL XGEC2Z( N, N, CVT, LDV, VT, LDV )
           CALL ZLACPY( 'A', M, N, A, LDA, A1, LDA )
           CGESVDQ_R = ZRESVD( M,N, NUMRANK, A1, LDA, ZU1,LDU, VT,LDV,
     $                 'C', DWORK, 'F', DWORK(N+1) )
           WRITE(*,1463) '  -> computed rezidual      =', CGESVDQ_R
       ELSE
           CGESVDQ_R = -DONE
       END IF       
*          
       CHKP1(it,3)  = CGESVDQ_R
       CHKP1(it,4)  = SVERR
       CHKP1(it,5)  = OFFMAXU
       CHKP1(it,6)  = DMINU
       CHKP1(it,7)  = DMAXU
       CHKP1(it,8)  = OFFMAXV
       CHKP1(it,9)  = DMINV
       CHKP1(it,10) = DMAXV
       CHKP1(it,11) = SVECERR_V
       CHKP1(it,12) = SVECERR_U
       CHKP1(it,13) = SVERR1 
       CHKP1(it,14) = NUMRANK
*       
*       CALL mkl_free_buffers
*       
 1234  CONTINUE  
 1225  CONTINUE      
 1224  CONTINUE           
 123   CONTINUE  
 2103  CONTINUE      
 1203  CONTINUE    
*      
 1462   FORMAT( A30, 1X, D15.8, 1X, E15.8 )       
 1463   FORMAT( A30,1X,D15.8 )       
 1464   FORMAT( A30,1X,D24.16 )  
 1465   FORMAT( A30, 1X, D24.16, 1X, D24.16 )  
 1466   FORMAT(E16.8, 3X, D24.16)         
*
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' ===================================================='
      WRITE(*,*) '| Grading the test-quick grade (check the graphs of |'
      WRITE(*,*) '| all measured errors (saved) for more details)     |'
      WRITE(*,*) '|===================================================='
      WRITE(*,*)
      WRITE(*,*)
      M100 = MAX(M,100)
      PASSED = 0 
      FAILED = 0 
*     .. residuals       
      IF ( LSVEC .AND. RSVEC ) THEN 
         WRITE(*,*) 'Residuals (||A * V - U*SIGMA||_F / ||A||_F):'    
         CGESVDQ_R  = DZERO
         ZGESVJ_R   = DZERO
         ZGESVDQ_R  = DZERO
         DO 700 j = 1, it
            ZGESVDQ_R = MAX(ZGESVDQ_R, ZCHKP2(j,3))               
            CGESVDQ_R = MAX(CGESVDQ_R,  CHKP1(j,3)) 
            ZGESVJ_R  = MAX( ZGESVJ_R, ZCHKP1(j,3))             
700      CONTINUE            
         WRITE(*,*) 'Maximal CGESVDQ residual was ', CGESVDQ_R
         WRITE(*,*) 'Maximal ZGESVDQ residual was ', ZGESVDQ_R
         IF ( (CGESVDQ_R  .LE. SLAMCH('E')*M100*SQRT(REAL(N))) .AND.
     $        (ZGESVDQ_R .LE. DLAMCH('E')*M100*SQRT(DBLE(N))) ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1 
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            WRITE(*,*) 'Check the graph of all residuals! Analyze!'
            FAILED = FAILED + 1 
         END IF     
         WRITE(*,*) 'Residual check for the reference routine.'
         WRITE(*,*) 'Maximal ZGESVJX (reference routine) residual was',
     $              ZGESVJ_R
         IF ( ZGESVJ_R .LE. DLAMCH('E')*M100*SQRT(FLOAT(N)) ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            WRITE(*,*) '<!!> Check the reference routine!'
            FAILED = FAILED + 1 
         END IF     
      END IF 
*     .. singular values, errors relative to ||A||_F        
      WRITE(*,*) 'Errors in the singular values:'
      SVERR1  = DZERO
      SVERRZ1 = DZERO
      DO 760 j = 1, it
         SVERR1   = MAX( SVERR1,   CHKP1(j,13))  
         SVERRZ1  = MAX( SVERRZ1, ZCHKP2(j,13))
760   CONTINUE    
      WRITE(*,*) 'Testing the quotient: max|d_sigma_i|/||A||_F'
      WRITE(*,*) '(Should be below some EPS x modest f(M,N); here use'
      WRITE(*,*) ' empirical value EPS x max(M,20). '
      WRITE(*,*) 'The maximal measured max|d_sigma_i|/||A||_F was ', 
     $            SVERR1, ' (CGESVDQ) and ', SVERRZ1, ' (ZGESVDQ)'
      
      IF ( (SVERR1 .LE. SLAMCH('E')*M100) .AND. 
     $     (SVERRZ1 .LE. DLAMCH('E')*M100) ) THEN
         WRITE(*,*) '..................................', ':) passed.'
         PASSED = PASSED + 1 
      ELSE
         WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ':( FAILED.'
         FAILED = FAILED + 1 
      END IF
*       
      IF ( ACCMODE .GE. 3 ) THEN 
         WRITE(*,*) 'Testing the quotient:', 
     $              ' max (|d_sigma_i|/sigma_i)/(EPS*SCOND)'
         WRITE(*,*) '(','Ideally, it is smaller than Max(M,20) = ', 
     $              M100, ';'
         WRITE(*,*) 'this empirical value is below the theoretical', 
     $              ' upper bound which'
         WRITE(*,*) 'is a modestly growing f(M,N) <= low degree', 
     $    ' polynomial in M, N.)'
         SVERR   = ZERO
         SVERRZ  = DZERO
         DO 770 j = 1, it
            SVERR  = MAX(SVERR,   CHKP1(j,4)/SLAMCH('E')/ CHKP1(j,1)) 
            SVERRZ = MAX(SVERRZ, ZCHKP2(j,4)/DLAMCH('E')/ZCHKP2(j,1))
770      CONTINUE           
         WRITE(*,*) 'Maximal CGESVDQ quotient', 
     $              ' (|d_sigma_i|/sigma_i)/(EPS*SCOND) was ', SVERR
         WRITE(*,*) 'Maximal ZGESVDQ quotient', 
     $              ' (|d_sigma_i|/sigma_i)/(EPS*SCOND) was ', SVERRZ
         IF ( (SVERR .LE. M100) .AND. ( SVERRZ .LE. M100 ) ) THEN
            WRITE(*,*) '...............................', ':) passed.'
            PASSED = PASSED + 1 
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ':( FAILED.'
            WRITE(*,*) '<!!> Check all measured errors.', 
     $               ' The failure could have'
            WRITE(*,*) 'occured at an ill-conditioned matrix where', 
     $                 ' small relative'
            WRITE(*,*) 'errors in smallest singular values were not', 
     $                 ' warranted.'
            WRITE(*,*) '<!!> Check the graph of all errors in the', 
     $                 ' sigmas.'
            FAILED = FAILED + 1 
         END IF   
      END IF
*               
      IF ( LSVEC ) THEN
       WRITE(*,*) 'Orthogonality of the left singular vectors'
       OFFMAXU   = DZERO
       OFFMAXUZ  = DZERO
       OFFMAXVZ  = DZERO
       DO 800 j = 1, it
          OFFMAXU  = MAX( OFFMAXU, CHKP1(j,5)  )
          OFFMAXUZ = MAX(OFFMAXUZ, ZCHKP1(j,5) )
          OFFMAXVZ = MAX(OFFMAXVZ, ZCHKP2(j,5) )
800    CONTINUE           
       WRITE(*,*) 'Maximal |u_i**u_j|_i.NE.j in  CGESVDQ  was ', 
     $            OFFMAXU
       IF ( OFFMAXU .LE. SLAMCH('E')*M*1.0E1 ) THEN
          WRITE(*,*) '..............................', ' :) passed.'
          PASSED = PASSED + 1 
       ELSE
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
          FAILED = FAILED + 1 
       END IF        
       WRITE(*,*) 'Maximal |u_i**u_j|_i.NE.j in ZGESVDQ  was ', 
     $            OFFMAXVZ
       IF ( OFFMAXVZ .LE. DLAMCH('E')*M*1.0D1 ) THEN
          WRITE(*,*) '..............................', ' :) passed.'
          PASSED = PASSED + 1
       ELSE
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
          FAILED = FAILED + 1
       END IF
         WRITE(*,*) 'Maximal |u_i**u_j|_i.NE.j in ZGESVJX  was ', 
     $              OFFMAXUZ
         IF ( OFFMAXUZ .LE. DLAMCH('E')*M*1.0D1 ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1 
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            FAILED = FAILED + 1
         END IF                  
      END IF    
*
      IF ( RSVEC ) THEN
         WRITE(*,*) 'Orthogonality of the right singular vectors'
         OFFMAXV   = DZERO
         OFFMAXVZ  = DZERO
         OFFMAXUZ  = DZERO
         DO 880 j = 1, it
            OFFMAXV  = MAX(  OFFMAXV, CHKP1(j,8)  )
            OFFMAXVZ = MAX( OFFMAXVZ, ZCHKP1(j,8) )
            OFFMAXUZ = MAX( OFFMAXUZ, ZCHKP2(j,8) )
880      CONTINUE      
         WRITE(*,*) 'Maximal |v_i**v_j|_i.NE.j in  CGESVDQ  was ', 
     $               OFFMAXV
         IF ( OFFMAXV .LE. SLAMCH('E')*M*1.0E1 ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1 
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            FAILED = FAILED + 1 
         END IF   
         WRITE(*,*) 'Maximal |v_i**v_j|_i.NE.j in ZGESVDQ  was ', 
     $               OFFMAXUZ
         IF ( OFFMAXUZ .LE. DLAMCH('E')*M*1.0D1 ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            FAILED = FAILED + 1
         END IF   
         WRITE(*,*) 'Maximal |v_i**v_j|_i.NE.j in ZGESVJX  was ', 
     $               OFFMAXVZ
         IF ( OFFMAXVZ .LE. DLAMCH('E')*M*1.0D1 ) THEN
            WRITE(*,*) '..............................', ' :) passed.'
            PASSED = PASSED + 1 
         ELSE
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ' :( FAILED.'
            FAILED = FAILED + 1
         END IF                    
      END IF       
*
      IF ( FAILED .EQ. 0 ) THEN 
          WRITE(*,*) '==============================================='
          WRITE(*,*) '          PASSED ALL TESTS                  '
          WRITE(*,*) '==============================================='
      ELSE
          WRITE(*,*) '==============================================='
          WRITE(*,*) ' PASSED ', PASSED, ' TESTS'
          WRITE(*,*) ' FAILED ', FAILED, ' TESTS'
          WRITE(*,*) '==============================================='
      END IF
      
      WRITE(*,*) 'Writing the test results to the files:             '
      WRITE(*,*) '1. CGESVDQ.table', 
     $           ' (all measured errros for CGESVDQ)     '
      WRITE(*,*) '2. ZGESVDQ.table', 
     $           ' (all measured errros for ZGESVDQ)     '
      WRITE(*,*) '3. ZGESVJ.table', 
     $           ' (basic check of the reference routine)'
      WRITE(*,*) '4. dimensions.data', 
     $  ' (the dimension of the test matices and test parameters) ' 
      WRITE(*,*) '5. SJEME.all',
     $  ' (SEED values - genes for all test matrices) ' 
*        
      OPEN (33, FILE = 'CGESVDQ.table')
      DO 2000 i = 1, it
          WRITE(33,3333) (CHKP1(i,j), j=1, 14)
 2000 CONTINUE
      CLOSE( 33 )
*        
      OPEN (37, FILE = 'ZGESVDQ.table')
      DO 2002 i = 1, it
         WRITE(37,5555) (ZCHKP2(i,j), j=1, 14)
 2002 CONTINUE
      CLOSE( 37 )   
*        
      OPEN (36, FILE = 'ZGESVJ.table')
      DO 5000 i = 1, it
         WRITE(36,5555) ( ZCHKP1(i,j), j=1, 14 )
 5000 CONTINUE
      CLOSE( 36 )  
*        
      OPEN( 34, FILE = 'dimensions.data')
      WRITE(34,3455) M, N, ICASE, ACCMODE, ICASEPIV, ICASER, NEXP 
      CLOSE(34)
*        
      OPEN (39, FILE = 'SJEME.all')
      DO 5050 i = 1, it
         WRITE(39,3945) ( SJEME(i,j), j=1, 10 )
 5050 CONTINUE
      CLOSE( 39 ) 
      WRITE(*,*) '... done.'
 3333 FORMAT( 14E18.8 ) 
 5555 FORMAT( 14D24.16)    
 3334 FORMAT( 2E18.8 )  
 3455 FORMAT( 2I5, 5I4 )
 3945 FORMAT( 4I7, 6I5)     
*        
      END PROGRAM Z_C_GESVDQ_test

