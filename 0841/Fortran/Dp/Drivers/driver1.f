      PROGRAM DRIVER
*-------------------------------------------------------
*     Update History:
*     July 9,  2004 by G. Howell
*     Nov. 21, 2002 by N. Diaa
*     Nov. 20, 2002 by N. Diaa
*-----------------------------------------------------------------
*
*     matrix entries are read from the file data
*
*-----------------------------------------------------------------
C     .. Parameters ..
      INTEGER LDX
      PARAMETER (LDX=21)
C     ..
C     .. Scalars in Common ..
      INTEGER SEED
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONDIT,SIZE2,TOL,XNORM1,XNORM2,XNORM3
      INTEGER AVGNU,I,IFLAG,INFO,J,M,MAXNU,NUMAX,SIZE
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BH(LDX,LDX),DUMV(LDX),MAT0(LDX,LDX),XIM(LDX),
     +                 XRE(LDX)
      INTEGER BHPIVT(LDX),KROW(LDX),KRVECT(LDX),NU(LDX),NUMON(LDX)
C     ..
C     .. External Subroutines ..
      EXTERNAL BHAP1,BHAP2,BHAP3,BHAPC,BHESS,BHRPCK,BHUNPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,INT,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON SEED
C     ..
      OPEN (33,FILE='data',STATUS='UNKNOWN',ERR=30)
      READ (33,FMT=*) TOL
      READ (33,FMT=*) SIZE2
      SIZE = INT(SIZE2)
      IF ((SIZE.GT.LDX) .OR. (SIZE.LE.2) .OR. (TOL.LE.0)) THEN
          WRITE(*,*)' ERROR: Invalid data in data '
          WRITE(*,*)
     +       '       Real number Tol or integer Size in data '
          WRITE(*,*)'        Tol must be real number > 0 '
          WRITE(*,*)
     +       '        Size must be integer: 2 < size <= ',LDX,'  '
          CLOSE (33)
          GO TO 30

      END IF

      DO I = 1,SIZE
          READ (33,FMT=*) (MAT0(I,J),J=1,SIZE)
      END DO
      CLOSE (33)

*--------------------------------------------------------------
*
*     A good general policy is to call a balancing routine
*     (using diagonal similarity transformations)
*     before BHESS (or any reduction to Hessenberg form).
*     A standard routine is balanc.f from EISPACK.
*
*--------------------------------------------------------------
*     Now we have the matrix MAT0,
*     either read from data or randomly generated
*--------------------------------------------------------------
   10 CONTINUE

      CALL BHESS(LDX,MAT0,SIZE,BHPIVT,NU,KRVECT,KROW,DUMV,TOL,INFO)

      IF (INFO.NE.0) GO TO 30
   20 CONTINUE

*--------------------------------------------------------------
*
*     The average bandwidth is computed as AVGNU
*
*--------------------------------------------------------------

      AVGNU = 0
      MAXNU = 0

      DO I = 1,SIZE
          AVGNU = AVGNU + NU(I)
          MAXNU = MAX(MAXNU,NU(I))
      END DO
      AVGNU = (AVGNU+SIZE/2)/SIZE
      WRITE(*,*)' Maximal upper bandwidth of returned matrix=',
     +          MAXNU
*----------------------------------------------------------------
*
*      Test BHAP1
*
      WRITE(*,*)' Testing BHAP1'
      DO I = 1,SIZE
          XRE(I) = 1.D0
          XIM(I) = 0.D0
      END DO
      CALL BHAP1(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,IFLAG)

      WRITE(*,*)'     XRE(I)      I  BHPIVT(I)   KROW(I)'
      XNORM1 = 0.0D0
      DO I = 1,SIZE
          WRITE(*,'(1x,e14.7, 2i5,7x,i5)')
     +            XRE(I),I,BHPIVT(I),KROW(I)
          XNORM1 = MAX(XNORM1,ABS(XRE(I)))
      END DO
      WRITE(*,*)' The Max Norm is = ',XNORM1
      WRITE(*,*)' XRE(2)= ',XRE(2)
      WRITE(*,*)' XRE(SIZE/2)= ',XRE(SIZE/2)
      WRITE(*,*)'                 '


*----------------------------------------------------------------
*
*      Test BHAP2
*
      WRITE(*,*)' Testing BHAP2'
      DO I = 1,SIZE
          XRE(I) = 1.D0
      END DO
      CALL BHAP2(SIZE,MAT0,LDX,XRE,BHPIVT,KROW,IFLAG)

      XNORM2 = 0.D0
      WRITE(*,*)'     XRE(I)      I  BHPIVT(I)   KROW(I)'
      DO I = 1,SIZE
          WRITE(*,'(1x,e14.7, 2i5,7x,i5)')
     +            XRE(I),I,BHPIVT(I),KROW(I)
          XNORM2 = MAX(XNORM2,ABS(XRE(I)))
      END DO
      WRITE(*,*)'                 '
      WRITE(*,*)' The Max Norm is  ',XNORM2
      WRITE(*,*)'               '
      WRITE(*,*)'               '

*----------------------------------------------------------------
*
*      Test BHAP3
*
*----------------------------------------------------------------

      WRITE(*,*) ' Testing BHAP3'
      DO I = 1,SIZE
          XRE(I) = 1.D0
          XIM(I) = 0.D0
      END DO
      CALL BHAP3(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,IFLAG)
      XNORM3 = 0.D0
      WRITE(*,*)'     XRE(I)      I  BHPIVT(I)   KROW(I)'
      DO I = 1,SIZE
          WRITE(*,'(1x,e14.7, 2i5,7x,i5)')
     +            XRE(I),I,BHPIVT(I),KROW(I)
          XNORM3 = MAX(XNORM3,ABS(XRE(I)))
      END DO
      WRITE(*,*)'The Max Norm is ',XNORM3
      WRITE(*,*)' XRE(2)= ',XRE(2)
      WRITE(*,*)' XRE(SIZE/2)= ',XRE(SIZE/2)
*     ---------------------------------------------------------------
*
*     The next statements call BHAPC -- a LINPACK style
*     estimator of the infinity norm of inv(Z)
*
*     ---------------------------------------------------------------
      WRITE(*,'(// ''Call to BHAPC'')')
      CALL BHAPC(SIZE,MAT0,LDX,XRE,XIM,BHPIVT,KROW,CONDIT,IFLAG)
      WRITE(*,*)
      WRITE(*,*)' Y vector'
      DO I = 1,SIZE
          WRITE(*,'(1x,a,i5,a,e14.7)')'Y(',I,')=',XRE(I)
      END DO
      WRITE(*,*)'                 '
      WRITE(*,*)' ESTIMATED CONDITION OF ACCUMULATED'
      WRITE(*,*)' TRANSFORMATIONS =',CONDIT
      WRITE(*,*)'                 '


*     ---------------------------------------------------------------
*
*      Test BHUNPK  or  BHRPCK
*

      WRITE(*,*) ' Testing BHUNPK'
      CALL BHUNPK(LDX,MAT0,BH,SIZE,NU,NUMON,NUMAX)
      WRITE(*,*) ' Relevant part of BH array'
      WRITE(*,*) 
      WRITE(*,*)'    I    J    BH(I,J)'
      DO I = 1,SIZE
          WRITE(*,'(1x,2i5,2x,e14.7)') (I,J,BH(I,J),J=1,NU(I)+2)
      END DO



*        --------------------------------------------
*         zero out multipliers in atobhes matrix MAT0
*        --------------------------------------------

      WRITE(*,*) ' Testing BHRPCK'
      CALL BHRPCK(LDX,SIZE,MAT0,NU)
      WRITE(*,*) ' Relevant part of MAT0 array'
      WRITE(*,*) 
      WRITE(*,*)'    I    J    MAT0(I,J)'
      DO I = 1,SIZE
          M = MIN(SIZE,NU(I)+2)
          DO J = 1,M
              IF (I+J-2.NE.0) THEN
                  WRITE(*,'(1x,2i5,2x,e14.7)')
     +                     I,I + J - 2,MAT0(I,I+J-2)
              END IF

          END DO
      END DO

   30 CONTINUE
      END

*     -----------------------------------------------
*
*====================================================
      SUBROUTINE BHBACK(LDX,N,A,B,W,V,TOL,NU,PIV,KRVECT,KROW)

*  This subroutine tests backward error of BHESS
*    by running BHESS to get a Hessenberg matrix
*    then using the stored multipliers to regenerate
*    the original matrix.
*
*  Outputs are the difference between the original matrix
*    and the reconstruction (computed backward error)
*
*  The upper bandwidth of the returned matrix and the estimated
*    condition number of the similariy transformation are
*    also returned.
*
*   Arguments
*-------------------
*
*     N -- Integer  Number of rows and columns of input matrices
*
*     LDX -- Integer Leading Dimension of Variable dimension arrays
*
*     A -- Double Precision Input matrix unchanged
*
*     B and W -- Double Precision Matrices, B must of at
*       least N+1 columns.   On return B is an approximation
*       of the input matrix A (to backward error) and W
*       holds the output from BHESS.
*
*     TOL -- Double precision parameter to control multiplier
*            size in BHESS.  Unchanged.
*
*     V  -- Double precision work vector of three columns.


*
*   Local Variables
*-------------------------------

*
*   Copy A to B
*
C     .. Scalar Arguments ..
      DOUBLE PRECISION TOL
      INTEGER LDX,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDX,*),B(LDX,*),V(LDX,3),W(LDX,*)
      INTEGER KROW(*),KRVECT(*),NU(*),PIV(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONDIT,FNORM,RNORM,SUM,SUM2
      REAL SUMA
      INTEGER I,IFLAG,INFO,J,MXBAND
C     ..
C     .. External Subroutines ..
      EXTERNAL BHAP1,BHAP3,BHAPC,BHESS,BHRPCK,DCOPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT,MAX
C     ..
      DO J = 1,N
          DO I = 1,N
              B(I,J) = A(I,J)
          END DO
      END DO
*--------------------------------------------------------
*
*   Call BHESS with B  Warning B requires N + 1 Columns
*
*--------------------------------------------------------

      CALL BHESS(LDX,B,N,PIV,NU,KRVECT,KROW,V,TOL,INFO)
*
*   Copy B to W
*
      DO J = 1,N
          DO I = 1,N
              W(I,J) = B(I,J)
          END DO
      END DO

*
*   Zero out all but upper Hessenberg in B
*
      CALL BHRPCK(LDX,N,B,NU)
*
*   Calculate backward A as B = Z * B * inv(Z)

*
*      First B <-- Z * B
      DO I = 1,N - 1,2
          CALL BHAP3(N,W,LDX,B(1,I),B(1,I+1),PIV,KROW,IFLAG)
      END DO
      IF (N.NE.2* (N/2)) THEN
          CALL BHAP3(N,W,LDX,B(1,N),V(1,2),PIV,KROW,IFLAG)
      END IF
*
*   Then B <-- B * inv(Z)       So B approximates A
*
      DO I = 1,N - 1,2
          CALL DCOPY(N,B(I,1),LDX,V(1,1),1)
          CALL DCOPY(N,B(I+1,1),LDX,V(1,2),1)
          CALL BHAP1(N,W,LDX,V(1,1),V(1,2),PIV,KROW,IFLAG)
          CALL DCOPY(N,V(1,1),1,B(I,1),LDX)
          CALL DCOPY(N,V(1,2),1,B(I+1,1),LDX)
      END DO
      IF (N.NE.2* (N/2)) THEN
          CALL DCOPY(N,B(N,1),LDX,V(1,1),1)
          CALL BHAP1(N,W,LDX,V(1,1),V(1,2),PIV,KROW,IFLAG)
          CALL DCOPY(N,V(1,1),1,B(N,1),LDX)
      END IF
*
*   Calculate || A - B ||
*
      SUM = 0.D0
      SUMA = 0.D0
      DO I = 1,N
          DO J = 1,N
              SUM = SUM + (A(I,J)-B(I,J))**2
              SUM2 = SUM2 + A(I,J)**2
          END DO
      END DO
* Frobenius norm normalized so that norm(eye(n)) = 1
      FNORM = DSQRT(SUM/N)
* Relative Frobenius norm
      RNORM = DSQRT(SUM/SUM2)
*
*   Estimate conditioning of Z
*
      CALL BHAPC(N,W,LDX,V(1,1),V(1,2),PIV,KROW,CONDIT,IFLAG)
*
*      Calculate Maximal Bandwidth
*
      MXBAND = 1
      DO I = 1,N
          MXBAND = MAX(NU(I),MXBAND)
      END DO
      WRITE(*,*)'    '
      WRITE(*,*)'    '
      WRITE(*,*)'                             n =  ',N
      WRITE(*,*)'                           tol = ',TOL
      WRITE(*,*)' Frobenius Backward Error Norm = ',FNORM
      WRITE(*,*)'  Relative Backward Error Norm = ',RNORM
      WRITE(*,*)'    Estimated Condition Number = ',CONDIT
      WRITE(*,*)'                     Bandwidth =  ',MXBAND

      RETURN

      END
