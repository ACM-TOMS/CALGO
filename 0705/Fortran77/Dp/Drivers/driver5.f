C--**--CH2913--705--P:LAP--28:10:1999
C--**--CH2897--705--C:SU--28:10:1999
C--**--CH2888--705--B:MA--28:10:1999
C--**--CH2866--705--A:H--28:10:1999
C--**--CH2860--705--P:RW--28:10:1999
C     -- PROGRAM TSYLGCF --
C
C     THIS PROGRAM TESTS THE SOFTWARE FOR SOLVING THE SYMMETRIC
C     SYLVESTER EQUATION
C
C        A*X*E' + E*X*A' + Q = 0    (' DENOTES TRANSPOSE)
C
C     READING THE COEFFICIENT MATRICES FROM A FILE.
C
C     THE SOLUTION IS CHECKED BY COMPUTING THE RESIDUAL AND ITS RELATIVE
C     1-NORM.  THE CONDITION ESTIMATE IS CHECKED FOR SMALL SYSTEMS ONLY
C     (N*N <= 40) BY FORMING THE KRONECKER PRODUCT MATRIX AND COMPUTING
C     ITS SINGULAR VALUE DECOMPOSITION.
C
C     INPUT FILE FORMAT (SEE EXAMPLE FILE TESTCD.IN):
C       INPUT IS FREE FORMAT; MATRICES ARE READ BY ROWS.
C       DIMENSION N (INTEGER)
C       A; EACH ROW MUST BEGIN ON A NEW LINE AND MAY TAKE AS MANY LINES
C       AS NECESSARY.
C       E; Q; SAME COMMENT.
C
C     SOLUTION MATRIX X IS WRITTEN TO AN OUTPUT FILE (SEE EXAMPLE FILE
C       TESTC.OUT).
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       SYLGC; (LINPACK) DSVDC; (MATU) MADD MULC TRNATA D1NRM MSAVE
C
C     WRITTEN -
C       14DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     -- MAXIMUM N IS 8 --
      INTEGER NAE, NQ, NW, N, IERR, NWX, IDMY
      INTEGER NMAX, WRKLEN, MAXSVD
      PARAMETER (NMAX=8, WRKLEN=2*NMAX*NMAX+3*NMAX,
     +           MAXSVD=40)
      DOUBLE PRECISION A(NMAX,NMAX), E(NMAX,NMAX), Q(NMAX,NMAX)
      DOUBLE PRECISION WKM1(NMAX,NMAX), WKM2(NMAX,NMAX),
     +                 WKV(WRKLEN)
C
      INTEGER I, J, NXN, II, JJ, INDX, JNDX
      DOUBLE PRECISION A1(NMAX,NMAX), E1(NMAX,NMAX), Q1(NMAX,NMAX), 
     +                 WKMX(MAXSVD,MAXSVD)

      DOUBLE PRECISION NRMR, NRMQ, RCOND, RCDSVD
      DOUBLE PRECISION RATIO, DUMY(1)
      DOUBLE PRECISION D1NRM
C     CHARACTER*20 INNAME, OUTNAM
      CHARACTER YESNO
C
      NAE = NMAX
      NQ = NMAX
      NW = NMAX
      NWX = MAXSVD
C
C     -- INITIALIZE --
      WRITE(*,90000)
      WRITE(*,90001)
90000 FORMAT("             T        T")
C
C
C      WRITE(*,91002)
C91002 FORMAT(/' ENTER NAME OF FILE CONTAINING COEFFICIENT MATRICES ')
C      READ(5,91003) INNAME
C91003 FORMAT(A)
C      OPEN(2,FILE=INNAME)
C      WRITE(*,91004)
C91004 FORMAT(/' ENTER NAME OF FILE FOR OUTPUT OF SOLUTION MATRIX')
C      READ(5,91003) OUTNAM
C      OPEN(3,FILE=OUTNAM)
C
C     -- READ COEFFICIENT MATRICES --
90001 FORMAT(" SOLVE  A*X*E  + E*X*A  + Q = 0  USING SYLGC")
      READ(*,*) N
      IF (N .GT. NAE) THEN
         WRITE(*,91005) N, NAE
91005    FORMAT(/' DIMENSION OF PROBLEM IS TOO LARGE:  N =',I2,         
     + '  MAXIMUM IS ',I2)
         STOP
      ENDIF
      WRITE(*,91008) N
C
91008 FORMAT(/' N =',I2)
      DO 201 I=1,N
         READ(*,*)(A(I,J),J=1,N)
 201  CONTINUE
      DO 202 I=1,N
         READ(*,*)(E(I,J),J=1,N)
 202  CONTINUE
      DO 203 I=1,N
         READ(*,*)(Q(I,J),J=1,N)
 203  CONTINUE
C
C        -- SAVE A COPY --
         CALL MSAVE(NAE, NAE, N, N, A, A1)
         CALL MSAVE(NAE, NAE, N, N, E, E1)
         CALL MSAVE(NQ, NQ, N, N, Q, Q1)
C
C        -- COMPUTE NORM OF Q --
         NRMQ = D1NRM(NQ, N, N, Q)
C
C        -- COMPUTE SOLUTION AND ESTIMATE CONDITION --
         IERR = 1
         CALL SYLGC(NAE, NQ, N, A, E, Q, WKV, IERR, RCOND)
         IF (IERR .NE. 0) THEN
            WRITE(*,90003) IERR
90003       FORMAT(" ERROR FROM SYLGC, IERR=", I2)
         ENDIF
C
C        -- SAVE SOLUTION --
         DO 301 I=1,N
            WRITE(*,92001) (Q(I,J),J=1,N)
92001       FORMAT(1X,6E12.4)
 301     CONTINUE
C
C        -- COMPUTE RESIDUAL --
         CALL MULC(NAE, NQ, NW, N, N, N, A1, Q, WKM1)
         CALL TRNATA(NAE, N, E1)
         CALL MULC(NW, NAE, NW, N, N, N, WKM1, E1, WKM2)
         CALL TRNATA(NAE, N, E1)
         CALL MADD(NQ, NW, NQ, N, N, Q1, WKM2, Q1)
         CALL MULC(NAE, NQ, NW, N, N, N, E1, Q, WKM1)
         CALL TRNATA(NAE, N, A1)
         CALL MULC(NW, NAE, NW, N, N, N, WKM1, A1, WKM2)
         CALL TRNATA(NAE, N, A1)
         CALL MADD(NQ, NW, NQ, N, N, Q1, WKM2, Q1)
C
C        -- COMPUTE NORM OF RESIDUAL --
         NRMR = D1NRM(NQ, N, N, Q1) / NRMQ
C
C        -- PRINT RESULTS --
         WRITE(*,90002) NRMR, RCOND
C
C     -- COMPUTE REAL CONDITION NUMBER IF REQUESTED --
90002    FORMAT(1X, '  (NRM RESID)/(NRM Q)=', E10.3,          '  RCOND(E
     +ST)=', E10.3)
      NXN = N*N
      IF (NXN .LE. MAXSVD) THEN
         WRITE(*,91006)
91006    FORMAT(/' DO YOU WISH TO HAVE THE ACTUAL CONDITION NUMBER',    
     +      ' COMPUTED? (Y OR N) ')
         READ(5,*)YESNO
         IF (YESNO .EQ. 'y'  .OR.  YESNO .EQ. 'Y') THEN
            DO 140 J = 1,N
               DO 130 JJ = 1,N
                  DO 120 I = 1,N
                     DO 110 II = 1,N
                        INDX = N*(I-1) + II
                        JNDX = N*(J-1) + JJ
                        WKMX(INDX,JNDX) = A1(II,JJ)*E1(I,J)
     *                                    + E1(II,JJ)*A1(I,J)
  110                CONTINUE
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
            IDMY = 1
            CALL DSVDC(WKMX, NWX, NXN, NXN, WKV, WKV(NXN+1), DUMY, IDMY,
     *                 DUMY, IDMY, WKM1, 00, IERR)
            IF (IERR .NE. 0) THEN
               WRITE(*,90006)
90006          FORMAT(" SVD FAILED")
               RCDSVD = 0.0D0
            ELSE
               RCDSVD = WKV(NXN) / WKV(1)
            ENDIF
            RATIO = RCOND / RCDSVD
C
            WRITE(*,90007) RCDSVD, RATIO
90007       FORMAT(1X, '  RCOND(TRUE)=', E10.3, ' EST/TRUE=', E10.3)
         ENDIF
      ENDIF
C
      CLOSE(2)
      CLOSE(3)
      STOP
C --- LAST LINE OF TSYLGCF ---
      END