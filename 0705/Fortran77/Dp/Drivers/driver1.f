C--**--CH2907--705--P:LAP--28:10:1999
C--**--CH2876--705--B:MA--28:10:1999
C--**--CH2862--705--A:H--28:10:1999
C     -- PROGRAM TSYLG --
C
C     THIS PROGRAM TESTS THE SOFTWARE FOR SOLVING THE GENERAL
C     SYLVESTER EQUATION
C
C        A * X * B' +  C * X * D' =  E    (' DENOTES TRANSPOSE)
C
C     USING A SET OF ILL-CONDITIONED EQUATIONS.  SEE "SOLUTION OF THE
C     SYLVESTER MATRIX EQUATION AXB'+CXD'=E" BY GARDINER, LAUB, AMATO,
C     AND MOLER FOR A DESCRIPTION OF THE TEST PROBLEMS.
C
C     THE SOLUTION IS CHECKED BY COMPUTING THE RESIDUAL AND ITS RELATIVE
C     1-NORM.  THE CONDITION ESTIMATE IS CHECKED FOR SMALL SYSTEMS ONLY
C     (M*N <= 40) BY FORMING THE KRONECKER PRODUCT MATRIX AND COMPUTING
C     ITS SINGULAR VALUE DECOMPOSITION.
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       SYLG; (LINPACK) DSVDC;
C       (MATU) MSAVE, D1NRM, MSUB, MULC, TRNATA
C
C     WRITTEN -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C     REVISED -
C       17DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     -- MAXIMUM M IS 9, MAXIMUM N IS 8 --
      INTEGER NAC, NBD, NE, NW1, NW2, NW3, M, N, IERR
      INTEGER MMAX, NMAX, MNMAX, WRKLEN, MAXSVD
      PARAMETER (MMAX=9, NMAX=8, MNMAX=9, MAXSVD=40,
     +  WRKLEN = 2*MMAX*MMAX+NMAX*NMAX+NMAX*MMAX+7*MMAX+MNMAX*MNMAX)
      DOUBLE PRECISION A0(MMAX,MMAX), B0(NMAX,NMAX), C0(MMAX,MMAX),
     +                 D0(NMAX,NMAX), E0(MMAX,NMAX)
      DOUBLE PRECISION A1(MMAX,MMAX), B1(NMAX,NMAX), C1(MMAX,MMAX),
     +                 D1(NMAX,NMAX), E1(MMAX,NMAX)
      DOUBLE PRECISION WKM1(MMAX,NMAX), WKM2(MMAX,NMAX), 
     +                 WKM3(MAXSVD,MAXSVD), WKV(WRKLEN)
      INTEGER IWKV(2*MMAX)
C
      INTEGER I, J, II, JJ, INDX, JNDX, MXN
      INTEGER IDMY, KT
      DOUBLE PRECISION TMP, NRMR, NRME, WORST, RCOND, RCDSVD
      DOUBLE PRECISION DUMY(1), WORSTL, WORSTH, RATIO, S
      DOUBLE PRECISION D1NRM
C     --- OUTPUT UNIT FOR WRITE ---
      NAC = MMAX
      NBD = NMAX
      NE = MMAX
      NW1 = MMAX
      NW2 = MMAX
      NW3 = MAXSVD
C
C     -- INITIALIZE --
      WORST = 0.0D0
      WORSTL = 1.0D100
      WORSTH = 0.0D0
C
      WRITE(*,90000)
      WRITE(*,90001)
90000 FORMAT("             T        T")
90001 FORMAT(" SOLVE  A*X*B  + C*X*D  = E  USING SYLG")
C
C     USE THESE TWO LINES TO RUN A COMPLETE TEST
C
       DO 320 M = 2,NAC
       DO 310 N = 2,NBD
C
C     USE THESE TWO LINES TO RUN A MINIMAL INSTALLATION TEST
C
C     DO 320 M = 3,3
C     DO 310 N = 2,2
C
      WRITE(*,90005) M, N
C
90005 FORMAT(/ " M =", I3, "   N =", I3)
      DO 230 KT = 0,40,10
         S = (2.0D0)**(-KT)
C
C GENERATE COEFFICIENT MATRICES  --  TRUE SOLUTION IS MATRIX OF ALL ONES
C
         DO 120 I=1,M
            DO 110 J=1,I
               A0(I,J) = 1.0D0
               A0(J,I) = 0.0D0
               C0(J,I) = S
               C0(I,J) = 0.0D0
  110       CONTINUE
            A0(I,I) = I
            C0(I,I) = 1.0D0
  120    CONTINUE
C
         DO 140 I=1,N
            DO 130 J=1,I
               B0(J,I) = S
               B0(I,J) = 0.0D0
               D0(I,J) = 1.0D0
               D0(J,I) = 0.0D0
  130       CONTINUE
            B0(I,I) = 1.0D0
            D0(I,I) = S - (N-I+1)
  140    CONTINUE
C
         DO 160 I=1,M
            WKV(I) = 0.0D0
            WKV(M+I) = 0.0D0
            DO 150 J=1,M
               WKV(I) = WKV(I) + A0(I,J)
               WKV(M+I) = WKV(M+I) + C0(I,J)
  150       CONTINUE
  160    CONTINUE
         DO 190 I=1,N
            TMP = 0.0D0
            DUMY(1) = 0.0D0
            DO 170 J=1,N
               TMP = TMP + B0(I,J)
               DUMY(1) = DUMY(1) + D0(I,J)
  170       CONTINUE
            DO 180 J=1,M
               E0(J,I) = TMP * WKV(J) + DUMY(1) * WKV(M+J)
  180       CONTINUE
  190    CONTINUE
C
C        -- SAVE A COPY --
         CALL MSAVE(NAC, NAC, M, M, A0, A1)
         CALL MSAVE(NAC, NAC, M, M, C0, C1)
         CALL MSAVE(NBD, NBD, N, N, B0, B1)
         CALL MSAVE(NBD, NBD, N, N, D0, D1)
         CALL MSAVE(NE, NE, M, N, E0, E1)
C
C        -- COMPUTE NORM OF E --
         NRME = D1NRM(NE, M, N, E0)
C
C        -- COMPUTE SOLUTION, ESTIMATE CONDITION IF SVD COMPUTABLE --
         MXN = M*N
         IF (MXN .LE. MAXSVD) THEN
            IERR = 1
         ELSE
            IERR = 0
         ENDIF
         CALL SYLG(NAC, NBD, NE, M, N, A0, B0, C0, D0, E0, WKV, IWKV,
     *             IERR, RCOND)
         IF (IERR .NE. 0) THEN
            WRITE(*,90003) KT, IERR
90003       FORMAT(" AT ITERATION", I2, " ERROR IN SYLG, IERR=", I2)
         ENDIF
C
C        -- COMPUTE RESIDUAL --
         CALL MULC(NAC, NE, NW1, M, M, N, A1, E0, WKM1)
         CALL TRNATA(NBD, N, B1)
         CALL MULC(NW1, NBD, NW2, M, N, N, WKM1, B1, WKM2)
         CALL TRNATA(NBD, N, B1)
         CALL MSUB(NE, NW2, NE, M, N, E1, WKM2, E1)
         CALL MULC(NAC, NE, NW1, M, M, N, C1, E0, WKM1)
         CALL TRNATA(NBD, N, D1)
         CALL MULC(NW1, NBD, NW2, M, N, N, WKM1, D1, WKM2)
         CALL TRNATA(NBD, N, D1)
         CALL MSUB(NE, NW2, NE, M, N, E1, WKM2, E1)
C
C        -- COMPUTE NORM OF RESIDUAL --
         NRMR = D1NRM(NE, M, N, E1)
         TMP = NRMR/NRME
         IF (TMP .GT. WORST) WORST = TMP
C
C        -- COMPUTE REAL CONDITION NUMBER --
         IF (MXN .LE. MAXSVD) THEN
            DO 40 J = 1,N
               DO 30 JJ = 1,M
                  DO 20 I = 1,N
                     DO 10 II = 1,M
                        INDX = M*(I-1) + II
                        JNDX = M*(J-1) + JJ
                        WKM3(INDX,JNDX) = A1(II,JJ)*B1(I,J)
     &                                    + C1(II,JJ)*D1(I,J)
   10                CONTINUE
   20             CONTINUE
   30          CONTINUE
   40       CONTINUE
            IDMY = 1
            CALL DSVDC(WKM3, NW3, MXN, MXN, WKV, WKV(MXN+1), DUMY, IDMY,
     *                 DUMY, IDMY, WKM1, 00, IERR)
            IF (IERR .NE. 0) THEN
               WRITE(*,90006) KT
90006          FORMAT(" ITERATION", I2, " SVD FAILED")
               RCDSVD = 0.0D0
            ELSE
               RCDSVD = WKV(MXN) / WKV(1)
            ENDIF
            RATIO = RCOND / RCDSVD
            IF (RATIO .LT. WORSTL) WORSTL = RATIO
            IF (RATIO .GT. WORSTH) WORSTH = RATIO
         ENDIF
C
C        -- PRINT RESULTS --
         IF (MXN .LE. MAXSVD) THEN
            WRITE(*,90002) KT, TMP, RCOND, RATIO
         ELSE
            WRITE(*,90007) KT, TMP
         ENDIF
90002    FORMAT(1X, I2, '  (NRM RESID)/(NRM E)=', E10.3,          '  RCO
     +ND(EST)=', E10.3, ' EST/TRUE=', E10.3)
C
90007    FORMAT(1X, I2, '  (NRM RESID)/(NRM E)=', E10.3)
  230 CONTINUE
C
  310 CONTINUE
  320 CONTINUE
C
      WRITE(*,90004) WORST
90004 FORMAT(/ " WORST CASE RESIDUAL IS", E10.3)
      WRITE(*,90008) WORSTL, WORSTH
C
90008 FORMAT(' WORST CASE RCOND RATIOS (EST/TRUE) ARE', E10.3,       ' (
     +LOW) AND', E10.3, ' (HIGH)')
      WRITE(*,90010)
90010 FORMAT(// ' THE RESIDUAL SHOULD BE ON THE ORDER OF 1E-14',        
     +  ' OR SMALLER ON MOST MACHINES.')
      WRITE(*,90011)
C
90011 FORMAT(' THE RCOND RATIOS SHOULD BE ON THE ORDER OF 1,',       ' I
     +.E., BETWEEN .1 AND 10.')
      STOP
C --- LAST LINE OF TSYLG ---
      END
