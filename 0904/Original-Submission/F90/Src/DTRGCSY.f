CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRGCSY( JOB, TRANZ, TRANAD, TRANBE, ISGN, I, J, M, MB, 
     $                    N, NB, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE,
     $                    F, LDF, DWORK, LDWORK, IWORK, LIWORK, ROWS, 
     $                    COLS, XYRST, XYCST, SCALE, INFO )
C
C  -- LAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     June 28, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          JOB, TRANAD, TRANBE, TRANZ
      INTEGER            I, J, MB, NB, INFO, ISGN, LDA, LDB, LDC, M, N,
     $                   ROWS, COLS, XYRST, XYCST, LDD, LDE, LDF, 
     $                   LIWORK, LDWORK
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), E( * ),
     $                   F( * ), DWORK( * )
C     ..
C
C     Purpose and description
C     =======================
C
C     Blocked GCSY solver, which solves for the (I,J)th
C     MB-by-NB subsystem or performs GEMM-updates with respect to
C     the already solved (I,J)th subsystem. What to do is
C     controlled by JOB (which takes the values 'S' for solve, 'U'
C     for update). Notice that if JOB = 'U', a previous call with 
C     JOB = 'S' for the same values of I and J is assumed such that 
C     the variables ROWS, COLS, XRST and XCST is input from that 
C     previuos call. Moreover, if this routine returned with
C     SCALE < 1.0 for JOB = 'S', the user must take care of the
C     scaling of the right hand side (excluding the (I,J)th 
C     subsolution) before re-calling this routine with JOB = 'U'. 
C
      LOGICAL SOLVE, UPDATE, TRANSAD, TRANSBE, ADEXT, BEEXT, ADDIM, 
     $        BEDIM, TRANSZ
      INTEGER IPA, IPB, IPC, IPD, IPE, IPF, K, IPX, IPY, ADRWS, THREADS,
     $        ADCLS, BERWS, BECLS, IPAKK, IPBKK, UPLOSIGN, II
      DOUBLE PRECISION USIGN, DIF, DUMMY, MACHINE(10)
C     
      LOGICAL  LSAME
      INTEGER  ICEIL, OMP_GET_NUM_THREADS
      EXTERNAL LSAME, ICEIL, OMP_GET_NUM_THREADS
C     
      SOLVE  = LSAME(JOB,'S') .OR. LSAME(JOB,'B')
      UPDATE = LSAME(JOB,'U') .OR. LSAME(JOB,'B')
      TRANSAD = LSAME(TRANAD,'T')
      TRANSBE = LSAME(TRANBE,'T')
      TRANSZ = LSAME(TRANZ,'T')
      USIGN = DBLE(-ISGN)
C
      IF( I.GT.ICEIL(M,MB) .OR. J.GT.ICEIL(N,NB) ) THEN
         ROWS = 0
         COLS = 0
         XYRST = 1 
         XYCST = 1
         RETURN
      END IF
C
      IF( SOLVE .OR. UPDATE ) THEN
         IF( MB.GT.M ) THEN
            ROWS = M
         ELSE
            ROWS = MIN( MB, M - (I-1)*MB )
         END IF
         IF( NB.GT.N ) THEN
            COLS = N
         ELSE
            COLS = MIN( NB, N - (J-1)*NB )
         END IF
         IPA = (I-1)*MB*LDA + (I-1)*MB + 1
         IPB = (J-1)*NB*LDB + (J-1)*NB + 1
         IPC = (J-1)*NB*LDC + (I-1)*MB + 1
         IPD = (I-1)*MB*LDD + (I-1)*MB + 1
         IPE = (J-1)*NB*LDE + (J-1)*NB + 1
         IPF = (J-1)*NB*LDF + (I-1)*MB + 1
         IF( I.LT.ICEIL(M,MB) ) THEN
            IF( A(IPA+(MB-1)*LDA+MB) .NE. 0.0D+00 ) THEN
               ROWS = ROWS + 1
               ADEXT = .TRUE.
            ELSE
               ADEXT = .FALSE.
            END IF
         ELSE
            ADEXT = .FALSE. 
         END IF
         IF( I.GT.1 ) THEN
            IF( A(IPA-LDA) .NE. 0.0D+00 ) THEN
               IPA = IPA + LDA + 1
               IPD = IPD + LDD + 1
               IPC = IPC + 1
               IPF = IPF + 1
               ROWS = ROWS - 1
               ADDIM = .TRUE.
            ELSE
               ADDIM = .FALSE.
            END IF
         ELSE
            ADDIM = .FALSE. 
         END IF
         IF( J.LT.ICEIL(N,NB) ) THEN
            IF( B(IPB+(NB-1)*LDB+NB) .NE. 0.0D+00 ) THEN
               COLS = COLS + 1
               BEEXT = .TRUE.
            ELSE
               BEEXT = .FALSE.
            END IF
         ELSE
            BEEXT = .FALSE.  
         END IF
         IF( J.GT.1 ) THEN
            IF( B(IPB-LDB) .NE. 0.0D+00 ) THEN
               IPB = IPB + LDB + 1
               IPE = IPE + LDE + 1
               IPC = IPC + LDC
               IPF = IPF + LDF
               COLS = COLS - 1
               BEDIM = .TRUE.
            ELSE
               BEDIM = .FALSE.
            END IF
         ELSE
            BEDIM = .FALSE.  
         END IF
      END IF
C
      IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) RETURN
C
      IF( SOLVE ) THEN
         IF( .NOT. TRANSZ ) THEN
            UPLOSIGN = (1-ISGN) / 2
            IF( TRANSAD ) UPLOSIGN = UPLOSIGN + 4
            IF( TRANSBE ) UPLOSIGN = UPLOSIGN + 2
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
            THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
            CALL RECSY_MACHINE( MACHINE )
            CALL RECGCSY_P( THREADS, UPLOSIGN, SCALE, ROWS, COLS, 
     $           A(IPA), LDA, B(IPB), LDB, C(IPC), LDC, D(IPD), LDD, 
     $           E(IPE), LDE, F(IPF), LDF, INFO, MACHINE, DWORK, 
     $           LDWORK )
#else
            CALL RECSY_MACHINE( MACHINE )
            CALL RECGCSY( UPLOSIGN, SCALE, ROWS, COLS, A(IPA), LDA, 
     $           B(IPB), LDB, C(IPC), LDC, D(IPD), LDD, E(IPE), LDE, 
     $           F(IPF), LDF, INFO, MACHINE, DWORK, LDWORK )
#endif 
         ELSE
            IF( .NOT. TRANSAD .AND. .NOT. TRANSBE ) THEN
               IF( ISGN.EQ.1 ) THEN
                  DO 2 II = 1, COLS
                     CALL DSCAL( ROWS, -1.0D+00, F(IPF+(II-1)*LDF), 1 ) 
 2                CONTINUE
               END IF
               CALL DTGSYL( 'Transpose', 0, ROWS, COLS,
     $              A(IPA), LDA, B(IPB), LDB, C(IPC), LDC, 
     $              D(IPD), LDD, E(IPE), LDE, F(IPF), LDF, SCALE, DIF, 
     $              DUMMY, 1, IWORK, INFO )
            ELSEIF( TRANSAD .AND. TRANSBE ) THEN
               DO 6 II = 1, COLS
                  CALL DSCAL( ROWS, -1.0D+00, C(IPC+(II-1)*LDC), 1 ) 
 6             CONTINUE
               CALL DLATCPY( 'All', ROWS, COLS, C(IPC), LDC, DWORK, 
     $              COLS )
               IF( ISGN.EQ.-1 ) THEN
                  DO 8 II = 1, COLS
                     CALL DSCAL( ROWS, -1.0D+00, F(IPF+(II-1)*LDF), 1 ) 
 8                CONTINUE
               END IF
               CALL DLATCPY( 'All', ROWS, COLS, F(IPF), LDF, 
     $              DWORK(ROWS*COLS+1), COLS )
               CALL DTGSYL( 'Transpose', 0, COLS, ROWS, B(IPB), LDB,
     $              A(IPA), LDA, DWORK(ROWS*COLS+1), COLS, E(IPE), LDE, 
     $              D(IPD), LDD, DWORK, COLS, SCALE, DIF, DUMMY, 1, 
     $              IWORK, INFO )
               CALL DLATCPY( 'All', COLS, ROWS, DWORK, COLS, F(IPF), 
     $                       LDF )
               CALL DLATCPY( 'All', COLS, ROWS, DWORK(ROWS*COLS+1), 
     $                       COLS, C(IPC), LDC )
            END IF
         END IF
         XYRST = MOD(IPC-1,LDC) + 1
         XYCST = ICEIL(IPC, LDC)
      END IF
C     
      IF( UPDATE ) THEN
         IF( .NOT. TRANSZ ) THEN
            IPX = IPC
            IPY = IPF
            IF( MB.LT.M ) THEN
               IF( .NOT. TRANSAD ) THEN
                  DO 10 K = 1, I-1
                     ADRWS = MIN( MB, M - (K-1)*MB )
                     IPA = (I-1)*MB*LDA + (K-1)*MB + 1
                     IPD = (I-1)*MB*LDD + (K-1)*MB + 1
                     IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                     IPF = (J-1)*NB*LDF + (K-1)*MB + 1
                     IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                     IF( K.LT.ICEIL(M,MB) ) THEN
                        IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                           ADRWS = ADRWS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                           ADRWS = ADRWS - 1
                           IPA = IPA + 1
                           IPD = IPD + 1
                           IPC = IPC + 1
                           IPF = IPF + 1
                        END IF
                     END IF
                     IF( ADDIM ) THEN
                        IPA = IPA + LDA
                        IPD = IPD + LDD
                     END IF
                     IF( BEDIM ) THEN
                        IPC = IPC + LDC
                        IPF = IPF + LDF
                     END IF
                     IF( ADRWS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'N', ADRWS, COLS, ROWS, 
     $                       -1.0D+00, A(IPA), LDA, C(IPX), LDC, 
     $                       1.0D+00, C(IPC), LDC )
                        CALL DGEMM( 'N', 'N', ADRWS, COLS, ROWS, 
     $                       -1.0D+00, D(IPD), LDD, C(IPX), LDC, 
     $                       1.0D+00, F(IPF), LDF )
                     END IF
 10               CONTINUE
               ELSEIF( TRANSAD ) THEN
                  DO 20 K = I+1, ICEIL(M,MB)
                     ADCLS = MIN( MB, M - (K-1)*MB )
                     IPA = (K-1)*MB*LDA + (I-1)*MB + 1
                     IPD = (K-1)*MB*LDD + (I-1)*MB + 1
                     IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                     IPF = (J-1)*NB*LDF + (K-1)*MB + 1
                     IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                     IF( K.LT.ICEIL(M,MB) ) THEN
                        IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                           ADCLS = ADCLS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                           ADCLS = ADCLS - 1
                           IPA = IPA + LDA
                           IPD = IPD + LDD
                           IPC = IPC + 1
                           IPF = IPF + 1
                        END IF
                     END IF
                     IF( ADDIM ) THEN
                        IPA = IPA + 1
                        IPD = IPD + 1
                     END IF
                     IF( BEDIM ) THEN
                        IPC = IPC + LDC
                        IPF = IPF + LDF
                     END IF
                     IF( ADCLS.GT.0 ) THEN
                        CALL DGEMM( 'T', 'N', ADCLS, COLS, ROWS, 
     $                       -1.0D+00, A(IPA), LDA, C(IPX), LDC, 
     $                       1.0D+00, C(IPC), LDC )
                        CALL DGEMM( 'T', 'N', ADCLS, COLS, ROWS, 
     $                       -1.0D+00, D(IPD), LDD, C(IPX), LDC, 
     $                       1.0D+00, F(IPF), LDF )
                     END IF
 20               CONTINUE
               END IF
            END IF
C     
            IF( NB.LT.N ) THEN
               IF( .NOT. TRANSBE ) THEN
                  DO 30 K = J+1, ICEIL(N,NB)
                     BECLS = MIN( NB, N - (K-1)*NB )
                     IPB = (K-1)*NB*LDB + (J-1)*NB + 1
                     IPE = (K-1)*NB*LDE + (J-1)*NB + 1
                     IPC = (K-1)*NB*LDC + (I-1)*MB + 1
                     IPF = (K-1)*NB*LDF + (I-1)*MB + 1
                     IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                     IF( K.LT.ICEIL(N,NB) ) THEN
                        IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                           BECLS = BECLS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                           BECLS = BECLS - 1
                           IPB = IPB + LDB
                           IPE = IPE + LDE
                           IPC = IPC + LDC
                           IPF = IPF + LDF
                        END IF
                     END IF
                     IF( BEDIM ) THEN
                        IPB = IPB + 1
                        IPE = IPE + 1
                     END IF
                     IF( ADDIM ) THEN
                        IPC = IPC + 1
                        IPF = IPF + 1
                     END IF
                     IF( BECLS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'N', ROWS, BECLS, COLS, USIGN, 
     $                       F(IPY), LDF, B(IPB), LDB, 1.0D+00, C(IPC), 
     $                       LDC )
                        CALL DGEMM( 'N', 'N', ROWS, BECLS, COLS, USIGN, 
     $                       F(IPY), LDF, E(IPE), LDE, 1.0D+00, F(IPF), 
     $                       LDF )
                     END IF
 30               CONTINUE
               ELSEIF( TRANSBE ) THEN
                  DO 40 K = 1, J-1
                     BERWS = MIN( NB, N - (K-1)*NB )
                     IPB = (J-1)*NB*LDB + (K-1)*NB + 1
                     IPE = (J-1)*NB*LDE + (K-1)*NB + 1
                     IPC = (K-1)*NB*LDC + (I-1)*MB + 1
                     IPF = (K-1)*NB*LDF + (I-1)*MB + 1
                     IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                     IF( K.LT.ICEIL(N,NB) ) THEN
                        IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                           BERWS = BERWS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                           BERWS = BERWS - 1
                           IPB = IPB + 1
                           IPE = IPE + 1
                           IPC = IPC + LDC
                           IPF = IPF + LDF
                        END IF
                     END IF
                     IF( BEDIM ) THEN
                        IPB = IPB + LDB
                        IPE = IPE + LDE
                     END IF
                     IF( ADDIM ) THEN
                        IPC = IPC + 1
                        IPF = IPF + 1
                     END IF
                     IF( BERWS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'T', ROWS, BERWS, COLS, USIGN, 
     $                       F(IPY), LDF, B(IPB), LDB, 1.0D+00, C(IPC), 
     $                       LDC )
                        CALL DGEMM( 'N', 'T', ROWS, BERWS, COLS, USIGN, 
     $                       F(IPY), LDF, E(IPE), LDE, 1.0D+00, F(IPF), 
     $                       LDF )
                     END IF
 40               CONTINUE
               END IF
            END IF
         ELSE
            IPX = IPC
            IPY = IPF
            IF( MB.LT.M ) THEN
               IF( TRANSAD ) THEN
                  DO 50 K = 1, I-1
                     ADRWS = MIN( MB, M - (K-1)*MB )
                     IPA = (I-1)*MB*LDA + (K-1)*MB + 1
                     IPD = (I-1)*MB*LDD + (K-1)*MB + 1
                     IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                     IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                     IF( K.LT.ICEIL(M,MB) ) THEN
                        IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                           ADRWS = ADRWS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                           ADRWS = ADRWS - 1
                           IPA = IPA + 1
                           IPD = IPD + 1
                           IPC = IPC + 1
                        END IF
                     END IF
                     IF( ADDIM ) THEN
                        IPA = IPA + LDA
                        IPD = IPD + LDD
                     END IF
                     IF( BEDIM ) THEN
                        IPC = IPC + LDC
                     END IF
                     IF( ADRWS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'N', ADRWS, COLS, ROWS, 
     $                       -1.0D+00, A(IPA), LDA, C(IPX), LDC, 
     $                       1.0D+00, C(IPC), LDC )
                        CALL DGEMM( 'N', 'N', ADRWS, COLS, ROWS, 
     $                       -1.0D+00, D(IPD), LDD, F(IPY), LDF, 
     $                       1.0D+00, C(IPC), LDC )
                     END IF
 50               CONTINUE
               ELSEIF( .NOT. TRANSAD ) THEN
                  DO 60 K = I+1, ICEIL(M,MB)
                     ADCLS = MIN( MB, M - (K-1)*MB )
                     IPA = (K-1)*MB*LDA + (I-1)*MB + 1
                     IPD = (K-1)*MB*LDD + (I-1)*MB + 1
                     IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                     IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                     IF( K.LT.ICEIL(M,MB) ) THEN
                        IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                           ADCLS = ADCLS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                           ADCLS = ADCLS - 1
                           IPA = IPA + LDA
                           IPD = IPD + LDD
                           IPC = IPC + 1
                        END IF
                     END IF
                     IF( ADDIM ) THEN
                        IPA = IPA + 1
                        IPD = IPD + 1
                     END IF
                     IF( BEDIM ) THEN
                        IPC = IPC + LDC
                     END IF
                     IF( ADCLS.GT.0 ) THEN
                        CALL DGEMM( 'T', 'N', ADCLS, COLS, ROWS, 
     $                       -1.0D+00, A(IPA), LDA, C(IPX), LDC, 
     $                       1.0D+00, C(IPC), LDC )
                        CALL DGEMM( 'T', 'N', ADCLS, COLS, ROWS, 
     $                       -1.0D+00, D(IPD), LDD, F(IPY), LDF, 
     $                       1.0D+00, C(IPC), LDC )
                     END IF
 60               CONTINUE
               END IF
            END IF
C     
            IF( NB.LT.N ) THEN
               IF( TRANSBE ) THEN
                  DO 70 K = J+1, ICEIL(N,NB)
                     BECLS = MIN( NB, N - (K-1)*NB )
                     IPB = (K-1)*NB*LDB + (J-1)*NB + 1
                     IPE = (K-1)*NB*LDE + (J-1)*NB + 1
                     IPF = (K-1)*NB*LDF + (I-1)*MB + 1
                     IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                     IF( K.LT.ICEIL(N,NB) ) THEN
                        IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                           BECLS = BECLS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                           BECLS = BECLS - 1
                           IPB = IPB + LDB
                           IPE = IPE + LDE
                           IPF = IPF + LDF
                        END IF
                     END IF
                     IF( BEDIM ) THEN
                        IPB = IPB + 1
                        IPE = IPE + 1
                     END IF
                     IF( ADDIM ) THEN
                        IPF = IPF + 1
                     END IF
                     IF( BECLS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'N', ROWS, BECLS, COLS, USIGN, 
     $                       C(IPX), LDC, B(IPB), LDB, 1.0D+00, F(IPF), 
     $                       LDF )
                        CALL DGEMM( 'N', 'N', ROWS, BECLS, COLS, USIGN, 
     $                       F(IPY), LDF, E(IPE), LDE, 1.0D+00, F(IPF), 
     $                       LDF )
                     END IF
 70               CONTINUE
               ELSEIF( .NOT. TRANSBE ) THEN
                  DO 80 K = 1, J-1
                     BERWS = MIN( NB, N - (K-1)*NB )
                     IPB = (J-1)*NB*LDB + (K-1)*NB + 1
                     IPE = (J-1)*NB*LDE + (K-1)*NB + 1
                     IPF = (K-1)*NB*LDF + (I-1)*MB + 1
                     IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                     IF( K.LT.ICEIL(N,NB) ) THEN
                        IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                           BERWS = BERWS + 1
                        END IF
                     END IF
                     IF( K.GT.1 ) THEN
                        IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                           BERWS = BERWS - 1
                           IPB = IPB + 1
                           IPE = IPE + 1
                           IPF = IPF + LDF
                        END IF
                     END IF
                     IF( BEDIM ) THEN
                        IPB = IPB + LDB
                        IPE = IPE + LDE
                     END IF
                     IF( ADDIM ) THEN
                        IPF = IPF + 1
                     END IF
                     IF( BERWS.GT.0 ) THEN
                        CALL DGEMM( 'N', 'T', ROWS, BERWS, COLS, USIGN, 
     $                       C(IPX), LDC, B(IPB), LDB, 1.0D+00, F(IPF), 
     $                       LDF )
                        CALL DGEMM( 'N', 'T', ROWS, BERWS, COLS, USIGN, 
     $                       F(IPY), LDF, E(IPE), LDE, 1.0D+00, F(IPF), 
     $                       LDF )
                     END IF
 80               CONTINUE
               END IF
            END IF
         END IF
      END IF
C     
      END
C     
C     End of DTRGCSY
C     
C *** Last line of DTRGCSY ***
