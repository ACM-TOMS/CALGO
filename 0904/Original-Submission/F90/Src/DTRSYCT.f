CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRSYCT( JOB, TRANA, TRANB, ISGN, I, J, M, MB, N, NB, 
     $                    A, LDA, B, LDB, C, LDC, ROWS, COLS, XRST, 
     $                    XCST, SCALE, INFO )
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
      CHARACTER          JOB, TRANA, TRANB
      INTEGER            I, J, MB, NB, INFO, ISGN, LDA, LDB, LDC, M, N,
     $                   ROWS, COLS, XRST, XCST
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), B( * ), C( * )
C     ..
C
C     Purpose and description
C     =======================
C
C     Blocked SYCT solver, which solves for the (I,J)th
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
      LOGICAL SOLVE, UPDATE, TRANSA, TRANSB, AEXT, BEXT, ADIM, BDIM
      INTEGER IPA, IPB, IPC, K, IPX, ARWS, ACLS, BRWS, BCLS, IPAKK,
     $        IPBKK, UPLOSIGN, THREADS
      DOUBLE PRECISION USIGN, MACHINE(10)
C     
      LOGICAL  LSAME
      INTEGER  ICEIL, OMP_GET_NUM_THREADS
      EXTERNAL LSAME, ICEIL, OMP_GET_NUM_THREADS
C     
      SOLVE  = LSAME(JOB,'S') 
      UPDATE = LSAME(JOB,'U')
      TRANSA = LSAME(TRANA,'T')
      TRANSB = LSAME(TRANB,'T')
      USIGN = DBLE(-ISGN)
C
      IF( I.GT.ICEIL(M,MB) .OR. J.GT.ICEIL(N,NB) ) THEN
         ROWS = 0
         COLS = 0
         XRST = 1 
         XCST = 1
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
         IF( I.LT.ICEIL(M,MB) ) THEN
            IF( A(IPA+(MB-1)*LDA+MB) .NE. 0.0D+00 ) THEN
               ROWS = ROWS + 1
               AEXT = .TRUE.
            ELSE
               AEXT = .FALSE.
            END IF
         ELSE
            AEXT = .FALSE. 
         END IF
         IF( I.GT.1 ) THEN
            IF( A(IPA-LDA) .NE. 0.0D+00 ) THEN
               IPA = IPA + LDA + 1
               IPC = IPC + 1
               ROWS = ROWS - 1
               ADIM = .TRUE.
            ELSE
               ADIM = .FALSE.
            END IF
         ELSE
            ADIM = .FALSE. 
         END IF
         IF( J.LT.ICEIL(N,NB) ) THEN
            IF( B(IPB+(NB-1)*LDB+NB) .NE. 0.0D+00 ) THEN
               COLS = COLS + 1
               BEXT = .TRUE.
            ELSE
               BEXT = .FALSE.
            END IF
         ELSE
            BEXT = .FALSE.  
         END IF
         IF( J.GT.1 ) THEN
            IF( B(IPB-LDB) .NE. 0.0D+00 ) THEN
               IPB = IPB + LDB + 1
               IPC = IPC + LDC
               COLS = COLS - 1
               BDIM = .TRUE.
            ELSE
               BDIM = .FALSE.
            END IF
         ELSE
            BDIM = .FALSE.  
         END IF
      END IF
C
      IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) RETURN
C
      IF( SOLVE ) THEN
         UPLOSIGN = (1-ISGN) / 2
         IF ( TRANSA ) UPLOSIGN = UPLOSIGN + 4
         IF ( TRANSB ) UPLOSIGN = UPLOSIGN + 2
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
         THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
         CALL RECSY_MACHINE( MACHINE )
         CALL RECSYCT_P( THREADS, UPLOSIGN, SCALE, ROWS, COLS, A(IPA), 
     $                   LDA, B(IPB), LDB, C(IPC), LDC, INFO, MACHINE )
#else
         CALL RECSY_MACHINE( MACHINE )
         CALL RECSYCT( UPLOSIGN, SCALE, ROWS, COLS, A(IPA), LDA, 
     $                 B(IPB), LDB, C(IPC), LDC, INFO, MACHINE )
#endif
         XRST = MOD(IPC-1,LDC) + 1
         XCST = ICEIL(IPC, LDC)
      END IF
C
      IF( UPDATE ) THEN
         IPX = IPC
         IF( MB.LT.M ) THEN
            IF( .NOT. TRANSA ) THEN
               DO 10 K = 1, I-1
                  ARWS = MIN( MB, M - (K-1)*MB )
                  IPA = (I-1)*MB*LDA + (K-1)*MB + 1
                  IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.ICEIL(M,MB) ) THEN
                     IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS - 1
                        IPA = IPA + 1
                        IPC = IPC + 1  
                     END IF
                  END IF
                  IF( ADIM ) IPA = IPA + LDA
                  IF( BDIM ) IPC = IPC + LDC
                  IF( ARWS.GT.0 )
     $                 CALL DGEMM( 'N', 'N', ARWS, COLS, ROWS, -1.0D+00, 
     $                 A(IPA), LDA, C(IPX), LDC, 1.0D+00, C(IPC), LDC )
 10            CONTINUE
            ELSEIF( TRANSA ) THEN
               DO 20 K = I+1, ICEIL(M,MB)
                  ACLS = MIN( MB, M - (K-1)*MB )
                  IPA = (K-1)*MB*LDA + (I-1)*MB + 1
                  IPC = (J-1)*NB*LDC + (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.ICEIL(M,MB) ) THEN
                     IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS - 1
                        IPA = IPA + LDA
                        IPC = IPC + 1
                     END IF
                  END IF
                  IF( ADIM ) IPA = IPA + 1
                  IF( BDIM ) IPC = IPC + LDC
                  IF( ACLS.GT.0 )
     $                 CALL DGEMM( 'T', 'N', ACLS, COLS, ROWS, -1.0D+00, 
     $                 A(IPA), LDA, C(IPX), LDC, 1.0D+00, C(IPC), LDC )
 20            CONTINUE
            END IF
         END IF
C     
         IF( NB.LT.N ) THEN
            IF( .NOT. TRANSB ) THEN
               DO 30 K = J+1, ICEIL(N,NB)
                  BCLS = MIN( NB, N - (K-1)*NB )
                  IPB = (K-1)*NB*LDB + (J-1)*NB + 1
                  IPC = (K-1)*NB*LDC + (I-1)*MB + 1
                  IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                        BCLS = BCLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                        BCLS = BCLS - 1
                        IPB = IPB + LDB
                        IPC = IPC + LDC
                     END IF
                  END IF
                  IF( BDIM ) IPB = IPB + 1
                  IF( ADIM ) IPC = IPC + 1
                  IF( BCLS.GT.0 )
     $                 CALL DGEMM( 'N', 'N', ROWS, BCLS, COLS, USIGN, 
     $                 C(IPX), LDC, B(IPB), LDB, 1.0D+00, C(IPC), LDC )
 30            CONTINUE
            ELSEIF( TRANSB ) THEN
               DO 40 K = 1, J-1
                  BRWS = MIN( NB, N - (K-1)*NB )
                  IPB = (J-1)*NB*LDB + (K-1)*NB + 1
                  IPC = (K-1)*NB*LDC + (I-1)*MB + 1
                  IPBKK = (K-1)*NB*LDB + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( B( IPBKK+(NB-1)*LDB+NB ) .NE. 0.0D+00 ) THEN
                        BRWS = BRWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( B( IPBKK-LDB ) .NE. 0.0D+00 ) THEN
                        BRWS = BRWS - 1
                        IPB = IPB + 1
                        IPC = IPC + LDC
                     END IF
                  END IF
                  IF( BDIM ) IPB = IPB + LDB
                  IF( ADIM ) IPC = IPC + 1
                  IF( BRWS.GT.0 )
     $                 CALL DGEMM( 'N', 'T', ROWS, BRWS, COLS, USIGN, 
     $                 C(IPX), LDC, B(IPB), LDB, 1.0D+00, C(IPC), LDC )
 40            CONTINUE
            END IF
         END IF
      END IF
C
      END
C     
C     End of DTRSYCT
C     
C *** Last line of DTRSYCT ***
