CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRLYCT( JOB, OP, I, J, N, NB, A, LDA, C, LDC, DWORK, 
     $                    LDWORK, ROWS, COLS, XRST, XCST, SCALE, INFO )
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
      CHARACTER          JOB, OP
      INTEGER            I, J, NB, INFO, LDA, LDC, N, ROWS, COLS, XRST, 
     $                   XCST, LDWORK
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), C( * ), DWORK( * )
C     ..
C
C     Purpose and description
C     =======================
C
C     Blocked LYCT solver, which solves for the (I,J)th
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
C     It is assumed that the solution is computed explicitly in the lower 
C     triangular part of X. If this routine is called with I and J 
C     pointing to any part of X that is not in the explicitly computed 
C     part, this routine will only compute the output information ROWS, 
C     COLS, XRST and XCST and return immediately.
C
      LOGICAL SOLVE, UPDATE, RSIDE, AEXT1, ADIM1, AEXT2, ADIM2
      INTEGER IPA1, IPA2, IPC, K, IPX, ARWS, ACLS, IPAKK, UPLOSIGN, 
     $     IPXT, THREADS
      DOUBLE PRECISION USIGN, MACHINE(10)
C     
      LOGICAL  LSAME
      INTEGER  ICEIL, OMP_GET_NUM_THREADS
      EXTERNAL LSAME, ICEIL, OMP_GET_NUM_THREADS
C     
      SOLVE  = LSAME(JOB,'S') .OR. LSAME(JOB,'B')
      UPDATE = LSAME(JOB,'U') .OR. LSAME(JOB,'B')
      RSIDE = LSAME(OP,'N')
      USIGN = -1.0D+00
C
      IF( I.GT.ICEIL(N,NB) .OR. J.GT.ICEIL(N,NB) ) THEN
         ROWS = 0
         COLS = 0
         XRST = 1 
         XCST = 1
         RETURN
      END IF
      IF( LDWORK.LT.(NB+1)**2 ) THEN
         INFO = -12
         RETURN
      END IF
C
      IF( SOLVE .OR. UPDATE ) THEN
         IF( NB.GT.N ) THEN
            ROWS = N
         ELSE
            ROWS = MIN( NB, N - (I-1)*NB )
         END IF
         IF( NB.GT.N ) THEN
            COLS = N
         ELSE
            COLS = MIN( NB, N - (J-1)*NB )
         END IF
         IPA1 = (I-1)*NB*LDA + (I-1)*NB + 1
         IPA2 = (J-1)*NB*LDA + (J-1)*NB + 1
         IPC  = (J-1)*NB*LDC + (I-1)*NB + 1
         IF( I.LT.ICEIL(N,NB) ) THEN
            IF( A(IPA1+(NB-1)*LDA+NB) .NE. 0.0D+00 ) THEN
               ROWS = ROWS + 1
               AEXT1 = .TRUE.
            ELSE
               AEXT1 = .FALSE.
            END IF
         ELSE
            AEXT1 = .FALSE. 
         END IF
         IF( I.GT.1 ) THEN
            IF( A(IPA1-LDA) .NE. 0.0D+00 ) THEN
               IPA1 = IPA1 + LDA + 1
               IPC = IPC + 1
               ROWS = ROWS - 1
               ADIM1 = .TRUE.
            ELSE
               ADIM1 = .FALSE.
            END IF
         ELSE
            ADIM1 = .FALSE. 
         END IF
         IF( I.NE.J ) THEN
            IF( J.LT.ICEIL(N,NB) ) THEN
               IF( A(IPA2+(NB-1)*LDA+NB) .NE. 0.0D+00 ) THEN
                  COLS = COLS + 1
                  AEXT2 = .TRUE.
               ELSE
                  AEXT2 = .FALSE.
               END IF
            ELSE
               AEXT2 = .FALSE.  
            END IF
            IF( J.GT.1 ) THEN
               IF( A(IPA2-LDA) .NE. 0.0D+00 ) THEN
                  IPA2 = IPA2 + LDA + 1
                  IPC = IPC + LDC
                  COLS = COLS - 1
                  ADIM2 = .TRUE.
               ELSE
                  ADIM2 = .FALSE.
               END IF
            ELSE
               ADIM2 = .FALSE.  
            END IF
         ELSE
            COLS = ROWS
            AEXT2 = AEXT1
            ADIM2 = ADIM1
            IF( ADIM2 ) IPC = IPC + LDC
         END IF
      END IF
C
      IF( ROWS.EQ.0 .OR. COLS.EQ.0 ) RETURN
C
      XRST = MOD(IPC-1,LDC) + 1
      XCST = ICEIL(IPC, LDC)
C
      IF( I.LT.J ) RETURN
C
      IPXT = (XRST-1)*LDC + XCST
      IF( SOLVE ) THEN
         IF( I.EQ.J ) THEN
            IF( RSIDE ) THEN
               UPLOSIGN = 0
            ELSE
               UPLOSIGN = 1
            END IF
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
            THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
            CALL RECSY_MACHINE( MACHINE )
            CALL RECLYCT_P( THREADS, UPLOSIGN, SCALE, ROWS, A(IPA1), 
     $                      LDA, C(IPC), LDC, INFO, MACHINE )
#else
            CALL RECSY_MACHINE( MACHINE )
            CALL RECLYCT( UPLOSIGN, SCALE, ROWS, A(IPA1), LDA, C(IPC), 
     $                    LDC, INFO, MACHINE )
#endif
            CALL DLATCPY( 'Upper', ROWS-1, ROWS-1, C(IPC+LDC), LDC,
     $                    C(IPC+1), LDC )
         ELSE
            IF ( RSIDE ) THEN
               UPLOSIGN = 2
            ELSE
               UPLOSIGN = 4
            END IF
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
            THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
            CALL RECSY_MACHINE( MACHINE )
            CALL RECSYCT_P( THREADS, UPLOSIGN, SCALE, ROWS, COLS, 
     $                      A(IPA1), LDA, A(IPA2), LDA, C(IPC), LDC, 
     $                      INFO, MACHINE )
#else
            CALL RECSY_MACHINE( MACHINE )
            CALL RECSYCT( UPLOSIGN, SCALE, ROWS, COLS, A(IPA1), LDA, 
     $                    A(IPA2), LDA, C(IPC), LDC, INFO, MACHINE )
#endif
            CALL DLATCPY( 'All', ROWS, COLS, C(IPC), LDC, C(IPXT), LDC )
         END IF
      END IF
C
      IF( UPDATE ) THEN
         IPX = IPC
         IF( NB.LT.N ) THEN
            IF( RSIDE ) THEN
               DO 10 K = J, I-1
                  ARWS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (I-1)*NB*LDA + (K-1)*NB + 1
                  IPC = (J-1)*NB*LDC + (K-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS - 1
                        IPA1 = IPA1 + 1
                        IPC = IPC + 1  
                     END IF
                  END IF
                  IF( ADIM1 ) IPA1 = IPA1 + LDA
                  IF( ADIM2 ) IPC = IPC + LDC
                  IF( ARWS.GT.0 ) THEN
                     IF( K.EQ.J ) THEN
                        CALL DSYR2K( 'U', 'N', COLS, ROWS, -1.0D+00, 
     $                       A(IPA1), LDA, C(IPXT), LDC, 1.0D+00, 
     $                       C(IPC), LDC )
                        CALL DLATCPY( 'Upper', COLS-1, COLS-1, 
     $                       C(IPC+LDC), LDC, C(IPC+1), LDC )
                     ELSE
                        CALL DGEMM( 'N', 'N', ARWS, COLS, ROWS, 
     $                       -1.0D+00, A(IPA1), LDA, C(IPX), LDC, 
     $                       1.0D+00, C(IPC), LDC )
                     END IF
                  END IF
 10            CONTINUE
            ELSE
               DO 20 K = I+1, ICEIL(N,NB)
                  ACLS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (K-1)*NB*LDA + (I-1)*NB + 1
                  IPC = (J-1)*NB*LDC + (K-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS - 1
                        IPA1 = IPA1 + LDA
                        IPC = IPC + 1
                     END IF
                  END IF
                  IF( ADIM1 ) IPA1 = IPA1 + 1
                  IF( ADIM2 ) IPC = IPC + LDC
                  IF( ACLS.GT.0 ) 
     $                 CALL DGEMM( 'T', 'N', ACLS, COLS, ROWS, -1.0D+00, 
     $                 A(IPA1), LDA, C(IPX), LDC, 1.0D+00, C(IPC), LDC )
 20            CONTINUE
            END IF
         END IF
C     
         IF( NB.LT.N ) THEN
            IF( .NOT. RSIDE ) THEN
               DO 30 K = J+1, I
                  ACLS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (K-1)*NB*LDA + (J-1)*NB + 1
                  IPC = (K-1)*NB*LDC + (I-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS - 1
                        IPA1 = IPA1 + LDA
                        IPC = IPC + LDC
                     END IF
                  END IF
                  IF( ADIM2 ) IPA1 = IPA1 + 1
                  IF( ADIM1 ) IPC = IPC + 1
                  IF( ACLS.GT.0 ) THEN
                     IF( K.EQ.I ) THEN
                        CALL DLATCPY( 'All', COLS, ROWS, A(IPA1), LDA,
     $                       DWORK, ROWS )
                        CALL DSYR2K( 'U', 'N', ROWS, COLS, -1.0D+00, 
     $                       DWORK, ROWS, C(IPX), LDC, 1.0D+00, 
     $                       C(IPC), LDC )
                        CALL DLATCPY( 'Upper', ROWS-1, ROWS-1, 
     $                       C(IPC+LDC), LDC, C(IPC+1), LDC )
                     ELSE
                        CALL DGEMM( 'N', 'N', ROWS, ACLS, COLS, USIGN, 
     $                       C(IPX), LDC, A(IPA1), LDA, 1.0D+00, C(IPC), 
     $                       LDC )
                     END IF
                  END IF
 30            CONTINUE
            ELSE
               DO 40 K = 1, J-1
                  ARWS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (J-1)*NB*LDA + (K-1)*NB + 1
                  IPC = (K-1)*NB*LDC + (I-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS - 1
                        IPA1 = IPA1 + 1
                        IPC = IPC + LDC
                     END IF
                  END IF
                  IF( ADIM2 ) IPA1 = IPA1 + LDA
                  IF( ADIM1 ) IPC = IPC + 1
                  IF( ARWS.GT.0 )
     $                 CALL DGEMM( 'N', 'T', ROWS, ARWS, COLS, USIGN, 
     $                 C(IPX), LDC, A(IPA1), LDA, 1.0D+00, C(IPC), LDC )
 40            CONTINUE
            END IF
         END IF
C
         IF( NB.LT.N ) THEN
            IF( RSIDE .AND. I.GT.J ) THEN
               DO 50 K = 1, J-1
                  ARWS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (I-1)*NB*LDA + (K-1)*NB + 1
                  IPC = (K-1)*NB*LDC + (J-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ARWS = ARWS - 1
                        IPA1 = IPA1 + 1
                        IPC = IPC + LDC
                     END IF
                  END IF
                  IF( ADIM1 ) IPA1 = IPA1 + LDA
                  IF( ADIM2 ) IPC = IPC + 1
                  IF( ARWS.GT.0 )
     $                 CALL DGEMM( 'N', 'T', COLS, ARWS, ROWS, USIGN, 
     $                 C(IPXT), LDC, A(IPA1), LDA, 1.0D+00, C(IPC), LDC)
 50            CONTINUE
            ELSEIF( I.GT.J ) THEN
               DO 60 K = I+1, ICEIL(N,NB)
                  ACLS = MIN( NB, N - (K-1)*NB )
                  IPA1 = (K-1)*NB*LDA + (J-1)*NB + 1
                  IPC = (I-1)*NB*LDC + (K-1)*NB + 1
                  IPAKK = (K-1)*NB*LDA + (K-1)*NB + 1
                  IF( K.LT.ICEIL(N,NB) ) THEN
                     IF( A( IPAKK+(NB-1)*LDA+NB ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACLS = ACLS - 1
                        IPA1 = IPA1 + LDA
                        IPC = IPC + 1
                     END IF
                  END IF
                  IF( ADIM2 ) IPA1 = IPA1 + 1
                  IF( ADIM1 ) IPC = IPC + LDC
                  IF( ACLS.GT.0 )
     $                 CALL DGEMM( 'T', 'N', ACLS, ROWS, COLS, USIGN, 
     $                 A(IPA1), LDA, C(IPXT), LDC, 1.0D+00, C(IPC), LDC)
 60            CONTINUE
            END IF
         END IF
      END IF
C
      END
C     
C     End of DTRLYCT
C     
C *** Last line of DTRLYCT ***
