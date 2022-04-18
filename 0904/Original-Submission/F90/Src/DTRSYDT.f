CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRSYDT( JOB, TRANA, TRANB, ISGN, I, M, MB, N,  A, 
     $     LDA, B, LDB, C, LDC, WORK, LWORK, ROWS, XRST, SCALE, INFO )
C
C  -- LAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     October 18, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          JOB, TRANA, TRANB
      INTEGER            I, J, MB, INFO, ISGN, LDA, LDB, LDC, M, N,
     $                   ROWS, XRST, LWORK
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), B( * ), C( * ), WORK( * )
C     ..
C
C     Purpose and description
C     =======================
C
C     Blocked SYDT solver, which solves for the Ith
C     MB-by-N subsystem or performs GEMM-updates with respect to
C     the already solved Ith subsystem. What to do is
C     controlled by JOB (which takes the values 'S' for solve, 'U'
C     for update). Notice that if JOB = 'U', a previous call with 
C     JOB = 'S' for the same value of I is assumed such that 
C     the variables ROWS and XRST ais input from that 
C     previuos call. Moreover, if this routine returned with
C     SCALE < 1.0 for JOB = 'S', the user must take care of the
C     scaling of the right hand side (excluding the (I,J)th 
C     subsolution) before re-calling this routine with JOB = 'U'. 
C
      LOGICAL SOLVE, UPDATE, TRANSA, TRANSB, AEXT, ADIM 
      INTEGER IPA, IPC, K, IPX, ARWS, ACLS, IPAKK, LWRK1, DBA, LDE, 
     $     IEND, ISTART, IPBJJ, II, ROWS2, ERWS, UPLOSIGN, THREADS
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
      USIGN = -1.0D+00
C     
      IF( I.LT.1 .OR. I.GT.ICEIL( M,MB ) ) THEN
         ROWS = 0
         XRST = 1 
         RETURN
      END IF
      LWRK1 = N*(MB+1)
      IF( LWORK.LT.LWRK1  ) THEN
         INFO = -18
         RETURN
      END IF
C     
      IF( SOLVE .OR. UPDATE ) THEN
         IF( MB.GT.M ) THEN
            ROWS = M
         ELSE
            ROWS = MIN( MB, M - (I-1)*MB )
         END IF
         IPA = (I-1)*MB*LDA + (I-1)*MB + 1
         IPC = (I-1)*MB + 1
         IF( I.LT.ICEIL(M,MB) ) THEN
            IF( A(IPA+(MB-1)*LDA+MB) .NE. 0.0D+00 ) THEN
               ROWS = ROWS + 1
            END IF 
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
      END IF
C     
      XRST = IPC
C
      IF( ROWS.EQ.0 ) RETURN
C     
      IF( SOLVE ) THEN
         UPLOSIGN = (1-ISGN) / 2
         IF( TRANSA ) UPLOSIGN = UPLOSIGN + 4
         IF( TRANSB ) UPLOSIGN = UPLOSIGN + 2
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
         THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
         CALL RECSY_MACHINE( MACHINE )
         CALL RECSYDT_P( THREADS, UPLOSIGN, SCALE, ROWS, N, A(IPA), LDA, 
     $                   B, LDB, C(IPC), LDC, INFO, MACHINE, WORK, 
     $                   LWORK )
#else
         CALL RECSY_MACHINE( MACHINE )
         CALL RECSYDT( UPLOSIGN, SCALE, ROWS, N, A(IPA), LDA, 
     $                 B, LDB, C(IPC), LDC, INFO, MACHINE, WORK, LWORK )
#endif
      END IF
C     
      IF( UPDATE ) THEN
         IPX = IPC
         DBA = ICEIL(M,MB)
         IF( MB.LT.M ) THEN
            IF( .NOT. TRANSA ) THEN
               DO 10 K = 1, I-1
                  ARWS = MIN( MB, M - (K-1)*MB )
                  IPA = (I-1)*MB*LDA + (K-1)*MB + 1
                  IPC = (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.DBA ) THEN
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
                  IF( ARWS.GT.0 )
     $                 CALL DGEMM( 'N', 'N', ARWS, N, ROWS, 1.0D+00, 
     $                 A(IPA), LDA, C(IPX), LDC, 0.0D+00, WORK, ARWS )
                  IF( ARWS.GT.0 )
     $                 CALL DGEMM( 'N', TRANB, ARWS, N, N, -1.0D+00, 
     $                 WORK, ARWS, B, LDB, 1.0D+00, C(IPC), LDC )
 10            CONTINUE
            ELSEIF( TRANSA ) THEN
               DO 20 K = I+1, DBA
                  ACLS = MIN( MB, M - (K-1)*MB )
                  IPA = (K-1)*MB*LDA + (I-1)*MB + 1
                  IPC = (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.DBA ) THEN
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
                  IF( ACLS.GT.0 )
     $                 CALL DGEMM( 'T', 'N', ACLS, N, ROWS, 1.0D+00, 
     $                 A(IPA), LDA, C(IPX), LDC, 0.0D+00, WORK, ACLS )
                  IF( ACLS.GT.0 )
     $                 CALL DGEMM( 'N', TRANB, ACLS, N, N, -1.0D+00, 
     $                 WORK, ACLS, B, LDB, 1.0D+00, C(IPC), LDC )
 20            CONTINUE
            END IF
         END IF
      END IF
C     
      END
C     
C     End of DTRSYDT
C     
C *** Last line of DTRSYDT ***
