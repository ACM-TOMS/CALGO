CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         This file is a part of the SCASY software library.           C
C         See http://www.cs.umu.se/parallel/research/scasy             C
C         Contributors: Robert Granat and Bo Kågström.                 C
C         SCASY Version 0.10,    March 31, 2006.                       C
C         Copyright 2006, Umeå University, Sweden.                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DTRGSYL( JOB, TRANAC, TRANBD, ISGN, I, M, MB, N,  A, 
     $     LDA, B, LDB, C, LDC, D, LDD, E, LDE, WORK, LWORK, ROWS, XRST, 
     $     SCALE, INFO )
C
C  -- LAPACK-style routine (preliminary version)--
C     HPC2N and Department of Computing Science,
C     University of Umeå, Sweden.
C     Written by Robert Granat, (granat@cs.umu.se)
C     October 20, 2006.
C
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      CHARACTER          JOB, TRANAC, TRANBD
      INTEGER            I, J, MB, INFO, ISGN, LDA, LDB, LDC, LDD, LDE,
     $                   M, N, ROWS, XRST, LWORK
      DOUBLE PRECISION   SCALE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), E( * ),
     $                   WORK( * )
C     ..
C
C     Purpose and description
C     =======================
C
C     Blocked GSYL solver, which solves for the Ith
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
      LOGICAL SOLVE, UPDATE, TRANSAC, TRANSBD, ACDIM 
      INTEGER IPAC, IPE, K, IPX, ACRWS, ACCLS, IPAKK, LWRK1, DBAC,  
     $        UPLOSIGN, THREADS
      DOUBLE PRECISION USIGN, MACHINE(10)
C     
      LOGICAL  LSAME, OMO_GET_NUM_THREADS
      INTEGER  ICEIL, OMP_GET_NUM_THREADS
      EXTERNAL LSAME, ICEIL
C     
      SOLVE  = LSAME(JOB,'S') 
      UPDATE = LSAME(JOB,'U')
      TRANSAC = LSAME(TRANAC,'T')
      TRANSBD = LSAME(TRANBD,'T')
      USIGN = -DBLE(ISGN)
C     
      IF( I.LT.1 .OR. I.GT.ICEIL( M,MB ) ) THEN
         ROWS = 0
         XRST = 1 
         RETURN
      END IF
      LWRK1 = N*(MB+1)
      IF( LWORK.LT.LWRK1  ) THEN
         INFO = -20
         RETURN
      END IF
C     
      IF( SOLVE .OR. UPDATE ) THEN
         IF( MB.GT.M ) THEN
            ROWS = M
         ELSE
            ROWS = MIN( MB, M - (I-1)*MB )
         END IF
         IPAC = (I-1)*MB*LDA + (I-1)*MB + 1
         IPE = (I-1)*MB + 1
         IF( I.LT.ICEIL(M,MB) ) THEN
            IF( A(IPAC+(MB-1)*LDA+MB) .NE. 0.0D+00 ) THEN
               ROWS = ROWS + 1
            END IF 
         END IF
         IF( I.GT.1 ) THEN
            IF( A(IPAC-LDA) .NE. 0.0D+00 ) THEN
               IPAC = IPAC + LDA + 1
               IPE = IPE + 1
               ROWS = ROWS - 1
               ACDIM = .TRUE.
            ELSE
               ACDIM = .FALSE.
            END IF
         ELSE
            ACDIM = .FALSE. 
         END IF
      END IF
C     
      XRST = IPE
C
      IF( ROWS.EQ.0 ) RETURN
C     
      IF( SOLVE ) THEN
         UPLOSIGN = (1-ISGN) / 2
         IF( TRANSAC ) UPLOSIGN = UPLOSIGN + 4
         IF( TRANSBD ) UPLOSIGN = UPLOSIGN + 2
#ifdef USE_OMP
C$OMP PARALLEL SHARED(THREADS)
C$OMP MASTER
            THREADS = OMP_GET_NUM_THREADS()
C$OMP END MASTER
C$OMP END PARALLEL 
         CALL RECSY_MACHINE( MACHINE )
         CALL RECGSYL_P( THREADS, UPLOSIGN, SCALE, ROWS, N, A(IPAC), 
     $                   LDA, B, LDB, C(IPAC), LDC, D, LDD, E(IPE),
     $                   LDE, INFO, MACHINE, WORK, LWORK )
#else
         CALL RECSY_MACHINE( MACHINE )
         CALL RECGSYL( UPLOSIGN, SCALE, ROWS, N, A(IPAC), LDA, 
     $                 B, LDB, C(IPAC), LDC, D, LDD, E(IPE),
     $                 LDE, INFO, MACHINE, WORK, LWORK )
#endif
      END IF
C     
      IF( UPDATE ) THEN
         IPX = IPE
         DBAC = ICEIL(M,MB)
         IF( MB.LT.M ) THEN
            IF( .NOT. TRANSAC ) THEN
               DO 10 K = 1, I-1
                  ACRWS = MIN( MB, M - (K-1)*MB )
                  IPAC = (I-1)*MB*LDA + (K-1)*MB + 1
                  IPE = (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.DBAC ) THEN
                     IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                        ACRWS = ACRWS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACRWS = ACRWS - 1
                        IPAC = IPAC + 1
                        IPE = IPE + 1
                     END IF
                  END IF
                  IF( ACDIM ) IPAC = IPAC + LDA
                  IF( ACRWS.GT.0 )
     $                 CALL DGEMM( 'N', 'N', ACRWS, N, ROWS, 1.0D+00, 
     $                 A(IPAC), LDA, E(IPX), LDE, 0.0D+00, WORK, ACRWS )
                  IF( ACRWS.GT.0 )
     $                 CALL DGEMM( 'N', TRANBD, ACRWS, N, N, -1.0D+00, 
     $                 WORK, ACRWS, B, LDB, 1.0D+00, E(IPE), LDE )
                  IF( ACRWS.GT.0 )
     $                 CALL DGEMM( 'N', 'N', ACRWS, N, ROWS, 1.0D+00, 
     $                 C(IPAC), LDC, E(IPX), LDE, 0.0D+00, WORK, ACRWS )
                  IF( ACRWS.GT.0 )
     $                 CALL DGEMM( 'N', TRANBD, ACRWS, N, N, USIGN, 
     $                 WORK, ACRWS, D, LDD, 1.0D+00, E(IPE), LDE )
 10            CONTINUE
            ELSEIF( TRANSAC ) THEN
               DO 20 K = I+1, DBAC
                  ACCLS = MIN( MB, M - (K-1)*MB )
                  IPAC = (K-1)*MB*LDA + (I-1)*MB + 1
                  IPE = (K-1)*MB + 1
                  IPAKK = (K-1)*MB*LDA + (K-1)*MB + 1
                  IF( K.LT.DBAC ) THEN
                     IF( A( IPAKK+(MB-1)*LDA+MB ) .NE. 0.0D+00 ) THEN
                        ACCLS = ACCLS + 1
                     END IF
                  END IF
                  IF( K.GT.1 ) THEN
                     IF( A( IPAKK-LDA ) .NE. 0.0D+00 ) THEN
                        ACCLS = ACCLS - 1
                        IPAC = IPAC + LDA
                        IPE = IPE + 1
                     END IF
                  END IF
                  IF( ACDIM ) IPAC = IPAC + 1
                  IF( ACCLS.GT.0 )
     $                 CALL DGEMM( 'T', 'N', ACCLS, N, ROWS, 1.0D+00, 
     $                 A(IPAC), LDA, E(IPX), LDE, 0.0D+00, WORK, ACCLS )
                  IF( ACCLS.GT.0 )
     $                 CALL DGEMM( 'N', TRANBD, ACCLS, N, N, -1.0D+00, 
     $                 WORK, ACCLS, B, LDB, 1.0D+00, E(IPE), LDE )
                  IF( ACCLS.GT.0 )
     $                 CALL DGEMM( 'T', 'N', ACCLS, N, ROWS, 1.0D+00, 
     $                 C(IPAC), LDC, E(IPX), LDE, 0.0D+00, WORK, ACCLS )
                  IF( ACCLS.GT.0 )
     $                 CALL DGEMM( 'N', TRANBD, ACCLS, N, N, USIGN, 
     $                 WORK, ACCLS, D, LDD, 1.0D+00, E(IPE), LDE )
 20            CONTINUE
            END IF
         END IF
      END IF
C     
      END
C     
C     End of DTRGSYL
C     
C *** Last line of DTRGSYL ***
