*
*     Auto-tuning routine
*
*     Find suitable PILAENVX(ISPEC=14) based on experiments.
*
*     This program assume that the following experiments are performed:
*         2x2 grid, N = 4000, 8000
*         4x4 grid, N = 4000, 8000, 16000
*         8x8 grid, N = 4000, 8000, 16000, 32000
*     with NIBBLE=10,20,...,80,90 for each case.
*
      PROGRAM TUNE2
      IMPLICIT NONE
*
      INTEGER TOTAL
      INTEGER SOLVER, NIBBLE, N, NB, NPROW, NPCOL, NIBBLE0
      DOUBLE PRECISION AA, BB, CC, AA0, BB0, CC0
      INTEGER NMAT(4), NPROC(3)
      DOUBLE PRECISION TIME
      DOUBLE PRECISION DATASET(9,4,3)
      DOUBLE PRECISION A(9,3), B(9), WORK(1024)
      INTEGER STAT
      CHARACTER*(16) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7
      CHARACTER*(512) LINE
      INTEGER N2IND, NP2IND, MININD
      DOUBLE PRECISION MINPARAM
      DOUBLE PRECISION LOG2
      DOUBLE PRECISION DIFF, TT
      INTEGER I, J, INFO
*     Threshold for replacing the parameters.
      DOUBLE PRECISION THRESH
      PARAMETER ( THRESH = 1.0D-1 )
*
      INTRINSIC DBLE, NINT, EXP, LOG
*
      EXTERNAL DGELS
      EXTERNAL GETNIBBLE
*
*     Setting for the experiments.
*
      NMAT(1) = 4000
      NMAT(2) = 8000
      NMAT(3) = 16000
      NMAT(4) = 32000
      NPROC(1) = 2
      NPROC(2) = 4
      NPROC(3) = 8
*
      LOG2 = LOG(2.0D0)
      A(1,1) = 4.0
      A(2,1) = 4.0
      A(3,1) = 16.0
      A(4,1) = 16.0
      A(5,1) = 16.0
      A(6,1) = 64.0
      A(7,1) = 64.0
      A(8,1) = 64.0
      A(9,1) = 64.0
      A(1,2) = 4.0*LOG(4.0D3/2.0)
      A(2,2) = 4.0*LOG(8.0D3/2.0)
      A(3,2) = 16.0*LOG(4.0D3/4.0)
      A(4,2) = 16.0*LOG(8.0D3/4.0)
      A(5,2) = 16.0*LOG(1.6D4/4.0)
      A(6,2) = 64.0*LOG(4.0D3/8.0)
      A(7,2) = 64.0*LOG(8.0D3/8.0)
      A(8,2) = 64.0*LOG(1.6D4/8.0)
      A(9,2) = 64.0*LOG(3.2D4/8.0)
      A(1,3) = 4.0*LOG2
      A(2,3) = 4.0*LOG2
      A(3,3) = 16.0*2.0*LOG2
      A(4,3) = 16.0*2.0*LOG2
      A(5,3) = 16.0*2.0*LOG2
      A(6,3) = 64.0*3.0*LOG2
      A(7,3) = 64.0*3.0*LOG2
      A(8,3) = 64.0*3.0*LOG2
      A(9,3) = 64.0*3.0*LOG2
*
*     Read the experimental data.
*
      TOTAL = 0
 1    CONTINUE
 2       CONTINUE
*           Skip comments
            READ( *, '(A)', IOSTAT=STAT ) LINE
            IF( STAT.LT.0 ) GO TO 3
            IF( INDEX(LINE, '#').NE.0 )
     $   GO TO 2
         READ( LINE, FMT = * )
     $        DUM1, DUM2, DUM3, DUM4, DUM5, DUM6, DUM7,
     $        SOLVER, NIBBLE, N, NB, NPROW, NPCOL, TIME
         DATASET(NIBBLE/10, N2IND(N), NP2IND(MIN(NPROW, NPCOL))) = TIME
         TOTAL = TOTAL+1
         GO TO 1
 3    CONTINUE
      IF( TOTAL.LT.81 ) THEN
         WRITE(*,*) 'Warning: Too few experiments.'
         WRITE(*,*) 'The result might be unreliable.'
         WRITE(*,*)
      END IF
*
*     Find suitable parameters by least squares fitting.
*
      B(1) = 4.0*MINPARAM(9, DATASET(1, N2IND(4000), NP2IND(2)))
      B(2) = 4.0*MINPARAM(9, DATASET(1, N2IND(8000), NP2IND(2)))
      B(3) = 16.0*MINPARAM(9, DATASET(1, N2IND(4000), NP2IND(4)))
      B(4) = 16.0*MINPARAM(9, DATASET(1, N2IND(8000), NP2IND(4)))
      B(5) = 16.0*MINPARAM(9, DATASET(1, N2IND(16000), NP2IND(4)))
      B(6) = 64.0*MINPARAM(9, DATASET(1, N2IND(4000), NP2IND(8)))
      B(7) = 64.0*MINPARAM(9, DATASET(1, N2IND(8000), NP2IND(8)))
      B(8) = 64.0*MINPARAM(9, DATASET(1, N2IND(16000), NP2IND(8)))
      B(9) = 64.0*MINPARAM(9, DATASET(1, N2IND(32000), NP2IND(8)))
      CALL DGELS( 'N', 9, 3, 1, A, 9, B, 9, WORK, 1024, INFO)
      AA = exp(B(1))
      BB = B(2)
      CC = (B(3)-B(2))*0.5
*
*     Comparison between the current setting with the measured one.
*
      CALL GETNIBBLE(AA0, BB0, CC0)
      DIFF = 0.0
      TT = 0.0
      DO I = 1, 4
      DO J = 1, 3
         IF(NMAT(I)/NPROC(J).LE.4000) THEN
            NIBBLE = MIN(90, NINT(AA*(NMAT(I))**BB*NPROC(J)**CC))
            NIBBLE = MAX((NIBBLE+5)/10, 1)
            NIBBLE0 = MIN(90, NINT(AA0*(NMAT(I))**BB0*NPROC(J)**CC0))
            NIBBLE0 = MAX((NIBBLE0+5)/10, 1)
            DIFF = DIFF + (DATASET(NIBBLE0, I, J)-DATASET(NIBBLE, I, J))
            TOTAL = TOTAL + DATASET(MININD(9, DATASET(1, I, J)), I, J)
         END IF
      END DO
      END DO
      WRITE(*,*), 'Suggestions on PILAENVX(ISPEC=14):'
      WRITE(*,*), '    PIPARMQ = AA * NH**BB * MIN(NPROW, NPCOL)**CC'
      WRITE(*,*), 'Current values:'
      WRITE(*,*), '    AA =', AA0
      WRITE(*,*), '    BB =', BB0
      WRITE(*,*), '    CC =', CC0
      WRITE(*,*), 'Computed values based on experiments:'
      WRITE(*,*), '    AA =', AA
      WRITE(*,*), '    BB =', BB
      WRITE(*,*), '    CC =', CC
      IF( DIFF/TOTAL.LE.THRESH ) THEN
         WRITE(*,*), 'The current setting is reasonably good.'
         WRITE(*,*), 'We recommend to keep the current setting.'
      ELSE
         WRITE(*,*), 'The current setting is not good.'
         WRITE(*,*), 'We recommend to replace the parameters ',
     $        'by the computed ones.'
      END IF
*
      END
*
      INTEGER FUNCTION N2IND(N)
      INTEGER N
*
      IF( N.EQ.4000 ) THEN
         N2IND = 1
      ELSE IF( N.EQ.8000 ) THEN
         N2IND = 2
      ELSE IF( N.EQ.16000 ) THEN
         N2IND = 3
      ELSE
         N2IND = 4
      END IF
*
      RETURN
      END
*
      INTEGER FUNCTION NP2IND(NP)
      INTEGER NP
*
      IF( NP.EQ.2 ) THEN
         NP2IND = 1
      ELSE IF( NP.EQ.4 ) THEN
         NP2IND = 2
      ELSE
         NP2IND = 3
      END IF
*
      RETURN
      END
*
      INTEGER FUNCTION MININD(N, X)
      INTEGER N
      DOUBLE PRECISION X(N)
*
      INTEGER I, J
*
      J = 1
      DO I = 2, N
         IF( X(I).LT.X(J) ) J = I
      END DO
      MININD = J
*
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION MINPARAM(N, X)
      INTEGER N
      DOUBLE PRECISION X(N)
*
      INTRINSIC LOG
*
      INTEGER MININD
      INTEGER I, L
*
      L = 0
      DO I = 2, N-1
         IF( X(I).LT.X(I-1) .AND. X(I).LT.X(I+1) ) L = L+1
      END DO
      IF( L.LE.1 ) THEN
         MINPARAM = MININD(N, X)
      ELSE
         I = 1
         DO
            IF( I.GE.N ) EXIT
            IF( X(I).LE.X(I+1) ) EXIT
            I = I+1
         END DO
         MINPARAM = I
      END IF
      MINPARAM = LOG(10.0D0*MINPARAM)
*
      RETURN
      END
