*
*     Auto-tuning routine
*
*     Find suitable PILAENVX(ISPEC=12) and PILAENVX(ISPEC=23) based on
*     experiments.
*
*     This program assume that the following experiments are performed:
*         N = 96, 192, 384, 768, 1536, 3072, 6144
*     on processor grids 1x1, 2x2, 4x4, 8x8 for each case.
*
      PROGRAM TUNE1
      IMPLICIT NONE
*
      INTEGER TOTAL
      INTEGER SOLVER, N, NB, NPROW, NPCOL
      INTEGER NMIN, PMIN, CURRENT
      INTEGER NMINS(4), PMINS(7)
      INTEGER NMAT(7), NPROC(4)
      DOUBLE PRECISION TIME
      DOUBLE PRECISION DATASET(3,7,4)
      INTEGER STAT
      CHARACTER*(16) DUM1, DUM2, DUM3, DUM4, DUM5, DUM6
      CHARACTER*(512) LINE
      INTEGER N2IND, NP2IND
      INTEGER I, J
*
      INTRINSIC DBLE, NINT
*
      INTEGER ICTXT
      EXTERNAL BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER PILAENVX
      EXTERNAL PILAENVX
*
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, '2D', 1, 1 )
*
*     Setting for the experiments.
*
      NMAT(1) = 96
      NMAT(2) = 192
      NMAT(3) = 384
      NMAT(4) = 768
      NMAT(5) = 1536
      NMAT(6) = 3072
      NMAT(7) = 6144
      NPROC(1) = 1
      NPROC(2) = 2
      NPROC(3) = 4
      NPROC(4) = 8
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
     $        DUM1, DUM2, DUM3, DUM4, DUM5, DUM6,
     $        SOLVER, N, NB, NPROW, NPCOL, TIME
         DATASET(SOLVER, N2IND(N), NP2IND(MIN(NPROW, NPCOL))) = TIME
         TOTAL = TOTAL+1
         GO TO 1
 3    CONTINUE
      IF( TOTAL.LT.56 ) THEN
         WRITE(*,*) 'Warning: Too few experiments.'
         WRITE(*,*) 'The result might be unreliable.'
         WRITE(*,*)
      END IF
*
*     Seek for the crossover point between PDLAQR1/PDLAQR0.
*
      DO J = 1,4
         NMINS(J) = NMAT(7)
         DO I = 1,7
            IF( DATASET(2,I,J).LE.DATASET(1,I,J) ) EXIT
         END DO
         NMINS(J) = NMAT(I)/NPROC(J)
      END DO
      NMIN = (NPROC(1)*NMINS(1)+NPROC(2)*NMINS(2)+NPROC(3)*NMINS(3)
     $     +NPROC(4)*NMINS(4)) / (NPROC(1)+NPROC(2)+NPROC(3)+NPROC(4))
      WRITE(*,*), 'Suggestions on PILAENVX(ISPEC=12) and ',
     $     'PIPARMQ(ISPEC=12):'
      WRITE(*,*), '    PIPARMQ = NMIN * MIN( NPROW, NPCOL )'
      WRITE(*,*), 'with', MAX(NB, NMIN/2), '<= NMIN <=', NMIN
      CURRENT = PILAENVX(ICTXT, 12, '', '', 0, 0, 0, 0)
      IF( MAX(NB, NMIN/2).LE.CURRENT .AND. CURRENT.LE.NMIN ) THEN
         WRITE(*,*), 'The current value ( NMIN =', CURRENT,
     $        ') is acceptable.'
      ELSE
         WRITE(*,*), 'The current value ( NMIN =', CURRENT,
     $        ') needs to be replaced.'
      END IF
      WRITE(*,*)
*
*     Seek for the minimal memory load which maintains scalability.
*
      NMIN = NMIN/4*3
      DO I = 1,7
         DO J = 1,4
            IF( NMAT(I)>NMIN*NPROC(J) ) THEN
               DATASET(3,I,J) = DATASET(2,I,J)
            ELSE
               DATASET(3,I,J) = DATASET(1,I,J)
            END IF
         END DO
         J = 1
         DO
            IF( J.EQ.4 ) EXIT
            IF( DATASET(3,I,J).LT.DATASET(3,I,J+1) ) EXIT
            J = J+1
         END DO
         PMINS(I) = NMAT(I)/NPROC(J)
      END DO
      PMIN = (1*PMINS(1)+1*PMINS(2)+2*PMINS(3)+8*PMINS(4)+8*PMINS(5)
     $     +8*PMINS(6)+4*PMINS(7)) / (1+1+2+8+8+8+4)
      WRITE(*,*), 'Suggestions on PILAENVX(ISPEC=23):'
      WRITE(*,*), '    PILAENVX = ICEIL(N1, ICEIL(NTHRESH, N2)*N2)'
      WRITE(*,*), 'with', MAX(NB, PMIN/2), '< NTHRESH <', PMIN
      CURRENT = PILAENVX(ICTXT, 23, '', '', 1000000000, 1, 0, 0)
      CURRENT = NINT(1.0D+9/DBLE(CURRENT))
      IF( MAX(NB, PMIN/2).LE.CURRENT .AND. CURRENT.LE.PMIN ) THEN
         WRITE(*,*), 'The current value ( NTHRESH =', CURRENT,
     $        ') is acceptable.'
      ELSE
         WRITE(*,*), 'The current value ( NTHRESH =', CURRENT,
     $        ') needs to be replaced.'
      END IF
      WRITE(*,*)
*
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
*
      END
*
      INTEGER FUNCTION N2IND(N)
      INTEGER N
*
      IF( N.EQ.96 ) THEN
         N2IND = 1
      ELSE IF( N.EQ.192 ) THEN
         N2IND = 2
      ELSE IF( N.EQ.384 ) THEN
         N2IND = 3
      ELSE IF( N.EQ.768 ) THEN
         N2IND = 4
      ELSE IF( N.EQ.1536 ) THEN
         N2IND = 5
      ELSE IF( N.EQ.3072 ) THEN
         N2IND = 6
      ELSE
         N2IND = 7
      END IF
*
      RETURN
      END
*
      INTEGER FUNCTION NP2IND(NP)
      INTEGER NP
*
      IF( NP.EQ.1 ) THEN
         NP2IND = 1
      ELSE IF( NP.EQ.2 ) THEN
         NP2IND = 2
      ELSE IF( NP.EQ.4 ) THEN
         NP2IND = 3
      ELSE
         NP2IND = 4
      END IF
*
      RETURN
      END
