      PROGRAM MAIN
C

C     .. Parameters ..
      REAL TAKSHT,ADVAN,FSTURN,YUTURN
      PARAMETER (TAKSHT=1,ADVAN=2,FSTURN=3,YUTURN=4)
      REAL RESTRE,TRMATE,ERROR
      PARAMETER (RESTRE=5,TRMATE=6,ERROR=7)
C     ..
C     .. Local Scalars ..
      INTEGER CAPO,CHECK,FINE,INFO,SNAPS,STEPS,WHATDO
C     ..
C     .. External Functions ..
      INTEGER REVOLV
      EXTERNAL REVOLV
C     ..
      WRITE (*,FMT=*) 'ENTER:   STEPS, SNAPS, INFO'
      READ (*,FMT=*) STEPS,SNAPS,INFO
      WRITE (*,FMT=*) STEPS,SNAPS,INFO
      CAPO = 0
      FINE = STEPS + CAPO
      CHECK = -1
   10 CONTINUE
      WHATDO = REVOLV(CHECK,CAPO,FINE,SNAPS,INFO)
      IF ((WHATDO.EQ.TAKSHT) .AND. (INFO.GT.1)) THEN
          WRITE (*,FMT=9000) CAPO
      END IF

      IF ((WHATDO.EQ.ADVAN) .AND. (INFO.GT.2)) THEN
          WRITE (*,FMT=9010) CAPO
      END IF

      IF ((WHATDO.EQ.FSTURN) .AND. (INFO.GT.2)) THEN
          WRITE (*,FMT=9020) CAPO
      END IF

      IF ((WHATDO.EQ.YUTURN) .AND. (INFO.GT.2)) THEN
          WRITE (*,FMT=9030) CAPO
      END IF

      IF ((WHATDO.EQ.RESTRE) .AND. (INFO.GT.2)) THEN
          WRITE (*,FMT=9040) CAPO
      END IF

      IF (WHATDO.EQ.ERROR) THEN
          WRITE (*,FMT=*) ' irregular termination of treeverse'
          IF (INFO.EQ.10) THEN
              WRITE (*,FMT=*)
     +          ' number of checkpoints stored exceeds CHEKUP,'
              WRITE (*,FMT=*) ' increase constant CHEKUP and recompile'
          END IF

          IF (INFO.EQ.11) THEN
              WRITE (*,FMT=*) ' number of checkpoints stored = ',
     +          CHECK + 1,' exceeds SNAPS,'
              WRITE (*,FMT=*) ' ensure SNAPS > 0 and ',
     +          'increase initial FINE'
          END IF

          IF (INFO.EQ.12) WRITE (*,FMT=*) ' error occurs in NUMFRW'
          IF (INFO.EQ.13) THEN
              WRITE (*,FMT=*) ' enhancement of FINE, SNAPS = ',
     +          CHECK + 1,'checkpoints stored, increase SNAPS'
          END IF
          IF (INFO.EQ.14) THEN
              WRITE (*,FMT=*) ' number of SNAPS = ',SNAPS,
     +          ' exceeds CHEKUP,'
              WRITE (*,FMT=*) ' increase constant CHEKUP and recompile'
          END IF
          IF (INFO.EQ.15) THEN
              WRITE (*,FMT=*) ' number of reps exceeds REPSUP, '
              WRITE (*,FMT=*) ' increase constant REPSUP and recompile'
          END IF

      END IF

      IF ((WHATDO.EQ.TRMATE) .OR. (WHATDO.EQ.ERROR)) THEN
          GO TO 20

      ELSE
          GO TO 10

      END IF

   20 CONTINUE

 9000 FORMAT (' takeshot at',I6)
 9010 FORMAT (' advance to',I7)
 9020 FORMAT (' firsturn at',I6)
 9030 FORMAT (' youturn at',I7)
 9040 FORMAT (' restore at',I7)
      END
