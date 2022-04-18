C   The function REVOLV coded below is meant to be used as a         * 
C   "controller" for running a time-dependent applications program   *
C   in the reverse mode with checkpointing described in the paper    *
C   "Achieving logarithmic Growth in temporal and spatial complexity *
C   in reverse automatic differentiation", Optimization Methods and  *
C   Software,  Vol.1 pp. 35-54.                                      *
C   A postscript source of that paper can be found in the ftp sites  *
C        info.mcs.anl.gov and nbtf02.math.tu-dresden.de.             *
C   Apart from REVOLV this file contains five auxiliary routines     * 
C   NUMFRW, EXPNS, MAXRGE, and ADJUST.                               *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   To utilize REVOLV the user must have procedures for              *
C     - Advancing the state of the modeled system to a certain time. *
C     - Saving the current state onto a stack of snapshots.          *
C     - Restoring the the most recently saved snapshot and           *
C       restarting the forward simulation from there.                *
C     - Initializing the adjoints at the end of forward sweep.       *
C     - Performing one combined forward and adjoint step.            * 
C   Through an encoding of its return value REVOLV asks the calling  *
C   program to perform one of these 'actions', which we will         *
C   refer to as                                                      *
C                                                                    *
C       'advance', 'takeshot', 'restore', 'firsturn' and 'youturn'  .*
C   There are two other return values, namely                        *
C       'terminate'   and     'error'                                *
C   which indicate a regular or faulty termination of the calls      *
C   to REVOLV.                                                       *
C                                                                    *
C   The action 'firsturn' includes a 'youturn', in that it requires  *
C     -advancing through the last time-step with recording           *
C      of intermediates                                              *
C     -initializing the adjoint values (possibly after               *
C      performing some IO)                                           *
C     -reversing the last time step using the record just written    *
C   The action 'firsturn' is obtained when the difference FINE-CAPO  *
C   has been reduced to 1 for the first time.                        *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   The calling sequence is                                          *
C                                                                    *
C               REVOLV(CHECK,CAPO,FINE,SNAPS,INFO)                   *
C                                                                    *
C   with the return value being one of the actions to be taken. The  *
C   calling parameters are all integers with the following meaning   *
C                                                                    *
C         CHECK     number of checkpoint being written or retrieved  *
C         CAPO      beginning of subrange currently being processed  *
C         FINE      end of subrange currently being processed        *
C         SNAPS     upper bound on number of checkpoints taken       *
C         INFO      determines how much information will be printed  *
C                   and contains information about an error occured  *
C                                                                    *
C   Since REVOLV involves only a few integer operations its run-time *
C   is truly negligible within any nontrivial application.           *
C                                                                    *
C   The parameter SNAPS is selected by the user (possibly with the   *
C   help of the routines EXPNS and  ADJUST described below ) and     *
C   remains unchanged throughout.                                    *
C                                                                    *
C   The pair (CAPO,FINE) always represents the initial and final     *
C   state of the subsequence of time steps currently being traversed *
C   backwards.                                                       *
C                                                                    *
C   The conditions                                                   *
C                    CHECK >= -1      and     CAPO <= FINE           *
C   are necessary and sufficient for a regular response of REVOLV.   *
C   If either condition is violated the value 'error' is returned.   *
C                                                                    *
C   The first call to REVOLV must be with CHECK=-1 so that           * 
C   appropriate initializations can be performed internally.         *
C                                                                    *
C   When CHECK =-1 and CAPO = FINE  then 'terminate' is returned as  *
C   action value. This combination necessarily arises after a        *
C   sufficiently large number of calls to REVOLV, which depends only * 
C   on the initial difference FINE-CAPO.                             *
C                                                                    *
C   The last parameter INFO determines how much information about    *
C   the actions performed will be printed. When INFO =0 no           *
C   information is sent to standard output. When INFO > 0 REVOLV     *
C   produces an output that contains a prediction of the number of   *    
C   forward steps and of the factor by which the execution will slow *    
C   down. When an error occurs, the return value of INFO contains    *
C   information about the reason:                                    *
C     INFO = 10: number of checkpoints stored exceeds CHECKUP,       *
C                increase constant CHECKUP and recompile             *
C     INFO = 11: number of checkpoints stored exceeds SNAPS, ensure  * 
C                SNAPS greater than 0 and increase initial FINE      *
C     INFO = 12: error occurs in NUMFORW                             *
C     INFO = 13: enhancement of FINE, SNAPS checkpoints stored,      *
C                SNAPS must be increased                             *
C     INFO = 14: number of SNAPS exceeds CHECKUP, increase constant  *
C                CHECKUP and recompile                               *
C     INFO = 15: number of REPS exceeds REPSUP, increase constant    *
C                REPSUP and recompile                                *
C                                                                    *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   Some further explanations and motivations:                       *
C                                                                    *
C   There is an implicit bound on CHECK through the dimensioning of  *
C   the integer array CH[CHEKUP] with CHECKUP = 64 being the default.*
C   If anybody wants to have that even larger he must change the     *
C   source. Also for the variable REPS an upper bound REPSUP is      *
C   defined. The default value equals 64. If during a call to        *
C   TREEVERSE a (CHECKUP+1)-st checkpoint would normally be called   * 
C   for then control is returned after an appropriate error message. * 
C   When the calculated REPS exceeds REPSUP also an error message    *
C   occurs.                                                          *
C   During the forward sweep the user is free to change the last     *
C   three parameters from call to call, except that FINE may never   *
C   be less than the current value of CAPO. This may be useful when  *
C   the total number of time STEPS to be taken is not a priori       *
C   known. The choice FINE=CAPO+1 initiates the reverse sweep, which * 
C   happens automatically if is left constant as CAPO is eventually  * 
C   moved up to it. Once the first reverse or restore action has     *
C   been taken only the last two parameters should be changed.       *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   The necessary number of forward steps without recording is       *
C   calculated by the function                                       *
C                                                                    *
C                      NUMFRW(STEPS,SNAPS)                           *
C                                                                    *
C   STEPS denotes the total number of time steps, i.e. FINE-CAPO     *
C   during the first call of REVOLV. When SNAPS is less than 1 an    * 
C   error message will be given and -1 is returned as value.         *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   To choose an appropriated value of SNAPS the function            *
C                                                                    *
C                      EXPNS(STEPS,SNAPS)                            *
C                                                                    *
C   estimates the run-time factor  incurred by REVOLV for a          *
C   particular value of SNAPS. The ratio NUMFRW(STEPS,SNAPS)/STEPS   *
C   is returned. This ratio corresponds to the run-time factor of    *
C   the execution relative to the run-time of one forward time step. *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   The auxiliary function                                           *
C                                                                    *
C                      MAXRGE(SNAPS,REPS)                            *
C                                                                    *
C   returns the integer (SNAPS+REPS)!/(SNAPS!REPS!) provided         *
C   SNAPS >=0, REPS >= 0. Otherwise there will be appropriate error  *
C   messages and the value -1 will be returned. If the binomial      *
C   expression is not representable as a  signed 4 byte integer,     *
C   greater than 2^31-1, this maximal value is returned and a        *
C   warning message printed.                                         *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C   Furthermore, the function                                        *
C                                                                    *
C                      ADJUST(STEPS)                                 *
C                                                                    *
C   is provided. It can be used to determine a value of SNAPS so     *
C   that the increase in spatial complexity equals approximately the *
C   increase in temporal complexity. For that ADJUST computes a      *
C   return value satisfying SNAPS ~= log_4 (STEPS) because of the    *
C   theory developed in the paper mentioned above.                   *
C                                                                    *
C--------------------------------------------------------------------*


      INTEGER FUNCTION NUMFRW(STEPS,SNAPS)


C     .. Parameters ..
      INTEGER CHEKUP,REPSUP
      PARAMETER (CHEKUP=64,REPSUP=64)
C     ..
C     .. Scalar Arguments ..
      INTEGER SNAPS,STEPS
C     ..
C     .. Local Scalars ..
      INTEGER RANGE,REPS
C     ..
      IF (SNAPS.LT.1) THEN
          NUMFRW = -1

      ELSE
          REPS = 0
          RANGE = 1
   10     CONTINUE
          IF (RANGE.LT.STEPS) THEN
              REPS = REPS + 1
              RANGE = RANGE* (REPS+SNAPS)/REPS
              GO TO 10

          END IF

          IF (REPS.GT.REPSUP) THEN
              WRITE (*,FMT=*) ' number of reps exceeds REPSUP, '
              WRITE (*,FMT=*) ' increase constant REPSUP and recompile'
              NUMFRW = -1
              RETURN

          END IF

          IF (SNAPS.LE.CHEKUP) THEN
              NUMFRW = REPS*STEPS - RANGE*REPS/(SNAPS+1)
          END IF

      END IF

      END

C--------------------------------------------------------------------*

      DOUBLE PRECISION FUNCTION EXPNS(STEPS,SNAPS)


C     .. Scalar Arguments ..
      INTEGER SNAPS,STEPS
C     ..
C     .. External Functions ..
      INTEGER NUMFRW
      EXTERNAL NUMFRW
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
      IF (SNAPS.LT.1) THEN
          WRITE (*,FMT=*) 'error occurs in EXPNS: SNAPS < 1'
          EXPNS = -1

      ELSE
          IF (SNAPS.LT.1) THEN
              WRITE (*,FMT=*) 'error occurs in EXPNS: SNAPS < 1'
              EXPNS = -1

          ELSE
              EXPNS = DBLE(NUMFRW(STEPS,SNAPS))
              IF (EXPNS.NE.-1) THEN
                  EXPNS = EXPNS/STEPS
              END IF

          END IF

      END IF

      END

C--------------------------------------------------------------------*

      INTEGER FUNCTION MAXRGE(SS,TT)

C     .. Scalar Arguments ..
      INTEGER SS,TT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RES
      INTEGER I
C     ..
      RES = 1.
      IF (TT.LT.0 .OR. SS.LT.0) THEN
          WRITE (*,FMT=*) 'error in MAXRGE: negative parameter '
          MAXRGE = -1
          GO TO 30

      END IF

      DO 10 I = 1,TT
          RES = RES* (SS+I)
          RES = RES/I
          IF (RES.GE.2.0**31) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IF (RES.LT.2.0**31-2) THEN
          MAXRGE = RES

      ELSE
          MAXRGE = 2.0**31 - 3
          WRITE (*,FMT=*)
     +      'warning from  MAXRGE: returned maximal integer'
          WRITE (*,FMT=*) MAXRGE
      END IF

   30 CONTINUE
      RETURN

      END

C--------------------------------------------------------------------*

      INTEGER FUNCTION ADJUST(STEPS)


C     .. Scalar Arguments ..
      INTEGER STEPS
C     ..
C     .. Local Scalars ..
      INTEGER REPS,S,SNAPS
C     ..
C     .. External Functions ..
      INTEGER MAXRGE
      EXTERNAL MAXRGE
C     ..
      SNAPS = 1
      REPS = 1
      S = 0
   10 IF (MAXRGE(SNAPS+S,REPS+S).GT.STEPS) THEN
          S = S - 1
          GO TO 10

      END IF

   20 IF (MAXRGE(SNAPS+S,REPS+S).LT.STEPS) THEN
          S = S + 1
          GO TO 20

      END IF

      SNAPS = SNAPS + S
      REPS = REPS + S
      S = -1
   30 IF (MAXRGE(SNAPS,REPS).GE.STEPS) THEN
          IF (SNAPS.GT.REPS) THEN
              SNAPS = SNAPS - 1
              S = 0

          ELSE
              REPS = REPS - 1
              S = 1
          END IF

          GO TO 30

      END IF

      IF (S.EQ.0) SNAPS = SNAPS + 1
      IF (S.EQ.1) REPS = REPS + 1
      ADJUST = SNAPS
      RETURN

      END

C--------------------------------------------------------------------*

      INTEGER FUNCTION REVOLV(CHECK,CAPO,FINE,SNAPS,INFO)

C     (CAPO ,FINE) is the time range currently under consideration
C     .. Parameters ..
      INTEGER CHEKUP,REPSUP
      PARAMETER (CHEKUP=64,REPSUP=64)
      INTEGER TAKSHT,ADVAN,FSTURN,YUTURN
      PARAMETER (TAKSHT=1,ADVAN=2,FSTURN=3,YUTURN=4)
      INTEGER RESTRE,TRMATE,ERROR
      PARAMETER (RESTRE=5,TRMATE=6,ERROR=7)
C     ..
C     .. Scalar Arguments ..
      INTEGER CAPO,CHECK,FINE,INFO,REPS,SNAPS
C     ..
C     .. Scalars in Common ..
      INTEGER NUMADV,NUMCMD,NUMTKS,OLDFIN,OLDSNP,TURN
C     ..
C     .. Arrays in Common ..
      INTEGER CH(0:CHEKUP)
C     ..
C     .. Local Scalars ..
      INTEGER DS,I,NUM,OLDCPO,RANGE
      REAL BINO1, BINO2, BINO3, BINO4, BINO5
C     ..
C     .. External Functions ..
      INTEGER NUMFRW
      EXTERNAL NUMFRW
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
C     .. Common blocks ..
      COMMON /INFNUM/NUMADV,NUMTKS,NUMCMD
      COMMON /SNREP/CH,TURN,OLDSNP,OLDFIN
C     ..
      IF (INFO.GT.3) INFO = 3
      NUMCMD = NUMCMD + 1
      IF ((CHECK.LT.-1) .AND. (CAPO.GT.FINE)) THEN
          REVOLV = ERROR
          GO TO 30

      END IF

      IF ((CHECK.EQ.-1) .AND. (CAPO.LT.FINE)) THEN
          IF (CHECK.EQ.-1) TURN = 0
          DO 10 I = 0,CHEKUP
              CH(I) = 0
   10     CONTINUE
          CH(0) = CAPO - 1
      END IF
C
      IF ((FINE-CAPO).EQ.0) THEN
          IF ((CHECK.EQ. (-1)) .OR. (CAPO.EQ.CH(0))) THEN
              CHECK = CHECK - 1
              IF (INFO.GT.0) THEN
                  WRITE (*,FMT=*)
                  WRITE (*,FMT=9000) NUMADV
                  WRITE (*,FMT=9010) NUMTKS
                  WRITE (*,FMT=9020) NUMCMD
              END IF

              REVOLV = TRMATE
              GO TO 30

          ELSE
              CAPO = CH(CHECK)
              OLDFIN = FINE
              REVOLV = RESTRE
              GO TO 30

          END IF

      END IF

      IF ((FINE-CAPO).EQ.1) THEN
          FINE = FINE - 1
          IF ((CHECK.GE.0) .AND. (CH(CHECK).EQ.CAPO)) THEN
              CHECK = CHECK - 1
          END IF

          IF (TURN.EQ.0) THEN
              OLDFIN = FINE
              REVOLV = FSTURN
              TURN = 1

          ELSE
              OLDFIN = FINE
              REVOLV = YUTURN
          END IF

          GO TO 30

      END IF

      IF ((CHECK.EQ. (-1)) .OR. (CH(CHECK).NE.CAPO)) THEN
          CHECK = CHECK + 1
          IF (CHECK.GE.CHEKUP) THEN
              INFO = 10
              REVOLV = ERROR
              GO TO 30

          ELSE
              IF (CHECK+1.GT.SNAPS) THEN
                  INFO = 11
                  REVOLV = ERROR
                  GO TO 30

              END IF

              CH(CHECK) = CAPO
              IF (CHECK.EQ.0) THEN
                  NUMADV = 0
                  NUMTKS = 0
                  NUMCMD = 1
                  OLDSNP = SNAPS
                  IF (SNAPS.GT.CHEKUP) THEN
                      INFO = 14
                      REVOLV = ERROR
                      GO TO 30

                  END IF

                  IF (INFO.GT.0) THEN
                      NUM = NUMFRW(FINE-CAPO,SNAPS)
                      IF (NUM.EQ.-1) THEN
                          INFO = 12
                          REVOLV = ERROR
                          GO TO 30

                      END IF

                      WRITE (*,FMT=9030) NUM
                      WRITE (*,FMT=9040) DBLE(NUM)/ (FINE-CAPO)
                      WRITE (*,FMT=*)
                  END IF

              END IF

              NUMTKS = NUMTKS + 1
              OLDFIN = FINE
              REVOLV = TAKSHT
              GO TO 30

          END IF

      ELSE
          IF ((OLDFIN.LT.FINE) .AND. (SNAPS.EQ.CHECK+1)) THEN
              INFO = 13
              REVOLV = ERROR
              GO TO 30

          END IF

          OLDCPO = CAPO
          DS = SNAPS - CHECK
          IF (DS.LT.1) THEN
              INFO = 11
              REVOLV = ERROR
              GO TO 30

          END IF

          REPS = 0
          RANGE = 1
   20     CONTINUE
          IF (RANGE.LT.FINE-CAPO) THEN
              REPS = REPS + 1
              RANGE = RANGE* (REPS+DS)/REPS
              GO TO 20

          END IF

          IF (REPS.GT.REPSUP) THEN
              INFO = 15
              REVOLV = ERROR
              GO TO 30

          END IF

          IF (SNAPS.NE.OLDSNP) THEN
              IF (SNAPS.GT.CHEKUP) THEN
                  INFO = 14
                  REVOLV = ERROR
                  GO TO 30

              END IF

          OLDSNP = SNAPS
          END IF

          BINO1 = RANGE*REPS/(DS+REPS)
          IF (DS .GT. 1) THEN
            BINO2 = BINO1*DS/(DS+REPS-1)
          ELSE
            BINO2 = 1
          ENDIF
          IF (DS .EQ. 1) THEN
            BINO3 = 0
          ELSE
            IF (DS .GT. 2) THEN
              BINO3 = BINO2*(DS-1)/(DS+REPS-2)
            ELSE
              BINO3 = 1
            ENDIF
          ENDIF
          BINO4 = BINO2*(REPS-1)/DS
          IF (DS .LT. 3) THEN
            BINO5 = 0
          ELSE
            IF (DS .GT. 3) THEN
              BINO5 = BINO3*(DS-1)/REPS
            ELSE
              BINO5 = 1
            ENDIF
          ENDIF

          IF (FINE-CAPO.LE.BINO1+BINO3) THEN
              CAPO = CAPO + BINO4
          ELSE
              IF (FINE-CAPO.GE.RANGE-BINO5) THEN
                  CAPO = CAPO + BINO1
              ELSE
                  CAPO = FINE - BINO2 - BINO3
              END IF
          END IF

          IF (CAPO.EQ.OLDCPO) CAPO = OLDCPO + 1
          NUMADV = NUMADV + CAPO - OLDCPO
          REVOLV = ADVAN
      END IF

   30 CONTINUE
      RETURN

 9000 FORMAT (' advances:',I8)
 9010 FORMAT (' takeshots:',I7)
 9020 FORMAT (' commands:',I8)
 9030 FORMAT (' prediction of needed forward steps: ',I7,' => ')
 9040 FORMAT (' slowdown factor: ',F8.4)
      END
