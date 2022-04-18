C     ALGORITHM 595, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 1,
C     MAR., 1983, P. 131-138.
C SAMPLE DRIVER PROGRAM FOR HC.                                         MAN   10
C                                                                       MAN   20
C THIS PROGRAM FINDS ONE OR MORE HAMILTONIAN CIRCUITS IN A              MAN   30
C DIRECTED GRAPH OF  N  VERTICES AND  M  ARCS.                          MAN   40
C                                                                       MAN   50
C COMPLETE DETAILS OF THE PARAMETERS MAY BE FOUND IN THE                MAN   60
C DOCUMENTATION OF SUBROUTINE HC.                                       MAN   70
C                                                                       MAN   80
C THE INPUT UNIT NUMBER IS ASSUMED TO BE  5 .                           MAN   90
C THE OUTPUT UNIT NUMBER IS ASSUMED TO BE  6 .                          MAN  100
C THE ARRAYS ARE CURRENTLY DIMENSIONED TO ALLOW PROBLEMS FOR            MAN  110
C WHICH  N .LE. 250  AND  M .LE. 2000 .                                 MAN  120
C                                                                       MAN  130
C THE PROGRAM MAY BE TESTED ON THE FOLLOWING DATA                       MAN  140
C                                                                       MAN  150
C N  =     6                                                            MAN  160
C PR =     0    3    5    8   11   13   16                              MAN  170
C AR =     3    5    6    6    3    6    1    4    5    1               MAN  180
C          2    2    3    4    2    5                                   MAN  190
C                                                                       MAN  200
      INTEGER PR(251), PC(251), AR(2000), AC(2000), S(250), VR(250),    MAN  210
     * VC(250), P(250), SUBR(250), RBUS(250), TOR(250)                  MAN  220
      NIN = 5                                                           MAN  230
      NOUT = 6                                                          MAN  240
C INPUT DATA                                                            MAN  250
      READ (NIN,99999) N                                                MAN  260
      NP1 = N + 1                                                       MAN  270
      READ (NIN,99999) (PR(J),J=1,NP1)                                  MAN  280
      M = PR(NP1)                                                       MAN  290
      READ (NIN,99999) (AR(J),J=1,M)                                    MAN  300
      WRITE (NOUT,99998)                                                MAN  310
      WRITE (NOUT,99997) N                                              MAN  320
      WRITE (NOUT,99996) (PR(J),J=1,NP1)                                MAN  330
      WRITE (NOUT,99995) (AR(J),J=1,M)                                  MAN  340
      WRITE (NOUT,99994)                                                MAN  350
C CALL HC AS EXACT PROCEDURE TO FIND (AND PRINT) A SINGLE               MAN  360
C HAMILTONIAN CIRCUIT, IF ONE EXISTS.                                   MAN  370
      NC = 1                                                            MAN  380
      NB = -1                                                           MAN  390
      CALL HC(N, PR, AR, NOUT, NC, NB, S, N+1, PR(N+1), PC, AC, VR, VC, MAN  400
     * P, SUBR, RBUS, TOR)                                              MAN  410
      WRITE (NOUT,99993) NC, NB                                         MAN  420
C CALL HC AS EXACT PROCEDURE TO FIND (AND PRINT) ALL THE                MAN  430
C HAMILTONIAN CIRCUITS.                                                 MAN  440
      NC = -1                                                           MAN  450
      NB = -1                                                           MAN  460
      CALL HC(N, PR, AR, NOUT, NC, NB, S, N+1, PR(N+1), PC, AC, VR, VC, MAN  470
     * P, SUBR, RBUS, TOR)                                              MAN  480
      WRITE (NOUT,99993) NC, NB                                         MAN  490
C CALL HC AS HEURISTIC PROCEDURE TO FIND (AND PRINT) A                  MAN  500
C SINGLE HAMILTONIAN CIRCUIT, IF ONE EXISTS, WITHOUT                    MAN  510
C PERFORMING MORE THAN  2  BACKTRACKINGS.                               MAN  520
      NC = 1                                                            MAN  530
      NB = 2                                                            MAN  540
      CALL HC(N, PR, AR, NOUT, NC, NB, S, N+1, PR(N+1), PC, AC, VR, VC, MAN  550
     * P, SUBR, RBUS, TOR)                                              MAN  560
      WRITE (NOUT,99993) NC, NB                                         MAN  570
C CALL HC AS HEURISTIC PROCEDURE TO FIND (AND PRINT) A                  MAN  580
C SINGLE HAMILTONIAN CIRCUIT, IF ONE EXISTS, WITHOUT                    MAN  590
C PERFORMING MORE THAN  4  BACKTRACKINGS.                               MAN  600
      NC = 1                                                            MAN  610
      NB = 4                                                            MAN  620
      CALL HC(N, PR, AR, NOUT, NC, NB, S, N+1, PR(N+1), PC, AC, VR, VC, MAN  630
     * P, SUBR, RBUS, TOR)                                              MAN  640
      WRITE (NOUT,99993) NC, NB                                         MAN  650
C CALL HC AS HEURISTIC PROCEDURE TO FIND (WITHOUT PRINTING)             MAN  660
C AT MOST  2  HAMILTONIAN CIRCUITS, WITHOUT PERFORMING MORE             MAN  670
C THAN  5  BACKTRACKINGS.                                               MAN  680
      NC = 2                                                            MAN  690
      NB = 5                                                            MAN  700
      CALL HC(N, PR, AR, -1, NC, NB, S, N+1, PR(N+1), PC, AC, VR, VC,   MAN  710
     * P, SUBR, RBUS, TOR)                                              MAN  720
      WRITE (NOUT,99993) NC, NB                                         MAN  730
      IF (NC.EQ.0) STOP                                                 MAN  740
C PRINT THE LAST HAMILTONIAN CIRCUIT FOUND                              MAN  750
      WRITE (NOUT,99992) (S(J),J=1,N)                                   MAN  760
      STOP                                                              MAN  770
99999 FORMAT (10I5)                                                     MAN  780
99998 FORMAT (1H1////////)                                              MAN  790
99997 FORMAT (6H N  = , I5)                                             MAN  800
99996 FORMAT (6H PR = , 25I5)                                           MAN  810
99995 FORMAT (6H AR = , 25I5)                                           MAN  820
99994 FORMAT (//////)                                                   MAN  830
99993 FORMAT (/I5, 10H CIRCUITS , I5, 14H BACKTRACKINGS///)             MAN  840
99992 FORMAT (4X, 4HS = , 25I5)                                         MAN  850
      END                                                               MAN  860
      SUBROUTINE HC(N, PR, AR, KW, NC, NB, S, NP1, M, PC, AC, VR, VC,   HC    10
     * P, SUBR, RBUS, TOR)
C
C SUBROUTINE TO FIND ONE OR MORE HAMILTONIAN CIRCUITS IN A
C DIRECTED GRAPH OF  N  VERTICES ( N .GT. 1 ) REPRESENTED
C BY THE INTEGERS  1, 2, ..., N  AND  M  ARCS.
C
C HC IS BASED ON AN ENUMERATIVE ALGORITHM AND CAN BE USED
C EITHER AS AN EXACT PROCEDURE OR AS A HEURISTIC PROCEDURE
C (BY LIMITING THE NUMBER OF BACKTRACKINGS ALLOWED).
C
C ENTRANCE TO HC IS ACHIEVED BY USING THE STATEMENT
C     CALL HC(N,PR,AR,KW,NC,NB,S,N+1,PR(N+1),PC,AC,VR,VC,
C    *        P,SUBR,RBUS,TOR)
C
C THE VALUES OF THE FIRST SIX PARAMETERS MUST BE DEFINED
C BY THE USER PRIOR TO CALLING HC. HC NEEDS  2  ARRAYS ( PR
C AND  PC ) OF LENGTH  N + 1 ,  2  ARRAYS ( AR  AND  AC )
C OF LENGTH  M  AND  7  ARRAYS ( S ,  VR ,  VC ,  SUBR ,
C RBUS  AND  TOR ) OF LENGTH  N . THESE ARRAYS MUST BE
C DIMENSIONED BY THE USER IN THE CALLING PROGRAM.
C
C HC CALLS  5  SUBROUTINES: PATH, FUPD, BUPD, IUPD, RARC.
C THESE SUBROUTINES ARE COMPLETELY LOCAL, I.E. THE INFORMA-
C TION THEY NEED IS PASSED THROUGH THE PARAMETER LIST.
C THE WHOLE PACKAGE IS COMPLETELY SELF CONTAINED AND COMMU-
C NICATION TO IT IS ACHIEVED SOLELY THROUGH THE PARAMETER
C LIST OF HC. NO MACHINE DEPENDENT CONSTANTS ARE USED.
C THE PACKAGE IS WRITTEN IN AMERICAN NATIONAL STANDARD
C FORTRAN AND IS ACCEPTED BY THE  FTN(EL=A)  COMPILER OF THE
C CDC CYBER 76 (OPTION  EL=A  CHECKS PROGRAM AND SUBROUTINES
C FOR ADHERENCE TO ANSI) AS WELL AS BY THE PFORT VERIFIER.
C THE PACKAGE HAS BEEN TESTED ON A  CDC CYBER 76 , ON A  CDC
C 6600 , ON AN  IBM 370/158  AND ON A  DIGITAL VAX 11/780 .
C
C MEANING OF THE INPUT PARAMETERS:
C N     = NUMBER OF VERTICES.
C PR(I) = SUM OF THE OUT-DEGREES OF VERTICES  1, ..., I-1
C         ( PR(1) = 0 ,  PR(N+1) = M ).
C AR    = ADJACENCY LIST. THE ELEMENTS FROM  AR(PR(I)+1)  TO
C         AR(PR(I+1))  ARE A RECORD CONTAINING,IN ANY ORDER,
C         ALL THE VERTICES  J  SUCH THAT ARC  (I,J)  EXISTS.
C         THE GRAPH SHOULD NOT CONTAIN ARCS STARTING AND
C         ENDING AT THE SAME VERTEX.
C KW    = UNIT TAPE ON WHICH TO WRITE THE HAMILTONIAN CIR-
C         CUITS FOUND, ACCORDING TO FORMAT  20I5 .  KW = - 1
C         IF NO WRITING IS DESIRED. THE CIRCUITS ARE WRITTEN
C         AS ORDERED SEQUENCES OF  N  VERTICES.
C
C MEANING OF THE INPUT-OUTPUT PARAMETERS:
C NC(INPUT)  = UPPER BOUND ON THE NUMBER OF HAMILTONIAN
C              CIRCUITS TO BE FOUND ( NC = - 1 IF ALL THE
C              HAMILTONIAN CIRCUITS ARE TO BE FOUND).
C NC(OUTPUT) = NUMBER OF HAMILTONIAN CIRCUITS FOUND.
C NB(INPUT)  = - 1  IF HC MUST BE EXECUTED AS AN EXACT
C              PROCEDURE.
C            = UPPER BOUND ON THE NUMBER OF BACKTRACKINGS IF
C              HC MUST BE EXECUTED AS A HEURISTIC PROCEDURE.
C NB(OUTPUT) = NUMBER OF BACKTRACKINGS PERFORMED. WHEN HC
C              HAS BEEN EXECUTED AS A HEURISTIC PROCEDURE,
C              IF  NB(OUTPUT) .LT. NB(INPUT)  THEN THE
C              RESULT OBTAINED IS EXACT.
C
C MEANING OF THE OUTPUT PARAMETER:
C S(I) = I-TH  VERTEX IN THE LAST HAMILTONIAN CIRCUIT FOUND.
C
C ON RETURN OF HC  N, PR  AND  KW  ARE UNCHANGED, WHILE IN
C AR  THE ORDER OF THE ELEMENTS WITHIN EACH RECORD MAY BE
C ALTERED.
C
C MEANING OF THE WORK ARRAYS:
C PC(I)   = SUM OF THE IN-DEGREES OF VERTICES  1, ..., I-1
C           ( PC(1) = 0 ).
C AC      = ADJACENCY LIST (BACKWARD). THE ELEMENTS FROM
C           AC(PC(I)+1)  TO  AC(PC(I+1))  CONTAIN, IN ANY
C           ORDER, ALL THE VERTICES  J  SUCH THAT ARC  (J,I)
C           EXISTS.
C WHEN AN ARC IS REMOVED FROM THE GRAPH AT THE  K-TH  LEVEL
C OF THE BRANCH-DECISION TREE, THE CORRESPONDING ELEMENTS
C AR(Q)  AND  AC(T)  ARE SET TO  - (K*(N+1) + AR(Q))  AND
C TO  - (K*(N+1) + AC(T)) , RESPECTIVELY.
C VR(I)   = CURRENT OUT-DEGREE OF VERTEX  I .
C VC(I)   = CURRENT IN-DEGREE OF VERTEX  I .
C SUBR(I) = - (K*(N+1) + J)  IF ARC  (I,J)  WAS IMPLIED AT
C           THE  K-TH  LEVEL OF THE BRANCH-DECISION TREE.
C         = 0  OTHERWISE.
C RBUS(I) = - J  IF ARC  (J,I)  IS CURRENTLY IMPLIED.
C         = 0  OTHERWISE.
C TOR(K)  = Q*(M+1) + T  IF THE ARC GOING FROM  S(K)  TO THE
C           ROOT, CORRESPONDING TO  AR(Q)  AND TO  AC(T),
C           WAS REMOVED FROM THE GRAPH AT THE  K-TH  LEVEL
C           OF THE BRANCH-DECISION TREE.
C         = 0  OTHERWISE.
C P(I)    = POINTER FOR THE FORWARD STEP. THE NEXT ARC
C           STARTING FROM  I  TO BE CONSIDERED IN THE
C           BRANCH-DECISION TREE IS  (I,AR(PR(I)+P(I)).
C
C MEANING OF THE MAIN WORK SIMPLE VARIABLES:
C JR  = ROOT. THE HAMILTONIAN CIRCUITS ARE DETERMINED AS
C       PATHS STARTING AND ENDING AT  JR .
C K   = CURRENT LEVEL OF THE BRANCH-DECISION TREE.
C M   = NUMBER OF ARCS.
C MP1 = M + 1 (USED FOR PACKING  TOR ).
C NP1 = N + 1 (USED FOR PACKING  AR ,  AC  AND  SUBR )
C
      INTEGER PR(NP1), PC(NP1), AR(M), AC(M), S(N), VR(N), VC(N), P(N),
     * SUBR(N), RBUS(N), TOR(N)
C
C S T E P   0   (INITIALIZE).
C
      NCO = NC
      NC = 0
      NBO = NB
      NB = 0
      DO 10 I=1,N
        VC(I) = 0
        SUBR(I) = 0
        RBUS(I) = 0
        P(I) = 1
   10 CONTINUE
      DO 30 I=1,N
        J1 = PR(I) + 1
        J2 = PR(I+1)
        VR(I) = J2 - J1 + 1
        IF (VR(I).EQ.0) RETURN
        DO 20 J=J1,J2
          JA = AR(J)
          VC(JA) = VC(JA) + 1
   20   CONTINUE
   30 CONTINUE
      PC(1) = 0
      DO 40 I=1,N
        IF (VC(I).EQ.0) RETURN
        PC(I+1) = PC(I) + VC(I)
        VC(I) = 0
   40 CONTINUE
      DO 60 I=1,N
        J1 = PR(I) + 1
        J2 = PR(I+1)
        DO 50 J=J1,J2
          JJ = AR(J)
          VC(JJ) = VC(JJ) + 1
          JA = PC(JJ) + VC(JJ)
          AC(JA) = I
   50   CONTINUE
   60 CONTINUE
      MP1 = M + 1
C SELECT AS ROOT  JR  THE VERTEX  I  WITH MAXIMUM  VC(I)
C (BREAK TIES BY CHOOSING  I  WITH MINIMUM  VR(I) ).
      MAXE = VC(1)
      MINU = VR(1)
      JR = 1
      DO 100 I=2,N
        IF (VC(I)-MAXE) 100, 70, 80
   70   IF (VR(I).GE.MINU) GO TO 100
        GO TO 90
   80   MAXE = VC(I)
   90   MINU = VR(I)
        JR = I
  100 CONTINUE
      K1 = -NP1
      K = 1
      S(1) = JR
C
C S T E P   1   (SEARCH FOR IMPLIED ARCS).
C
  110 DO 120 J=1,N
        IF (VR(J).EQ.1) GO TO 130
        IF (VC(J).EQ.1) GO TO 170
  120 CONTINUE
C NO FURTHER ARC CAN BE IMPLIED.
      GO TO 220
C ARC  (J,JL)  IS IMPLIED BECAUSE  VR(J) = 1 .
  130 L1 = PR(J) + 1
      L2 = PR(J+1)
      DO 140 L=L1,L2
        IF (AR(L).GT.0) GO TO 150
  140 CONTINUE
  150 JL = AR(L)
C FIND THE STARTING VERTEX  IT1  AND THE ENDING VERTEX  IT2
C OF THE LARGEST PATH OF IMPLIED ARCS CONTAINING  (J,JL) .
      CALL PATH(J, JL, SUBR, RBUS, AR, PR, S, N, NP, IT1, IT2, K, JR,
     * M, NP1)
      IF (NP.EQ.0) GO TO 160
      IF (NP.EQ.(-1)) GO TO 340
C SUBROUTINE PATH FOUND A HAMILTONIAN CIRCUIT.
      K = K + 1
      GO TO 320
  160 SUBR(J) = K1 - JL
      RBUS(JL) = J
C REMOVE FROM THE GRAPH ALL ARCS TERMINATING AT  JL .
      CALL IUPD(J, JL, L, AC, AR, PC, PR, VC, VR, K1, N, M, NP1)
      IF (J.EQ.0) GO TO 340
      GO TO 210
C ARC  (JL,J)  IS IMPLIED BECAUSE  VC(J) = 1 .
  170 L1 = PC(J) + 1
      L2 = PC(J+1)
      DO 180 L=L1,L2
        IF (AC(L).GT.0) GO TO 190
  180 CONTINUE
  190 JL = AC(L)
C FIND THE STARTING VERTEX  IT1  AND THE ENDING VERTEX  IT2
C OF THE LARGEST PATH OF IMPLIED ARCS CONTAINING  (JL,J) .
      CALL PATH(JL, J, SUBR, RBUS, AR, PR, S, N, NP, IT1, IT2, K, JR,
     * M, NP1)
      IF (NP.EQ.0) GO TO 200
      IF (NP.EQ.(-1)) GO TO 340
C SUBROUTINE PATH FOUND A HAMILTONIAN CIRCUIT.
      I = S(K)
      K = K + 1
      GO TO 320
  200 SUBR(JL) = K1 - J
      RBUS(J) = JL
C REMOVE FROM THE GRAPH ALL ARCS EMANATING FROM  JL .
      CALL IUPD(J, JL, L, AR, AC, PR, PC, VR, VC, K1, N, M, NP1)
      IF (J.EQ.0) GO TO 340
C IF ARC  (IT2,IT1)  IS IN THE GRAPH, REMOVE IT.
  210 CALL RARC(IT2, IT1, AR, AC, PR, PC, VR, VC, K1, JJ, LL, N, M, NP1)
      IF (JJ.EQ.(-1)) GO TO 340
      GO TO 110
C
C S T E P   2   (ADD IMPLIED ARCS TO  S ).
C
  220 I = S(K)
      IF (SUBR(I).EQ.0) GO TO 230
      JSUBR = -SUBR(I) + SUBR(I)/NP1*NP1
      IF (JSUBR.EQ.JR) GO TO 340
      K = K + 1
      S(K) = JSUBR
      IF (K.NE.N) GO TO 220
      IF (SUBR(JSUBR).LT.0) GO TO 320
      GO TO 340
C
C S T E P   3   (BRANCH).
C
  230 L1 = PR(I) + P(I)
      L2 = PR(I+1)
      IF (L1.GT.L2) GO TO 340
C FIND THE NEXT ARC  (I,JL)  TO BE ADDED TO  S .
      DENS = N**3
      J1 = 0
      J2 = 0
      DO 310 J=L1,L2
        JL = AR(J)
        IF (JL.LT.0) GO TO 310
        IF (VR(JL).GT.0) GO TO 260
        IF (SUBR(JL).EQ.0) GO TO 310
        IF (JL.EQ.JR) GO TO 310
        IEND = JL
  240   IEND = -SUBR(IEND) + SUBR(IEND)/NP1*NP1
        IF (SUBR(IEND).NE.0) GO TO 240
        IF (VC(JL).LT.VR(IEND)) GO TO 250
        SCORE = VR(IEND)*N + VC(JL)
        GO TO 280
  250   SCORE = VC(JL)*N + VR(IEND)
        GO TO 280
  260   IF (VC(JL).LT.VR(JL)) GO TO 270
        SCORE = VR(JL)*N + VC(JL)
        GO TO 280
  270   SCORE = VC(JL)*N + VR(JL)
  280   IF (DENS.LE.SCORE) GO TO 290
        DENS = SCORE
        IPI = J
  290   IF (J1.EQ.0) GO TO 300
        IF (J2.EQ.0) J2 = J
        GO TO 310
  300   J1 = J
  310 CONTINUE
      IF (J1.EQ.0) GO TO 340
      JL = AR(IPI)
      AR(IPI) = AR(J1)
      AR(J1) = JL
      IF (J2.EQ.0) J2 = PR(I+1) + 1
      P(I) = J2 - PR(I)
      K = K + 1
      S(K) = JL
      K1 = -K*NP1
C REMOVE FROM THE GRAPH ALL ARCS EMANATING FROM  I .
      CALL FUPD(AR, AC, PR, PC, VR, VC, I, K1, N, M, NP1)
C REMOVE FROM THE GRAPH ALL ARCS TERMINATING AT JL .
      CALL FUPD(AC, AR, PC, PR, VC, VR, JL, K1, N, M, NP1)
      TOR(K) = 0
C IF ARC  (JL,JR)  IS IN THE GRAPH, REMOVE IT.
      CALL RARC(JL, JR, AR, AC, PR, PC, VR, VC, K1, JJ, LL, N, M, NP1)
      IF (JJ.EQ.0) GO TO 110
      IF (JJ.EQ.(-1)) GO TO 340
      TOR(K) = JJ*MP1 + LL
      GO TO 110
C
C S T E P   4   (HAMILTONIAN CIRCUIT FOUND).
C
  320 NC = NC + 1
      IF (KW.EQ.(-1)) GO TO 330
      WRITE (KW,99999) (S(KJ),KJ=1,N)
99999 FORMAT (20I5)
  330 IF (NC.EQ.NCO) GO TO 430
      K = K - 1
C
C S T E P   5   (BACKTRACK).
C
  340 IF (K.LE.1) GO TO 430
      JA = S(K)
      P(JA) = 1
      JA = S(K-1)
      IF (SUBR(JA).EQ.0) GO TO 350
C BACKTRACKING FOR AN IMPLIED ARC.
      K = K - 1
      GO TO 340
  350 IF (NB.EQ.NBO) GO TO 430
      NB = NB + 1
      K1 = -K*NP1
      K2 = -(K+1)*NP1
      I = S(K-1)
C BACKTRACKING FOR THE ARCS IMPLIED AT LEVEL  K .
      IFF = 0
      DO 360 J=1,N
        IF (SUBR(J).GT.K1) GO TO 360
        IF (SUBR(J).LT.K2) GO TO 360
        JA = K1 - SUBR(J)
        RBUS(JA) = 0
        SUBR(J) = 0
        IFF = 1
  360 CONTINUE
      IF (IFF.EQ.1) GO TO 370
C NO ARC WAS IMPLIED AT LEVEL  K .
      CALL BUPD(AR, AC, PR, PC, VR, VC, I, K1, K2, N, M, NP1)
      CALL BUPD(AC, AR, PC, PR, VC, VR, S(K), K1, K2, N, M, NP1)
      IF (TOR(K).EQ.0) GO TO 420
      J1 = TOR(K)/MP1
      J2 = TOR(K) - J1*MP1
      AR(J1) = JR
      JA = S(K)
      VR(JA) = VR(JA) + 1
      AC(J2) = S(K)
      VC(JR) = VC(JR) + 1
      GO TO 420
C AT LEAST ONE ARC WAS IMPLIED AT LEVEL  K .
  370 DO 410 J=1,N
        L1 = PR(J) + 1
        L2 = PR(J+1)
        DO 400 L=L1,L2
          JL = AR(L)
          IF (JL.GT.K1) GO TO 400
          IF (JL.LT.K2) GO TO 400
          JL = K1 - JL
          AR(L) = JL
          VR(J) = VR(J) + 1
          LL1 = PC(JL) + 1
          LL2 = PC(JL+1)
          DO 380 LL=LL1,LL2
            IF (K1-AC(LL).EQ.J) GO TO 390
  380     CONTINUE
  390     AC(LL) = J
          VC(JL) = VC(JL) + 1
  400   CONTINUE
  410 CONTINUE
  420 K = K - 1
      GO TO 230
C
C RE-STORE THE ORIGINAL VECTOR  AR .
C
  430 DO 440 J=1,M
        IF (AR(J).GT.0) GO TO 440
        AR(J) = -AR(J) + AR(J)/NP1*NP1
  440 CONTINUE
      RETURN
      END
      SUBROUTINE PATH(I, J, SUBR, RBUS, AR, PR, S, N, NP, I1, I2, K,    PAT   10
     * JR, M, NP1)
C SUBROUTINE TO FIND THE STARTING VERTEX  I1  AND THE ENDING
C VERTEX  I2  OF THE LARGEST PATH OF IMPLIED ARCS CONTAINING
C ARC  (I,J) .
C MEANING OF THE OUTPUT PARAMETER  NP :
C NP =  0  IF THE PATH CONTAINS  L .LT. N  VERTICES.
C    =  1  IF THE PATH CONTAINS  N  VERTICES AND ARC (I2,I1)
C          EXISTS (THE HAMILTONIAN CIRCUIT IS STORED IN  S )
C    = -1  IF THE PATH CONTAINS  N  VERTICES BUT ARC (I2,I1)
C          DOES NOT EXIST.
      INTEGER SUBR(N), RBUS(N), AR(M), PR(NP1), S(N)
      NP = 0
      L = 1
      I1 = I
   10 IF (RBUS(I1).EQ.0) GO TO 20
      I1 = RBUS(I1)
      L = L + 1
      GO TO 10
   20 I2 = J
      L = L + 1
   30 IF (SUBR(I2).EQ.0) GO TO 40
      I2 = -SUBR(I2) + SUBR(I2)/NP1*NP1
      L = L + 1
      GO TO 30
   40 CONTINUE
      IF (L.LT.N) RETURN
C THE PATH CONTAINS  N  VERTICES.
      K1 = -K*NP1
      L1 = PR(I2) + 1
      L2 = PR(I2+1)
      DO 60 L=L1,L2
        IF (AR(L).LT.0) GO TO 50
        IF (AR(L).EQ.I1) GO TO 70
        GO TO 60
   50   IF (K1-AR(L).EQ.I1) GO TO 70
   60 CONTINUE
C NO HAMILTONIAN CIRCUIT CAN BE DETERMINED.
      NP = -1
      RETURN
C A HAMILTONIAN CIRCUIT EXISTS. STORE IT IN  S .
   70 NP = 1
      RBUS(J) = I
      RBUS(I1) = I2
      S(N) = RBUS(JR)
      L = N - 1
   80 IF (L.EQ.K) GO TO 90
      JA = S(L+1)
      S(L) = RBUS(JA)
      L = L - 1
      GO TO 80
   90 RBUS(I1) = 0
      RBUS(J) = 0
      RETURN
      END
      SUBROUTINE FUPD(A1, A2, P1, P2, V1, V2, I1, K1, N, M, NP1)        FUP   10
C FORWARD STEP UPDATING
      INTEGER A1(M), A2(M), P1(NP1), P2(NP1), V1(N), V2(N)
      J1 = P1(I1) + 1
      J2 = P1(I1+1)
      DO 30 J=J1,J2
        IF (A1(J).LT.0) GO TO 30
        IA = A1(J)
        L1 = P2(IA) + 1
        L2 = P2(IA+1)
        DO 10 L=L1,L2
          IF (A2(L).EQ.I1) GO TO 20
   10   CONTINUE
   20   V2(IA) = V2(IA) - 1
        A2(L) = K1 - A2(L)
        A1(J) = K1 - IA
   30 CONTINUE
      V1(I1) = 0
      RETURN
      END
      SUBROUTINE BUPD(A1, A2, P1, P2, V1, V2, II, K1, K2, N, M, NP1)    BUP   10
C BACKTRACKING STEP UPDATING
      INTEGER A1(M), A2(M), P1(NP1), P2(NP1), V1(N), V2(N)
      L1 = P1(II) + 1
      L2 = P1(II+1)
      DO 30 L=L1,L2
        IF (A1(L).GT.K1) GO TO 30
        IF (A1(L).LT.K2) GO TO 30
        IA = K1 - A1(L)
        A1(L) = IA
        V1(II) = V1(II) + 1
        LL1 = P2(IA) + 1
        LL2 = P2(IA+1)
        DO 10 LL=LL1,LL2
          IF (K1-A2(LL).EQ.II) GO TO 20
   10   CONTINUE
   20   A2(LL) = II
        V2(IA) = V2(IA) + 1
   30 CONTINUE
      RETURN
      END
      SUBROUTINE IUPD(IA, IB, L, A1, A2, P1, P2, V1, V2, K1, N, M, NP1) IUP   10
      INTEGER A1(M), A2(M), P1(NP1), P2(NP1), V1(N), V2(N)
C UPDATING FOR IMPLIED ARC
      M1 = P1(IB) + 1
      M2 = P1(IB+1)
      DO 40 MM=M1,M2
        IARC = A1(MM)
        IF (IARC.LT.0) GO TO 40
        IF (V2(IARC).NE.1) GO TO 10
        IF (IARC.NE.IA) GO TO 50
        JJ = L
        GO TO 30
   10   J1 = P2(IARC) + 1
        J2 = P2(IARC+1)
        DO 20 JJ=J1,J2
          IF (A2(JJ).EQ.IB) GO TO 30
   20   CONTINUE
   30   A2(JJ) = K1 - A2(JJ)
        V2(IARC) = V2(IARC) - 1
        A1(MM) = K1 - IARC
        V1(IB) = V1(IB) - 1
   40 CONTINUE
      RETURN
   50 IA = 0
      RETURN
      END
      SUBROUTINE RARC(IA, IB, AR, AC, PR, PC, VR, VC, K1, JJ, LL, N, M, RAR   10
     * NP1)
C SUBROUTINE TO REMOVE ARC  (IA,IB)  FROM THE GRAPH.
C MEANING OF THE OUTPUT PARAMETERS  JJ  AND  LL :
C JJ =  LOCATION OF THE ELEMENT OF  AR  CORRESPONDING TO THE
C       REMOVED ARC.
C    =  0  IF ARC  (IA,IB)  IS NOT IN THE GRAPH.
C    = -1  IF, AFTER THE REMOVAL OF ARC  (IA,IB) , THE GRAPH
C          WOULD ADMIT NO HAMILTONIAN CIRCUIT.
C LL =  LOCATION OF THE ELEMENT OF  AC  CORRESPONDING TO THE
C       REMOVED ARC (DEFINED ONLY IF  JJ .GT. 0 ).
      INTEGER AR(M), AC(M), PR(NP1), PC(NP1), VR(N), VC(N)
      J1 = PR(IA) + 1
      J2 = PR(IA+1)
      DO 20 JJ=J1,J2
        IF (AR(JJ).LT.0) GO TO 20
        IF (AR(JJ).NE.IB) GO TO 20
        L1 = PC(IB) + 1
        L2 = PC(IB+1)
        DO 10 LL=L1,L2
          IF (AC(LL).EQ.IA) GO TO 30
   10   CONTINUE
   20 CONTINUE
C ARC  (IA,IB)  IS NOT IN THE GRAPH.
      JJ = 0
      RETURN
   30 IF (VR(IA).EQ.1) GO TO 40
      IF (VC(IB).EQ.1) GO TO 40
      AR(JJ) = K1 - IB
      VR(IA) = VR(IA) - 1
      AC(LL) = K1 - IA
      VC(IB) = VC(IB) - 1
      RETURN
C ARC  (IA,IB)  CANNOT BE REMOVED FROM THE GRAPH.
   40 JJ = -1
      RETURN
      END
    6
    0    3    5    8   11   13   16
    3    5    6    6    3    6    1    4    5    1
    2    2    3    4    2    5
