      SUBROUTINE BBUPDT ( N, G, S, XX, Y, H, ETA, CT, CNTRST, LASTPT,
     -                   IDENTY, NUPS, STEEPD, RSTEP, QNPART, ALPHA,
     -                    STG, UTG, NU, UPDATT, INNER, IW, RW, DW      )

C## A R G U M E N T S:
                      INTEGER N, CT, CNTRST, LASTPT, NUPS, UPDATT, IW(*)
                      LOGICAL STEEPD, RSTEP, QNPART, IDENTY

      DOUBLE PRECISION G(N), S(N), XX(N), Y(N), H(*)
C!!!! REAL             G(N), S(N), XX(N), Y(N), H(*)

      DOUBLE PRECISION ALPHA, ETA, STG, UTG, NU
C!!!! REAL             ALPHA, ETA, STG, UTG, NU

      EXTERNAL                INNER
      DOUBLE PRECISION DW(*), INNER
      REAL             RW(*)

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!           THIS VERSION IS IN   S I N G L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.

C>RCS $HEADER: UPDT.F,V 2.3 91/12/31 14:53:08 BUCKLEY EXP $
C>RCS $LOG:     UPDT.F,V $
C>RCS REVISION 2.3  91/12/31  14:53:08  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.2  91/11/22  11:33:52  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.1  90/08/02  13:38:57  BUCKLEY
C>RCS REMOVED DOUBLE DECLARATIONS; MINOR FIX
C>RCS
C>RCS REVISION 2.0  90/07/31  12:04:18  BUCKLEY
C>RCS FIXED LONG LINES.
C>RCS
C>RCS REVISION 1.10  90/07/31  10:50:16  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:12:52  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:39:27  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:55:39  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:54:30  BUCKLEY
C>RCS INITIAL REVISION
C>RCS

C## D E S C R I P T I O N:
C
C     THE BASIC PURPOSE OF THIS ROUTINE IS TO COMPUTE THE VALUE OF
C     THE UPDATE MATRIX  H  FOR THE NEW POINT.
C
C     NOTE THAT THERE ARE SEVERAL VARIABLES DEFINED IN THE MAIN
C     ROUTINE BBLNIR WHICH AFFECT THIS ROUTINE. HOWEVER, SINCE
C     THEY ARE INVARIANT BETWEEN CALLS TO BBUPDT, THEY ARE SET ONCE
C     WITH A CALL TO THE ENTRY POINT BBSUPT AND THEY ARE RETAINED FROM
C     CALL TO CALL WITH A NUMBER OF SAVE VARIABLES.
C
C     ON ENTRY THE FOLLOWING VALUES ARE REQUIRED.
C
C       N      THE DIMENSION OF THE PROBLEM.
C       G      THE GRADIENT AT THE NEW POINT X.
C       S      THE STEP TAKEN ON THE ITERATION JUST COMPLETED.
C       XX     MAY BE USED AS A SCRATCH VECTOR; ALSO SEE BELOW.
C       Y      THE CHANGE IN GRADIENT FROM THE PREVIOUS POINT.
C              THIS MAY ALSO BE USED AS A SCRATCH VECTOR.
C       ETA    = S'Y OVER LAST STEP
C       H      THE CURRENT MATRIX H.
C
C       INNER  SEE THE DISCUSSION IN BBMULT.
C
C     IN ADDITION, IF CG IS TRUE (SEE THE ENTRY POINT BBSUPT),
C     THE FOLLOWING VALUES MUST BE DEFINED.
C
C       CT     THE ITERATION NUMBER
C       CNTRST THE NUMBER OF RESTARTS DONE.
C       NUPS   THE NUMBER OF TERMS DEFINING THE CURRENT UPDATE MATRIX.
C       ALPHA  STEP LENGTH; NEEDED BY POWELL
C       STEEPD A FLAG WHICH IS TRUE WHEN THE LAST SEARCH DIRECTION
C              WAS ALONG THE DIRECTION OF STEEPEST DESCENT.
C       RSTEP  A FLAG WHICH IS TRUE WHEN THIS IS A RESTART POINT.
C              THIS FLAG WILL ALWAYS BE FALSE WHEN CG IS FALSE.
C       IDENTY TRUE IF H[0] = I.
C       UPDATT TO INDICATE WHAT TYPE OF QN UPDATES ARE BEING STORED,
C              I.E. SUM FORM (BBLNIR) OR PRODUCT FORM (NOCEDAL)
C              OR FACTORED FORM (POWELL).
C       QNPART  A FLAG WHICH IS TRUE WHEN WE ARE IN THE QUASI-NEWTON
C               PART OF THE CODE.
C
C     ON EXIT FROM BBUPDT, THE MATRIX H MUST HAVE BEEN UPDATED.
C
C     IN THE QUASI-NEWTON CASE, THAT UPDATE WILL HAVE BEEN DONE IN
C     PLACE, I.E. THE NEW MATRIX H WILL JUST HAVE OVERWRITTEN THE OLD.
C
C     IN THE CONJUGATE GRADIENT CASE (I.E. WHEN CG IS TRUE), ANOTHER
C     TERM WILL HAVE BEEN ADDED TO THE SUM FORM OF H (UNLESS ALL
C     THE SPACE FOR UPDATES HAS BEEN USED), OR ELSE, IN THE EVENT
C     OF A RESTART, H WILL HAVE REDEFINED BY A SINGLE UPDATE TERM.
C     THUS, THE FOLLOWING VALUES MUST BE SET BEFORE RETURNING:
C
C       CNTRST MUST HAVE BEEN INCREMENTED BY 1 IF A RESTART WAS DONE.
C       IDENTY WILL BE SET TO TRUE WHENEVER H0 IS THE IDENTITY, EVEN
C              IF SCDIAG IS TRUE.
C       NUPS   MUST HAVE BEEN REVISED TO THE NUMBER OF SUM TERMS
C              DEFINED BY THE NEW H, WHETHER 1 OR AN INCREMENT OF THE
C              PREVIOUS VALUE.
C       CT     MUST BE RESET IF THE UPDATE WAS A RESTART. THIS IS THE
C              ACTUAL ITERATION COUNTER, WHICH STARTS FROM 1 AT A
C              RESTART POINT AND IS INCREMENTED FOR EACH NEW POINT.
C       LASTPT MUST BE SET, IF A RESTART, TO INDICATE THE NEXT POINT
C              AT WHICH A RESTART MUST BE FORCED, REGARDLESS OF THE
C              TESTING MECHANISM.
C       XX, STG, UTG, NU
C
C    IN THE CASE THAT PRODUCT FORM UPDATES ARE BEING STORED, THESE
C    VALUES MUST ALSO BE UPDATED, BUT THERE ARE SOME NOTABLE DIFF-
C    ERENCES. THERE ARE NO RESTARTS, AND WHEN THE MEMORY LIMIT IS
C    REACHED, EARLIER UPDATE TERMS ARE SIMPLY OVERWRITTEN IN A CIRCU-
C    LAR FASHION.
C
C    THE TRACE VARIABLES TR7, TR8 AND TRV ARE EXPLAINED WITHIN THE
C    DESCRIPTION OF BBLNIR.
C
C## E N T R Y   P O I N T S: BBUPDT   THE NATURAL ENTRY POINT.
C                            BBSUPD   INITIALIZE FIXED ARGUMENTS.

C## S U B R O U T I N E S: BBMULT    TO MULTIPLY BY A SUM FORM H.
C                          MOD, SQRT INTRINSICS.
C                          INNER     AN EXTERNAL ARGUMENT.

C## P A R A M E T E R S:

       INTEGER     SUMFRM,     PRDFRM,     MJDFRM
       PARAMETER ( SUMFRM = 1, PRDFRM = 2, MJDFRM = 3 )
      DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
C!!!! REAL              ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
C!!!! REAL              FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      DOUBLE PRECISION  EIGHT,         NINE,          TEN
C!!!! REAL              EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

      LOGICAL     DONORM,          NONORM
      PARAMETER ( DONORM = .TRUE., NONORM = .FALSE. )

      LOGICAL     T,          F
      PARAMETER ( T = .TRUE., F = .FALSE. )

      CHARACTER*(*) TRUE,          QT,       FALSE,           QF
      PARAMETER (   TRUE = 'TRUE', QT = 'T', FALSE = 'FALSE', QF = 'F' )

      INTEGER     ITRUE,     IFALSE
      PARAMETER ( ITRUE = 1, IFALSE = 0 )

      DOUBLE PRECISION  RTRUE,        RFALSE
C!!!! REAL              RTRUE,        RFALSE
      PARAMETER      (  RTRUE = 1.D0, RFALSE = 0.D0 )

      INTEGER     GAMOFF,     GAMONE,     GAMALL
      PARAMETER ( GAMOFF = 0, GAMONE = 1, GAMALL = 2 )

C## L O C A L   D E C L:

C-----CONTROL PARAMETERS FOR ENTRY BBSUPD.

      INTEGER   M,  BASE,  INCR, SCGAMM,  TRU
      INTEGER  SM, SBASE, SINCR, SSGAMM, STRU

      LOGICAL   CG, SCDIAG, USESHN, FROMRS,  TR7,  TR10,  TRV
      LOGICAL  SCG,  SDIAG, SUSEHN, SFRMRS, STR7, STR10, STRV

C-----GENERAL DECLARATIONS.

      INTEGER          IETA, INU, IS, IU, J, K, KJ, IY, ID, IG

      DOUBLE PRECISION TP1, TP2, R
C!!!! REAL             TP1, TP2, R

C## S A V E:
                       SAVE  M, BASE, INCR, INU, SCGAMM, TRU, CG,
     -                       SCDIAG, USESHN, FROMRS, TR7, TR10, TRV

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NO DATA VALUES SET.
C##                                            E X E C U T I O N
C##                                            E X E C U T I O N
C
C     IN THIS DESCRIPTION, "H" WILL DENOTE THE UPDATE MATRIX DEFINED
C     WHEN THE CURRENT POINT IS REACHED; "H^" WILL DENOTE THE UPDATE
C     MATRIX TO BE COMPUTED AND USED IN FORMING THE NEXT SEARCH
C     DIRECTION.

      IF ( TR7 ) WRITE(TRU,*) ' ***[ENTERING UPDT]***'
      IF ( CG ) THEN

          IF ( RSTEP ) THEN

C            P H A S E  X - A : R E S T A R T.<<<<<<<<<<<<<<<<<<<

             IF (TR7 .OR. TR10) WRITE(TRU,*) ' [UPDT] RESTART! NUPS->1'
C            COUNT NUMBER OF RESTARTS
             CNTRST = CNTRST + 1

C                                   SET POINT TO FORCE THE NEXT RESTART.
             IF ( FROMRS .AND. .NOT. USESHN ) THEN
                LASTPT = 1 + N
             ELSE
                LASTPT = M + 1 + N
             ENDIF

C            AFTER A RESTART, THE DIAGONAL SCALING MATRIX IS ALWAYS JUST  I.

             IDENTY = T

C            IF M = 0, WE CAN NOT SAVE UPDATES SO REVERT TO A STEEPEST DESCENT
C            RESTART.

             IF ( M .EQ. 0 ) THEN
                CT   = 0
                NUPS = 0
             ELSE
                CT   = 1
                NUPS = 1
                QNPART = T
             ENDIF

C            SINCE A RESTART IS INDICATED, SAVE THE CURRENT S AND U = H*Y
C            = I*Y = Y (THE BEALE RESTART VECTORS) AND SAVE  NU = Y'*H*Y
C            = Y'Y AND ETA = S'Y IN H(INU) AND H(IETA), I.E. DEFINE H[1].

             IF (  M .NE. 0  ) THEN
                INU   = BASE + 1
                IETA  = INU  + 1
                IU    = IETA
                IS    = IU   + N
                CALL ZZCOPY ( N, Y, 1, H(IU+1), 1 )
                CALL ZZCOPY ( N, S, 1, H(IS+1), 1 )

                H(INU)  = INNER ( N, H(IU+1), H(IU+1), NONORM,IW,RW,DW)
                H(IETA) = ETA
                IF(TR7)
     -            WRITE(TRU,*) ' [UPDT] SAVED NU, ETA->',H(INU),H(IETA)
             ENDIF

          ELSE IF       ( CG .AND.  UPDATT .EQ. SUMFRM ) THEN

C            P H A S E  X - B: CG, SUM FORM UPDATE. <<<<<<<<<

C            COMPUTE U = H*Y (INTO XX).
C            ALSO COMPUTE  S'G (IN STG) AND U'G (IN UTG), AS WELL
C            AS ACCUMULATING NU = Y'*H*Y  AND  ETA = S'Y.  NOTE THAT THE
C            COMPUTATION IS THE SAME FOR THE CG OR QN PARTS.

             IF ( TR7 ) WRITE(TRU,*) ' ***DOING SUM FORM UPDATE'

             CALL BBMULT ( H, Y, XX, N, NUPS, T, IDENTY, INNER,IW,RW,DW)

             NU  = INNER( N,  Y, XX, NONORM, IW, RW, DW )
             STG = INNER( N,  S,  G, NONORM, IW, RW, DW )
             UTG = INNER( N, XX,  G, NONORM, IW, RW, DW )

             IF ( QNPART ) THEN
C               SAVE UPDATE TERMS: PUT  NU,ETA,U AND S IN THE ARRAY H.
                NUPS = NUPS + 1
                INU  = INU + INCR
                IETA = INU + 1
                IU   = IETA
                IS   = IU + N
                CALL ZZCOPY ( N, XX, 1, H(IU+1), 1 )
                CALL ZZCOPY ( N, S,  1, H(IS+1), 1 )
                H(INU)  = NU
                H(IETA) = ETA
                IF ( TR7 .OR. TR10 ) WRITE(TRU,*) ' [UPDT] NO RESTART;'
     -                                       //' NUPS->',NUPS
                IF ( TR7 ) WRITE(TRU,*) ' [UPDT] SAVED NU, ETA->',NU,ETA
             ELSE
                IF ( TR7 .OR. TR10 ) WRITE(TRU,*) ' [UPDT] NO RESTART;'
     -                                 //' NUPS->',NUPS+1,'(NOT STORED)'
             ENDIF
C                  ...FOR THE  "IF QNPART SO SAVE...".
          ELSE IF       ( CG .AND.  UPDATT .EQ. PRDFRM ) THEN

C            P H A S E  X - C: CG, PRODUCT FORM UPDATE. <<<<<<<<<

             IF ( TR7 ) WRITE(TRU,*) ' ***DOING PRODUCT FORM UPDATE'
             NUPS = MOD(NUPS,M)
             IETA = BASE + NUPS*INCR + 1
             IS   = IETA
             IY   = IS + N
             NUPS = NUPS + 1

             CALL ZZCOPY ( N, S, 1, H(IS+1), 1 )
             CALL ZZCOPY ( N, Y, 1, H(IY+1), 1 )

             H(IETA) = ETA

             IF ( TR7 .OR. TR10 ) WRITE(TRU,*) ' [UPDT] SAVED NOCEDAL'
     -                                       //' UPDATE TERM.'

          ELSEIF        ( CG .AND.  UPDATT .EQ. MJDFRM ) THEN

C            P H A S E  X - D: CG, FACTORED FORM UPDATE. <<<<<<<<

             IF ( TR7 ) WRITE(TRU,*) ' ***DOING FACTORED FORM UPDATE'
             IY   = BASE + NUPS*N
             NUPS = NUPS + 1
             ID   = IY   + N*M
             IG   = ID   + N*M

             R = SQRT ( ETA )
             CALL ZZSCAL ( N, -ALPHA,  H(IY+1), 1 )
             CALL ZZCOPY ( N, S, 1,    H(ID+1), 1 )
             CALL ZZSCAL ( N, (ONE/R), H(ID+1), 1 )
             CALL ZZCOPY ( N, Y, 1,    H(IG+1), 1 )
             CALL ZZSCAL ( N, (ONE/R), H(IG+1), 1 )

                IF ( TR7 .OR. TR10 ) WRITE(TRU,*) ' [UPDT] SAVED POWELL'
     -                                       //' UPDATE TERM.'
         ENDIF
      ELSE

C        P H A S E  X - E:   Q N   C A S E.<<<<<<<<<<<<<<<<<<

C        A VARIABLE METRIC ALGORITHM IS BEING USED. CALCULATE GRADIENT
C        DIFFERENCE  Y  AND ETA = S'Y.  S  IS THE STEP .

C        IF STEEPD IS .TRUE., SET UP THE INITIAL SCALED APPROXIMATE
C        HESSIAN. THIS IS THE INITIAL STEP.

         IF ( STEEPD ) THEN

C           CALCULATE  NU = Y'*H*Y, WHICH HERE IS  NU = Y'*Y.

            NU = INNER ( N, Y, Y, NONORM, IW, RW, DW )

C           STORE THE INITIAL HESSIAN, WHICH IS  H = (S'Y/Y'Y)*I =
C           (ETA/NU)*I. THEN WE NEED TO RECALCULATE THE INITIAL NU =
C           Y'*H*Y = (ETA/NU)*(NU ABOVE) = ETA, AND TO FIND XX = H*Y.

            KJ  = 1
            TP1 = ETA/NU

            DO 6000 K=1,N

C              NOTE: INNER LOOP IS FROM K TO N SO ONLY HALF OF H.

               H(KJ) = TP1
               KJ    = KJ + 1

               DO 5900 J = K+1,N
                  H(KJ) = ZERO
                  KJ    = KJ + 1
 5900          CONTINUE
               XX(K) = TP1*Y(K)
 6000       CONTINUE
            NU = ETA
         ELSE

C           CALCULATE XX = H*Y AND NU = Y'*H*Y.  REMEMBER THAT ONLY
C           THE SYMMETRIC UPPER HALF OF H IS STORED (IN ROW ORDER).

            NU = ZERO
            DO 6500 K = 1,N
               TP1 = ZERO
               KJ  = K
               DO 6200 J=1,K-1
                  TP1 = TP1 + H(KJ)*Y(J)
                  KJ  = KJ  + (N-J)
 6200          CONTINUE
               DO 6400 J=K,N
                  TP1 = TP1 + H(KJ)*Y(J)
                  KJ  = KJ+1
 6400          CONTINUE
               NU    = NU + TP1*Y(K)
               XX(K) = TP1
 6500        CONTINUE
          ENDIF
C               ...FOR " IF STEEPD".

C         NOW CALCULATE THE UPDATED APPROXIMATE HESSIAN H^.
C         USE THE BFGS UPDATE. NU, ETA AND H*Y (IN XX) ARE KNOWN.

          TP1 = ONE + NU/ETA

          CALL ZZCOPY ( N, XX,     1, Y, 1 )
          CALL ZZSCAL ( N, -ONE,      Y, 1 )
          CALL ZZAXPY ( N, TP1, S, 1, Y, 1 )

          KJ = 1

          DO 7400 K=1,N

             TP2 =  S(K)/ETA
             TP1 = XX(K)/ETA

             DO 7200 J=K,N
                H(KJ) = H(KJ) + TP2*Y(J) - TP1*S(J)
                KJ    = KJ+1
 7200        CONTINUE

 7400     CONTINUE

      E N D I F
C               ...FOR THE UPDATE CHOICES.

      GOTO 90000

C## E N T R Y  BBSUPD:
                      ENTRY BBSUPD ( SM, SBASE, SINCR, SSGAMM,
     -           SCG, SDIAG, SUSEHN, SFRMRS, STR7, STR10, STRV, STRU )
      M      = SM
      BASE   = SBASE
      INCR   = SINCR
      SCGAMM = SSGAMM
      CG     = SCG
      SCDIAG = SDIAG
      USESHN = SUSEHN
      FROMRS = SFRMRS
      TR7    = STR7
      TR10   = STR10
      TRV    = STRV
      TRU    = STRU
      RETURN

C## E X I T
90000        IF ( TR7 ) WRITE(TRU,*) ' ===[LEAVING UPDT].'
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF BBUPDT.
                    END
