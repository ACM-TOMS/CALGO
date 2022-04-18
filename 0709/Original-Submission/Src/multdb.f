      SUBROUTINE BBMULT ( H, V, HV, N, NUPS, CMPALL, IDENTY, INNER,
     -                                              IW, RW, DW    )
C## A R G U M E N T S:
                      INTEGER          N, NUPS, IW(*)
                      LOGICAL          IDENTY, CMPALL

                      DOUBLE PRECISION H(*), V(N), HV(N)
C!!!!                 REAL             H(*), V(N), HV(N)

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
C
C>RCS $HEADER: MULT.F,V 1.13 91/12/31 14:53:05 BUCKLEY EXP $
C>RCS $LOG:     MULT.F,V $
C>RCS REVISION 1.13  91/12/31  14:53:05  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.12  91/12/16  11:37:19  BUCKLEY
C>RCS MINOR FIX FOR TOMS.
C>RCS
C>RCS REVISION 1.11  91/11/19  15:33:35  BUCKLEY
C>RCS FINAL SUBMISSION FOR TOMS
C>RCS
C>RCS REVISION 1.10  90/07/31  10:50:01  BUCKLEY
C>RCS ADDED REVISED BLAS.
C>RCS
C>RCS REVISION 1.9  89/06/30  13:12:49  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3  89/05/18  12:39:23  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:55:34  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:54:29  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     GIVEN THE QUASI NEWTON UPDATE MATRIX  H (IN SUM FORM) AND
C     GIVEN THE VECTOR V, THIS ROUTINE COMPUTES
C
C               HV = H * V  .
C
C-----NOTE THAT A NUMBER OF PARAMETERS WHICH WILL BE THE SAME
C     FOR EACH CALL TO BBMULT DURING ANY PARTICULAR MINIMIZATION
C     PROBLEM ARE SET JUST ONCE THROUGH AN ENTRY POINT.
C
C------EACH UPDATE TERM OF H REQUIRES 2N+2 ENTRIES OF H. THE ORDER IS
C
C                   NU(I), ETA(I),  U(I) AND  S(I).
C
C      EACH BLOCK OF 2N+2 ENTRIES IS CALLED A "TERM" OF THE UPDATE.
C
C      H IS THE MATRIX DEFINED AT THE CURRENT POINT I AND
C      H^ DENOTES THE MATRIX OBTAINED BY UPDATING H.
C
C      HERE    N    = THE DIMENSION OF THE PROBLEM
C              S    = X[I] - X[I-1] = ALPHA * D
C              Y    = G[I] - G[I-1]
C              U    = H * Y
C              NU   = Y' * H * Y
C              ETA  = S' * Y
C
C      NOTE THAT THIS ROUTINE USES  S' TO DENOTE THE TRANSPOSE OF THE
C      COLUMN VECTOR  S, SO THAT S'*Y, FOR EXAMPLE, IS A SCALAR. THE
C      SUBSCRIPT [I] IS DROPPED IN MOST OF THE SUBSEQUENT DESCRIPTION.
C
C      IN FACT, ALL INNER PRODUCTS (AND 2-NORMS) ARE COMPUTED BY CALLING
C      THE PROCEDURE INNER, WHICH IS PASSED AS AN ARGUMENT TO BBMULT.
C      BY DEFAULT THEN, IF ZZINNR IS PASSED IN FOR INNER, NORMAL
C      EUCLIDEAN INNER PRODUCTS AND NORMS ARE OBTAINED FOR S'*Y=(S,Y)
C      AND OTHER INNER PRODUCTS. HOWEVER THE USER MAY REPLACE ZZINNR
C      WITH ANY SUITABLE ROUTINE OF HIS CHOICE. IN THIS CASE, INNER
C      MAY MAKE USE OF THE DATA IN THE VECTORS IW, RW AND DW.
C
C--NUPS = NUMBER OF TERMS IN THE UPDATE MATRIX H.
C
C         IF NUPS = 0,  H IS JUST H0  AND H * V IS JUST H0*V.
C         IN PARTICULAR, IF H0 = I THIS GIVES H*V = V.
C
C--CMPALL IS A FLAG WHICH CONTROLS THE COMPUTATION TO BE DONE.
C      TRUE   COMPUTE H*V   USING ALL THE TERMS WHICH DEFINE H.
C      FALSE  COMPUTE BY ADDING JUST ONE LAST TERM; I.E. WE COMPUTE
C             (H^)*V, ASSUMING THAT H*V WAS DONE EARLIER AND IS IN HV,
C             AND THAT H^ IS THE UPDATE OF H DEFINED BY THE LAST TERM.
C
C--IDENTY IS TRUE TO INDICATE THAT H0 = I; THIS MAY BE TRUE
C         EVEN IF SCDIAG IS TRUE.
C
C-------- IN THE ENTRY POINT...
C
C--BETA   IS THE PARAMETER DEFINING THE BROYDEN FAMILY OF UPDATES.
C         THE FORM USED AT EACH POINT IS
C              H^ = H(DFP) + BETA * NU * W'W
C         SO THAT BETA = 1 GIVES THE BFGS UPDATE.
C
C--SCDIAG IF .TRUE., H0 IS TAKEN TO BE A DIAGONAL MATRIX WHICH IS
C         AVAILABLE IN THE FIRST  N LOCATIONS OF THE ARRAY  H.
C         OTHERWISE  H0 = I, AND IT IS OF COURSE NOT STORED.
C
C--SCGAMM
C  =GAMALL  THEN THE SO-CALLED GAMMA SCALING OF OREN AND
C         SPEDICATO, WHICH IS DESCRIBED BY SHANNO, IS USED AT EACH
C         UPDATE STEP. THIS CAN IN FACT BE DONE ONLY IF THE BFGS
C         UPDATE IS BEING USED, I.E. IF BETA = 1. NO EXTRA STORAGE
C         IS NEEDED TO IMPLEMENT THIS SCALING.
C
C  =GAMONE  THEN SCALING IS DONE, AS JUST DESCRIBED FOR SCGAMM=GAMALL,
C         BUT IT ONLY APPLIES TO THE FIRST UPDATE TERM.
C
C  =GAMOFF  THEN NO SCALING IS DONE.
C
C--INCR   IS THE CONSTANT 2N+2, THE LENGTH OF EACH TERM.
C  BASE   IS THE NO. OF LOCATIONS FOR THE DIAG. H0, EITHER 0 OR N.
C
C--TRACES TURN ON TR TO SEE NU, ETA, GAMMA, HV AND S'V.
C         THESE WILL BE ON THE UNIT TRACUN.
C         VECTORS ARE TRACED ONLY IF TRV IS TRUE AS WELL.
C
C## E N T R Y   P O I N T S: BBMULT  THE NATURAL ENTRY POINT.
C                            BBSMLT  SET FIXED PARAMETERS.
C## S U B R O U T I N E S:   INNER   (ARGUMENT)
C## P A R A M E T E R S:
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

      INTEGER     GAMOFF,     GAMONE,     GAMALL
      PARAMETER ( GAMOFF = 0, GAMONE = 1, GAMALL = 2 )

C## L O C A L   D E C L:

      INTEGER          COUNT, IU, K, PTNU, IS, I, II

      DOUBLE PRECISION NU, ETA, UV, SV, GAMMA, MU, SIGMA
C!!!! REAL             NU, ETA, UV, SV, GAMMA, MU, SIGMA

C-----VARIABLES FOR THE ENTRY POINT.

      LOGICAL     TR,  TRV , SCDIAG
      LOGICAL    STR , STRV, SSCDAG

      INTEGER     SCGAMM,  INCR,  BASE, TRACUN
      INTEGER     SSCGAM, SINCR, SBASE, STRACN

      DOUBLE PRECISION BETA, SBETA
C!!!! REAL             BETA, SBETA

C## S A V E:
            SAVE   TR, TRV, SCDIAG, SCGAMM, INCR, BASE, TRACUN, BETA

C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
      IF ( TR ) WRITE (TRACUN,*) ' ***[MULT ENTERED]***'

C     INITIALIZE COUNTERS AND INITIALIZE FOR VARIOUS  H0.

      IF (.NOT. CMPALL) THEN
         COUNT = NUPS
         PTNU  = BASE + 1 + INCR*(NUPS-1)
      ELSE

C        SET HV = H0 * V, WHERE H0 IS THE INITIAL POSITIVE
C        DEFINITE MATRIX, WHICH MAY BE EITHER THE IDENTITY
C        OR A DIAGONAL SCALING MATRIX.

         IF (SCDIAG .AND. .NOT. IDENTY) THEN
            DO 200 K=1,N
               HV(K) = H(K) * V(K)
  200       CONTINUE
         ELSE
            CALL ZZCOPY ( N, V, 1, HV, 1 )
         ENDIF

         COUNT = 1
         PTNU  = BASE + 1
      ENDIF

C     COMPUTE THE TERMS OF THE PRODUCT.

      DO  4000 I= COUNT, NUPS

         NU  = H(PTNU)
         ETA = H(PTNU+1)
         IU  = PTNU + 1
         IS  = IU   + N

         IF ( TR ) WRITE (TRACUN,*) ' [MULT] NU,ETA,PTNU,NUPS->',
     -                                       NU,ETA,PTNU,NUPS

C         COMPUTE  UV = U' * V  AND  SV = S' * V.

         UV = ZERO
         SV = ZERO
C$DOACROSS SHARE(H,V,IU,IS), REDUCTION(UV,SV), CHUNK=2000
C        DO 3000 II = 1,N
C           UV = UV + H(IU+II)*V(II)
C           SV = SV + H(IS+II)*V(II)
C3000    CONTINUE
         UV = INNER ( N, H(IU+1), V, NONORM, IW, RW, DW )
         SV = INNER ( N, H(IS+1), V, NONORM, IW, RW, DW )

         IF ( TR ) WRITE ( TRACUN, * ) ' [MULT] SV->', SV
         IF ( TR ) WRITE ( TRACUN, * ) ' [MULT] UV->', UV

C        COMPUTE NEXT TERM AND ADD INTO HV.  USE GENERAL FORM
C        H(DFP) + BETA* NU*W'*W. BETA = 1 GIVES A BFGS UPDATE.

C        IF GAMMA-SCALING IS REQUIRED, SET GAMMA = ETA/NU, AND USE THE
C        MODIFIED UPDATE FORMULA WHICH CAN BE DERIVED FROM SHANNO'S
C        WORK. AGAIN, THIS ONLY APPLIES TO THE BFGS UPDATE.

         IF ( (BETA .EQ. ONE)
     -         .AND. ( ( SCGAMM .EQ. GAMALL )
     -             .OR. (SCGAMM .EQ. GAMONE .AND. I .EQ. 1)) ) THEN

            GAMMA = ETA/NU
            IF ( TR ) WRITE (TRACUN, * ) ' [MULT] GAMMA->',GAMMA

            CALL ZZSCAL ( N, GAMMA, HV, 1 )

            MU    = - SV/NU
            SIGMA = (TWO*SV/ETA)  - (UV/NU)

         ELSEIF ( BETA .EQ. ONE ) THEN
            MU    = - SV/ETA
            SIGMA = - ( ONE + NU/ETA )*MU - UV/ETA
         ELSE
            MU    = ( (BETA - ONE)*UV/NU ) - ( BETA*SV/ETA )
            SIGMA = SV* (ETA + BETA*NU)/(ETA*ETA) - (BETA*UV/ETA)
         ENDIF
         IF ( TR ) WRITE (TRACUN, * ) ' [MULT] MU,SIGMA->',MU,SIGMA

C$DOACROSS SHARE(MU,SIGMA, IU,IS,HV,H), CHUNK=2000
C        DO II = 1,N
C            HV(II) = HV(II) + MU*H(IU+II) + SIGMA*H(IS+II)
C        ENDDO
         CALL ZZAXPY ( N, MU, H(IU+1), 1, HV, 1 )
         CALL ZZAXPY ( N, SIGMA, H(IS+1), 1, HV, 1 )

         IF ( TRV .AND. TR ) WRITE (TRACUN, * ) ' [MULT] H*V->',HV
         PTNU = PTNU + INCR
 4000 CONTINUE
      GOTO 90000

C## E N T R Y  BBSMLT:
                      ENTRY  BBSMLT ( STR, STRV, SSCDAG, SSCGAM,
     -                             STRACN, SBASE, SINCR, SBETA )
      TR     = STR
      TRV    = STRV
      SCDIAG = SSCDAG
      SCGAMM = SSCGAM
      TRACUN = STRACN
      BASE   = SBASE
      INCR   = SINCR
      BETA   = SBETA
      RETURN

C## E X I T

90000       IF ( TR ) WRITE (TRACUN,*) ' ===[LEAVING MULT].'
      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF BBMULT.
                    END