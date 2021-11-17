! SUBROUTINE PEQNU              ALL SYSTEMS                   97/01/22
! PURPOSE :
! EASY TO USE SUBROUTINE FOR SOLUTION OF SPARSE SYSTEMS OF NONLINEAR
! EQUATIONS USING THE DISCRETE NEWTON METHOD.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  AF(N)  VECTOR CONTAINING VALUES OF THE APPROXIMATED
!         FUNCTIONS.
!  II  IAG(N+1)  POSITION OF THE FIRST ROWS ELEMENTS IN THE JACOBIAN
!         MATRIX.
!  II  JAG(MA) COLUMN INDICES OF ELEMENTS IN THE JACOBIAN MATRIX.
!  II  IPAR(7)  INTEGER PAREMETERS:
!      IPAR(1)  MAXIMUM NUMBER OF ITERATIONS.
!      IPAR(2)  MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
!      IPAR(3)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!      IPAR(4)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!      IPAR(5)  CHOICE OF THE SMOOTHING STRATEGY FOR THE CONJUGATE
!         GRADIENT SQUARED METHOD. IPAR(5)=1-SMOOTHING IS NOT USED.
!         IPAR(5)=2-SINGLE SMOOTHING STRATEGY IS USED. IPAR(5)=3-DOUBLE
!         SMOOTHING STRATEGY IS USED.
!      IPAR(6)  CHOICE OF PRECONDITIONING. IPAR(6)=1-PRECONDITIONING
!         IS NOT USED. IPAR(6)=2-PRECONDITIONING BY THE INCOMPLETE
!         GILL-MURRAY DECOMPOSITION. IPAR(6)=3-PRECONDITIONING BY THE
!         INCOMPLETE GILL-MURRAY DECOMPOSITION WITH A PRELIMINARY
!         SOLUTION OF THE PRECONDITIONED SYSTEM WHICH IS USED IF IT
!         SATISFIES THE TERMINATION CRITERION.
!      IPAR(7)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!  RI  RPAR(9)  REAL PARAMETERS:
!      RPAR(1)  MAXIMUM STEPSIZE.
!      RPAR(2)  TOLERANCE FOR THE CHANGE OF VARIABLES.
!      RPAR(3)  TOLERANCE FOR THE CHANGE OF FUNCTION VALUES.
!      RPAR(4)  TOLERANCE FOR THE FUNCTION FALUE.
!      RPAR(5)  TOLERANCE FOR THE GRADIENT NORM.
!      RPAR(6)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!      RPAR(7)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!      RPAR(8)  DAMPING PARAMETER FOR AN INCOMPLETE LU PRECONDITIONER.
!      RPAR(9)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PEQN.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE.
!  II  IDER  DEGREE OF ANALYTICALLY COMPUTED DERIVATIVES (0 OR 1).
!  II  ISPAS  INPUT SPARSE STRUCTURE. ISPAS=1-STANDARD COORDINATE
!         FORM. ISPAS=2-SPARSE STRUCTURE COMPRESSED BY ROWS.
!  II  IPRNT  PRINT SPECIFICATION. IPRNT=0-NO PRINT.
!         ABS(IPRNT)=1-PRINT OF FINAL RESULTS.
!         ABS(IPRNT)=2-PRINT OF FINAL RESULTS AND ITERATIONS.
!         IPRNT>0-BASIC FINAL RESULTS. IPRNT<0-EXTENDED FINAL
!         RESULTS.
!  IO  ITERM  VARIABLE THAT INDICATES THE CAUSE OF TERMINATION.
!         ITERM=1-IF ABS(X-XO) WAS LESS THAN OR EQUAL TO TOLX IN
!                   MTESX (USUALLY TWO) SUBSEQUENT ITERATIONS.
!         ITERM=2-IF ABS(F-FO) WAS LESS THAN OR EQUAL TO TOLF IN
!                   MTESF (USUALLY TWO) SUBSEQUENT ITERATIONS.
!         ITERM=3-IF F IS LESS THAN OR EQUAL TO TOLB.
!         ITERM=4-IF GMAX IS LESS THAN OR EQUAL TO TOLG.
!         ITERM=6-IF THE TERMINATION CRITERION WAS NOT SATISFIED,
!                   BUT THE SOLUTION OBTAINED IS PROBABLY ACCEPTABLE.
!         ITERM=11-IF NIT EXCEEDED MIT. ITERM=12-IF NFV EXCEEDED MFV.
!         ITERM=13-IF NFG EXCEEDED MFG. ITERM<0-IF THE METHOD FAILED.
!
! VARIABLES IN COMMON /STAT/ (STATISTICS) :
!  IO  NRES  NUMBER OF RESTARTS.
!  IO  NDEC  NUMBER OF MATRIX DECOMPOSITIONS.
!  IO  NIN  NUMBER OF INNER ITERATIONS.
!  IO  NIT  NUMBER OF ITERATIONS.
!  IO  NFV  NUMBER OF FUNCTION EVALUATIONS.
!  IO  NFG  NUMBER OF GRADIENT EVALUATIONS.
!  IO  NFH  NUMBER OF HESSIAN EVALUATIONS.
!
! SUBPROGRAMS USED :
!  S   PEQN  SOLUTION OF SPARSE NONLINEAR SYSTEMS OF EQUATIONS BY THE
!         NEWTON METHOD USING THE PRECONDITIONED SMOOTHED CGS METHOD
!         FOR ITERATIVE SOLUTION OF THE LINEARIZED SYSTEM.
!
! EXTERNAL SUBROUTINES :
!  SE  FUN  COMPUTATION OF THE VALUE OF THE APPROXIMATED FUNCTION.
!         CALLING SEQUENCE: CALL FUN(N,KA,X,FA) WHERE N IS A NUMBER
!         OF VARIALES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(N) IS A VECTOR OF VARIABLES AND FA IS THE VALUE OF THE
!         APPROXIMATED FUNCTION.
!
      SUBROUTINE PEQNU (N, MA, X, AF, IAG, JAG, IPAR, RPAR, F, GMAX,
     &IDER, ISPAS, IPRNT, ITERM)
      DOUBLE PRECISION F,GMAX
      INTEGER IDER,ISPAS,IPRNT,ITERM,MA,N
      DOUBLE PRECISION AF(*),RPAR(9),X(*)
      INTEGER IAG(*),IPAR(7),JAG(*)
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      INTEGER LAFD,LAFO,LAG,LG,LGA,LGO,LGP,LGS,LIB,LIW1,LIW2,LIW3,LIW4,
     &LS,LXO,LXP,LXS,IER
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      INTEGER IA(:)
      DOUBLE PRECISION RA(:)
      ALLOCATABLE IA,RA
      IF (ISPAS.LE.1) THEN
        CALL PASED3 (N, N, MA, IAG, JAG, IER)
        IF (IER.NE.0) THEN
          WRITE (6,'(''INPUT ERROR : IER = '',I3)') IER
          STOP
        END IF
      ELSE
        MA=IAG(N+1)-1
      END IF
      ALLOCATE(IA(5*N),RA(11*N+2*MA))
!
!     POINTERS FOR AUXILIARY ARRAYS
!
      LGA=1
      LAG=LGA+N
      LG=LAG+2*MA
      LS=LG+N
      LXO=LS+N
      LGO=LXO+N
      LXS=LGO+N
      LGS=LXS+N
      LXP=LGS+N
      LGP=LXP+N
      LAFO=LGP+N
      LAFD=LAFO+N
      LIB=1
      LIW1=LIB+N
      LIW2=LIW1+N
      LIW3=LIW2+N
      LIW4=LIW3+N
      CALL PEQN (N, X, RA(LGA), RA(LAG), IAG, JAG, IA(LIB), IA(LIW1),
     &IA(LIW2), IA(LIW3), IA(LIW4), RA(LG), RA(LS), RA(LXO), RA(LGO),
     &RA(LXS), RA(LGS), RA(LXP), RA(LGP), AF, RA(LAFO), RA(LAFD),
     &RPAR(1), RPAR(2), RPAR(3), RPAR(4), RPAR(5), RPAR(8), GMAX, F,
     &IPAR(1), IPAR(2), IPAR(5), IPAR(6), IDER, IPRNT, ITERM)
      DEALLOCATE(IA,RA)
      RETURN
      END
! SUBROUTINE PEQN                   ALL SYSTEMS                95/12/01
! PORTABILITY : ALL SYSTEMS
! 95/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
! SOLUTION OF SPARSE NONLINEAR SYSTEMS OF EQUATIONS BY THE NEWTON
! METHOD USING THE PRECONDITIONED SMOOTHED CGS SUBALGORITHM FOR
! ITERATIVE SOLUTION OF LINEARIZED SYSTEMS.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RA  GA(N)  GRADIENT OF THE APPROXIMATED FUNCTION.
!  RA  AG(IAG(N+1)-1)  SPARSE RECTANGULAR MATRIX WHICH IS USED FOR THE
!         DIRECTION VECTOR DETERMINATION.
!  II  IAG(N+1)  POSITION OF THE FIRST ROWS ELEMENTS IN THE FIELD AG.
!  II  JAG(IAG(N+1)-1) COLUMN INDICES OF ELEMENTS IN THE FIELD AG.
!  IA  IB(N)  PERMUTATION VECTOR.
!  IA  IW1(N)  AUXILIARY VECTOR.
!  IA  IW2(N)  AUXILIARY VECTOR.
!  IA  IW3(N)  AUXILIARY VECTOR.
!  IA  IW4(N)  AUXILIARY VECTOR.
!  RA  G(N)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  S(N)  DIRECTION VECTOR.
!  RA  XO(N)  AUXILIARY VECTOR.
!  RA  GO(N)  AUXILIARY VECTOR.
!  RA  XS(N)  AUXILIARY VECTOR.
!  RA  GS(N)  AUXILIARY VECTOR.
!  RA  XP(N)  AUXILIARY VECTOR.
!  RA  GP(N)  AUXILIARY VECTOR.
!  RO  AF(N)  VECTOR WHOSE ELEMENTS ARE VALUES OF THE APPROXIMATED
!         FUNCTIONS.
!  RA  AFO(N)  AUXILIARY VECTOR.
!  RA  AFD(N)  AUXILIARY VECTOR.
!  RI  XMAX  MAXIMUM STEPSIZE.
!  RI  TOLX  TOLERANCE FOR CHANGE OF VARIABLES.
!  RI  TOLF  TOLERANCE FOR CHANGE OF FUNCTION VALUES.
!  RI  TOLB  TOLERANCE FOR THE FUNCTION VALUE.
!  RI  TOLG  TOLERANCE FOR THE GRADIENT OF THE LAGRANGIAN FUNCTION.
!  RI  ETA2  DAMPING PARAMETER FOR AN INCOMPLETE LU PRECONDITIONER.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  II  MIT  MAXIMUM NUMBER OF ITERATIONS.
!  II  MFV  MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
!  II  MOS1  CHOICE OF SMOOTHING STRATEGY FOR THE CGS METHOD.
!         MOS1=1-NO SMOOTHING. MOS1=2-SINGLE SMOOTHING STRATEGY
!         IS USED. MOS1=3-DOUBLE SMOOTHING STRATEGY IS USED.
!  II  MOS2  TYPE OF PRECONDITIONING. MOS2=1-PRECONDITIONING IS NOT
!         USED. MOS2=2-PRECONDITIONING BY THE INCOMPLETE GILL-MURRAY
!         DECOMPOSITION. MOS2=3-PRECONDITIONING BY THE INCOMPLETE
!         GILL-MURRAY DECOMPOSITION WITH A PRELIMINARY SOLUTION OF
!         THE PRECONDITIONED SYSTEM WHICH IS USED IF IT SATISFIES
!         THE TERMINATION CRITERION.
!  II  IDER  DEGREE OF ANALYTICALLY COMPUTED DERIVATIVES (0 OR 1).
!  II  IPRNT  PRINT SPECIFICATION. IPRNT=0-NO PRINT.
!         ABS(IPRNT)=1-PRINT OF FINAL RESULTS.
!         ABS(IPRNT)=2-PRINT OF FINAL RESULTS AND ITERATIONS.
!         IPRNT>0-BASIC FINAL RESULTS. IPRNT<0-EXTENDED FINAL
!         RESULTS.
!  IO  ITERM  VARIABLE THAT INDICATES THE CAUSE OF TERMINATION.
!         ITERM=1-IF ABS(X-XO) WAS LESS THAN OR EQUAL TO TOLX IN
!                   MTESX (USUALLY TWO) SUBSEQUENT ITERATIONS.
!         ITERM=2-IF ABS(F-FO) WAS LESS THAN OR EQUAL TO TOLF IN
!                   MTESF (USUALLY TWO) SUBSEQUENT ITERATIONS.
!         ITERM=3-IF F IS LESS THAN OR EQUAL TO TOLB.
!         ITERM=4-IF GMAX IS LESS THAN OR EQUAL TO TOLG.
!         ITERM=6-IF THE TERMINATION CRITERION WAS NOT SATISFIED,
!                   BUT THE SOLUTION OBTAINED IS PROBABLY ACCEPTABLE.
!         ITERM=11-IF NIT EXCEEDED MIT. ITERM=12-IF NFV EXCEEDED MFV.
!         ITERM=13-IF NFG EXCEEDED MFG. ITERM<0-IF THE METHOD FAILED.
!
! VARIABLES IN COMMON /STAT/ (STATISTICS) :
!  IO  NRES  NUMBER OF RESTARTS.
!  IO  NDEC  NUMBER OF MATRIX DECOMPOSITIONS.
!  IO  NIN  NUMBER OF INNER ITERATIONS.
!  IO  NIT  NUMBER OF ITERATIONS.
!  IO  NFV  NUMBER OF FUNCTION EVALUATIONS.
!  IO  NFG  NUMBER OF GRADIENT EVALUATIONS.
!  IO  NFH  NUMBER OF HESSIAN EVALUATIONS.
!
! SUBPROGRAMS USED :
!  S   PA0SQ3  COMPUTATION OF THE VALUE AND THE GRADIENT OF THE
!         OBJECTIVE FUNCTION WHICH IS DEFINED AS A SUM OF SQUARES
!         OF THE APPROXIMATED FUNCTIONS (THE SPARSE CASE).
!  S   PS0L02  LINE SEARCH USING ONLY FUNCTION VALUES.
!  S   PYFUT1  TEST ON TERMINATION.
!  S   MXDPGB  BACK SUBSTITUTION USING THE GILL-MURRAY DECOMPOSITION
!         OBTAINED BY MXDPGF.
!  S   MXDPGF  GILL-MURRAY DECOMPOSITION OF A DENSE SYMMETRIC MATRIX.
!  S   MXSCMM  MATRIX-VECTOR PRODUCT. SPARSE RECTANGULAR MATRIX IS
!         STORED COLUMNWISE.
!  S   MXSGIB  BACK SUBSTITUTION USING THE INCOMPLETE LU DECOMPOSITION
!         OBTAINED BY MXSGIF.
!  S   MXSGIF  INCOMPLETE LU DECOMPOSITION OF A SPARSE NONSYMMETRIC
!         MATRIX.
!  S   MXSRMD  MATRIX-VECTOR PRODUCT FOLLOWED BY THE ADDITION OF A
!         SCALED VECTOR. SPARSE RECTANGULAR MATRIX IS STORED ROWWISE.
!  S   MXSRMM  MATRIX-VECTOR PRODUCT. SPARSE RECTANGULAR MATRIX IS
!         STORED ROWWISE.
!  S   MXSRSP  ROW PERMUTATIONS FOR OBTAINING DIAGONAL NONZEROS.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIF  DIFFERENCE OF TWO VECTORS.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  RF  MXVMAX  L-INFINITY NORM OF A VECTOR
!  S   MXVNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
!  RF  MXVNOR  EUCLIDEAN NORM OF A VECTOR.
!  S   MXVSAV  DIFFERENCE OF TWO VECTORS WITH COPYING AND SAVING THE
!         SUBSTRACTED ONE.
!  S   MXVSET  INITIATION OF A VECTOR.
!  S   MXVSUM  SUM OF TWO VECTORS.
!
! EXTERNAL SUBROUTINES :
!  SE  FUN  COMPUTATION OF THE VALUE OF THE APPROXIMATED FUNCTION.
!         CALLING SEQUENCE: CALL FUN(N,KA,X,FA) WHERE N IS A NUMBER
!         OF VARIABLES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(N) IS A VECTOR OF VARIABLES AND FA IS THE VALUE OF THE
!         APPROXIMATED FUNCTION.
!
! METHOD :
! PRECONDITIONED SMOOTHED CGS METHOD WITH INEXACT TERMINATION.
!
      SUBROUTINE PEQN (N, X, GA, AG, IAG, JAG, IB, IW1, IW2, IW3, IW4,
     &G, S, XO, GO, XS, GS, XP, GP, AF, AFO, AFD, XMAX, TOLX, TOLF,
     &TOLB, TOLG, ETA2, GMAX, F, MIT, MFV, MOS1, MOS2, IDER, IPRNT,
     &ITERM)
      DOUBLE PRECISION ETA2,F,GMAX,TOLB,TOLD,TOLF,TOLG,TOLS,TOLX,XMAX
      INTEGER IDER,IPRNT,ITERM,MES,MFV,MIT,MOS,MOS1,MOS2,N
      DOUBLE PRECISION AF(*),AFD(*),AFO(*),AG(*),G(*),GA(*),GO(*),GP(*),
     &GS(*),S(*),X(*),XO(*),XP(*),XS(*)
      INTEGER IAG(*),IB(*),IW1(*),IW2(*),IW3(*),IW4(*),JAG(*)
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      DOUBLE PRECISION ALF1,ALF2,DMAX,EPS6,ETA0,ETA9,FMAX,FMIN,FO,FP,
     &GNORM,P,PO,PP,R,RMAX,RMIN,RO,RP,SNORM,BTB(3),BTR(2),RMU,RNU,ALF,
     &BET,SIG,RHO,RHO1,RHO2,PAR,UMAX
      INTEGER I,INF,IPOM1,IPOM2,IREST,ITERD,ITERS,NRED,KD,KIT,LD,MA,
     &MRED,MTESF,MTESX,NTESF,NTESX,LDS,IDECA,INITS,KTERS,IEST,ITES,
     &MAXST,IRES1,IRES2,ISYS,MFG
      DOUBLE PRECISION MXVDOT,MXVNOR,MXVMAX
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      IF (ABS(IPRNT).GT.1) WRITE (6,'(1X,''ENTRY TO PEQN :'')')
!
!     INITIATION
!
      NRES=0
      NDEC=0
      NIN=0
      NIT=0
      NFV=0
      NFG=0
      NFH=0
      ISYS=0
      IEST=1
      ITES=1
      MTESX=2
      MTESF=2
      INITS=1
      ITERM=0
      ITERD=0
      ITERS=2
      KTERS=5
      IREST=1
      IRES1=999
      IRES2=0
      MRED=20
      IDECA=0
      IPOM1=0
      IPOM2=0
      ETA0=1.0D-15
      IF (ETA2.LE.0.0D0.OR.ETA2.GE.1.0D0) ETA2=0.0D0
      ETA9=1.0D120
      EPS6=2.5D-1
      ALF1=1.0D-15
      ALF2=1.0D10
      FMAX=1.0D60
      FMIN=0.0D0
      IF (XMAX.LE.0.0D0) XMAX=1.0D16
      IF (TOLX.LE.0.0D0) TOLX=1.0D-16
      IF (TOLF.LE.0.0D0) TOLF=1.0D-16
      IF (TOLB.LE.0.0D0) TOLB=1.0D-16
      IF (TOLG.LE.0.0D0) TOLG=1.0D-6
      TOLD=1.0D-12
      TOLS=1.0D-4
      MES=1
      MOS=1
      IF (MOS1.LE.0) MOS1=3
      IF (MOS2.EQ.0) MOS2=3
      IDER=MAX(IDER,0)
      IF (MIT.LE.0) MIT=1000
      IF (MFV.LE.0) MFV=1000
      MFG=MFV
      KD=0
      LD=-1
      KIT=0
      FO=FMIN
      GMAX=ETA9
      DMAX=ETA9
!
!     SYMBOLIC PREPATION OF INCOMPLETE LU DECOMPOSITION
!
      MA=IAG(N+1)-1
      IF (MOS2.GT.1) CALL MXSRSP (N, IAG, JAG, IB, INF, IW1, IW2, IW3,
     &IW4)
!
!     COMPUTATION OF THE VALUE OF THE OBJECTIVE FUNCTION
!
      CALL PA0SQ3 (N, X, F, AF, GA, AG, IAG, JAG, G, ETA0, KD, LD, NFV,
     &NFG, IDER)
   10 IF (ABS(IPRNT).GT.1) WRITE (6,'(1X,''NIT='',I5,2X,''NFV='',I5,2X,'
     &'NFG='',I5,2X,       ''F='', G13.6,2X,''G='',G13.6)') NIT,NFV,NFG,
     &F,GNORM
!
!     START OF THE ITERATION WITH TESTS FOR TERMINATION.
!
      CALL PYFUT1 (N, F, FO, UMAX, GMAX, DMAX, TOLX, TOLF, TOLB, TOLG,
     &KD, NIT, KIT, MIT, NFV, MFV, NFG, MFG, NTESX, MTESX, NTESF, MTESF,
     & ITES, IRES1, IRES2, IREST, ITERS, ITERM)
      IF (ITERM.NE.0) GO TO 100
   20 IF (IREST.LE.0) GO TO 30
!
!     RESTART
!
      KD=1
      CALL PA0SQ3 (N, X, F, AF, GA, AG, IAG, JAG, G, ETA0, KD, LD, NFV,
     &NFG, IDER)
      IDECA=0
      IF (KIT.LT.NIT) THEN
        NRES=NRES+1
        KIT=NIT
      ELSE
        ITERM=-10
        IF (ITERS.LT.0) ITERM=ITERS-5
      END IF
   30 CONTINUE
      IF (ITERM.NE.0) GO TO 100
!
!     DIRECTION DETERMINATION USING PRECONDITIONED SMOOTHED CGS
!     ALGORITHM
!
      IF (IDECA.LT.0) IDECA=0
      IF (IDECA.EQ.2) THEN
      ELSE IF (IDECA.NE.0) THEN
        ITERD=-1
        GO TO 60
      ELSE IF (MOS2.GT.1) THEN
!
!     CONSTRUCTION OF PRECONDITIONER
!
        INF=0
        CALL MXVCOP (MA, AG, AG(MA+1))
        CALL MXSGIF (N, AG(MA+1), IAG, JAG, IB, IW1, IW2, ETA2, INF)
        IF (INF.LT.0) THEN
          ITERD=INF
          GO TO 60
        ELSE
          NDEC=NDEC+1
          IDECA=2
        END IF
      END IF
      IF (MOS.EQ.1) THEN
        IF (LD.LE.0) CALL MXSCMM (N, N, AG, IAG, JAG, AF, G)
      ELSE
        CALL MXVCOP (N, AF, G)
      END IF
      GNORM=SQRT(MXVDOT(N,G,G))
      PAR=SQRT(F/FO)**1.618D0
      PAR=MAX(PAR,SQRT(SQRT(2.0D0*F)))
      PAR=MIN(EPS6,PAR)
      IF (PAR.GT.1.0D1*1.0D-3) THEN
        PAR=MIN(PAR,1.0D0/DBLE(NIT))
      END IF
      PAR=PAR*PAR
      RHO2=MXVDOT(N,AF,AF)
      IF (MOS2.GT.2) THEN
!
!     PRELIMINARY INEXACT SOLUTION
!
        CALL MXVNEG (N, AF, S)
        CALL MXSGIB (N, AG(MA+1), IAG, JAG, IB, IW1, S, XO, 0)
        CALL MXSRMD (N, AG, IAG, JAG, S, 1.0D0, AF, AFO)
        RHO1=MXVDOT(N,AFO,AFO)
        IF (RHO1.LE.PAR*RHO2) THEN
          SNORM=SQRT(MXVDOT(N,S,S))
          ITERD=1
          GO TO 50
        END IF
      END IF
      ITERD=2
!
!     CGS INITIATION
!
      SNORM=0.0D0
      RHO=1.0D0
      CALL MXVNEG (N, AF, AFO)
      CALL MXVNEG (N, AF, AFD)
      CALL MXVSET (N, 0.0D0, S)
      CALL MXVSET (N, 0.0D0, XO)
      CALL MXVSET (N, 0.0D0, GO)
      CALL MXVSET (N, 0.0D0, XS)
      SIG=MXVNOR(N,AFD)
      NRED=0
!
!    CGS ITERATIONS
!
      DO 40 NRED=1,2*N
        RHO1=RHO
        IF (RHO1.EQ.0.0D0) THEN
          ITERD=-4
          GO TO 60
        END IF
        RHO=MXVDOT(N,G,AFD)
        BET=RHO/RHO1
        CALL MXVDIR (N, BET, XS, AFD, GS)
        CALL MXVDIR (N, BET, GO, XS, GO)
        CALL MXVDIR (N, BET, GO, GS, GO)
!
!     CGS PRECONDITIONING
!
        CALL MXVCOP (N, GO, GA)
        IF (MOS2.GT.1) CALL MXSGIB (N, AG(MA+1), IAG, JAG, IB, IW1, GA,
     &   XP, 0)
        CALL MXSRMM (N, AG, IAG, JAG, GA, XP)
        SIG=MXVDOT(N,G,XP)
        IF (SIG.EQ.0.0D0) THEN
          ITERD=-5
          GO TO 60
        END IF
        ALF=RHO/SIG
!
!     CGS STEP
!
        CALL MXVDIR (N, -ALF, XP, GS, XS)
        CALL MXVSUM (N, GS, XS, GS)
!
!     CGS PRECONDITIONING
!
        IF (MOS2.GT.1) CALL MXSGIB (N, AG(MA+1), IAG, JAG, IB, IW1, GS,
     &   GP, 0)
        CALL MXSRMM (N, AG, IAG, JAG, GS, GP)
!
!     CGS STEP
!
        CALL MXVDIR (N, -ALF, GP, AFD, AFD)
        CALL MXVDIR (N, ALF, GS, XO, XO)
        NIN=NIN+1
!
!     CGS SMOOTHING
!
        IF (MOS1.EQ.1) THEN
          CALL MXVCOP (N, AFD, AFO)
          CALL MXVCOP (N, XO, S)
        ELSE
          RMU=ETA0**2
          CALL MXVDIF (N, AFO, AFD, GP)
          BTB(1)=MXVDOT(N,GP,GP)
          BTR(1)=MXVDOT(N,GP,AFD)
          IF (MOS1.EQ.3) THEN
            BTB(2)=MXVDOT(N,GP,XP)
            BTB(3)=MXVDOT(N,XP,XP)
            BTR(2)=MXVDOT(N,XP,AFD)
            CALL MXDPGF (2, BTB, INF, RMU, RNU)
            CALL MXDPGB (2, BTB, BTR, 0)
            RMU=-BTR(1)
            RNU=-BTR(2)
          ELSE
            RMU=-BTR(1)/MAX(BTB(1),RMU)
          END IF
          CALL MXVDIR (N, RMU, GP, AFD, AFO)
          CALL MXVDIF (N, S, XO, GP)
          CALL MXVDIR (N, RMU, GP, XO, S)
          IF (MOS1.EQ.3) THEN
            CALL MXVDIR (N, RNU, XP, AFO, AFO)
            CALL MXVDIR (N, -RNU, GA, S, S)
          END IF
        END IF
        SNORM=MXVNOR(N,S)
        IF (SNORM.GE.XMAX) GO TO 50
        RHO1=MXVDOT(N,AFO,AFO)
        IF (RHO1.LE.PAR*RHO2) GO TO 50
   40 CONTINUE
!
!     AN INEXACT SOLUTION IS OBTAINED
!
   50 CONTINUE
      P=-F
   60 CONTINUE
      IF (ITERD.LT.0) THEN
        ITERM=ITERD
      ELSE
!
!     TEST FOR SUFFICIENT DESCENT
!
        IF (SNORM.LE.0.0D0) THEN
          IREST=MAX(IREST,1)
        ELSE IF (P+TOLD*GNORM*SNORM.LE.0.0D0) THEN
          IREST=0
        ELSE
!
!     UNIFORM DESCENT CRITERION
!
          IREST=MAX(IREST,1)
        END IF
        IF (IREST.EQ.0) THEN
!
!     PREPARATION OF LINE SEARCH
!
          NRED=0
          RMIN=ALF1*GNORM/SNORM
          RMAX=MIN(ALF2*GNORM/SNORM,XMAX/SNORM)
        END IF
      END IF
      IF (ITERM.NE.0) GO TO 100
      IF (IREST.NE.0) GO TO 20
      LDS=LD
      FP=FO
      RO=0.0D0
      FO=F
      PO=P
      CALL MXVCOP (N, X, XO)
      CALL MXVCOP (N, AF, AFO)
   70 CALL PS0L02 (R, RO, RP, F, FO, FP, PO, PP, FMIN, FMAX, RMIN, RMAX,
     & TOLS, KD, LD, NIT, KIT, NRED, MRED, MAXST, IEST, INITS, ITERS,
     &KTERS, MES, ISYS)
      IF (ISYS.EQ.0) GO TO 80
      CALL MXVDIR (N, R, S, XO, X)
      CALL PA0SQ3 (N, X, F, AF, GA, AG, IAG, JAG, G, ETA0, KD, LD, NFV,
     &NFG, IDER)
      GO TO 70
   80 CONTINUE
      IF (ITERS.LE.0) THEN
        R=0.0D0
        F=FO
        P=PO
        CALL MXVCOP (N, XO, X)
        CALL MXVCOP (N, AFO, AF)
        IREST=MAX(IREST,1)
        LD=LDS
        GO TO 20
      END IF
      IF (IPOM1.EQ.1) KD=1
      IF (KD.GT.LD) THEN
        CALL PA0SQ3 (N, X, F, AF, GA, AG, IAG, JAG, G, ETA0, KD, LD,
     &   NFV, NFG, IDER)
      END IF
      KD=0
      IREST=1
      IF (ITERS.GT.0) THEN
        CALL MXVDIF (N, X, XO, XO)
        PO=R*PO
        P=R*P
      ELSE
        F=FO
        P=PO
        CALL MXVSAV (N, X, XO)
        LD=KD
      END IF
      DMAX=0.0D0
      DO 90 I=1,N
        DMAX=MAX(DMAX,ABS(XO(I))/MAX(ABS(X(I)),1.0D0))
   90 CONTINUE
      GO TO 10
  100 CONTINUE
      IF (IPRNT.GT.1.OR.IPRNT.LT.0) WRITE (6,'(1X,''EXIT FROM PEQN :'')'
     &)
      IF (IPRNT.NE.0) THEN
        GMAX=MXVMAX(N,G)
        WRITE (6,'(1X,''NIT='',I5,2X,''NFV='',I5,2X,''NFG='',I5,2X,
     &   ''F='', G13.6,2X,''G='',G13.6,2X,''ITERM='',I3)') NIT,NFV,NFG,
     &   F,GMAX,ITERM
      END IF
      IF (IPRNT.LT.0) WRITE (6,'(1X,''X='',5(G14.7,1X):/(3X,5(G14.7,1X))
     &)') (X(I),I=1,N)
      RETURN
      END
