! SUBROUTINE PSENU              ALL SYSTEMS                   97/01/22
! PURPOSE :
! EASY TO USE SUBROUTINE FOR LARGE-SCALE UNCONSTRAINED MINIMIZATION
! OF NONSMOOTH PARTIALLY SEPARABLE FUNCTIONS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NA  NUMBER OF PARTIAL FUNCTIONS.
!  IU  MA  NUMBER OF NONZERO ELEMENTS IN THE JACOBIAN MATRIX.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  RO  AF(NA)  VECTOR CONTAINING VALUES OF THE APPROXIMATED
!         FUNCTIONS.
!  RI  IAG(NA+1)  POSITION OF THE FIRST ROWS ELEMENTS IN THE JACOBIAN
!         MATRIX.
!  RI  JAG(MA)  COLUMN INDICES OF ELEMENTS IN THE JACOBIAN MATRIX.
!  II  IPAR(7)  INTEGER PAREMETERS:
!      IPAR(1)  MAXIMUM NUMBER OF ITERATIONS.
!      IPAR(2)  MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
!      IPAR(3)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PMAX.
!      IPAR(4)  ESTIMATION INDICATOR. IPAR(4)=0-MINIMUM IS NOT
!         ESTIMATED. IPAR(4)=1-MINIMUM IS ESTIMATED BY THE VALUE
!         RPAR(6).
!      IPAR(5)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PSEN.
!      IPAR(6)  DIMENSION OF A BUNDLE USED IN THE LINE SEARCH.
!      IPAR(7)  NUMBER DEFINING THE SPACE FOR FILL-IN (THE SIZE OF THIS
!         SPACE IS IFIL TIMES THE STANDARD SIZE). THE DEFAULT VALUE IS
!         IFIL=1. THE DEFAULT VALUE HAS TO BE INCREASED IF ITERM IS
!         LESS OR EQUAL TO -40.
!  RI  RPAR(9)  REAL PARAMETERS:
!      RPAR(1)  MAXIMUM STEPSIZE.
!      RPAR(2)  TOLERANCE FOR THE CHANGE OF VARIABLES.
!      RPAR(3)  TOLERANCE FOR THE CHANGE OF FUNCTION VALUES.
!      RPAR(4)  TOLERANCE FOR THE FUNCTION FALUE.
!      RPAR(5)  TOLERANCE FOR THE GRADIENT NORM.
!      RPAR(6)  ESTIMATION OF THE MINIMUM FUNCTION VALUE.
!      RPAR(7)  THIS PARAMETER IS NOT USED IN THE SUBROUTINE PSEN.
!      RPAR(8)  CORRECTION PARAMETER.
!      RPAR(9)  PARAMETER FOR SUBGRADIENT LOCALITY MEASURE.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE.
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
!         VALUES ITERM<=-40 DETECT A LACK OF SPACE. IN THIS CASE,
!         PARAMETER IPAR(7) HAS TO BE INCREASED.
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
!  S   PSEN  BUNDLE VARIABLE METRIC METHOD FOR NONSMOOTH PARTIALLY
!         SEPARABLE FUNCTION.
!  S   PASED3  COMPRESSED SPARSE STRUCTURE OF THE JACOBIAN MATRIX IS
!         COMPUTED FROM THE COORDINATE FORM.
!  S   PFSET2  NUMBER OF NONZERO ELEMENTS IN THE PARTITIONED HESSIAN
!         MATRIX.
!
! EXTERNAL SUBROUTINES :
!  SE  FUN  COMPUTATION OF THE VALUE OF THE PARTIAL FUNCTION.
!         CALLING SEQUENCE: CALL FUN(NF,KA,X,FA) WHERE NF IS A NUMBER
!         OF VARIABLES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(NF) IS A VECTOR OF VARIABLES AND FA IS THE VALUE OF THE
!         APPROXIMATED FUNCTION.
!  SE  DFUN  COMPUTATION OF THE SUBGRADIENT OF THE PARTIAL FUNCTION.
!         CALLING SEQUENCE: CALL DFUN(NF,KA,X,GA) WHERE NF IS A NUMBER
!         OF VARIABLES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(NF) IS A VECTOR OF VARIABLES AND GA(NF) IS THE SUBGRADIENT
!         OF THE APPROXIMATED FUNCTION.
!
      SUBROUTINE PSENU (NF, NA, MA, X, AF, IAG, JAG, IPAR, RPAR, F,
     &GMAX, ISPAS, IPRNT, ITERM)
      INTEGER NF,NA,MA,IAG(*),JAG(*),IPAR(7),ISPAS,IPRNT,ITERM
      DOUBLE PRECISION X(*),AF(*),RPAR(9),F,GMAX
      INTEGER LAG,LAGO,LAH,LGA,LG,LH,LS,LXO,LGO,LXS,LGS,LGP,LIH,LJH,LAX,
     &LAY,LAZ,LPSL,LPERM,LINVP,LWN11,LWN12,LWN13,LWN14,MB,MC,MH,IFIL,
     &IER
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      INTEGER IA(:)
      DOUBLE PRECISION RA(:)
      ALLOCATABLE IA,RA
      MB=IPAR(6)
      IF (MB.LE.0) MB=20
      IFIL=IPAR(7)
      IF (IFIL.LE.0) IFIL=1
      IF (ISPAS.LE.1) THEN
        CALL PASED3 (NF, NA, MA, IAG, JAG, IER)
        IF (IER.NE.0) THEN
          WRITE (6,'(''INPUT ERROR : IER = '',I3)') IER
          STOP
        END IF
      ELSE
        MA=IAG(NA+1)-1
      END IF
      CALL PFSET2 (NA, MH, MC, IAG)
      ALLOCATE(IA(8*NF+6+(IFIL+3)*MH),RA(2*MA+8*NF+2*(NF+2)*MB+(IFIL+4)*
     &MH))
!
!     POINTERS FOR AUXILIARY ARRAYS
!
      LAG=1
      LAGO=LAG+MA
      LAH=LAGO+MA
      LGA=LAH+MH
      LG=LGA+NF
      LS=LG+NF
      LXO=LS+NF
      LGO=LXO+NF
      LXS=LGO+NF
      LGS=LXS+NF
      LGP=LGS+NF
      LAX=LGP+NF
      LAY=LAX+NF*MB
      LAZ=LAY+NF*MB
      LH=LAZ+4*MB
      LPSL=1
      LPERM=LPSL+NF+1
      LINVP=LPERM+NF
      LWN11=LINVP+NF
      LWN12=LWN11+NF+1
      LWN13=LWN12+NF+1
      LWN14=LWN13+NF+1
      LIH=LWN14+NF+1
      LJH=LIH+NF+1
      CALL PSEN (NF, NA, MB, (IFIL+3)*MH, X, IA, AF, RA(LAG), RA(LAGO),
     &RA(LAH), RA(LGA), RA(LG), RA(LH), IA(LIH), IA(LJH), IAG, JAG,
     &RA(LS), RA(LXO), RA(LGO), RA(LXS), RA(LGS), RA(LGP), RA(LAX),
     &RA(LAY), RA(LAZ), IA(LPSL), IA(LPERM), IA(LINVP), IA(LWN11),
     &IA(LWN12), IA(LWN13), IA(LWN14), RPAR(1), RPAR(2), RPAR(3),
     &RPAR(4), RPAR(5), RPAR(6), RPAR(8), RPAR(9), GMAX, F, IPAR(1),
     &IPAR(2), IPAR(4), IPRNT, ITERM)
      DEALLOCATE(IA,RA)
      RETURN
      END
! SUBROUTINE PSEN               ALL SYSTEMS                   01/09/22
! PURPOSE :
! GENERAL SUBROUTINE FOR LARGE-SCALE UNCONSTRAINED MINIMIZATION
! OF NONSMOOTH PARTIALLY SEPARABLE FUNCTIONS.
!
! PARAMETERS :
!  II  NF  NUMBER OF VARIABLES.
!  II  NA  NUMBER OF LINEAR APPROXIMATED FUNCTIONS.
!  II  MB  DIMENSION OF A BUNDLE USED IN THE LINE SEARCH.
!  II  MMAX  MAXIMUM DIMENSION OF THE SPARSE TABLEAU.
!  RI  X(NF)  VECTOR OF VARIABLES.
!  IA  IX(NF)  AUXILIARY VECTOR.
!  RI  AF(NA)  VECTOR CONTAINING VALUES OF THE APPROXIMATED
!         FUNCTIONS.
!  RA  AG(MA)  JACOBIAN MATRIX OF THE PARTITIONED FUNCTION.
!  RA  AGO(MA)  OLD JACOBIAN MATRIX OF THE PARTITIONED FUNCTION,
!  RA  AH(MB)  ELEMENTS OF THE PARTITIONED HESSIAN MATRIX.
!  RA  GA(NF)  GRADIENT OF THE SELECTED APPROXIMATED FUNCTION.
!  RO  G(NF)  GRADIENT OF THE OBJECTIVE FUNCTION.
!  RA  H(MMAX)  NONZERO ELEMENTS OF THE APPROXIMATION OF THE SPARSE
!         HESSIAN MATRIX TOGETHER WITH AN ADDITIONAL SPACE USED FOR
!         THE NUMERICAL DIFFERENTIATION.
!  II  IH(NF+1)  POINTERS OF DIAGONAL ELEMENTS OF THE MATRIX H.
!  IU  JH(MMAX)  INDICES OF NONZERO ELEMENTS OF THE MATRIX H
!         TOGETHER WITH AN ADDITIONAL SPACE USED FOR THE NUMERICAL
!         DIFFERENTIATION.
!  RI  IAG(NA+1)  POSITION OF THE FIRST ROWS ELEMENTS IN THE FIELD AG.
!  RI  JAG(MA)  COLUMN INDICES OF ELEMENTS IN THE FIELD AG.
!  RA  S(NF)  DIRECTION VECTOR.
!  RA  XO(NF)  VECTORS OF VARIABLES DIFFERENCE.
!  RA  GO(NF)  GRADIENTS DIFFERENCE.
!  RA  XS(NF)  AUXILIARY VECTOR.
!  RA  GS(NF)  AUXILIARY VECTOR.
!  RA  GP(NF)  AUXILIARY VECTOR.
!  RA  AX(NF*MB)  AUXILIARY VECTOR.
!  RA  AY(NF*MB)  AUXILIARY VECTOR.
!  RA  AZ(4*MB)  AUXILIARY VECTOR.
!  IA  PSL(NF+1)  POINTER VECTOR OF THE COMPACT FORM OF THE TRIANGULAR
!         FACTOR OF THE HESSIAN APPROXIMATION.
!  IA  PERM(NF)  PERMUTATION VECTOR.
!  IA  INVP(NF)  INVERSE PERMUTATION VECTOR.
!  IA  WN11(NF+1) AUXILIARY VECTOR.
!  IA  WN12(NF+1) AUXILIARY VECTOR.
!  IA  WN13(NF+1) AUXILIARY VECTOR.
!  IA  WN14(NF+1) AUXILIARY VECTOR.
!  RI  XMAX  MAXIMUM STEPSIZE.
!  RI  TOLX  TOLERANCE FOR CHANGE OF VARIABLES.
!  RI  TOLF  TOLERANCE FOR CHANGE OF FUNCTION VALUES.
!  RI  TOLB  TOLERANCE FOR THE FUNCTION VALUE.
!  RI  TOLG  TOLERANCE FOR THE GRADIENT NORM.
!  RI  FMIN  ESTIMATION OF THE MINIMUM FUNCTION VALUE.
!  RI  ETA3  CORRECTION PARAMETER.
!  RI  ETA5  PARAMETER FOR SUBGRADIENT LOCALITY MEASURE.
!  RO  GMAX  MAXIMUM PARTIAL DERIVATIVE.
!  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
!  II  MIT  MAXIMUM NUMBER OF ITERATIONS.
!  II  MFV  MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
!  II  IEST  ESTIMATION INDICATOR. IEST=0-MINIMUM IS NOT ESTIMATED.
!         IEST=1-MINIMUM IS ESTIMATED BY THE VALUE FMIN.
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
!         VALUES ITERM<=-40 DETECT A LACK OF SPACE.
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
!  S   PA1SF3  COMPUTATION OF THE VALUE AND THE SUBGRADIENT OF A
!         PARTIALLY SEPARABLE OBJECTIVE FUNCTION.
!  S   PFSET3  PREPARATION OF THE SPARSE HESSIAN MATRIX
!  S   PFSET4  PREPARATION OF THE PARTITIONED HESSIAN MATRIX
!  S   PS1L18  SPECIAL NONSMOOTH LINE SEARCH.
!  S   PUBVI2  VARIABLE METRIC UPDATE OF THE PARTITIONED HESSIAN MATRIX.
!  S   PYABU1  DETERMINATION OF THE AGGREGATE SUBGRADIENT AS A SOLUTION
!         OF THE THREE-TERM QUADRATIC PROGRAMMING SUBPROBLEM.
!  S   PYABU2  DETERMINATION OF THE AGGREGATE SUBGRADIENT AS A SOLUTION
!         OF THE TWO-TERM QUADRATIC PROGRAMMING SUBPROBLEM.
!  S   PYBUN1  BUNDLE SELECTION.
!  S   PYTRCD  COMPUTATION OF PROJECTED DIFFERENCES FOR THE VARIABLE MET
!         UPDATE.
!  S   PYTRCG  COMPUTATION OF THE PROJECTED GRADIENT.
!  S   PYTRCS  COMPUTATION OF THE PROJECTED DIRECTION VECTOR.
!  S   PYTSCH  CORRECTION OF THE HESSIAN MATRIX.
!  S   MXBSMI  INITIATION OF THE PARTITIONED HESSIAN MATRIX.
!  S   MXSPCB  BACK SUBSTITUTION USING THE SPARSE DECOMPOSITION
!         OBTAINED BY MXSPCF.
!  S   MXSPCC  SPARSE MATRIX REORDERING, SYMBOLIC FACTORIZATION, DATA
!         STRUCTURES TRANSFORMATION. INITIATION OF THE DIRECT SPARSE
!         SOLVER.
!  S   MXSPCF  GILL-MURRAY DECOMPOSITION OD A SPARSE SYMMETRIC MATRIX.
!  S   MXSPCT  COPYING A SPARSE SYMMETRIC MATRIX INTO THE PERMUTED
!         FACTORIZED COMPACT SCHEME.
!  RF  MXUDOT  DOT PRODUCT OF TWO VECTORS.
!  S   MXUNEG  COPYING OF A VECTOR WITH CHANGE OF THE SIGN.
!  S   MXVCOP  COPYING OF A VECTOR.
!  S   MXVDIR  VECTOR AUGMENTED BY THE SCALED VECTOR.
!  RF  MXVDOT  DOT PRODUCT OF TWO VECTORS.
!  S   MXVINE  RESTORATION OF A SPARSE SYMMETRIC MATRIX OBTAINED BY
!         MXVINB
!  S   MXVINS  INITIATION OF THE INTEGER VECTOR.
!  RF  MXVNOR  EUCLIDEAN NORM OF A VECTOR.
!  S   MXVSAB  SUM OF ABSOLUTE VALUES.
!  S   MXVSBP  INVERSE PERMUTATION OF A VECTOR
!  S   MXVSET  INITIATION OF A VECTOR.
!  S   MXVSFP  PERMUTATION OF A VECTOR.
!
! EXTERNAL SUBROUTINES :
!  SE  FUN  COMPUTATION OF THE VALUE OF THE APPROXIMATED FUNCTION.
!         CALLING SEQUENCE: CALL FUN(NF,KA,X,FA) WHERE NF IS A NUMBER
!         OF VARIABLES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(NF) IS A VECTOR OF VARIABLES AND FA IS THE VALUE OF THE
!         APPROXIMATED FUNCTION.
!  SE  DFUN  COMPUTATION OF THE SUBGRADIENT OF THE APPROXIMATED FUNCTION
!         CALLING SEQUENCE: CALL DFUN(NF,KA,X,GA) WHERE NF IS A NUMBER
!         OF VARIABLES, KA IS THE INDEX OF THE APPROXIMATED FUNCTION,
!         X(NF) IS A VECTOR OF VARIABLES AND GA(NF) IS THE GRADIENT OF
!         THE APPROXIMATED FUNCTION.
!
! METHOD :
! BUNDLE VARIABLE METRIC METHOD FOR MINIMIZATION OF NONSMOOTH PARTIALLY
! SEPARABLE FUNCTIONS.
!
      SUBROUTINE PSEN (NF, NA, MB, MMAX, X, IX, AF, AG, AGO, AH, GA, G,
     &H, IH, JH, IAG, JAG, S, XO, GO, XS, GS, GP, AX, AY, AZ, PSL, PERM,
     & INVP, WN11, WN12, WN13, WN14, XMAX, TOLX, TOLF, TOLB, TOLG, FMIN,
     & ETA3, ETA5, GMAX, F, MIT, MFV, IEST, IPRNT, ITERM)
      INTEGER NA,NF,MB,MMAX,IX(*),IH(*),JH(*),IAG(*),JAG(*),PSL(*),
     &PERM(*),INVP(*),WN11(*),WN12(*),WN13(*),WN14(*),MIT,MFV,MFG,IEST,
     &IPRNT,ITERM
      DOUBLE PRECISION X(*),AF(*),AG(*),AGO(*),AH(*),GA(*),G(*),H(*),S(*
     &),XO(*),GO(*),XS(*),GS(*),GP(*),AX(*),AY(*),AZ(*),XMAX,TOLX,TOLF,
     &TOLB,TOLG,FMIN,ETA3,ETA5,GMAX,F
      INTEGER IDECF,ITERD,ITERS,KD,LD,NTESX,NTESF,MTESX,MTESF,MRED,KIT,
     &IREST,KBF,IOLD,INITD,IER,N,NB,MA,ISYS,KTERS,IRES1,IRES2,NRED,I,
     &INITS,JC,JR,JE,NNK,ITERH,NNV,M,MM,MH,ISNA,MAM,MOS3
      DOUBLE PRECISION R,RO,RP,FO,FP,P,PO,GNORM,SNORM,XNORM,RMIN,RMAX,
     &FMAX,DMAX,UMAX,ETA0,ETA2,ETA9,EPS0,EPS1,EPS2,ALF1,ALF2,POM,ALFN,
     &ALFV,DF
      DOUBLE PRECISION MXVDOT,MXVNOR,MXVSAB
      INTEGER NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      COMMON /STAT/ NRES,NDEC,NIN,NIT,NFV,NFG,NFH
      INTEGER NEXT
      DOUBLE PRECISION FB
      COMMON /PROB/ FB,NEXT
      IF (ABS(IPRNT).GT.1) WRITE (6,'(1X,''ENTRY TO PSEN :'')')
      KBF=0
      NRES=0
      NDEC=0
      NIN=0
      NIT=0
      NFV=0
      NFG=0
      NFH=0
      ISYS=0
      NTESX=0
      NTESF=0
      MTESX=10
      MTESF=10
      INITS=1
      INITD=1
      ITERM=0
      ITERD=0
      ITERS=2
      ITERH=0
      KTERS=5
      IREST=0
      IRES1=111
      IRES2=999
      IDECF=0
      MRED=20
      MOS3=2
      ETA0=1.0D-15
      ETA2=1.0D-12
      ETA9=1.0D120
      EPS0=1.0D-6
      EPS1=1.0D-4
      EPS2=2.5D-1
      ALF1=1.0D-10
      ALF2=1.0D10
      RMAX=ETA9
      DMAX=ETA9
      FMAX=1.0D20
      IF (IEST.LE.0) FMIN=-1.0D60
      IF (IEST.GT.0) IEST=1
      IF (XMAX.LE.0.0D0) XMAX=1.0D16
      IF (TOLX.LE.0.0D0) TOLX=1.0D-16
      IF (TOLF.LE.0.0D0) TOLF=1.0D-12
      IF (TOLB.LE.0.0D0) TOLB=FMIN+1.0D-12
      IF (TOLG.LE.0.0D0) TOLG=1.0D-8
      IF (ETA3.LE.0.0D0) ETA3=1.0D-12
      IF (ETA5.LE.0.0D0) ETA5=2.5D-1
      IF (MIT.LE.0) MIT=20000
      IF (MFV.LE.0) MFV=20000
      MFG=MFV
      MA=IAG(NA+1)-1
      CALL MXVINP (NF+1, IH)
      CALL MXVINP (NF, JH)
      CALL PFSET3 (NF, NA, M, MMAX, IH, JH, IAG, JAG, ITERM)
      IF (ITERM.LT.0) STOP
      IF (ITERM.NE.0) GO TO 80
      MH=0
      CALL MXVINE (IH(NF+1)-1, JH)
      CALL MXSPCC (NF, M, MH, MMAX, H, IH, JH, PSL, PERM, INVP, WN11,
     &WN12, WN13, WN14, IER)
      IF (IER.NE.0) THEN
        ITERM=IER
      END IF
      IF (ITERM.NE.0) GO TO 70
!
!     BUNDLE VARIABLE METRIC METHOD
!
      LD=-1
      KD=1
      ISNA=2
      KIT=-(IRES1*NF+IRES2)
      FO=FMIN
      NB=0
      JR=0
      JE=0
      NNK=0
      NNV=0
!
!     MODEL DESCRIPTION
!
      CALL PA1SF3 (NF, NA, X, GA, G, AG, IAG, JAG, F, AF, KD, LD, ISNA,
     &NFV, NFG)
!
!     END OF MODEL DESCRIPTION
!
      DF=ABS(F)+1.0D0
      CALL PYBUN1 (NF, MB, NB, X, G, F, AX, AY, AZ, ITERS)
   10 IF (ITERS.GT.0) THEN
        JC=0
        ALFN=0.0D0
        ALFV=0.0D0
        CALL MXVCOP (NF, G, GP)
      END IF
      CALL PYTRCG (NF, N, IX, GP, UMAX, GMAX, KBF, IOLD)
      IF (ITERM.LT.0) GO TO 70
      IF (ABS(IPRNT).GT.1) WRITE (6,'(1X,''NIT='',I5,2X,''NFV='',I5,2X,'
     &'NFG='',I5,2X,       ''F='', G16.9,2X,''G='',E10.3)') NIT,NFV,NFG,
     &F,GMAX
      IF (F.LE.TOLB) THEN
        ITERM=3
        GO TO 70
      END IF
      IF (NIT.GE.MIT) THEN
        ITERM=11
        GO TO 70
      END IF
      IF (NFV.GE.MFV) THEN
        ITERM=12
        GO TO 70
      END IF
      IF (NFG.GE.MFG) THEN
        ITERM=13
        GO TO 70
      END IF
      ITERM=0
      IF (NF.NE.0.AND.NIT-KIT.GE.IRES1*NF+IRES2) THEN
        IREST=MAX(IREST,1)
      END IF
      NIT=NIT+1
   20 IF (IREST.GT.0) THEN
        CALL MXBSMI (NA, AH, IAG)
!      LD=MIN(LD,1)
        IDECF=0
        IF (KIT.LT.NIT) THEN
          NRES=NRES+1
          KIT=NIT
        ELSE
          ITERM=-10
          IF (ITERS.LT.0) ITERM=ITERS-5
        END IF
      END IF
      IF (ITERM.NE.0) GO TO 70
      CALL MXVSET (IH(NF+1)-1, 0.0D0, H)
      CALL PFSET4 (NA, H, IH, JH, AH, IAG, JAG)
      CALL PYTSCH (NF, IX, H, IH, JH, KBF)
      IF (JE.GT.0) GO TO 30
!
!     DIRECTION DETERMINATION
!
      MM=IH(NF+1)-1
      CALL PDSLM1 (NF, MMAX, MH, IX, GP, H, IH, JH, S, XO, PSL, PERM,
     &WN11, WN12, GNORM, SNORM, ETA2, KBF, IDECF, NDEC, ITERD, ITERM)
      JC=1
      IF (JC.EQ.1) THEN
        CALL MXVDIR (NF, -ETA3, GP, S, S)
        SNORM=SQRT(MXVDOT(NF,S,S))
      END IF
!
!     END OF DIRECTION DETERMINATION
!
      IF (KD.GT.0) P=MXVDOT(NF,GP,S)
!
!     TEST ON DESCENT DIRECTION AND PREPARATION OF LINE SEARCH
!
      IF (ITERD.LT.0) THEN
        ITERM=ITERD
      ELSE
!
!     TEST ON DESCENT DIRECTION
!
        IF (SNORM.LE.0.0D0) THEN
          IREST=MAX(IREST,1)
        ELSE IF (P+EPS0*GNORM*SNORM.LE.0.0D0) THEN
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
      IF (ITERS.EQ.0) IREST=0
      IF (SNORM.LE.0.0D0) IREST=0
      IF (ITERM.NE.0) GO TO 70
      IF (IREST.NE.0) GO TO 20
      XNORM=-P+2.0D0*ALFV
   30 CONTINUE
      IF (XNORM.LE.TOLG) THEN
        IF (SNORM.LE.0.0D0) ITERM=4
        NTESX=NTESX+1
        IF (ITERS.GT.0.AND.DF.LT.1.0D2*TOLF*MAX(ABS(F),1.0D0)) ITERM=4
        IF (NTESX.GE.2.AND.NNK.GT.1) ITERM=4
      ELSE
        NTESX=0
      END IF
      IF (ITERM.NE.0) GO TO 70
      IF (SNORM.GT.0.0D0) RMAX=XMAX/SNORM
      RMIN=MIN(1.0D-10,RMAX/1.0D1)
      CALL PYTRCS (NF, X, IX, XO, X, X, G, GO, S, RO, FP, FO, F, PO, P,
     &RMAX, ETA9, KBF)
      CALL MXVCOP (MA, AG, AGO)
      IF (RMAX.LE.RMIN) THEN
        R=0.0D0
        GO TO 60
      END IF
   40 CALL PS1L18 (NF, MB, NB, X, G, S, XO, AZ, AY, AX, R, RP, FO, F,
     &PO, P, RMIN, RMAX, SNORM, XNORM, EPS1, EPS2, ETA5, ETA9, KD, LD,
     &JE, MOS3, ITERS, ISYS)
      IF (ISYS.EQ.0) GO TO 50
      CALL MXVDIR (NF, R, S, XO, X)
      CALL PA1SF3 (NF, NA, X, GA, G, AG, IAG, JAG, F, AF, KD, LD, ISNA,
     &NFV, NFG)
      P=MXVDOT(NF,G,S)
      GO TO 40
   50 CONTINUE
      NNV=NNV+1
      POM=DF
      IF (ABS(FO-F).GE.DF*1.0D-5) POM=ABS(FO-F)
      IF (ITERS.GT.0) DF=POM
      IF (POM.LE.TOLF*MAX(ABS(F),1.0D0).OR.FO.EQ.F) THEN
        NTESF=NTESF+1
        IF (NTESF.GE.MTESF) THEN
          F=FO
          ITERM=2
          GO TO 70
        END IF
      ELSE
        NTESF=0
      END IF
      CALL PYBUN1 (NF, MB, NB, X, G, F, AX, AY, AZ, ITERS)
      IF (ITERS.EQ.0) THEN
        NNK=NNK+1
        ALFN=MAX(ABS(FO-F+P*R),ETA5*(SNORM*R)**MOS3)
        MAM=MA
        IF (NNK.EQ.1) THEN
          CALL PYABU2 (NF, H(MM+1), JH(MM+1), PSL, PERM, G, GP, S, GS,
     &     ALFN, ALFV, ETA3, JC)
        ELSE
          CALL PYABU1 (NF, H(MM+1), JH(MM+1), PSL, PERM, G, GO, GP, S,
     &     GS, XS, ALFN, ALFV, ETA3, JC)
        END IF
        F=FO
      ELSE
        NNK=0
      END IF
      POM=P
      CALL PYTRCD (NF, X, IX, XO, G, GO, R, F, FO, P, PO, DMAX, KBF, KD,
     & LD, ITERS)
      P=POM
      POM=MXVSAB(NF,GO)
      IF (POM.EQ.0.0D0.AND.ITERS.GT.0) THEN
        ITERM=6
        GO TO 70
      ELSE
        JE=0
      END IF
      POM=MXVDOT(NF,XO,GO)
      IF (POM.GT.R*0.0D0.AND.ABS(POM).GT.1.0D-6*MXVNOR(NF,XO)*MXVNOR(NF,
     &GO)) THEN
        IDECF=0
        CALL PUBVI2 (NA, AH, IAG, JAG, AG, AGO, XO, S, GS, ETA9, NNK,
     &   NIT, ITERH)
      END IF
   60 CONTINUE
      GO TO 10
   70 CONTINUE
      GMAX=XNORM
   80 CONTINUE
      IF (IPRNT.GT.1.OR.IPRNT.LT.0) WRITE (6,'(1X,''EXIT FROM PSUM :'')'
     &)
      IF (IPRNT.NE.0) WRITE (6,'(1X,''NIT='',I5,2X,''NFV='',I5,2X,''NFG=
     &'',I5,2X,       ''F='', G16.9,2X,''G='',E10.3,2X,''ITERM='',I3)')
     &NIT,NFV,NFG,F,GMAX,ITERM
      IF (IPRNT.LT.0) WRITE (6,'(1X,''X='',5(G14.7,1X):/(3X,5(G14.7,1X))
     &)') (X(I),I=1,NF)
      END