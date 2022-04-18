      PROGRAM EXMPL
C
C
Ccc   INCLUDE 'PARNEWTON'
C
C PARNEWTON
C
C Parameters for Newton process
C MAXNIT : Max. number of Newton iterations
C MAXJAC : Max. number of Jacobian / preconditioner evaluations during
C          a Newton process
C TOLNEW : Tolerance for Newton process:
C          rho/(1-rho)*|| corr.||_w < TOLNEW
      INTEGER MAXNIT, MAXJAC
      REAL TOLNEW
      PARAMETER (MAXNIT = 10, MAXJAC = 2, TOLNEW = 1.0)
C
C end INCLUDE 'PARNEWTON'
C
C
Ccc   INCLUDE 'PARGCRO'
C
C PARGCRO
C
C Parameters for linear system solver GCRO + (block-)diagonal
C    preconditioner
C IDIAGP : 0: block-diagonal + first order derivatives
C          1: block-diagonal neglecting first order derivatives
C          2: diagonal + first order derivatives
C          3: diagonal neglecting first order derivatives
C NRRMAX : Max. number of restarts of outer loop
C MAXLR  : Max. number of iterations in outer loop
C MAXL   : Max. number of iterations in GMRES inner loop
C TOLLSC : Tolerance for linear system solver
      INTEGER IDIAGP, NRRMAX, MAXLR, MAXL
      REAL TOLLSC
      PARAMETER (NRRMAX = 1, MAXLR = 5, MAXL = 20)
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 10)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
C
C
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK, LENLWK
      PARAMETER (MXLEV=2, NPD=3, NPTS=5000)
      PARAMETER (LENIWK=NPTS*(5*MXLEV+14),
     +           LENRWK=NPTS*NPD*(5*MXLEV+9 +
     +                  9*NPD+(2*MAXLR+MAXL+6+NPD)),
     +           LENLWK=2*NPTS)
C
C-----------------------------------------------------------------------
C
      INTEGER LUNDMP
      PARAMETER (LUNDMP = 89)
C
      CHARACTER FILE*7
      INTEGER NPDE, INFO(7), IWK(LENIWK), MNTR, I
      LOGICAL LWK(LENLWK)
      REAL T, TOUT(4), DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, RINFO(2+3*NPD), RWK(LENRWK)
C
C First call of VLUGR2
      MNTR = 0
      NPDE = 3
      T    = 0.0
      TOUT(1) = 500.0
      TOUT(2) = 5000.0
      TOUT(3) = 10000.0
      TOUT(4) = 20000.0
      DT   = 0.1
      XL = 0.0
      XR = 1.0
      YL = 0.0
      YU = 1.0
      DX = 0.05
      DY = 0.05
      TOLS = 0.1
      TOLT = 0.1
      INFO(1) = 1
C MAXLEV
      INFO(2) = 3
C Domain is a rectangle
      INFO(3) = 0
C Linear system solver
      PRINT *, 'Lin.sys.solver; BiCGStab, GCRO or matrix-free GCRO ?'
      PRINT *, '   (0 / 10,11,12,13 / 20,21,22,23 ) ?'
      READ *, INFO(4)
      OPEN (UNIT=61,FILE='RunInfo')
C Write integration history to unit # 61
      INFO(5) = 61
C Write Newton info to unit # 61
      INFO(6) = 61
C Write Linear system solver info to unit # 61
      INFO(7) = 61
C DTMIN = 1E-3
      RINFO(1) = 1.0E-3
C DTMAX = 1.0
      RINFO(2) = 10000.0
C UMAX
      RINFO(3) = 1.1E+5
      RINFO(4) = 0.25
      RINFO(5) = 292.0
C SPCWGT = 1.0
      RINFO(6) = 1.0
      RINFO(7) = 1.0
      RINFO(8) = 1.0
C TIMWGT = 1.0
      RINFO( 9) = 1.0
      RINFO(10) = 1.0
      RINFO(11) = 1.0
C
C Call main routine
      FILE='DUMP'
      DO 10 I = 1, 4
         CALL  VLUGR2 (NPDE, T, TOUT(I), DT, XL, YL, XR, YU, DX, DY,
     +      TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK,
     +      LWK, LENLWK, MNTR)
C
C Save info on file
         WRITE(FILE(5:7),'(I3.3)') I
         OPEN(UNIT=LUNDMP,FILE=FILE,FORM='UNFORMATTED')
         CALL DUMP (LUNDMP, RWK, IWK)
         CLOSE(LUNDMP)
C
Check MNTR value
         IF (MNTR .NE. 1) THEN
            PRINT *, 'VLUGR2 returned with MNTR=', MNTR
            STOP
         ENDIF
   10 CONTINUE
      END
      SUBROUTINE PDEIV (T, X, Y, U, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), U(NPTS,NPDE)
C
Ccc PURPOSE:
C Define (initial) solution of PDE.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which (initial) solution should be given
C X      : IN. Array of X-coordinates for the gridpoints
C Y      : IN. Array of Y-coordinates for the gridpoints
C U      : OUT. Array of PDE component values for the gridpoints.
C NPTS   : IN. Number of gridpoints
C NPDE   : IN. # PDE components
C
C-----------------------------------------------------------------------
C
      INTEGER I
C
      REAL N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
Ccc Problem parameters
      N      = 0.4
      KAPPA  = 1.0E-10
      G      = 9.81
      DM     = 0.0
      AT     = 0.002
      AL     = 0.01
      CF     = 4182.0
      TKAPPA = 4.0
      LT     = 0.001
      LL     = 0.01
      CS     = 840.0
      RHOS   = 2500.0
      RHO0   = 1.0E+3
      T0     = 290.0
      P0     = 1.0E+5
      ALPHA  = -3.0E-4
      BETA   = 4.45E-10
      GAMMA  = LOG(1.2)
      MU0    = 1.0E-3
      W0     = 0.25
      QC     = 1.0E-4
      TC     = 292.0
C
Ccc Initial solution
      DO 10 I = 1, NPTS
         U(I,1) = P0 + (1.0 - Y(I))*RHO0*G
         U(I,2) = 0.0
         U(I,3) = T0
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PDEF (T, X, Y, U, UT, UX, UY, UXX, UXY, UYY, RES,
     +   NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), U(NPTS,NPDE),
     +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE),
     +     UXX(NPTS,NPDE), UXY(NPTS,NPDE), UYY(NPTS,NPDE),
     +     RES(NPTS,NPDE)
C
Ccc PURPOSE:
C Define residual of PDE on interior of domain. Boundary values will be
C overwritten later on.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which residual should be evaluated
C X      : IN. Array of X-coordinates for the gridpoints
C Y      : IN. Array of Y-coordinates for the gridpoints
C U      : IN. Array of PDE components for the gridpoints.
C UT     : IN. Array of time derivative of PDE components
C UX     : IN. -I
C UY     : IN.  I
C UXX    : IN.  I Arrays containing space derivatives of PDE components
C UXY    : IN.  I
C UYY    : IN. -I
C RES    : OUT. Array containg PDE residual at gridpoints in interior of
C              domain. The residual values at boundary points will be
C              overwritten by a call to PDEBC.
C NPTS   : IN. Number of gridpoints
C NPDE   : IN. Number of PDE components
C
C
C
Ccc   INCLUDE 'CMNCMMACH'
C
C CMNCMMACH
C
C COMMON with `machine numbers'
C LUNOUT : Logical unit # standard output         -I
C LUNERR : Logical unit # standard error           I  Set in the routine
C UROUND : Smallest machine number such that       I      MACNUM
C          1.0+UROUND > 1.0 and 1.0-UROUND < 1.0   I
C XMIN   : Smallest floating-point number         -I
      INTEGER LUNOUT, LUNERR
      REAL UROUND, XMIN
      COMMON /IMACH/ LUNOUT, LUNERR
      COMMON /RMACH/ UROUND, XMIN
      SAVE /IMACH/, /RMACH/
C
C end INCLUDE 'CMNCMMACH'
C
C-----------------------------------------------------------------------
C
      REAL N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
      INTEGER I
      REAL P, PT, PX, PY, W, WT, WX, WY, T1, TT, TX, TY,
     +   RHO, RHOX, RHOY,
     +   MU, MUX, MUY, KAPMU, KAPMU2, KAPMUX, KAPMUY, Q1, Q2, QL,
     +   ND11, ND12, ND22, H11, H12, H22,
     +   PXX, PXY, PYY, WXX, WXY, WYY, TXX, TXY, TYY,
     +   ND11Q1, ND11Q2, ND12Q1, ND12Q2, ND22Q1, ND22Q2,
     +   H11Q1, H11Q2, H12Q1, H12Q2, H22Q1, H22Q2, Q1X, Q1Y, Q2X, Q2Y,
     +   ND11X, ND12X, ND12Y, ND22Y, JW1, JW2, JW1X, JW2Y,
     +   H11X, H12X, H12Y, H22Y, JT1X, JT2Y
C
      DO 10 I = 1, NPTS
         P      = U(I,1)
         PT     = UT(I,1)
         PX     = UX(I,1)
         PY     = UY(I,1)
         W      = U(I,2)
         WT     = UT(I,2)
         WX     = UX(I,2)
         WY     = UY(I,2)
         T1     = U(I,3)
         TT     = UT(I,3)
         TX     = UX(I,3)
         TY     = UY(I,3)
         RHO    = RHO0*EXP(ALPHA*(T1-T0)+BETA*(P-P0)+GAMMA*W)
         RHOX   = RHO*(ALPHA*TX+BETA*PX+GAMMA*WX)
         RHOY   = RHO*(ALPHA*TY+BETA*PY+GAMMA*WY)
         MU     = MU0*(1+1.85*W-4.0*W*W)
         MUX    = MU0*(1.85*WX-8.0*W*WX)
         MUY    = MU0*(1.85*WY-8.0*W*WY)
         KAPMU  = KAPPA/MU
         KAPMU2 = -KAPMU/MU
         KAPMUX = KAPMU2*MUX
         KAPMUY = KAPMU2*MUY
         Q1     = -KAPMU*PX
         Q2     = -KAPMU*(PY+RHO*G)
         QL     = MAX(SQRT(Q1*Q1+Q2*Q2),UROUND)
         ND11   = N*DM + AT*QL + (AL-AT)*Q1*Q1/QL
         ND12   = (AL-AT)*Q1*Q2/QL
         ND22   = N*DM + AT*QL + (AL-AT)*Q2*Q2/QL
         H11    = TKAPPA + LT*QL + (LL-LT)*Q1*Q1/QL
         H12    = (LL-LT)*Q1*Q2/QL
         H22    = TKAPPA + LT*QL + (LL-LT)*Q2*Q2/QL
         PXX    = UXX(I,1)
         PXY    = UXY(I,1)
         PYY    = UYY(I,1)
         WXX    = UXX(I,2)
         WXY    = UXY(I,2)
         WYY    = UYY(I,2)
         TXX    = UXX(I,3)
         TXY    = UXY(I,3)
         TYY    = UYY(I,3)
         ND11Q1 = (AT + (AL-AT)*(2-(Q1/QL)**2))*Q1/QL
         ND11Q2 = (AT - (AL-AT)*((Q1/QL)**2))*Q2/QL
         ND12Q1 = (AL-AT)*(Q2/QL)**3
         ND12Q2 = (AL-AT)*(Q1/QL)**3
         ND22Q1 = (AT - (AL-AT)*((Q2/QL)**2))*Q1/QL
         ND22Q2 = (AT + (AL-AT)*(2-(Q2/QL)**2))*Q2/QL
         H11Q1  = (LT + (LL-LT)*(2-(Q1/QL)**2))*Q1/QL
         H11Q2  = (LT - (LL-LT)*((Q1/QL)**2))*Q2/QL
         H12Q1  = (LL-LT)*(Q2/QL)**3
         H12Q2  = (LL-LT)*(Q1/QL)**3
         H22Q1  = (LT - (LL-LT)*((Q2/QL)**2))*Q1/QL
         H22Q2  = (LT + (LL-LT)*(2-(Q2/QL)**2))*Q2/QL
         Q1X    = -(KAPMUX*PX+KAPMU*PXX)
         Q1Y    = -(KAPMUY*PX+KAPMU*PXY)
         Q2X    = -(KAPMUX*(PY+RHO*G)+KAPMU*(PXY+RHOX*G))
         Q2Y    = -(KAPMUY*(PY+RHO*G)+KAPMU*(PYY+RHOY*G))
         ND11X  = ND11Q1*Q1X + ND11Q2*Q2X
         ND12X  = ND12Q1*Q1X + ND12Q2*Q2X
         ND12Y  = ND12Q1*Q1Y + ND12Q2*Q2Y
         ND22Y  = ND22Q1*Q1Y + ND22Q2*Q2Y
         JW1    = -(ND11*WX + ND12*WY)
         JW2    = -(ND12*WX + ND22*WY)
         JW1X   = -(ND11X*WX+ND11*WXX + ND12X*WY+ND12*WXY)
         JW2Y   = -(ND12Y*WX+ND12*WXY + ND22Y*WY+ND22*WYY)
         H11X   = H11Q1*Q1X + H11Q2*Q2X
         H12X   = H12Q1*Q1X + H12Q2*Q2X
         H12Y   = H12Q1*Q1Y + H12Q2*Q2Y
         H22Y   = H22Q1*Q1Y + H22Q2*Q2Y
         JT1X   = -(H11X*TX+H11*TXX + H12X*TY+H12*TXY)
         JT2Y   = -(H12Y*TX+H12*TXY + H22Y*TY+H22*TYY)
C
         RES(I,1) = N*RHO*(BETA*PT+GAMMA*WT+ALPHA*TT) +
     +      RHOX*Q1+RHO*Q1X + RHOY*Q2+RHO*Q2Y
         RES(I,2) = N*RHO*WT +
     +      RHO*Q1*WX + RHO*Q2*WY +
     +      RHOX*JW1+RHO*JW1X + RHOY*JW2+RHO*JW2Y
         RES(I,3) = (N*CF*RHO+(1-N)*RHOS*CS)*TT +
     +      RHO*CF*Q1*TX + RHO*CF*Q2*TY +
     +      JT1X + JT2Y
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PDEBC (T, X, Y, U, UT, UX, UY, RES,
     +   NPTS, NPDE, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      REAL T, X(NPTS), Y(NPTS), U(NPTS,NPDE),
     +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE),
     +     RES(NPTS,NPDE)
C
Ccc PURPOSE:
C Define residual of boundary equations of PDE. The residual on interior
C points has already been stored in RES.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which BC's should be evaluated
C X      : IN. Array of X-coordinates for the gridpoints
C Y      : IN. Array of Y-coordinates for the gridpoints
C U      : IN. Array of PDE components for the gridpoints.
C UT     : IN. Array of time derivative of PDE components
C UX     : IN. -I
C UY     : IN. -I Arrays containing space derivatives of PDE components
C RES    : INOUT.
C          IN: PDE residual for interior points (set by PDEF)
C          OUT: Array with PDE residual at physical boundary points
C               inserted
C NPTS   : IN. Number of grid components
C NPDE   : IN. Number of PDE components
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical boundaries and corners in
C                     actual domain.
C             NB. corners should be stored as an independent boundary
C             (cf. ILBND). The order in LLBND should be first the
C             boundaries and then the corners.
C          LLBND(1:NBNDS): pointers to a specific boundary or corner in
C                          LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Lower boundary -I
C                     2: Left  boundary  I
C                     3: Upper boundary  I max. first order derivative
C                     4: Right boundary -I
C                    12: Lowerleft corner  -I
C                    23: Leftupper corner   I corners of 90 degrees
C                    34: Upperright corner  I (external corners)
C                    41: Rightlower corner -I max. first order deriv.
C                    21: Leftlower corner  -I
C                    32: Upperleft corner   I corners of 270 degrees
C                    43: Rightupper corner  I (internal corners)
C                    14: Lowerright corner -I max. first order deriv.
C LBND   : IN. (NBDPTS)
C          LBND(LB): pointer to boundary point in actual grid
C                    structure (as in X, Y, and U)
C
C
Ccc   INCLUDE 'CMNCMMACH'
C
C CMNCMMACH
C
C COMMON with `machine numbers'
C LUNOUT : Logical unit # standard output         -I
C LUNERR : Logical unit # standard error           I  Set in the routine
C UROUND : Smallest machine number such that       I      MACNUM
C          1.0+UROUND > 1.0 and 1.0-UROUND < 1.0   I
C XMIN   : Smallest floating-point number         -I
      INTEGER LUNOUT, LUNERR
      REAL UROUND, XMIN
      COMMON /IMACH/ LUNOUT, LUNERR
      COMMON /RMACH/ UROUND, XMIN
      SAVE /IMACH/, /RMACH/
C
C end INCLUDE 'CMNCMMACH'
C
C
C-----------------------------------------------------------------------
C
C
      REAL N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
      INTEGER I, J, K, NBNDS
      REAL P, PY, W, T1, RHO, MU, KAPMU, Q2
C
      NBNDS = LLBND(0)
      DO 10 J = 1, NBNDS
         IF (ILBND(J) .EQ. 1) THEN
C
C yL boundary: dp/dy - rho.g2 = 0
C              dw/dy = 0                0<x<1/11, 2/11<x<9/11, 10/11<x<1
C              dT/dy = 0
C              q2 = qc
C              dw/dy = 0                1/11<=x<=2/11, 9/11<=x<=10/11
C              dT/dy = 0
C
            DO 20 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF ((     0 .LT. X(I) .AND. X(I) .LT. 1.0/11) .OR.
     +             (2.0/11 .LT. X(I) .AND. X(I) .LT. 9.0/11) .OR.
     +            (10.0/11 .LT. X(I) .AND. X(I) .LT. 1)) THEN
                  P      = U(I,1)
                  PY     = UY(I,1)
                  W      = U(I,2)
                  T1     = U(I,3)
                  RHO    = RHO0*EXP(ALPHA*(T1-T0)+BETA*(P-P0)+GAMMA*W)
                  MU     = MU0*(1+1.85*W-4.0*W*W)
                  KAPMU  = KAPPA/MU
                  Q2     = -KAPMU*(PY+RHO*G)
                  RES(I,1) = Q2
                  RES(I,2) = UY(I,2)
                  RES(I,3) = UY(I,3)
               ELSE
                  P      = U(I,1)
                  PY     = UY(I,1)
                  W      = U(I,2)
                  T1     = U(I,3)
                  RHO    = RHO0*EXP(ALPHA*(T1-T0)+BETA*(P-P0)+GAMMA*W)
                  MU     = MU0*(1+1.85*W-4.0*W*W)
                  KAPMU  = KAPPA/MU
                  Q2     = -KAPMU*(PY+RHO*G)

                  RES(I,1) = Q2 - QC
                  RES(I,2) = W - W0
                  RES(I,3) = T1 - TC
               ENDIF
   20       CONTINUE
         ELSE IF (ILBND(J) .EQ. 2 .OR. ILBND(J) .EQ. 4) THEN
C
C xL boundary. dp/dx = 0 dT/dx = 0
C xR boundary` dw/dx = 0
C
            DO 30 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = UX(I,1)
               RES(I,2) = UX(I,2)
               RES(I,3) = UX(I,3)
   30       CONTINUE
         ELSE IF (ILBND(J) .EQ. 3) THEN
C
C yU boundary: p(x,0) = p0
C              dw/dy = 0                xL < x < xR
C              dT/dy = 0
C
            DO 40 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = U(I,1) - P0
               RES(I,2) = UY(I,2)
               RES(I,3) = UY(I,3)
   40       CONTINUE
         ELSE
CDIR$ NOVECTOR
C
C Corners: dp/dx = 0 dT/dx = 0
C          dw/dx = 0
C
            DO 50 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = UX(I,1)
               RES(I,2) = UX(I,2)
               RES(I,3) = UX(I,3)
   50       CONTINUE
         ENDIF
   10 CONTINUE

      RETURN
      END
