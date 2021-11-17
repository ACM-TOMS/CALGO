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
      DOUBLE PRECISION TOLNEW
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
      DOUBLE PRECISION TOLLSC
      PARAMETER (NRRMAX = 1, MAXLR = 5, MAXL = 20)
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 10)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
C
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK, LENLWK
      PARAMETER (MXLEV=2, NPD=2, NPTS=10000)
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
      DOUBLE PRECISION T, TOUT(4), DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, RINFO(2+3*NPD), RWK(LENRWK)
C
C First call of VLUGR2
      MNTR = 0
      NPDE = 2
      T    = 0.0
      TOUT(1) = 10000.0
      TOUT(2) = 20000.0
      TOUT(3) = 30000.0
      TOUT(4) = 100000.0
      DT   = 0.1
C Since domain is not a rectangle the domain parameters have not to be
C specified here
      TOLS = 0.1
      TOLT = 0.1
      INFO(1) = 1
C MAXLEV
      INFO(2) = 3
C Domain is not a rectangle
      INFO(3) = 1
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
C DTMIN
      RINFO(1) = 1.0D-3
C DTMAX
      RINFO(2) = 50000.0
C UMAX
      RINFO(3) = 1.1D+5
      RINFO(4) = 0.25
C SPCWGT = 1.0
      RINFO(5) = 1.0
      RINFO(6) = 1.0
C TIMWGT = 1.0
      RINFO(7) = 1.0
      RINFO(8) = 1.0
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
      LOGICAL FUNCTION INIDOM (MAXPTS, XL, YL, XR, YU, DX, DY,
     +   LROW, IROW, ICOL, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER MAXPTS, LROW(0:*), IROW(*), ICOL(*),
     +   LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION XL, YL, XR, YU, DX, DY
C
Ccc PURPOSE:
C Define initial domain. NB. Boundaries should consist of as many points
C as are necessary to employ second order space discretization, i.e.,
C a boundary enclosing the internal part of the domain should not
C include less than 3 grid points including the corners. If Neumann
C boundaries are used the minimum is 4 since otherwise the Jacobian
C matrix will be singular.
C
C A (virtual) rectangle is placed upon the (irregular) domain. The
C lowerleft point of this rectangle is (XL,YL) in physical coordinates
C and (0,0) in column, resp. row coordinates. The upperright point is
C (XR,YU) resp. (Nx, Ny), where Nx = (XR-XL)/DX and Ny = (YU-YL)/DY.
C Only real grid points are stored.
C The coordinate values of the initial grid should be stored rowwise
C in LROW, IROW, ICOL.
C Pointers to the boundary points should be stored in a list together
C with the type of the boundary. (LLBND, ILBND, LBND)
C
C On exit INIDOM = .FALSE. if the # grid points required is larger
C than MAXPTS and MAXPTS is set to the required # points.
C
Ccc PARAMETER DESCRIPTION:
C MAXPTS : INOUT.
C          IN:  Max. # grid points allowed by the available workspace
C          OUT: # grid points required, if larger than # points allowed
C XL     : OUT. X-coordinate of lower-left point of virtual rectangle
C YL     : OUT. Y-coordinate of lower-left point of virtual rectangle
C XR     : OUT. X-coordinate of upper-right point of virtual rectangle
C YU     : OUT. Y-coordinate of upper-right point of virtual rectangle
C DX     : OUT. Grid width in X-direction
C DY     : OUT. Grid width in Y-direction
C LROW   : OUT. INTEGER array of dimension (0:LROW(0)+1)
C          LROW(0) = NROWS: Actual # rows in grid
C          LROW(1:NROWS): pointers to the start of a row in the grid
C                         structure
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : OUT. INTEGER array of dimension (NROWS)
C          IROW(IR): row number of row IR in virtual rectangle
C ICOL   : OUT. INTEGER array of dimension (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual
C                     rectangle
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
C LBND   : OUT. INTEGER array of dimension (NBDPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C                      structure
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
C NX should be even and NY a quintuple
      INTEGER NX, NY
      PARAMETER (NX = 20, NY = 20)
C
      INTEGER I, IPT, J, NROWS, NPTS, NBNDS, NX1, NY1, NY2
C
Ccc Make initial grid, check MAXPTS against rough estimate of NPTS
      IF (MAXPTS .LT. (NX+1)*(NY+1)) THEN
         INIDOM = .FALSE.
         MAXPTS = (NX+1)*(NY+1)
         RETURN
      ELSE
         INIDOM = .TRUE.
      ENDIF
      NROWS = NY+1

      XL = 0.0
      YL = 0.0
      XR = 1.0
      YU = 1.0
      DX = (XR-XL)/NX
      DY = (YU-YL)/NY

      NX1 = NX/2
      NY1 = NINT(NY*0.4)
      NY2 = NINT(NY*0.6)
C
C Make grid structure
      LROW(0) = NROWS
      IPT = 1
      DO 10 I = 0, NY1
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 20 J = 0, NX
            ICOL(IPT) = J
            IPT = IPT + 1
   20    CONTINUE
   10 CONTINUE
      DO 30 I = NY1+1, NY2-1
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 40 J = NX/2, NX
            ICOL(IPT) = J
            IPT = IPT + 1
   40    CONTINUE
   30 CONTINUE
      DO 50 I = NY2, NY
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 60 J = 0, NX
            ICOL(IPT) = J
            IPT = IPT + 1
   60    CONTINUE
   50 CONTINUE
      LROW(NROWS+1) = IPT
      NPTS = IPT-1
C
C Boundaries
      NBNDS = 16
      ILBND(1) = 1
      ILBND(2) = 2
      ILBND(3) = 3
      ILBND(4) = 2
      ILBND(5) = 1
      ILBND(6) = 2
      ILBND(7) = 3
      ILBND(8) = 4
      ILBND( 9) = 12
      ILBND(10) = 23
      ILBND(11) = 32
      ILBND(12) = 21
      ILBND(13) = 12
      ILBND(14) = 23
      ILBND(15) = 34
      ILBND(16) = 41
      LLBND(0) = NBNDS
      LLBND(1) = 1
      LLBND(2) = LLBND(1) + (NX-1)
      LLBND(3) = LLBND(2) + (NY1-1)
      LLBND(4) = LLBND(3) + (NX1-1)
      LLBND(5) = LLBND(4) + (NY2-NY1-1)
      LLBND(6) = LLBND(5) + (NX1-1)
      LLBND(7) = LLBND(6) + (NY1-1)
      LLBND(8) = LLBND(7) + (NX-1)
      LLBND( 9) = LLBND( 8) + (NY-1)
      LLBND(10) = LLBND( 9) + 1
      LLBND(11) = LLBND(10) + 1
      LLBND(12) = LLBND(11) + 1
      LLBND(13) = LLBND(12) + 1
      LLBND(14) = LLBND(13) + 1
      LLBND(15) = LLBND(14) + 1
      LLBND(16) = LLBND(15) + 1
      LLBND(17) = LLBND(16) + 1
C Lower and upper boundary pointers
      DO 100 J = 1, NX-1
         LBND(LLBND(1)+J-1) = J + 1
         LBND(LLBND(7)+J-1) = NPTS - J
  100 CONTINUE
C Left boundary pointers
      DO 120 I = 1, NY1-1
         LBND(LLBND(2)+I-1) = I*(NX+1) + 1
         LBND(LLBND(6)+I-1) = NPTS - (I+1)*(NX+1) + 1
  120 CONTINUE
      DO 130 I = 1, NY2-NY1-1
         LBND(LLBND(4)+I-1) = NY1*(NX+1) + (I+1)*(NX1+1)
  130 CONTINUE
      DO 140 I = 1, NX1-1
         LBND(LLBND(3)+I-1) = NY1*(NX+1) + 1 + I
         LBND(LLBND(5)+I-1) = NPTS - (NY1+1)*(NX+1)  + 1 + I
  140 CONTINUE
C Right boundary pointers
      DO 110 I = 1, NY1
         LBND(LLBND(8)+I-1) = (I+1)*(NX+1)
  110 CONTINUE
      J = LLBND(8)+NY1-1
      DO 113 I = 1, NY2-NY1-1
         LBND(J+I) = LBND(J) + I*(NX1+1)
  113 CONTINUE
      J = LLBND(8)+NY2-1
      DO 116 I = 0, NY1-1
         LBND(J+I) = LBND(J-1) + (I+1)*(NX+1)
  116 CONTINUE
C Corners
      LBND(LLBND( 9)) = 1
      LBND(LLBND(16)) = NX+1
      LBND(LLBND(10)) = NY1*(NX+1) + 1
      LBND(LLBND(11)) = LBND(LLBND(10)) + NX1
      LBND(LLBND(13)) = (NY1+1)*(NX+1) + (NY2-NY1-1)*(NX1+1) + 1
      LBND(LLBND(12)) = LBND(LLBND(13)) + NX1
      LBND(LLBND(14)) = NPTS - NX
      LBND(LLBND(15)) = NPTS

      RETURN
      END
      SUBROUTINE PDEIV (T, X, Y, U, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), U(NPTS,NPDE)
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
      DOUBLE PRECISION N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL,
     +   AT, CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
Ccc Problem parameters
      N      = 0.4
      KAPPA  = 1.0D-10
      G      = 9.81
      DM     = 0.0
      AT     = 0.002
      AL     = 0.01
      RHO0   = 1.0D+3
      P0     = 1.0D+5
      BETA   = 0.0
      GAMMA  = LOG(1.2)
      MU0    = 1.0D-3
      W0     = 0.25
      QC     = 1.0D-4
C
Ccc Initial solution
      DO 10 I = 1, NPTS
         U(I,1) = P0 + (1.0 - Y(I))*RHO0*G
         U(I,2) = 0.0
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
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), U(NPTS,NPDE),
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
      DOUBLE PRECISION UROUND, XMIN
      COMMON /IMACH/ LUNOUT, LUNERR
      COMMON /RMACH/ UROUND, XMIN
      SAVE /IMACH/, /RMACH/
C
C end INCLUDE 'CMNCMMACH'
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL,
     +   AT, CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
      INTEGER I
      DOUBLE PRECISION P, PT, PX, PY, W, WT, WX, WY,
     +   RHO, RHOX, RHOY,
     +   MU, MUX, MUY, KAPMU, KAPMU2, KAPMUX, KAPMUY, Q1, Q2, QL,
     +   ND11, ND12, ND22,
     +   PXX, PXY, PYY, WXX, WXY, WYY,
     +   ND11Q1, ND11Q2, ND12Q1, ND12Q2, ND22Q1, ND22Q2,
     +   Q1X, Q1Y, Q2X, Q2Y,
     +   ND11X, ND12X, ND12Y, ND22Y, JW1, JW2, JW1X, JW2Y
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
         RHO    = RHO0*EXP(BETA*(P-P0)+GAMMA*W)
         RHOX   = RHO*(BETA*PX+GAMMA*WX)
         RHOY   = RHO*(BETA*PY+GAMMA*WY)
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
         PXX    = UXX(I,1)
         PXY    = UXY(I,1)
         PYY    = UYY(I,1)
         WXX    = UXX(I,2)
         WXY    = UXY(I,2)
         WYY    = UYY(I,2)
         ND11Q1 = (AT + (AL-AT)*(2-(Q1/QL)**2))*Q1/QL
         ND11Q2 = (AT - (AL-AT)*((Q1/QL)**2))*Q2/QL
         ND12Q1 = (AL-AT)*(Q2/QL)**3
         ND12Q2 = (AL-AT)*(Q1/QL)**3
         ND22Q1 = (AT - (AL-AT)*((Q2/QL)**2))*Q1/QL
         ND22Q2 = (AT + (AL-AT)*(2-(Q2/QL)**2))*Q2/QL
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
C
         RES(I,1) = N*RHO*(BETA*PT+GAMMA*WT) +
     +      RHOX*Q1+RHO*Q1X + RHOY*Q2+RHO*Q2Y
         RES(I,2) = N*RHO*WT +
     +      RHO*Q1*WX + RHO*Q2*WY +
     +      RHOX*JW1+RHO*JW1X + RHOY*JW2+RHO*JW2Y
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
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), U(NPTS,NPDE),
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
C LLBND  : IN. (0:LLBND(0)+1)
C          LLBND(0) = NBNDS: total # physical boundaries in actual grid
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in the list + 1
C          NB. corners with 2 different types of boundaries should be
C          pointed at twice.
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Dirichlet
C                     2: Lower boundary -I
C                     3: Left  boundary  I
C                     4: Upper boundary  I max. first order derivative
C                     5: Right boundary -I
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
      DOUBLE PRECISION UROUND, XMIN
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
      DOUBLE PRECISION N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL,
     +   AT, CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      COMMON /PROBLM/ N, GAMMA, MU0, RHO0, P0, W0, G, DM, KAPPA, AL, AT,
     +   CF, TKAPPA, LT, LL, CS, RHOS, T0, ALPHA, BETA, QC, TC
      SAVE /PROBLM/
C
      INTEGER I, J, K
      DOUBLE PRECISION P, PY, W, RHO, MU, KAPMU, Q2
C
      J = 1
C
C yL boundary: q2 = 0
C              dw/dy = 0                0<x<1/11, 2/11<x<1
C              q2 = qc
C              dw/dy = 0                1/11<=x<=2/11
C
            DO 10 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF ((     0 .LT. X(I) .AND. X(I) .LT. 1.0/11) .OR.
     +             (2.0/11 .LT. X(I) .AND. X(I) .LT. 1)) THEN
                  P      = U(I,1)
                  PY     = UY(I,1)
                  W      = U(I,2)
                  RHO    = RHO0*EXP(BETA*(P-P0)+GAMMA*W)
                  MU     = MU0*(1+1.85*W-4.0*W*W)
                  KAPMU  = KAPPA/MU
                  Q2     = -KAPMU*(PY+RHO*G)
                  RES(I,1) = Q2
                  RES(I,2) = UY(I,2)
               ELSE
                  P      = U(I,1)
                  PY     = UY(I,1)
                  W      = U(I,2)
                  RHO    = RHO0*EXP(BETA*(P-P0)+GAMMA*W)
                  MU     = MU0*(1+1.85*W-4.0*W*W)
                  KAPMU  = KAPPA/MU
                  Q2     = -KAPMU*(PY+RHO*G)

                  RES(I,1) = Q2 - QC
                  RES(I,2) = W - W0
               ENDIF
   10       CONTINUE
C
      J = 7
C
C yU boundary: p(x,0) = p0
C              dw/dy = 0                xL < x < xR
C
            DO 20 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = U(I,1) - P0
               RES(I,2) = UY(I,2)
   20       CONTINUE
C
      DO 30 J = 2, 8, 2
C
C xL boundary. dp/dx = 0
C xR boundary` dw/dx = 0
C
            DO 40 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = UX(I,1)
               RES(I,2) = UX(I,2)
   40       CONTINUE
   30 CONTINUE
C
      DO 50 J = 3, 5, 2
C
C `Interface conditions' : dp/dy = -rho.g      0 <= x <= 0.5
C                          dw/dy = 0            y = 0.4, 0.6
C
            DO 60 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               P      = U(I,1)
               PY     = UY(I,1)
               W      = U(I,2)
               RHO    = RHO0*EXP(BETA*(P-P0)+GAMMA*W)
               RES(I,1) = UY(I,1) + RHO*G
               RES(I,2) = UY(I,2)
   60       CONTINUE
   50 CONTINUE
C
      DO 70 J = 9, 16
C
C Corners: dp/dx = 0
C          dw/dx = 0
C
            DO 80 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               RES(I,1) = UX(I,1)
               RES(I,2) = UX(I,2)
   80       CONTINUE
   70 CONTINUE
C
      RETURN
      END
