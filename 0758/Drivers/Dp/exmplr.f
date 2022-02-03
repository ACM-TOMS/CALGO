      PROGRAM EXMPLR
C
C Restart of EXMPL, default values, Jacobian derivatives exact
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK, LENLWK
      PARAMETER (MXLEV=2, NPD=2, NPTS=2500)
      PARAMETER (LENIWK=NPTS*(5*MXLEV+14),
     +           LENRWK=NPTS*NPD*(5*MXLEV + 9+18*NPD),
     +           LENLWK=2*NPTS)
C
C-----------------------------------------------------------------------
C
      INTEGER LUNDMP
      PARAMETER (LUNDMP = 89)
C
      INTEGER NPDE, INFO(1), IWK(LENIWK), MNTR
      LOGICAL LWK(LENLWK)
      DOUBLE PRECISION T, TOUT, DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, RINFO(1), RWK(LENRWK)

C Continuation call of VLUGR2
      MNTR = 1
      TOUT = 3.0
      TOLS = 0.1
      TOLT = 0.05
C Default choices
      INFO(1) = 0
C
      OPEN(UNIT=LUNDMP,FILE='DUMP',FORM='UNFORMATTED')
      CALL RDDUMP (LUNDMP, RWK, LENRWK, IWK, LENIWK)
      CLOSE(LUNDMP)
C
C call main routine
      CALL  VLUGR2 (NPDE, T, TOUT, DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK, LWK, LENLWK,
     +   MNTR)
      PRINT *, 'VLUGR2 returned with MNTR=', MNTR
C
      OPEN(UNIT=LUNDMP,FILE='DUMP2',FORM='UNFORMATTED')
      CALL DUMP (LUNDMP, RWK, IWK)
      CLOSE(LUNDMP)
      END
      SUBROUTINE DERIVF (F, T, X, Y, NPTS, NPDE, U, A0, DT, DX, DY,
     +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UXX, UXY, UYY,
     +   ABSTOL, DEL, WORK,
     +   FU, FUX, FUY, FUXX, FUXY, FUYY)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION F(NPTS,NPDE), T, X(NPTS), Y(NPTS), U(NPTS,NPDE),
     +   A0, DT, DX, DY, UIB(*),
     +   UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE),
     +   UXX(NPTS,NPDE), UXY(NPTS,NPDE), UYY(NPTS,NPDE),
     +   ABSTOL(NPDE), DEL(NPTS), WORK(2*NPTS*NPDE),
     +   FU(NPTS,NPDE,NPDE), FUX(NPTS,NPDE,NPDE), FUY(NPTS,NPDE,NPDE),
     +   FUXX(NPTS,NPDE,NPDE),FUXY(NPTS,NPDE,NPDE),FUYY(NPTS,NPDE,NPDE)
C
Ccc PURPOSE:
C Compute derivatives of residual wrt (derivatives of) U
C
C PARAMETER DESCRIPTION:
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y    : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C U      : IN. Solution at T on current grid
C A0     : IN. Coefficient of U_n+1 in time derivative
C DT     : IN. Current time step size
C DX     : IN. Cell width in X-direction for current grid
C DY     : IN. Cell width in Y-direction for current grid
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
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UXX    : IN.  I Space derivatives of U on current grid
C UXY    : IN.  I
C UYY    : IN. -I
C ABSTOL : IN. Absolute tolerance for Newton process
C DEL    : WORK. (NPTS)
C WORK   : WORK. (2.LENU)
C FU     : OUT. dF(U,Ut)dU
C FUX    : OUT. dF(Ux)dUx
C FUY    : OUT. dF(Uy)dUy
C FUXX   : OUT. dF(Uxx)dUxx
C FUXY   : OUT. dF(Uxy)dUxy
C FUYY   : OUT. dF(Uyy)dUyy
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-3)
C
      INTEGER IC, IPT, LB, NBNDS
C
Ccc Loop over the components of the (derivatives of) U
      IC = 1
C
C dF(U,Ut)/dU_ic
         DO 20 IPT = 1, NPTS
            FU(IPT,1,IC) = A0 - (-UX(IPT,1))
            FU(IPT,2,IC) = - (-UX(IPT,2))
   20    CONTINUE
C
C dF(Ux)/dUx_ic
         DO 40 IPT = 1, NPTS
            FUX(IPT,1,IC) = - (-U(IPT,1))
            FUX(IPT,2,IC) = 0.0
   40    CONTINUE
C
C dF(Uy)/dUy_ic
         DO 50 IPT = 1, NPTS
            FUY(IPT,1,IC) = - (-U(IPT,2))
            FUY(IPT,2,IC) = 0.0
   50    CONTINUE
C
C dF(Uxx)/dUxx_ic
         DO 60 IPT = 1, NPTS
            FUXX(IPT,1,IC) = - (EPS)
            FUXX(IPT,2,IC) = 0.0
   60    CONTINUE
C
C dF(Uxy)/dUxy_ic
         DO 70 IPT = 1, NPTS
            FUXY(IPT,1,IC) = 0.0
            FUXY(IPT,2,IC) = 0.0
   70    CONTINUE
C
C dF(Uyy)/dUyy_ic
         DO 80 IPT = 1, NPTS
            FUYY(IPT,1,IC) = - (EPS)
            FUYY(IPT,2,IC) = 0.0
   80    CONTINUE
      
      IC = 2
C
C dF(U,Ut)/dU_ic
         DO 120 IPT = 1, NPTS
            FU(IPT,1,IC) = - (-UY(IPT,1))
            FU(IPT,2,IC) = A0 - (-UY(IPT,2))
  120    CONTINUE
C
C dF(Ux)/dUx_ic
         DO 140 IPT = 1, NPTS
            FUX(IPT,1,IC) = 0.0
            FUX(IPT,2,IC) = - (-U(IPT,1))
  140    CONTINUE
C
C dF(Uy)/dUy_ic
         DO 150 IPT = 1, NPTS
            FUY(IPT,1,IC) = 0.0
            FUY(IPT,2,IC) = - (-U(IPT,2))
  150    CONTINUE
C
C dF(Uxx)/dUxx_ic
         DO 160 IPT = 1, NPTS
            FUXX(IPT,1,IC) = 0.0
            FUXX(IPT,2,IC) = - (EPS)
  160    CONTINUE
C
C dF(Uxy)/dUxy_ic
         DO 170 IPT = 1, NPTS
            FUXY(IPT,1,IC) = 0.0
            FUXY(IPT,2,IC) = 0.0
  170    CONTINUE
C
C dF(Uyy)/dUyy_ic
         DO 180 IPT = 1, NPTS
            FUYY(IPT,1,IC) = 0.0
            FUYY(IPT,2,IC) = - (EPS)
  180    CONTINUE
      
C
C Correct boundaries (incl. the internal)
      NBNDS = LLBND(0)
      DO 100 LB = LLBND(1), LLBND(NBNDS+2)-1
         IPT = LBND(LB)
         FU(IPT,1,1) = 1.0
         FU(IPT,1,2) = 0.0
         FU(IPT,2,1) = 0.0
         FU(IPT,2,2) = 1.0
         FUX(IPT,1,1) = 0.0
         FUX(IPT,1,2) = 0.0
         FUX(IPT,2,1) = 0.0
         FUX(IPT,2,2) = 0.0
         FUY(IPT,1,1) = 0.0
         FUY(IPT,1,2) = 0.0
         FUY(IPT,2,1) = 0.0
         FUY(IPT,2,2) = 0.0
         FUXX(IPT,1,1) = 0.0
         FUXX(IPT,1,2) = 0.0
         FUXX(IPT,2,1) = 0.0
         FUXX(IPT,2,2) = 0.0
         FUXY(IPT,1,1) = 0.0
         FUXY(IPT,1,2) = 0.0
         FUXY(IPT,2,1) = 0.0
         FUXY(IPT,2,2) = 0.0
         FUYY(IPT,1,1) = 0.0
         FUYY(IPT,1,2) = 0.0
         FUYY(IPT,2,1) = 0.0
         FUYY(IPT,2,2) = 0.0
  100 CONTINUE

      RETURN
      END
      SUBROUTINE MONITR (T, DT, DTNEW, XL, YL, DXB, DYB,
     +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
      DOUBLE PRECISION T, DT, DTNEW, XL, YL, DXB, DYB, SOL(*)
C
Ccc PURPOSE:
C Control after a successful time step. The solution can be printed,
C plotted or compared with the exact solution.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN.  Current value of time variable
C DT     : IN.  Current time step size
C DTNEW  : IN.  Time step size for next time step
C XL     : IN.  X-coordinate of lowerleft corner of (virtual) rectangle
C YL     : IN.  Y-coordinate of lowerleft corner of (virtual) rectangle
C DXB    : IN.  Cell width in X-direction of base grid
C DYB    : IN.  Cell width in Y-direction of base grid
C LGRID  : IN.  (0:*)
C          LGRID(0) = max. grid level used at T
C          LGRID(1): pointer to base grid structure ISTRUC
C          LGRID(LEVEL): pointer to grid structure (LROW, IROW, ICOL)
C                        of refinement level LEVEL for time T
C ISTRUC : IN.  (*)
C          ISTRUC(LGRID(LEVEL):.) contains (LROW,IROW,ICOL) of grid
C                                 level LEVEL,
C          LROW   : (0:LROW(0)+1)
C             LROW(0) = NROWS: Actual # rows in grid
C             LROW(1:NROWS): pointers to the start of a row in the grid
C             LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C          IROW   : (NROWS)
C             IROW(IR): row number of row IR in virtual rectangle
C          ICOL   : (NPTS)
C             ICOL(IPT): column number of grid point IPT in virtual
C                        rectangle
C LSOL   : IN.  (*)
C          LSOL(LEVEL): pointer to (injected) solution at grid
C                       of refinement level LEVEL for time T
C SOL    : IN.  (*)
C          SOL(LSOL(LEVEL)+1:LSOL(LEVEL)+NPTS(LEVEL)*NPDE) contains
C          U_LEVEL(NPTS,NPDE)
C
C
C Local arrays:
      INTEGER MAXPTS, NPDE
      PARAMETER (MAXPTS=10000, NPDE=2)
      DOUBLE PRECISION X(MAXPTS), Y(MAXPTS), UEX(MAXPTS*NPDE)
C
C-----------------------------------------------------------------------
C
      INTEGER MAXLEV, LEVEL, LLROW, NROWS, NPTS, LIROW, LICOL
      DOUBLE PRECISION DX, DY
C
C Loop over the grid levels from coarse to fine.
C Get physical coordinates of grid points
C Compute ||err||_max
      MAXLEV = LGRID(0)
      DX = DXB
      DY = DYB
      DO 10 LEVEL = 1, MAXLEV
         LLROW  = LGRID(LEVEL)
         NROWS  = ISTRUC(LLROW)
         NPTS   = ISTRUC(LLROW+NROWS+1)-1
         LIROW  = LLROW+NROWS+2
         LICOL  = LIROW+NROWS
         CALL SETXY (XL, YL, DX, DY,
     +      ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL), X, Y)
         DX = DX/2
         DY = DY/2
         CALL PRERR (LEVEL, T, NPTS, NPDE, X, Y, SOL(LSOL(LEVEL)+1),
     +      UEX)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE PRERR (LEVEL, T, NPTS, NPDE, X, Y, U, UEX)
      INTEGER LEVEL, NPTS, NPDE
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), U(NPTS,NPDE), UEX(NPTS,NPDE)

      INTEGER I,J
      DOUBLE PRECISION RMAX(2)

      CALL PDEIV (T, X, Y, UEX, NPTS, NPDE)
      RMAX(1) = 0.0
      RMAX(2) = 0.0
      DO 10 I = 1, NPTS
         J = 1
         RMAX(J) = MAX(RMAX(J),ABS(UEX(I,J)-U(I,J)))
         J = 2
         RMAX(J) = MAX(RMAX(J),ABS(UEX(I,J)-U(I,J)))
   10 CONTINUE
      PRINT '(''Error at T='',E9.3,'', level='',I1,'' :'',2E11.3,I10)',
     +   T, LEVEL, RMAX(1), RMAX(2), NPTS

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
C Burgers' equation, solution wave front at y = x+0.25t, speed of
C propagation sqrt(2)/8 perpendicular to wave front.
C    U = 3/4 - 1/4/(1+exp((-4x+4y-t)/(32*eps)))
C    V = 3/4 + 1/4/(1+exp((-4x+4y-t)/(32*eps)))
C
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-3)
      INTEGER I

      DO 10 I = 1, NPTS
         U(I,1) = 0.75 - 0.25/(1+EXP((-4*X(I)+4*Y(I)-T)/(32*EPS)))
         U(I,2) = 0.75 + 0.25/(1+EXP((-4*X(I)+4*Y(I)-T)/(32*EPS)))
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
C-----------------------------------------------------------------------
C
C Burgers' equation Ut = - U.Ux - V.Uy + eps.(Uxx + Uyy)
C                   Vt = - U.Vx - V.Vy + eps.(Vxx + Vyy)
C
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-3)
      INTEGER I

      DO 10 I = 2, NPTS-1
         RES(I,1) = UT(I,1) -
     +      (-U(I,1)*UX(I,1) - U(I,2)*UY(I,1) + EPS*(UXX(I,1)+UYY(I,1)))
         RES(I,2) = UT(I,2) -
     +      (-U(I,1)*UX(I,2) - U(I,2)*UY(I,2) + EPS*(UXX(I,2)+UYY(I,2)))
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
C-----------------------------------------------------------------------
C
C Burgers' equation, Dirichlet boundaries.
C    U = 3/4 - 1/4/(1+exp((-4x+4y-t)/(32*eps)))
C    V = 3/4 + 1/4/(1+exp((-4x+4y-t)/(32*eps)))
C
      DOUBLE PRECISION EPS
      PARAMETER (EPS = 1D-3)
      INTEGER I, K, NBNDS

      NBNDS = LLBND(0)
      DO 10 K = LLBND(1), LLBND(NBNDS+1)-1
         I = LBND(K)
         RES(I,1) = U(I,1) -
     +            (0.75 - 0.25/(1+EXP((-4*X(I)+4*Y(I)-T)/(32*EPS))))
         RES(I,2) = U(I,2) -
     +            (0.75 + 0.25/(1+EXP((-4*X(I)+4*Y(I)-T)/(32*EPS))))
   10 CONTINUE

      RETURN
      END
