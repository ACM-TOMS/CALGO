      SUBROUTINE VLUGR3 (NPDE, T, TOUT, DT,
     +   XL,YF,ZD, XR,YB,ZU, DX, DY, DZ,
     +   TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK, LWK, LENLWK,
     +   MNTR)
C
C=======================================================================
C
Ccc PURPOSE:
C==========
C This code solves systems of PDEs of the type
C          F(t,x,y,z,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz)=0
C with boundary conditions
C          B(t,x,y,z,U,Ut,Ux,Uy,Uz)=0
C and initial values
C          U(t0,x,y,z)=U0
C on a 3D domain bounded by right-angled polyhedrons.
C
C In space Local Uniform Grid Refinement is applied to resolve local
C sharp gradients in the solution. For the time integration the
C implicit BDF2 method is used with variable stepsizes.
C Although time-independent and hyperbolic PDEs fit into the problem
C class, it should be observed that VLUGR3 is tuned for time-dependent
C parabolic PDEs (see below `HOW TO REPLACE MODULES' and the part on
C INCLUDEd files for the (non)linear solvers).
C
C
C
Ccc PARAMETER SPECIFICATION:
C==========================
      INTEGER LENIWK
      INTEGER NPDE, INFO(*), LENRWK, IWK(LENIWK), LENLWK, MNTR
      LOGICAL LWK(LENLWK)
      DOUBLE PRECISION T, TOUT, DT,
     +   XL,YF,ZD, XR,YB,ZU, DX, DY, DZ,
     +   TOLS, TOLT, RINFO(*), RWK(LENRWK)
C
Ccc LANGUAGE: FORTRAN 77
C===========
C
Ccc TYPE: Single precision
C=======
C
Ccc REFERENCE:
C============
C VLUGR3: A Vectorizable Adaptive Grid Solver for PDEs in 3D
C    Part I. Algorithmic Aspects and Applications
C J.G. Blom and J.G. Verwer, Applied Numerical Mathematics, Vol.16
C    pp.129-156 (1994).
C
C VLUGR3: A Vectorizable Adaptive Grid Solver for PDEs in 3D
C    Part II. Code Description
C J.G. Blom and J.G. Verwer, Report NM-R9405, CWI, Amsterdam.
C (to appear in ACM TOMS)
C
C
C
Ccc PARAMETER DESCRIPTION:
C========================
C NPDE   : IN.  # PDE components.
C T      : INOUT. Current value of time variable
C          IN:  If this is the first call the initial time
C          OUT: Time to which PDE has been integrated
C TOUT   : IN.  Time point at which solution is desired
C DT     : INOUT.
C          IN:  If this is the first call the initial time stepsize
C          OUT: Stepsize for next time step
C XL     : IN.  If this is the first call and INFO(3) = 0
C               X-coordinate of left/front/down point of rectangular
C               prism
C YF     : IN.  If this is the first call and INFO(3) = 0
C               Y-coordinate of left/front/down point of rectangular
C               prism
C ZD     : IN.  If this is the first call and INFO(3) = 0
C               Z-coordinate of left/front/down point of rectangular
C               prism
C XR     : IN.  If this is the first call and INFO(3) = 0
C               X-coordinate of right/back/upper point of rectangular
C               prism
C YB     : IN.  If this is the first call and INFO(3) = 0
C               Y-coordinate of right/back/upper point of rectangular
C               prism
C ZU     : IN.  If this is the first call and INFO(3) = 0
C               Z-coordinate of right/back/upper point of rectangular
C               prism
C DX     : IN.  If this is the first call and INFO(3) = 0
C               Cell width in X-direction of base grid
C DY     : IN.  If this is the first call and INFO(3) = 0
C               Cell width in Y-direction of base grid
C DZ     : IN.  If this is the first call and INFO(3) = 0
C               Cell width in Z-direction of base grid
C TOLS   : IN.  Space tolerance
C TOLT   : IN.  Time tolerance
C INFO   : IN.  If INFO(1)=0, default parameters are used, otherwise
C RINFO  : IN.  they should be specified in INFO and RINFO array
C               (for description see below)
C RWK    : WORK. (LENRWK)
C LENRWK : IN.  Dimension of RWK. (6.NPDE for VLUGR3)+:
C          Let NPTS be the max. # points on a grid level and
C          NPTSA the average # points over all grid levels.
C          Then LENRWK should be:
C          MAXLEV=1: 3.NPTS.NPDE+3.NPTS+13.NPTS.NPDE + LSSWRK
C             LSSWRK:
C                ( INFO(4)=0
C                | 38.NPDE.NPTS.NPDE
C                !:INFO(4)=10
C                | (MAX(NPDE.7+3,2.MAXLR+MAXL+6)+NPDE).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=11
C                | (MAX(NPDE.4+3,2.MAXLR+MAXL+6)+NPDE).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=12
C                | (2.MAXLR+MAXL+7).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=13
C                | (2.MAXLR+MAXL+7).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                )
C             (default: MAXLR = 5, MAXL = 20)
C          Indication of the length for a maximum grid level
C          MAXLEV (default value MAXLEV=3):
C             5.NPTSA.NPDE.MAXLEV+(3+13.NPDE).NPTS + LSSWRK
C IWK    : WORK. (LENIWK)
C LENIWK : IN.  Dimension of IWK. (8.MAXLEV+3 for VLUGR3)+:
C          MAXLEV=1: 28.NPTS
C          Indication of the length for a maximum grid level MAXLEV,
C          7.NPTSA.MAXLEV+7.NPTS + ( INFO(4)=0| 19.NPTS )
C LWK    : WORK. (LENLWK)
C LENLWK : IN.  Dimension of LWK. Indication of the length
C          2.NPTS
C MNTR   : INOUT. Monitor of VLUGR3
C          IN: State of integration
C           0. First call
C           1. Continuation call
C          OUT: Error return flag
C           1. OK
C          -1. Workspace too small
C          -2. Time step size too small
C         -10. COMMON to keep the statistics is too small
C
C
C
Ccc HOW TO USE: Default case
C===========================
C
C 3 problem defining routines should be specified
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE PDEIV (T, X, Y, Z, U, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER NPTS, NPDE
C     REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE)
C
Ccc PURPOSE:
C Define (initial) solution of PDE.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which (initial) solution should be given
C X,Y,Z  : IN. Arrays of physical coordinates for the gridpoints
C U      : OUT. Array of PDE component values for the gridpoints.
C NPTS   : IN. Number of gridpoints
C NPDE   : IN. # PDE components
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE PDEF (T, X, Y, Z, U,
C    +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, RES, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER NPTS, NPDE
C     REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
C    +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE), UZ(NPTS,NPDE),
C    +     UXX(NPTS,NPDE), UYY(NPTS,NPDE), UZZ(NPTS,NPDE),
C    +     UXY(NPTS,NPDE), UXZ(NPTS,NPDE), UYZ(NPTS,NPDE),
C    +     RES(NPTS,NPDE)
C
Ccc PURPOSE:
C Define residual of PDE on interior of domain. Boundary values will be
C overwritten later on.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which residual should be evaluated
C X,Y,Z  : IN. Physical coordinates of gridpoints
C U      : IN. Array of PDE components for the gridpoints.
C UT     : IN. Array of time derivative of PDE components
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C RES    : OUT. Array containg PDE residual at gridpoints in interior of
C              domain. The residual values at boundary points will be
C              overwritten by a call to PDEBC.
C NPTS   : IN. Number of gridpoints
C NPDE   : IN. Number of PDE components
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE PDEBC (T, X, Y, Z, U, UT, UX, UY, UZ, RES, NPTS, NPDE,
C    +   LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
C     REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
C    +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE), UZ(NPTS,NPDE),
C    +     RES(NPTS,NPDE)
C
Ccc PURPOSE:
C Define residual of boundary equations of PDE. The residual on interior
C points has already been stored in RES.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Time at which BC's should be evaluated
C X,Y,Z  : IN. Physical coordinates of gridpoints
C U      : IN. Array of PDE components for the gridpoints.
C UT     : IN. Array of time derivative of PDE components
C UX     : IN. -I
C UY     : IN.  I Arrays containing space derivatives of PDE components
C UZ     : IN. -I
C RES    : INOUT.
C          IN: PDE residual for interior points (set by PDEF)
C          OUT: Array with PDE residual at physical boundary points
C               inserted
C NPTS   : IN. Number of grid components
C NPDE   : IN. Number of PDE components
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical horizontal planes in
C             actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBDPTS)
C          LBND(LB): pointer to boundary point in actual grid
C                    structure (as in X, Y, Z, and U)
C
C-----------------------------------------------------------------------
C
C
C
Ccc HOW TO USE: Extra's
C======================
C
C If INFO(1) <> 0 a number of parameters can be specified in INFO and
C RINFO that are described below. The parenthesized value is the
C default value.
C
C INFO(2)         : MAXLEV (3)
C                   maximum # grid levels allowed
C INFO(3)         : RCTDOM (0)
C                   If RCTDOM=0 the initial domain is a rectangular
C                   prism otherwise the user should specify a subroutine
C                   INIDOM to define the initial grid (see below)
C INFO(4)         : LINSYS (0)
C                   Linear system solver in use
C                    0: BiCGStab + ILU
C                   10: GCRO + Block-diagonal preconditioning
C                   11: GCRO + Block-diagonal preconditioning
C                              (neglecting first-order derivatives
C                               at the boundaries)
C                   12: GCRO + Diagonal preconditioning
C                   13: GCRO + Diagonal preconditioning
C                              (neglecting first-order derivatives
C                               at the boundaries)
C                   NB. 10-13 are matrix-free solvers
C INFO(5)         : LUNPDS (0)
C                   Logical Unit # of file for information on the
C                   integration history. If 0, only global information
C                   will be written on standard output.
C INFO(6)         : LUNNLS (0)
C                   Logical Unit # of file for information on the
C                   Newton process. If 0, no information will be
C                   written.
C INFO(7)         : LUNLSS (0)
C                   Logical Unit # of file for information on the
C                   linear system solver. If 0, no information will be
C                   written.
C
C RINFO(1)        : DTMIN (0.0)
C                   minimum time stepsize allowed
C RINFO(2)        : DTMAX (TOUT-T)
C                   maximum time stepsize allowed
C RINFO(3)        : UMAX ((1.0))
C                   approx. max. value of the PDE solution components.
C                   Used for scaling purposes
C RINFO(3+NPDE)   : SPCWGT ((1.0))
C                   weighting factor used in the space monitor to
C                   indicate the relative importance of a PDE
C                   component on the space monitor
C RINFO(3+2.NPDE) : TIMWGT ((1.0))
C                   weighting factor used in the time monitor to
C                   indicate the relative importance of a PDE
C                   component on the time monitor
C
C
C
C After each successful time step a subroutine MONITR is called.
C Default is an empty body, but it can be overloaded with
C-----------------------------------------------------------------------
C
C     SUBROUTINE MONITR (T, DT, DTNEW, XL, YL, ZD, DXB, DYB, DZB,
C    +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
C     REAL T, DT, DTNEW, XL, YL, ZD, DXB, DYB, DZB, SOL(*)
C
Ccc PURPOSE:
C Control after a successful time step. The solution can be printed,
C plotted or compared with the exact solution.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN.  Current value of time variable
C DT     : IN.  Current time step size
C DTNEW  : IN.  Time step size for next time step
C XL     : IN.  X-coordinate of left/front/down point of (virtual) box
C YF     : IN.  Y-coordinate of left/front/down point of (virtual) box
C ZD     : IN.  Z-coordinate of left/front/down point of (virtual) box
C DXB    : IN.  Cell width in X-direction of base grid
C DYB    : IN.  Cell width in Y-direction of base grid
C DZB    : IN.  Cell width in Z-direction of base grid
C LGRID  : IN.  (0:*)
C          LGRID(0) = max. grid level used at T
C          LGRID(1): pointer to base grid structure ISTRUC
C          LGRID(LEVEL): pointer to grid structure
C                        (LPLN,IPLN,LROW,IROW,ICOL) of refinement
C                        level LEVEL for time T
C ISTRUC : IN.  (*)
C          ISTRUC(LGRID(LEVEL):.) contains (LPLN,IPLN,LROW,IROW,ICOL)
C                                 of grid level LEVEL,
C          LPLN   : (0:LPLN(0)+1)
C             LPLN(0) = NPLNS: Actual # horizontal planes in LROW
C             LPLN(1:NPLNS): pointers to the start of a plane in LROW
C             LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C          IPLN   : (NPLNS)
C             IPLN(IP): plane number of plane IP in virtual box
C          LROW   : (NROWS+1)
C             LROW(1:NROWS): pointers to the start of a row in the grid
C             LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C          IROW   : (NROWS)
C             IROW(IR): row number of row IR in virtual box
C          ICOL   : (NPTS)
C             ICOL(IPT): column number of grid point IPT in virtual box
C LSOL   : IN.  (*)
C          LSOL(LEVEL): pointer to (injected) solution at grid
C                       of refinement level LEVEL for time T
C SOL    : IN.  (*)
C          SOL(LSOL(LEVEL)+1:LSOL(LEVEL)+NPTS(LEVEL)*NPDE) contains
C          U_LEVEL(NPTS,NPDE)
C
C-----------------------------------------------------------------------
C
C
C
C To force grid refinement at a specific point in space and time and
C on a specific level, one can overload the routine CHSPCM with
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE CHSPCM (T, LEVEL, NPTS, X, Y, Z, NPDE, U, SPCMON, TOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LEVEL, NPTS, NPDE
C     REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE), SPCMON(NPTS), TOL
C
Ccc PURPOSE:
C Force grid refinement.
C If for a node IPT SPCMON(IPT) > TOL the 64 surrounding cells will be
C refined.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN.  Current value of time variable
C LEVEL  : IN.  Current grid level
C NPTS   : IN.  Number of grid points at this level
C X,Y,Z  : IN.  Arrays of physical coordinates for the gridpoints
C NPDE   : IN.  Number of PDE components
C U      : IN.  Array of PDE components for the gridpoints
C SPCMON : INOUT.
C          IN:  Space monitor values as determined by VLUGR3
C          OUT: Changed to a value > TOL where refinement is required
C TOL    : IN.  Tolerance with which SPCMON will be compared
C
C-----------------------------------------------------------------------
C
C
C
C If the initial domain is not a rectangular prism one should specify
C the initial grid via the function INIDOM
C
C-----------------------------------------------------------------------
C
C     LOGICAL FUNCTION INIDOM (MAXPTS,
C    +   XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
C    +   LPLN, IPLN, LROW, IROW, ICOL, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER MAXPTS, LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
C    +   LLBND(0:*), ILBND(*), LBND(*)
C     REAL XL, YF, ZD, XR, YB, ZU, DX, DY, DZ
C
Ccc PURPOSE:
C Define initial domain. NB. Boundaries should consist of as many points
C as are necessary to employ second-order space discretization, i.e.,
C a boundary enclosing the internal part of the domain should not
C include less than 3 grid points in any coordinate direction including
C the corners. If Neumann boundaries are used the minimum is 4 since
C otherwise the Jacobian matrix will be singular.
C
C A (virtual) box is placed around the (irregular) domain.
C The left/front/down point of this box is (XL,YF,ZD) in physical
C coordinates and (0,0,0) in column, row, plane coordinates, resp..
C The right/back/upper point is (XR,YB,ZU) resp. (Nx,Ny,Nz), where
C Nx = (XR-XL)/DX, Ny = (YB-YF)/DY, and Nz = (ZU-ZD)/DZ.
C Only real grid points are stored.
C The coordinate values of the initial grid should be stored plane
C after plane and rowwise in LPLN, IPLN, LROW, IROW, ICOL.
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
C XL     : OUT. X-coordinate of left/front/down point of virtual box
C YF     : OUT. Y-coordinate of left/front/down point of virtual box
C ZD     : OUT. Z-coordinate of left/front/down point of virtual box
C XR     : OUT. X-coordinate of right/back/upper point of virtual box
C YB     : OUT. Y-coordinate of right/back/upper point of virtual box
C ZU     : OUT. Z-coordinate of right/back/upper point of virtual box
C DX     : OUT. Grid width in X-direction
C DY     : OUT. Grid width in Y-direction
C DZ     : OUT. Grid width in Z-direction
C LPLN   : OUT. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # horizontal planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : OUT. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : OUT. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : OUT. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : OUT. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LLBND  : OUT. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : OUT. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C                      structure
C
C-----------------------------------------------------------------------
C
C
C
C To store the exact partial derivatives of the residual F with respect
C to (the derivatives of) U.
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE DERIVF (F, T, X, Y, Z, NPTS, NPDE, U,
C    +   A0, DT, DX, DY, DZ,
C    +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
C    +   UXY, UXZ, UYZ, ATOL, DEL, WORK,
C    +   FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ, FUXY, FUXZ, FUYZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
C     REAL F(NPTS,NPDE), T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
C    +   A0, DT, DX, DY, DZ, UIB(*),
C    +   UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE), UZ(NPTS,NPDE),
C    +   UXX(NPTS,NPDE), UYY(NPTS,NPDE), UZZ(NPTS,NPDE),
C    +   UXY(NPTS,NPDE), UXZ(NPTS,NPDE), UYZ(NPTS,NPDE),
C    +   ATOL(NPDE), DEL(NPTS), WORK(2*NPTS*NPDE),
C    +   FUX(NPTS,NPDE,NPDE), FUY(NPTS,NPDE,NPDE), FUZ(NPTS,NPDE,NPDE),
C    +   FUXX(NPTS,NPDE,NPDE),FUYY(NPTS,NPDE,NPDE),FUZZ(NPTS,NPDE,NPDE),
C    +   FUXY(NPTS,NPDE,NPDE),FUXZ(NPTS,NPDE,NPDE),FUYZ(NPTS,NPDE,NPDE)
C
Ccc PURPOSE:
C Compute derivatives of residual wrt (derivatives of) U
C
C PARAMETER DESCRIPTION:
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C U      : IN. Solution at T on current grid
C A0     : IN. Coefficient of U_n+1 in time derivative
C DT     : IN. Current time step size
C DX     : IN. Cell width in X-direction for current grid
C DY     : IN. Cell width in Y-direction for current grid
C DZ     : IN. Cell width in Z-direction for current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C ATOL   : IN. Absolute tolerance for Newton process
C DEL    : WORK. (NPTS)
C WORK   : WORK. (2.NPTS.NPDE)
C FU     : OUT. dF(U,Ut)dU
C FUX    : OUT. dF(Ux)dUx
C FUY    : OUT. dF(Uy)dUy
C FUZ    : OUT. dF(Uz)dUz
C FUXX   : OUT. dF(Uxx)dUxx
C FUYY   : OUT. dF(Uyy)dUyy
C FUZZ   : OUT. dF(Uzz)dUzz
C FUXY   : OUT. dF(Uxy)dUxy
C FUXZ   : OUT. dF(Uxz)dUxz
C FUYZ   : OUT. dF(Uyz)dUyz
C
C-----------------------------------------------------------------------
C
C
C
Ccc `HANDY' ROUTINES:
C===================
C
C VLUGR3 contains some routines that facilitate the use of the
C data structure.
C
C
C
C To make a printout of the domain one has defined with INIDOM one
C can call PRDOM
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE PRDOM (LPLN, IPLN, LROW, IROW, ICOL,
C    +   LLBND, ILBND, LBND, IDOM, NX, NY, NZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
C    +   LLBND(0:*), ILBND(*), LBND(*), IDOM(0:*), NX, NY, NZ
C
Ccc PURPOSE:
C Print domain plane-wise. Internal points are .., external points XX,
C physical plane-boundary points their ILBND value. Edges are given
C both ILBND values, corners an explicated 2-character value, and
C internal boundary values II.
C
Ccc PARAMETER DESCRIPTION:
C See INIDOM
C
C-----------------------------------------------------------------------
C
C
C
C To get the X-,Y- and Z-coordinates corresponding with the grid points
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE SETXYZ (XL, YF, ZD, DX, DY, DZ,
C    +   LPLN, IPLN, LROW, IROW, ICOL, X, Y, Z)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*)
C     REAL XL, YF, ZD, DX, DY, DZ, X(*), Y(*), Z(*)
C
Ccc PURPOSE:
C Store X-, Y- and Z-coordinates of the grid points.
C
Ccc PARAMETER DESCRIPTION:
C See MONITR.
C NB. DX = DXB.2^(1-LEVEL); the same for DY and DZ.
C
C-----------------------------------------------------------------------
C
C
C
C To print the solution and the corresponding coordinate values at all
C grid levels
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE PRSOL (LUN, T, NPDE, XL, YF, ZD, DXB, DYB, DZB,
C    +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LUN, NPDE, LGRID(0:*), ISTRUC(*), LSOL(*)
C     REAL T, XL, YF, ZD, DXB, DYB, DZB, SOL(*)
C
Ccc PURPOSE:
C Print solution and coordinate values at all grid levels.
C
Ccc PARAMETER DESCRIPTION:
C LUN    : IN.  Logical unit number of print file
C NPDE   : IN.  # PDE components
C Others see MONITR.
C
C-----------------------------------------------------------------------
C
C
C
C To write to file the (interpolated) solution values on a uniform grid
C of a specified grid level and the maximum grid level used in each
C point
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE WRUNI (LUNS, LUNG, UNILEV,
C    +   T, NPDE, XL, YF, ZD, DXB, DYB, DZB, NXB, NYB, NZB,
C    +   LGRID, ISTRUC, LSOL, SOL, UNIFRM, NX, NY, NZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C     INTEGER LUNS, LUNG, UNILEV,
C    +   NPDE, NXB, NYB, NZB, LGRID(0:*), ISTRUC(*), LSOL(*), NX, NY, NZ
C     DOUBLE PRECISION T, XL, YF, ZD, DXB, DYB, DZB, SOL(*),
C    +   UNIFRM(0:NX,0:NY,0:NZ,NPDE)
C
Ccc PURPOSE:
C Write (interpolated) solution values at grid level UNILEV to file
C LUNS.
C Write maximum gridlevel used in each point to file LUNG.
C NB. The data will not be correct for a domain with holes in it with
C a size of the width of the base grid.
C
Ccc PARAMETER DESCRIPTION:
C LUNS   : IN.  Logical unit number of solution file
C LUNG   : IN.  Logical unit number of grid level file
C UNILEV : IN.  Maximum grid level to be used to generate the data
C NPDE   : IN.  # PDE components
C NXB,NYB,NZB: IN. # gridcells in X-, Y- and Z-direction, resp., on grid
C          of base level
C UNIFRM : WORK. (Interpolated) solution on level UNILEV / max. grid
C          level used.
C NX,NY,NZ: IN. # gridcells in X-, Y- and Z-direction, resp., on grid
C          of level UNILEV
C Others see MONITR.
C
C-----------------------------------------------------------------------
C
C
C
C To dump all necessary information for a restart on file.
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE DUMP (LUNDMP, RWK, IWK)
C
C-----------------------------------------------------------------------
C
C
C
C
C To read all necessary information for a restart from the dump file.
C
C-----------------------------------------------------------------------
C
C     SUBROUTINE RDDUMP (LUNDMP, RWK, LENRWK, IWK, LENIWK)
C
C-----------------------------------------------------------------------
C
C
C
Ccc HOW TO REPLACE MODULES:
C=========================
C
Ccc Space discretization.
C Replace the computation of the derivatives in subroutine DERIVS by the
C desired discretization.
C If the new space discretization uses a larger stencil than the
C implemented one (internally a central 19-point stencil and at the
C boundary a 3-point one sided), one should use as linear system
C solver the matrix-free GCRO variant (INFO(4)=10,11,12 or 13).
C Moreover, one should check whether the required grid points are
C available on the current grid level, e.g. using the x-, y-
C and z-coordinates of the grid points (see SETXYZ above).
C Note that the refinement strategy results in subgrids of at least
C 5 points in every coordinate direction.
C
Ccc Linear system solver.
C If the new solver is matrix-free:
C    rewrite the body of subroutine GCRO using the routines
C    MVDIFF to compute y=Ax
C    If the (block-)diagonal preconditioner is wanted, use the routine
C    BCKBDI to compute w=P^(-1).v
C    (copy the call used in GCRO and replace the vector arguments for
C    x, y, v, w, and, optionally for the workspace needed)
C otherwise, if the ILU preconditioner is to be used:
C    rewrite the body of subroutine BICGST using the routines
C    BCKSLV to compute v=P^(-1).v and
C    MVDIAG to compute y=Ax
C    (copy the call used in BICGST and replace the vector argument(s)).
C
C If a user-made preconditioner is wanted, one should adapt INTGRB
C (when the Jacobian is used) or INTGRC (for a matrix-free solver).
C The calls to JACPB and PINIT, resp. should be replaced by calls to
C the routine that computes the preconditioner. In BICGST and GCRO,
C resp., one should call one's own routine to compute w=P^(-1).v
C instead of BCKSLV and BCKBDI, resp..
C
C If extra workspace is needed, the easiest way is to declare it in
C the subroutine.
C
C
C
Ccc DESCRIPTION OF THE SETUP IN THE WORKARRAYS:
C=============================================
C
Ccc Datastructure for the solution at a grid level
C The solution is stored plane after plane, rowwise, one component
C vector after the other in
C     REAL U(0:NPTS*NPDE)
C The element U(0) is added because pointers to non-existing nodes point
C to 0.
C
Ccc Solutions from 3 different time levels have to be saved. For Tn-1
C only the injected one (U); for Tn the original solution (S) belonging
C to a specific grid, the injected solution (U), and the injected
C solution at the Tn+1 grid; and for Tn+1 the solution (S) and when
C finished the injected solution (U).
C
C The real work storage is set up as follows:
C First some method related arrays of length NPDE each: SPCTOL, TIMWGT,
C RELTOL, ABSTOL, RTOL, ATOL.
C From 6*NPDE+1 work storage for PDESOL where the array RWK starts with
C index 1. From there it will contain the following items:
C First the X-, Y- and Z- coordinates for the base grid: X(NPTSB),
C Y(NPTSB), and Z(NPTSB).
C From 3*NPTSB+1 the solutions are stored:
C First for Tn-1: U_i for i=LSGNM1(0),(-1),1
C Next for Tn: S_i for i=1,...,LSGN(0)
C              U_i for i=LSGN(0)-1,(-1),1
C Next for Tn+1: S_1
C                U_i(Tn) at grid LSGNP1(i) I
C                S_i(Tn+1)                 I for i=2,...,LEVEL
C                when refinement is finished:
C                U_i(Tn+1) for i=LSGNP1(0)-1,(-1),1
C After the solutions work storage is available for the (interpolated)
C solutions from Tn-1 at the current grid, the current X- and
C Y-coordinates, if necessary the (interpolated) solution values at the
C internal boundary, the initial solution at Tn+1 at the current grid
C (since the not updated solution of the old time level has to be used),
C and for the derivatives and the linear solver.
C
Ccc If the linear solver uses a Jacobian and an ILU preconditioner
C (INFO(4)=0) the Jacobian is stored as a block 19-diagonal matrix.
C If a second-order discretization is used at the boundary the extra
C information will be stored in one of the `mixed-derivative blocks'.
C Addressing is done with the use of pointers to off 3-diagonal blocks
C (cf. LLDG and LUDG below).
C For the incomplete LU the second-order discretization at the
C boundaries is replaced by a first order discretization, since a true
C block 19-diagonal matrix is required to apply the hyperplane method.
C The same block structure will be used as for the Jacobian itself.
C
C
C
Ccc Datastructure for the grid at the current grid level
C A (virtual) box is placed around the irregular domain.
C The left/front/down point of this box is (XL,YF,ZD) in physical
C coordinates and (0,0,0) in column, row, plane coordinates, resp..
C The right/back/upper point is (XR,YB,ZU) resp. (Nx,Ny,Nz), where
C Nx = (XR-XL)/DX, Ny = (YB-YF)/DY, and Nz = (ZU-ZD)/DZ.
C Only real grid points are stored, plane after plane, rowwise.
C
C     INTEGER ISTRUC(0:*)
C
Ccc ISTRUC contains the following arrays:
C LPLN   : (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # horizontal planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C                      structure
C LBLWY  : (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C
C   The next 6 arrays are only stored if INFO(4) = 0
C   The next 2 arrays are used for the Jacobian structure and its ILU
C LLDG   : (NPTS,-9:-2)
C          LLDG(IPT,-9): pointer to node Y-below Z-below
C                        or to node Z-below Z-below
C          LLDG(IPT,-8): pointer to node left of Z-below
C          LLDG(IPT,-7): pointer to node Z-below
C          LLDG(IPT,-6): pointer to node right of Z-below
C          LLDG(IPT,-5): pointer to node Y-above Z-below
C          LLDG(IPT,-4): pointer to node left of Y-below
C                        or to node Y-below Y-below
C          LLDG(IPT,-3): pointer to node Y-below
C          LLDG(IPT,-2): pointer to node right of Y-below
C                        or to node left of the node left
C LUDG   : (NPTS,2:9)
C          LUDG(IPT,2): pointer to node left of Y-above
C                       or to node right of the node right
C          LUDG(IPT,3): pointer to node Y-above
C          LUDG(IPT,4): pointer to node right of node Y-above
C                       or to node Y-above Y-above
C          LUDG(IPT,5): pointer to node Y-below Z-above
C          LUDG(IPT,6): pointer to node left of Z-above
C          LUDG(IPT,7): pointer to node Z-above
C          LUDG(IPT,8): pointer to node right of Z-above
C          LUDG(IPT,9): pointer to node Y-above Z-above
C                       or to node Z-above Z-above
C
C   The next 4 arrays are used to hold the data dependency lists
C   for the ILU factorization and the forward, resp. backward
C   sweep of the backsolve
C LSL    : (NPTS)
C          LSL(ISLPT): pointer to node in actual grid
C LLSL   : (0:LLSL(0))
C          LLSL(0) = # independent data dependency lists in ILU
C                    factorization and forward sweep
C          LLSL(1:LLSL(0)): pointers to the start of a list in LSL
C LSU    : (NPTS)
C          LSU(ISLPT): pointer to node in actual grid
C LLSU   : (0:LLSU(0))
C          LLSU(0) = # independent data dependency lists in backward
C                    sweep
C          LLSU(1:LLSU(0)): pointers to the start of a list in LSU
C
C For the base grid the complete datastructure is saved (including
C the last 6 arrays because of restart), for higher level grids
C only the first 5 arrays LPLN, IPLN, LROW, IROW and ICOL.
C
C Pointers to the specific arrays in ISTRUC are obtained by
C     LLPLN  = 0
C     NPLNS  = ISTRUC(LLPLN)
C     LIPLN  = LLPLN+NPLNS+2
C     LLROW  = LIPLN+NPLNS
C     NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
C     NPTS   = ISTRUC(LLROW+NROWS)-1
C     LIROW  = LLROW+NROWS+1
C     LICOL  = LIROW+NROWS
C     LLLBND = LICOL+NPTS
C     NBNDS  = ISTRUC(LLLBND)
C     NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
C     LILBND = LLLBND+NBNDS+3
C     LLBNDP = LILBND+NBNDS
C     LLBLWY = LLBNDP+NBIPTS
C     LLABVY = LLBLWY+NPTS+1
C     LLBLWZ = LLABVY+NPTS+1
C     LLABVZ = LLBLWZ+NPTS+1
C     LIWK   = LLABVZ+NPTS
C
C     LLLDG  = LIWK
C     LLUDG  = LLLDG+NPTS*8
C     LLSLP  = LLUDG+NPTS*8
C     LLLSL  = LLSLP+NPTS
C     LLSUP  = LLLSL+ISTRUC(LLLSL)+1
C     LLLSU  = LLSUP+NPTS
C     LIWK   = LLLSU+ISTRUC(LLLSU)+1
C
C
Ccc All grids from 3 different time levels have to be saved
C The integer work storage is set up as follows:
C LSGNM1 : (0:MAXLEV)
C          LSGNM1(0) = max. grid level used at Tn-1
C          LSGNM1(1): pointer to base grid structure ISTRUC
C          LSGNM1(LEVEL): pointer to grid structure
C                     (LPLN, IPLN, LROW, IROW, ICOL)
C                     of refinement level LEVEL for time Tn-1
C LSGN   : (0:MAXLEV)
C          LSGN(0) = max. grid level used at Tn
C          LSGN(1): pointer to base grid structure ISTRUC
C          LSGN(LEVEL): pointer to grid structure
C                     (LPLN, IPLN, LROW, IROW, ICOL)
C                     of refinement level LEVEL for time Tn
C LSGNP1 : (0:MAXLEV)
C          LSGNP1(0) = max. grid level used at Tn+1
C          LSGNP1(1): pointer to base grid structure ISTRUC
C          LSGNP1(2): pointer after grid structure of max. refinement
C                     level for time Tn
C          LSGNP1(LEVEL): pointer to augmented grid structure
C                     (LPLN, IPLN, LROW, IROW, ICOL, LLBND, ILBND, LBND)
C                     of refinement level LEVEL for time Tn+1
C          LSGNP1(LEVEL+1): pointer to grid structure ISTRUC of
C                     refinement level LEVEL+1 for time Tn+1
C LSUNM1 : (MAXLEV)
C          LSUNM1(LEVEL): pointer to (injected) solution at grid
C                     of refinement level LEVEL for time Tn-1
C LSSN   : (MAXLEV)
C          LSSN(LEVEL): pointer to original solution belonging
C                     to refinement level LEVEL for time Tn
C LSUN   : (MAXLEV)
C          LSUN(LEVEL): pointer to (injected) solution at grid
C                     of refinement level LEVEL for time Tn
C LSSNP1 : (MAXLEV)
C          LSSNP1(LEVEL): pointer to original solution belonging
C                     to refinement level LEVEL for time Tn+1
C LSUNP1 : (MAXLEV)
C          LSUNP1(LEVEL): pointer to (injected) solution at grid
C                     of refinement level LEVEL for time Tn+1
C From 8*MAXLEV+4 the grids are stored, in PDESOL the array IWK starts
C with the grids at index 1.
C Storage order:
C First ISTRUC for the base grid
C Next for Tn-1: (LPLN, IPLN, LROW, IROW, ICOL)_i for i=2,...,LSGNM1(0)
C Next for Tn:   (LPLN, IPLN, LROW, IROW, ICOL)_i for i=2,...,LSGN(0)
C Next for Tn+1: (LPLN, IPLN, LROW, IROW, ICOL, LLBND, ILBND, LBND)_i
C                                                 for i=2,...,LEVEL
C                ISTRUC_i for i=LEVEL+1
C After the grids work storage is available for domain flags and
C the linear solver
C
C=======================================================================
C
C IMPORTANT:
C ========= 
C
C The INCLUDEd file CMNCMMACH contains machine numbers that
C are set in the routine MACNUM by calling the appropriate functions
C of the BLAS library. If I1MACH and D1MACH of the file blas.f are used,
C the functions should be altered for the particular machine used (cf.
C comment in I1MACH and D1MACH).
C
Ccc CMNCMMACH
C
C COMMON with `machine numbers'
C LUNOUT : Logical unit # standard output         -I
C LUNERR : Logical unit # standard error           I  Set in the routine
C UROUND : Smallest machine number such that       I      MACNUM
C          1.0+UROUND > 1.0 and 1.0-UROUND < 1.0   I
C XMIN   : Smallest floating-point number > 0.0   -I
C     INTEGER LUNOUT, LUNERR
C     REAL UROUND, XMIN
C     COMMON /IMACH/ LUNOUT, LUNERR
C     COMMON /RMACH/ UROUND, XMIN
C     SAVE /IMACH/, /RMACH/
C
C
C
C The INCLUDE files PARNEWTON, PARBICGSTAB, and PARGCRO contain the
C method parameters for the corresponding (non)linear solvers. These
C parameters may be changed by the user.
C
Ccc PARNEWTON
C
C Parameters for Newton process
C MAXNIT : Max. number of Newton iterations
C          NB. If MAXNIT > 20 the include file CMNSTATS
C          ==  should also be changed.
C MAXJAC : Max. number of Jacobian / preconditioner evaluations during
C          a Newton process
C TOLNEW : Tolerance for Newton process:
C          rho/(1-rho)*|| corr.||_w < TOLNEW
C     INTEGER MAXNIT, MAXJAC
C     REAL TOLNEW
C     PARAMETER (MAXNIT = 10, MAXJAC = 2, TOLNEW = 1.0)
C
Ccc PARBICGSTAB
C
C Parameters for linear system solver BiCGStab
C MAXLIT : Max. number of BiCGStab iterations
C TOLLSB : Tolerance for linear system solver:
C          || P^(-1).residual ||_w < TOLLSB/2^NIT
C     INTEGER MAXLIT
C     REAL TOLLSB
C     PARAMETER (MAXLIT = 100, TOLLSB = TOLNEW/10)
C
Ccc PARGCRO
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
C          || P^(-1).residual ||_w < TOLLSC/2^NIT
C     INTEGER IDIAGP, NRRMAX, MAXLR, MAXL
C     REAL TOLLSC
C     PARAMETER (NRRMAX = 1, MAXLR = 5, MAXL = 20)
C     PARAMETER (TOLLSC = TOLNEW/10)
C     COMMON /IGCRO/ IDIAGP
C     SAVE /IGCRO/
C
C Note, that in the actual code the INCLUDE statements have been
C replaced by
C
CCcc   INCLUDE 'file'
C  ...  code in file
CC end INCLUDE 'file'
C
C So if one wishes to change the method parameters care should be taken
C that it is done for all occurrences.
C
C=======================================================================
C
Ccc EXTERNALS USED:
      EXTERNAL ICOPY, INTGRB, INTGRC, IYPOC, PDESOL, RCOPY
C
C
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
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
Ccc   INCLUDE 'CMNWRITEF'
C
C CMNWRITEF
C
C COMMON needed for continuation calls
      INTEGER MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      LOGICAL FIRST, SECOND
      DOUBLE PRECISION T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,
     +   DZB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW,
     +   DXB,DYB,DZB, DTO
      SAVE /WRITIF/, /WRITLF/, /WRITRF/
C
C end INCLUDE 'CMNWRITEF'
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
Ccc   INCLUDE 'PARBICGSTAB'
C
C PARBICGSTAB
C
C Parameters for linear system solver BiCGStab
C MAXLIT : Max. number of BiCGStab iterations
C TOLLSB : Tolerance for linear system solver
      INTEGER MAXLIT
      DOUBLE PRECISION TOLLSB
      PARAMETER (MAXLIT = 100, TOLLSB = TOLNEW/10)
C
C end INCLUDE 'PARBICGSTAB'
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
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 15)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
C
C
C-----------------------------------------------------------------------
C
      INTEGER LSGNM1, LSGN, LSGNP1, LSUNM1, LSUN, LSUNP1, LSSN, LSSNP1,
     +        LGNM1, LGN, LGNP1, LUNM1, LUN, LSN, LSNP1,
     +   MAXLEV, LSPCTL, LTIMWT, LRELTL, LABSTL, LRTOL, LATOL, LSPCWT,
     +   LUMAX, RCTDOM, LINSYS, I, I1, I2, IC, J, LRWK, LIWK, LIWKPN
      DOUBLE PRECISION DTMIN, DTMAX, TOL
C
CDIR$ NOVECTOR
      IF (MNTR .EQ. 1) THEN
         NPDE = NPDEW
         T    = TW
         DT   = DTW
         XL   = XLW
         YF   = YFW
         ZD   = ZDW
         XR   = XRW
         YB   = YBW
         ZU   = ZUW
      ENDIF
C
C Set machine numbers in /CMMACH/
      CALL MACNUM
C
C Setup real work storage
      LSPCTL = 1
      LTIMWT = LSPCTL+NPDE
      LRELTL = LTIMWT+NPDE
      LABSTL = LRELTL+NPDE
      LRTOL  = LABSTL+NPDE
      LATOL  = LRTOL +NPDE
      LRWKPS = LATOL +NPDE
      LSPCWT = LSPCTL
      LUMAX  = LATOL
C
C Get User info
      IF (INFO(1) .EQ. 0) THEN
         MAXLEV = 3
         RCTDOM = 0
         LINSYS = 0
         LUNPDS = 0
         LUNNLS = 0
         LUNLSS = 0
         DTMIN  = 0.0
         DTMAX  = TOUT - T
         DO 10 IC = 1, NPDE
            RWK(LUMAX-1 +IC) = 1.0
            RWK(LSPCWT-1+IC) = 1.0
            RWK(LTIMWT-1+IC) = 1.0
   10    CONTINUE
      ELSE
         MAXLEV = INFO(2)
         IF (MAXLEV .EQ. 0) MAXLEV = 3
         RCTDOM = INFO(3)
         IDIAGP = MOD(INFO(4),10)
         LINSYS = INFO(4)/10
         LUNPDS = INFO(5)
         LUNNLS = INFO(6)
         LUNLSS = INFO(7)
         DTMIN  = RINFO(1)
         DTMAX  = RINFO(2)
         IF (DTMAX .EQ. 0.0) DTMAX  = TOUT - T
         DO 20 IC = 1, NPDE
            RWK(LUMAX-1 +IC) = RINFO(2+IC)
            RWK(LSPCWT-1+IC) = RINFO(2+NPDE+IC)
            RWK(LTIMWT-1+IC) = RINFO(2+2*NPDE+IC)
   20    CONTINUE
      ENDIF
C
C Store method arrays
      TOL = 1D-1*MIN(TOLT*TOLT,TOLS)
      DO 30 IC = 1, NPDE
         RWK(LSPCTL-1+IC) = RWK(LSPCWT-1+IC)/(RWK(LUMAX-1+IC)*TOLS)
         RWK(LRELTL-1+IC) = TOLT
         RWK(LABSTL-1+IC) = 0.01*RWK(LUMAX-1+IC)*TOLT
         RWK(LRTOL-1+IC)  = TOL
         RWK(LATOL-1+IC)  = 0.01*RWK(LUMAX-1+IC)*TOL
   30 CONTINUE
C
C Setup integer work storage
      IF (MXCLEV .LT. MAXLEV) THEN
         WRITE(LUNERR,*) 'Arrays for the statistic are too small'
         WRITE(LUNERR,*) 'Either MAXLEV > 10 or MAXNIT > 20'
         WRITE(LUNERR,*) 'Adapt the parameter statements for /STATS/'
         MNTR = -10
         RETURN
      ENDIF
      LSGNM1 = 1
      LSGN   = LSGNM1 + MAXLEV+1
      LSGNP1 = LSGN   + MAXLEV+1
      LSUNM1 = LSGNP1 + MAXLEV+1
      LSSN   = LSUNM1 + MAXLEV
      LSUN   = LSSN   + MAXLEV
      LSSNP1 = LSUN   + MAXLEV
      LSUNP1 = LSSNP1 + MAXLEV
      LIWKPN = LSUNP1 + MAXLEV
      IF (MNTR .EQ. 0) THEN
C This is the first call, initialize pointer arrays and STATS common
         DO 50 I = 1, LIWKPN-1
            IWK(I) = 1
   50    CONTINUE
         NSTEPS = 0
         NREJS  = 0
         DO 60 I = 1, MXCLEV
            NJACS(I)  = 0
            NRESID(I) = 0
            NNIT(I)   = 0
            DO 70 J = 1, MXCNIT
               NLSIT(I,J) = 0
   70       CONTINUE
   60    CONTINUE
      ELSE IF (MAXLEV .GT. MAXLVW) THEN
C MAXLEV larger than previous call; shift info in IWK array backwards
         IF (LENIWK .LT. LIWKPN+LIWKB) THEN
            WRITE(LUNERR,*) 'Integer work space too small, required:',
     +         LIWKPN+LIWKB
            MNTR = -1
            RETURN
         ENDIF
         CALL IYPOC (LIWKB, IWK(LIWKPS), IWK(LIWKPN))
         LIWK   = LIWKPS - MAXLVW
         CALL IYPOC (MAXLVW, IWK(LIWK), IWK(LSUNP1))
         LIWK   = LIWK - MAXLVW
         CALL IYPOC (MAXLVW, IWK(LIWK), IWK(LSSNP1))
         LIWK   = LIWK - MAXLVW
         CALL IYPOC (MAXLVW, IWK(LIWK), IWK(LSUN))
         LIWK   = LIWK - MAXLVW
         CALL IYPOC (MAXLVW, IWK(LIWK), IWK(LSSN))
         LIWK   = LIWK - MAXLVW
         CALL IYPOC (MAXLVW, IWK(LIWK), IWK(LSUNM1))
         LIWK   = LIWK - MAXLVW-1
         CALL IYPOC (MAXLVW+1, IWK(LIWK), IWK(LSGNP1))
         LIWK   = LIWK - MAXLVW-1
         CALL IYPOC (MAXLVW+1, IWK(LIWK), IWK(LSGN))
      ELSE IF (MAXLEV .LT. MAXLVW) THEN
C MAXLEV smaller than previous call; shift info in IWK array forwards
         LGNM1 = 1
         LGN   = LGNM1 + MAXLVW+1
         LGNP1 = LGN   + MAXLVW+1
         LUNM1 = LGNP1 + MAXLVW+1
         LSN   = LUNM1 + MAXLVW
         LUN   = LSN   + MAXLVW
         LSNP1 = LUN   + MAXLVW
         IF (IWK(LGNM1) .GT. MAXLEV) THEN
C       Shift grid_n forwards to LGNM1(MAXLEV+1)
            I1  = IWK(LGN+2)
            I2  = IWK(LGNM1+MAXLEV+1)
            CALL ICOPY (LIWKB-I1, IWK(I1), IWK(I2))
            DO 110 I = 2, IWK(LGN)
               IWK(LGN+I) = IWK(LGN+I) - (I1-I2)
  110       CONTINUE
            LIWKB = LIWKB - (I1-I2)
C       Shift info from U_n-1(MAXLEV) forwards to LUNM1(LGNM1(0))
            I1 = IWK(LUNM1-1+MAXLEV)
            I2 = IWK(LUNM1-1+IWK(LGNM1))
            CALL RCOPY (LRWKB-I1, RWK(I1), RWK(I2))
            DO 120 I = 1, MAXLEV
               IWK(LUNM1-1+I) = IWK(LUNM1-1+I) - (I1-I2)
  120       CONTINUE
            DO 130 I = 1, IWK(LGN)
               IWK(LSN-1+I) = IWK(LSN-1+I) - (I1-I2)
               IWK(LUN-1+I) = IWK(LUN-1+I) - (I1-I2)
  130       CONTINUE
            IWK(LSNP1) = IWK(LSNP1) - (I1-I2)
            LRWKB = LRWKB - (I1-I2)
            IWK(LGNM1) = MAXLEV
         ENDIF
         IF (IWK(LGN) .GT. MAXLEV) THEN
            LIWKB = IWK(LGN+MAXLEV+1)
C       Shift info from U_n(MAXLEV) forwards to LUN(LGN(0))
            I1 = IWK(LUN-1+MAXLEV)
            I2 = IWK(LUN-1+IWK(LGN))
            CALL RCOPY (LRWKB-I1, RWK(I1), RWK(I2))
            DO 140 I = 1, MAXLEV
               IWK(LUN-1+I) = IWK(LUN-1+I) - (I1-I2)
  140       CONTINUE
            IWK(LSNP1) = IWK(LSNP1) - (I1-I2)
            LRWKB = LRWKB - (I1-I2)
            IWK(LGN) = MAXLEV
         ENDIF
C Shift pointer arrays and grids forwards
         CALL ICOPY (MAXLEV+1, IWK(LGN),   IWK(LSGN))
         CALL ICOPY (MAXLEV+1, IWK(LGNP1), IWK(LSGNP1))
         CALL ICOPY (MAXLEV,   IWK(LUNM1), IWK(LSUNM1))
         CALL ICOPY (MAXLEV,   IWK(LSN),   IWK(LSSN))
         CALL ICOPY (MAXLEV,   IWK(LUN),   IWK(LSUN))
         IWK(LSSNP1) = IWK(LSNP1)
         CALL ICOPY (LIWKB, IWK(LIWKPS), IWK(LIWKPN))
      ENDIF
      LIWKPS = LIWKPN

      IF (LUNPDS .NE. 0) THEN
         LUN = LUNPDS
      ELSE
         LUN = LUNOUT
      ENDIF
C
C Call main routine
      LRWK = LENRWK - LRWKPS+1
      LIWK = LENIWK - LIWKPS+1
      WRITE(LUN,*) 'Newton: MAXNIT, MAXJAC, TOLNEW:',
     +   MAXNIT, MAXJAC, TOLNEW
      IF (LINSYS .EQ. 0) THEN
C Lin. sys. solver = BiCGStab
      WRITE(LUN,*) 'Lin. solver BiCGStab + ILU: MAXLIT, TOLLSB:',
     +   MAXLIT, TOLLSB
      CALL PDESOL (MAXLEV, NPDE, IWK(LSGNM1), IWK(LSGN), IWK(LSGNP1),
     +   IWK(LSUNM1), IWK(LSSN), IWK(LSUN), IWK(LSSNP1), IWK(LSUNP1),
     +   T, TOUT, DT, DTMIN, DTMAX, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   RWK(LRTOL), RWK(LATOL), RWK(LSPCTL), RWK(LTIMWT),
     +   RWK(LRELTL), RWK(LABSTL),
     +   LINSYS, INTGRB,
     +   RWK(LRWKPS), LRWK, IWK(LIWKPS), LIWK, LWK, LENLWK, MNTR)
      ELSE IF (LINSYS .EQ. 1) THEN
C Lin. sys. solver = matrix-free GCRO
      IF (IDIAGP .LE. 1) THEN
         WRITE(LUN,*) 'Lin. solver matrix-free GCRO + Block-diag:',
     +      'NRRMAX, MAXLR, MAXL, TOLLSC:', NRRMAX, MAXLR, MAXL, TOLLSC
      ELSE
         WRITE(LUN,*) 'Lin. solver matrix-free GCRO + Diag:',
     +      'NRRMAX, MAXLR, MAXL, TOLLSC:', NRRMAX, MAXLR, MAXL, TOLLSC
      ENDIF
      CALL PDESOL (MAXLEV, NPDE, IWK(LSGNM1), IWK(LSGN), IWK(LSGNP1),
     +   IWK(LSUNM1), IWK(LSSN), IWK(LSUN), IWK(LSSNP1), IWK(LSUNP1),
     +   T, TOUT, DT, DTMIN, DTMAX, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   RWK(LRTOL), RWK(LATOL), RWK(LSPCTL), RWK(LTIMWT),
     +   RWK(LRELTL), RWK(LABSTL),
     +   LINSYS, INTGRC,
     +   RWK(LRWKPS), LRWK, IWK(LIWKPS), LIWK, LWK, LENLWK, MNTR)
       ENDIF
C
C Give final statistics
      IF (MNTR .NE. 0) THEN
         WRITE(LUN,'(''Error exit PDESOL, MNTR='',I4)') MNTR
      ELSE
         MNTR = 1
      ENDIF
      WRITE(LUN,*)
      WRITE(LUN,'(''Statistics:'')')
      WRITE(LUN,'('' # accepted timesteps ='', I5,
     +   '', # rejected timesteps ='', I5)') NSTEPS, NREJS
      WRITE(LUN,'('' Level # Nit  # Jacs   # Res'')')
      DO 200 I = 1, MXCLEV
      IF (NNIT(I) .NE. 0)
     +   WRITE(LUN,'(2I6,2I8)') I, NNIT(I), NJACS(I), NRESID(I)
  200 CONTINUE
      WRITE(LUN,'(''   Nit Level  # Lin. sys. it'')')
      DO 210 J = 1, MXCNIT
      DO 210 I = 1, MXCLEV
      IF (NLSIT(I,J) .NE. 0)
     +   WRITE(LUN,'(2I6,I12)') J, I, NLSIT(I,J)
  210 CONTINUE
C
C Take care of all information needed to dump info to file
      MAXLVW = MAXLEV
      NPDEW  = NPDE
      LRWKB  = IWK(LSSNP1)
      TW     = T
      TEW    = TOUT
      DTW    = DT
      XLW    = XL
      YFW    = YF
      ZDW    = ZD
      XRW    = XR
      YBW    = YB
      ZUW    = ZU

      RETURN
      END
      SUBROUTINE PDESOL (MAXLEV, NPDE, LSGNM1, LSGN, LSGNP1,
     +   LSUNM1, LSSN, LSUN, LSSNP1, LSUNP1,
     +   TN, TE, DT, DTMIN, DTMAX, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   RTOL, ATOL, SPCTOL, TIMWGT, RELTOL, ABSTOL,
     +   LINSYS, INTGRT,
     +   RWK, LENRWK, IWK, LENIWK, LWK, LENLWK, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LENIWK
      INTEGER MAXLEV, NPDE,
     +   LSGNM1(0:MAXLEV), LSGN(0:MAXLEV), LSGNP1(0:MAXLEV),
     +   LSUNM1(MAXLEV), LSSN(MAXLEV), LSUN(MAXLEV),
     +   LSSNP1(MAXLEV), LSUNP1(MAXLEV), LINSYS, 
     +   LENRWK, IWK(LENIWK), LENLWK, IERR
      LOGICAL LWK(LENLWK)
      DOUBLE PRECISION TN, TE, DT, DTMIN, DTMAX, XL, YF, ZD, XR, YB, ZU,
     +   DX, DY, DZ,
     +   RTOL(NPDE), ATOL(NPDE), SPCTOL(NPDE), TIMWGT(NPDE),
     +   RELTOL(NPDE), ABSTOL(NPDE), RWK(LENRWK)
      EXTERNAL INTGRT
C
Ccc PARAMETER DESCRIPTION:
C MAXLEV : IN.  Max. # grid levels allowed
C NPDE   : IN.  # PDE components.
C LSGNM1 : IN.  (0:MAXLEV)
C          LSGNM1(0) = max. grid level used at Tn-1
C          LSGNM1(1): pointer to base grid structure ISTRUC
C          LSGNM1(LEVEL): pointer to grid structure
C                         (LPLN,IPLN,LROW,IROW,ICOL)
C                         of refinement level LEVEL for time Tn-1
C LSGN   : IN.  (0:MAXLEV)
C          LSGN(0) = max. grid level used at Tn
C          LSGN(1): pointer to base grid structure ISTRUC
C          LSGN(LEVEL): pointer to grid structure
C                       (LPLN,IPLN,LROW,IROW,ICOL)
C                       of refinement level LEVEL for time Tn
C LSGNP1 : IN.  (0:MAXLEV)
C          LSGNP1(0) = max. grid level used at Tn+1
C          LSGNP1(1): pointer to base grid structure ISTRUC
C          LSGNP1(2): pointer after grid structure of max. refinement
C                     level for time Tn
C          LSGNP1(LEVEL): pointer to augmented grid structure
C                         (LPLN,IPLN,LROW,IROW,ICOL,LLBND,ILBND,LBND)
C                         of refinement level LEVEL for time Tn+1
C          LSGNP1(LEVEL+1): pointer to grid structure ISTRUC of
C                           refinement level LEVEL+1 for time Tn+1
C LSUNM1 : IN.  (MAXLEV)
C          LSUNM1(LEVEL): pointer to (injected) solution at grid
C                         of refinement level LEVEL for time Tn-1
C LSSN   : IN.  (MAXLEV)
C          LSSN(LEVEL): pointer to original solution belonging
C                       to refinement level LEVEL for time Tn
C LSUN   : IN.  (MAXLEV)
C          LSUN(LEVEL): pointer to (injected) solution at grid
C                       of refinement level LEVEL for time Tn
C LSSNP1 : IN.  (MAXLEV)
C          LSSNP1(LEVEL): pointer to original solution belonging
C                         to refinement level LEVEL for time Tn+1
C LSUNP1 : IN.  (MAXLEV)
C          LSUNP1(LEVEL): pointer to (injected) solution at grid
C                         of refinement level LEVEL for time Tn+1
C NB. All the above pointers should be initialized on 1
C TN     : INOUT. Current value of time variable
C          IN:  Initial time
C          OUT: Time to which PDE has been integrated
C TE     : IN.  Time point at which solution is desired
C DT     : INOUT.
C          IN:  The initial time stepsize
C          OUT: Stepsize for next time step
C DTMIN  : IN.  Minimum time stepsize allowed
C          If IERR=0 and domain a rectangular prism:
C DTMAX  : IN.  Maximum time stepsize allowed
C XL     : IN. X-coordinate of left/front/down point of rectangular
C              prism
C YF     : IN. Y-coordinate of left/front/down point of rectangular
C              prism
C ZD     : IN. Z-coordinate of left/front/down point of rectangular
C              prism
C XR     : IN. X-coordinate of right/back/upper point of rectangular
C              prism
C YB     : IN. Y-coordinate of right/back/upper point of rectangular
C              prism
C ZU     : IN. Z-coordinate of right/back/upper point of rectangular
C              prism
C DX     : IN. Cell width in X-direction of base grid
C DY     : IN. Cell width in Y-direction of base grid
C DZ     : IN. Cell width in Z-direction of base grid
C
C RTOL   : IN.  (NPDE)
C          Relative tolerance for the Newton iteration process
C ATOL   : IN.  (NPDE)
C          Absolute tolerance for the Newton iteration process
C SPCTOL : IN.  (NPDE)
C          Space tolerance used to determine if resolution of grid
C          is large enough
C TIMWGT : IN.  (NPDE)
C          Time weights used in check if time stepsize can be accepted
C RELTOL : IN.  (NPDE)
C          Relative time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C ABSTOL : IN.  (NPDE)
C          Absolute time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C LINSYS : IN.  Linear system solver in use
C               0: BiCGStab + ILU
C               1: GCRO + (Block-)diagonally preconditioning
C INTGRT : IN.  Name of external routine that performs the integration
C          If LINSYS=0: INTGRB, otherwise INTGRC
C RWK    : WORK. (LENRWK)
C LENRWK : IN.  Dimension of RWK.
C          Let NPTS be the max. # points on a grid level and
C          NPTSA the average # points over all grid levels.
C          Then LENRWK should be:
C          MAXLEV=1: 3.NPTS.NPDE+3.NPTS+13.NPTS.NPDE + LSSWRK
C             LSSWRK:
C                ( INFO(4)=0
C                | 38.NPDE.NPTS.NPDE
C                !:INFO(4)=10
C                | (MAX(NPDE.7+3,2.MAXLR+MAXL+6)+NPDE).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=11
C                | (MAX(NPDE.4+3,2.MAXLR+MAXL+6)+NPDE).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=12
C                | (MAX(10,      2.MAXLR+MAXL+6)   +1).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                !:INFO(4)=13
C                | (MAX( 7,      2.MAXLR+MAXL+6)   +1).NPTS.NPDE+
C                  MAXLR*MAXLR+(MAXL+3).MAXL+1
C                )
C             (default: MAXLR = 5, MAXL = 20)
C          Indication of the length for a maximum grid level
C          MAXLEV (default value MAXLEV=3):
C             5.NPTSA.NPDE.MAXLEV+(3+13.NPDE).NPTS + LSSWRK
C IWK    : WORK. (LENIWK)
C LENIWK : IN.  Dimension of IWK.
C          MAXLEV=1: 28.NPTS
C          Indication of the length for a maximum grid level MAXLEV,
C          7.NPTSA.MAXLEV+7.NPTS + ( INFO(4)=0| 19.NPTS )
C LWK    : WORK. (LENLWK)
C LENLWK : IN.  Dimension of LWK >= NPTS+1
C IERR   : INOUT.
C          IN: 0: First call of PDESOL
C              1: Continuation call
C          OUT: 0: OK
C              -1: Workspace too small for required # gridpoints in
C                  base grid. No continuation possible
C              -2: Stepsize too small
C
Ccc EXTERNALS USED:
      LOGICAL CHKWRK, CHKGRD, CHKTIM
      EXTERNAL CHKWRK, CHKGRD, CHKTIM, GETSOL, GETINI, ICOPY, INIGRD,
     +   MKFGRD, MONITR, PDEIV, PUTSOL, RCOPY, SETXYZ
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
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 15)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
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
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
C
C
Ccc   INCLUDE 'CMNWRITEF'
C
C CMNWRITEF
C
C COMMON needed for continuation calls
      INTEGER MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      LOGICAL FIRST, SECOND
      DOUBLE PRECISION T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,
     +   DZB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW,
     +   DXB,DYB,DZB, DTO
      SAVE /WRITIF/, /WRITLF/, /WRITRF/
C
C end INCLUDE 'CMNWRITEF'
C
C
C-----------------------------------------------------------------------
C
CDIR$ NOVECTOR
C
      INTEGER NPTSB, LENUB, LXB, LYB, LZB,
     +   LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP, NIBPTS,
     +   LGNP1, LX, LY, LZ, LUNM1, LUN, LUNP1, LUNP1I, LENU, LUIB,
     +   LGNP1C, LXC, LUNM1C, LUNC, LUNP1C, LENUC,
     +   LGNM1, LGN, LSN, LUNM1T, LENS, LENG, LENGN, LENUN, LEUNP1,
     +   LWT, LF, LCORR, LEV, MLVNM1, MLVN, MLVNP1,
     +   LISTRF, LIWK, LRWK, LLWKN, LIWKN, LRWKN, MAXPTS,
     +   LENPRE, LENLSW, LU, LUO
      LOGICAL LEVN, LEVNM1, OK
      DOUBLE PRECISION DTNEW, DTRAT, TOLWGT, SPCMON, TIMMON
C
      IF (IERR .EQ. 0) THEN
C
Ccc This is the first call of PDESOL.
         T0 = TN
C
C Initialize datastructure, X- and Y-coordinates for base grid
         IF (LINSYS .EQ. 0) THEN
            LRWKN = (3 +16*NPDE + 38*NPDE*NPDE)
         ELSE
            IF (IDIAGP .EQ. 0) THEN
               LENPRE = NPDE
               LENLSW = NPDE*7+3
            ELSE IF (IDIAGP .EQ. 1) THEN
               LENPRE = NPDE
               LENLSW = NPDE*4+3
            ELSE IF (IDIAGP .EQ. 2) THEN
               LENPRE = 1
               LENLSW = 10
            ELSE
               LENPRE = 1
               LENLSW = 7
            ENDIF
            LENLSW = MAX(LENLSW,2*MAXLR+MAXL+6)
            LRWKN = (3 + 16*NPDE + (LENLSW + LENPRE)*NPDE)
         ENDIF
         LIWKN = 28
         MAXPTS = MIN((LENRWK-3)/LRWKN, LENIWK/LIWKN, LENLWK-1)
         CALL INIGRD (MAXPTS, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +      RWK, IWK, NPTSB, LIWKB, IERR)
         DXB = DX
         DYB = DY
         DZB = DZ
         IF (IERR .EQ. 1) THEN
            LRWKN = LRWKN*NPTSB+3
            LIWKN = LIWKN*NPTSB
            LLWKN = NPTSB+1
            OK = CHKWRK (LRWKN+6*NPDE, LENRWK+6*NPDE,
     +              LIWKN+8*MAXLEV+3, LENIWK+8*MAXLEV+3, LLWKN, LENLWK)
            IERR = -1
            RETURN
         ELSE IF (IERR .NE. 0) THEN
            STOP 'Return from INIGRD with unknown IERR'
         ENDIF
         LXB    = 1
         LYB    = LXB + NPTSB
         LZB    = LYB + NPTSB
C
C Set max. grid levels for Tn and Tn-1 at 1 
         LSGNM1(0) = 1
         LSGN  (0) = 1
C
C Set pointers to base grid data structures for Tn-1, Tn and Tn+1
C and to solution for Tn-1 and Tn
         LSGNM1(1) = 1
         LSGN  (1) = 1
         LSGNP1(1) = 1
         LSUNM1(1) = LZB + NPTSB
         LSSN  (1) = LSUNM1(1)
         LSUN  (1) = LSUNM1(1)
C
C Initialize solution values at base grid at Tn = T0
         RWK(LSUN(1)) = 0.0
         CALL PDEIV (T0, RWK(LXB), RWK(LYB), RWK(LZB), RWK(LSUN(1)+1),
     +      NPTSB, NPDE)
         LENUB  = NPTSB*NPDE+1
C
C Set pointer to not updated base grid solution at Tn+1
         LSSNP1(1) = LSUN(1) + LENUB
C
C Initialize time integration variables
         FIRST  = .TRUE.
         SECOND = .FALSE.
C
      ELSE IF (IERR .EQ. 1) THEN
C
Ccc This is a continuation call of PDESOL.
C Set all required variables that were not saved in COMMON
         IF (LINSYS .EQ. 1) THEN
            IF (IDIAGP .EQ. 0) THEN
               LENPRE = NPDE
               LENLSW = NPDE*7+3
            ELSE IF (IDIAGP .EQ. 1) THEN
               LENPRE = NPDE
               LENLSW = NPDE*4+3
            ELSE IF (IDIAGP .EQ. 2) THEN
               LENPRE = 1
               LENLSW = 10
            ELSE
               LENPRE = 1
               LENLSW = 7
            ENDIF
            LENLSW = MAX(LENLSW,2*MAXLR+MAXL+6)
         ENDIF
         NPTSB = IWK(IWK(IWK(1)+2)+2*IWK(1)+2)-1
         LENUB = NPTSB*NPDE+1
         LXB   = 1
         LYB   = LXB + NPTSB
         LZB   = LYB + NPTSB
C
      ELSE
C      This shouldn't happen
         STOP 'PDESOL called with unknown IERR'
      ENDIF
C
Ccccc Time integration loop
   10 CONTINUE
C Adjust time stepsize such that interval TE-TN takes an integer # of
C time steps of this size
      DT = (TE-TN)/INT((TE-TN)/DT+0.95)
      DT = (TN+DT)-TN
C Check if time stepsize is acceptable
      IF (DT .LT. DTMIN) THEN
         WRITE(LUNERR,'(''Time step size too small, DT ='',E16.7)') DT
         IERR = -2
         RETURN
      ENDIF
C
C Time integration method: BE in first time step, BDF2 in following.
C DTRAT = DT / DT_old; 0 => BE
         IF (FIRST) THEN
            DTRAT = 0
         ELSE
            DTRAT = DT / DTO
         ENDIF
C
      LEVEL  = 1
C
Ccc Set pointer to first free element after grid structure of max.
C refinement level for Tn
      IF (MAXLEV .GT. 1) LSGNP1(2) = LIWKB
C
      LGNP1  = LSGNP1(1)
      LX     = LXB
      LY     = LYB
      LZ     = LZB
      DX     = DXB
      DY     = DYB
      DZ     = DZB
      LUNM1  = LSUNM1(1)
      LUN    = LSUN  (1)
      LUNP1  = LSSNP1(1)
      LENU   = LENUB
      LUIB   = LUNP1+LENU
C Pointer to space for eventual refined grid structure
      LISTRF = LIWKB
      LIWK   = LIWKB
      LRWK   = LUNP1 + LENU
C
Ccc Initial solution at coarse grid is coarse grid solution of previous
C time level
      LUNP1I = LSSN(1)
      CALL RCOPY (LENU, RWK(LUNP1I), RWK(LUNP1))
C
Ccccc Grid refinement Loop
C
  100 CONTINUE
      IF (LUNPDS .NE. 0) THEN
         NPLNS = IWK(LGNP1)
         NROWS = IWK(LGNP1+NPLNS+1)-1
         NPTS  = IWK(LGNP1+2*NPLNS+1+NROWS+1)-1
         WRITE(LUNPDS,
     +      '(''Time integration at T='',E10.2,'', Grid level='',I3,
     +        '', NPTS='',I6)') TN+DT, LEVEL, NPTS
      ENDIF
C
Ccc Timestep on current level
      LWT    = LRWK
      LF     = LWT + LENU-1
      LCORR  = LF  + LENU-1
      LRWK   = LCORR + LENU-1
      CALL INTGRT (IWK(LGNP1), RWK(LX), RWK(LY), RWK(LZ), NPDE,
     +   RWK(LUIB), RWK(LUNP1), RWK(LUN), RWK(LUNM1), RTOL, ATOL,
     +   TN, DT, DTRAT, DX, DY, DZ, RWK(LWT), RWK(LF), RWK(LCORR),
     +   RWK(LRWK), IERR)
      LRWK   = LWT
      IF (IERR .EQ. 10) THEN
C      If Newton failure redo time step with stepsize quartered
         NREJS = NREJS+1
         IF (LUNPDS .NE. 0) THEN
            WRITE(LUNPDS,
     +         '(''Newton failure at T='',E10.2,'', Grid level'',I3)')
     +         TN+DT, LEVEL
         ENDIF
         IERR = 0
         DT = DT/4
         GOTO 10
      ELSE IF (IERR .NE. 0) THEN
C      This shouldn't happen
         STOP 'Return from INTGRT with unknown IERR'
      ENDIF
C
Ccc Compute space monitor and if necessary determine new grid
      IF (LSGN(0) .GT. LEVEL) THEN
C    More severe tolerance on grid monitor if max.grid level at Tn 
C    exceeded current level
         TOLWGT = 0.9
      ELSE
         TOLWGT = 1.0
      ENDIF
      OK = CHKGRD (TN+DT, LEVEL, RWK(LUNP1), NPDE,
     +   RWK(LX), RWK(LY), RWK(LZ),
     +   SPCTOL, TOLWGT, IWK(LGNP1), RWK(LRWK), LWK, SPCMON)
C      If no grid refinement needed, check time error
      IF (LUNPDS .NE. 0) THEN
         WRITE(LUNPDS,'(''T='',E10.2,'', LEVEL='',I3,
     +      '' ,TOLWGT='',F3.1,'', SPCMON='',E10.2)')
     +      TN+DT, LEVEL, TOLWGT, SPCMON
      ENDIF
      IF (OK) GOTO 200
      IF (LEVEL .EQ. MAXLEV) THEN
         WRITE(LUNERR,'(''Max. grid level exceeded at T='',E16.7)')
     +      TN+DT
         GOTO 200
      ENDIF
C
Ccc Create refined grid
C Save coarse grid pointers
      LGNP1C = LGNP1
      LUNP1C = LUNP1
      LUNC   = LUN
      LUNM1C = LUNM1
      LXC    = LX
      LENUC  = LENU
C
C Make fine grid structure
      LGNP1  = LISTRF
      CALL MKFGRD (LWK, IWK, LENIWK, LGNP1C, LGNP1, LINSYS,
     +   NPTS, LIWK, IERR)
      LENU   = NPTS*NPDE+1
C
C Check on workspace needed
      IF (LINSYS .EQ. 0) THEN
         LRWKN  = LUNP1C+LENUC+8*LENU+3*NPTS+10*LENU+38*NPDE*LENU
      ELSE
         LRWKN  = LUNP1C+LENUC+8*LENU+3*NPTS+10*LENU+
     +            (LENLSW+LENPRE)*LENU+MAXLR*MAXLR+(MAXL+3)*MAXL
      ENDIF
      LIWKN  = LIWK+NPTS+1
      LLWKN  = NPTS+1
      OK = CHKWRK (LRWKN+6*NPDE, LENRWK+6*NPDE,
     +        LIWKN+8*MAXLEV+3, LENIWK+8*MAXLEV+3, LLWKN, LENLWK)
      IF (.NOT. OK) THEN
         IERR = -1
         RETURN
      ENDIF
C
C Set fine grid pointers and values
      LLPLN  = LGNP1
      NPLNS  = IWK(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      NROWS  = IWK(LLPLN+NPLNS+1)-1
      LLROW  = LIPLN+NPLNS
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = IWK(LLLBND)
      NBDPTS = IWK(LLLBND+NBNDS+1)-1
      NBIPTS = IWK(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LISTRF = LLBNDP+NBIPTS
      LUN    = LUNP1C+LENUC
      LUNP1  = LUN+LENU
      LUNM1  = LUNP1+LENU
      LX     = LUNM1+LENU
      LY     = LX+NPTS
      LZ     = LY+NPTS
      LUIB   = LZ+NPTS
      NIBPTS = NBIPTS-NBDPTS
      LRWK   = LUIB+NIBPTS*NPDE
      LSGNP1(LEVEL+1) = LGNP1
      LSSNP1(LEVEL+1) = LUNP1
      DX = DX/2
      DY = DY/2
      DZ = DZ/2
C
C Save initial solution at current grid level at end of workspace to
C prevent overwriting
      CALL RCOPY (LENUC, RWK(LUNP1I), RWK(LENRWK-LENUC))
      LUNP1I = LENRWK-LENUC
C
C Store grid values at Tn and Tn-1 in temporary storage
      LEVN   = LSGN(0) .GE. LEVEL+1
      LEVNM1 = LSGNM1(0) .GE. LEVEL+1
      IF (FIRST) THEN
C      Store X- and Y- coordinates, and initial solution in Un = Un-1
         LUNM1 = LUN
         LX    = LUNP1+LENU
         LY    = LX+NPTS
         LZ    = LY+NPTS
         LUIB  = LZ+NPTS
         LRWK  = LUIB+NIBPTS*NPDE
         CALL SETXYZ (XL,YF,ZD, DX, DY, DZ,
     +      IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW), IWK(LICOL),
     +      RWK(LX), RWK(LY), RWK(LZ))
         RWK(LUN) = 0.0
         CALL PDEIV (T0, RWK(LX), RWK(LY), RWK(LZ), RWK(LUN+1),
     +      NPTS, NPDE)
C
      ELSE IF (SECOND) THEN
C      Get Un on refined grid
         CALL GETSOL (NPDE, RWK(LUNC), IWK(LGNP1C),
     +      LEVN, RWK(LSUN(LEVEL+1)), IWK(LSGN(LEVEL+1)),
     +      RWK(LUN), IWK(LGNP1), IWK(LIWK), RWK(LRWK))
C      Store X- and Y- coordinates and initial solution in Un-1
         CALL SETXYZ (XL,YF,ZD, DX, DY, DZ,
     +      IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW), IWK(LICOL),
     +      RWK(LX), RWK(LY), RWK(LZ))
         RWK(LUNM1) = 0.0
         CALL PDEIV (T0, RWK(LX), RWK(LY), RWK(LZ), RWK(LUNM1+1),
     +      NPTS, NPDE)
C
      ELSE
C      Get Un-1 and Un on refined grid
         LUNM1T = MAX(LUNM1,LXC)
         CALL GETSOL (NPDE, RWK(LUNM1C), IWK(LGNP1C),
     +      LEVNM1, RWK(LSUNM1(LEVEL+1)), IWK(LSGNM1(LEVEL+1)),
     +      RWK(LUNM1T), IWK(LGNP1), IWK(LIWK), RWK(LRWK))
         IF (LUNM1T .GT. LUNM1)
     +      CALL RCOPY (LENU, RWK(LUNM1T), RWK(LUNM1))
         CALL GETSOL (NPDE, RWK(LUNC), IWK(LGNP1C),
     +      LEVN, RWK(LSUN(LEVEL+1)), IWK(LSGN(LEVEL+1)),
     +      RWK(LUN), IWK(LGNP1), IWK(LIWK), RWK(LRWK))
C      Store X- and Y- coordinates
         CALL SETXYZ (XL,YF,ZD, DX, DY, DZ,
     +      IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW), IWK(LICOL),
     +      RWK(LX), RWK(LY), RWK(LZ))
C
      ENDIF
C Get initial solution Un+1, store internal boundary values also in UIB
C list
      CALL GETINI (NPDE, RWK(LUNP1I), RWK(LUNP1C), IWK(LGNP1C),
     +   LEVN, RWK(LSSN(LEVEL+1)), IWK(LSGN(LEVEL+1)),
     +   RWK(LUNP1), IWK(LGNP1), RWK(LUIB), IWK(LIWK), RWK(LRWK))
      LUNP1I = LRWK
      LRWK   = LUNP1I + LENU
      CALL RCOPY (LENU, RWK(LUNP1), RWK(LUNP1I))
      LEVEL = LEVEL+1
      GOTO 100
Ccc End Refinement Loop
C
  200 CONTINUE
C
Ccc Time step finished
C Inject values from finest level
      LSGNP1(0) = LEVEL
      LSUNP1(LEVEL) = LSSNP1(LEVEL)
      DO 210 LEV = LEVEL, 2, -1
         LSUNP1(LEV-1) = LSUNP1(LEV) + LENU
         CALL PUTSOL (NPDE, RWK(LSUNP1(LEV)), IWK(LSGNP1(LEV)),
     +      RWK(LSSNP1(LEV-1)), IWK(LSGNP1(LEV-1)),
     +      RWK(LSUNP1(LEV-1)), LENU)
  210 CONTINUE
      LRWK = LSUNP1(1) + LENU
C
Ccc Check time-error
      LU  = LSUNP1(1)+LENUB
      LUO = LSSNP1(1)-LENUB
      OK = CHKTIM (RWK, LU, LUO, NPDE, IWK,
     +   LSGNP1, TIMWGT, RELTOL, ABSTOL, RWK(LRWK), DT, DTNEW, TIMMON)
      IF (LUNPDS .NE. 0) THEN
         WRITE(LUNPDS,'(''TN='',E10.2,'', DT='',E10.2,
     +      '', DTNEW='',E10.2, '', TIMMON='',E10.2)')
     +      TN, DT, DTNEW, TIMMON
      ENDIF
C Restrict stepsize
      DTNEW = MIN(DTNEW, DTMAX)
      IF (.NOT. OK) THEN
C
Ccc Time step rejected
         NREJS = NREJS+1
         IF (LUNPDS .NE. 0) THEN
            WRITE(LUNPDS,'(''Time step rejected'')')
         ENDIF
         DT = DTNEW
         GOTO 10
      ELSE
C
Ccc Time step accepted
         NSTEPS = NSTEPS+1
C
Ccc Time step accepted; move data saved for Tn to nm1-save and
C data at Tn+1 to n-save.
C
C Move updated solution at Tn (Un) to Unm1 save, and gridstructure at Tn
C to Gnm1 save
C NB. For first step this is not necessary, but harmless
         MLVNM1 = LSGNM1(0)
         MLVN   = LSGN(0)
C       Start of Unm1 data (= 2*NPTSB+1)
         LUNM1  = LSUNM1(MLVNM1)
C       Start of updated Un data
         LUN    = LSUN(MLVN)
C       LSSNP1(1)-1: end of updated Un data
         LENUN  = LSSNP1(1) - LUN
         CALL RCOPY (LENUN, RWK(LUN), RWK(LUNM1))
C       Adjust pointers to Unm1 data
         DO 220 LEV = MLVN, 1, -1
            LSUNM1(LEV) = LSUN(LEV) - (LUN-LUNM1)
  220    CONTINUE
C       New start of not-updated Un data
         LSN   = LUNM1 + LENUN
C
C       New max. Gnm1-level
         LSGNM1(0) = MLVN
         IF (MLVNM1 .EQ. 1) THEN
C          Grids already in place, adjust pointers
            DO 230 LEV = 2, MLVN
               LSGNM1(LEV) = LSGN(LEV)
  230       CONTINUE
C          New start of Gn data is old one
            LGN    = LIWKB
         ELSE IF (MLVN .GT. 1) THEN
C       Both Gnm1 and Gn have more than 1 level, move Gn
C       Start of Gnm1 data (after base grid)
            LGNM1  = LSGNM1(2)
C       Start of Gn data
            LGN    = LSGN(2)
C       LSGNP1(2)-1: end of Gn data
            LENGN  = LSGNP1(2) - LGN
            CALL ICOPY (LENGN, IWK(LGN), IWK(LGNM1))
C       Adjust pointers to Gnm1 data
            DO 240 LEV = 2, MLVN
               LSGNM1(LEV) = LSGN(LEV) - (LGN-LGNM1)
  240       CONTINUE
C       New start of Gn data
            LGN    = LSGNM1(2) + LENGN
         ELSE
C       At Tn only base grid, new start of Gn data is after base grid
            LGN    = LSGNM1(2)
         ENDIF
C
C Move Tn+1 data, not_updated solution (Snp1) to Sn save, gridstructure
C to Gn save, and injected solution to Un save
         MLVNP1 = LSGNP1(0)
         LSGN(0) = MLVNP1
C   Move not-updated solution Snp1 on base grid
         CALL RCOPY (LENUB, RWK(LSSNP1(1)), RWK(LSN))
         LSSN(1) = LSN
         LSN = LSN + LENUB
C   Move Snp1 and (LROW,IROW,ICOL) of higher levels, adjust pointers to
C   Sn and Gn data
         DO 250 LEV = 2, MLVNP1
            LLPLN  = LSGNP1(LEV)
            NPLNS  = IWK(LLPLN)
            LIPLN  = LLPLN+NPLNS+2
            NROWS  = IWK(LLPLN+NPLNS+1)-1
            LLROW  = LIPLN+NPLNS
            NPTS   = IWK(LLROW+NROWS)-1
            LENS   = NPTS*NPDE+1
            LENG   = NPLNS+2 + NPLNS + NROWS+1 + NROWS + NPTS
            CALL RCOPY (LENS, RWK(LSSNP1(LEV)), RWK(LSN))
            LSSN(LEV) = LSN
            LSN = LSN + LENS
            CALL ICOPY (LENG, IWK(LSGNP1(LEV)), IWK(LGN))
            LSGN(LEV) = LGN
            LGN = LGN + LENG
  250    CONTINUE
C
C   Adjust pointer to solution on highest grid level
         LSUN(MLVNP1) = LSSN(MLVNP1)
         IF (MLVNP1 .GT. 1) THEN
C       Move updated solutions on grids (max.lev-1),...,2 and adjust
C       pointers to Un data
            LUNP1  = LSUNP1(MLVNP1-1)
            LEUNP1 = LSUNP1(1)+LENUB - LUNP1
            CALL RCOPY (LEUNP1, RWK(LUNP1), RWK(LSN))
            DO 260 LEV = 1, MLVNP1-1
               LSUN(LEV) = LSUNP1(LEV) - (LUNP1-LSN)
  260       CONTINUE
         ENDIF
C
Ccc Set pointer to not updated base grid solution at Tn+1
         LSSNP1(1) = LSUN(1) + LENUB
Ccc Set pointer to first free element after grid structure of max.
C refinement level for Tn
         LIWKB = LGN
C
Ccc Adapt time variables
         CALL MONITR (TN+DT, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB,
     +      LSGN, IWK, LSUN, RWK)
         TN = TN + DT
         DTO = DT
         DT  = DTNEW
         IF (FIRST) THEN
            FIRST  = .FALSE.
            SECOND = .TRUE.
         ELSE IF (SECOND) THEN
            SECOND = .FALSE.
         ENDIF
         IF (TN .GE. TE) THEN
            IF (LUNPDS .NE. 0) THEN
               WRITE(LUNPDS,'(''# steps accepted:'',I5,
     +            '', # steps rejected:'',I5)') NSTEPS, NREJS
            ENDIF
            RETURN
         ELSE
            GOTO 10
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE INIGRD (MAXPTS, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   XYZ, IWK, NPTS, LIWK, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER MAXPTS, IWK(*), NPTS, LIWK, IERR
      DOUBLE PRECISION XL, YF, ZD, XR, YB, ZU, DX, DY, DZ, XYZ(*)
C
Ccc PURPOSE:
C Stores datastructure and coordinate values of initial grid (rowwise).
C A (virtual) rectangular box is placed around the irregular domain. The
C intersection point of the left, front, and down plane of this box is
C (XL,YF,ZD) in physical coordinates and (0,0,0) in column, row, resp.
C plane coordinates. The intersection point of the right, back, and
C upper plane is (XR,YB,ZU) resp. (Nx,Ny,Nz), where Nx = (XR-XL)/DX,
C Ny = (YB-YF)/DY, and Nz = (ZU-ZD)/DZ. Only real grid points are stored
C in the order: plane, row, column.
C In the default case the domain is a rectangular prism and the user
C has to specify only the (XL,YF,ZD)- and the (XR,YB,ZU)-point, and
C the gridwidth in each direction. Otherwise the user has to write
C the domain specifying routine INIDOM in which also the coordinate
C values and the cellwidths can be specified.
C
Ccc PARAMETER DESCRIPTION:
C MAXPTS : IN.  Max. # grid points allowed by the available workspace
C XL     : INOUT. X-coordinate of left/front/down point of virtual box
C YF     : INOUT. Y-coordinate of left/front/down point of virtual box
C ZD     : INOUT. Z-coordinate of left/front/down point of virtual box
C XR     : INOUT. X-coordinate of right/back/upper point of virtual box
C YB     : INOUT. Y-coordinate of right/back/upper point of virtual box
C ZU     : INOUT. Z-coordinate of right/back/upper point of virtual box
C DX     : INOUT. Grid width in X-direction
C DY     : INOUT. Grid width in Y-direction
C DZ     : INOUT. Grid width in Z-direction
C XYZ    : OUT. Contains the X-, Y- and Z-coordinates for the base grid
C IWK    : OUT. Contains the following arrays:
CcLPLN   : (0:LPLN(0)+1)
Cc         LPLN(0) = NPLNS: Actual # planes in LROW
Cc         LPLN(1:NPLNS): pointers to the start of a plane in LROW
Cc         LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
CcIPLN   : (NPLNS)
Cc         IPLN(IP): plane number of plane IP in virtual box
CcLROW   : (NROWS+1)
Cc         LROW(1:NROWS): pointers to the start of a row in the grid
Cc         LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
CcIROW   : (NROWS)
Cc         IROW(IR): row number of row IR in virtual box
CcICOL   : (NPTS)
Cc         ICOL(IPT): column number of grid point IPT in virtual box
CcLLBND  : (0:LLBND(0)+2)
Cc         LLBND(0) = NBNDS: total # physical planes in actual domain.
Cc            NB. edges and corners are stored for each plane they
Cc            belong to.
Cc         LLBND(1:NBNDS): pointers to a specific boundary in LBND
Cc         LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
Cc                                    in LBND + 1
Cc         LLBND(NBNDS+1): pointer to internal boundary in LBND
Cc         LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
CcILBND  : (NBNDS)
Cc         ILBND(IB): type of boundary:
Cc                    1: Left  plane -I
Cc                    2: Down  plane  I
Cc                    3: Right plane  I max. first order derivative
Cc                    4: Up    plane  I
Cc                    5: Front plane  I
Cc                    6: Back  plane -I
CcLBND   : (NBIPTS)
Cc         LBND(IBPT): pointer to boundary point in actual grid
CcLBLWY  : IN. (0:NPTS)
Cc         LBLWY(IPT): pointer to node below in Y-direction in
Cc                     actual grid
Cc                      0, if index node is front-plane boundary point
CcLABVY  : IN. (0:NPTS)
Cc         LABVY(IPT): pointer to node above in Y-direction in
Cc                     actual grid
Cc                      0, if index node is back-plane boundary point
CcLBLWZ  : IN. (0:NPTS)
Cc         LBLWZ(IPT): pointer to node below in Z-direction in
Cc                     actual grid
Cc                      0, if index node is down-plane boundary point
CcLABVZ  : IN. (0:NPTS)
Cc         LABVZ(IPT): pointer to node above in Z-direction in
Cc                     actual grid
Cc                      0, if index node is up-plane boundary point
Cc         (Even if LINSYS/=0, because of restart:)
Cc         The next 2 arrays are used for the Jacobian structure and its
Cc         ILU
CcLLDG   : (NPTS,-9:-2)
Cc         LLDG(IPT,-9): pointer to node Y-below Z-below
Cc                       or to node Z-below Z-below
Cc         LLDG(IPT,-8): pointer to node left of Z-below
Cc         LLDG(IPT,-7): pointer to node Z-below
Cc         LLDG(IPT,-6): pointer to node right of Z-below
Cc         LLDG(IPT,-5): pointer to node Y-above Z-below
Cc         LLDG(IPT,-4): pointer to node left of Y-below
Cc                       or to node Y-below Y-below
Cc         LLDG(IPT,-3): pointer to node Y-below
Cc         LLDG(IPT,-2): pointer to node right of Y-below
Cc                       or to node left of the node left
CcLUDG   : (NPTS,2:9)
Cc         LUDG(IPT,2): pointer to node left of Y-above
Cc                      or to node right of the node right
Cc         LUDG(IPT,3): pointer to node Y-above
Cc         LUDG(IPT,4): pointer to node right of node Y-above
Cc                      or to node Y-above Y-above
Cc         LUDG(IPT,5): pointer to node Y-below Z-above
Cc         LUDG(IPT,6): pointer to node left of Z-above
Cc         LUDG(IPT,7): pointer to node Z-above
Cc         LUDG(IPT,8): pointer to node right of Z-above
Cc         LUDG(IPT,9): pointer to node Y-above Z-above
Cc                      or to node Z-above Z-above
Cc         the next 4 arrays are used to hold the data dependency lists
Cc         for the ILU factorization and the forward, resp. backward
Cc         sweep of the backsolve
CcLSL    : LSL(NPTS)
Cc         LSL(ISLPT): pointer to node in actual grid
CcLLSL   : LLSL(0:LLSL(0))
Cc         LLSL(0) = # independent data dependency lists in ILU
Cc                   factorization and forward sweep
Cc         LLSL(1:LLSL(0)): pointers to the start of a list in LSL
CcLSU    : LSU(NPTS)
Cc         LSU(ISLPT): pointer to node in actual grid
CcLLSU   : LLSU(0:LLSU(0))
Cc         LLSU(0) = # independent data dependency lists in backward
C                    sweep
Cc         LLSU(1:LLSU(0)): pointers to the start of a list in LSU
C NPTS   : OUT. # grid points in base grid
C LIWK   : OUT. Pointer to first free element in IWK
C IERR   : OUT. Error return flag
C          0: OK.
C          1: workspace too small for required # gridpoints
C
Ccc EXTERNALS USED:
      LOGICAL INIDOM
      EXTERNAL ICOPY, INIDOM, JACSDP, JACSLP, JACSUP, SETBA, SETXYZ
C
C-----------------------------------------------------------------------
C
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8,
     +   LLPLN, LIPLN, LLROW, LIROW, LICOL,
     +   LLLBND, LILBND, LLBNDP, LLBLWY, LLABVY, LLBLWZ, LLABVZ,
     +   NPLNS, NROWS, NBNDS, NBDPTS,
     +   LLLDG, LLUDG, LLSLP, LLLSL, LLSUP, LLLSU
      LOGICAL OK
C
      IERR = 0
C
Ccc Get (user defined) initial domain
      I1 = 1
      I2 = I1 + MAXPTS
      I3 = I2 + MAXPTS
      I4 = I3 + MAXPTS
      I5 = I4 + MAXPTS
      I6 = I5 + MAXPTS
      I7 = I6 + MAXPTS
      I8 = I7 + MAXPTS
      OK = INIDOM (MAXPTS, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   IWK(I1), IWK(I2), IWK(I3), IWK(I4), IWK(I5), IWK(I6), IWK(I7),
     +   IWK(I8))
      IF (.NOT. OK) THEN
         IERR = 1
         NPTS = MAXPTS
         RETURN
      ENDIF
C
Ccc Copy integer arrays to their correct position in the IWK array
      NPLNS  = IWK(I1)
      NROWS  = IWK(I1+NPLNS+1)-1
      NPTS   = IWK(I3+NROWS)-1
      NBNDS  = IWK(I6)
      NBDPTS = IWK(I6+NBNDS+1)-1
C LPLN at correct position
      LLPLN = 1
C Copy IPLN
      LIPLN = LLPLN+NPLNS+2
      CALL ICOPY (NPLNS, IWK(I2), IWK(LIPLN))
C Copy LROW
      LLROW = LIPLN+NPLNS
      CALL ICOPY (NROWS+1, IWK(I3), IWK(LLROW))
C Copy IROW
      LIROW = LLROW+NROWS+1
      CALL ICOPY (NROWS, IWK(I4), IWK(LIROW))
C Copy ICOL
      LICOL = LIROW+NROWS
      CALL ICOPY (NPTS, IWK(I5), IWK(LICOL))
C Copy LLBND
      LLLBND = LICOL+NPTS
      CALL ICOPY (NBNDS+2, IWK(I6), IWK(LLLBND))
C No internal boundaries
      IWK(LLLBND+NBNDS+2) = NBDPTS+1
C Copy ILBND
      LILBND = LLLBND+NBNDS+3
      CALL ICOPY (NBNDS, IWK(I7), IWK(LILBND))
C Copy LBND
      LLBNDP = LILBND+NBNDS
      CALL ICOPY (NBDPTS, IWK(I8), IWK(LLBNDP))
C
Ccc Store X-, Y-, and Z-coordinates
      CALL SETXYZ (XL, YF, ZD, DX, DY, DZ,
     +   IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW), IWK(LICOL),
     +   XYZ(1), XYZ(1+NPTS), XYZ(1+2*NPTS))
C
Ccc Set pointers to nodes below and above a grid point
      LLBLWY = LLBNDP+NBDPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
      CALL SETBA (IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW),
     +   IWK(LICOL), IWK(LLBLWY), IWK(LLABVY), IWK(LLBLWZ), IWK(LLABVZ))
      LIWK = LLABVZ+NPTS+1
C
Ccc Set pointers to lower and upper diagonals in Jacobian for base grid
      LLLDG  = LIWK
      LLUDG  = LLLDG + NPTS*8 
      CALL JACSDP (NPTS, IWK(LLLBND), IWK(LILBND), IWK(LLBNDP),
     +   IWK(LLBLWY), IWK(LLABVY), IWK(LLBLWZ), IWK(LLABVZ),
     +   IWK(LLLDG), IWK(LLUDG))
      LIWK = LLUDG + NPTS*8 
C
Ccc Make data-dependency lists for ILU on base-grid Jacobian
      LLSLP  = LIWK
      LLLSL  = LLSLP  + NPTS
      LIWK   = LLLSL  + NPTS
      CALL JACSLP (NPTS, IWK(LLLBND), IWK(LILBND), IWK(LLBNDP),
     +   IWK(LLLDG), IWK(LIWK), IWK(LLLSL), IWK(LLSLP))
      LLSUP  = LLLSL + IWK(LLLSL)+1
      LLLSU  = LLSUP + NPTS
      LIWK   = LLLSU + NPTS
      CALL JACSUP (NPTS, IWK(LLLBND), IWK(LILBND), IWK(LLBNDP),
     +   IWK(LLUDG), IWK(LIWK), IWK(LLLSU), IWK(LLSUP))
      LIWK   = LLLSU + IWK(LLLSU)+1

      RETURN
      END
      SUBROUTINE SETBA (LPLN, IPLN, LROW, IROW, ICOL,
     +   LBLWY, LABVY, LBLWZ, LABVZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
     +   LBLWY(0:*), LABVY(0:*), LBLWZ(0:*), LABVZ(0:*)
C
Ccc PURPOSE:
C Set pointers to nodes below and above a grid point, if such a node
C exists, otherwise the pointer is set to zero.
C
Ccc PARAMETER DESCRIPTION:
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LBLWY  : OUT. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : OUT. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : OUT. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : OUT. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IP, IPT, IPTA, IR, IRA, NPLNS, NROWS, NPTS
C
      NPLNS = LPLN(0)
      NROWS = LPLN(NPLNS+1)-1
      NPTS  = LROW(NROWS+1)-1
      DO 10 IPT = 0, NPTS
         LBLWY(IPT) = 0
         LABVY(IPT) = 0
         LBLWZ(IPT) = 0
         LABVZ(IPT) = 0
   10 CONTINUE
      DO 20 IP = 1, NPLNS
         DO 30 IR = LPLN(IP), LPLN(IP+1)-2
C Check if next row in datastructure is next row in virtual plane
            IF (IROW(IR) .EQ. IROW(IR+1)-1) THEN
C            Compare column indices in row with those in row above,
C            until the two match or one of both rows is finished
               IPTA = LROW(IR+1)
               DO 40 IPT = LROW(IR), LROW(IR+1)-1
   50             IF (ICOL(IPT) .LT. ICOL(IPTA)) THEN
                     GOTO 40
                  ELSE IF (ICOL(IPT) .EQ. ICOL(IPTA)) THEN
C               Set above pointer in current row and below pointer in
C               next row
                     LBLWY(IPTA) = IPT
                     LABVY(IPT)  = IPTA
                  ELSE
                     IPTA = IPTA + 1
                     IF (IPTA .GT. LROW(IR+2)-1) GOTO 30
                     GOTO 50
                  ENDIF
   40          CONTINUE
            ENDIF
   30    CONTINUE
   20 CONTINUE
      DO 100 IP = 1, NPLNS-1
C Check if next plane in datastructure is next plane in virtual box
         IF (IPLN(IP) .EQ. IPLN(IP+1)-1) THEN
C         Compare row indices in plane with those in plane above, until
C         the two match or one of both planes is finished
            IRA = LPLN(IP+1)
            DO 110 IR = LPLN(IP), LPLN(IP+1)-1
  120          IF (IROW(IR) .LT. IROW(IRA)) THEN
                  GOTO 110
               ELSE IF (IROW(IR) .EQ. IROW(IRA)) THEN
C               Compare column indices in row with those in row above,
C               until the two match or one of both rows is finished
                  IPTA = LROW(IRA)
                  DO 130 IPT = LROW(IR), LROW(IR+1)-1
  140                IF (ICOL(IPT) .LT. ICOL(IPTA)) THEN
                        GOTO 130
                     ELSE IF (ICOL(IPT) .EQ. ICOL(IPTA)) THEN
C                     Set above pointer in current row and below pointer
C                     in next row
                        LBLWZ(IPTA) = IPT
                        LABVZ(IPT)  = IPTA
                     ELSE
                        IPTA = IPTA + 1
                        IF (IPTA .GT. LROW(IRA+1)-1) GOTO 110
                        GOTO 140
                     ENDIF
  130             CONTINUE
               ELSE
                  IRA = IRA + 1
                  IF (IRA .GT. LPLN(IP+2)-1) GOTO 100
                  GOTO 120
               ENDIF
  110       CONTINUE
         ENDIF
  100 CONTINUE

      RETURN
      END
      LOGICAL FUNCTION CHKGRD (T, LEVEL, U, NPDE, X, Y, Z, SPCTOL,
     +   TOLWGT, ISTRUC, WORK, REFFLG, SPCMON)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LEVEL, NPDE, ISTRUC(0:*)
      LOGICAL REFFLG(0:*)
      DOUBLE PRECISION T, U(0:*), X(*), Y(*), Z(*), 
     +   SPCTOL(NPDE), TOLWGT, WORK(*), SPCMON
C
Ccc PURPOSE:
C Check if grid needs to be refined. If so, CHKGRD = .FALSE. and
C flags are set where the refinement is needed.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Current time level
C LEVEL  : IN. Current grid level
C U      : IN. Solution on current grid
C NPDE   : IN. # PDE components
C X,Y,Z  : IN. Physical coordinates of grid points
C SPCTOL : IN. User defined space tolerance for each PDE component
C TOLWGT : IN. Weight factor for tolerance. If new level at previous
C          time existed TOLWGT < 1, else 1
C ISTRUC : IN. Datastructure for current grid
C WORK   : WORK. (3*NPTS*NPDE)
C REFFLG : OUT. If one of the corners of a cell is flagged the cell
C          needs to be refined
C SPCMON : OUT. Value of space monitor
C
Ccc EXTERNALS USED:
      LOGICAL CHKREF
      EXTERNAL CHKREF
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ
C
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
C
Ccc Compute space monitor and check if grid needs to be refined
      CHKGRD = CHKREF (T, LEVEL, U, NPTS, NPDE, X, Y, Z,
     +   ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP), 
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY), ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +   SPCTOL, TOLWGT, WORK(1),WORK(1+NPTS*NPDE),WORK(1+2*NPTS*NPDE),
     +   REFFLG, SPCMON)

      RETURN
      END
      SUBROUTINE MKFGRD (REFFLG, IWK, LENIWK, LISTRC, LISTRF, LINSYS,
     +   NPTSF, LIWK, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LENIWK
      INTEGER IWK(LENIWK), LISTRC, LISTRF, LINSYS, NPTSF, LIWK, IERR
      LOGICAL REFFLG(0:*)
C
Ccc PURPOSE:
C Make fine grid structure and store in IWK(LISTRF+...)
C
Ccc PARAMETER DESCRIPTION:
C REFFLG : IN. If one of the corners of a cell is flagged the cell
C          needs to be refined
C IWK    : INOUTWORK. Integer workspace.
C          IN: Coarse grid structure
C          OUT: If coarse grid is base grid: coarse grid structure,
C               otherwise (LPLN,...,LBND) of coarse grid structure, and
C               fine grid structure (LPLN,...,LLLSU)
C          WORK: (NPTS+1+NPTSF+1) at end of IWK for domain flags
C LENIWK : IN. Length of IWK
C LISTRC : IN. Pointer to coarse grid structure in IWK
C LISTRF : IN. Pointer to place where fine grid structure should be
C          stored in IWK
C LINSYS : IN.  Linear system solver in use
C               0: BiCGStab + ILU
C               1: Diagonally scaled GCRO
C NPTSF  : OUT. # grid points in fine grid
C LIWK   : OUT. Pointer to first free element in IWK after fine grid
C          structure
C IERR   : OUT. Error return flag
C          0: OK
C          1: workspace too small for required # fine gridpoints
C
Ccc EXTERNALS USED:
      EXTERNAL DOMFLG, ICOPY, JACSDP, JACSLP, JACSUP, MKBND, REFDOM,
     +   SETBA
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ, LLLDG,
     +   LLPLNF, NPLNSF, LIPLNF, LLROWF, NROWSF, LIROWF, LICOLF,
     +   LLLBDF, NBIPTF, LILBDF, LLBNDF,
     +   LLBLYF, LLABYF, LLBLZF, LLABZF,
     +   LLLDGF, LLUDGF, LLSLF, LLLSLF, LLSUF, LLLSUF,
     +   LIDOM, LIDOMF, MAXPTS
C
      IERR = 0
C
      LLPLN  = LISTRC
      NPLNS  = IWK(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NROWS  = IWK(LLPLN+NPLNS+1)-1
      NPTS   = IWK(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = IWK(LLLBND)
      NBIPTS = IWK(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
      LLLDG  = LLABVZ+NPTS+1
      LIDOM  = LENIWK-NPTS
C
Ccc Make data structure fine grid
      MAXPTS = (LIDOM-MAX(LISTRF,LLLDG))/6
      LIDOMF = LIDOM -MAXPTS-1
      LICOLF = LIDOMF-MAXPTS
      LIROWF = LICOLF-MAXPTS+1
      LLROWF = LIROWF-MAXPTS
      LIPLNF = LLROWF-MAXPTS
      LLPLNF = LIPLNF-MAXPTS
C Set domain flags for coarse grid
      CALL DOMFLG (NPTS, IWK(LLLBND), IWK(LLBNDP), IWK(LIDOM))
      CALL REFDOM (MAXPTS, REFFLG, NBNDS,
     +   IWK(LLPLN), IWK(LIPLN), IWK(LLROW), IWK(LIROW), IWK(LICOL),
     +   IWK(LLBLWY), IWK(LLABVY), IWK(LLBLWZ), IWK(LLABVZ), IWK(LIDOM),
     +   IWK(LLPLNF), IWK(LIPLNF), IWK(LLROWF), IWK(LIROWF),
     +   IWK(LICOLF), IWK(LIDOMF), NPTSF, IERR)
      IF (IERR .EQ. 1) THEN
         LIWK = LENIWK+1
         RETURN
      ENDIF
C
Ccc Move fine grid structure to their correct position
      NPLNSF = IWK(LLPLNF)
      CALL ICOPY (NPLNSF+2, IWK(LLPLNF), IWK(LISTRF))
      LLPLNF = LISTRF
      CALL ICOPY (NPLNSF, IWK(LIPLNF), IWK(LLPLNF+NPLNSF+2))
      LIPLNF = LLPLNF+NPLNSF+2
      NROWSF = IWK(LLPLNF+NPLNSF+1)-1
      CALL ICOPY (NROWSF+1, IWK(LLROWF), IWK(LIPLNF+NPLNSF))
      LLROWF = LIPLNF+NPLNSF
      CALL ICOPY (NROWSF, IWK(LIROWF), IWK(LLROWF+NROWSF+1))
      LIROWF = LLROWF+NROWSF+1
      CALL ICOPY (NPTSF, IWK(LICOLF), IWK(LIROWF+NROWSF))
      LICOLF = LIROWF+NROWSF
C
Ccc Copy # physical boundaries and boundary types from coarse grid
      LLLBDF = LICOLF+NPTSF
      LILBDF = LLLBDF+NBNDS+3
      LLBNDF = LILBDF+NBNDS
      IWK(LLLBDF) = NBNDS
      CALL ICOPY (NBNDS, IWK(LILBND), IWK(LILBDF))
C
Ccc Set pointers below and above and new boundary lists
      LLBLYF = LLBNDF+NPTSF
      LLABYF = LLBLYF+NPTSF+1
      LLBLZF = LLABYF+NPTSF+1
      LLABZF = LLBLZF+NPTSF+1
      LIWK   = LLABZF+NPTSF+1
      IF (LIWK .GT. LENIWK) THEN
         IERR = 1
         RETURN
      ENDIF
      CALL SETBA (IWK(LLPLNF), IWK(LIPLNF), IWK(LLROWF), IWK(LIROWF),
     +   IWK(LICOLF), IWK(LLBLYF),IWK(LLABYF), IWK(LLBLZF),IWK(LLABZF))
      CALL MKBND (NPTSF,
     +   IWK(LLPLNF),IWK(LIPLNF), IWK(LLROWF),IWK(LIROWF), IWK(LICOLF),
     +   IWK(LLBLYF), IWK(LLABYF), IWK(LLBLZF),IWK(LLABZF),
     +   IWK(LIDOMF), IWK(LLLBDF), IWK(LILBDF), IWK(LLBNDF))
C
Ccc Move below/above pointers to their correct position
      NBIPTF = IWK(LLLBDF+NBNDS+2)-1
      CALL ICOPY (NPTSF+1, IWK(LLBLYF), IWK(LLBNDF+NBIPTF))
      LLBLYF = LLBNDF+NBIPTF
      CALL ICOPY (NPTSF+1, IWK(LLABYF), IWK(LLBLYF+NPTSF+1))
      LLABYF = LLBLYF+NPTSF+1
      CALL ICOPY (NPTSF+1, IWK(LLBLZF), IWK(LLABYF+NPTSF+1))
      LLBLZF = LLABYF+NPTSF+1
      CALL ICOPY (NPTSF+1, IWK(LLABZF), IWK(LLBLZF+NPTSF+1))
      LLABZF = LLBLZF+NPTSF+1
      LIWK   = LLABZF+NPTSF+1
      IF (LINSYS .NE. 0) RETURN
C
Ccc Set pointers to lower and upper diagonals in Jacobian for fine grid
      LLLDGF = LIWK
      LLUDGF = LLLDGF+NPTSF*8
      LIWK   = LLUDGF+NPTSF*8
      IF (LIWK .GT. LENIWK) THEN
         IERR = 1
         RETURN
      ENDIF
      CALL JACSDP (NPTSF, IWK(LLLBDF), IWK(LILBDF), IWK(LLBNDF),
     +   IWK(LLBLYF), IWK(LLABYF), IWK(LLBLZF), IWK(LLABZF),
     +   IWK(LLLDGF), IWK(LLUDGF))
C
Ccc Make data-dependency lists for ILU on fine-grid Jacobian
      LLSLF  = LIWK
      LLLSLF = LLSLF +NPTSF
      LIWK   = LLLSLF+NPTSF
      IF (LIWK+NPTSF .GT. LENIWK) THEN
         IERR = 1
         RETURN
      ENDIF
      CALL JACSLP (NPTSF, IWK(LLLBDF), IWK(LILBDF), IWK(LLBNDF),
     +   IWK(LLLDGF), IWK(LIWK), IWK(LLLSLF), IWK(LLSLF))
      LLSUF  = LLLSLF+IWK(LLLSLF)+1
      LLLSUF = LLSUF +NPTSF
      LIWK   = LLLSUF+NPTSF
      IF (LIWK+NPTSF .GT. LENIWK) THEN
         LIWK = LIWK+NPTSF
         IERR = 1
         RETURN
      ENDIF
      CALL JACSUP (NPTSF, IWK(LLLBDF), IWK(LILBDF), IWK(LLBNDF),
     +   IWK(LLUDGF), IWK(LIWK), IWK(LLLSUF), IWK(LLSUF))
      LIWK   = LLLSUF+ IWK(LLLSUF)+1
 
      RETURN
      END
      LOGICAL FUNCTION CHKREF (T, LEVEL, U, NPTS, NPDE, X, Y, Z,
     +   LLBND, ILBND, LBND, LBLWY, LABVY, LBLWZ, LABVZ,
     +   SPCTOL, TOLWGT, W1, W2, W3, REFFLG, SPCMON)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LEVEL, NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*),
     +   LBLWY(0:NPTS), LABVY(0:NPTS), LBLWZ(0:NPTS), LABVZ(0:NPTS)
      LOGICAL REFFLG(0:NPTS)
      DOUBLE PRECISION T, U(0:NPTS*NPDE), X(NPTS), Y(NPTS), Z(NPTS),
     +   SPCTOL(NPDE), TOLWGT,
     +   W1(NPTS*NPDE), W2(NPTS*NPDE), W3(NPTS*NPDE), SPCMON
C
Ccc PURPOSE:
C Check if grid needs to be refined. If so, CHKREF = .FALSE. and
C flags are set where the refinement is needed.
C
C Space monitor:
C SPCMON(ipt) = max SPCTOL(ic).(|(dx)^2.Uxx(ipt)| + |(dy)^2.Uyy(ipt)|
C            (ic = 1, NPDE)                       + |(dz)^2.Uzz(ipt)|)
C A user routine is called to eventually enforce refinement by setting
C SPCMON.
C If max SPCMON(ipt) < TOLWGT then no refinement is needed,
C    (ipt = 1, NPTS)
C otherwise all gridpoints for which SPCMON(ipt) > 1/4 are flagged
C plus their 26 neighbours.
C On exit CHKREF = .TRUE. if no refinement is required
C
Ccc PARAMETER DESCRIPTION:
C T      : IN. Current time level
C LEVEL  : IN. Current grid level
C U      : IN. Array of solution values.
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C X,Y,Z  : IN. Physical coordinates of grid points
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LBLWY  : OUT. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LABVY  : OUT. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LBLWZ  : OUT. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : OUT. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C SPCTOL : IN. User defined space tolerance for the monitor values of
C          different components
C TOLWGT : IN. Weight factor for tolerance. If new level at previous
C          time existed TOLWGT < 1, else 1
C W1,W2,W3 : WORK.
C REFFLG : OUT. If the solution in a grid point violates the space
C          monitor condition, the gridpoint and its 26 neighbours are
C          flagged
C SPCMON : OUT. Max SPCMON(ipt)
C
Ccc EXTERNALS USED:
      EXTERNAL CHSPCM
C
C-----------------------------------------------------------------------
C
      INTEGER I, IB, IC, IPT, IM1, IM2, IP1, IP2, LB, NBNDS, NBIPTS
C
      NBNDS  = LLBND(0)
      NBIPTS = LLBND(NBNDS+2)-1
C
Ccc Store (dx)^2.Uxx in W1, (dy)^2.Uyy in W2, and (dz)^2.Uzz in W3
C First interior points, boundary values will be rubbish
      DO 10 I = 2, NPTS*NPDE-1
         W1(I) = U(I-1) - 2*U(I) + U(I+1)
   10 CONTINUE
      DO 20 IC = 1, NPDE
      DO 20 IPT = 1, NPTS
         I = IPT + (IC-1)*NPTS
         IM1 = LBLWY(IPT) + (IC-1)*NPTS
         IP1 = LABVY(IPT) + (IC-1)*NPTS
         W2(I) = U(IM1) - 2*U(I) + U(IP1)
         IM1 = LBLWZ(IPT) + (IC-1)*NPTS
         IP1 = LABVZ(IPT) + (IC-1)*NPTS
         W3(I) = U(IM1) - 2*U(I) + U(IP1)
   20 CONTINUE
C
C Correct boundaries, first the physical boundaries then the internal
C ones
      DO 30 IB = 1, NBNDS
         IF (ILBND(IB) .EQ. 1) THEN
C       Left boundary plane, correct (dx)^2.Uxx in W1
            DO 40 IC = 1, NPDE
            DO 40 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               W1(I) = U(I) - 2*U(I+1) + U(I+2)
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down boundary plane, correct (dz)^2.Uzz in W3
            DO 50 IC = 1, NPDE
            DO 50 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               IP1 = LABVZ(IPT)
               IP2 = LABVZ(IP1) + (IC-1)*NPTS
               IP1 = IP1 + (IC-1)*NPTS
               W3(I) = U(I) - 2*U(IP1) + U(IP2)
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right boundary plane, correct (dx)^2.Uxx in W1
            DO 60 IC = 1, NPDE
            DO 60 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               W1(I) = U(I) - 2*U(I-1) + U(I-2)
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up boundary plane, correct (dz)^2.Uzz in W3
            DO 70 IC = 1, NPDE
            DO 70 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               IM1 = LBLWZ(IPT)
               IM2 = LBLWZ(IM1) + (IC-1)*NPTS
               IM1 = IM1 + (IC-1)*NPTS
               W3(I) = U(I) - 2*U(IM1) + U(IM2)
   70       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front boundary plane, correct (dy)^2.Uyy in W2
            DO 80 IC = 1, NPDE
            DO 80 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               IP1 = LABVY(IPT)
               IP2 = LABVY(IP1) + (IC-1)*NPTS
               IP1 = IP1 + (IC-1)*NPTS
               W2(I) = U(I) - 2*U(IP1) + U(IP2)
   80       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back boundary plane, correct (dy)^2.Uyy in W2
            DO 90 IC = 1, NPDE
            DO 90 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               IM1 = LBLWY(IPT)
               IM2 = LBLWY(IM1) + (IC-1)*NPTS
               IM1 = IM1 + (IC-1)*NPTS
               W2(I) = U(I) - 2*U(IM1) + U(IM2)
   90       CONTINUE
         ENDIF
   30 CONTINUE
      IB = NBNDS + 1
C       Internal boundary, Dirichlet condition, space error = 0
            DO 210 IC = 1, NPDE
            DO 210 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I = IPT + (IC-1)*NPTS
               W1(I) = 0.0
               W2(I) = 0.0
               W3(I) = 0.0
  210       CONTINUE
C
Ccc Compute space monitor SPCMON, and its maximum
      IC = 1
      DO 300 IPT = 1, NPTS
         W1(IPT) = SPCTOL(IC)*(ABS(W1(IPT)) + ABS(W2(IPT)) +
     +                         ABS(W3(IPT)))
  300 CONTINUE
      DO 310 IC = 2, NPDE
      DO 310 IPT = 1, NPTS
         I = IPT + (IC-1)*NPTS
         W1(IPT) = MAX(W1(IPT),SPCTOL(IC)*(ABS(W1(I)) + ABS(W2(I)) +
     +                                     ABS(W3(I))))
  310 CONTINUE
C
C Call user routine to possibly force refinement
      CALL CHSPCM (T, LEVEL, NPTS, X, Y, Z, NPDE, U(1), W1, TOLWGT)
C
C Compute maximum
      SPCMON = W1(1)
      DO 320 IPT = 2, NPTS
         SPCMON = MAX(SPCMON,W1(IPT))
  320 CONTINUE
C
Ccc Check if grid refinement is needed
      IF (SPCMON .LT. TOLWGT) THEN
C      No refinement needed
         CHKREF = .TRUE.
         RETURN
      ENDIF
C
Ccc Flag each node where space monitor is too large + its 26 neighbors.
C Cells will be refined if a flag is set on 1 corner
      CHKREF = .FALSE.
      DO 400 IPT = 0, NPTS
         REFFLG(IPT) = .FALSE.
  400 CONTINUE
C
C If neighbors in the grid datastructure are not physical neighbors in
C the grid the former are wrongly flagged but since those points
C are boundary points the flags will be unset later on
      DO 410 IPT = 1, NPTS
         IF (W1(IPT) .GE. 0.25) THEN
            IM1 = LBLWZ(IPT)
            IF (IM1 .GT. 0) THEN
               REFFLG(LABVY(IM1-1)) = .TRUE.
               REFFLG(      IM1-1)  = .TRUE.
               REFFLG(LBLWY(IM1-1)) = .TRUE.
            ENDIF
            REFFLG(LABVY(IM1))   = .TRUE.
            REFFLG(      IM1)    = .TRUE.
            REFFLG(LBLWY(IM1))   = .TRUE.
            REFFLG(LABVY(IM1+1)) = .TRUE.
            REFFLG(      IM1+1)  = .TRUE.
            REFFLG(LBLWY(IM1+1)) = .TRUE.
C
            REFFLG(LABVY(IPT-1)) = .TRUE.
            REFFLG(      IPT-1)  = .TRUE.
            REFFLG(LBLWY(IPT-1)) = .TRUE.
            REFFLG(LABVY(IPT))   = .TRUE.
            REFFLG(      IPT)    = .TRUE.
            REFFLG(LBLWY(IPT))   = .TRUE.
            IF (IPT .LT. NPTS) THEN
               REFFLG(LABVY(IPT+1)) = .TRUE.
               REFFLG(      IPT+1)  = .TRUE.
               REFFLG(LBLWY(IPT+1)) = .TRUE.
            ENDIF
C
            IP1 = LABVZ(IPT)
            IF (IP1 .GT. 0) THEN
               REFFLG(LABVY(IP1-1)) = .TRUE.
               REFFLG(      IP1-1)  = .TRUE.
               REFFLG(LBLWY(IP1-1)) = .TRUE.
            ENDIF
            REFFLG(LABVY(IP1))   = .TRUE.
            REFFLG(      IP1)    = .TRUE.
            REFFLG(LBLWY(IP1))   = .TRUE.
            IF (IP1 .LT. NPTS) THEN
               REFFLG(LABVY(IP1+1)) = .TRUE.
               REFFLG(      IP1+1)  = .TRUE.
               REFFLG(LBLWY(IP1+1)) = .TRUE.
            ENDIF
         ENDIF
  410 CONTINUE
      REFFLG(0) = .FALSE.
C Unset errorflags at (physical and internal) boundary
      DO 430 LB = 1, NBIPTS
         IPT = LBND(LB)
         REFFLG(IPT) = .FALSE.
  430 CONTINUE

      RETURN
      END
      SUBROUTINE DOMFLG (NPTS, LLBND, LBND, IDOM)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, LLBND(0:*), LBND(*), IDOM(0:NPTS)
C
Ccc PURPOSE:
C Set domain flags for determination of location of grid point in grid
C 
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary or corner in
C                          LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C IDOM   : OUT. IDOM(IPT): location in domain of node IPT
C          0: interior point
C          IB: on boundary
C              if not on edge or corner, or if internal boundary
C              then IB, IB = 1, NBNDS+1
C              else | IB4 | IB3 | IB2 | IBYTE | if point is part of
C                   the boundary planes IBYTE, IB2, IB3, and IB4
C
Ccc EXTERNALS USED: NONE
C
      INTEGER BYTE
      PARAMETER (BYTE = 2**8)
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, IB, LB, NBNDS
C
      NBNDS = LLBND(0)
C
C Set domain flags
      DO 10 IPT = 0, NPTS
         IDOM(IPT) = 0
   10 CONTINUE
      DO 20 IB = 1, NBNDS
         DO 30 LB = LLBND(IB), LLBND(IB+1)-1
            IPT = LBND(LB)
            IF (IDOM(IPT) .EQ. 0) THEN
               IDOM(IPT) = IB
            ELSE
               IDOM(IPT) = IB+BYTE*IDOM(IPT)
            ENDIF
   30    CONTINUE
   20 CONTINUE
      IB = NBNDS+1
         DO 40 LB = LLBND(IB), LLBND(IB+1)-1
            IPT = LBND(LB)
            IDOM(IPT) = IB
   40    CONTINUE

      RETURN
      END
      SUBROUTINE REFDOM (MAXPTS, REFFLG, NBNDS,
     +   LPLNC, IPLNC, LROWC, IROWC, ICOLC,
     +   LBLWYC, LABVYC, LBLWZC, LABVZC, IDOMC,
     +   LPLN, IPLN, LROW,  IROW,  ICOL, IDOM, NPTS, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER MAXPTS, NBNDS,
     +   LPLNC(0:*), IPLNC(*), LROWC(*), IROWC(*), ICOLC(*),
     +   LABVYC(0:*), LBLWYC(0:*), LABVZC(0:*), LBLWZC(0:*), IDOMC(0:*),
     +   LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*), IDOM(0:*),
     +   NPTS, IERR
      LOGICAL REFFLG(0:*)
C
Ccc PURPOSE:
C Create refined grid. If one of the corners of a cell is flagged, the
C cell is divided in 8. The (LPLN,IPLN,LROW,IROW,ICOL) structure of the
C fine grid will be stored and IDOM will contain domainflags (only for
C rows corresponding with the coarse grid) to indicate the location of a
C node in the fine grid.
C
Ccc PARAMETER DESCRIPTION:
C MAXPTS : IN. Max. # grid points allowed on fine grid
C REFFLG : IN. (0:NPTSC)
C          If the solution in a grid point violates the space monitor
C          condition, the gridpoint and its 26 neighbors are flagged.
C          Gridpoints at the boundaries are not flagged
C LPLNC  : IN. (0:LPLNC(0)+1)
C          LPLNC(0) = NPLNSC: Actual # planes in LROWC
C          LPLNC(1:NPLNSC): pointers to the start of a plane in LROWC
C          LPLNC(NPLNSC+1) = NROWSC+1: Total # rows in coarse grid + 1
C IPLNC  : IN. (NPLNSC)
C          IPLNC(IP): plane number of plane IP in virtual box
C LROWC  : IN. (NROWSC+1)
C          LROWC(1:NROWSC): pointers to the start of a row in the grid
C          LROWC(NROWSC+1) = NPTSC+1: Actual # nodes in grid + 1
C IROWC  : IN. (NROWSC)
C          IROWC(IR): row number of row IR in virtual box
C ICOLC  : IN. (NPTSC)
C          ICOLC(IPT): column number of grid point IPT in virtual box
C LBLWYC : IN. (0:NPTSC)
C          LBLWYC(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVYC : IN. (0:NPTSC)
C          LABVYC(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZC : IN. (0:NPTSC)
C          LBLWZC(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZC : IN. (0:NPTSC)
C          LABVZC(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C IDOMC  : IN. (0:NPTSC)
C          IDOMC(IPT): location in coarse grid of node IPT
C          0: interior point
C          IB: on boundary
C              if not on edge or corner, or if internal boundary
C              then IB, IB = 1, NBNDS+1
C              else | IB4 | IB3 | IB2 | IBYTE | if point is part of
C                   the boundary planes IBYTE, IB2, IB3, and IB4
C LPLN   : OUT. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : OUT. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : OUT. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : OUT. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : OUT. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C IDOM   : OUT. (0:NPTS)
C          IDOM(IPT): location in coarse grid of node IPT (only set for
C          rows corresponding with coarse grid rows)
C          0: interior point or new horizontal internal boundary
C          IB: on boundary
C              if not on edge or corner, or if internal boundary
C              then IB, IB = 1, NBNDS+1
C              else | IB4 | IB3 | IB2 | IBYTE | if point is part of
C                   the boundary planes IBYTE, IB2, IB3, and IB4
C NPTS   : OUT. # grid points on fine grid or MAXPTS if IERR=1
C IERR   : OUT. Error return flag
C          0: OK.
C          1: workspace too small for required # fine gridpoints
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPC, IP, IPTC, IPT, IPTO, IPTOLD, IRC, IR, IROLD,
     +   NPLNSC, NPLNS, NROWS
      LOGICAL LEFT, MIDDLE, RIGHT
C
      IERR = 0
      IDOM(0) = 0
      DO 1 IPT = 1, MAXPTS
         IDOM(IPT) = -1
    1 CONTINUE
C
      NPLNSC = LPLNC(0)
C
Ccc Create new grid level, all cells with at least one flagged corner
C will be refined
      IP  = 0
      IR  = 0
      IPT = 0
      DO 10 IPC = 1, NPLNSC
C
Ccc Make coarse grid plane of fine grid
         IROLD = IR
         DO 100 IRC = LPLNC(IPC), LPLNC(IPC+1)-1
C
Ccc Make coarse grid row of fine grid
            IPTOLD = IPT
            IPTO   = IPT
            LEFT = .FALSE.
            IPTC = LROWC(IRC)
            MIDDLE = REFFLG(LBLWYC(LBLWZC(IPTC))) .OR.
     +         REFFLG(LBLWZC(IPTC)) .OR. REFFLG(LABVYC(LBLWZC(IPTC)))
     +               .OR. REFFLG(LBLWYC(IPTC)) .OR.
     +         REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +               REFFLG(LBLWYC(LABVZC(IPTC))) .OR.
     +         REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
            DO 110 IPTC = LROWC(IRC)+1, LROWC(IRC+1)-1
               RIGHT = REFFLG(LBLWYC(LBLWZC(IPTC))) .OR.
     +            REFFLG(LBLWZC(IPTC)) .OR. REFFLG(LABVYC(LBLWZC(IPTC)))
     +                 .OR. REFFLG(LBLWYC(IPTC)) .OR.
     +            REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +                 REFFLG(LBLWYC(LABVZC(IPTC))) .OR.
     +            REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
               IF (MIDDLE .OR. RIGHT) THEN
C          Refine cell
                  IF (IPT+2 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC),IPLNC(IPC)) +
C                   its right neighbor
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
                  IF (IDOMC(IPTC-1) .NE. 0) THEN
C             Coarse grid point is at (physical) boundary, so is new
                     IDOM(IPT) = IDOMC(IPTC-1)
                  ELSE IF (IPT .EQ. IPTO+1) THEN
C             First new point in this (sub)row: internal boundary
                     IDOM(IPT) = NBNDS+1
                  ELSE
C             Internal point, or internal X-, or Z- boundary
                     IDOM(IPT) = 0
                  ENDIF
                  IPT = IPT + 1
                  ICOL(IPT) = ICOL(IPT-1)+1
C             If one of both (coarse) neighbors is an internal point, so
C             is IPT; otherwise it lies on a X-, or Z- boundary and
C             only one of the neighbors can be a physical corner
                  IF (IDOMC(IPTC-1) .EQ. 0 .OR.
     +                IDOMC(IPTC) .EQ. 0) THEN
                     IDOM(IPT) = 0
                  ELSE IF (IDOMC(IPTC-1) .EQ. NBNDS+1 .OR.
     +                     IDOMC(IPTC) .EQ. NBNDS+1) THEN
                     IDOM(IPT) = NBNDS+1
                  ELSE
                     IDOM(IPT) = MIN(IDOMC(IPTC-1),IDOMC(IPTC))
                  ENDIF
               ELSE IF (LEFT) THEN
C          Previous cell is refined, current not
                  IF (IPT+1 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC),IPLNC(IPC))
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
                  IF (IDOMC(IPTC-1) .NE. 0) THEN
C             Coarse grid point is at (physical) boundary, so is new
                     IDOM(IPT) = IDOMC(IPTC-1)
                  ELSE
C             Internal boundary
                     IDOM(IPT) = NBNDS+1
                  ENDIF
                  IPTO = IPT
               ENDIF
               LEFT = MIDDLE
               MIDDLE = RIGHT
  110       CONTINUE
            IF (LEFT) THEN
C       Last cell in row has been refined
               IF (IPT+1 .GT. MAXPTS) GOTO 900
C       Add last coarse node
               IPTC = LROWC(IRC+1)-1
               IPT = IPT + 1
               ICOL(IPT) = 2*ICOLC(IPTC)
C       Coarse grid point is at physical or internal boundary, so is new
               IDOM(IPT) = IDOMC(IPTC)
            ENDIF
C
            IF (IPT .GT. IPTOLD) THEN
C       Current coarse grid row has been refined
               IR = IR + 1
               LROW(IR) = IPTOLD + 1
               IROW(IR) = 2*IROWC(IRC)
            ENDIF
C
            IF (IRC .EQ. LPLNC(IPC+1)-1) GOTO 100
C
Ccc Make intermediate row of fine grid
            IPTOLD = IPT
            LEFT = .FALSE.
            IPTC = LROWC(IRC)
            MIDDLE =
     +         REFFLG(LBLWZC(IPTC)) .OR. REFFLG(LABVYC(LBLWZC(IPTC)))
     +         .OR. REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +         REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
            DO 120 IPTC = LROWC(IRC)+1, LROWC(IRC+1)-1
               RIGHT =
     +            REFFLG(LBLWZC(IPTC)) .OR. REFFLG(LABVYC(LBLWZC(IPTC)))
     +            .OR. REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +            REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
               IF (MIDDLE .OR. RIGHT) THEN
C          Refine cell
                  IF (IPT+2 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC)+1/2,IPLNC(IPC)) +
C                   its right neighbor
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
                  IPT = IPT + 1
                  ICOL(IPT) = ICOL(IPT-1)+1
               ELSE IF (LEFT) THEN
C          Previous cell is refined, current not
                  IF (IPT+1 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC)+1/2,IPLNC(IPC))
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
               ENDIF
               LEFT = MIDDLE
               MIDDLE = RIGHT
  120       CONTINUE
            IF (LEFT) THEN
C       Last cell in row has been refined
               IF (IPT+1 .GT. MAXPTS) GOTO 900
C       Add last coarse node
               IPTC = LROWC(IRC+1)-1
               IPT = IPT + 1
               ICOL(IPT) = 2*ICOLC(IPTC)
            ENDIF
C
            IF (IPT .GT. IPTOLD) THEN
C       Current intermediate row has been refined
               IR = IR + 1
               LROW(IR) = IPTOLD + 1
               IROW(IR) = 2*IROWC(IRC)+1
            ENDIF
  100    CONTINUE
C
         IF (IR .GT. IROLD) THEN
C       Current coarse grid plane has been refined
            IP = IP + 1
            LPLN(IP) = IROLD + 1
            IPLN(IP) = 2*IPLNC(IPC)
         ENDIF
         IF (IPC .EQ. NPLNSC) GOTO 10
C
Ccc Make intermediate grid plane of fine grid
         IROLD = IR
         DO 200 IRC = LPLNC(IPC), LPLNC(IPC+1)-1
C
Ccc Make coarse grid row of fine grid
            IPTOLD = IPT
            LEFT = .FALSE.
            IPTC = LROWC(IRC)
            MIDDLE = REFFLG(LBLWYC(IPTC)) .OR.
     +         REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +               REFFLG(LBLWYC(LABVZC(IPTC))) .OR.
     +         REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
            DO 210 IPTC = LROWC(IRC)+1, LROWC(IRC+1)-1
               RIGHT = REFFLG(LBLWYC(IPTC)) .OR.
     +            REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +                  REFFLG(LBLWYC(LABVZC(IPTC))) .OR.
     +            REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
               IF (MIDDLE .OR. RIGHT) THEN
C          Refine cell
                  IF (IPT+2 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC),IPNC(IPC)+1/2) +
C                   its right neighbor
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
                  IPT = IPT + 1
                  ICOL(IPT) = ICOL(IPT-1)+1
               ELSE IF (LEFT) THEN
C          Previous cell is refined, current not
                  IF (IPT+1 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC),IPNC(IPC)+1/2)
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
               ENDIF
               LEFT = MIDDLE
               MIDDLE = RIGHT
  210       CONTINUE
            IF (LEFT) THEN
C       Last cell in row has been refined
               IF (IPT+1 .GT. MAXPTS) GOTO 900
C       Add last coarse node
               IPTC = LROWC(IRC+1)-1
               IPT = IPT + 1
               ICOL(IPT) = 2*ICOLC(IPTC)
            ENDIF
C
            IF (IPT .GT. IPTOLD) THEN
C       Current coarse grid row has been refined
               IR = IR + 1
               LROW(IR) = IPTOLD + 1
               IROW(IR) = 2*IROWC(IRC)
            ENDIF
C
            IF (IRC .EQ. LPLNC(IPC+1)-1) GOTO 200
C
Ccc Make intermediate row of fine grid
            IPTOLD = IPT
            LEFT = .FALSE.
            IPTC = LROWC(IRC)
            MIDDLE =
     +         REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +         REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
            DO 220 IPTC = LROWC(IRC)+1, LROWC(IRC+1)-1
               RIGHT =
     +            REFFLG(IPTC) .OR. REFFLG(LABVYC(IPTC)) .OR.
     +            REFFLG(LABVZC(IPTC)) .OR. REFFLG(LABVYC(LABVZC(IPTC)))
               IF (MIDDLE .OR. RIGHT) THEN
C          Refine cell
                  IF (IPT+2 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC)+1/2,IPNC(IPC)+1/2) +
C                   its right neighbor
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
                  IPT = IPT + 1
                  ICOL(IPT) = ICOL(IPT-1)+1
               ELSE IF (LEFT) THEN
C          Previous cell is refined, current not
                  IF (IPT+1 .GT. MAXPTS) GOTO 900
C          Add node (ICOLC(IPTC-1),IROWC(IRC)+1/2,IPNC(IPC)+1/2)
                  IPT = IPT + 1
                  ICOL(IPT) = 2*ICOLC(IPTC-1)
               ENDIF
               LEFT = MIDDLE
               MIDDLE = RIGHT
  220       CONTINUE
            IF (LEFT) THEN
C       Last cell in row has been refined
               IF (IPT+1 .GT. MAXPTS) GOTO 900
C       Add last coarse node
               IPTC = LROWC(IRC+1)-1
               IPT = IPT + 1
               ICOL(IPT) = 2*ICOLC(IPTC)
            ENDIF
C
            IF (IPT .GT. IPTOLD) THEN
C       Current intermediate row has been refined
               IR = IR + 1
               LROW(IR) = IPTOLD + 1
               IROW(IR) = 2*IROWC(IRC)+1
            ENDIF
  200    CONTINUE
         IF (IR .GT. IROLD) THEN
C       Current intermediate grid plane has been refined
            IP = IP + 1
            LPLN(IP) = IROLD + 1
            IPLN(IP) = 2*IPLNC(IPC)+1
         ENDIF
   10 CONTINUE
C
Ccc Store # find grid planes in LPLN(0) and # fine grid points in NPTS
C and LROW(NROWS+1)
      NPLNS = IP
      NROWS = IR
      NPTS  = IPT
      LPLN(0) = NPLNS
      LPLN(NPLNS+1) = NROWS + 1
      LROW(NROWS+1) = NPTS + 1
C
      RETURN
C
Ccc Error return
  900 CONTINUE
      NPTS = MAXPTS
      IERR = 1
C
      RETURN
      END
      SUBROUTINE MKBND (NPTSF, LPLN, IPLN, LROW, IROW, ICOL,
     +   LBLWY, LABVY, LBLWZ, LABVZ,
     +   IDOM, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTSF
      INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(NPTSF),
     +   LBLWY(0:NPTSF), LABVY(0:NPTSF), LBLWZ(0:NPTSF), LABVZ(0:NPTSF),
     +   IDOM(0:NPTSF),  LLBND(0:*), ILBND(*), LBND(*)
C
Ccc PURPOSE:
C Make boundary list for refined grid using domain flags set on grid
C points corresponding to coarse grid points
C
Ccc PARAMETER DESCRIPTION:
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C IDOM   : INWORK. (0:NPTS)
C          IDOM(IPT): location in coarse grid of node IPT (only set for
C          rows corresponding with coarse grid rows)
C          0: interior point or new horizontal internal boundary
C          IB: on boundary
C              if not on edge or corner, or if internal boundary
C              then IB, IB = 1, NBNDS+1
C              else | IB4 | IB3 | IB2 | IBYTE | if point is part of
C                   the boundary planes IBYTE, IB2, IB3, and IB4
C LLBND  : OUT. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : OUT. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : OUT. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C
Ccc EXTERNALS USED: NONE
C
      INTEGER BYTE
      PARAMETER (BYTE = 2**8)
C
C-----------------------------------------------------------------------
C
      INTEGER IB, ID, ID1, IDA, IDB, IP, IPT, IR, NPLNS, NROWS,
     +   NPTS, NBNDS
C
      NPLNS = LPLN(0)
      NROWS = LPLN(NPLNS+1)-1
      NPTS  = LROW(NROWS+1)-1
      NBNDS = LLBND(0)
C
Ccc Domain flags have been set at nodes corresponding with coarse grid
C nodes and their X-neighbors, but some internal boundaries parallel to
C the X- or Z-axis can still be marked as internal points. Correct these
C by checking if any neighbor is missing. These points have X-neighbors.
C Correct the IDOM value for planes that do not longer exist in the
C refinement
      DO 10 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 10
         DO 20 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .NE. IROW(IR)) GOTO 20
            DO 30 IPT = LROW(IR), LROW(IR+1)-1
               IF (IDOM(IPT) .EQ. 0) THEN
                  IF (LBLWY(IPT-1)*LBLWY(IPT)*LBLWY(IPT+1) .EQ. 0 .OR.
     +                LABVY(IPT-1)*LABVY(IPT)*LABVY(IPT+1) .EQ. 0 .OR.
     +                LBLWZ(IPT-1)*LBLWZ(IPT)*LBLWZ(IPT+1) .EQ. 0 .OR.
     +                LABVZ(IPT-1)*LABVZ(IPT)*LABVZ(IPT+1) .EQ. 0 .OR.
     +                LBLWZ(LBLWY(IPT-1))*LBLWZ(LBLWY(IPT))*
     +                                 LBLWZ(LBLWY(IPT+1)) .EQ. 0 .OR.
     +                LABVZ(LBLWY(IPT-1))*LABVZ(LBLWY(IPT))*
     +                                 LABVZ(LBLWY(IPT+1)) .EQ. 0 .OR.
     +                LBLWZ(LABVY(IPT-1))*LBLWZ(LABVY(IPT))*
     +                                 LBLWZ(LABVY(IPT+1)) .EQ. 0 .OR.
     +                LABVZ(LABVY(IPT-1))*LABVZ(LABVY(IPT))*
     +                                 LABVZ(LABVY(IPT+1)) .EQ. 0)
     +                IDOM(IPT) = NBNDS+1
               ELSE IF (IDOM(IPT) .GT. BYTE) THEN
                  ID  = 0
                  ID1 = MOD(IDOM(IPT),BYTE)
                  IB  = ILBND(ID1)
                  IF (IB .EQ. 1 .AND. ICOL(IPT)+1 .EQ. ICOL(IPT+1)) THEN
                     ID = ID1
                  ELSE IF (IB .EQ. 2 .AND. LABVZ(IPT) .NE. 0) THEN
                            ID = ID1
                  ELSE IF (IB .EQ. 3 .AND. ICOL(IPT)-1 .EQ. ICOL(IPT-1))
     +            THEN
                     ID = ID1
                  ELSE IF (IB .EQ. 4 .AND. LBLWZ(IPT) .NE. 0) THEN
                     ID = ID1
                  ELSE IF (IB .EQ. 5 .AND. LABVY(IPT) .NE. 0) THEN
                     ID = ID1
                  ELSE IF (IB .EQ. 6 .AND. LBLWY(IPT) .NE. 0) THEN
                     ID = ID1
                  ENDIF
                  IDOM(IPT) = IDOM(IPT)/BYTE
                  ID1 = MOD(IDOM(IPT),BYTE)
                  IB  = ILBND(ID1)
                  IF (IB .EQ. 1 .AND. ICOL(IPT)+1 .EQ. ICOL(IPT+1)) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 2 .AND. LABVZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 3 .AND. ICOL(IPT)-1 .EQ. ICOL(IPT-1))
     +            THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 4 .AND. LBLWZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 5 .AND. LABVY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 6 .AND. LBLWY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ENDIF
                  IDOM(IPT) = IDOM(IPT)/BYTE
                  IF (IDOM(IPT) .EQ. 0) GOTO 40
                  ID1 = MOD(IDOM(IPT),BYTE)
                  IB  = ILBND(ID1)
                  IF (IB .EQ. 1 .AND. ICOL(IPT)+1 .EQ. ICOL(IPT+1)) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 2 .AND. LABVZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 3 .AND. ICOL(IPT)-1 .EQ. ICOL(IPT-1))
     +            THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 4 .AND. LBLWZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 5 .AND. LABVY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 6 .AND. LBLWY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ENDIF
                  IDOM(IPT) = IDOM(IPT)/BYTE
                  IF (IDOM(IPT) .EQ. 0) GOTO 40
                  ID1 = MOD(IDOM(IPT),BYTE)
                  IB  = ILBND(ID1)
                  IF (IB .EQ. 1 .AND. ICOL(IPT)+1 .EQ. ICOL(IPT+1)) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 2 .AND. LABVZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 3 .AND. ICOL(IPT)-1 .EQ. ICOL(IPT-1))
     +            THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 4 .AND. LBLWZ(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 5 .AND. LABVY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ELSE IF (IB .EQ. 6 .AND. LBLWY(IPT) .NE. 0) THEN
                     ID = ID*BYTE+ID1
                  ENDIF
   40             IDOM(IPT) = ID
               ENDIF
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
Ccc Set domain flags in intermediate rows of coarse grid planes
      DO 100 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 100
         DO 110 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .EQ. IROW(IR)) GOTO 110
            DO 120 IPT = LROW(IR), LROW(IR+1)-1
               IDA = IDOM(LABVY(IPT))
               IDB = IDOM(LBLWY(IPT))
C          If one of both neighbors is an internal point, so
C          is IPT; otherwise it lies on a Y- or Z-boundary,
C          if both Y-neighbors are lying on a plane only one can be a
C          physical edge, if both are lying on an edge only one can be a
C          physical corner
               IF (IDA .EQ. 0 .OR. IDB .EQ. 0) THEN
                  IDOM(IPT) = 0
               ELSE IF (IDA .EQ. NBNDS+1 .OR. IDB .EQ. NBNDS+1) THEN
                  IDOM(IPT) = NBNDS+1
               ELSE
                  IDOM(IPT) = MIN(IDA,IDB)
               ENDIF
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
Ccc Set domain flags in intermediate planes
      DO 150 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .EQ. IPLN(IP)) GOTO 150
         DO 160 IR = LPLN(IP), LPLN(IP+1)-1
            DO 170 IPT = LROW(IR), LROW(IR+1)-1
               IDA = IDOM(LABVZ(IPT))
               IDB = IDOM(LBLWZ(IPT))
C          If one of both neighbors is an internal point, so
C          is IPT; otherwise it lies on a X- or Y-boundary,
C          if both Z-neighbors are lying on a plane only one can be a
C          physical edge, if both are lying on an edge only one can be a
C          physical corner
               IF (IDA .EQ. 0 .OR. IDB .EQ. 0) THEN
                  IDOM(IPT) = 0
               ELSE IF (IDA .EQ. NBNDS+1 .OR. IDB .EQ. NBNDS+1) THEN
                  IDOM(IPT) = NBNDS+1
               ELSE
                  IDOM(IPT) = MIN(IDA,IDB)
               ENDIF
  170       CONTINUE
  160    CONTINUE
  150 CONTINUE
C
Ccc Edges between physical and internal boundaries can still wrongly be
C marked physical
      DO 200 IP = 1, NPLNS
         DO 210 IR = LPLN(IP), LPLN(IP+1)-1
            DO 220 IPT = LROW(IR), LROW(IR+1)-1
               ID = IDOM(IPT)
               IF (ID .EQ. 0 .OR. ID .EQ. NBNDS+1 .OR.
     +             ID .GT. BYTE) GOTO 220
               IB = ILBND(ID)
               IF (IB .EQ. 1) THEN
                  ID1 = IDOM(IPT+1)
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ELSE IF (IB .EQ. 2) THEN
                  ID1 = IDOM(LABVZ(IPT))
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ELSE IF (IB .EQ. 3) THEN
                  ID1 = IDOM(IPT-1)
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ELSE IF (IB .EQ. 4) THEN
                  ID1 = IDOM(LBLWZ(IPT))
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ELSE IF (IB .EQ. 5) THEN
                  ID1 = IDOM(LABVY(IPT))
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ELSE IF (IB .EQ. 6) THEN
                  ID1 = IDOM(LBLWY(IPT))
                  IF (ID1 .EQ. NBNDS+1) IDOM(IPT) = NBNDS+1
               ENDIF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
Ccc Corners can still wrongly be marked physical. If a point has at
C least 2 neighbors that are internal boundary points, so is the
C point itself
      DO 230 IP = 1, NPLNS
         DO 240 IR = LPLN(IP), LPLN(IP+1)-1
            IPT = LROW(IR)
               IF (IDOM(IPT) .LT. BYTE) GOTO 249
               IB = 0
               IF (IDOM(IPT+1)      .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IB .GE. 2) IDOM(IPT) = NBNDS+1
  249       DO 250 IPT = LROW(IR)+1, LROW(IR+1)-2
               IF (IDOM(IPT) .LT. BYTE) GOTO 250
               IB = 0
               IF (ICOL(IPT-1) .EQ. ICOL(IPT)-1 .AND.
     +             IDOM(IPT-1) .EQ. NBNDS+1) IB = IB+1
               IF (ICOL(IPT+1) .EQ. ICOL(IPT)+1 .AND.
     +             IDOM(IPT+1) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IB .GE. 2) IDOM(IPT) = NBNDS+1
  250       CONTINUE
            IPT = LROW(IR+1)-1
               IF (IDOM(IPT) .LT. BYTE) GOTO 240
               IB = 0
               IF (IDOM(IPT-1)      .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVY(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LBLWZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IDOM(LABVZ(IPT)) .EQ. NBNDS+1) IB = IB+1
               IF (IB .GE. 2) IDOM(IPT) = NBNDS+1
  240    CONTINUE
  230 CONTINUE

Ccc Make boundary lists
      DO 300 IB = 0, NBNDS+1
         LLBND(IB) = 0
  300 CONTINUE
      DO 310 IPT = 1, NPTS
         ID = IDOM(IPT)
         ID1 = MOD(ID,BYTE)
         LLBND(ID1) = LLBND(ID1) + 1
         IF (ID .EQ. ID1) GOTO 310
         ID = ID/BYTE
         ID1 = MOD(ID,BYTE)
         LLBND(ID1) = LLBND(ID1) + 1
         IF (ID .EQ. ID1) GOTO 310
         ID = ID/BYTE
         ID1 = MOD(ID,BYTE)
         LLBND(ID1) = LLBND(ID1) + 1
         IF (ID .EQ. ID1) GOTO 310
         ID = ID/BYTE
         LLBND(ID) = LLBND(ID) + 1
  310 CONTINUE
      LLBND(0) = 1
      DO 320 IB = 1, NBNDS
         LLBND(IB) = LLBND(IB-1) + LLBND(IB)
  320 CONTINUE
      DO 330 IB = NBNDS, 0, -1
         LLBND(IB+2) = LLBND(IB)
  330 CONTINUE
      DO 340 IPT = 1, NPTS
         IF (IDOM(IPT) .EQ. 0) GOTO 340
         ID = IDOM(IPT)
         ID1 = MOD(ID,BYTE)
         LBND(LLBND(ID1+1)) = IPT
         LLBND(ID1+1) = LLBND(ID1+1) + 1
         IF (ID .EQ. ID1) GOTO 340
         ID = ID/BYTE
         ID1 = MOD(ID,BYTE)
         LBND(LLBND(ID1+1)) = IPT
         LLBND(ID1+1) = LLBND(ID1+1) + 1
         IF (ID .EQ. ID1) GOTO 340
         ID = ID/BYTE
         ID1 = MOD(ID,BYTE)
         LBND(LLBND(ID1+1)) = IPT
         LLBND(ID1+1) = LLBND(ID1+1) + 1
         IF (ID .EQ. ID1) GOTO 340
         ID = ID/BYTE
         LBND(LLBND(ID+1)) = IPT
         LLBND(ID+1) = LLBND(ID+1) + 1
  340 CONTINUE
      LLBND(0) = NBNDS
      LLBND(1) = 1

      RETURN
      END
      SUBROUTINE GETSOL (NPDE, UC, ISTRCN, LEVO, UO, ISTRFO, U, ISTRUC,
     +   IWORK, RWORK)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, ISTRCN(0:*), ISTRFO(0:*), ISTRUC(0:*), IWORK(*)
      LOGICAL LEVO
      DOUBLE PRECISION UC(0:*), UO(0:*), U(0:*), RWORK(*)
C
Ccc PURPOSE:
C Store solution at a previous time level on a grid of the current time
C level.
C First copy the (embedded) coarser grid solution, then copy all
C available values from the old time grid of the same grid level, and
C finally calculate all non-initialized values by interpolation.
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C UC     : IN. Solution at embedded coarse grid
C ISTRCN : IN. Datastructure of the embedded coarse grid
C LEVO   : IN. .TRUE. if new grid level existed at previous time level
C UO     : IN. If LEVO the solution at grid ISTRFO
C ISTRFO : IN. If LEVO the datastructure of the grid with the same level
C          as the new grid level but on a previous time level
C U      : OUT. Solution of previous time level on new grid level
C ISTRUC : IN. Data structure of new grid level at current time level
C IWORK  : WORK. (NPTS)
C RWORK  : WORK. (0)
C
Ccc EXTERNALS USED:
      EXTERNAL INJCF, INJON, INTPOL
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLNC,NPLNSC,LIPLNC, LLROWC,NROWSC,NPTSC,LIROWC, LICOLC,
     +   LLPLNO,NPLNSO,LIPLNO, LLROWO, NROWSO, NPTSO, LIROWO, LICOLO,
     +   LLPLN,NPLNS,LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ
C
      LLPLNC = 0
      NPLNSC = ISTRCN(LLPLNC)
      NROWSC = ISTRCN(LLPLNC+NPLNSC+1)-1
      LIPLNC = LLPLNC+NPLNSC+2
      LLROWC = LIPLNC+NPLNSC
      NPTSC  = ISTRCN(LLROWC+NROWSC)-1
      LIROWC = LLROWC+NROWSC+1
      LICOLC = LIROWC+NROWSC
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
C
      U(0) = 0.0
C
C Copy embedded coarse grid solution
      CALL INJCF (NPDE, UC, U, IWORK,
     +   ISTRCN(LLPLNC), ISTRCN(LIPLNC),
     +   ISTRCN(LLROWC), ISTRCN(LIROWC), ISTRCN(LICOLC),
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL))
C
      IF (LEVO) THEN
C Copy solution of grid with same level but on previous time level
         LLPLNO = 0
         NPLNSO = ISTRFO(LLPLNO)
         NROWSO = ISTRFO(LLPLNO+NPLNSO+1)-1
         LIPLNO = LLPLNO+NPLNSO+2
         LLROWO = LIPLNO+NPLNSO
         NPTSO  = ISTRFO(LLROWO+NROWSO)-1
         LIROWO = LLROWO+NROWSO+1
         LICOLO = LIROWO+NROWSO
         CALL INJON (NPDE, UO, U, IWORK,
     +      ISTRFO(LLPLNO), ISTRFO(LIPLNO),
     +      ISTRFO(LLROWO), ISTRFO(LIROWO), ISTRFO(LICOLO),
     +      ISTRUC(LLPLN), ISTRUC(LIPLN),
     +      ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL))
      ENDIF
C
C Calculate all uninitialized values by interpolation
      CALL INTPOL (NPDE, U, IWORK,
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +   ISTRUC(LLBLWZ), ISTRUC(LLABVZ), RWORK)
C
      RETURN
      END
      SUBROUTINE GETINI (NPDE, UIC, UC, ISTRCC, LEVO, UO, ISTRFO,
     +   U, ISTRUC, UIB, IWORK, RWORK)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, ISTRCC(0:*), ISTRFO(0:*), ISTRUC(0:*), IWORK(*)
      LOGICAL LEVO
      DOUBLE PRECISION UIC(0:*), UC(0:*), UO(0:*), U(0:*), UIB(*),
     +   RWORK(*)
C
Ccc PURPOSE:
C Initialize solution at current time level on the next finer grid
C level. Store (Dirichlet) internal boundary values in UIB.
C First get the internal bounday values from interpolation of the
C solution at the (embedded) coarser grid.
C The initial solution is obtained by first copying the (embedded)
C coarser grid initial solution, and then all available values from the
C solution at the previous time level on the grid of the required grid
C level. Finally all as yet non-initialized values are calculated by
C interpolation.
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C UIC    : IN. Initial solution at embedded coarse grid
C UC     : IN. Solution at embedded coarse grid
C ISTRCC : IN. Datastructure of the embedded coarse grid
C LEVO   : IN. .TRUE. if new grid level existed at previous time level
C UO     : IN. If LEVO the solution at grid ISTRFO
C ISTRFO : IN. If LEVO the datastructure of the grid with the same level
C          as the new grid level but on a previous time level
C U      : OUT. Solution of current time level on new grid level
C ISTRUC : IN. Data structure of new grid level at current time level
C UIB    : OUT. List of internal boundary values
C IWORK  : WORK. (NPTS)
C RWORK  : WORK. (0)
C
Ccc EXTERNALS USED:
      EXTERNAL INJCF, INJCFB, INJON, INTPOL
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLNC,NPLNSC,LIPLNC, LLROWC,NROWSC,NPTSC,LIROWC, LICOLC,
     +   LLPLNO, NPLNSO, LIPLNO, LLROWO, NROWSO, NPTSO, LIROWO, LICOLO,
     +   LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ,
     +   NIBPTS, I, IB, IC, LB
C
      LLPLNC = 0
      NPLNSC = ISTRCC(LLPLNC)
      NROWSC = ISTRCC(LLPLNC+NPLNSC+1)-1
      LIPLNC = LLPLNC+NPLNSC+2
      LLROWC = LIPLNC+NPLNSC
      NPTSC  = ISTRCC(LLROWC+NROWSC)-1
      LIROWC = LLROWC+NROWSC+1
      LICOLC = LIROWC+NROWSC
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBDPTS = ISTRUC(LLLBND+NBNDS+1)-1
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
      NIBPTS = NBIPTS - NBDPTS
C
      U(0) = 0.0
C
Ccc Get internal boundary values
C Copy embedded coarse grid solution
      CALL INJCFB (NPDE, UC, U, IWORK,
     +   ISTRCC(LLPLNC), ISTRCC(LIPLNC),
     +   ISTRCC(LLROWC), ISTRCC(LIROWC), ISTRCC(LICOLC),
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL),
     +   ISTRUC(LLLBND), ISTRUC(LLBNDP))
C
C Calculate all uninitialized values at the internal boundary by
C interpolation
      CALL INTPOL (NPDE, U, IWORK,
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +   ISTRUC(LLBLWZ), ISTRUC(LLABVZ), RWORK)
C
C Store internal boundary values in list
      DO 10 IC = 1, NPDE
      DO 10 LB = 1, NIBPTS
          I  = ISTRUC(LLBNDP+NBDPTS-1+LB) + (IC-1)*NPTS
         IB = LB + (IC-1)*NIBPTS
         UIB(IB) = U(I)
   10 CONTINUE
C
Ccc Get initial solution
C Copy embedded coarse grid initial solution
      CALL INJCF (NPDE, UIC, U, IWORK,
     +   ISTRCC(LLPLNC), ISTRCC(LIPLNC),
     +   ISTRCC(LLROWC), ISTRCC(LIROWC), ISTRCC(LICOLC),
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL))
C
      IF (LEVO) THEN
C
C Copy solution of grid with same level but on previous time level
         LLPLNO = 0
         NPLNSO = ISTRFO(LLPLNO)
         NROWSO = ISTRFO(LLPLNO+NPLNSO+1)-1
         LIPLNO = LLPLNO+NPLNSO+2
         LLROWO = LIPLNO+NPLNSO
         NPTSO  = ISTRFO(LLROWO+NROWSO)-1
         LIROWO = LLROWO+NROWSO+1
         LICOLO = LIROWO+NROWSO
         CALL INJON (NPDE, UO, U, IWORK,
     +      ISTRFO(LLPLNO), ISTRFO(LIPLNO),
     +      ISTRFO(LLROWO), ISTRFO(LIROWO), ISTRFO(LICOLO),
     +      ISTRUC(LLPLN), ISTRUC(LIPLN),
     +      ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL))
C
      ENDIF
C
C Calculate all uninitialized values by interpolation
      CALL INTPOL (NPDE, U, IWORK,
     +   ISTRUC(LLPLN), ISTRUC(LIPLN),
     +   ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +   ISTRUC(LLBLWZ), ISTRUC(LLABVZ), RWORK)
C
      RETURN
      END
      SUBROUTINE INJCFB (NPDE, UC, U, IPDOM,
     +   LPLNC, IPLNC, LROWC, IROWC, ICOLC,
     +   LPLN, IPLN, LROW, IROW, ICOL, LLBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, IPDOM(*),
     +        LPLNC(0:*), IPLNC(*), LROWC(*), IROWC(*), ICOLC(*),
     +        LPLN(0:*), IPLN(*), LROW(*),  IROW(*),  ICOL(*),
     +        LLBND(0:*), LBND(*)
      DOUBLE PRECISION UC(0:*), U(0:*)
C
Ccc PURPOSE:
C Inject solution from coarse grid into (embedded) fine grid
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C UC     : IN. Solution at coarse grid
C U      : OUT. Solution at coarse gridpoints at internal boundary of
C               refined grid
C IPDOM  : OUT. Domain flags wrt to interpolation
C                0: Injected point
C               -1: Otherwise
C LPLNC  : IN. -I
C IPLNC  : IN.  I
C LROWC  : IN.  I Data structure of the coarse grid
C IROWC  : IN.  I see description for fine grid below
C ICOLC  : IN. -I
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LLBND  : (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary or corner in
C                          LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C LBND   : OUT. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IB, IC, IPC, IP, IPTC, IPT, IRC, IR,
     +   NBNDS, NPLNSC, NPLNS, NPTSC, NPTS, NROWSC, NROWS
C
      NPLNS  = LPLN(0)
      NROWS  = LPLN(NPLNS+1)-1
      NPTS   = LROW(NROWS+1)-1
      NBNDS  = LLBND(0)
      NPLNSC = LPLNC(0)
      NROWSC = LPLNC(NPLNSC+1)-1
      NPTSC  = LROWC(NROWSC+1)-1
C
Ccc Initialize interpolation flags
      DO 10 IPT = 1, NPTS
         IPDOM(IPT) = 0
   10 CONTINUE
      DO 20 IB = LLBND(1), LLBND(NBNDS+2)-1
         IPT = LBND(IB)
         IPDOM(IPT) = -1
   20 CONTINUE
C
Ccc Inject values from coarse level into fine grid solution
CDIR$ NOVECTOR
      IPC = 0
      DO 30 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 30
         IPC = IPC + 1
   40    IF (2*IPLNC(IPC) .NE. IPLN(IP)) THEN
            IPC = IPC + 1
            GOTO 40
         ENDIF
         IRC = LPLNC(IPC)-1
         DO 50 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .NE. IROW(IR)) GOTO 50
            IRC = IRC + 1
   60       IF (2*IROWC(IRC) .NE. IROW(IR)) THEN
               IRC = IRC + 1
               GOTO 60
            ENDIF
            IPTC = LROWC(IRC)-1
            DO 70 IPT = LROW(IR), LROW(IR+1)-1
               IF (ICOL(IPT)/2*2 .NE. ICOL(IPT)) GOTO 70
               IPTC = IPTC + 1
   80          IF (2*ICOLC(IPTC) .NE. ICOL(IPT)) THEN
                  IPTC = IPTC + 1
                  GOTO 80
               ENDIF
               DO 90 IC = 1, NPDE
                  U(IPT+(IC-1)*NPTS) = UC(IPTC+(IC-1)*NPTSC)
   90          CONTINUE
               IPDOM(IPT) = 0
   70       CONTINUE
   50    CONTINUE
   30 CONTINUE

      RETURN
      END
      SUBROUTINE INJCF (NPDE, UC, U, IPDOM,
     +   LPLNC, IPLNC, LROWC, IROWC, ICOLC,
     +   LPLN, IPLN, LROW, IROW, ICOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, IPDOM(*),
     +        LPLNC(0:*), IPLNC(*), LROWC(*), IROWC(*), ICOLC(*),
     +        LPLN(0:*), IPLN(*), LROW(*),  IROW(*),  ICOL(*)
      DOUBLE PRECISION UC(0:*), U(0:*)
C
Ccc PURPOSE:
C Inject solution from coarse grid into (embedded) fine grid
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C UC     : IN. Solution at coarse grid
C U      : OUT. Solution at coarse gridpoints in refined grid
C IPDOM  : OUT. Domain flags wrt to interpolation
C                0: Injected point
C               -1: Otherwise
C LPLNC  : IN. -I
C IPLNC  : IN.  I
C LROWC  : IN.  I Data structure of the coarse grid
C IROWC  : IN.  I see description for fine grid below
C ICOLC  : IN. -I
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, IPC, IP, IPTC, IPT, IRC, IR,
     +   NPLNSC, NPLNS, NPTSC, NPTS, NROWSC, NROWS
C
      NPLNS  = LPLN(0)
      NROWS  = LPLN(NPLNS+1)-1
      NPTS   = LROW(NROWS+1)-1
      NPLNSC = LPLNC(0)
      NROWSC = LPLNC(NPLNSC+1)-1
      NPTSC  = LROWC(NROWSC+1)-1
C
Ccc Initialize interpolation flags
      DO 10 IPT = 1, NPTS
         IPDOM(IPT) = -1
   10 CONTINUE
C
Ccc Inject values from coarse level into fine grid solution
CDIR$ NOVECTOR
      IPC = 0
      DO 30 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 30
         IPC = IPC + 1
   40    IF (2*IPLNC(IPC) .NE. IPLN(IP)) THEN
            IPC = IPC + 1
            GOTO 40
         ENDIF
         IRC = LPLNC(IPC)-1
         DO 50 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .NE. IROW(IR)) GOTO 50
            IRC = IRC + 1
   60       IF (2*IROWC(IRC) .NE. IROW(IR)) THEN
               IRC = IRC + 1
               GOTO 60
            ENDIF
            IPTC = LROWC(IRC)-1
            DO 70 IPT = LROW(IR), LROW(IR+1)-1
               IF (ICOL(IPT)/2*2 .NE. ICOL(IPT)) GOTO 70
               IPTC = IPTC + 1
   80          IF (2*ICOLC(IPTC) .NE. ICOL(IPT)) THEN
                  IPTC = IPTC + 1
                  GOTO 80
               ENDIF
               DO 90 IC = 1, NPDE
                  U(IPT+(IC-1)*NPTS) = UC(IPTC+(IC-1)*NPTSC)
   90          CONTINUE
               IPDOM(IPT) = 0
   70       CONTINUE
   50    CONTINUE
   30 CONTINUE

      RETURN
      END
      SUBROUTINE INJON (NPDE, UO, U, IPDOM,
     +   LPLNO, IPLNO, LROWO, IROWO, ICOLO,
     +   LPLN, IPLN, LROW, IROW, ICOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, IPDOM(*),
     +        LPLNO(0:*), IPLNO(*), LROWO(*), IROWO(*), ICOLO(*),
     +        LPLN(0:*), IPLN(*), LROW(*),  IROW(*),  ICOL(*)
      DOUBLE PRECISION UO(0:*), U(0:*)
C
Ccc PURPOSE:
C Inject solution from previous time at grid from same level into
C solution at current time grid
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C UO     : IN. Solution at previous time
C U      : INOUT. Solution at corresponding gridpoints injected from UO
C IPDOM  : INOUT. Domain flags wrt to interpolation
C                0: Injected point
C               -1: Otherwise
C LPLNO  : IN. -I
C IPLNO  : IN.  I
C LROWO  : IN.  I Data structure of the old-time grid
C IROWO  : IN.  I see description for current time grid below
C ICOLO  : IN. -I
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
CDIR$ NOVECTOR
C
      INTEGER IC, IPO, IP, IPTO, IPT, IRO, IR,
     +   NPLNSO, NPLNS, NPTSO, NPTS, NROWSO, NROWS
C
      NPLNS  = LPLN(0)
      NROWS  = LPLN(NPLNS+1)-1
      NPTS   = LROW(NROWS+1)-1
      NPLNSO = LPLNO(0)
      NROWSO = LPLNO(NPLNSO+1)-1
      NPTSO  = LROWO(NROWSO+1)-1
C
Ccc Inject values from old time level into current-time solution
      IPO = 1
      DO 10 IP = 1, NPLNS
   20    IF (IPLNO(IPO) .LT. IPLN(IP)) THEN
            IPO = IPO + 1
            IF (IPO .LE. NPLNSO) GOTO 20
            RETURN
         ELSE IF (IPLNO(IPO) .GT. IPLN(IP)) THEN
            GOTO 10
         ENDIF
         IRO = LPLNO(IPO)
         DO 30 IR = LPLN(IP), LPLN(IP+1)-1
   40       IF (IROWO(IRO) .LT. IROW(IR)) THEN
               IRO = IRO + 1
               IF (IRO .LE. NROWSO) GOTO 40
               GOTO 10
            ELSE IF (IROWO(IRO) .GT. IROW(IR)) THEN
               GOTO 30
            ENDIF
            IPTO = LROWO(IRO)
            DO 50 IPT = LROW(IR), LROW(IR+1)-1
   60          IF (ICOLO(IPTO) .LT. ICOL(IPT)) THEN
                  IPTO = IPTO + 1
                  IF (IPTO .LE. LROWO(IRO+1)-1) GOTO 60
                  GOTO 30
               ELSE IF (ICOLO(IPTO) .GT. ICOL(IPT)) THEN
                  GOTO 50
               ENDIF
               DO 70 IC = 1, NPDE
                  U(IPT+(IC-1)*NPTS) = UO(IPTO+(IC-1)*NPTSO)
   70          CONTINUE
               IPDOM(IPT) = 0
   50       CONTINUE
   30    CONTINUE
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PUTSOL (NPDE, U, ISTRUC, UC, ISTRCC, UI, LENUC)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, ISTRUC(0:*), ISTRCC(0:*), LENUC
      DOUBLE PRECISION U(0:*), UC(0:*), UI(0:*)
C
Ccc PURPOSE:
C Copy coarse grid solution UC to UI and inject fine grid solution U
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C U      : IN. Solution at fine grid
C ISTRUC : IN. Datastructure of the fine grid
C UC     : IN. Solution at coarse grid
C ISTRCC : IN. Datastructure of the coarse grid
C UI     : OUT.Injected solution at coarse grid
C LENUC  : OUT.Dimension of coarse grid solution array
C
Ccc EXTERNALS USED:
      EXTERNAL INJFC, RCOPY
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, NROWS, LLROW, NPTS, LIROW, LICOL,
     +   LLPLNC, NPLNSC, LIPLNC, NROWSC, LLROWC, NPTSC, LIROWC, LICOLC
C
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      LLROW  = LIPLN+NPLNS
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLPLNC = 0
      NPLNSC = ISTRCC(LLPLNC)
      LIPLNC = LLPLNC+NPLNSC+2
      NROWSC = ISTRCC(LLPLNC+NPLNSC+1)-1
      LLROWC = LIPLNC+NPLNSC
      NPTSC  = ISTRCC(LLROWC+NROWSC)-1
      LIROWC = LLROWC+NROWSC+1
      LICOLC = LIROWC+NROWSC
      LENUC  = NPTSC*NPDE + 1
C
C Copy coarse grid solution to UI
      CALL RCOPY (LENUC, UC, UI)
C
C Inject fine grid solution into UI
      CALL INJFC (NPDE, U, UI,
     +   ISTRUC(LLPLN),  ISTRUC(LIPLN),
     +   ISTRUC(LLROW),  ISTRUC(LIROW),  ISTRUC(LICOL),
     +   ISTRCC(LLPLNC), ISTRCC(LIPLNC),
     +   ISTRCC(LLROWC), ISTRCC(LIROWC), ISTRCC(LICOLC))
      
      RETURN
      END
      SUBROUTINE INJFC (NPDE, U, UC,
     +   LPLN, IPLN, LROW, IROW, ICOL,
     +   LPLNC, IPLNC, LROWC, IROWC, ICOLC)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
     +   LPLNC(0:*), IPLNC(*), LROWC(*), IROWC(*), ICOLC(*)
      DOUBLE PRECISION U(0:*), UC(0:*)
C
Ccc PURPOSE:
C Inject solution from (embedded) fine grid into coarser grid
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C U      : IN. Fine grid solution
C UC     : INOUT.
C          IN: Coarse grid solution
C          OUT: Injected solution at coarse grid
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LPLNC  : IN. -I
C IPLNC  : IN.  I
C LROWC  : IN.  I Data structure of the coarse grid
C IROWC  : IN.  I see description for fine grid above
C ICOLC  : IN. -I
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
CDIR$ NOVECTOR
C
      INTEGER IC, IPC, IP, IPTC, IPT, IRC, IR,
     +   NPLNSC, NPLNS, NPTSC, NPTS, NROWSC, NROWS
C
      NPLNS  = LPLN(0)
      NROWS  = LPLN(NPLNS+1)-1
      NPTS   = LROW(NROWS+1)-1
      NPLNSC = LPLNC(0)
      NROWSC = LPLNC(NPLNSC+1)-1
      NPTSC  = LROWC(NROWSC+1)-1
C
Ccc Inject values from fine level into coarse grid solution
      IPC = 0
      DO 10 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 10
         IPC = IPC + 1
   20    IF (2*IPLNC(IPC) .NE. IPLN(IP)) THEN
            IPC = IPC + 1
            GOTO 20
         ENDIF
         IRC = LPLNC(IPC)-1
         DO 30 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .NE. IROW(IR)) GOTO 30
            IRC = IRC + 1
   40       IF (2*IROWC(IRC) .NE. IROW(IR)) THEN
               IRC = IRC + 1
               GOTO 40
            ENDIF
            IPTC = LROWC(IRC)-1
            DO 50 IPT = LROW(IR), LROW(IR+1)-1
               IF (ICOL(IPT)/2*2 .NE. ICOL(IPT)) GOTO 50
               IPTC = IPTC + 1
   60          IF (2*ICOLC(IPTC) .NE. ICOL(IPT)) THEN
                  IPTC = IPTC + 1
                  GOTO 60
               ENDIF
               DO 70 IC = 1, NPDE
                  UC(IPTC+(IC-1)*NPTSC) = U(IPT+(IC-1)*NPTS)
   70          CONTINUE
   50       CONTINUE
   30    CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE INTPOL (NPDE, U, IPDOM,
     +   LPLN, IPLN, LROW, IROW, ICOL, LBLWY, LABVY, LBLWZ, LABVZ, WORK)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPDE, IPDOM(*),
     +   LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
     +   LBLWY(0:*), LABVY(0:*), LBLWZ(0:*), LABVZ(0:*)
      DOUBLE PRECISION U(0:*), WORK(*)
C
Ccc PURPOSE:
C Interpolate where necessary solution
C
Ccc PARAMETER DESCRIPTION:
C NPDE   : IN. # PDE components
C U      : INOUT.
C          IN: Solution values at injected points
C          OUT: Interpolated solution values at other points
C IPDOM  : IN. Domain flags wrt to interpolation
C                0: Injected point
C               -1: Otherwise
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
CDIR$ NOVECTOR
C
      INTEGER IC, IP, IPT, IR, NPLNS, NPTS, NROWS
C
      NPLNS = LPLN(0)
      NROWS = LPLN(NPLNS+1)-1
      NPTS  = LROW(NROWS+1)-1
C
      DO 10 IP = 1, NPLNS
         IF (IPLN(IP)/2*2 .NE. IPLN(IP)) GOTO 10
C
Ccc Interpolation in coarse grid planes
C
C       Interpolation in X-direction
         DO 20 IR = LPLN(IP), LPLN(IP+1)-1
            IF (IROW(IR)/2*2 .NE. IROW(IR)) GOTO 20
            DO 30 IPT = LROW(IR)+1, LROW(IR+1)-2
               IF (IPDOM(IPT) .NE. 0) THEN
                  DO 40 IC = 1, NPDE
                     U(IPT+(IC-1)*NPTS) =
     +                  (U(IPT-1+(IC-1)*NPTS) + U(IPT+1+(IC-1)*NPTS))/2
   40             CONTINUE
               ENDIF
   30       CONTINUE
   20    CONTINUE
C
C       Interpolation in Y-direction
         DO 50 IR = LPLN(IP)+1, LPLN(IP+1)-2
            IF (IROW(IR)/2*2 .EQ. IROW(IR)) GOTO 50
            DO 60 IPT = LROW(IR), LROW(IR+1)-1
               IF (IPDOM(IPT) .NE. 0) THEN
                  DO 70 IC = 1, NPDE
                     U(IPT+(IC-1)*NPTS) =
     +                  (U(LBLWY(IPT)+(IC-1)*NPTS) +
     +                   U(LABVY(IPT)+(IC-1)*NPTS))/2
   70             CONTINUE
               ENDIF
   60       CONTINUE
   50    CONTINUE
   10 CONTINUE

      DO 100 IP = 2, NPLNS-1
         IF (IPLN(IP)/2*2 .EQ. IPLN(IP)) GOTO 100
C
Ccc Interpolation in other then coarse grid planes
C
C       Interpolation in Z-direction
         DO 110 IR = LPLN(IP), LPLN(IP+1)-1
            DO 120 IPT = LROW(IR), LROW(IR+1)-1
               IF (IPDOM(IPT) .NE. 0) THEN
                  DO 130 IC = 1, NPDE
                     U(IPT+(IC-1)*NPTS) =
     +                  (U(LBLWZ(IPT)+(IC-1)*NPTS) +
     +                   U(LABVZ(IPT)+(IC-1)*NPTS))/2
  130             CONTINUE
               ENDIF
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

      RETURN
      END
      SUBROUTINE RESID (T, X, Y, Z, NPTS, NPDE, UNP1, UN, UNM1,
     +   DT, DTRAT, UIB, LLBND, ILBND, LBND, LBLWY, LABVY, LBLWZ, LABVZ,
     +   DX, DY, DZ, UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, F)
C
C-----------------------------------------------------------------------
C
C PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*),
     +   LBLWY(0:*), LABVY(0:*), LBLWZ(0:*), LABVZ(0:*)
      DOUBLE PRECISION T, X(*), Y(*), Z(*), UNP1(0:*), UN(0:*),
     +   UNM1(0:*),
     +   DT, DTRAT, UIB(*), DX, DY, DZ, UT(*), UX(*), UY(*), UZ(*),
     +   UXX(*), UYY(*), UZZ(*), UXY(*), UXZ(*), UYZ(*), F(*)
C
C PURPOSE:
C Compute time and space derivatives of U and residual F(t,Un+1,Udot)
C
C PARAMETER DESCRIPTION:
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C UNP1   : IN. Solution at Tn+1 on current grid
C UN     : IN. Solution at Tn   on current grid
C UNM1   : IN. Solution at Tn-1 on current grid
C DT     : IN. Current time stepsize
C DTRAT  : IN. 0 or DT/DT_old
C UIB    : IN. Solution at T on internal boundaries
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C UT     : OUT. Time derivative of U on current grid
C UX     : OUT. -I
C UY     : OUT.  I
C UZ     : OUT.  I
C UXX    : OUT.  I Space derivatives of U on current grid
C UYY    : OUT.  I
C UZZ    : OUT.  I
C UXY    : OUT.  I
C UXZ    : OUT.  I
C UYZ    : OUT. -I
C F      : OUT. Residual
C
Ccc EXTERNALS USED:
      EXTERNAL DERIVS, DERIVT, RES
C
C-----------------------------------------------------------------------
C
Ccc Compute derivatives
      CALL DERIVT (NPTS, NPDE, UNP1(1), UN(1), UNM1(1), DT, DTRAT, UT)
      CALL DERIVS (NPTS, NPDE, UNP1, LLBND, ILBND, LBND,
     +   LBLWY, LABVY, LBLWZ, LABVZ, DX, DY, DZ,
     +   UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ)
C
Ccc Compute residual
      CALL RES (T, X, Y, Z, NPTS, NPDE, UNP1(1), LLBND, ILBND, LBND,
     +   UIB, UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, F)
C
      RETURN
      END
      SUBROUTINE DERIVT (NPTS, NPDE, UNP1, UN, UNM1, DT, DTRAT, UT)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION UNP1(NPTS,NPDE), UN(NPTS,NPDE), UNM1(NPTS,NPDE),
     +   DT, DTRAT,
     +   UT(NPTS,NPDE)
C
Ccc PURPOSE:
C Compute time derivative. If DTRAT = 0 first order results,
C if DTRAT = DT/DT_old second order.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points of current grid
C NPDE   : IN. # PDE components
C UNP1   : IN. Solution at Tn+1 on current grid
C UN     : IN. Solution at Tn   on current grid
C UNM1   : IN. Solution at Tn-1 on current grid
C DT     : IN. Current time stepsize
C DTRAT  : IN. 0 or DT/DT_old
C UT     : OUT. Time derivative of U on current grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER I, IC
      DOUBLE PRECISION A0, A1, A2

      A0 = (1+2*DTRAT)   / ((1+DTRAT)*DT)
      A1 = -(1+DTRAT)**2 / ((1+DTRAT)*DT)
      A2 = DTRAT**2      / ((1+DTRAT)*DT)
      DO 10 IC = 1, NPDE
      DO 10 I = 1, NPTS
         UT(I,IC) = A0*UNP1(I,IC) + A1*UN(I,IC) + A2*UNM1(I,IC)
   10 CONTINUE
      
      RETURN
      END
      SUBROUTINE DERIVS (NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, LBLWY, LABVY, LBLWZ, LABVZ,
     +   DX, DY, DZ, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*),
     +   LBLWY(0:NPTS), LABVY(0:NPTS), LBLWZ(0:NPTS), LABVZ(0:NPTS)
      DOUBLE PRECISION U(0:NPTS*NPDE), DX, DY, DZ,
     +   UX(NPTS*NPDE), UY(NPTS*NPDE), UZ(NPTS*NPDE),
     +   UXX(NPTS*NPDE), UYY(NPTS*NPDE), UZZ(NPTS*NPDE),
     +   UXY(NPTS*NPDE), UXZ(NPTS*NPDE), UYZ(NPTS*NPDE)
C
Ccc PURPOSE:
C Compute space derivatives with second order approximation. Second
C order derivatives are required only in the interior of the domain.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points of current grid
C NPDE   : IN. # PDE components
C U      : IN. Solution on current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C UX     : OUT. -I
C UY     : OUT.  I
C UZ     : OUT.  I
C UXX    : OUT.  I Space derivatives of U on current grid
C UYY    : OUT.  I
C UZZ    : OUT.  I
C UXY    : OUT.  I
C UXZ    : OUT.  I
C UYZ    : OUT. -I
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER I, IB, IC, IPT, LB, IM1, IM2, IP1, IP2, I1, I2, I3, I4
      DOUBLE PRECISION FACX, FACY, FACZ, FACXX, FACYY, FACZZ, FACXY,
     +   FACXZ, FACYZ
C
Ccc Zero derivative arrays in first and last point
C (will possibly not be initialized)
CDIR$ NEXTSCALAR
      DO 10 IC = 1, NPDE
         IPT = 1
         I = IPT + (IC-1)*NPTS
         UX (I) = 0.0
         UY (I) = 0.0
         UZ (I) = 0.0
         UXX(I) = 0.0
         UYY(I) = 0.0
         UZZ(I) = 0.0
         UXY(I) = 0.0
         UXZ(I) = 0.0
         UYZ(I) = 0.0
         IPT = NPTS
         I = IPT + (IC-1)*NPTS
         UX (I) = 0.0
         UY (I) = 0.0
         UZ (I) = 0.0
         UXX(I) = 0.0
         UYY(I) = 0.0
         UZZ(I) = 0.0
         UXY(I) = 0.0
         UXZ(I) = 0.0
         UYZ(I) = 0.0
   10 CONTINUE
C
Ccc Compute derivatives in interior points, boundary values will be
C rubbish
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACZ  = 1/(2*DZ)
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
      FACXY = 1/(2*DX*2*DY)
      FACXZ = 1/(2*DX*2*DZ)
      FACYZ = 1/(2*DY*2*DZ)
C
      DO 20 I = 2, NPTS*NPDE-1
         UX (I) = (U(I+1)-U(I-1))*FACX
         UXX(I) = (U(I+1)-2*U(I)+U(I-1))*FACXX
   20 CONTINUE
C
      DO 30 IC = 1, NPDE
      DO 30 IPT = 2, NPTS-1
         IM1 = LBLWY(IPT) + (IC-1)*NPTS
         I   =       IPT  + (IC-1)*NPTS
         IP1 = LABVY(IPT) + (IC-1)*NPTS
         UY (I) = (U(IP1)-U(IM1))*FACY
         UYY(I) = (U(IP1)-2*U(I)+U(IM1))*FACYY
         IM1 = LBLWZ(IPT) + (IC-1)*NPTS
         I   =       IPT  + (IC-1)*NPTS
         IP1 = LABVZ(IPT) + (IC-1)*NPTS
         UZ (I) = (U(IP1)-U(IM1))*FACZ
         UZZ(I) = (U(IP1)-2*U(I)+U(IM1))*FACZZ
C
         I1 = LABVY(IPT-1) + (IC-1)*NPTS
         I2 = LABVY(IPT+1) + (IC-1)*NPTS
         I3 = LBLWY(IPT-1) + (IC-1)*NPTS
         I4 = LBLWY(IPT+1) + (IC-1)*NPTS
         UXY(I) = (U(I2)-U(I1)-U(I4)+U(I3))*FACXY
         I1 = LABVZ(IPT-1) + (IC-1)*NPTS
         I2 = LABVZ(IPT+1) + (IC-1)*NPTS
         I3 = LBLWZ(IPT-1) + (IC-1)*NPTS
         I4 = LBLWZ(IPT+1) + (IC-1)*NPTS
         UXZ(I) = (U(I2)-U(I1)-U(I4)+U(I3))*FACXZ
         I1 = LABVZ(LBLWY(IPT)) + (IC-1)*NPTS
         I2 = LABVZ(LABVY(IPT)) + (IC-1)*NPTS
         I3 = LBLWZ(LBLWY(IPT)) + (IC-1)*NPTS
         I4 = LBLWZ(LABVY(IPT)) + (IC-1)*NPTS
         UYZ(I) = (U(I2)-U(I1)-U(I4)+U(I3))*FACYZ
   30 CONTINUE
C
Ccc Correct physical boundaries
      DO 40 IB = 1, LLBND(0)
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane, correct Ux
            DO 50 IC = 1, NPDE
            DO 50 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   = IPT + (IC-1)*NPTS
               UX(I) = (-3*U(I)+4*U(I+1)-U(I+2))*FACX
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down plane, correct Uz
            DO 60 IC = 1, NPDE
            DO 60 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   =       IPT  + (IC-1)*NPTS
               IP1 = LABVZ(IPT)
               IP2 = LABVZ(IP1) + (IC-1)*NPTS
               IP1 =       IP1  + (IC-1)*NPTS
               UZ(I) = (-3*U(I)+4*U(IP1)-U(IP2))*FACZ
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane, correct Ux
            DO 70 IC = 1, NPDE
            DO 70 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   = IPT + (IC-1)*NPTS
               UX(I) = (+3*U(I)-4*U(I-1)+U(I-2))*FACX
   70       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane, correct Uz
            DO 80 IC = 1, NPDE
            DO 80 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   =       IPT  + (IC-1)*NPTS
               IM1 = LBLWZ(IPT)
               IM2 = LBLWZ(IM1) + (IC-1)*NPTS
               IM1 =       IM1  + (IC-1)*NPTS
               UZ(I) = (+3*U(I)-4*U(IM1)+U(IM2))*FACZ
   80       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front plane, correct Uy
            DO 90 IC = 1, NPDE
            DO 90 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   =       IPT  + (IC-1)*NPTS
               IP1 = LABVY(IPT)
               IP2 = LABVY(IP1) + (IC-1)*NPTS
               IP1 =       IP1  + (IC-1)*NPTS
               UY(I) = (-3*U(I)+4*U(IP1)-U(IP2))*FACY
   90       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back boundary, correct Uy
            DO 100 IC = 1, NPDE
            DO 100 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               I   =       IPT  + (IC-1)*NPTS
               IM1 = LBLWY(IPT)
               IM2 = LBLWY(IM1) + (IC-1)*NPTS
               IM1 =       IM1  + (IC-1)*NPTS
               UY(I) = (+3*U(I)-4*U(IM1)+U(IM2))*FACY
  100       CONTINUE
         ENDIF
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, F)
C
C-----------------------------------------------------------------------
C
C PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION T, X(*), Y(*), Z(*), U(NPTS*NPDE), UIB(*),
     +   UT(*), UX(*), UY(*), UZ(*),
     +   UXX(*), UYY(*), UZZ(*), UXY(*), UXZ(*), UYZ(*), F(NPTS*NPDE)
C
Ccc PURPOSE:
C Compute residual F(t,U,Ut)
C
C PARAMETER DESCRIPTION:
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C U      : IN. Solution at T on current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C F      : OUT. Residual
C
Ccc EXTERNALS USED:
      EXTERNAL PDEBC, PDEF
C
C-----------------------------------------------------------------------
C
      INTEGER I, IB, IC, LB, NBNDS, NIBPTS, IBS, IBE
C
Ccc Get residual on internal domain
      CALL PDEF (T, X, Y, Z, U,
     +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, F, NPTS, NPDE)
C
Ccc Correct residual on physical boundaries
      CALL PDEBC (T, X, Y, Z, U, UT, UX, UY, UZ, F, NPTS, NPDE,
     +   LLBND, ILBND, LBND)
C
Ccc Correct residual on internal boundaries
      NBNDS  = LLBND(0)
      IBS    = LLBND(NBNDS+1)
      IBE    = LLBND(NBNDS+2)-1
      NIBPTS = IBE-IBS+1
      DO 10 IC = 1, NPDE
      DO 10 LB = IBS, IBE
         I  = LBND(LB) + (IC-1)*NPTS
         IB = LB-IBS+1 + (IC-1)*NIBPTS
         F(I) = U(I) - UIB(IB)
   10 CONTINUE

      RETURN
      END
      LOGICAL FUNCTION CHKTIM (RWK, LU, LUO, NPDE, IWK, LSG,
     +    TIMWGT, RELTOL, ABSTOL, WORK, DT, DTNEW, TIMON)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LU, LUO, NPDE, IWK(*), LSG(0:*)
      DOUBLE PRECISION RWK(*), TIMWGT(NPDE), RELTOL(NPDE), ABSTOL(NPDE),
     +   WORK(*), DT, DTNEW, TIMON
C
Ccc PURPOSE:
C Check if time step was OK. If so, CHKTIM = .TRUE. and DTNEW is set
C to the stepsize for the next time step. If not CHKTIM = .FALSE. and
C DTNEW is the stepsize for the next try.
C
Ccc PARAMETER DESCRIPTION:
C RWK    : IN. Work array containing both U and U_old on all grids of
C          this time level
C LU     : IN. Pointer after last element of U on base grid
C LUO    : IN. Pointer to first element of U_old on base grid
C NPDE   : IN. Number of PDE components
C IWK    : IN. Work array containing the datastructures for the
C          different grids on this level
C LSG    : IN. (0:LSG(0))
C          LSG(0): # grid levels for this time step
C          LSG(I): pointer in IWK to datastructure for grid of level I
C TIMWGT : IN. User defined time weight for each PDE component
C          used in check if time stepsize can be accepted
C RELTOL : IN.  (NPDE)
C          Relative time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C ABSTOL : IN.  (NPDE)
C          Absolute time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C WORK   : WORK. (max. # grid points on a level)
C DT     : IN. Current time stepsize
C DTNEW  : OUT. Stepsize for, new or retry of, timestep
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION TIMMON
      EXTERNAL TIMMON
C
C-----------------------------------------------------------------------
C
      INTEGER LEVEL, LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW,
     +   LICOL, LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP, LENU
      DOUBLE PRECISION TIMONL
C
      TIMON = 0.0
      DO 10 LEVEL = 1, LSG(0)
         LLPLN  = LSG(LEVEL)
         NPLNS  = IWK(LLPLN)
         LIPLN  = LLPLN+NPLNS+2
         LLROW  = LIPLN+NPLNS
         NROWS  = IWK(LLPLN+NPLNS+1)-1
         NPTS   = IWK(LLROW+NROWS)-1
         LIROW  = LLROW+NROWS+1
         LICOL  = LIROW+NROWS
         LLLBND = LICOL+NPTS
         NBNDS  = IWK(LLLBND)
         NBDPTS = IWK(LLLBND+NBNDS+1)-1
         NBIPTS = IWK(LLLBND+NBNDS+2)-1
         LILBND = LLLBND+NBNDS+3
         LLBNDP = LILBND+NBNDS
         LENU   = NPTS*NPDE+1
         LU     = LU - LENU
C Compute time monitor for this grid level
         TIMONL = TIMMON (RWK(LU+1), RWK(LUO+1), NPTS, NPDE,
     +      IWK(LLBNDP), NBIPTS, TIMWGT, RELTOL, ABSTOL, WORK)
         LUO    = LUO + LENU*2
C Compute maximum of monitor values for all levels
         TIMON  = MAX(TIMON,TIMONL)
   10 CONTINUE
C
Ccc Compute new stepsize and check if current step can be accepted
C Compute new stepsize such that prediction of next time monitor is 0.5
      IF (TIMON .GT. 1.0) THEN
C Reject step
         CHKTIM = .FALSE.
         DTNEW = 0.5 / TIMON * DT
C       restrict time step variance
         DTNEW = MAX(DTNEW, DT/4)
      ELSE
C Accept step
         CHKTIM = .TRUE.
C       restrict time step variance
         DTNEW = 2*DT
         IF (TIMON .GT. 0.25) DTNEW = 0.5 / TIMON * DT
      ENDIF

      RETURN
      END
      DOUBLE PRECISION FUNCTION TIMMON (U, UO, NPTS, NPDE, LBND, NBIPTS,
     +   TIMWGT,
     +   RELTOL, ABSTOL, DTUT)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NBIPTS
      INTEGER NPTS, NPDE, LBND(NBIPTS)
      DOUBLE PRECISION U(NPTS,NPDE), UO(NPTS,NPDE), TIMWGT(NPDE),
     +   RELTOL(NPDE), ABSTOL(NPDE), DTUT(NPTS)
C
Ccc PURPOSE:
C Compute time monitor for a specific grid level.
C
C Time monitor:
C TIMMON =
C    sqrt{ sum TIMWGT(ic)/N sum [dt.Ut(ipt,ic) / w(ipt,ic)] ** 2 }
C      (ic=1,NPDE)     (ipt=1,NPTS)
C with N = NPTS*NPDE and w(ipt,ic) = ABSTOL(ic) + RELTOL(ic).|U(ipt,ic)|
C On the boundaries Ut is set to zero.
C
Ccc PARAMETER DESCRIPTION:
C U      : IN. Array of solution values at Tn+1 on current grid
C UO     : IN. Array of solution values at Tn on current grid
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C LBND   : IN. Array containing pointers to boundary points in the grid
C NBIPTS : IN. Total # boundary points
C TIMWGT : IN. User defined time weight for each PDE component
C          used in check if time stepsize can be accepted
C RELTOL : IN.  (NPDE)
C          Relative time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C ABSTOL : IN.  (NPDE)
C          Absolute time tolerance used to determine if time stepsize
C          can be accepted and to determine the new step size
C DTUT   : WORK. (NPTS)
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, IPT, LB, N
      DOUBLE PRECISION TMIC, W2
C
      N = NPTS*NPDE
C
      TIMMON = 0.0
      DO 10 IC = 1, NPDE
         IF (TIMWGT(IC) .EQ. 0.0) GOTO 10
         DO 20 IPT = 1, NPTS
            DTUT(IPT) = U(IPT,IC)-UO(IPT,IC)
   20    CONTINUE
         DO 30 LB = 1, NBIPTS
            IPT = LBND(LB)
            DTUT(IPT) = 0.0
   30    CONTINUE
         TMIC = 0.0
         DO 40 IPT = 1, NPTS
            W2 = ABSTOL(IC) + RELTOL(IC)*ABS(U(IPT,IC))
            TMIC = TMIC + (DTUT(IPT) / W2) ** 2
   40    CONTINUE
         TIMMON = TIMMON + TIMWGT(IC)*TMIC/N
   10 CONTINUE
      TIMMON = SQRT(TIMMON)
      
      RETURN
      END
      SUBROUTINE INTGRB (ISTRUC, X, Y, Z, NPDE, UIB, UNP1, UN, UNM1,
     +   RELTOL, ABSTOL, TN, DT, DTRAT, DX, DY, DZ, WT, F, CORR, RWORK,
     +   IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER ISTRUC(0:*), NPDE, IERR
      DOUBLE PRECISION X(*), Y(*), Z(*), UIB(*), UNP1(0:*), UN(0:*),
     +   UNM1(0:*),
     +   RELTOL(NPDE), ABSTOL(NPDE),
     +   TN, DT, DTRAT, DX, DY, DZ,
     +   WT(*), F(*), CORR(*), RWORK(*)
C
Ccc PURPOSE:
C Integration in time with BDF2 (first timestep BE).
C Solve nonlinear system F(Tn+1, Un+1, Udot) = 0 with modified Newton.
C Solve linear systems with ILU-preconditioned BiCGStab.
C
Ccc PARAMETER DESCRIPTION:
C ISTRUC : IN. Data structure Un+1 grid.
C X,Y,Z  : IN. Physical coordinates grid.
C NPDE   : IN. # PDE components
C UIB    : IN. Dirichlet boundary values on internal boundary.
C UNP1   : INOUT. On entry: Initial solution, on exit final solution
C          Newton converged
C UN     : IN. Solution at Tn on Un+1 grid
C UNM1   : IN. Solution at Tn-1 on Un+1 grid
C RELTOL : IN. Relative tolerance for Newton process
C ABSTOL : IN. Absolute tolerance for Newton process
C TN     : IN. Previous time
C DT     : IN. Current time step
C DX     : IN. Current grid spacing in X-direction
C DY     : IN. Current grid spacing in Y-direction
C DZ     : IN. Current grid spacing in Z-direction
C DTRAT  : IN. If BE: 0, if BDF2: DT/DT_old
C WT     : WORK. (NPTS*NPDE)
C          Weight function for norm computation
C F      : WORK. (NPTS*NPDE)
C          Residual
C CORR   : WORK. (NPTS*NPDE)
C          Correction in Newton iteration
C RWORK  : WORK. (JACILU+max(RESWRK,LSSWRK))
C          JACILU: 2.19.NPDE.LENU
C          RESWRK: LENU.10
C          LSSWRK: LENU.5
C          LENU  : NPTS*NPDE
C IERR   : OUT.
C           0: OK.
C          10: Newton process did not converge
C
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION MAXNRM, WDNRM2
      EXTERNAL ERRWGT, BICGST, JAC, JACPB, MAXNRM, RESID, WDNRM2
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
Ccc   INCLUDE 'PARBICGSTAB'
C
C PARBICGSTAB
C
C Parameters for linear system solver BiCGStab
C MAXLIT : Max. number of BiCGStab iterations
C TOLLSB : Tolerance for linear system solver
      INTEGER MAXLIT
      DOUBLE PRECISION TOLLSB
      PARAMETER (MAXLIT = 100, TOLLSB = TOLNEW/10)
C
C end INCLUDE 'PARBICGSTAB'
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
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
C
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ,
     +   LLLDG, LLUDG, LLSLP, LLLSL, LLSUP, LLLSU, LENU, LENGLU,
     +   LUT, LUX, LUY, LUZ, LUXX, LUYY, LUZZ, LUXY, LUXZ, LUYZ,
     +   LBCG1, LBCG2, LBCG3, LBCG4, LBCG5, LG, LGLU, LJACWK,
     +   NJAC, NRES, I, NIT, ITER
      LOGICAL NEWJAC
      DOUBLE PRECISION ERR, CORNRM, OLDNRM, RATE, TOL, UNRM
C
      IERR = 0
C
      IF (LUNNLS .NE. 0) THEN
         WRITE(LUNNLS,'(''Nonlinear system solver at T ='',E16.7)')
     +      TN+DT
      ENDIF
C
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBDPTS = ISTRUC(LLLBND+NBNDS+1)-1
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
C
      LLLDG  = LLABVZ+NPTS+1
      LLUDG  = LLLDG+NPTS*8
      LLSLP  = LLUDG+NPTS*8
      LLLSL  = LLSLP+NPTS
      LLSUP  = LLLSL+ISTRUC(LLLSL)+1
      LLLSU  = LLSUP+NPTS
C
      LENU   = NPTS*NPDE
      LENGLU = LENU*NPDE*19
C
      LUT    = 1
      LUX    = LUT  + LENU
      LUY    = LUX  + LENU
      LUZ    = LUY  + LENU
      LUXX   = LUZ  + LENU
      LUYY   = LUXX + LENU
      LUZZ   = LUYY + LENU
      LUXY   = LUZZ + LENU
      LUXZ   = LUXY + LENU
      LUYZ   = LUXZ + LENU
C
      LBCG1  = LUX
      LBCG2  = LBCG1 + LENU
      LBCG3  = LBCG2 + LENU
      LBCG4  = LBCG3 + LENU
      LBCG5  = LBCG4 + LENU
C
      LG     = MAX (LUYZ+LENU, LBCG5+LENU)
      LGLU   = LG+LENGLU
      LJACWK = LGLU
C
Ccc Set error weights for use in Newton process
      CALL ERRWGT (NPTS, NPDE, UNP1(1), RELTOL, ABSTOL, WT)
C
Ccc Compute weighted norm of initial solution for convergence check
      UNRM = WDNRM2 (LENU, UNP1(1), WT)
C
Ccc Compute derivatives and residual
      CALL RESID (TN+DT, X, Y, Z, NPTS, NPDE, UNP1, UN, UNM1, DT, DTRAT,
     +   UIB, ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY), ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +   DX, DY, DZ,
     +   RWORK(LUT),  RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ), F)
      NRES = 1
      IF (LUNNLS .NE. 0) THEN
         WRITE(LUNNLS,'('' Max. and WRMS norm residual='',2E16.7)')
     +      MAXNRM(LENU, F), WDNRM2 (LENU, F, WT)
      ENDIF
C
Ccc Compute Jacobian G = dF/dU and its incomplete factorization GLU
      CALL JAC (NPTS, NPDE, F, TN+DT, X, Y, Z, DT, DTRAT, DX, DY, DZ,
     +   UNP1(1), ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP), UIB,
     +   RWORK(LUT), RWORK(LUX), RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ),
     +   ABSTOL, RWORK(LJACWK), RWORK(LG))
C Copy Jacobian for factorization
      CALL RCOPY (LENGLU, RWORK(LG), RWORK(LGLU))
C Compute ILU
      CALL JACPB (NPTS, NPDE, RWORK(LGLU), ISTRUC(LLLDG),
     +   ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +   ISTRUC(LLSLP), ISTRUC(LLLSL))
      NEWJAC = .TRUE.
      NJAC = 1
C
Ccc Newton iteration loop
    9 CONTINUE
      DO 10 NIT = 1, MAXNIT
C
Cccccc Solve G.corr = F. Store the residual in F.
         TOL = TOLLSB / (2**NIT)
         CALL BICGST (NPTS, NPDE, RWORK(LG), CORR, F, WT, TOL,
     +      MAXLIT, RWORK(LGLU), ISTRUC(LLLDG), ISTRUC(LLUDG),
     +      ISTRUC(LLSLP), ISTRUC(LLLSL), ISTRUC(LLSUP), ISTRUC(LLLSU),
     +      LUNLSS, RWORK(LBCG1), RWORK(LBCG2), RWORK(LBCG3),
     +      RWORK(LBCG4), RWORK(LBCG5), ITER, ERR, IERR)
         NLSIT(LEVEL,NIT) = NLSIT(LEVEL,NIT)+ ITER
         IF (IERR .NE. 0) GOTO 100
C
Cccccc Test for convergence
         CORNRM = WDNRM2 (LENU, CORR, WT)
         IF (LUNNLS .NE. 0) THEN
            WRITE(LUNNLS,'('' NI:'',I3,'', NLI:'',I4,'', ERLI'':,E16.7,
     +         '', ERNI:'',E16.7)') NIT, ITER, ERR, CORNRM
         ENDIF
         IF (CORNRM .LE. 100*UROUND*UNRM) GOTO 900
         IF (.NOT. NEWJAC) THEN
            RATE = SQRT(CORNRM/OLDNRM)
            IF (RATE .GT. 0.9) THEN
C          Divergence
               GOTO 100
            ELSE IF (RATE/(1-RATE)*CORNRM .LE. TOLNEW) THEN
C          Convergence
               GOTO 900
            ENDIF
         ENDIF
         OLDNRM = CORNRM
C
Ccccc Update solution
         DO 20 I = 1, LENU
            UNP1(I) = UNP1(I) - CORR(I)
   20    CONTINUE
C
Ccc Compute derivatives and residual and start next iteration
         CALL RESID (TN+DT, X, Y, Z, NPTS, NPDE, UNP1, UN, UNM1,
     +      DT, DTRAT, UIB,
     +      ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +      ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +      ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +      DX, DY, DZ,
     +      RWORK(LUT),  RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +      RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +      RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ), F)
         NRES   = NRES+1
         NEWJAC = .FALSE.
         IF (LUNNLS .NE. 0) THEN
            WRITE(LUNNLS,'('' Max. and WRMS norm residual='',2E16.7)')
     +         MAXNRM(LENU, F), WDNRM2 (LENU, F, WT)
         ENDIF
   10 CONTINUE
Ccc End Newton iteration loop
C
Ccc No convergence in max. # iterations
C
Ccccc Check if Jacobian is recent
  100 CONTINUE
         IF (.NOT. NEWJAC .AND. NJAC .LT. MAXJAC) THEN
            IF (LUNNLS .NE. 0) THEN
               WRITE(LUNNLS,'('' New Jacobian, NIT='',I4)') NIT
            ENDIF
C       Compute new Jacobian and retry
            CALL DERIVS (NPTS, NPDE, UNP1,
     +         ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +         ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +         ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +         DX, DY, DZ,
     +         RWORK(LUX), RWORK(LUY), RWORK(LUZ),
     +         RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +         RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ))
            CALL JAC (NPTS, NPDE, F, TN+DT, X, Y, Z, DT, DTRAT,
     +         DX, DY, DZ, UNP1(1),
     +         ISTRUC(LLLBND),ISTRUC(LILBND),ISTRUC(LLBNDP), UIB,
     +         RWORK(LUT), RWORK(LUX), RWORK(LUY), RWORK(LUZ),
     +         RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +         RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ),
     +         ABSTOL, RWORK(LJACWK), RWORK(LG))
C       Copy Jacobian for factorization
            CALL RCOPY (LENGLU, RWORK(LG), RWORK(LGLU))
C       Compute ILU
            CALL JACPB (NPTS, NPDE, RWORK(LGLU), ISTRUC(LLLDG),
     +         ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +         ISTRUC(LLSLP), ISTRUC(LLLSL))
C
            NEWJAC = .TRUE.
            NJAC = NJAC + 1
            GOTO 9
         ELSE
C       Newton failure
            IERR = 10
            NNIT(LEVEL)   = NNIT(LEVEL)+NIT
            NRESID(LEVEL) = NRESID(LEVEL)+NRES
            NJACS(LEVEL)  = NJACS(LEVEL)+NJAC
            IF (LUNNLS .NE. 0) THEN
               WRITE(LUNNLS,'(''Newton failure, NIT='',I4)') NIT
            ENDIF
            RETURN
         ENDIF
C
Ccc Nonlinear proces has been solved
  900 CONTINUE
C Update solution
      DO 30 I = 1, LENU
         UNP1(I) = UNP1(I) - CORR(I)
   30 CONTINUE
C
      NNIT(LEVEL)   = NNIT(LEVEL)+NIT
      NRESID(LEVEL) = NRESID(LEVEL)+NRES
      NJACS(LEVEL)  = NJACS(LEVEL)+NJAC
C
      RETURN
      END
      SUBROUTINE JAC (NPTS, NPDE, F, T, X, Y, Z, DT, DTRAT,
     +   DX, DY, DZ, U, LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, ABSTOL, WORK, G)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION F(*), T, X(*), Y(*), Z(*), DT, DTRAT, DX, DY, DZ,
     +   U(*), UIB(*), UT(*), UX(*), UY(*), UZ(*),
     +   UXX(*), UYY(*), UZZ(*), UXY(*), UXZ(*), UYZ(*),
     +   ABSTOL(*), WORK(*), G(*)
C
Ccc PURPOSE:
C Compute Jacobian G = dF/dU and store in block 19-diagonal mode.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C DT     : IN. Current time stepsize
C DTRAT  : IN. 0 or DT/DT_old
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C U      : IN. Solution at T on current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : OUT. Time derivative of U on current grid
C UX     : OUT. -I
C UY     : OUT.  I
C UZ     : OUT.  I
C UXX    : OUT.  I Space derivatives of U on current grid
C UYY    : OUT.  I
C UZZ    : OUT.  I
C UXY    : OUT.  I
C UXZ    : OUT.  I
C UYZ    : OUT. -
C ABSTOL : IN. Absolute tolerance for Newton process
C WORK   : WORK. (10*LENFU+2*LENU+NPTS)
C G      : OUT. Jacobian stored in block 19-diagonal mode
C
Ccc EXTERNALS USED:
      EXTERNAL DERIVF, JACG
C
C-----------------------------------------------------------------------
C
      INTEGER LENU, LENFU, LFU, LFUX, LFUY, LFUZ, LFUXX, LFUYY, LFUZZ,
     +   LFUXY, LFUXZ, LFUYZ, LDEL, LRWK
      DOUBLE PRECISION A0
C
      LENU  = NPTS*NPDE
      LENFU = LENU*NPDE
C
      LFU   = 1
      LFUX  = LFU   + LENFU
      LFUY  = LFUX  + LENFU
      LFUZ  = LFUY  + LENFU
      LFUXX = LFUZ  + LENFU
      LFUYY = LFUXX + LENFU
      LFUZZ = LFUYY + LENFU
      LFUXY = LFUZZ + LENFU
      LFUXZ = LFUXY + LENFU
      LFUYZ = LFUXZ + LENFU
      LDEL  = LFUYZ + LENFU
      LRWK  = LDEL  + NPTS
C
Ccc Compute dF/dU, dF/dUx, dF/dUy, dF/dUz, dF/dUxx, dF/dUyy, dF/dUzz,
C dF/dUxy, dF/dUxz, dF/dUyz 
      A0 = (1+2*DTRAT) / ((1+DTRAT)*DT)
      CALL DERIVF (F, T, X, Y, Z, NPTS, NPDE, U, A0, DT, DX, DY, DZ,
     +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +   UXY, UXZ, UYZ, ABSTOL, WORK(LDEL), WORK(LRWK),
     +   WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +   WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ),
     +   WORK(LFUXY), WORK(LFUXZ), WORK(LFUYZ))
C
Ccc Compute G = dF/dU + dF/dUx.dUx/dU + ...
      CALL JACG (NPTS, NPDE, DX, DY, DZ, LLBND, ILBND, LBND,
     +   WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +   WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ),
     +   WORK(LFUXY), WORK(LFUXZ), WORK(LFUYZ), G)
C
      RETURN
      END
      SUBROUTINE PRTRBU (ICPTB, NPTS, NPDE, U, A0, DT, UT, TOL, DEL,
     +   UBAR, UTBAR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER ICPTB, NPTS, NPDE
      DOUBLE PRECISION U(NPTS,NPDE), A0, DT, UT(NPTS,NPDE), TOL,
     +   DEL(NPTS),
     +   UBAR(NPTS,NPDE), UTBAR(NPTS,NPDE)
C
Ccc PURPOSE:
C Perturb the ICPTB-th component of U. Return perturbance in DEL and
C perturbed U in UBAR.
C
Ccc PARAMETER DESCRIPTION:
C ICPTB  : IN. Component of U to be perturbed
C NPTS   : IN. # gridpoints
C NPDE   : IN. # PDE components
C U      : IN. Solution or derivative of solution to be perturbed
C TOL    : IN. Threshold for perturbation
C DEL    : OUT. Perturbation values
C UBAR   : OUT. Perturbed values of U
C
Ccc EXTERNALS USED:
      EXTERNAL RCOPY
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
      INTEGER IPT
      DOUBLE PRECISION DELI, EPS

      CALL RCOPY (NPTS*NPDE, U, UBAR)
      CALL RCOPY (NPTS*NPDE, UT, UTBAR)

      EPS = SQRT(UROUND)
      DO 10 IPT = 1, NPTS
C Compute perturbance, if U=0, U(T+dt)=dtUt, if both are zero take
C threshold
         DELI = EPS*MAX(ABS(U(IPT,ICPTB)),ABS(DT*UT(IPT,ICPTB)),TOL)
         DELI = SIGN(DELI,DT*UT(IPT,ICPTB))
C To ensure that the perturbance is the same machine number as the
C denominator
         DEL(IPT) = (U(IPT,ICPTB)+DELI)-U(IPT,ICPTB)
         UBAR(IPT,ICPTB) = U(IPT,ICPTB) + DEL(IPT)
         UTBAR(IPT,ICPTB) = UT(IPT,ICPTB) + A0*DEL(IPT)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PERTRB (ICPTB, NPTS, NPDE, U, TOL, DEL, UBAR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER ICPTB, NPTS, NPDE
      DOUBLE PRECISION U(NPTS,NPDE), TOL, DEL(NPTS), UBAR(NPTS,NPDE)
C
Ccc PURPOSE:
C Perturb the ICPTB-th component of U. Return perturbance in DEL and
C perturbed U in UBAR.
C
Ccc PARAMETER DESCRIPTION:
C ICPTB  : IN. Component of U to be perturbed
C NPTS   : IN. # gridpoints
C NPDE   : IN. # PDE components
C U      : IN. Solution or derivative of solution to be perturbed
C TOL    : IN. Threshold for perturbation
C DEL    : OUT. Perturbation values
C UBAR   : OUT. Perturbed values of U
C
Ccc EXTERNALS USED:
      EXTERNAL RCOPY
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
      INTEGER IPT
      DOUBLE PRECISION DELI, EPS

      CALL RCOPY (NPTS*NPDE, U, UBAR)

      EPS = SQRT(UROUND)
      DO 10 IPT = 1, NPTS
C Compute perturbance
         DELI = EPS*MAX(ABS(U(IPT,ICPTB)),TOL)
C To ensure that UBAR has the same sign as U
         DELI = SIGN(DELI,U(IPT,ICPTB))
C To ensure that the perturbance is the same machine number as the
C denominator
         DEL(IPT) = (U(IPT,ICPTB)+DELI)-U(IPT,ICPTB)
         UBAR(IPT,ICPTB) = U(IPT,ICPTB) + DEL(IPT)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE JACG (NPTS, NPDE, DX, DY, DZ, LLBND, ILBND, LBND,
     +   FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ, FUXY, FUXZ, FUYZ, G)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION DX, DY, DZ,
     +   FU(NPTS*NPDE,NPDE),
     +   FUX(NPTS*NPDE,NPDE), FUY(NPTS*NPDE,NPDE), FUZ(NPTS*NPDE,NPDE),
     +   FUXX(NPTS*NPDE,NPDE),FUYY(NPTS*NPDE,NPDE),FUZZ(NPTS*NPDE,NPDE),
     +   FUXY(NPTS*NPDE,NPDE),FUXZ(NPTS*NPDE,NPDE),FUYZ(NPTS*NPDE,NPDE),
     +   G(NPTS*NPDE,NPDE,-9:9)
C
Ccc PURPOSE:
C Compute Jacobian G = dF/dU using derivatives of residual wrt
C (derivatives of) U
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FU     : IN. Derivative residual F(.,U,Ut,.) wrt U
C FUX    : IN. Derivative residual F(.,Ux,.)   wrt Ux
C FUY    : IN. Derivative residual F(.,Uy,.)   wrt Uy
C FUZ    : IN. Derivative residual F(.,Uz,.)   wrt Uz
C FUXX   : IN. Derivative residual F(.,Uxx,.)  wrt Uxx
C FUYY   : IN. Derivative residual F(.,Uyy,.)  wrt Uyy
C FUZZ   : IN. Derivative residual F(.,Uzz,.)  wrt Uzz
C FUXY   : IN. Derivative residual F(.,Uxy,.)  wrt Uxy
C FUXZ   : IN. Derivative residual F(.,Uxz,.)  wrt Uxz
C FUYZ   : IN. Derivative residual F(.,Uyz,.)  wrt Uyz
C G      : OUT. Jacobian stored in block 19-diagonal mode
C
Ccc EXTERNALS USED:
      EXTERNAL JACGBD
C
C-----------------------------------------------------------------------
C
      INTEGER I, JC, LENU
      DOUBLE PRECISION FACX, FACY, FACZ, FACXX, FACYY, FACZZ, FACXY,
     +   FACXZ, FACYZ
C
      LENU  = NPTS*NPDE
C
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACZ  = 1/(2*DZ)
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
      FACXY = 1/(2*DX*2*DY)
      FACXZ = 1/(2*DX*2*DZ)
      FACYZ = 1/(2*DY*2*DZ)
C
Ccc First internal domain
      DO 10 JC = 1, NPDE
      DO 10 I = 1, LENU
C dF(ipt,ic)/dU(blwz(blwy(ipt)),jc)
         G(I,JC,-9) = FUYZ(I,JC)*(+FACYZ)
C dF(ipt,ic)/dU(blwz(ipt-1),jc)
         G(I,JC,-8) = FUXZ(I,JC)*(+FACXZ)
C dF(ipt,ic)/dU(blwz(ipt),jc)
         G(I,JC,-7) =
     +      FUZ(I,JC)*(-FACZ) + FUZZ(I,JC)*(+FACZZ)
C dF(ipt,ic)/dU(blwz(ipt+1),jc)
         G(I,JC,-6) = FUXZ(I,JC)*(-FACXZ)
C dF(ipt,ic)/dU(blwz(abvy(ipt)),jc)
         G(I,JC,-5) = FUYZ(I,JC)*(-FACYZ)
C dF(ipt,ic)/dU(blwy(ipt)-1,jc)
         G(I,JC,-4) = FUXY(I,JC)*(+FACXY)
C dF(ipt,ic)/dU(blwy(ipt),jc)
         G(I,JC,-3) =
     +      FUY(I,JC)*(-FACY) + FUYY(I,JC)*(+FACYY)
C dF(ipt,ic)/dU(blwy(ipt)+1,jc)
         G(I,JC,-2) = FUXY(I,JC)*(-FACXY)
C dF(ipt,ic)/dU(ipt-1,jc)
         G(I,JC,-1) =
     +      FUX(I,JC)*(-FACX) + FUXX(I,JC)*(+FACXX)
C dF(ipt,ic)/dU(ipt,jc)
         G(I,JC, 0) = FU(I,JC) + 
     +      FUXX(I,JC)*(-2*FACXX) + FUYY(I,JC)*(-2*FACYY) +
     +      FUZZ(I,JC)*(-2*FACZZ)
C dF(ipt,ic)/dU(ipt+1,jc)
         G(I,JC,+1) =
     +      FUX(I,JC)*(+FACX) + FUXX(I,JC)*(+FACXX)
C dF(ipt,ic)/dU(abvy(ipt)-1,jc)
         G(I,JC,+2) = FUXY(I,JC)*(-FACXY)
C dF(ipt,ic)/dU(abvy(ipt),jc)
         G(I,JC,+3) =
     +      FUY(I,JC)*(+FACY) + FUYY(I,JC)*(+FACYY)
C dF(ipt,ic)/dU(abvy(ipt)+1,jc)
         G(I,JC,+4) = FUXY(I,JC)*(+FACXY)
C dF(ipt,ic)/dU(abvz(blwy(ipt)),jc)
         G(I,JC,+5) = FUYZ(I,JC)*(-FACYZ)
C dF(ipt,ic)/dU(abvz(ipt-1),jc)
         G(I,JC,+6) = FUXZ(I,JC)*(-FACXZ)
C dF(ipt,ic)/dU(abvz(ipt),jc)
         G(I,JC,+7) =
     +      FUZ(I,JC)*(+FACZ) + FUZZ(I,JC)*(+FACZZ)
C dF(ipt,ic)/dU(abvz(ipt+1),jc)
         G(I,JC,+8) = FUXZ(I,JC)*(+FACXZ)
C dF(ipt,ic)/dU(abvz(abvy(ipt)),jc)
         G(I,JC,+9) = FUYZ(I,JC)*(+FACYZ)
   10 CONTINUE
C
C Correct boundaries
      CALL JACGBD (NPTS, NPDE, FACX, FACY, FACZ,
     +   LLBND, ILBND, LBND, FUX, FUY, FUZ, G)
C
      RETURN
      END
      SUBROUTINE JACGBD (NPTS, NPDE, FACX, FACY, FACZ,
     +   LLBND, ILBND, LBND, FUX, FUY, FUZ, G)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION FACX, FACY, FACZ,
     +   FUX(NPTS,NPDE,NPDE), FUY(NPTS,NPDE,NPDE), FUZ(NPTS,NPDE,NPDE),
     +   G(NPTS,NPDE,NPDE,-9:9)
C
Ccc PURPOSE:
C Correct Jacobian G = dF/dU for second order approximation of
C first order derivatives at boundaries
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C FACX   : IN. 1/(2*DX)
C FACY   : IN. 1/(2*DY)
C FACZ   : IN. 1/(2*DZ)
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FUX    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Ux
C FUY    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uy
C FUZ    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uz
C G      : INOUT.
C          IN: Jacobian stored in block 19-diagonal mode
C          OUT: Jacobian corrected for first order derivatives at
C               boundaries
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, IC, JC, IB, LB
C
Ccc Boundary corrections, no corrections needed for internal boundaries
      DO 10 IB = 1, LLBND(0)
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane, correction needed for dF/dUx
            DO 20 JC = 1, NPDE
            DO 20 IC = 1, NPDE
CDIR$ IVDEP
            DO 25 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC,-1) = 0.0
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUX(IPT,IC,JC)*(-3*FACX)
               G(IPT,IC,JC,+1) = FUX(IPT,IC,JC)*(+4*FACX)
C             dF(ipt,ic)/dU(ipt+2),jc)
               G(IPT,IC,JC,+2) = FUX(IPT,IC,JC)*(-FACX)
   25       CONTINUE
   20       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down plane, correction needed for dF/dUz
            DO 30 JC = 1, NPDE
            DO 30 IC = 1, NPDE
CDIR$ IVDEP
            DO 35 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC,-7) = 0.0
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUZ(IPT,IC,JC)*(-3*FACZ)
               G(IPT,IC,JC,+7) = FUZ(IPT,IC,JC)*(+4*FACZ)
C             dF(ipt,ic)/dU(above(above(ipt)),jc)
               G(IPT,IC,JC,+9) = FUZ(IPT,IC,JC)*(-FACZ)
   35       CONTINUE
   30       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane, correction needed for dF/dUx
            DO 40 JC = 1, NPDE
            DO 40 IC = 1, NPDE
CDIR$ IVDEP
            DO 45 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
C             dF(ipt,ic)/dU(ipt-2),jc)
               G(IPT,IC,JC,-2) = FUX(IPT,IC,JC)*(+FACX)
               G(IPT,IC,JC,-1) = FUX(IPT,IC,JC)*(-4*FACX)
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUX(IPT,IC,JC)*(+3*FACX)
               G(IPT,IC,JC,+1) = 0.0
   45       CONTINUE
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane, correction needed for dF/dUz
            DO 50 JC = 1, NPDE
            DO 50 IC = 1, NPDE
CDIR$ IVDEP
            DO 55 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
C             dF(ipt,ic)/dU(below(below(ipt)),jc)
               G(IPT,IC,JC,-9) = FUZ(IPT,IC,JC)*(+FACZ)
               G(IPT,IC,JC,-7) = FUZ(IPT,IC,JC)*(-4*FACZ)
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUZ(IPT,IC,JC)*(+3*FACZ)
               G(IPT,IC,JC,+7) = 0.0
   55       CONTINUE
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front plane, correction needed for dF/dUy
            DO 60 JC = 1, NPDE
            DO 60 IC = 1, NPDE
CDIR$ IVDEP
            DO 65 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC,-3) = 0.0
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUY(IPT,IC,JC)*(-3*FACY)
               G(IPT,IC,JC,+3) = FUY(IPT,IC,JC)*(+4*FACY)
C             dF(ipt,ic)/dU(above(above(ipt)),jc)
               G(IPT,IC,JC,+4) = FUY(IPT,IC,JC)*(-FACY)
   65       CONTINUE
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back plane, correction needed for dF/dUy
            DO 70 JC = 1, NPDE
            DO 70 IC = 1, NPDE
CDIR$ IVDEP
            DO 75 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
C             dF(ipt,ic)/dU(below(below(ipt)),jc)
               G(IPT,IC,JC,-4) = FUY(IPT,IC,JC)*(+FACY)
               G(IPT,IC,JC,-3) = FUY(IPT,IC,JC)*(-4*FACY)
               G(IPT,IC,JC, 0) = G(IPT,IC,JC,0) +
     +            FUY(IPT,IC,JC)*(+3*FACY)
               G(IPT,IC,JC,+3) = 0.0
   75       CONTINUE
   70       CONTINUE
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE JACSDP (NPTS, LLBND, ILBND, LBND,
     +   LBLWY, LABVY, LBLWZ, LABVZ, LLDG, LUDG)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, LLBND(0:*), ILBND(*), LBND(*),
     +   LBLWY(0:NPTS), LABVY(0:NPTS), LBLWZ(0:NPTS), LABVZ(0:NPTS),
     +   LLDG(NPTS,-9:-2), LUDG(NPTS,2:9)
C
Ccc PURPOSE:
C Set pointers to nodes of lower 8 subdiagonals of Jacobian in LLDG and
C to nodes of upper 8 superdiagonals in LUDG. All nonexisting diagonals
C should point to the main diagonal nodes.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C LLDG   : OUT. (NPTS,-9:-2)
C          LLDG(IPT,-9): pointer to node Y-below Z-below
C                        or to node Z-below Z-below
C          LLDG(IPT,-8): pointer to node left of Z-below
C          LLDG(IPT,-7): pointer to node Z-below
C          LLDG(IPT,-6): pointer to node right of Z-below
C          LLDG(IPT,-5): pointer to node Y-above Z-below
C          LLDG(IPT,-4): pointer to node left of Y-below
C                        or to node Y-below Y-below
C          LLDG(IPT,-3): pointer to node Y-below
C          LLDG(IPT,-2): pointer to node right of Y-below
C                        or to node left of the node left
C LUDG   : OUT. (NPTS,2:9)
C          LUDG(IPT,2): pointer to node left of Y-above
C                       or to node right of the node right
C          LUDG(IPT,3): pointer to node Y-above
C          LUDG(IPT,4): pointer to node right of node Y-above
C                       or to node Y-above Y-above
C          LUDG(IPT,5): pointer to node Y-below Z-above
C          LUDG(IPT,6): pointer to node left of Z-above
C          LUDG(IPT,7): pointer to node Z-above
C          LUDG(IPT,8): pointer to node right of Z-above
C          LUDG(IPT,9): pointer to node Y-above Z-above
C                       or to node Z-above Z-above
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, NBNDS, IB, LB
C
Ccc First internal domain
      DO 10 IPT = 1, NPTS-1
         LLDG(IPT,-9) = LBLWZ(LBLWY(IPT))
         LLDG(IPT,-8) = LBLWZ(IPT-1)
         LLDG(IPT,-7) = LBLWZ(IPT)
         LLDG(IPT,-6) = LBLWZ(IPT+1)
         LLDG(IPT,-5) = LBLWZ(LABVY(IPT))
         LLDG(IPT,-4) = LBLWY(IPT-1)
         LLDG(IPT,-3) = LBLWY(IPT)
         LLDG(IPT,-2) = LBLWY(IPT+1)
         LUDG(IPT,+2) = LABVY(IPT-1)
         LUDG(IPT,+3) = LABVY(IPT)
         LUDG(IPT,+4) = LABVY(IPT+1)
         LUDG(IPT,+5) = LABVZ(LBLWY(IPT))
         LUDG(IPT,+6) = LABVZ(IPT-1)
         LUDG(IPT,+7) = LABVZ(IPT)
         LUDG(IPT,+8) = LABVZ(IPT+1)
         LUDG(IPT,+9) = LABVZ(LABVY(IPT))
   10 CONTINUE
      IPT = NPTS
         LLDG(IPT,-9) = LBLWZ(LBLWY(IPT))
         LLDG(IPT,-8) = LBLWZ(IPT-1)
         LLDG(IPT,-7) = LBLWZ(IPT)
         LLDG(IPT,-5) = LBLWZ(LABVY(IPT))
         LLDG(IPT,-4) = LBLWY(IPT-1)
         LLDG(IPT,-3) = LBLWY(IPT)
         LUDG(IPT,+2) = LABVY(IPT-1)
         LUDG(IPT,+3) = LABVY(IPT)
         LUDG(IPT,+5) = LABVZ(LBLWY(IPT))
         LUDG(IPT,+6) = LABVZ(IPT-1)
         LUDG(IPT,+7) = LABVZ(IPT)
         LUDG(IPT,+9) = LABVZ(LABVY(IPT))
C
Ccc Correct boundaries
      NBNDS  = LLBND(0)
      DO 20 IB = 1, NBNDS
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane
C          I /
C          O - -
C        / I
            DO 30 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-8) = IPT
               IF (LLDG(IPT,-4) .EQ. LBLWY(IPT-1)) LLDG(IPT,-4) = IPT
               LUDG(IPT,+2) = IPT+2
               LUDG(IPT,+6) = IPT
   30       CONTINUE
        ELSE IF (ILBND(IB) .EQ. 2) THEN
C      Down plane
C            I
C            I /
C          - O -
C          /
            DO 40 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-9) = IPT
               LLDG(IPT,-8) = IPT
               LLDG(IPT,-7) = IPT
               LLDG(IPT,-6) = IPT
               LLDG(IPT,-5) = IPT
               LUDG(IPT,+9) = LABVZ(LABVZ(IPT))
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane
C              I /
C          - - O
C            / I
            DO 50 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-6) = IPT
               LLDG(IPT,-2) = IPT-2
	       IF (IPT .EQ. NPTS) THEN
		  LUDG(IPT,+4) = IPT
	       ELSE IF (LUDG(IPT,+4) .EQ. LABVY(IPT+1)) THEN
		  LUDG(IPT,+4) = IPT
	       ENDIF
               LUDG(IPT,+8) = IPT
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane
C              /
C          - O -
C          / I
C            I
            DO 60 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-9) = LBLWZ(LBLWZ(IPT))
               LUDG(IPT,+5) = IPT
               LUDG(IPT,+6) = IPT
               LUDG(IPT,+7) = IPT
               LUDG(IPT,+8) = IPT
               LUDG(IPT,+9) = IPT
   60       CONTINUE
        ELSE IF (ILBND(IB) .EQ. 5) THEN
C      Front plane
C                /
C            I /
C          - O -
C            I
            DO 70 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               IF (LLDG(IPT,-9) .EQ. LBLWZ(LBLWY(IPT)))
     +            LLDG(IPT,-9) = IPT
               LLDG(IPT,-4) = IPT
               LLDG(IPT,-3) = IPT
               IF (LLDG(IPT,-2) .EQ. LBLWY(IPT+1)) LLDG(IPT,-2) = IPT
               LUDG(IPT,+4) = LABVY(LABVY(IPT))
               LUDG(IPT,+5) = IPT
   70       CONTINUE
        ELSE IF (ILBND(IB) .EQ. 6) THEN
C      Back plane
C            I
C          - O -
C          / I
C        /
            DO 80 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-5) = IPT
               LLDG(IPT,-4) = LBLWY(LBLWY(IPT))
               IF (LUDG(IPT,+2) .EQ. LABVY(IPT-1)) LUDG(IPT,+2) = IPT
               LUDG(IPT,+3) = IPT
               LUDG(IPT,+4) = IPT
               IF (LUDG(IPT,+9) .EQ. LABVZ(LABVY(IPT)))
     +            LUDG(IPT,+9) = IPT
   80       CONTINUE
        ENDIF
   20 CONTINUE
C
      IB = NBNDS+1
CDIR$ VECTOR
C       Internal boundary, Dirichlet condition, no off diagonals
C          . . .
C          . O .
C          . . .
            DO 200 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLDG(IPT,-9) = IPT
               LLDG(IPT,-8) = IPT
               LLDG(IPT,-7) = IPT
               LLDG(IPT,-6) = IPT
               LLDG(IPT,-5) = IPT
               LLDG(IPT,-4) = IPT
               LLDG(IPT,-3) = IPT
               LLDG(IPT,-2) = IPT
               LUDG(IPT,+2) = IPT
               LUDG(IPT,+3) = IPT
               LUDG(IPT,+4) = IPT
               LUDG(IPT,+5) = IPT
               LUDG(IPT,+6) = IPT
               LUDG(IPT,+7) = IPT
               LUDG(IPT,+8) = IPT
               LUDG(IPT,+9) = IPT
  200       CONTINUE
C
      RETURN
      END
      SUBROUTINE JACPB (NPTS, NPDE, GLU, LLDG,
     +   LLBND, ILBND, LBND, LSL, LLSL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLDG(*),
     +   LLBND(0:*), ILBND(*), LBND(*), LSL(*), LLSL(0:*)
      DOUBLE PRECISION GLU(*)

C
Ccc PURPOSE:
C Compute ILU factorization of the Jacobian in block 19-diagonal mode.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C GLU    : INOUT.
C          IN: Jacobian stored in block 19-diagonal mode
C          OUT: ILU factorization of Jacobian stored in block
C               19-diagonal mode
C LLDG   : IN. (NPTS,-9:-2)
C          LLDG(IPT,-9): pointer to node Y-below Z-below
C                        or to node Z-below Z-below
C          LLDG(IPT,-8): pointer to node left of Z-below
C          LLDG(IPT,-7): pointer to node Z-below
C          LLDG(IPT,-6): pointer to node right of Z-below
C          LLDG(IPT,-5): pointer to node Y-above Z-below
C          LLDG(IPT,-4): pointer to node left of Y-below
C                        or to node Y-below Y-below
C          LLDG(IPT,-3): pointer to node Y-below
C          LLDG(IPT,-2): pointer to node right of Y-below
C                        or to node left of the node left
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LSL    : IN. (NPTS)
C          LSL(ISLPT): pointer to node in actual grid
C LLSL   : IN. (0:LLSL(0))
C          LLSL(0) = # independent data dependency lists in ILU
C                    factorization and forward sweep
C          LLSL(1:LLSL(0)): pointers to the start of a list in LSL
C
Ccc EXTERNALS USED:
      EXTERNAL ILU, JAC19
C
C-----------------------------------------------------------------------
C
Ccc Adapt Jacobian to real block 19-diagonal structure by replacing
C second-order boundary discretization by first-order
      CALL JAC19 (NPTS, NPDE, GLU, LLBND, ILBND, LBND)
C
Ccc Incomplete LU factorization
      CALL ILU (NPTS, NPDE, GLU, LLDG, LSL, LLSL)
C
      RETURN
      END
      SUBROUTINE JAC19 (NPTS, NPDE, A, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION A(NPTS,NPDE,NPDE,-9:9)
C
Ccc PURPOSE:
C Replace second-order boundary discretization by first-order in
C Jacobian to get real block 19-diagonal structure
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C A      : INOUT.
C          IN: Jacobian
C          OUT: Jacobian with second-order boundary discretization
C               replaced by first-order
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, IC, JC, IB, LB
C
      DO 10 IB = 1, LLBND(0)
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane, correction needed for dF/dUx
            DO 20 IC = 1, NPDE
            DO 20 JC = 1, NPDE
CDIR$ IVDEP
            DO 25 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,2)
               A(IPT,IC,JC,2) = 0.0
   25       CONTINUE
   20       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down plane, correction needed for dF/dUz
            DO 30 IC = 1, NPDE
            DO 30 JC = 1, NPDE
CDIR$ IVDEP
            DO 35 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,9)
               A(IPT,IC,JC,9) = 0.0
   35       CONTINUE
   30       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane, correction needed for dF/dUx
            DO 40 IC = 1, NPDE
            DO 40 JC = 1, NPDE
CDIR$ IVDEP
            DO 45 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,-2)
               A(IPT,IC,JC,-2) = 0.0
   45       CONTINUE
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane, correction needed for dF/dUz
            DO 50 IC = 1, NPDE
            DO 50 JC = 1, NPDE
CDIR$ IVDEP
            DO 55 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,-9)
               A(IPT,IC,JC,-9) = 0.0
   55       CONTINUE
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front plane, correction needed for dF/dUy
            DO 60 IC = 1, NPDE
            DO 60 JC = 1, NPDE
            DO 65 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,4)
               A(IPT,IC,JC,4) = 0.0
   65       CONTINUE
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back plane, correction needed for dF/dUy
            DO 70 IC = 1, NPDE
            DO 70 JC = 1, NPDE
            DO 75 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               A(IPT,IC,JC,0) = A(IPT,IC,JC,0) + A(IPT,IC,JC,-4)
               A(IPT,IC,JC,-4) = 0.0
   75       CONTINUE
   70       CONTINUE
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE JACSLP (NPTS, LLBND, ILBND, LBND, LLDG, M, LLS, LS)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, LLBND(0:*), ILBND(*), LBND(*), LLDG(NPTS,-9:-2),
     +   M(NPTS), LLS(0:*), LS(NPTS)
C
Ccc PURPOSE:
C Make data-dependency list for ILU factorization and forward sweep of
C backsolve.
C
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LLDG   : IN. (NPTS,-9:-2)
C          LLDG(IPT,-9): pointer to node Y-below Z-below
C                        or to node Z-below Z-below
C          LLDG(IPT,-8): pointer to node left of Z-below
C          LLDG(IPT,-7): pointer to node Z-below
C          LLDG(IPT,-6): pointer to node right of Z-below
C          LLDG(IPT,-5): pointer to node Y-above Z-below
C          LLDG(IPT,-4): pointer to node left of Y-below
C                        or to node Y-below Y-below
C          LLDG(IPT,-3): pointer to node Y-below
C          LLDG(IPT,-2): pointer to node right of Y-below
C                        or to node left of the node left
C M      : WORK. (NPTS)
C          M(IPT) contains list # of node IPT
C LLS    : OUT. (0:LLS(0))
C          LLS(0) = # independent data dependency lists in ILU
C                   factorization and forward sweep
C          LLS(1:LLS(0)): pointers to the start of a list in LS
C LS     : OUT. (NPTS)
C          LS(ISPT): pointer to node in actual grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IBP, IPT, NBNDS, IB, LB, IDS, IDW, ID, IDE, IDN,
     +   IW, ISW, IS, ISE, LDF, MAXM, MI
C
Ccc Determine for each grid point the # of its data dependency list.
C LLS(MI) contains # nodes in list MI
C M(IPT) contains list # of node IPT
C
C Initialize LLS and M
      DO 1 IPT = 1, NPTS
         LLS(IPT-1) = 0
         M(IPT) = 0
    1 CONTINUE
C
C First list contains independent points, i.e., left/down/front corners
C and internal boundary points.
C For first list the pointers to the nodes in the grid can already be
C stored in LS
      NBNDS = LLBND(0)
      DO 10 IB = 1, NBNDS
C       Store boundary info into work array M
            IBP = 8**ILBND(IB)
            DO 20 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               M(IPT) = M(IPT)+IBP
   20       CONTINUE
   10 CONTINUE
      LDF = 8**1+8**2+8**5
      DO 30 IPT = 1, NPTS
         IF (M(IPT) .EQ. LDF) THEN
C       left/down/front corner, node in starting list
            LLS(1) = LLS(1)+1
            LS(LLS(1)) = IPT
            M(IPT) = 1
         ELSE IF (MOD(INT(M(IPT)/8),8) .EQ. 1) THEN
C       Left boundary, mark node
            M(IPT) = -1
         ELSE
            M(IPT) = 0
         ENDIF
   30 CONTINUE
      IB = NBNDS+1
C       Internal boundary, Dirichlet condition, node in starting list
            DO 40 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLS(1) = LLS(1)+1 
               LS(LLS(1)) = IPT
               M(IPT) = 1
   40       CONTINUE
C
C
C Compute for rest of nodes their list #; a node is dependent from
C its neighbors at compass points W, SW, S and SE in the same plane and
C in the plane below the point directly below and the points
C S, W, E, and N of it.
      MAXM = 0
      DO 50 IPT = 1, NPTS
         IF (M(IPT) .GT. 0) THEN
C       Node already in list
            GOTO 50
         ELSE IF (M(IPT) .LT. 0) THEN
C       Left boundary
            IW = IPT
            M(IPT) = 0
         ELSE
            IW = IPT-1
         ENDIF
         IDS = LLDG(IPT,-9)
         IDW = LLDG(IPT,-8)
         ID  = LLDG(IPT,-7)
         IDE = LLDG(IPT,-6)
         IDN = LLDG(IPT,-5)
         ISW = LLDG(IPT,-4)
         IS  = LLDG(IPT,-3)
         ISE = LLDG(IPT,-2)
         MI = MAX(M(IDS),M(IDW),M(ID),M(IDE),M(IDN),
     +            M(IW),M(ISW),M(IS),M(ISE)) + 1
         M(IPT) = MI
         LLS(MI) = LLS(MI) + 1
         MAXM = MAX(MAXM,MI)
   50 CONTINUE
C
Ccc Store list pointers in LLS and grid pointers in LS
C
C LLS(i):=SUM (# nodes in list_j)
C          j=1,i
      DO 60 IS = 2, MAXM
         LLS(IS) = LLS(IS) + LLS(IS-1)
   60 CONTINUE
C
C Store grid pointers
C LLS(i-1) is pointer to next free place in list i-1 in LS
      LLS(0) = LLS(1)
      DO 70 IPT = 2, NPTS
         IF (M(IPT) .NE. 1) THEN
            MI = M(IPT)
            LLS(MI-1) = LLS(MI-1) + 1
            LS(LLS(MI-1)) = IPT
         ENDIF
   70 CONTINUE
C LLS(i-1) points to list i in LS, should be i-1
      DO 80 IS = MAXM, 1, -1
         LLS(IS) = LLS(IS-1)
   80 CONTINUE
C
Ccc Store # lists in LLS(0)
      LLS(0) = MAXM

      RETURN
      END
      SUBROUTINE JACSUP (NPTS, LLBND, ILBND, LBND, LUDG, M, LLS, LS)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, LLBND(0:*), ILBND(*), LBND(*), LUDG(NPTS,2:9),
     +   M(NPTS), LLS(0:*), LS(NPTS)
C
Ccc PURPOSE:
C Make data-dependency list for backward sweep of backsolve.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C LUDG   : IN. (NPTS,2:9)
C          LUDG(IPT,2): pointer to node left of Y-above
C                       or to node right of the node right
C          LUDG(IPT,3): pointer to node Y-above
C          LUDG(IPT,4): pointer to node right of node Y-above
C                       or to node Y-above Y-above
C          LUDG(IPT,5): pointer to node Y-below Z-above
C          LUDG(IPT,6): pointer to node left of Z-above
C          LUDG(IPT,7): pointer to node Z-above
C          LUDG(IPT,8): pointer to node right of Z-above
C          LUDG(IPT,9): pointer to node Y-above Z-above
C                       or to node Z-above Z-above
C M      : WORK. (NPTS)
C          M(IPT) contains list # of node IPT
C LLS    : OUT. (0:LLS(0))
C          LLS(0) = # independent data dependency lists in
C                   backward sweep
C          LLS(1:LLS(0)): pointers to the start of a list in LS
C LS     : OUT. (NPTS)
C          LS(ISPT): pointer to node in actual grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IBP, IPT, NBNDS, IB, LB, IUS, IUW, IU, IUE, IUN,
     +   IE, INW, IN, INE, IS, RUB, MAXM, MI
C
Ccc Determine for each grid point the # of its data dependency list.
C LLS(MI) contains # nodes in list MI
C M(IPT) contains list # of node IPT
C
C Initialize LLS and M
      DO 1 IPT = 1, NPTS
         LLS(IPT-1) = 0
         M(IPT) = 0
    1 CONTINUE
C
C First list contains independent points, i.e., right/up/back corners
C and internal boundary points.
C For first list the pointers to the nodes in the grid can already be
C stored in LS
      NBNDS = LLBND(0)
      DO 10 IB = 1, NBNDS
C       Store boundary info into work array M
            IBP = 8**ILBND(IB)
            DO 20 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               M(IPT) = M(IPT)+IBP
   20       CONTINUE
   10 CONTINUE
      RUB = 8**3+8**4+8**6
      DO 30 IPT = 1, NPTS
         IF (M(IPT) .EQ. RUB) THEN
C       right/up/back corner, node in starting list
            LLS(1) = LLS(1)+1
            LS(LLS(1)) = IPT
            M(IPT) = 1
         ELSE IF (MOD(INT(M(IPT)/8**3),8) .EQ. 1) THEN
C       Right plane, mark node
            M(IPT) = -1
         ELSE
            M(IPT) = 0
         ENDIF
   30 CONTINUE
      IB = NBNDS+1
C       Internal boundary, Dirichlet condition, node in starting list
            DO 40 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               LLS(1) = LLS(1)+1 
               LS(LLS(1)) = IPT
               M(IPT) = 1
   40       CONTINUE
C
C
C Compute for rest of nodes their list #; a node is dependent from
C its neighbors at compass points E, NW, N and NE in the same plane and
C in the plane above the point directly below and the points
C S, W, E, and N of it.
      MAXM = 0
      DO 50 IPT = NPTS, 1, -1
         IF (M(IPT) .GT. 0) THEN
C       Node already in list
            GOTO 50
         ELSE IF (M(IPT) .LT. 0) THEN
C       Right boundary
            IE = IPT
            M(IPT) = 0
         ELSE
            IE = IPT+1
         ENDIF
         INW = LUDG(IPT,2)
         IN  = LUDG(IPT,3)
         INE = LUDG(IPT,4)
         IUS = LUDG(IPT,5)
         IUW = LUDG(IPT,6)
         IU  = LUDG(IPT,7)
         IUE = LUDG(IPT,8)
         IUN = LUDG(IPT,9)
         MI = MAX(M(IE),M(INW),M(IN),M(INE),
     +            M(IUS),M(IUW),M(IU),M(IUE),M(IUN)) + 1
         M(IPT) = MI
         LLS(MI) = LLS(MI) + 1
         MAXM = MAX(MAXM,MI)
   50 CONTINUE
C
Ccc Store list pointers in LLS and grid pointers in LS
C
C LLS(i):=SUM (# nodes in list_j)
C          j=1,i
      DO 60 IS = 2, MAXM
         LLS(IS) = LLS(IS) + LLS(IS-1)
   60 CONTINUE
C
C Store grid pointers
C LLS(i-1) is pointer to next free place in list i-1 in LS
      LLS(0) = LLS(1)
      DO 70 IPT = NPTS-1, 1, -1
         IF (M(IPT) .NE. 1) THEN
            MI = M(IPT)
            LLS(MI-1) = LLS(MI-1) + 1
            LS(LLS(MI-1)) = IPT
         ENDIF
   70 CONTINUE
C LLS(i-1) points to list i in LS, should be i-1
      DO 80 IS = MAXM, 1, -1
         LLS(IS) = LLS(IS-1)
   80 CONTINUE
C
Ccc Store # lists in LLS(0)
      LLS(0) = MAXM

      RETURN
      END
      SUBROUTINE BICGST (NPTS, NPDE, A, X, B, WT, TOL, ITMAX,
     +   ALU, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   LUN, R, R0, P, T, V, ITER, ERR, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, ITMAX,
     +   LLDG(*), LUDG(*), LSL(*), LLSL(*), LSU(*), LLSU(*),
     +   LUN, ITER, IERR
      DOUBLE PRECISION A(*), X(NPTS*NPDE), B(NPTS*NPDE), WT(NPTS*NPDE),
     +   TOL, ALU(*),
     +   R(NPTS*NPDE), R0(NPTS*NPDE), P(NPTS*NPDE),
     +   T(NPTS*NPDE), V(NPTS*NPDE), ERR
C
Ccc PURPOSE:
C Solve a Non-Symmetric linear system Ax = b using the Preconditioned
C BiConjugate Gradient STAB method. Preconditioning is done with an
C Incomplete LU factorization of A.
C Actually solved is the system [P^(-1).A].x = [P^(-1).b]
C until ||residual||_WRMS < TOL.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C A      : IN. -I
C LLDG   : IN   I These arrays hold the matrix A in block diagonal
C LUDG   : IN  -I storage mode (see description in MVDIAG)
C X      : OUT. Final approximate solution.
C B      : IN. Right-hand side vector.
C WT     : IN. Contains weight factors to compute weighted norm.
C TOL    : IN. System is considered to be solved if
C          weighted max. norm < TOL
C ITMAX  : IN. Maximum number of iterations.
C ALU    : IN. -I
C LSL    : IN   I These arrays should hold the ILU factorization of A in
C LLSL   : IN   I diagonal storage mode and the data dependency lists
C LSU    : IN   I for the forward and the backward solve
C LLSU   : IN. -I (see description in BCKSLV)
C LUN    : IN. Logical unit # of file on which to write the error at
C          each iteration, if this is desired for monitoring convergence
C          If LUN = 0, no writing will occur.
C R      : WORK. (NPTS*NPDE)
C R0     : WORK. (NPTS*NPDE)
C P      : WORK. (NPTS*NPDE)
C T      : WORK. (NPTS*NPDE)
C V      : WORK. (NPTS*NPDE)
C ITER   : OUT. Number of iterations required to reach convergence, or 
C          ITMAX+1 if convergence criterion could not be achieved in 
C          ITMAX iterations.
C ERR    : OUT. Weighted max. norm of error estimate in final
C          approximate solution
C IERR   : OUT. Error return flag
C          0: OK
C          1: Method failed to converge in ITMAX steps
C          2: Breakdown of the method detected (<R0,R> ~ 0.0)
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION WDNRM2, DDOT
      EXTERNAL MVDIAG, BCKSLV, RCOPY, WDNRM2, DDOT, DAXPY
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
      INTEGER I, N
      DOUBLE PRECISION ALPHA, BETA, OMEGA, RHOIM1, RHOI, SXMIN, TNRM2
C
      N = NPTS*NPDE
      ITER = 0
      IERR = 0
      SXMIN = SQRT(XMIN)
C         
Ccc Initialize X and set initial residual to B
      CALL ZERO (N, X)
      DO 10 I = 1, N
         R0(I)  = B(I)
   10 CONTINUE
      CALL BCKSLV (NPTS, NPDE, ALU, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   R0)
C
Ccc Check stopping criterion         
      ERR = WDNRM2 (N, R0, WT)
      IF (LUN .NE. 0)  THEN
         WRITE(LUN,'(''ILU preconditioned BiCGStab for N ='',I6)') N
         WRITE(LUN,'('' ITER  Error Estimate         Alpha'',
     +            ''            Beta            Omega'')')
         WRITE(LUN,'(I5,E16.7)') ITER, ERR
      ENDIF
      IF (ERR .LT. TOL) RETURN
C         
Ccc BiCGStab loop
      CALL RCOPY (N, R0, R)
      DO 100 ITER = 1, ITMAX
C Compute innerproduct original residual with previous residual
         RHOI = DDOT(N, R0, 1, R, 1)
C Calculate coefficient Beta and direction vector Pi
         IF( ITER.EQ.1 ) THEN
            DO 110 I = 1, N
               BETA = 0.0
               P(I) = R(I)
  110       CONTINUE
         ELSE
            BETA = RHOI/RHOIM1*ALPHA/OMEGA
            DO 120 I = 1, N
               P(I) = R(I) + BETA*(P(I)-OMEGA*V(I))
  120       CONTINUE
         ENDIF
C Calculate Vi and coefficient Alfa
         CALL MVDIAG (NPTS, NPDE, A, P, LLDG, LUDG, V)
         CALL BCKSLV (NPTS, NPDE, ALU, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +      V)
         ALPHA = RHOI / DDOT(N, R0, 1, V, 1)
C Calculate polynomial coefficient Omega_i
C store intermediate vector S in R
         DO 130 I = 1, N
            R(I) = R(I) - ALPHA*V(I)
  130    CONTINUE
         CALL MVDIAG (NPTS, NPDE, A, R, LLDG, LUDG, T)
         CALL BCKSLV (NPTS, NPDE, ALU, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +      T)
         TNRM2 = DDOT(N,T,1,T,1)
         IF (TNRM2 .LT. SXMIN) THEN
C    Lucky breakdown
            OMEGA = 0.0
         ELSE
            OMEGA = DDOT(N,T,1,R,1) / TNRM2
         ENDIF
C Adapt Xi = Xi-1 + Alfa*Pi + Omega_i*S.
         CALL DAXPY (N, ALPHA, P, 1, X, 1)
         CALL DAXPY (N, OMEGA, R, 1, X, 1)
C Compute residual R = S - Omega_i*T
         DO 140 I = 1, N
            R(I) = R(I) - OMEGA*T(I)
  140    CONTINUE
C         
C Check stopping criterion.
         ERR = WDNRM2 (N, R, WT)
         IF(LUN .NE. 0)
     +      WRITE(LUN,'(I5,4E16.7)') ITER, ERR, ALPHA, BETA, OMEGA
         IF (ERR .LT. TOL) RETURN
C
C Check if last residual is not parallel to original residual
         IF (ABS(RHOI) .LT. SXMIN) GOTO 990
         RHOIM1 = RHOI
  100 CONTINUE
C         
Ccc end of BiCGStab loop
C
Ccc Stopping criterion not satisfied
      ITER = ITMAX + 1
      IERR = 1
      RETURN
C
Ccc Breakdown of method detected.
  990 IERR = 2
      RETURN
C
      END
      SUBROUTINE MVDIAG (NPTS, NPDE, AD, X, LLDG, LUDG, Y)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLDG(NPTS,-9:-2), LUDG(NPTS,2:9)
      DOUBLE PRECISION AD(NPTS,NPDE,NPDE,-9:9), X(NPTS,NPDE),
     +   Y(NPTS,NPDE)
C
Ccc PURPOSE:
C Compute y = Ax where A is stored in block 19-diagonal mode.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C AD     : IN. A(.,1:NPDE,1:NPDE,-9:-1)  : lower block diagonals
C              A(.,1:NPDE,1:NPDE,0)      :  main block diagonal
C              A(.,1:NPDE,1:NPDE, 1: 9)  : upper block diagonals
C X      : IN. Multiplying vector
C LLDG   : IN. (NPTS,-9:-2)
C          LLDG(IPT,-4): pointer to node left of node below
C                        or to node below the node below
C          LLDG(IPT,-3): pointer to node below
C          LLDG(IPT,-2): pointer to node right of node below
C                        or to node left of the node left
C LUDG   : IN. (NPTS,2:4)
C          LUDG(IPT,2): pointer to node left of node above
C                       or to node right of the node right
C          LUDG(IPT,3): pointer to node above
C          LUDG(IPT,4): pointer to node right of node above
C                       or to node above the node above
C          If one of the above nodes does not exist, the pointer is
C          to the node itself.
C Y      : OUT. Result vector
C
Ccc EXTERNALS USED:
      EXTERNAL ZERO
C
C ----------------------------------------------------------------------
C
      INTEGER IC, JC, IPT, JD
C
      CALL ZERO (NPTS*NPDE, Y)
C
      DO 10 JC = 1, NPDE
      DO 10 IC = 1, NPDE
         DO 20 IPT = 1, NPTS
            Y(IPT,IC) = Y(IPT,IC) + AD(IPT,IC,JC, 0)*X(IPT,JC)
   20    CONTINUE
         DO 30 IPT = 1, NPTS-1
            Y(IPT,IC) = Y(IPT,IC) + AD(IPT,IC,JC,+1)*X(IPT+1,JC)
   30    CONTINUE
         DO 40 IPT = 2, NPTS
            Y(IPT,IC) = Y(IPT,IC) + AD(IPT,IC,JC,-1)*X(IPT-1,JC)
   40    CONTINUE
C
C The next loops can be done for all points, because if an entry
C in the Jacobian does not exist in reality the value in AD is zero
C and the pointer in LUDG or LLDG points to the node itself.
         DO 60 JD = 2, 9
         DO 60 IPT = 1, NPTS
            Y(IPT,IC) = Y(IPT,IC) + AD(IPT,IC,JC,JD)*X(LUDG(IPT,JD),JC)
   60    CONTINUE
         DO 70 JD = -2, -9, -1
         DO 70 IPT = 1, NPTS
            Y(IPT,IC) = Y(IPT,IC) + AD(IPT,IC,JC,JD)*X(LLDG(IPT,JD),JC)
   70    CONTINUE
   10 CONTINUE

      RETURN
      END
      SUBROUTINE INTGRC (ISTRUC, X, Y, Z, NPDE, UIB, UNP1, UN, UNM1,
     +   RELTOL, ABSTOL, TN, DT, DTRAT, DX, DY, DZ, WT, F, CORR, RWORK,
     +   IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER ISTRUC(0:*), NPDE, IERR
      DOUBLE PRECISION X(*), Y(*), Z(*), UIB(*), UNP1(0:*), UN(0:*),
     +   UNM1(0:*),
     +   RELTOL(NPDE), ABSTOL(NPDE),
     +   TN, DT, DTRAT, DX, DY, DZ,
     +   WT(*), F(*), CORR(*), RWORK(*)
C
Ccc PURPOSE:
C Integration in time with BDF2 (first timestep BE).
C Solve nonlinear system F(Tn+1, Un+1, Udot) = 0 with matrix-free
C Newton.
C Solve linear systems with (block) diagonally scaled GCRO.
C
Ccc PARAMETER DESCRIPTION:
C ISTRUC : IN. Data structure Un+1 grid.
C X,Y,Z  : IN. Physical coordinates grid.
C NPDE   : IN. # PDE components
C UIB    : IN. Dirichlet boundary values on internal boundary.
C UNP1   : INOUT. On entry: Initial solution, on exit final solution
C          Newton converged
C UN     : IN. Solution at Tn on Un+1 grid
C UNM1   : IN. Solution at Tn-1 on Un+1 grid
C RELTOL : IN. Relative tolerance for Newton process
C ABSTOL : IN. Absolute tolerance for Newton process
C TN     : IN. Previous time
C DT     : IN. Current time step
C DX     : IN. Current grid spacing in X-direction
C DY     : IN. Current grid spacing in Y-direction
C DZ     : IN. Current grid spacing in Z-direction
C DTRAT  : IN. If BE: 0, if BDF2: DT/DT_old
C WT     : WORK. (NPTS*NPDE)
C          Weight function for norm computation
C F      : WORK. (NPTS*NPDE)
C          Residual
C CORR   : WORK. (NPTS*NPDE)
C          Correction in Newton iteration
C RWORK  : WORK. (RESWRK+LSSWRK)
C          RESWRK: LENU.10
C          LSSWRK: MAX(LENPWK,LENU.(2.MAXLR+MAXL+6))+LENPRE+
C                  MAXLR.MAXLR+(MAXL+3).MAXL+1
C          LENPRE: ( IDIAGP <= 1 ! LENU.NPDE ! LENU )
C          LENPWK: ( IDIAGP = 0 ! LENU.(NPDE.7+2)+NPTS
C                  |:IDIAGP = 1 ! LENU.(NPDE.4+2)+NPTS
C                  |:IDIAGP = 2 ! LENU.10
C                  |:IDIAGP = 3 ! LENU.7 )
C          LENU  : NPTS.NPDE
C IERR   : OUT.
C           0: OK.
C          10: Newton process did not converge
C
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION MAXNRM, WDNRM2
      EXTERNAL ERRWGT, GCRO, PINIT, MAXNRM, RESID, WDNRM2
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
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 15)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
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
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
C
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ, 
     +   LIWK, LENU,
     +   LUT, LUX, LUY, LUZ, LUXX, LUYY, LUZZ, LUXY, LUXZ, LUYZ, 
     +   LPREC, LR, LU, LC, LZW, LRWK,
     +   NPRE, NRES, I, NIT, ITER
      LOGICAL BDPREC, NEWPRE
      DOUBLE PRECISION A0, ERR, CORNRM, OLDNRM, RATE, TOL, UNRM
C
      IERR = 0
C
      A0 = (1+2*DTRAT) / ((1+DTRAT)*DT)
      BDPREC = IDIAGP .LE. 1
C
      IF (LUNNLS .NE. 0) THEN
         WRITE(LUNNLS,'(''Nonlinear system solver at T ='',E16.7)')
     +      TN+DT
      ENDIF
C
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBDPTS = ISTRUC(LLLBND+NBNDS+1)-1
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
C
      LIWK   = LLABVZ+NPTS+1
C
      LENU   = NPTS*NPDE
C
      LUT    = 1
      LUX    = LUT  + LENU
      LUY    = LUX  + LENU
      LUZ    = LUY  + LENU
      LUXX   = LUZ  + LENU
      LUYY   = LUXX + LENU
      LUZZ   = LUYY + LENU
      LUXY   = LUZZ + LENU
      LUXZ   = LUXY + LENU
      LUYZ   = LUXZ + LENU
C
      LPREC  = LUYZ+LENU
      IF (BDPREC) THEN
C Block-diagonal preconditioner
         LR     = LPREC + LENU*NPDE
      ELSE
C Diagonal preconditioner
         LR     = LPREC + LENU
      ENDIF
      LU     = LR + LENU
      LC     = LU + (LENU*MAXLR)
      LZW    = LC + (LENU*MAXLR)
      LRWK   = LZW+ (MAXLR*MAXLR)
C
Ccc Set error weights for use in Newton process
      CALL ERRWGT (NPTS, NPDE, UNP1(1), RELTOL, ABSTOL, WT)
C
Ccc Compute weighted norm of initial solution for convergence check
      UNRM = WDNRM2 (LENU, UNP1(1), WT)
C
Ccc Compute derivatives and residual
      CALL RESID (TN+DT, X, Y, Z, NPTS, NPDE, UNP1, UN, UNM1, DT, DTRAT,
     +   UIB, ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY), ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +   DX, DY, DZ,
     +   RWORK(LUT),  RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ), F)
      NRES = 1
      IF (LUNNLS .NE. 0) THEN
         WRITE(LUNNLS,'('' Max. and WRMS norm residual='',2E16.7)')
     +      MAXNRM(LENU, F), WDNRM2 (LENU, F, WT)
      ENDIF
C
Ccc Compute preconditioner: (block-)diagonal of Jacobian G = dF/dU.
C Store LU-decomposition in PREC, main diagonal inverted.
      CALL PINIT (NPTS, NPDE, F, TN+DT, X, Y, Z, DT, DTRAT, DX, DY, DZ,
     +   UNP1(1), ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP), UIB,
     +   RWORK(LUT), RWORK(LUX), RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ),
     +   ABSTOL, RWORK(LR), IDIAGP, RWORK(LPREC))
      NEWPRE = .TRUE.
      NPRE = 1
C
Ccc Newton iteration loop
    9 CONTINUE
      DO 10 NIT = 1, MAXNIT
C
Cccccc Solve G.corr = F. Store the residual in F.
         TOL = TOLLSC / (2**NIT)
         CALL GCRO (LENU, CORR, F, WT, TOL, BDPREC, RWORK(LPREC),
     +      ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +      TN+DT, A0, DX, DY, DZ, RWORK,
     +      NRRMAX, MAXLR, MAXL, LUNLSS,
     +      RWORK(LR), RWORK(LU), RWORK(LC), RWORK(LZW), RWORK(LRWK),
     +      ITER, ERR, IERR)
         NLSIT(LEVEL,NIT) = NLSIT(LEVEL,NIT)+ ITER
         IF (IERR .NE. 0) GOTO 100
C
Cccccc Test for convergence
         CORNRM = WDNRM2 (LENU, CORR, WT)
         IF (LUNNLS .NE. 0) THEN
            WRITE(LUNNLS,'('' NI:'',I3,'', NLI:'',I4,'', ERLI'':,E16.7,
     +         '', ERNI:'',E16.7)') NIT, ITER, ERR, CORNRM
         ENDIF
         IF (CORNRM .LE. 100*UROUND*UNRM) GOTO 900
         IF (.NOT. NEWPRE) THEN
            RATE = SQRT(CORNRM/OLDNRM)
            IF (RATE .GT. 0.9) THEN
C          Divergence
               GOTO 100
            ELSE IF (RATE/(1-RATE)*CORNRM .LE. TOLNEW) THEN
C          Convergence
               GOTO 900
            ENDIF
         ENDIF
         OLDNRM = CORNRM
C
Ccccc Update solution
         DO 20 I = 1, LENU
            UNP1(I) = UNP1(I) - CORR(I)
   20    CONTINUE
C
Ccc Compute derivatives and residual and start next iteration
         CALL RESID (TN+DT, X, Y, Z, NPTS, NPDE, UNP1, UN, UNM1,
     +      DT, DTRAT, UIB,
     +      ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +      ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +      ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +      DX, DY, DZ,
     +      RWORK(LUT),  RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +      RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +      RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ), F)
         NRES   = NRES+1
         NEWPRE = .FALSE.
         IF (LUNNLS .NE. 0) THEN
            WRITE(LUNNLS,'('' Max. and WRMS norm residual='',2E16.7)')
     +         MAXNRM(LENU, F), WDNRM2 (LENU, F, WT)
         ENDIF
   10 CONTINUE
Ccc End Newton iteration loop
C
Ccc No convergence in max. # iterations
C
Ccccc Check if preconditioner is recent
  100 CONTINUE
         IF (.NOT. NEWPRE .AND. NPRE .LT. MAXJAC) THEN
            IF (LUNNLS .NE. 0) THEN
               WRITE(LUNNLS,'('' New preconditioner, NIT='',I4)') NIT
            ENDIF
C       Compute new preconditioner and retry
C       Compute space derivatives anew since they are disturbed by
C       MVDIFF
            CALL DERIVS (NPTS, NPDE, UNP1,
     +         ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +         ISTRUC(LLBLWY), ISTRUC(LLABVY),
     +         ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +         DX, DY, DZ,
     +         RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +         RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +         RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ))
            CALL PINIT (NPTS, NPDE, F, TN+DT, X, Y, Z, DT, DTRAT,
     +         DX, DY, DZ, UNP1(1),
     +         ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP), UIB,
     +         RWORK(LUT), RWORK(LUX), RWORK(LUY), RWORK(LUZ),
     +         RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +         RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ),
     +         ABSTOL, RWORK(LR), IDIAGP, RWORK(LPREC))
            NEWPRE = .TRUE.
            NPRE = NPRE + 1
            GOTO 9
         ELSE
C       Newton failure
            IERR = 10
            NNIT(LEVEL)   = NNIT(LEVEL)+NIT
            NRESID(LEVEL) = NRESID(LEVEL)+NRES
            NJACS(LEVEL)  = NJACS(LEVEL)+NPRE
            IF (LUNNLS .NE. 0) THEN
               WRITE(LUNNLS,'(''Newton failure, NIT='',I4)') NIT
            ENDIF
            RETURN
         ENDIF
C
Ccc Nonlinear proces has been solved
  900 CONTINUE
C Update solution
      DO 30 I = 1, LENU
         UNP1(I) = UNP1(I) - CORR(I)
   30 CONTINUE
C
      NNIT(LEVEL)   = NNIT(LEVEL)+NIT
      NRESID(LEVEL) = NRESID(LEVEL)+NRES
      NJACS(LEVEL)  = NJACS(LEVEL)+NPRE
C
      RETURN
      END
      SUBROUTINE PINIT (NPTS, NPDE, F, T, X, Y, Z, DT, DTRAT,
     +   DX, DY, DZ, U, LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, ABSTOL, WORK,
     +   IDIAGP, PREC)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*), IDIAGP
      DOUBLE PRECISION F(*), T, X(*), Y(*), Z(*), DT, DTRAT, DX, DY, DZ,
     +   U(*), UIB(*), UT(*), UX(*), UY(*), UZ(*),
     +   UXX(*), UYY(*), UZZ(*), UXY(*), UXZ(*), UYZ(*),
     +   ABSTOL(*), WORK(*), PREC(NPTS,NPDE,*)
C 
Ccc PURPOSE:
C Store the LU-decomposition of the (block-)diagonal of the Jacobian
C G = dF/dU in PREC, main diagonal inverted.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C DT     : IN. Current time stepsize
C DTRAT  : IN. 0 or DT/DT_old
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C U      : IN. Solution at T on current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C LBLWY  : IN. (0:NPTS)
C          LBLWY(IPT): pointer to node below in Y-direction in
C                      actual grid
C                       0, if index node is front-plane boundary point
C LABVY  : IN. (0:NPTS)
C          LABVY(IPT): pointer to node above in Y-direction in
C                      actual grid
C                       0, if index node is back-plane boundary point
C LBLWZ  : IN. (0:NPTS)
C          LBLWZ(IPT): pointer to node below in Z-direction in
C                      actual grid
C                       0, if index node is down-plane boundary point
C LABVZ  : IN. (0:NPTS)
C          LABVZ(IPT): pointer to node above in Z-direction in
C                      actual grid
C                       0, if index node is up-plane boundary point
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C ABSTOL : IN. Absolute tolerance for Newton process
C WORK   : WORK. ( ( IDIAGP = 0 ! LENU.(NPDE.7+2)+NPTS
C                  |:IDIAGP = 1 ! LENU.(NPDE.4+2)+NPTS
C                  |:IDIAGP = 2 ! LENU.10
C                  |:IDIAGP = 3 ! LENU.7 ) )
C IDIAGP : IN. Type of preconditioner
C             0: block-diagonal + first order derivatives
C             1: block-diagonal neglecting first order derivatives
C             2: diagonal + first order derivatives
C             3: diagonal neglecting first order derivatives
C PREC   : OUT. LU-decomposition of the (block-)diagonal of the Jacobian
C             G = dF/dU in PREC, main diagonal inverted.
C
Ccc EXTERNALS USED:
      EXTERNAL BLU, DERVF, DERVFB, PREG, PREGB
C
C-----------------------------------------------------------------------
C
      INTEGER LENDEL, LENU, LENFU, LENFU1,
     +   LFU, LFUX, LFUY, LFUZ, LFUXX, LFUYY, LFUZZ, LDEL, LRWK
      LOGICAL PRECFO
      DOUBLE PRECISION A0
C
      PRECFO = IDIAGP .EQ. 0 .OR. IDIAGP .EQ. 2
      LENU  = NPTS*NPDE
      IF (IDIAGP .LE. 1) THEN
         LENDEL = NPTS
         LENFU  = LENU*NPDE
      ELSE
         LENDEL = LENU
         LENFU  = LENU
      ENDIF
      IF (PRECFO) THEN
         LENFU1 = LENFU
      ELSE
         LENFU1 = 0
      ENDIF
C
      LFU   = 1
      LFUX  = LFU   + LENFU
      LFUY  = LFUX  + LENFU1
      LFUZ  = LFUY  + LENFU1
      LFUXX = LFUZ  + LENFU1
      LFUYY = LFUXX + LENFU
      LFUZZ = LFUYY + LENFU
      LDEL  = LFUZZ + LENFU
      LRWK  = LDEL  + LENDEL
      A0 = (1+2*DTRAT) / ((1+DTRAT)*DT)
      IF (IDIAGP .LE. 1) THEN
C
Ccc Compute dF/dU, (dF/dUx, dF/dUy, dF/dUz,) dF/dUxx, dF/dUyy, dF/dUzz
         CALL DERVFB (F, T, X, Y, Z, NPTS, NPDE, U, A0, DT, DX, DY, DZ,
     +      LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +      UXY, UXZ, UYZ, ABSTOL, WORK(LDEL), WORK(LRWK),
     +      PRECFO, WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +      WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ))
C
Ccc Compute block-diagonal
C G = dF/dU + (dF/dUx.dUx/dU + ...) + dF/dUxx.dUxx/dU + ...
         CALL PREGB (NPTS, NPDE, DX, DY, DZ,
     +      LLBND, ILBND, LBND,
     +      WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +      WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ),
     +      PRECFO, PREC)
C
Ccc Store LU of G in PREC, invert main diagonal
         CALL BLU (NPTS, NPDE, PREC)
      ELSE
C
Ccc Compute dF/dU, (dF/dUx, dF/dUy, dF/dUz,) dF/dUxx, dF/dUyy, dF/dUzz
         CALL DERVF (F, T, X, Y, Z, NPTS, NPDE, U, A0, DT, DX, DY, DZ,
     +      LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +      UXY, UXZ, UYZ, ABSTOL, WORK(LDEL), WORK(LRWK),
     +      PRECFO, WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +      WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ))
C
Ccc Compute G = dF/dU + (dF/dUx.dUx/dU + ...) + dF/dUxx.dUxx/dU + ...
C   Store inverted in PREC
         CALL PREG (NPTS, NPDE, DX, DY, DZ,
     +      LLBND, ILBND, LBND,
     +      WORK(LFU), WORK(LFUX), WORK(LFUY), WORK(LFUZ),
     +      WORK(LFUXX), WORK(LFUYY), WORK(LFUZZ),
     +      PRECFO, PREC)
      ENDIF
C
      RETURN
      END
      SUBROUTINE DERVFB (F, T, X, Y, Z, NPTS, NPDE, U,
     +   A0, DT, DX, DY, DZ,
     +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +   UXY, UXZ, UYZ, ABSTOL, DEL, WORK,
     +   PRECFO, FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      LOGICAL PRECFO
      DOUBLE PRECISION F(NPTS*NPDE), T, X(*), Y(*), Z(*), U(*), A0, DT,
     +   DX, DY, DZ,
     +   UIB(*), UT(*), UX(*), UY(*), UZ(*), UXX(*), UYY(*), UZZ(*),
     +   UXY(*), UXZ(*), UYZ(*), ABSTOL(*),
     +   DEL(NPTS), WORK(2*NPTS*NPDE),
     +   FU(NPTS*NPDE,NPDE),
     +   FUX(NPTS*NPDE,NPDE), FUY(NPTS*NPDE,NPDE), FUZ(NPTS*NPDE,NPDE),
     +   FUXX(NPTS*NPDE,NPDE), FUYY(NPTS*NPDE,NPDE),
     +   FUZZ(NPTS*NPDE,NPDE)
C
Ccc PURPOSE:
C Compute derivatives of residual wrt (derivatives of) U by numerical
C differencing
C
C PARAMETER DESCRIPTION:
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C U      : IN. Solution at T on current grid
C A0     : IN. Coefficient of U_n+1 in time derivative
C DT     : IN. Current time step size
C DX     : IN. Cell width in X-direction for current grid
C DY     : IN. Cell width in Y-direction for current grid
C DZ     : IN. Cell width in Z-direction for current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C ABSTOL : IN. Absolute tolerance for Newton process
C DEL    : WORK. (NPTS)
C WORK   : WORK. (2.LENU)
C PRECFO : IN. If FALSE first order derivatives may be neglected
C FU     : OUT. dF(U,Ut)dU
C FUX    : OUT. dF(Ux)dUx
C FUY    : OUT. dF(Uy)dUy
C FUZ    : OUT. dF(Uz)dUz
C FUXX   : OUT. dF(Uxx)dUxx
C FUYY   : OUT. dF(Uyy)dUyy
C FUZZ   : OUT. dF(Uzz)dUzz
C
Ccc EXTERNALS USED:
      EXTERNAL PERTRB, PRTRBU, RES
C
C-----------------------------------------------------------------------
C
      INTEGER I, IC, ICPTB, IPT, LUTBAR
      DOUBLE PRECISION FACX, FACY, FACZ, FACXX, FACYY, FACZZ, TOL

      LUTBAR = 1 + NPTS*NPDE
C
Ccc How to decide if derivatives are `zero'?
C Take `zero'-value of U divided by the grid width
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACZ  = 1/(2*DZ)
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
C
Ccc Loop over the components of the (derivatives of) U
      DO 10 ICPTB = 1, NPDE
C
C dF(U,Ut)/dU
         TOL = ABSTOL(ICPTB)
         CALL PRTRBU (ICPTB, NPTS, NPDE, U, A0, DT, UT, TOL, DEL,
     +      WORK, WORK(LUTBAR))
         CALL RES (T, X, Y, Z, NPTS, NPDE, WORK,
     +      LLBND, ILBND, LBND, UIB,
     +      WORK(LUTBAR), UX, UY, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FU(1,ICPTB))
         DO 20 IC = 1, NPDE
         DO 20 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FU(I,ICPTB) = (FU(I,ICPTB) - F(I)) / DEL(IPT)
   20    CONTINUE
         
         IF (PRECFO) THEN
C
C dF(Ux)/dUx
         TOL = ABSTOL(ICPTB)*FACX
         CALL PERTRB (ICPTB, NPTS, NPDE, UX, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, WORK, UY, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUX(1,ICPTB))
         DO 21 IC = 1, NPDE
         DO 21 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUX(I,ICPTB) = (FUX(I,ICPTB) - F(I)) / DEL(IPT)
   21    CONTINUE
C
C dF(Uy)/dUy
         TOL = ABSTOL(ICPTB)*FACY
         CALL PERTRB (ICPTB, NPTS, NPDE, UY, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, WORK, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUY(1,ICPTB))
         DO 22 IC = 1, NPDE
         DO 22 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUY(I,ICPTB) = (FUY(I,ICPTB) - F(I)) / DEL(IPT)
   22    CONTINUE
C
C dF(Uz)/dUz
         TOL = ABSTOL(ICPTB)*FACZ
         CALL PERTRB (ICPTB, NPTS, NPDE, UZ, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, WORK,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUZ(1,ICPTB))
         DO 23 IC = 1, NPDE
         DO 23 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUZ(I,ICPTB) = (FUZ(I,ICPTB) - F(I)) / DEL(IPT)
   23    CONTINUE
         
         ENDIF
C
C dF(Uxx)/dUxx
         TOL = ABSTOL(ICPTB)*FACXX
         CALL PERTRB (ICPTB, NPTS, NPDE, UXX, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      WORK, UYY, UZZ, UXY, UXZ, UYZ, FUXX(1,ICPTB))
         DO 70 IC = 1, NPDE
         DO 70 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUXX(I,ICPTB) = (FUXX(I,ICPTB) - F(I)) / DEL(IPT)
   70    CONTINUE
C
C dF(Uyy)/dUyy
         TOL = ABSTOL(ICPTB)*FACYY
         CALL PERTRB (ICPTB, NPTS, NPDE, UYY, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      UXX, WORK, UZZ, UXY, UXZ, UYZ, FUYY(1,ICPTB))
         DO 80 IC = 1, NPDE
         DO 80 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUYY(I,ICPTB) = (FUYY(I,ICPTB) - F(I)) / DEL(IPT)
   80    CONTINUE
C
C dF(Uzz)/dUzz
         TOL = ABSTOL(ICPTB)*FACZZ
         CALL PERTRB (ICPTB, NPTS, NPDE, UZZ, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      UXX, UYY, WORK, UXY, UXZ, UYZ, FUZZ(1,ICPTB))
         DO 90 IC = 1, NPDE
         DO 90 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUZZ(I,ICPTB) = (FUZZ(I,ICPTB) - F(I)) / DEL(IPT)
   90    CONTINUE
   10 CONTINUE
      
      RETURN
      END
      SUBROUTINE PREGB (NPTS, NPDE, DX, DY, DZ,
     +   LLBND, ILBND, LBND,
     +   FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ, PRECFO, PREC)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      LOGICAL PRECFO
      DOUBLE PRECISION DX, DY, DZ,
     +   FU(NPTS*NPDE,NPDE),
     +   FUX(NPTS*NPDE,NPDE), FUY(NPTS*NPDE,NPDE), FUZ(NPTS*NPDE,NPDE),
     +   FUXX(NPTS*NPDE,NPDE), FUYY(NPTS*NPDE,NPDE),
     +   FUZZ(NPTS*NPDE,NPDE), PREC(NPTS*NPDE,NPDE)
C
Ccc PURPOSE:
C Compute block-diagonal of Jacobian G = dF/dU using derivatives
C of residual wrt (derivatives of) U.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FU     : IN. Derivative residual F(.,U,Ut,.) wrt U
C FUX    : IN. Derivative residual F(.,Ux,.)  wrt Ux
C FUY    : IN. Derivative residual F(.,Uy,.)  wrt Uy
C FUZ    : IN. Derivative residual F(.,Uz,.)  wrt Uz
C FUXX   : IN. Derivative residual F(.,Uxx,.)  wrt Uxx
C FUYY   : IN. Derivative residual F(.,Uyy,.)  wrt Uyy
C FUZZ   : IN. Derivative residual F(.,Uzz,.)  wrt Uzz
C PRECFO : IN. If FALSE first order derivatives may be neglected
C PREC   : OUT. Block-diagonal of Jacobian.
C
Ccc EXTERNALS USED:
      EXTERNAL PRGBBD
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
      INTEGER I, JC, LENU
      DOUBLE PRECISION FACX, FACY, FACZ, FACXX, FACYY, FACZZ
C
      LENU  = NPTS*NPDE
C
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
C
      DO 10 JC = 1, NPDE
      DO 10 I = 1, LENU
C dF(ipt,ic)/dU(ipt,ic)
         PREC(I,JC) = FU(I,JC) + 
     +      FUXX(I,JC)*(-2*FACXX) + FUYY(I,JC)*(-2*FACYY) +
     +      FUZZ(I,JC)*(-2*FACZZ)
   10 CONTINUE
C
      IF (PRECFO) THEN
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACZ  = 1/(2*DZ)
      CALL PRGBBD (NPTS, NPDE, FACX, FACY, FACZ,
     +   LLBND, ILBND, LBND, FUX, FUY, FUZ, PREC)
      ENDIF
C
      RETURN
      END
      SUBROUTINE PRGBBD (NPTS, NPDE, FACX, FACY, FACZ,
     +   LLBND, ILBND, LBND, FUX, FUY, FUZ, G)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION FACX, FACY, FACZ,
     +   FUX(NPTS,NPDE,NPDE), FUY(NPTS,NPDE,NPDE), FUZ(NPTS,NPDE,NPDE),
     +   G(NPTS,NPDE,NPDE)
C
Ccc PURPOSE:
C Correct Jacobian G = dF/dU for second order approximation of
C first order derivatives at boundaries
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C FACX   : IN. 1/(2*DX)
C FACY   : IN. 1/(2*DY)
C FACZ   : IN. 1/(2*DZ)
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FUX    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Ux
C FUY    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uy
C FUZ    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uz
C G      : INOUT.
C          IN: block-diagonal of Jacobian
C          OUT: corrected for first order derivatives at boundaries
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, IC, JC, IB, LB
C
Ccc Boundary corrections, no corrections needed for internal boundaries
      DO 10 IB = 1, LLBND(0)
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane, correction needed for dF/dUx
            DO 20 JC = 1, NPDE
            DO 20 IC = 1, NPDE
CDIR$ IVDEP
            DO 25 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUX(IPT,IC,JC)*(-3*FACX)
   25       CONTINUE
   20       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down plane, correction needed for dF/dUz
            DO 30 JC = 1, NPDE
            DO 30 IC = 1, NPDE
CDIR$ IVDEP
            DO 35 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUZ(IPT,IC,JC)*(-3*FACZ)
   35       CONTINUE
   30       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane, correction needed for dF/dUx
            DO 40 JC = 1, NPDE
            DO 40 IC = 1, NPDE
CDIR$ IVDEP
            DO 45 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUX(IPT,IC,JC)*(+3*FACX)
   45       CONTINUE
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane, correction needed for dF/dUz
            DO 50 JC = 1, NPDE
            DO 50 IC = 1, NPDE
CDIR$ IVDEP
            DO 55 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
C             dF(ipt,ic)/dU(below(below(ipt)),jc)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUZ(IPT,IC,JC)*(+3*FACZ)
   55       CONTINUE
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front plane, correction needed for dF/dUy
            DO 60 JC = 1, NPDE
            DO 60 IC = 1, NPDE
CDIR$ IVDEP
            DO 65 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUY(IPT,IC,JC)*(-3*FACY)
   65       CONTINUE
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back plane, correction needed for dF/dUy
            DO 70 JC = 1, NPDE
            DO 70 IC = 1, NPDE
CDIR$ IVDEP
            DO 75 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC,JC) = G(IPT,IC,JC) +
     +            FUY(IPT,IC,JC)*(+3*FACY)
   75       CONTINUE
   70       CONTINUE
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BLU (NPTS, NPDE, A)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION A(NPTS,NPDE,NPDE)
C
Ccc PURPOSE:
C LU decomposition of block-diagonal A
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : INOUT.
C          IN: main block diagonal
C          OUT: A(.,ic,jc):      jc < ic: block diagonal of L
C                                         diagonal L == I
C                                jc >=ic: block diagonal of U
C                                         diagonal U inverted
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, JC, LC, N
      DOUBLE PRECISION D
C
         DO 550 IC = 1, NPDE
            DO 554 LC = 1, IC-1
            DO 555 JC = IC, NPDE
CDIR$ IVDEP
            DO 551 N = 1, NPTS
               A(N,IC,JC) = A(N,IC,JC)
     +                        - A(N,IC,LC)*A(N,LC,JC)
  551       CONTINUE
  555       CONTINUE
            DO 556 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 552 N = 1, NPTS
               A(N,JC,IC) = A(N,JC,IC)
     +                        - A(N,JC,LC)*A(N,LC,IC)
  552       CONTINUE
  556       CONTINUE
  554       CONTINUE
CDIR$ IVDEP
            DO 553 N = 1, NPTS
               D = A(N,IC,IC)
               IF (ABS(D) .LT. 1D-7) THEN
                  A(N,IC,IC) = 1.0
               ELSE
                  A(N,IC,IC) = 1.0 / D
               ENDIF
  553       CONTINUE
            DO 557 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 559 N = 1, NPTS
               A(N,JC,IC) = A(N,JC,IC) * A(N,IC,IC)
  559       CONTINUE
  557       CONTINUE
  550    CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKBDI (NPTS, NPDE, A, B, X)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION A(NPTS,NPDE,NPDE), X(NPTS,NPDE), B(NPTS,NPDE)
C
Ccc PURPOSE:
C Solve LUx = b
C    A is a block-diagonal matrix
C A((i,j,k),1:NPDE,1:NPDE) contains a block of NPDE.NPDE elements
C    corresponding with node (i,j,k)
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : IN.  A(.,ic,jc):      jc < ic: block diagonal of L
C                                         diagonal L == I
C                                jc >=ic: block diagonal of U
C                                         diagonal U inverted
C X      : OUT: solution vector x
C B      : IN: right-hand side vector b
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, JC, N
C
      CALL RCOPY (NPTS*NPDE, B, X)
C
CCC Ly = b
C
      DO 100 IC = 2, NPDE
      DO 101 JC = 1, IC-1
CDIR$ IVDEP
      DO 1 N = 1, NPTS
         X(N,IC) = X(N,IC) - A(N,IC,JC)*X(N,JC)
    1 CONTINUE
  101 CONTINUE
  100 CONTINUE
C
CCC Ux = y
C
      DO 130 IC = NPDE, 1, -1
      DO 131 JC = NPDE, IC+1, -1
CDIR$ IVDEP
      DO 132 N = 1, NPTS
         X(N,IC) = X(N,IC) - A(N,IC,JC)*X(N,JC)
  132 CONTINUE
  131 CONTINUE
CDIR$ IVDEP
      DO 133 N = 1, NPTS
         X(N,IC) = X(N,IC) * A(N,IC,IC)
  133 CONTINUE
  130 CONTINUE
C
      RETURN
      END
      SUBROUTINE DERVF (F, T, X, Y, Z, NPTS, NPDE, U,
     +   A0, DT, DX, DY, DZ,
     +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +   UXY, UXZ, UYZ, ABSTOL, DEL, WORK,
     +   PRECFO, FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      LOGICAL PRECFO
      DOUBLE PRECISION F(NPTS*NPDE), T, X(*), Y(*), Z(*), U(*), A0, DT,
     +   DX, DY, DZ,
     +   UIB(*), UT(*), UX(*), UY(*), UZ(*), UXX(*), UYY(*), UZZ(*),
     +   UXY(*), UXZ(*), UYZ(*), ABSTOL(*),
     +   DEL(NPTS*NPDE), WORK(2*NPTS*NPDE),
     +   FU(NPTS*NPDE),
     +   FUX(NPTS*NPDE), FUY(NPTS*NPDE), FUZ(NPTS*NPDE),
     +   FUXX(NPTS*NPDE), FUYY(NPTS*NPDE), FUZZ(NPTS*NPDE)
C
Ccc PURPOSE:
C Compute derivatives of residual wrt (derivatives of) U by numerical
C differencing
C
C PARAMETER DESCRIPTION:
C F      : IN. Residual F(t,U,Ut)
C T      : IN. Current time
C X,Y,Z  : IN. Physical coordinates of gridpoints
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C U      : IN. Solution at T on current grid
C A0     : IN. Coefficient of U_n+1 in time derivative
C DT     : IN. Current time step size
C DX     : IN. Cell width in X-direction for current grid
C DY     : IN. Cell width in Y-direction for current grid
C DZ     : IN. Cell width in Z-direction for current grid
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C UIB    : IN. Solution at T on internal boundaries
C UT     : IN. Time derivative of U on current grid
C UX     : IN. -I
C UY     : IN.  I
C UZ     : IN.  I
C UXX    : IN.  I
C UYY    : IN.  I Space derivatives of U on current grid
C UZZ    : IN.  I
C UXY    : IN.  I
C UXZ    : IN.  I
C UYZ    : IN. -I
C ABSTOL : IN. Absolute tolerance for Newton process
C DEL    : WORK. (LENU)
C WORK   : WORK. (2.LENU)
C PRECFO : IN. If FALSE first order derivatives may be neglected
C FU     : OUT. dF(U,Ut)dU
C FUX    : OUT. dF(Ux)dUx
C FUY    : OUT. dF(Uy)dUy
C FUZ    : OUT. dF(Uz)dUz
C FUXX   : OUT. dF(Uxx)dUxx
C FUYY   : OUT. dF(Uyy)dUyy
C FUZZ   : OUT. dF(Uzz)dUzz
C
Ccc EXTERNALS USED:
      EXTERNAL PERTRG, PRTRGU, RES
C
C-----------------------------------------------------------------------
C
      INTEGER I, LENU, LUTBAR
      DOUBLE PRECISION FACX, FACY, FACZ, FACXX, FACYY, FACZZ

      LENU = NPTS*NPDE
      LUTBAR = 1 + LENU
C
Ccc How to decide if derivatives are `zero'?
C Take `zero'-value of U divided by the grid width
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
C
C dF(U,Ut)/dU
      CALL PRTRGU (NPTS, NPDE, U, A0, DT, UT, ABSTOL, DEL,
     +   WORK, WORK(LUTBAR))
      CALL RES (T, X, Y, Z, NPTS, NPDE, WORK,
     +   LLBND, ILBND, LBND, UIB,
     +   WORK(LUTBAR), UX, UY, UZ,
     +   UXX, UYY, UZZ, UXY, UXZ, UYZ, FU)
      DO 10 I = 1, LENU
         FU(I) = (FU(I) - F(I)) / DEL(I)
   10 CONTINUE

      IF (PRECFO) THEN
C
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACZ  = 1/(2*DZ)
C
C dF(Ux)/dUx
      CALL PERTRG (NPTS, NPDE, UX, ABSTOL, FACX, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, WORK, UY, UZ,
     +   UXX, UYY, UZZ, UXY, UXZ, UYZ, FUX)
      DO 11 I = 1, LENU
         FUX(I) = (FUX(I) - F(I)) / DEL(I)
   11 CONTINUE
C
C dF(Uy)/dUy
      CALL PERTRG (NPTS, NPDE, UY, ABSTOL, FACY, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, WORK, UZ,
     +   UXX, UYY, UZZ, UXY, UXZ, UYZ, FUY)
      DO 12 I = 1, LENU
         FUY(I) = (FUY(I) - F(I)) / DEL(I)
   12 CONTINUE
C
C dF(Uz)/dUz
      CALL PERTRG (NPTS, NPDE, UZ, ABSTOL, FACZ, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, WORK,
     +   UXX, UYY, UZZ, UXY, UXZ, UYZ, FUZ)
      DO 13 I = 1, LENU
         FUZ(I) = (FUZ(I) - F(I)) / DEL(I)
   13 CONTINUE

      ENDIF
C
C dF(Uxx)/dUxx
      CALL PERTRG (NPTS, NPDE, UXX, ABSTOL, FACXX, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ,
     +   WORK, UYY, UZZ, UXY, UXZ, UYZ, FUXX)
      DO 20 I = 1, LENU
         FUXX(I) = (FUXX(I) - F(I)) / DEL(I)
   20 CONTINUE
C
C dF(Uyy)/dUyy
      CALL PERTRG (NPTS, NPDE, UYY, ABSTOL, FACYY, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ,
     +   UXX, WORK, UZZ, UXY, UXZ, UYZ, FUYY)
      DO 30 I = 1, LENU
         FUYY(I) = (FUYY(I) - F(I)) / DEL(I)
   30 CONTINUE
C
C dF(Uzz)/dUzz
      CALL PERTRG (NPTS, NPDE, UZZ, ABSTOL, FACZZ, DEL, WORK)
      CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +   LLBND, ILBND, LBND, UIB,
     +   UT, UX, UY, UZ,
     +   UXX, UYY, WORK, UXY, UXZ, UYZ, FUZZ)
      DO 40 I = 1, LENU
         FUZZ(I) = (FUZZ(I) - F(I)) / DEL(I)
   40 CONTINUE

      RETURN
      END
      SUBROUTINE PRTRGU (NPTS, NPDE, U, A0, DT, UT, ABSTOL, DEL,
     +   UBAR, UTBAR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION U(NPTS,NPDE), A0, DT, UT(NPTS,NPDE),
     +   ABSTOL(NPDE),
     +   DEL(NPTS,NPDE), UBAR(NPTS,NPDE), UTBAR(NPTS,NPDE)
C
Ccc PURPOSE:
C Perturb U. Return perturbance in DEL and perturbed U in UBAR.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # PDE components
C U      : IN. Solution or derivative of solution to be perturbed
C A0     : IN. Coefficient of U_n+1 in time derivative
C DT     : IN. Current time step size
C UT     : IN. Time derivative of U on current grid
C ABSTOL : IN. Absolute tolerance for Newton process
C DEL    : OUT. Perturbation values
C UBAR   : OUT. Perturbed values of U
C UTBAR  : OUT. Perturbed values of UT
C
Ccc EXTERNALS USED:
      EXTERNAL RCOPY
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
      INTEGER IC, IPT
      DOUBLE PRECISION DELI, EPS, TOL

      CALL RCOPY (NPTS*NPDE, U, UBAR)
      CALL RCOPY (NPTS*NPDE, UT, UTBAR)

      EPS = SQRT(UROUND)
      DO 10 IC = 1, NPDE
         TOL = ABSTOL(IC)
         DO 20 IPT = 1, NPTS
C Compute perturbance, if U=0, U(T+dt)=dtUt, if both are zero take
C threshold
            DELI = EPS*MAX(ABS(U(IPT,IC)),ABS(DT*UT(IPT,IC)),TOL)
            DELI = SIGN(DELI,DT*UT(IPT,IC))
C To ensure that the perturbance is the same machine number as the
C denominator
            DEL(IPT,IC) = (U(IPT,IC)+DELI)-U(IPT,IC)
            UBAR(IPT,IC) = U(IPT,IC) + DEL(IPT,IC)
            UTBAR(IPT,IC) = UT(IPT,IC) + A0*DEL(IPT,IC)
   20    CONTINUE
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PERTRG (NPTS, NPDE, U, ABSTOL, FAC, DEL, UBAR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION U(NPTS,NPDE), ABSTOL(NPDE), FAC, DEL(NPTS,NPDE),
     +   UBAR(NPTS,NPDE)
C
Ccc PURPOSE:
C Perturb U. Return perturbance in DEL and perturbed U in UBAR.
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # PDE components
C U      : IN. Derivative of solution to be perturbed
C ABSTOL : IN. Absolute tolerance for Newton process
C FAC    : IN. Grid factor for tolerance to get threshold
C DEL    : OUT. Perturbation values
C UBAR   : OUT. Perturbed values of U
C
Ccc EXTERNALS USED:
      EXTERNAL RCOPY
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
      INTEGER IC, IPT
      DOUBLE PRECISION DELI, EPS, TOL

      CALL RCOPY (NPTS*NPDE, U, UBAR)

      EPS = SQRT(UROUND)
      DO 10 IC = 1, NPDE
         TOL = ABSTOL(IC)*FAC
         DO 20 IPT = 1, NPTS
C Compute perturbance
            DELI = EPS*MAX(ABS(U(IPT,IC)),TOL)
C To ensure that UBAR has the same sign as U
            DELI = SIGN(DELI,U(IPT,IC))
C To ensure that the perturbance is the same machine number as the
C denominator
            DEL(IPT,IC) = (U(IPT,IC)+DELI)-U(IPT,IC)
            UBAR(IPT,IC) = U(IPT,IC) + DEL(IPT,IC)
   20    CONTINUE
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PREG (NPTS, NPDE, DX, DY, DZ,
     +   LLBND, ILBND, LBND,
     +   FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ, PRECFO, DGINV)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      LOGICAL PRECFO
      DOUBLE PRECISION DX, DY, DZ,
     +   FU(NPTS*NPDE),
     +   FUX(NPTS*NPDE), FUY(NPTS*NPDE), FUZ(NPTS*NPDE),
     +   FUXX(NPTS*NPDE), FUYY(NPTS*NPDE), FUZZ(NPTS*NPDE),
     +   DGINV(NPTS*NPDE)
C
Ccc PURPOSE:
C Compute inverse of diagonal of Jacobian G = dF/dU using derivatives
C of residual wrt (derivatives of) U.
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C DX     : IN. Current grid width in X-direction
C DY     : IN. Current grid width in Y-direction
C DZ     : IN. Current grid width in Z-direction
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FU     : IN. Derivative residual F(.,U,Ut,.) wrt U
C FUXX   : IN. Derivative residual F(.,Uxx,.)  wrt Uxx
C FUYY   : IN. Derivative residual F(.,Uyy,.)  wrt Uyy
C FUZZ   : IN. Derivative residual F(.,Uzz,.)  wrt Uzz
C PRECFO : IN. If FALSE first order derivatives may be neglected
C DGINV  : OUT. Inverse of diagonal of Jacobian.
C
Ccc EXTERNALS USED:
      EXTERNAL PREGBD
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
      INTEGER I, LENU
      DOUBLE PRECISION DG, EPS, FACX, FACY, FACZ, FACXX, FACYY, FACZZ
C
      EPS = SQRT(UROUND)
C
      LENU  = NPTS*NPDE
C
      FACXX = 1/DX**2
      FACYY = 1/DY**2
      FACZZ = 1/DZ**2
C
      DO 10 I = 1, LENU
C dF(ipt,ic)/dU(ipt,ic)
         DGINV(I) = FU(I) +
     +      FUXX(I)*(-2*FACXX) + FUYY(I)*(-2*FACYY) +
     +      FUZZ(I)*(-2*FACZZ)
   10 CONTINUE

      IF (PRECFO) THEN
         FACX = 1/(2*DX)
         FACY = 1/(2*DY)
         FACZ = 1/(2*DZ)
         CALL PREGBD (NPTS, NPDE, FACX, FACY, FACZ,
     +      LLBND, ILBND, LBND, FUX, FUY, FUZ, DGINV)
      ENDIF

      DO 20 I = 1, LENU
         DG = DGINV(I)
         IF (ABS(DG) .LT. EPS) THEN
            DGINV(I) = 1.0
         ELSE
            DGINV(I) = 1.0/DG
         ENDIF
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE PREGBD (NPTS, NPDE, FACX, FACY, FACZ,
     +   LLBND, ILBND, LBND, FUX, FUY, FUZ, G)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      DOUBLE PRECISION FACX, FACY, FACZ,
     +   FUX(NPTS,NPDE), FUY(NPTS,NPDE), FUZ(NPTS,NPDE),
     +   G(NPTS,NPDE)
C
Ccc PURPOSE:
C Correct Jacobian G = dF/dU for second order approximation of
C first order derivatives at boundaries
C
C PARAMETER DESCRIPTION:
C NPTS   : IN. # grid points
C NPDE   : IN. # PDE components
C FACX   : IN. 1/(2*DX)
C FACY   : IN. 1/(2*DY)
C FACZ   : IN. 1/(2*DZ)
C LLBND  : IN. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : IN. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : IN. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C FUX    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Ux
C FUY    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uy
C FUZ    : IN. Derivative residual
C              F(t,U,Ut,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz) wrt Uz
C G      : INOUT.
C          IN: Main diagonal of Jacobian
C          OUT: G corrected for first order derivatives at
C               boundaries
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IPT, IC, IB, LB
C
Ccc Boundary corrections, no corrections needed for internal boundaries
      DO 10 IB = 1, LLBND(0)
         IF (ILBND(IB) .EQ. 1) THEN
C       Left plane, correction needed for dF/dUx
            DO 20 IC = 1, NPDE
CDIR$ IVDEP
            DO 25 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUX(IPT,IC)*(-3*FACX)
   25       CONTINUE
   20       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 2) THEN
C       Down plane, correction needed for dF/dUz
            DO 30 IC = 1, NPDE
CDIR$ IVDEP
            DO 35 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUZ(IPT,IC)*(-3*FACZ)
   35       CONTINUE
   30       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 3) THEN
C       Right plane, correction needed for dF/dUx
            DO 40 IC = 1, NPDE
CDIR$ IVDEP
            DO 45 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUX(IPT,IC)*(+3*FACX)
   45       CONTINUE
   40       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 4) THEN
C       Up plane, correction needed for dF/dUz
            DO 50 IC = 1, NPDE
CDIR$ IVDEP
            DO 55 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUZ(IPT,IC)*(+3*FACZ)
   55       CONTINUE
   50       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 5) THEN
C       Front plane, correction needed for dF/dUy
            DO 60 IC = 1, NPDE
CDIR$ IVDEP
            DO 65 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUY(IPT,IC)*(-3*FACY)
   65       CONTINUE
   60       CONTINUE
         ELSE IF (ILBND(IB) .EQ. 6) THEN
C       Back plane, correction needed for dF/dUy
            DO 70 IC = 1, NPDE
CDIR$ IVDEP
            DO 75 LB = LLBND(IB), LLBND(IB+1)-1
               IPT = LBND(LB)
               G(IPT,IC) = G(IPT,IC) + FUY(IPT,IC)*(+3*FACY)
   75       CONTINUE
   70       CONTINUE
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE GCRO (N, XV, F, WT, TOL, BDPREC, PREC,
     +   ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +   TNP1, A0, DX, DY, DZ, RWORK,
     +   NRRMAX, MAXLR, MAXL, LUN,
     +   R, U, C, ZW, WORK, ITER, ERR, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, ISTRUC(0:*), NPDE, NRRMAX, MAXLR, MAXL, LUN,
     +   ITER, IERR
      LOGICAL BDPREC
      DOUBLE PRECISION XV(N), F(N), WT(N), TOL, PREC(N),
     +   X(*), Y(*), Z(*), UIB(*), UNP1(0:*),
     +   TNP1, A0, DX, DY, DZ, RWORK(*),
     +   R(N), U(N,0:MAXLR-1), C(N,0:MAXLR-1), ZW(0:MAXLR-1,0:MAXLR-1),
     +   WORK(*), ERR
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
C
Ccc PURPOSE:
C Solve a Non-Symmetric linear system Gx = F using the matrix-free
C (block)-diagonally scaled GCRO(NRRMAX+1,MAXLR,(MAXL)) method.
C Actually solved is the system [W.D^(-1).G.W^(-1)].[W.x] = [W.D^(-1).F]
C where W = diag(WT(i)) and D is the (block) diagonal of G
C until max(||residual||_2,||GM^(-1).residual||_2) < TOL,
C with GM the projection of the matrix unto the Krylov base obtained
C with the GMRES inner iteration.
C
Ccc PARAMETER DESCRIPTION:
C N      : IN. Dimension of the system
C XV     : OUT. Final approximate solution.
C F      : IN. Right-hand side vector.
C WT     : IN. Contains weight factors to compute weighted norm.
C TOL    : IN. System is considered to be solved if
C          weighted 2-norm < TOL
C BDPREC : IN. if true a block-diagonal preconditioner is used
C PREC   : IN. LU decomposition of (block-)diagonal of G. main diagonal
C          inverted
C ISTRUC : IN. -I Parameters
C  ...          I    for
C RWORK  : IN. -I residual evaluation
C NRRMAX : IN. # restarts outer loop
C MAXLR  : IN. max. iterations outer loop
C MAXL   : IN. max. iterations GMRES (no restarts)
C LUN    : IN. Logical unit # of file on which to write the error at
C          each iteration, if this is desired for monitoring convergence
C          If LUN = 0, no writing will occur.
C R      : WORK.
C U      : WORK.
C C      : WORK.
C ZW     : WORK.
C WORK   : WORK. (N.(MAXL+1)+(MAXL+3).MAXL+4.N+1)
C ITER   : OUT. Number of iterations required to reach convergence, or
C          until (NRRMAX+1).MAXLR. outer loop iterations have been
C          performed. ITER is the sum of the number of outerloop
C          iterations + number of GMRES (preconditioner) iterations.
C ERR    : OUT. Weighted 2-norm of error estimate in final
C          approximate solution
C IERR   : OUT. Error return flag
C          0: OK
C          1: Method failed to converge in (NRRMAX+1).MAXLR. outer loop
C             iterations
C          2: Break down in outer loop
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL BCKBDI, GMRESO, MVDIFF, RCOPY, DAXPY, DDOT, DNRM2, ZERO
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
      INTEGER I, IR, J, K, GRITER, PRITER, GMITER,
     +   LV, LHES, LQ, LRWRK
      DOUBLE PRECISION RNRM, UNRM
C
Ccc Distribute workspace for GMRES
      LV   = 1
      LHES = LV + N*(MAXL+1)
      LQ   = LHES + (MAXL+1)*MAXL
      LRWRK= LQ + 2*MAXL
C
      ITER   = 0
      IERR   = 0
      GRITER = 0
      PRITER = 0
C
Ccc Initialize X and set initial residual to r_0 = W.D^(-1).F
      CALL ZERO (N, XV)
      IF (BDPREC) THEN
         CALL BCKBDI (N/NPDE, NPDE, PREC, F, R)
         DO 10 I = 1, N
            R(I) = WT(I)*R(I)
   10    CONTINUE
      ELSE
         DO 11 I = 1, N
            R(I) = WT(I)*PREC(I)*F(I)
   11    CONTINUE
      ENDIF
C
C
Ccc Check stopping criterion
      ERR = DNRM2 (N, R, 1)
      IF (LUN .NE. 0)  THEN
         WRITE(LUN,*)
     +      'Diag. scaled GCRO(NRRMAX,MAXLR))'
         WRITE(LUN,'(''NRRMAX, MAXLR, N:'',3I10)')
     +      NRRMAX, MAXLR, N
         WRITE(LUN,*) '# it. GCRO  # it.GMRES   Error Estimate'
         WRITE(LUN,'(2I10,E20.7)') GRITER, PRITER, ERR
      ENDIF
      IF (ERR .LT. SQRT(UROUND)*TOL) RETURN
C
Ccc Restart loop
      DO 150 IR = 0, NRRMAX
         IERR   = 0
C
Ccc Outer loop
      DO 200 K = 0, MAXLR-1
C
C Perform the diagonally scaled GMRES algorithm to solve
C    (I-C_k-1C_k-1^T).[W.D^(-1).G.W^(-1)].u_k = A_k-1.u_k = r_k-1,
C    r_k = r_k-1 - A_k-1.u_k
C    c_k = (r_k-1 - r_k)/<c_k,c_k>
C to compute the initial preconditioner for the E-N process
C If u_k is solved accurately enough ||u_k = H_k.r_k||_2 is a good
C measure for the error when solving  DAE systems.
         CALL GMRESO (N, U(1,K), R, WT, BDPREC, PREC, C, K, TOL, MAXL,
     +      F, ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +      TNP1, A0, DX, DY, DZ, RWORK,
     +      WORK(LV), WORK(LHES), WORK(LQ), WORK(LRWRK),
     +      LUN, GMITER, ERR, IERR)
         PRITER = PRITER + GMITER
         IF (LUN .GT. 0) THEN
            WRITE(LUN,*) 'Result GMRES:', GMITER, TOL, ERR, IERR
         ENDIF
         IF (IERR .GT. 1) THEN
            PRINT *, 'wat nu?'
            STOP
         ENDIF
         IERR = 0
C
Ccc Check stopping criterion
         UNRM = DNRM2 (N, U(1,K), 1)
C
C Compute v = [W.D^(-1).G.W^(-1)].u_k
         DO 210 I = 1, N
            WORK(LV+N-1+I) = U(I,K)/WT(I)/UNRM
  210    CONTINUE
         CALL MVDIFF (N, F, WORK(LV+N),
     +      ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +      TNP1, A0, DX, DY, DZ, RWORK,
     +      WORK(LHES), WORK(LV))
         IF (BDPREC) THEN
            CALL BCKBDI (N/NPDE, NPDE, PREC, WORK(LV), WORK(LV+N))
            DO 220 I = 1, N
               WORK(LV-1+I) = WT(I)*WORK(LV+N-1+I)*UNRM
  220       CONTINUE
         ELSE
            DO 221 I = 1, N
               WORK(LV-1+I) = WT(I)*PREC(I)*WORK(LV-1+I)*UNRM
  221       CONTINUE
         ENDIF
C
C
C Compute ZW[0:k-1,k] = C_k^T.v
         DO 300 I = 0, K-1
            ZW(I,K) = DDOT(N, C(1,I),1, WORK(LV),1)
  300    CONTINUE
C
         GRITER = GRITER + 1
C
Ccc Check stopping criterion
         RNRM = DNRM2 (N, R, 1)
         IF (LUN .NE. 0)  THEN
            WRITE(LUN,'(2I10,2E20.7)') GRITER, PRITER, RNRM, UNRM
         ENDIF
         IF (MAX(RNRM,UNRM) .LT. TOL) THEN
C Compute x = x + U_k.Z_k^(-1).1
            DO 310 I = K, 0, -1
               WORK(LV+I) = 1
               DO 320 J = I+1, K
                  WORK(LV+I) = WORK(LV+I) - ZW(I,J)*WORK(LV+J)
  320          CONTINUE
               CALL DAXPY (N, WORK(LV+I), U(1,I), 1, XV, 1)
  310       CONTINUE
            ITER = GRITER + PRITER
            GOTO 900
         ENDIF
  200 CONTINUE
Ccc End outer loop
C Compute x = x + U_k.Z_k^(-1).1
      K = MAXLR-1
      DO 330 I = K, 0, -1
         WORK(LV+I) = 1
         DO 340 J = I+1, K
            WORK(LV+I) = WORK(LV+I) - ZW(I,J)*WORK(LV+J)
  340    CONTINUE
         CALL DAXPY (N, WORK(LV+I), U(1,I), 1, XV, 1)
  330    CONTINUE
C
  150 CONTINUE
Ccc End Restart loop
C
      IERR = 1
      ITER = GRITER + PRITER
C
  900 CONTINUE
C Unscale x
      DO 910 I = 1, N
         XV(I) = XV(I) / WT(I)
  910 CONTINUE

      RETURN
      END
      SUBROUTINE GMRESO (N, XV, BV, WT, BDPREC, PREC, CO, M, TOL, MAXL,
     +   F, ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +   TNP1, A0, DX, DY, DZ, RWORK,
     +   V, HES, Q, WORK,
     +   LUN, ITER, ERR, IERR)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER  N, M, ISTRUC(0:*), NPDE, MAXL, LUN,
     +   ITER, IERR
      LOGICAL BDPREC
      DOUBLE PRECISION XV(N), BV(N), WT(N), PREC(N), CO(N,0:M), TOL,
     +   F(*), X(*), Y(*), Z(*), UIB(*), UNP1(0:*),
     +   TNP1, A0, DX, DY, DZ, RWORK(*),
     +   V(N,MAXL+1), HES(MAXL+1,MAXL),
     +   Q(2*MAXL), WORK(*), ERR
C
Ccc PURPOSE:
C Solve a Non-Symmetric linear system
C    [W.D^(-1).G.W^(-1)].[Wx] = [W.D^(-1).b]
C using the (block)-diagonally scaled GMRES(MAXL) method, orthogonalize
C new V_k not only against previous ones but also against C's from
C outer iteration.
C W = diag(WT(i)) and D is the (block) diagonal of G.
C The right hand-side W.D^(-1).b is stored in B,
C the matrix G and the preconditioner are stored in RWORK and IWORK.
C (Dx) is returned in X.
C The routine MVDIFF (N, RWORK, IWORK, X, Y) should perform y = Gx
C
Ccc PARAMETER DESCRIPTION:
C N      : IN. # grid points
C XV     : OUT. Final approximate solution.
C BV     : IN. Preconditioned right-hand side vector.
C          OUT. Residual vector.
C WT     : IN. Contains weight factors to compute weighted norm.
C BDPREC : IN. if true a block-diagonal preconditioner is used
C PREC   : IN. LU decomposition of (block-)diagonal of G. main diagonal
C          inverted
C CO     : IN. (.,0:M-1): vectors from outer iteration against which
C          V's should be orthogonalized.
C          OUT. (.,M) = Residual_outer_old - Residual
C M      : IN. Outer loop iteration count
C TOL    : IN. System is considered to be solved if
C          2-norm < TOL
C F      : IN. -I Parameters
C  ...          I    for
C RWORK  : IN. -I residual evaluation
C MAXL   : IN. max. iterations GMRES (no restarts)
C V      : WORK.
C HES    : WORK.
C Q      : WORK.
C WORK   : WORK. 4.N+1
C LUN    : IN. Logical unit # of file on which to write the error at
C          each iteration, if this is desired for monitoring convergence
C          If LUN = 0, no writing will occur.
C ITER   : OUT. Number of iterations required to reach convergence, or
C          MAXL+1 if convergence criterion could not be achieved in
C          MAXL iterations.
C ERR    : OUT. Weighted max. norm of error estimate in final
C          approximate solution
C IERR   : OUT. Error return flag
C          0: OK
C          1: Method failed to converge in MAXL iterations
C
Ccc EXTERNALS USED:
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL BCKBDI, MVDIFF, RCOPY, DAXPY, DDOT, DNRM2, ZERO
C
C-----------------------------------------------------------------------
C
      INTEGER I, J, K
      DOUBLE PRECISION C, CNM2, R0NRM, PROD, RHO, S, TEM, T1, T2, T,
     +   VNRM
C
      IERR = 0
      ITER = 0
C
Ccc Initialize solution on zero, the initial residual R0 is the
C left preconditioned vector B
      CALL ZERO (N, XV)
      CALL RCOPY (N, BV, V(1,1))
      R0NRM = DNRM2(N, V(1,1),1)
C
Ccc Check stopping criterion
      ERR = R0NRM
      IF (LUN .NE. 0)  THEN
         WRITE(LUN,'(''Diagonally scaled GMRESO(MAXL)'',I5)')
     +      MAXL
         WRITE(LUN,
     +      '('' ITER      Error Estimate'')')
         WRITE(LUN,'(I5,E20.7)') ITER, ERR
      ENDIF
C
Ccc Rescale so that the norm of V(1,1) is one
      DO 80 I = 1, N
         V(I,1) = V(I,1)/R0NRM
   80 CONTINUE
C
Ccc Initialize HES array.
      CALL ZERO (MAXL*(MAXL+1), HES)
C
Ccc Main loop to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
      PROD = 1.0
      DO 10 K = 1, MAXL
         ITER = K
C
C V(*,K+1) = [W.D^(-1).G.W^(-1)] . V(*,K)
         DO 11 I = 1, N
            WORK(I) = V(I,K)/WT(I)
   11    CONTINUE
         CALL MVDIFF (N, F, WORK,
     +      ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +      TNP1, A0, DX, DY, DZ, RWORK,
     +      WORK(1+N), V(1,K+1))
         IF (BDPREC) THEN
            CALL BCKBDI (N/NPDE, NPDE, PREC, V(1,K+1), WORK)
            DO 12 I = 1, N
               V(I,K+1) = WT(I)*WORK(I)
   12       CONTINUE
         ELSE
            DO 13 I = 1, N
               V(I,K+1) = WT(I)*PREC(I)*V(I,K+1)
   13       CONTINUE
         ENDIF
C
C Orthogonalize V(*,K+1) first against the previous C using
C modified Gram-Schmidt
         DO 801 I = 0, M-1
            TEM = DDOT (N, CO(1,I), 1, V(1,K+1), 1) /
     +            DDOT (N, CO(1,I), 1, CO(1,I), 1)
            CALL DAXPY (N, -TEM, CO(1,I), 1, V(1,K+1), 1)
  801    CONTINUE
C
C Orthogonalize V(*,K+1) against the previous V using
C modified Gram-Schmidt
         DO 81 I = 1, K
            HES(I,K) = DDOT (N, V(1,I), 1, V(1,K+1), 1)
            CALL DAXPY (N, -HES(I,K), V(1,I), 1, V(1,K+1), 1)
   81    CONTINUE
         VNRM = DNRM2(N, V(1,K+1), 1)
         HES(K+1,K) = VNRM
C
C Update the QR factors of HES (Q.HES = R) using Givens rotations
C First, multiply new column by previous Givens rotations
         DO 82 I = 1, K-1
            T1 = HES(I,K)
            T2 = HES(I+1,K)
            C = Q(2*I-1)
            S = Q(2*I)
            HES(I  ,K) = C*T1 - S*T2
            HES(I+1,K) = S*T1 + C*T2
   82    CONTINUE
C Form last Givens rotation and multiply it with last 2 elements of HES
         T1 = HES(K,K)
         T2 = HES(K+1,K)
         IF (T2 .EQ. 0.0) THEN
            C = 1.0
            S = 0.0
         ELSE IF (ABS(T2) .GE. ABS(T1)) THEN
            T = T1/T2
            S = -1.0/SQRT(1.0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0/SQRT(1.0+T*T)
            S = -C*T
         ENDIF
         Q(2*K-1) = C
         Q(2*K  ) = S
         HES(K,K) = C*T1 - S*T2
         IF (HES(K,K) .EQ. 0.0) THEN
            IERR = 2
            RETURN
         ENDIF
C
C Update RHO, the estimate of the norm of the residual R0-A*XL.
         PROD = PROD*Q(2*K)
         RHO  = ABS(PROD*R0NRM)
C
Ccc Check stopping criterion
         ERR = RHO
         IF (LUN .NE. 0)  THEN
            WRITE(LUN,'(I5,2E20.7)') ITER, ERR, ERR/R0NRM
         ENDIF
         IF (ERR/R0NRM .LT. 0.001 .AND. ERR .LT. TOL) GOTO 100
         IF (K .EQ. MAXL) GOTO 20
C
C Rescale so that the norm of V(1,K+1) is one.
         DO 83 I = 1, N
            V(I,K+1) = V(I,K+1)/VNRM
   83    CONTINUE
   10 CONTINUE
C
   20 CONTINUE
      IF (RHO .GT. R0NRM) THEN
         IERR = 2
         RETURN
      ELSE
         IERR = 1
      ENDIF
C
Ccc Compute the approximation XL to the solution.
C Min. ||beta.e1 - Hk+1k.y||_2
C X = X + Vk.y
 100  CONTINUE
      K = ITER
      WORK(1) = R0NRM
      DO 110 I = 2, K+1
         WORK(I) = 0.0
 110  CONTINUE
C Q.beta.e1
      DO 84 I = 1, K
         C = Q(2*I-1)
         S = Q(2*I)
         T1 = WORK(I)
         T2 = WORK(I+1)
         WORK(I  ) = C*T1 - S*T2
         WORK(I+1) = S*T1 + C*T2
   84 CONTINUE
C Solve R.y = Q.beta.e1
      DO 85 I = 1, K
         J = K+1-I
         WORK(J) = WORK(J) / HES(J,J)
         CALL DAXPY (J-1, -WORK(J), HES(1,J),1, WORK,1)
   85 CONTINUE
C
C X = X + Vk.y
      DO 120 I = 1,K
         CALL DAXPY(N, WORK(I), V(1,I), 1, XV, 1)
 120  CONTINUE
C
C Calculate the residual vector
      CALL RCOPY (N, V(1,1), WORK(K+1))
      DO 86 I = 1, K-1
         S = Q(2*I)
         C = Q(2*I-1)
         DO 87 J = 1, N
            WORK(K+J) = S*WORK(K+J) + C*V(J,I+1)
   87    CONTINUE
   86 CONTINUE
      I = K
         S = Q(2*I)
         C = Q(2*I-1)/VNRM
         DO 88 J = 1, N
            WORK(K+J) = S*WORK(K+J) + C*V(J,I+1)
   88    CONTINUE
      DO 89 J = 1, N
         WORK(K+J) = WORK(K+J)*R0NRM*PROD
   89 CONTINUE
C
C Compute c_m = (b - r) / <c_m,c_m>
      DO 130 J = 1, N
         CO(J,M) = BV(J) - WORK(K+J)
  130 CONTINUE
      CNM2 = 1 / DDOT (N, CO(1,M), 1, CO(1,M), 1)
      DO 140 J = 1, N
         CO(J,M) = CO(J,M) * CNM2
  140 CONTINUE
C
C Inner residual = outer residual
      CALL RCOPY (N, WORK(K+1), BV)
C
      RETURN
      END
      SUBROUTINE MVDIFF (N, F, XV,
     +   ISTRUC, X, Y, Z, NPDE, UIB, UNP1,
     +   TNP1, A0, DX, DY, DZ, RWORK,
     +   WORK, YV)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, ISTRUC(0:*), NPDE
      DOUBLE PRECISION F(N), XV(N),
     +   X(*), Y(*), Z(*), UIB(*), UNP1(0:*),
     +   TNP1, A0, DX, DY, DZ, RWORK(*),
     +   WORK(*), YV(N)
C
Ccc PURPOSE:
C Compute y = Gx where Gx ~ F(t,U+x,(U+x)_t)-F(t,U,(U+x)_t)
C
Ccc PARAMETER DESCRIPTION:
C N      : IN. Dimension of x
C F      : IN. Residual F(t,U,Udot), U=UNP1, Udot = A0.U+UH
C XV     : IN. Multiplying vector
C ISTRUC : IN. -I Parameters
C  ...          I    for
C RWORK  : IN. -I residual evaluation
C WORK   : WORK. (N+1 + 2N)
C YV     : OUT. Result vector
C
Ccc EXTERNALS USED:
      EXTERNAL RESID
C
C-----------------------------------------------------------------------
C
      INTEGER LLPLN, NPLNS, LIPLN, LLROW, NROWS, NPTS, LIROW, LICOL,
     +   LLLBND, NBNDS, NBDPTS, NBIPTS, LILBND, LLBNDP,
     +   LLBLWY, LLABVY, LLBLWZ, LLABVZ,
     +   LUT, LUX, LUY, LUZ, LUXX, LUYY, LUZZ, LUXY, LUXZ, LUYZ,
     +   LUBAR, LUTBAR, LFBAR, I
C
      LLPLN  = 0
      NPLNS  = ISTRUC(LLPLN)
      LIPLN  = LLPLN+NPLNS+2
      LLROW  = LIPLN+NPLNS
      NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
      NPTS   = ISTRUC(LLROW+NROWS)-1
      LIROW  = LLROW+NROWS+1
      LICOL  = LIROW+NROWS
      LLLBND = LICOL+NPTS
      NBNDS  = ISTRUC(LLLBND)
      NBDPTS = ISTRUC(LLLBND+NBNDS+1)-1
      NBIPTS = ISTRUC(LLLBND+NBNDS+2)-1
      LILBND = LLLBND+NBNDS+3
      LLBNDP = LILBND+NBNDS
      LLBLWY = LLBNDP+NBIPTS
      LLABVY = LLBLWY+NPTS+1
      LLBLWZ = LLABVY+NPTS+1
      LLABVZ = LLBLWZ+NPTS+1
C
      LUT    = 1
      LUX    = LUT  + N
      LUY    = LUX  + N
      LUZ    = LUY  + N
      LUXX   = LUZ  + N
      LUYY   = LUXX + N
      LUZZ   = LUYY + N
      LUXY   = LUZZ + N
      LUXZ   = LUXY + N
      LUYZ   = LUXZ + N
C
      LUBAR  = 1
      LUTBAR = LUBAR + 1+N
      LFBAR  = LUTBAR + N
C
Ccc Store U+x in WORK(LUBAR), and d(U+x)/dt in WORK(LUTBAR)
      WORK(LUBAR) = 0.0
      DO 10 I = 1, N
         WORK(LUBAR+I)    = UNP1(I) + XV(I)
         WORK(LUTBAR-1+I) = RWORK(LUT-1+I) + A0*XV(I)
   10 CONTINUE

C
Ccc Compute space derivatives and residual
      CALL DERIVS (NPTS, NPDE, WORK(LUBAR),
     +   ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP),
     +   ISTRUC(LLBLWY), ISTRUC(LLABVY), ISTRUC(LLBLWZ), ISTRUC(LLABVZ),
     +   DX, DY, DZ,
     +   RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ))
      CALL RES (TNP1, X, Y, Z, NPTS, NPDE, WORK(LUBAR+1),
     +   ISTRUC(LLLBND), ISTRUC(LILBND), ISTRUC(LLBNDP), UIB,
     +   WORK(LUTBAR),  RWORK(LUX),  RWORK(LUY), RWORK(LUZ),
     +   RWORK(LUXX), RWORK(LUYY), RWORK(LUZZ),
     +   RWORK(LUXY), RWORK(LUXZ), RWORK(LUYZ), WORK(LFBAR))
C
      DO 20 I = 1, N
         YV(I) = WORK(LFBAR-1+I) - F(I)
   20 CONTINUE

      RETURN
      END
      SUBROUTINE PRDOM (LPLN, IPLN, LROW, IROW, ICOL,
     +   LLBND, ILBND, LBND, IDOM, NX, NY, NZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
     +   LLBND(0:*), ILBND(*), LBND(*), IDOM(0:*), NX, NY, NZ
C
Ccc PURPOSE:
C Print domain plane-wise. Internal points are .., external points XX,
C physical plane-boundary points their ILBND value. Edges are given
C both ILBND values, corners an explicated 2-character value, and
C internal boundary values II.
C 
Ccc PARAMETER DESCRIPTION:
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C LLBND  : OUT. (0:LLBND(0)+2)
C          LLBND(0) = NBNDS: total # physical planes in actual domain.
C             NB. edges and corners are stored for each plane they
C             belong to.
C          LLBND(1:NBNDS): pointers to a specific boundary in LBND
C          LLBND(NBNDS+1) = NBDPTS+1: total # physical boundary points
C                                     in LBND + 1
C          LLBND(NBNDS+1): pointer to internal boundary in LBND
C          LLBND(NBNDS+2) = NBIPTS+1: total # points in LBND + 1
C ILBND  : OUT. (NBNDS)
C          ILBND(IB): type of boundary:
C                     1: Left  plane -I
C                     2: Down  plane  I
C                     3: Right plane  I max. first order derivative
C                     4: Up    plane  I
C                     5: Front plane  I
C                     6: Back  plane -I
C LBND   : (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C IDOM   : OUT. IDOM(IPT): location in domain of node IPT
C          0: interior point
C          IB: on boundary
C              if not on edge or corner, or if internal boundary
C              then IB, IB = 1, NBNDS+1
C              else | IB4 | IB3 | IB2 | IBYTE | if point is part of
C                   the boundary planes IBYTE, IB2, IB3, and IB4
C
Ccc EXTERNALS USED:
      EXTERNAL DOMFLG
C
      INTEGER BYTE1, BYTE2
      PARAMETER (BYTE1 = 2**8, BYTE2 = 2**16)
C
C-----------------------------------------------------------------------
C
      INTEGER MAXC, MAXN
      PARAMETER (MAXC = 100, MAXN = 25)
      CHARACTER*80 LINE(0:MAXN)
      INTEGER I, J, K, L, IB, IP, IR, IPT, NPLNS, NROWS, NPTS, NBNDS
      INTEGER IA, IE, IEVAL(0:MAXC)
C
      IF (NX .GE. MAXN .OR. NY .GE. MAXN) THEN
         PRINT *, 'Sorry, Nx, Ny should be < MAXN; adapt PRDOM'
         RETURN
      ENDIF
C
      NPLNS = LPLN(0)
      NROWS = LPLN(NPLNS+1)-1
      NPTS  = LROW(NROWS+1)-1
      NBNDS = LLBND(0)
C
C Set domain values
      CALL DOMFLG (NPTS, LLBND, LBND, IDOM)

      IA = ICHAR('a')
      IE = -1
      DO 10 IP = 1, NPLNS
         K = IPLN(IP)
         LINE(0) = ' '
         DO 20 I = 1, NX+1
            WRITE(LINE(0)(3*I-2:3*I),'('' XX'')')
   20    CONTINUE
         DO 30 J = 1, NY
            LINE(J) = LINE(0)
   30    CONTINUE
         DO 40 IR = LPLN(IP), LPLN(IP+1)-1
            J = IROW(IR)
            DO 50 IPT = LROW(IR), LROW(IR+1)-1
               I = ICOL(IPT)+1
               IF (IDOM(IPT) .EQ. 0) THEN
                  WRITE(LINE(J)(3*I-2:3*I),'('' ..'')')
               ELSE IF (IDOM(IPT) .EQ. NBNDS+1) THEN
                  WRITE(LINE(J)(3*I-2:3*I),'('' II'')')
               ELSE
                  IB = IDOM(IPT)
                  IF (IB .LT. BYTE1) THEN
                     WRITE(LINE(J)(3*I-2:3*I),'(I3)') ILBND(IB)
                  ELSE IF (IB .LT. BYTE2) THEN
                     WRITE(LINE(J)(3*I-1:3*I),'(2I1)')
     +                  ILBND(IB/BYTE1), ILBND(MOD(IB,BYTE1))
                  ELSE
                     DO 55 L = 0, IE
                        IF (IB .EQ. IEVAL(L)) GOTO 56
   55                CONTINUE
                     IE = IE+1
                     IF (IE .GT. MAXC) THEN
                        PRINT *, 'Sorry, # corners > MAXC; adapt PRDOM'
                        RETURN
                     ENDIF
                     IEVAL(IE) = IB
                     L = IE
   56                WRITE(LINE(J)(3*I-1:3*I),'(2A)')
     +                  CHAR(MOD(L,26)+IA), CHAR(L/26+IA)
                  ENDIF
               ENDIF
   50       CONTINUE
   40    CONTINUE
         PRINT '(/,''Plane:'',I3)', K
         DO 60 J = NY, 0, -1
            PRINT '(A)', LINE(J)
            PRINT *
   60    CONTINUE
   10 CONTINUE

      IF (IE .EQ. -1) RETURN
      PRINT '(///,''Legenda corners:'')'
      DO 100 L = 0, IE
         LINE(0) = ' '
         IB = MOD(IEVAL(L),BYTE1)
         WRITE(LINE(0)(1:2),'(I2.1)') ILBND(IB)
         IEVAL(L) = IEVAL(L)/BYTE1
         IB = MOD(IEVAL(L),BYTE1)
         WRITE(LINE(0)(3:4),'(I2.1)') ILBND(IB)
         IF (IB .EQ. IEVAL(L)) GOTO 110
         IEVAL(L) = IEVAL(L)/BYTE1
         IB = MOD(IEVAL(L),BYTE1)
         WRITE(LINE(0)(5:6),'(I2.1)') ILBND(IB)
         IF (IB .EQ. IEVAL(L)) GOTO 110
         IEVAL(L) = IEVAL(L)/BYTE1
         IB = MOD(IEVAL(L),BYTE1)
         WRITE(LINE(0)(7:8),'(I2.1)') ILBND(IB)
  110    PRINT '(X,2A,'':'',A8)', CHAR(MOD(L,26)+IA), CHAR(L/26+IA),
     +      LINE(0)(1:8)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SETXYZ (XL, YF, ZD, DX, DY, DZ,
     +   LPLN, IPLN, LROW, IROW, ICOL, X, Y, Z)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*)
      DOUBLE PRECISION XL, YF, ZD, DX, DY, DZ, X(*), Y(*), Z(*)
C
Ccc PURPOSE:
C Store X-, Y- and Z-coordinates of the grid points.
C
Ccc PARAMETER DESCRIPTION:
C XL     : IN. X-coordinate of left/front/down point of virtual box
C YF     : IN. Y-coordinate of left/front/down point of virtual box
C ZD     : IN. Z-coordinate of left/front/down point of virtual box
C DX     : IN. Grid width in X-direction
C DY     : IN. Grid width in Y-direction
C DZ     : IN. Grid width in Z-direction
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C X      : OUT. Contains the X-coordinates for the grid
C Y      : OUT. Contains the Y-coordinates for the grid
C Z      : OUT. Contains the Z-coordinates for the grid
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IP, IPT, IR, NPLNS
      DOUBLE PRECISION YIR, ZIP
C
      NPLNS = LPLN(0)
      DO 10 IP = 1, NPLNS
         ZIP = ZD + IPLN(IP)*DZ
         DO 20 IR = LPLN(IP), LPLN(IP+1)-1
            YIR = YF + IROW(IR)*DY
            DO 30 IPT = LROW(IR), LROW(IR+1)-1
               X(IPT) = XL + ICOL(IPT)*DX
               Y(IPT) = YIR
               Z(IPT) = ZIP
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      
      RETURN
      END
      SUBROUTINE PRSOL (LUN, T, NPDE, XL, YF, ZD, DXB, DYB, DZB,
     +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUN, NPDE, LGRID(0:*), ISTRUC(*), LSOL(*)
      DOUBLE PRECISION T, XL, YF, ZD, DXB, DYB, DZB, SOL(*)
C
Ccc PURPOSE:
C Print solution and coordinate values at all grid levels.
C
Ccc PARAMETER DESCRIPTION:
C LUN    : IN.  Logical unit number of print file
C T      : IN.  Current value of time variable
C NPDE   : IN.  # PDE components
C XL     : IN. X-coordinate of left/front/down point of virtual box
C YF     : IN. Y-coordinate of left/front/down point of virtual box
C ZD     : IN. Z-coordinate of left/front/down point of virtual box
C DXB    : IN.  Cell width in X-direction of base grid
C DYB    : IN.  Cell width in Y-direction of base grid
C DZB    : IN.  Cell width in Z-direction of base grid
C LGRID  : IN.  (0:*)
C          LGRID(0) = max. grid level used at T
C          LGRID(1): pointer to base grid structure ISTRUC
C          LGRID(LEVEL): pointer to grid structure
C                        (LPLN,IPLN,LROW,IROW,ICOL)
C                        of refinement level LEVEL for time T
C ISTRUC : IN.  (*)
C          ISTRUC(LGRID(LEVEL):.) contains (LPLN,IPLN,LROW,IROW,ICOL)
C                                 of grid level LEVEL,
C          LPLN   : (0:LPLN(0)+1)
C             LPLN(0) = NPLNS: Actual # planes in grid
C             LPLN(1:NPLNS): pointers to the start of a plane in LROW
C             LPLN(NPLNS+1) = NROWS+1: Actual # rows in grid + 1
C          IPLN   : (NPLNS)
C             IPLN(IP): plane number of plane IP in virtual box
C          LROW   : (NROWS+1)
C             LROW(1:NROWS): pointers to the start of a row in the grid
C             LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C          IROW   : (NROWS)
C             IROW(IR): row number of row IR in virtual box
C          ICOL   : (NPTS)
C             ICOL(IPT): column number of grid point IPT in virtual box
C LSOL   : IN.  (*)
C          LSOL(LEVEL): pointer to (injected) solution at grid
C                       of refinement level LEVEL for time T
C SOL    : IN.  (*)
C          SOL(LSOL(LEVEL)+1:LSOL(LEVEL)+NPTS(LEVEL)*NPDE) contains
C          U_LEVEL(NPTS,NPDE)
C
Ccc EXTERNALS USED:
      EXTERNAL PRSOLL
C
C-----------------------------------------------------------------------
C
      INTEGER MAXLEV, LEVEL, LLPLN, LIPLN, LLROW, LIROW, LICOL,
     +   NPLNS, NROWS, NPTS
      DOUBLE PRECISION DX, DY, DZ

      MAXLEV = LGRID(0)
      DX = DXB
      DY = DYB
      DZ = DZB
      DO 10 LEVEL = 1, MAXLEV
         LLPLN  = LGRID(LEVEL)
         NPLNS  = ISTRUC(LLPLN)
         NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
         LIPLN  = LLPLN+NPLNS+2
         LLROW  = LIPLN+NPLNS
         NPTS   = ISTRUC(LLROW+NROWS)-1
         LIROW  = LLROW+NROWS+1
         LICOL  = LIROW+NROWS
         CALL PRSOLL (LUN, LEVEL, T, NPTS, NPDE, XL, YF, ZD, DX, DY, DZ,
     +      ISTRUC(LLPLN), ISTRUC(LIPLN), ISTRUC(LLROW), ISTRUC(LIROW),
     +      ISTRUC(LICOL), SOL(LSOL(LEVEL)+1))
         DX = DX/2
         DY = DY/2
         DZ = DZ/2
   10 CONTINUE
      RETURN
      END
      SUBROUTINE PRSOLL (LUN, LEVEL, T, NPTS, NPDE, XL, YF, ZD,
     +   DX, DY, DZ, LPLN, IPLN, LROW, IROW, ICOL, U)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUN, LEVEL, NPTS, NPDE, LPLN(0:*), IPLN(*),
     +   LROW(*), IROW(*), ICOL(*)
      DOUBLE PRECISION T, XL, YF, ZD, DX, DY, DZ, U(NPTS,NPDE)
C
Ccc PURPOSE:
C Print solution and  X-, Y- and Z-coordinates of gridlevel LEVEL.
C
Ccc PARAMETER DESCRIPTION:
C LUN    : IN.  Logical unit number of print file
C LEVEL  : IN.  Grid level corresponding with solution U.
C T      : IN.  Current value of time variable
C NPTS   : IN.  # grid points at this level
C NPDE   : IN.  # PDE components
C XL     : IN. X-coordinate of left/front/down point of virtual box
C YF     : IN. Y-coordinate of left/front/down point of virtual box
C ZD     : IN. Z-coordinate of left/front/down point of virtual box
C DX     : IN. Grid width in X-direction
C DY     : IN. Grid width in Y-direction
C DZ     : IN. Grid width in Z-direction
C LPLN   : IN. (0:LPLN(0)+1)
C          LPLN(0) = NPLNS: Actual # planes in LROW
C          LPLN(1:NPLNS): pointers to the start of a plane in LROW
C          LPLN(NPLNS+1) = NROWS+1: Total # rows in grid + 1
C IPLN   : IN. (NPLNS)
C          IPLN(IP): plane number of plane IP in virtual box
C LROW   : IN. (NROWS+1)
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN. (NROWS)
C          IROW(IR): row number of row IR in virtual box
C ICOL   : IN. (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual box
C U      : IN.  Solution on this grid level
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, IP, IPT, IR, NPLNS
      DOUBLE PRECISION X, Y, Z
C
      NPLNS = LPLN(0)

      WRITE(LUN,'(//// A,T14,A,T30,A,T46,A,T62,A,T71,A //)')
     +   'Lev', 't', 'Z', 'Y', 'X', 'Solution'
      IP = 1
         Z = ZD + IPLN(IP)*DZ
         IR = LPLN(IP)
            Y = YF + IROW(IR)*DY
            IPT = LROW(IR)
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(I3,T5,E12.5,T21,E12.5,T37,E12.5,T53,E12.5,T69,E12.5)')
     +         LEVEL, T, Z, Y, X, U(IPT,1)
               DO 10 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   10          CONTINUE
               DO 14 IPT = LROW(IR)+1, LROW(IR+1)-1
                  X = XL + ICOL(IPT)*DX
                  WRITE(LUN,
     +            '(T53,E12.5,T69,E12.5)')
     +            X, U(IPT,1)
                  DO 15 IC = 2, NPDE
                     WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   15             CONTINUE
   14          CONTINUE
         DO 20 IR = LPLN(IP)+1, LPLN(IP+1)-1
            Y = YF + IROW(IR)*DY
            IPT = LROW(IR)
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(T37,E12.5,T53,E12.5,T69,E12.5)')
     +         Y, X, U(IPT,1)
               DO 30 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   30          CONTINUE
            DO 40 IPT = LROW(IR)+1, LROW(IR+1)-1
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(T53,E12.5,T69,E12.5)')
     +         X, U(IPT,1)
               DO 50 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   50          CONTINUE
   40       CONTINUE
   20    CONTINUE
      DO 100 IP = 2, NPLNS
         Z = ZD + IPLN(IP)*DZ
         IR = LPLN(IP)
            Y = YF + IROW(IR)*DY
            IPT = LROW(IR)
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(T21,E12.5,T37,E12.5,T53,E12.5,T69,E12.5)')
     +         Z, Y, X, U(IPT,1)
               DO 110 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
  110          CONTINUE
               DO 114 IPT = LROW(IR)+1, LROW(IR+1)-1
                  X = XL + ICOL(IPT)*DX
                  WRITE(LUN,
     +            '(T53,E12.5,T69,E12.5)')
     +            X, U(IPT,1)
                  DO 115 IC = 2, NPDE
                     WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
  115             CONTINUE
  114          CONTINUE
         DO 120 IR = LPLN(IP)+1, LPLN(IP+1)-1
            Y = YF + IROW(IR)*DY
            IPT = LROW(IR)
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(T37,E12.5,T53,E12.5,T69,E12.5)')
     +         Y, X, U(IPT,1)
               DO 130 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
  130          CONTINUE
            DO 140 IPT = LROW(IR)+1, LROW(IR+1)-1
               X = XL + ICOL(IPT)*DX
               WRITE(LUN,
     +         '(T53,E12.5,T69,E12.5)')
     +         X, U(IPT,1)
               DO 150 IC = 2, NPDE
                  WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
  150          CONTINUE
  140       CONTINUE
  120    CONTINUE
  100 CONTINUE

      RETURN
      END
      SUBROUTINE WRUNI (LUNS, LUNG, UNILEV,
     +   T, NPDE, XL, YF, ZD, DXB, DYB, DZB, NXB, NYB, NZB,
     +   LGRID, ISTRUC, LSOL, SOL, UNIFRM, NX, NY, NZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUNS, LUNG, UNILEV,
     +   NPDE, NXB, NYB, NZB, LGRID(0:*), ISTRUC(*), LSOL(*), NX, NY, NZ
      DOUBLE PRECISION T, XL, YF, ZD, DXB, DYB, DZB, SOL(*),
     +   UNIFRM(0:NX,0:NY,0:NZ,NPDE)
C
Ccc PURPOSE:
C Write (interpolated) solution values at grid level UNILEV to file
C LUNS.
C Write maximum gridlevel used in each point to file LUNG.
C NB. The data will not be correct for a domain with holes in it with
C a size of the width of the base grid.
C
Ccc PARAMETER DESCRIPTION:
C LUNS   : IN.  Logical unit number of solution file
C LUNG   : IN.  Logical unit number of grid level file
C UNILEV : IN.  Maximum grid level to be used to generate the data
C T      : IN.  Value of time variable
C NPDE   : IN.  # PDE components
C XL     : IN. X-coordinate of left/front/down point of virtual box
C YF     : IN. Y-coordinate of left/front/down point of virtual box
C ZD     : IN. Z-coordinate of left/front/down point of virtual box
C DXB    : IN.  Cell width in X-direction of base grid
C DYB    : IN.  Cell width in Y-direction of base grid
C DZB    : IN.  Cell width in Z-direction of base grid
C NXB,NYB,NZB: IN. # gridcells in X-, Y- and Z-direction, resp., on grid
C          of base level
C LGRID  : IN.  (0:*)
C          LGRID(0) = max. grid level used at T
C          LGRID(1): pointer to base grid structure ISTRUC
C          LGRID(LEVEL): pointer to grid structure
C                        (LPLN,IPLN,LROW,IROW,ICOL)
C                        of refinement level LEVEL for time T
C ISTRUC : IN.  (*)
C          ISTRUC(LGRID(LEVEL):.) contains (LPLN,IPLN,LROW,IROW,ICOL)
C                                 of grid level LEVEL,
C          LPLN   : (0:LPLN(0)+1)
C             LPLN(0) = NPLNS: Actual # planes in grid
C             LPLN(1:NPLNS): pointers to the start of a plane in LROW
C             LPLN(NPLNS+1) = NROWS+1: Actual # rows in grid + 1
C          IPLN   : (NPLNS)
C             IPLN(IP): plane number of plane IP in virtual box
C          LROW   : (NROWS+1)
C             LROW(1:NROWS): pointers to the start of a row in the grid
C             LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C          IROW   : (NROWS)
C             IROW(IR): row number of row IR in virtual box
C          ICOL   : (NPTS)
C             ICOL(IPT): column number of grid point IPT in virtual box
C LSOL   : IN.  (*)
C          LSOL(LEVEL): pointer to (injected) solution at grid
C                       of refinement level LEVEL for time T
C SOL    : IN.  (*)
C          SOL(LSOL(LEVEL)+1:LSOL(LEVEL)+NPTS(LEVEL)*NPDE) contains
C          U_LEVEL(NPTS,NPDE)
C UNIFRM : WORK. (Interpolated) solution on level UNILEV / max. grid
C          level used.
C NX,NY,NZ: IN. # gridcells in X-, Y- and Z-direction, resp., on grid
C          of level UNILEV
C
C-----------------------------------------------------------------------
C
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!                                                                !!!
C !!! In subroutine WRUNI the constant NONVAL should be adjusted to  !!!
C !!! the data (NONVAL = impossible value for the first componenent) !!!
C !!!                                                                !!!
      DOUBLE PRECISION NONVAL
      PARAMETER (NONVAL = -999.999)
C !!!                                                                !!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C-----------------------------------------------------------------------
C

      INTEGER I, IC, ICOL, IMUL, IP, IPLN, IPT, IR, IROW, J, K,
     +   LEVEL, LLPLN, LIPLN, LLROW, LIROW, LICOL, MAXLEV,
     +   NPLNS, NROWS, NPTS

      DO 1 IC = 1, NPDE
      DO 1 IPLN = 0, NZ
      DO 1 IROW = 0, NY
      DO 1 ICOL = 0, NX
         UNIFRM(ICOL,IROW,IPLN,IC) = NONVAL
    1 CONTINUE

      MAXLEV = LGRID(0)
      DO 10 LEVEL = 1, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         LLPLN  = LGRID(LEVEL)
         NPLNS  = ISTRUC(LLPLN)
         NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
         LIPLN  = LLPLN+NPLNS+2
         LLROW  = LIPLN+NPLNS
         NPTS   = ISTRUC(LLROW+NROWS)-1
         LIROW  = LLROW+NROWS+1
         LICOL  = LIROW+NROWS
         DO 20 IP = 1, NPLNS
            IPLN = ISTRUC(LIPLN-1+IP)*IMUL
            DO 30 IR = ISTRUC(LLPLN+IP), ISTRUC(LLPLN+IP+1)-1
               IROW = ISTRUC(LIROW-1+IR)*IMUL
               DO 40 IPT = ISTRUC(LLROW-1+IR), ISTRUC(LLROW+IR)-1
                  ICOL = ISTRUC(LICOL-1+IPT)*IMUL
                  DO 50 IC = 1, NPDE
                     UNIFRM(ICOL,IROW,IPLN,IC) =
     +                  SOL(LSOL(LEVEL)+(IC-1)*NPTS+IPT)
   50             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      DO 100 LEVEL = 2, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         DO 110 K = IMUL, NZ, IMUL*2
         DO 110 J = 0, NY, IMUL*2
         DO 110 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .EQ. NONVAL) THEN
               DO 120 IC = 1, NPDE
                  UNIFRM(I,J,K,IC) = 
     +               (UNIFRM(I,J,K-IMUL,IC)+UNIFRM(I,J,K+IMUL,IC))/2
  120          CONTINUE
            ENDIF
  110    CONTINUE
         DO 130 K = 0, NZ, IMUL
         DO 130 J = IMUL, NY, IMUL*2
         DO 130 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .EQ. NONVAL) THEN
               DO 140 IC = 1, NPDE
                  UNIFRM(I,J,K,IC) = 
     +               (UNIFRM(I,J-IMUL,K,IC)+UNIFRM(I,J+IMUL,K,IC))/2
  140          CONTINUE
            ENDIF
  130    CONTINUE
         DO 150 K = 0, NZ, IMUL
         DO 150 J = 0, NY, IMUL
         DO 150 I = IMUL, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .EQ. NONVAL) THEN
               DO 160 IC = 1, NPDE
                  UNIFRM(I,J,K,IC) = 
     +               (UNIFRM(I-IMUL,J,K,IC)+UNIFRM(I+IMUL,J,K,IC))/2
  160          CONTINUE
            ENDIF
  150    CONTINUE
  100 CONTINUE

      DO 170 K = 0, NZ
      DO 170 J = 0, NY
      DO 170 I = 0, NX
         WRITE(LUNS,'(100E13.3)') (UNIFRM(I,J,K,IC), IC = 1, NPDE)
  170 CONTINUE
C
C Grids
      DO 201 IPLN = 0, NZ
      DO 201 IROW = 0, NY
      DO 201 ICOL = 0, NX
         UNIFRM(ICOL,IROW,IPLN,1) = 0
  201 CONTINUE
      DO 210 LEVEL = 1, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         LLPLN  = LGRID(LEVEL)
         NPLNS  = ISTRUC(LLPLN)
         NROWS  = ISTRUC(LLPLN+NPLNS+1)-1
         LIPLN  = LLPLN+NPLNS+2
         LLROW  = LIPLN+NPLNS
         NPTS   = ISTRUC(LLROW+NROWS)-1
         LIROW  = LLROW+NROWS+1
         LICOL  = LIROW+NROWS
         DO 220 IP = 1, NPLNS
            IPLN = ISTRUC(LIPLN-1+IP)*IMUL
            DO 230 IR = ISTRUC(LLPLN+IP), ISTRUC(LLPLN+IP+1)-1
               IROW = ISTRUC(LIROW-1+IR)*IMUL
               DO 240 IPT = ISTRUC(LLROW-1+IR), ISTRUC(LLROW+IR)-1
                  ICOL = ISTRUC(LICOL-1+IPT)*IMUL
                  UNIFRM(ICOL,IROW,IPLN,1) = LEVEL
  240          CONTINUE
  230       CONTINUE
  220    CONTINUE
  210 CONTINUE
      DO 300 LEVEL = 2, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         DO 310 K = IMUL, NZ, IMUL*2
         DO 310 J = 0, NY, IMUL*2
         DO 310 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .LT. LEVEL) THEN
               UNIFRM(I,J,K,1) =
     +            MIN(UNIFRM(I,J,K-IMUL,1),UNIFRM(I,J,K+IMUL,1))
            ENDIF
  310    CONTINUE
         DO 320 K = 0, NZ, IMUL
         DO 320 J = IMUL, NY, IMUL*2
         DO 320 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .LT. LEVEL) THEN
               UNIFRM(I,J,K,1) =
     +            MIN(UNIFRM(I,J-IMUL,K,1),UNIFRM(I,J+IMUL,K,1))
            ENDIF
  320    CONTINUE
         DO 330 K = 0, NZ, IMUL
         DO 330 J = 0, NY, IMUL
         DO 330 I = IMUL, NX, IMUL*2
            IF (UNIFRM(I,J,K,1) .LT. LEVEL) THEN
               UNIFRM(I,J,K,1) =
     +            MIN(UNIFRM(I-IMUL,J,K,1),UNIFRM(I+IMUL,J,K,1))
            ENDIF
  330    CONTINUE
  300 CONTINUE

      DO 350 K = 0, NZ
      DO 350 J = 0, NY
      DO 350 I = 0, NX
         WRITE(LUNG,'(I2)') NINT(UNIFRM(I,J,K,1))
  350 CONTINUE
      RETURN
      END
      SUBROUTINE DUMP (LUNDMP, RWK, IWK)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUNDMP, IWK(*)
      DOUBLE PRECISION RWK(*)
C
Ccc PURPOSE:
C Dump all information necessary for a restart of VLUGR3 on file
C
Ccc PARAMETER DESCRIPTION:
C LUNDMP : IN.  Logical unit number of dumpfile. Should be opened as an
C          unformatted file.
C RWK    : IN.  Real workstorage as returned from VLUGR3
C IWK    : IN.  Integer workstorage as returned from VLUGR3
C
Ccc EXTERNALS USED: NONE
C
C
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
C
C
Ccc   INCLUDE 'CMNWRITEF'
C
C CMNWRITEF
C
C COMMON needed for continuation calls
      INTEGER MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      LOGICAL FIRST, SECOND
      DOUBLE PRECISION T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,
     +   DZB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW,
     +   DXB,DYB,DZB, DTO
      SAVE /WRITIF/, /WRITLF/, /WRITRF/
C
C end INCLUDE 'CMNWRITEF'
C
C
C-----------------------------------------------------------------------
C
      INTEGER I, J

      WRITE(LUNDMP) MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB,
     +   FIRST, SECOND,
     +   T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,DZB, DTO
      WRITE(LUNDMP) LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   (NJACS(I), I=1,MXCLEV), (NRESID(I), I=1,MXCLEV),
     +   (NNIT(I), I=1,MXCLEV), ((NLSIT(I,J), I=1,MXCLEV), J=1,MXCNIT)
      WRITE(LUNDMP) (RWK(I), I=1,LRWKPS+LRWKB)
      WRITE(LUNDMP) (IWK(I), I=1,LIWKPS+LIWKB)

      RETURN
      END
      SUBROUTINE RDDUMP (LUNDMP, RWK, LENRWK, IWK, LENIWK)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LENIWK
      INTEGER LUNDMP, LENRWK, IWK(LENIWK)
      DOUBLE PRECISION RWK(LENRWK)
C
Ccc PURPOSE:
C Read all information necessary for a restart of VLUGR3 from file
C
Ccc PARAMETER DESCRIPTION:
C LUNDMP : IN.  Logical unit number of dumpfile. Should be opened as an
C          unformatted file.
C RWK    : OUT. Real workstorage intended to pass to VLUGR3
C LENRWK : IN.  Dimension of RWK.
C IWK    : OUT. Integer workstorage intended to pass to VLUGR3
C LENIWK : IN.  Dimension of IWK.
C
Ccc EXTERNALS USED: NONE
C
C
Ccc   INCLUDE 'CMNSTATS'
C
C CMNSTATS
C
C COMMON with integration statistics
      INTEGER MXCLEV, MXCNIT
      PARAMETER (MXCLEV = 10, MXCNIT = 20)
      INTEGER LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS(MXCLEV), NRESID(MXCLEV), NNIT(MXCLEV),
     +   NLSIT(MXCLEV,MXCNIT)
      COMMON /STATS/ LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   NJACS, NRESID, NNIT, NLSIT
      SAVE /STATS/
C
C end INCLUDE 'CMNSTATS'
C
C
Ccc   INCLUDE 'CMNWRITEF'
C
C CMNWRITEF
C
C COMMON needed for continuation calls
      INTEGER MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      LOGICAL FIRST, SECOND
      DOUBLE PRECISION T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,
     +   DZB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW,
     +   DXB,DYB,DZB, DTO
      SAVE /WRITIF/, /WRITLF/, /WRITRF/
C
C end INCLUDE 'CMNWRITEF'
C
C
C-----------------------------------------------------------------------
C
      INTEGER I, J

      READ(LUNDMP) MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB,
     +   FIRST, SECOND,
     +   T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,DZB, DTO
      IF (LENRWK .LT. LRWKPS+LRWKB .OR. LENIWK .LT. LIWKPS+LIWKB) THEN
         PRINT *, LENRWK, LRWKPS+LRWKB, LENIWK, LIWKPS+LIWKB
         STOP 'work space too small'
      ENDIF
      READ(LUNDMP) LUNPDS, LUNNLS, LUNLSS, LEVEL, NSTEPS, NREJS,
     +   (NJACS(I), I=1,MXCLEV), (NRESID(I), I=1,MXCLEV),
     +   (NNIT(I), I=1,MXCLEV), ((NLSIT(I,J), I=1,MXCLEV), J=1,MXCNIT)
      READ(LUNDMP) (RWK(I), I=1,LRWKPS+LRWKB)
      READ(LUNDMP) (IWK(I), I=1,LIWKPS+LIWKB)
C
      RETURN
      END
      LOGICAL FUNCTION CHKWRK (LRWKN, LENRWK, LIWKN, LENIWK,
     +   LLWKN, LENLWK)
C-----------------------------------------------------------------------
      INTEGER LRWKN, LENRWK, LIWKN, LENIWK, LLWKN, LENLWK
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
      CHKWRK = .TRUE.
      IF (LRWKN .GT. LENRWK) THEN
         WRITE(LUNERR,
     +      '(''Real workspace too small, required at least:'',I10)')
     +      LRWKN
         CHKWRK = .FALSE.
      ENDIF
      IF (LIWKN .GT. LENIWK) THEN
         WRITE(LUNERR,
     +      '(''Integer workspace too small, required at least:'',I10)')
     +      LIWKN
         CHKWRK = .FALSE.
      ENDIF
      IF (LLWKN .GT. LENLWK) THEN
         WRITE(LUNERR,
     +      '(''Logical workspace too small, required at least:'',I10)')
     +      LLWKN
         CHKWRK = .FALSE.
      ENDIF

      RETURN
      END
      SUBROUTINE ERRWGT (NPTS, NPDE, U, RELTOL, ABSTOL, WT)
C-----------------------------------------------------------------------
      INTEGER NPTS, NPDE
      DOUBLE PRECISION U(NPTS,NPDE), RELTOL(NPDE), ABSTOL(NPDE),
     +   WT(NPTS,NPDE)
C-----------------------------------------------------------------------

      INTEGER IC, IPT
      DOUBLE PRECISION SN

      SN = 1.0/SQRT(DBLE(NPTS*NPDE))
      DO 10 IC = 1, NPDE
         DO 20 IPT = 1, NPTS
            WT(IPT,IC) = SN/(RELTOL(IC)*ABS(U(IPT,IC)) + ABSTOL(IC))
   20    CONTINUE
   10 CONTINUE
      
      RETURN
      END
      DOUBLE PRECISION FUNCTION MAXNRM (N, V)
C-----------------------------------------------------------------------
      INTEGER N
      DOUBLE PRECISION V(N)
C-----------------------------------------------------------------------

      INTEGER I

      MAXNRM = 0.0
      DO 10 I = 1, N
         MAXNRM = MAX(MAXNRM, ABS(V(I)))
   10 CONTINUE

      RETURN
      END
      DOUBLE PRECISION FUNCTION WMXNRM (N, V, W)
C-----------------------------------------------------------------------
      INTEGER N
      DOUBLE PRECISION V(N), W(N)
C-----------------------------------------------------------------------

      INTEGER I

      WMXNRM = 0.0
      DO 10 I = 1, N
         WMXNRM = MAX(WMXNRM, ABS(V(I)*W(I)))
   10 CONTINUE
      WMXNRM = WMXNRM*SQRT(DBLE(N))

      RETURN
      END
      DOUBLE PRECISION FUNCTION WDNRM2 (N, V, W)
C-----------------------------------------------------------------------
      INTEGER N
      DOUBLE PRECISION V(N), W(N)
C-----------------------------------------------------------------------

      INTEGER I

      WDNRM2 = 0.0
      DO 10 I = 1, N
         WDNRM2 = WDNRM2 + (V(I)*W(I))**2
   10 CONTINUE
      WDNRM2 = SQRT(WDNRM2)

      RETURN
      END
      SUBROUTINE ICOPY (LEN, A, B)
C-----------------------------------------------------------------------
      INTEGER LEN
      INTEGER A(LEN), B(LEN)
C-----------------------------------------------------------------------

      INTEGER I

      DO 10 I = 1, LEN
         B(I) = A(I)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE IYPOC (LEN, A, B)
C-----------------------------------------------------------------------
      INTEGER LEN
      INTEGER A(LEN), B(LEN)
C-----------------------------------------------------------------------

      INTEGER I

      DO 10 I = LEN, 1, -1
         B(I) = A(I)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE RCOPY (LEN, A, B)
C-----------------------------------------------------------------------
      INTEGER LEN
      DOUBLE PRECISION A(LEN), B(LEN)
C-----------------------------------------------------------------------

      INTEGER I

      DO 10 I = 1, LEN
         B(I) = A(I)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE ZERO (LEN, A)
C-----------------------------------------------------------------------
      INTEGER LEN
      DOUBLE PRECISION A(LEN)
C-----------------------------------------------------------------------

      INTEGER I

      DO 10 I = 1, LEN
         A(I) = 0.0
   10 CONTINUE

      RETURN
      END
      SUBROUTINE MACNUM
C-----------------------------------------------------------------------
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
Ccc EXTERNALS USED:
      INTEGER I1MACH
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH, I1MACH
C-----------------------------------------------------------------------
C
      LUNOUT = I1MACH(2)
      LUNERR = I1MACH(4)
      UROUND = D1MACH(4)
      XMIN   = D1MACH(1)

      RETURN
      END
