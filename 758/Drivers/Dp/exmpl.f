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
      PARAMETER (MXLEV=2, NPD=2, NPTS=2500)
      PARAMETER (LENIWK=NPTS*(5*MXLEV+11),
     +           LENRWK=NPTS*NPD*(5*MXLEV+9 +
     +                            (9*NPD+2*MAXLR+MAXL+7)),
     +           LENLWK=2*NPTS)
C
C-----------------------------------------------------------------------
C
      INTEGER LUNDMP
      PARAMETER (LUNDMP = 89)
C
      INTEGER NPDE, INFO(7), IWK(LENIWK), MNTR
      LOGICAL LWK(LENLWK)
      DOUBLE PRECISION T, TOUT, DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, RINFO(2+3*NPD), RWK(LENRWK)

C First call of VLUGR2
      MNTR = 0
      NPDE = 2
      T    = 0.0
      TOUT = 1.0
      DT   = 0.001
C Since domain is not a rectangle the grid parameters need not to be
C specified here (cf. INIDOM)
      TOLS = 0.1
      TOLT = 0.05
      INFO(1) = 1
C MAXLEV
      INFO(2) = 5
C Domain not a rectangle
      INFO(3) = 1
C Linear system solver: GCRO + Diagonal scaling
C (no first order derivatives at the boundaries)
      INFO(4) = 13
      OPEN (UNIT=61,FILE='RunInfo')
C Write integration history to unit # 61
      INFO(5) = 61
C Write Newton info to unit # 61
      INFO(6) = 61
C Write GCRO info to unit # 61
      INFO(7) = 61
C DTMIN = 1D-7
      RINFO(1) = 1.0D-7
C DTMAX = 1.0
      RINFO(2) = 1.0
C UMAX = 1.0
      RINFO(3) = 1.0
      RINFO(4) = 1.0
C SPCWGT = 1.0
      RINFO(5) = 1.0
      RINFO(6) = 1.0
C TIMWGT = 1.0
      RINFO(7) = 1.0
      RINFO(8) = 1.0
C
C Call main routine
      CALL  VLUGR2 (NPDE, T, TOUT, DT, XL, YL, XR, YU, DX, DY,
     +   TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK, LWK, LENLWK,
     +   MNTR)
      PRINT *, 'VLUGR2 returned with MNTR=', MNTR
C
C Save info on file
      OPEN(UNIT=LUNDMP,FILE='DUMP',FORM='UNFORMATTED')
      CALL DUMP (LUNDMP, RWK, IWK)
      CLOSE(LUNDMP)
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
C XL     : OUT. X-coordinate of lowerleft point of virtual rectangle
C YL     : OUT. Y-coordinate of lowerleft point of virtual rectangle
C XR     : OUT. X-coordinate of upperright point of virtual rectangle
C YU     : OUT. Y-coordinate of upperright point of virtual rectangle
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
C LBND   : (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C                      structure
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
C Square domain [0,1]x[0,1] with holes. Dirichlet boundaries.
C
      INTEGER NX, NY
      PARAMETER (NX = 10, NY = 10)
      INTEGER IDOM((NX+1)*(NY+1))
C
      INTEGER I, IPT, J, NROWS, NPTS, NBNDS

      NPTS = (NX+1)*(NY+1) - (NX-2) - 2*3 - 2
      IF (MAXPTS .LT. NPTS) THEN
         INIDOM = .FALSE.
         MAXPTS = NPTS
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
C
C Make grid structure
      LROW(0) = NROWS
      IPT = 1
      DO 10 I = 0, 0
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 15 J = 0, 2
            ICOL(IPT) = J
            IPT = IPT + 1
   15    CONTINUE
   10 CONTINUE
      DO 20 I = 1, 3
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 23 J = 0, 2
            ICOL(IPT) = J
            IPT = IPT + 1
   23    CONTINUE
         DO 26 J = 3, 10
            ICOL(IPT) = J
            IPT = IPT + 1
   26    CONTINUE
   20 CONTINUE
      DO 30 I = 4, 4
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 33 J = 0, 2
            ICOL(IPT) = J
            IPT = IPT + 1
   33    CONTINUE
         DO 36 J = 3, 5
            ICOL(IPT) = J
            IPT = IPT + 1
   36    CONTINUE
         DO 39 J = 8, 10
            ICOL(IPT) = J
            IPT = IPT + 1
   39    CONTINUE
   30 CONTINUE
      DO 40 I = 5, 7
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 43 J = 0, 2
            ICOL(IPT) = J
            IPT = IPT + 1
   43    CONTINUE
         DO 46 J = 3, 10
            ICOL(IPT) = J
            IPT = IPT + 1
   46    CONTINUE
   40 CONTINUE
      DO 50 I = 8, 10
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 56 J = 0, 8
            ICOL(IPT) = J
            IPT = IPT + 1
   56    CONTINUE
   50 CONTINUE
      LROW(NROWS+1) = IPT
C
C Boundaries
      NBNDS = 28
      ILBND(1) = 1
      LLBND(1) = 1
         IPT = 2
         LBND(LLBND(1)) = IPT
      ILBND(2) = 2
      LLBND(2) = LLBND(1) + 1
         IPT = 4
         LBND(LLBND(2)  ) = IPT
         IPT = 15
         LBND(LLBND(2)+1) = IPT
         IPT = 26
         LBND(LLBND(2)+2) = IPT
         IPT = 37
         LBND(LLBND(2)+3) = IPT
         IPT = 46
         LBND(LLBND(2)+4) = IPT
         IPT = 57
         LBND(LLBND(2)+5) = IPT
         IPT = 68
         LBND(LLBND(2)+6) = IPT
         IPT = 79
         LBND(LLBND(2)+7) = IPT
         IPT = 88
         LBND(LLBND(2)+8) = IPT
      ILBND(3) = 3
      LLBND(3) = LLBND(2) + 9
      DO 130 J = 0, 6 
         IPT = 98+J
         LBND(LLBND(3)+J) = IPT
  130 CONTINUE
      ILBND(4) = 4
      LLBND(4) = LLBND(3) + 7
         IPT = 96
         LBND(LLBND(4)) = IPT
      ILBND(5) = 1
      LLBND(5) = LLBND(4) + 1
      DO 150 J = 0, 4
         IPT = 86-J
         LBND(LLBND(5)+J) = IPT
  150 CONTINUE
      ILBND(6) = 4
      LLBND(6) = LLBND(5) + 5
      DO 160 J = 0, 6
         IPT = LBND(LLBND(2)+J) + 2
         LBND(LLBND(6)+J) = IPT
  160 CONTINUE
      ILBND(7) = 1
      LLBND(7) = LLBND(6) + 7
      DO 170 J = 0, 5 
         IPT = 8+J
         LBND(LLBND(7)+J) = IPT
  170 CONTINUE
      ILBND(8) = 2
      LLBND(8) = LLBND(7) + 6
      DO 180 J = 0, 4
         IPT = LBND(LLBND(6)+J+1) + 1
         LBND(LLBND(8)+J) = IPT
  180 CONTINUE
      ILBND(9) = 3
      LLBND(9) = LLBND(8) + 5
      DO 190 J = 0, 5 
         IPT = 72+J
         LBND(LLBND(9)+J) = IPT
  190 CONTINUE
      ILBND(10) = 4
      LLBND(10) = LLBND(9) + 6
         IPT = 67
         LBND(LLBND(10)  ) = IPT
         IPT = 56
         LBND(LLBND(10)+1) = IPT
         IPT = 45
         LBND(LLBND(10)+2) = IPT
         IPT = 36
         LBND(LLBND(10)+3) = IPT
         IPT = 25
         LBND(LLBND(10)+4) = IPT
      ILBND(11) = 1
      LLBND(11) = LLBND(10) + 5
         IPT = 52
         LBND(LLBND(11)  ) = IPT
         IPT = 53
         LBND(LLBND(11)+1) = IPT
      ILBND(12) = 2
      LLBND(12) = LLBND(11) + 2
         IPT = 43
         LBND(LLBND(12)  ) = IPT
      ILBND(13) = 3
      LLBND(13) = LLBND(12) + 1
         IPT = 32
         LBND(LLBND(13)  ) = IPT
         IPT = 33
         LBND(LLBND(13)+1) = IPT
      ILBND(14) = 4
      LLBND(14) = LLBND(13) + 2
         IPT = 42
         LBND(LLBND(14)  ) = IPT
C
      ILBND(15) = 12
      LLBND(15) = LLBND(14) + 1
         IPT = 1
         LBND(LLBND(15)) = IPT
      ILBND(16) = 23
      LLBND(16) = LLBND(15) + 1
         IPT = 97
         LBND(LLBND(16)) = IPT
      ILBND(17) = 34
      LLBND(17) = LLBND(16) + 1
         IPT = 105
         LBND(LLBND(17)) = IPT
      ILBND(18) = 41
      LLBND(18) = LLBND(17) + 1
         IPT = 87
         LBND(LLBND(18)) = IPT
      ILBND(19) = 14
      LLBND(19) = LLBND(18) + 1
         IPT = 81
         LBND(LLBND(19)) = IPT
      ILBND(20) = 41
      LLBND(20) = LLBND(19) + 1
         IPT = 3
         LBND(LLBND(20)) = IPT
      ILBND(21) = 12
      LLBND(21) = LLBND(20) + 1
         IPT = 7
         LBND(LLBND(21)) = IPT
      ILBND(22) = 23
      LLBND(22) = LLBND(21) + 1
         IPT = 71
         LBND(LLBND(22)) = IPT
      ILBND(23) = 34
      LLBND(23) = LLBND(22) + 1
         IPT = 78
         LBND(LLBND(23)) = IPT
      ILBND(24) = 41
      LLBND(24) = LLBND(23) + 1
         IPT = 14
         LBND(LLBND(24)) = IPT
      ILBND(25) = 14
      LLBND(25) = LLBND(24) + 1
         IPT = 51
         LBND(LLBND(25)) = IPT
      ILBND(26) = 43
      LLBND(26) = LLBND(25) + 1
         IPT = 31
         LBND(LLBND(26)) = IPT
      ILBND(27) = 32
      LLBND(27) = LLBND(26) + 1
         IPT = 34
         LBND(LLBND(27)) = IPT
      ILBND(28) = 21
      LLBND(28) = LLBND(27) + 1
         IPT = 54
         LBND(LLBND(28)) = IPT
C
      LLBND(29) = LLBND(28) + 1
      LLBND( 0) = NBNDS

C No internal boundaries
C (only necessary because we want to print the domain)
      LLBND(NBNDS+2) = LLBND(NBNDS+1)
      PRINT *, 'Input domain:'
      CALL PRDOM (LROW, IROW, ICOL, LLBND, ILBND, LBND,
     +   IDOM, NX, NY)

      RETURN
      END
      SUBROUTINE CHSPCM (T, LEVEL, NPTS, X, Y, NPDE, U, SPCMON, TOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LEVEL, NPTS, NPDE
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), U(NPTS,NPDE), SPCMON(NPTS),
     +   TOL
C
Ccc PURPOSE:
C Force grid refinement.
C If for a node IPT SPCMON(IPT) > TOL the 16 surrounding cells will be
C refined.
C
Ccc PARAMETER DESCRIPTION:
C T      : IN.  Current value of time variable
C LEVEL  : IN.  Current grid level
C NPTS   : IN.  Number of grid points at this level
C X      : IN.  Array of X-coordinates for the gridpoints
C Y      : IN.  Array of Y-coordinates for the gridpoints
C NPDE   : IN.  Number of PDE components
C U      : IN.  Array of PDE components for the gridpoints
C SPCMON : INOUT.
C          IN:  Space monitor values as determined by VLUGR2
C          OUT: Changed to a value > TOL where refinement is required
C TOL    : IN.  Tolerance with which SPCMON will be compared
C
C-----------------------------------------------------------------------
C
      INTEGER I
C
      IF (LEVEL .GE. 3) RETURN
      DO 10 I = 1, NPTS
         IF (ABS(X(I)-1.0) .LT. 0.0001 .AND. 
     +       ABS(Y(I)-0.1) .LT. 0.0001) THEN
            SPCMON(I) = 2*TOL
         ENDIF
   10 CONTINUE
C
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
