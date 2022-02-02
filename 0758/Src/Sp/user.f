      LOGICAL FUNCTION INIDOM (MAXPTS, XL, YL, XR, YU, DX, DY,
     +   LROW, IROW, ICOL, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER MAXPTS, LROW(0:*), IROW(*), ICOL(*),
     +   LLBND(0:*), ILBND(*), LBND(*)
      REAL XL, YL, XR, YU, DX, DY
C
Ccc PURPOSE:
C Define grid for initial rectangular domain ((XL,YL),(XR,YU)) in
C in physical coordinates and ((0,0),(Nx,Ny) in column, resp. row
C coordinates, where Nx = (XR-XL)/DX and Ny = (YU-YL)/DY.
C Only real grid points are stored.
C The coordinate values of the initial grid should be stored rowwise,
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
C XL     : IN.  X-coordinate of lowerleft point of rectangle
C YL     : IN.  Y-coordinate of lowerleft point of rectangle
C XR     : IN.  X-coordinate of upperright point of rectangle
C YU     : IN.  Y-coordinate of upperright point of rectangle
C DX     : IN.  Grid width in X-direction
C DY     : IN.  Grid width in Y-direction
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
      INTEGER NX, NY, I, IPT, J, NROWS, NPTS, NBNDS

      NX = NINT((XR-XL)/DX)
      NY = NINT((YU-YL)/DY)
C
Ccc Make initial grid
      NPTS = (NX+1)*(NY+1)
      IF (MAXPTS .LT. NPTS) THEN
         INIDOM = .FALSE.
         MAXPTS = NPTS
         RETURN
      ELSE
         INIDOM = .TRUE.
      ENDIF
      NROWS = NY+1
C
C Make grid structure
      LROW(0) = NROWS
      IPT = 1
      DO 10 I = 0, NY
         LROW(I+1) = IPT
         IROW(I+1) = I
         DO 20 J = 0, NX
            ICOL(IPT) = J
            IPT = IPT + 1
   20    CONTINUE
   10 CONTINUE
      LROW(NROWS+1) = NPTS+1
C
C Boundaries
      NBNDS = 8
      ILBND(1) = 1
      ILBND(2) = 2
      ILBND(3) = 3
      ILBND(4) = 4
      ILBND(5) = 12
      ILBND(6) = 23
      ILBND(7) = 34
      ILBND(8) = 41
      LLBND(0) = NBNDS
      LLBND(1) = 1
      LLBND(2) = LLBND(1) + (NX-1)
      LLBND(3) = LLBND(2) + (NY-1)
      LLBND(4) = LLBND(3) + (NX-1)
      LLBND(5) = LLBND(4) + (NY-1)
      LLBND(6) = LLBND(5) + 1
      LLBND(7) = LLBND(6) + 1
      LLBND(8) = LLBND(7) + 1
      LLBND(9) = LLBND(8) + 1
C Lower and upper boundary pointers
      DO 50 J = 1, NX-1
         LBND(LLBND(1)+J-1) = J + 1
         LBND(LLBND(3)+J-1) = NPTS - J
   50 CONTINUE
C Left and right boundary pointers
      DO 60 I = 1, NY-1
         LBND(LLBND(2)+I-1) = I*(NX+1) + 1
         LBND(LLBND(4)+I-1) = (I+1)*(NX+1)
   60 CONTINUE
C Corners
      LBND(LLBND(5)) = 1
      LBND(LLBND(6)) = NPTS - NX
      LBND(LLBND(7)) = NPTS
      LBND(LLBND(8)) = NX+1

      RETURN
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
      REAL F(NPTS*NPDE), T, X(*), Y(*), U(*), A0, DT, DX, DY, UIB(*),
     +   UT(*), UX(*), UY(*), UXX(*), UXY(*), UYY(*),
     +   ABSTOL(*), DEL(NPTS), WORK(2*NPTS*NPDE),
     +   FU(NPTS*NPDE,NPDE), FUX(NPTS*NPDE,NPDE), FUY(NPTS*NPDE,NPDE),
     +   FUXX(NPTS*NPDE,NPDE),FUXY(NPTS*NPDE,NPDE),FUYY(NPTS*NPDE,NPDE)
C
Ccc PURPOSE:
C Compute derivatives of residual wrt (derivatives of) U by numerical
C differencing
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
Ccc EXTERNALS USED:
      EXTERNAL PERTRB, PRTRBU, RES
C
C-----------------------------------------------------------------------
C
      INTEGER I, IC, ICPTB, IPT, LUTBAR
      REAL FACX, FACY, FACXX, FACXY, FACYY, TOL

      LUTBAR = 1 + NPTS*NPDE
C
Ccc How to decide if derivatives are `zero'?
C Take `zero'-value of U divided by the grid width
      FACX  = 1/(2*DX)
      FACY  = 1/(2*DY)
      FACXX = 1/DX**2
      FACXY = 1/(2*DX*2*DY)
      FACYY = 1/DY**2
C
Ccc Loop over the components of the (derivatives of) U
      DO 10 ICPTB = 1, NPDE
C
C dF(U,Ut)/dU
         TOL = ABSTOL(ICPTB)
         CALL PRTRBU (ICPTB, NPTS, NPDE, U, A0, DT, UT, TOL, DEL,
     +      WORK, WORK(LUTBAR))
         CALL RES (T, X, Y, NPTS, NPDE, WORK, LLBND, ILBND, LBND, UIB,
     +      WORK(LUTBAR), UX, UY, UXX, UXY, UYY, FU(1,ICPTB))
         DO 20 IC = 1, NPDE
         DO 20 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FU(I,ICPTB) = (FU(I,ICPTB) - F(I)) / DEL(IPT)
   20    CONTINUE
C
C dF(Ux)/dUx
         TOL = ABSTOL(ICPTB)*FACX
         CALL PERTRB (ICPTB, NPTS, NPDE, UX, TOL, DEL, WORK)
         CALL RES (T, X, Y, NPTS, NPDE, U, LLBND, ILBND, LBND, UIB,
     +      UT, WORK, UY, UXX, UXY, UYY, FUX(1,ICPTB))
         DO 40 IC = 1, NPDE
         DO 40 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUX(I,ICPTB) = (FUX(I,ICPTB) - F(I)) / DEL(IPT)
   40    CONTINUE
C
C dF(Uy)/dUy
         TOL = ABSTOL(ICPTB)*FACY
         CALL PERTRB (ICPTB, NPTS, NPDE, UY, TOL, DEL, WORK)
         CALL RES (T, X, Y, NPTS, NPDE, U, LLBND, ILBND, LBND, UIB,
     +      UT, UX, WORK, UXX, UXY, UYY, FUY(1,ICPTB))
         DO 50 IC = 1, NPDE
         DO 50 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUY(I,ICPTB) = (FUY(I,ICPTB) - F(I)) / DEL(IPT)
   50    CONTINUE
C
C dF(Uxx)/dUxx
         TOL = ABSTOL(ICPTB)*FACXX
         CALL PERTRB (ICPTB, NPTS, NPDE, UXX, TOL, DEL, WORK)
         CALL RES (T, X, Y, NPTS, NPDE, U, LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, WORK, UXY, UYY, FUXX(1,ICPTB))
         DO 60 IC = 1, NPDE
         DO 60 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUXX(I,ICPTB) = (FUXX(I,ICPTB) - F(I)) / DEL(IPT)
   60    CONTINUE
C
C dF(Uxy)/dUxy
         TOL = ABSTOL(ICPTB)*FACXY
         CALL PERTRB (ICPTB, NPTS, NPDE, UXY, TOL, DEL, WORK)
         CALL RES (T, X, Y, NPTS, NPDE, U, LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UXX, WORK, UYY, FUXY(1,ICPTB))
         DO 70 IC = 1, NPDE
         DO 70 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUXY(I,ICPTB) = (FUXY(I,ICPTB) - F(I)) / DEL(IPT)
   70    CONTINUE
C
C dF(Uyy)/dUyy
         TOL = ABSTOL(ICPTB)*FACYY
         CALL PERTRB (ICPTB, NPTS, NPDE, UYY, TOL, DEL, WORK)
         CALL RES (T, X, Y, NPTS, NPDE, U, LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UXX, UXY, WORK, FUYY(1,ICPTB))
         DO 80 IC = 1, NPDE
         DO 80 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUYY(I,ICPTB) = (FUYY(I,ICPTB) - F(I)) / DEL(IPT)
   80    CONTINUE
   10 CONTINUE
      
      RETURN
      END
      SUBROUTINE MONITR (T, DT, DTNEW, XL, YL, DXB, DYB,
     +   LGRID, ISTRUC, LSOL, SOL)
      INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
      REAL T, DT, DTNEW, XL, YL, DXB, DYB, SOL(*)
      RETURN
      END
      SUBROUTINE CHSPCM (T, LEVEL, NPTS, X, Y, NPDE, U, SPCMON, TOL)
      INTEGER LEVEL, NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), U(NPTS,NPDE), SPCMON(NPTS), TOL
      RETURN
      END
