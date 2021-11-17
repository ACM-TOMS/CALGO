      LOGICAL FUNCTION INIDOM (MAXPTS,
     +   XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   LPLN, IPLN, LROW, IROW, ICOL, LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER MAXPTS, LPLN(0:*), IPLN(*), LROW(*), IROW(*), ICOL(*),
     +   LLBND(0:*), ILBND(*), LBND(*)
      REAL XL, YF, ZD, XR, YB, ZU, DX, DY, DZ
C
Ccc PURPOSE:
C Define grid for initial rectangular-prism domain
C ((XL,YF,ZD),(XR,YB,ZU)) in physical coordinates and
C (( 0, 0, 0),(Nx,Ny,Nz)) in computational grid coordinates,
C where Nx = (XR-XL)/DX, Ny = (YB-YF)/DY, and Nz = (ZU-ZD)/DZ.
C Only real grid points are stored.
C The coordinate values of the initial grid should be stored rowwise,
C in LPLN, IPLN, LROW, IROW, ICOL.
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
C DX     : IN. Grid width in X-direction
C DY     : IN. Grid width in Y-direction
C DZ     : IN. Grid width in Z-direction
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
C                      structure
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER NX, NY, NZ, I, IPT, IR, J, K, NPLNS, NROWS, NPTS, NBNDS,
     +   NPTSPL

      NX = NINT((XR-XL)/DX)
      NY = NINT((YB-YF)/DY)
      NZ = NINT((ZU-ZD)/DZ)
C
Ccc Make initial grid
      NPLNS = NZ+1
      NROWS = (NY+1)*NPLNS
      NPTS  = (NX+1)*NROWS
      IF (MAXPTS .LT. NPTS) THEN
         INIDOM = .FALSE.
         MAXPTS = NPTS
         RETURN
      ELSE
         INIDOM = .TRUE.
      ENDIF
C
C Make grid structure
      LPLN(0) = NPLNS
      IPT = 1
      IR = 1
      DO 10 K = 0, NZ
         LPLN(K+1) = IR
         IPLN(K+1) = K
         DO 20 I = 0, NY
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR + 1
            DO 30 J = 0, NX
               ICOL(IPT) = J
               IPT = IPT + 1
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      LROW(NROWS+1) = NPTS+1
      LPLN(NPLNS+1) = NROWS+1
C
C Boundaries
      NPTSPL = (NX+1)*(NY+1)
      NBNDS = 6
      ILBND(1) = 1
      ILBND(2) = 2
      ILBND(3) = 3
      ILBND(4) = 4
      ILBND(5) = 5
      ILBND(6) = 6
      LLBND(0) = NBNDS
      LLBND(1) = 1
      LLBND(2) = LLBND(1) + (NY+1)*(NZ+1)
      LLBND(3) = LLBND(2) + (NX+1)*(NY+1)
      LLBND(4) = LLBND(3) + (NY+1)*(NZ+1)
      LLBND(5) = LLBND(4) + (NX+1)*(NY+1)
      LLBND(6) = LLBND(5) + (NX+1)*(NZ+1)
      LLBND(7) = LLBND(6) + (NX+1)*(NZ+1)
C Left and right boundary plane pointers
      DO 100 K = 0, NZ
         DO 110 I = 0, NY
            LBND(LLBND(1)+K*(NY+1)+I) = K*NPTSPL + I*(NX+1) + 1
            LBND(LLBND(3)+K*(NY+1)+I) = (K+1)*NPTSPL - I*(NX+1)
  110    CONTINUE
  100 CONTINUE
C Down and up boundary plane pointers
      DO 120 I = 0, NY
         DO 130 J = 0, NX
            LBND(LLBND(2)+I*(NX+1)+J) = I*(NX+1) + J + 1
            LBND(LLBND(4)+I*(NX+1)+J) = NPTS - (I*(NX+1)+J)
  130    CONTINUE
  120 CONTINUE
C Front and back boundary plane pointers
      DO 140 K = 0, NZ
         DO 150 J = 0, NX
            LBND(LLBND(5)+K*(NX+1)+J) = K*NPTSPL + J + 1
            LBND(LLBND(6)+K*(NX+1)+J) = NPTS - (K*NPTSPL+J)
  150    CONTINUE
  140 CONTINUE
C
      RETURN
      END
      SUBROUTINE DERIVF (F, T, X, Y, Z, NPTS, NPDE, U,
     +   A0, DT, DX, DY, DZ,
     +   LLBND, ILBND, LBND, UIB, UT, UX, UY, UZ, UXX, UYY, UZZ,
     +   UXY, UXZ, UYZ, ABSTOL, DEL, WORK,
     +   FU, FUX, FUY, FUZ, FUXX, FUYY, FUZZ, FUXY, FUXZ, FUYZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      REAL F(NPTS*NPDE), T, X(*), Y(*), Z(*), U(*), A0, DT, DX, DY, DZ,
     +   UIB(*), UT(*), UX(*), UY(*), UZ(*), UXX(*), UYY(*), UZZ(*),
     +   UXY(*), UXZ(*), UYZ(*),
     +   ABSTOL(*), DEL(NPTS), WORK(2*NPTS*NPDE),
     +   FU(NPTS*NPDE,NPDE),
     +   FUX(NPTS*NPDE,NPDE), FUY(NPTS*NPDE,NPDE), FUZ(NPTS*NPDE,NPDE),
     +   FUXX(NPTS*NPDE,NPDE),FUYY(NPTS*NPDE,NPDE),FUZZ(NPTS*NPDE,NPDE),
     +   FUXY(NPTS*NPDE,NPDE),FUXZ(NPTS*NPDE,NPDE),FUYZ(NPTS*NPDE,NPDE)
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
Ccc EXTERNALS USED:
      EXTERNAL PERTRB, PRTRBU, RES
C
C-----------------------------------------------------------------------
C
      INTEGER I, IC, ICPTB, IPT, LUTBAR
      REAL FACX, FACY, FACZ, FACXX, FACYY, FACZZ, FACXY, FACXZ, FACYZ,
     +   TOL

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
      FACXY = 1/(2*DX*2*DY)
      FACXZ = 1/(2*DX*2*DZ)
      FACYZ = 1/(2*DY*2*DZ)
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
C
C dF(Ux)/dUx
         TOL = ABSTOL(ICPTB)*FACX
         CALL PERTRB (ICPTB, NPTS, NPDE, UX, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, WORK, UY, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUX(1,ICPTB))
         DO 40 IC = 1, NPDE
         DO 40 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUX(I,ICPTB) = (FUX(I,ICPTB) - F(I)) / DEL(IPT)
   40    CONTINUE
C
C dF(Uy)/dUy
         TOL = ABSTOL(ICPTB)*FACY
         CALL PERTRB (ICPTB, NPTS, NPDE, UY, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, WORK, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUY(1,ICPTB))
         DO 50 IC = 1, NPDE
         DO 50 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUY(I,ICPTB) = (FUY(I,ICPTB) - F(I)) / DEL(IPT)
   50    CONTINUE
C
C dF(Uz)/dUz
         TOL = ABSTOL(ICPTB)*FACZ
         CALL PERTRB (ICPTB, NPTS, NPDE, UZ, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, WORK,
     +      UXX, UYY, UZZ, UXY, UXZ, UYZ, FUZ(1,ICPTB))
         DO 60 IC = 1, NPDE
         DO 60 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUZ(I,ICPTB) = (FUZ(I,ICPTB) - F(I)) / DEL(IPT)
   60    CONTINUE
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
C
C dF(Uxy)/dUxy
         TOL = ABSTOL(ICPTB)*FACXY
         CALL PERTRB (ICPTB, NPTS, NPDE, UXY, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      UXX, UYY, UZZ, WORK, UXZ, UYZ, FUXY(1,ICPTB))
         DO 100 IC = 1, NPDE
         DO 100 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUXY(I,ICPTB) = (FUXY(I,ICPTB) - F(I)) / DEL(IPT)
  100    CONTINUE
C
C dF(Uxz)/dUxz
         TOL = ABSTOL(ICPTB)*FACXZ
         CALL PERTRB (ICPTB, NPTS, NPDE, UXZ, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      UXX, UYY, UZZ, UXY, WORK, UYZ, FUXZ(1,ICPTB))
         DO 110 IC = 1, NPDE
         DO 110 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUXZ(I,ICPTB) = (FUXZ(I,ICPTB) - F(I)) / DEL(IPT)
  110    CONTINUE
C
C dF(Uyz)/dUyz
         TOL = ABSTOL(ICPTB)*FACYZ
         CALL PERTRB (ICPTB, NPTS, NPDE, UYZ, TOL, DEL, WORK)
         CALL RES (T, X, Y, Z, NPTS, NPDE, U,
     +      LLBND, ILBND, LBND, UIB,
     +      UT, UX, UY, UZ,
     +      UXX, UYY, UZZ, UXY, UXZ, WORK, FUYZ(1,ICPTB))
         DO 120 IC = 1, NPDE
         DO 120 IPT = 1, NPTS
            I = IPT + (IC-1)*NPTS
            FUYZ(I,ICPTB) = (FUYZ(I,ICPTB) - F(I)) / DEL(IPT)
  120    CONTINUE
   10 CONTINUE
      
      RETURN
      END
      SUBROUTINE MONITR (T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB,
     +   LGRID, ISTRUC, LSOL, SOL)
      INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
      REAL T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB, SOL(*)
      RETURN
      END
      SUBROUTINE CHSPCM (T, LEVEL, NPTS, X, Y, Z, NPDE, U, SPCMON, TOL)
      INTEGER LEVEL, NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE), SPCMON(NPTS), TOL
      RETURN
      END
