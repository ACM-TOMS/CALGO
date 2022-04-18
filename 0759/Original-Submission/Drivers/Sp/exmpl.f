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
C     PARAMETER (NRRMAX = 1, MAXLR = 3, MAXL = 15)
      PARAMETER (TOLLSC = TOLNEW/10)
      COMMON /IGCRO/ IDIAGP
      SAVE /IGCRO/
C
C end INCLUDE 'PARGCRO'
C
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK, LENLWK
      PARAMETER (MXLEV=2, NPD=3, NPTS=61000)
      PARAMETER (LENIWK=NPTS*(7*MXLEV+7),
     +           LENRWK=NPTS*NPD*(5*MXLEV+13 + (2*MAXLR+MAXL+7)),
     +           LENLWK=2*NPTS)
C
C-----------------------------------------------------------------------
C
      INTEGER LUNDMP
      PARAMETER (LUNDMP = 89)
C
      INTEGER NPDE, INFO(7), IWK(LENIWK), MNTR
      LOGICAL LWK(LENLWK)
      REAL T, TOUT, DT, XL, YF, ZD, XR, YB, ZU, DX, DY, DZ,
     +   TOLS, TOLT, RINFO(2+3*NPD), RWK(LENRWK)

C First call of VLUGR3
      MNTR = 0
      NPDE = 3
      T    = 0.0
      TOUT = 1.0
      DT   = 0.001
C Since domain is not a rectangular prism the grid parameters need not
C to be specified here (cf. INIDOM)
      TOLS = 0.1
      TOLT = 0.1
      INFO(1) = 1
C MAXLEV
      INFO(2) = 4
C Domain not a rectangular prism
      INFO(3) = 1
C Linear system solver: matrix-free GCRO + Diagonal scaling
C (no first order derivatives at the boundaries)
      INFO(4) = 13
      OPEN (UNIT=61,FILE='RunInfo')
C Write integration history to unit # 61
      INFO(5) = 61
C Write Newton info to unit # 61
      INFO(6) = 61
C Write GCRO info to unit # 61
      INFO(7) = 61
C DTMIN = 1E-7
      RINFO(1) = 1.0E-7
C DTMAX = 1.0
      RINFO(2) = 1.0
C UMAX = 1.0
      RINFO(3) = 1.0
      RINFO(4) = 1.0
      RINFO(5) = 1.0
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
      CALL  VLUGR3 (NPDE, T, TOUT, DT, XL,YF,ZD, XR,YB,ZU, DX, DY, DZ,
     +   TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK, LWK, LENLWK,
     +   MNTR)
      PRINT *, 'VLUGR3 returned with MNTR=', MNTR
C
C Save info on file
      OPEN(UNIT=LUNDMP,FILE='DUMP',FORM='UNFORMATTED')
      CALL DUMP (LUNDMP, RWK, IWK)
      CLOSE(LUNDMP)
      END
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
C Define initial domain. NB. Boundaries should consist of as many points
C as are necessary to employ second order space discretization, i.e.,
C a boundary enclosing the internal part of the domain should not
C include less than 3 grid points in any coordinate direction including
C the corners. If Neumann boundaries are used the minimum is 4 since
C otherwise the Jacobian matrix will be singular.
C
C A (virtual) box is placed upon the (irregular) domain.
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
C LBND   : OUT. (NBIPTS)
C          LBND(IBPT): pointer to boundary point in actual grid
C                      structure
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
C Domain [0,1]x[0,1]x[0,1] with hole in the middle and projection
C at [1,1.333]x[0,1]x[0.666,1].
C Virtual box: [0,1.333]x[0,1]x[0,1]
C
      INTEGER NX, NY, NZ
      PARAMETER (NX = 8, NY = 6, NZ = 6)
      INTEGER IDOM(0:(NX+1)*(NY+1)*(NZ+1))
C
      INTEGER I, IPT, IR, J, K, NPLNS, NROWS, NPTS, NBNDS,
     +   NPTSP1, NPTSP2

      NPLNS = NZ+1
      NROWS = (NY+1)*NPLNS
      NPTS  = (NX+1)*NROWS-1-2*(NY+1)*(NZ-2)
      IF (MAXPTS .LT. NPTS) THEN
         INIDOM = .FALSE.
         MAXPTS = NPTS
         RETURN
      ELSE
         INIDOM = .TRUE.
      ENDIF

      XL = 0.0
      YF = 0.0
      ZD = 0.0
      XR = 4.0/3.0
      YB = 1.0
      ZU = 1.0
      DX = (XR-XL)/NX
      DY = (YB-YF)/NY
      DZ = (ZU-ZD)/NZ
C
C Make grid structure
      LPLN(0) = NPLNS
      IPT = 1
      IR = 1
      DO 10 K = 0, 2
         LPLN(K+1) = IR
         IPLN(K+1) = K
         DO 20 I = 0, NY
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR+1
            DO 30 J = 0, NX-2
               ICOL(IPT) = J
               IPT = IPT + 1
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      K = 3
         LPLN(K+1) = IR
         IPLN(K+1) = K
         DO 40 I = 0, 2
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR+1
            DO 50 J = 0, NX-2
               ICOL(IPT) = J
               IPT = IPT + 1
   50       CONTINUE
   40    CONTINUE
         I = 3
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR+1
            DO 60 J = 0, 2
               ICOL(IPT) = J
               IPT = IPT + 1
   60       CONTINUE
            DO 70 J = 4, NX-2
               ICOL(IPT) = J
               IPT = IPT + 1
   70       CONTINUE
         DO 80 I = 4, NY
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR+1
            DO 90 J = 0, NX-2
               ICOL(IPT) = J
               IPT = IPT + 1
   90       CONTINUE
   80    CONTINUE
      DO 100 K = 4, NZ
         LPLN(K+1) = IR
         IPLN(K+1) = K
         DO 110 I = 0, NY
            LROW(IR) = IPT
            IROW(IR) = I
            IR = IR+1
            DO 120 J = 0, NX
               ICOL(IPT) = J
               IPT = IPT + 1
  120       CONTINUE
  110    CONTINUE
  100    CONTINUE
      LROW(NROWS+1) = NPTS+1
      LPLN(NPLNS+1) = NROWS+1
C
Ccc Boundaries
      NBNDS = 14
      ILBND( 1) = 1
      ILBND( 2) = 2
      ILBND( 3) = 3
      ILBND( 4) = 2
      ILBND( 5) = 3
      ILBND( 6) = 4
      ILBND( 7) = 5
      ILBND( 8) = 6
      ILBND( 9) = 1
      ILBND(10) = 2
      ILBND(11) = 3
      ILBND(12) = 4
      ILBND(13) = 5
      ILBND(14) = 6
      LLBND( 0) = NBNDS
      LLBND( 1) = 1
      LLBND( 2) = LLBND( 1) + (NY+1)*(NZ+1)
      LLBND( 3) = LLBND( 2) + (NX-1)*(NY+1)
      LLBND( 4) = LLBND( 3) + (NY+1)*(NZ-1)
      LLBND( 5) = LLBND( 4) +     3 *(NY+1)
      LLBND( 6) = LLBND( 5) + (NY+1)* 3
      LLBND( 7) = LLBND( 6) + (NX+1)*(NY+1)
      LLBND( 8) = LLBND( 7) + (NX-1)*(NZ-2)+(NX+1)*3
      LLBND( 9) = LLBND( 8) + (NX-1)*(NZ-2)+(NX+1)*3
      LLBND(10) = LLBND( 9) + 9
      LLBND(11) = LLBND(10) + 9
      LLBND(12) = LLBND(11) + 9
      LLBND(13) = LLBND(12) + 9
      LLBND(14) = LLBND(13) + 9
      LLBND(15) = LLBND(14) + 9
C
Ccc Outer planes
C Left boundary plane pointers
      NPTSP1 = (NX-1)*(NY+1)
      NPTSP2 = (NX+1)*(NY+1)
      DO 200 K = 0, 3
         DO 201 I = 0, NY
            IPT = K*NPTSP1 + I*(NX-1) + 1
            IF (K .EQ. 3 .AND. I .GT. 3) IPT = IPT-1
            LBND(LLBND(1)+K*(NY+1)+I) = IPT
  201    CONTINUE
  200 CONTINUE
      DO 202 K = NZ-2, NZ
         DO 203 I = 0, NY
            IPT = (NZ-2)*NPTSP1+(K-NZ+2)*NPTSP2 + I*(NX+1)
            LBND(LLBND(1)+K*(NY+1)+I) = IPT
  203    CONTINUE
  202 CONTINUE
C Right boundary plane pointers
      DO 210 K = 0, 3
         DO 211 I = 0, NY
            IPT = (K+1)*NPTSP1 - I*(NX-1)
            IF (K .EQ. 3 .AND. I .LE. 3) IPT = IPT-1
            LBND(LLBND(3)+K*(NY+1)+I) = IPT
  211    CONTINUE
  210 CONTINUE
      K = NZ-2
         DO 209 I = 0, NY
            IPT = 4*NPTSP1 + NPTSP2-3 - I*(NX+1)
            LBND(LLBND(3)+K*(NY+1)+I) = IPT
  209    CONTINUE
      DO 212 I = 0, NY
         DO 213 J = NX-2, NX
            IPT = NPTSP1*(NZ-2)+I*(NX+1) + J
            LBND(LLBND(4)+I*3+J-NX+2) = IPT
  213    CONTINUE
  212 CONTINUE
      DO 214 K = NZ-2, NZ
         DO 215 I = 0, NY
            IPT = NPTSP1*(NZ-2)+(K-NZ+3)*NPTSP2 - I*(NX+1) - 1
            LBND(LLBND(5)+(K-NZ+2)*(NY+1)+I) = IPT
  215    CONTINUE
  214 CONTINUE
C Down and up boundary plane pointers
      DO 220 I = 0, NY
         DO 221 J = 0, NX-2
            IPT = I*(NX-1) + J + 1
            LBND(LLBND(2)+I*(NX-1)+J) = IPT
  221    CONTINUE
  220 CONTINUE
      DO 230 I = 0, NY
         DO 231 J = 0, NX
            IPT = (NPTS - (I*(NX+1)+J))
            LBND(LLBND(6)+I*(NX+1)+J) = IPT
  231    CONTINUE
  230 CONTINUE
C Front and back boundary plane pointers
      DO 240 K = 0, 3
         DO 241 J = 0, NX-2
            IPT = K*NPTSP1 + J + 1
            LBND(LLBND(7)+K*(NX-1)+J) = IPT
            IPT = (K+1)*NPTSP1 - J
            IF (K .EQ. 3) IPT = IPT-1
            LBND(LLBND(8)+K*(NX-1)+J) = IPT
  241    CONTINUE
  240 CONTINUE
      DO 242 K = NZ-2, NZ
         DO 243 J = 0, NX
            IPT = (NZ-2)*NPTSP1+(K-NZ+2)*NPTSP2 + J
            LBND(LLBND(7)+4*(NX-1)+(K-NZ+2)*(NX+1)+J) = IPT
            IPT = 4*NPTSP1+(K-NZ+3)*NPTSP2 - J - 1
            LBND(LLBND(8)+4*(NX-1)+(K-NZ+2)*(NX+1)+J) = IPT
  243    CONTINUE
  242 CONTINUE
C
Ccc Inner planes
C Left and right boundary plane pointers
      LBND(LLBND( 9)  ) = 117
      LBND(LLBND( 9)+1) = 124
      LBND(LLBND( 9)+2) = 131
      LBND(LLBND( 9)+3) = 166
      LBND(LLBND( 9)+4) = 172
      LBND(LLBND( 9)+5) = 179
      LBND(LLBND( 9)+6) = 218
      LBND(LLBND( 9)+7) = 227
      LBND(LLBND( 9)+8) = 236
      LBND(LLBND(11)  ) = 115
      LBND(LLBND(11)+1) = 122
      LBND(LLBND(11)+2) = 129
      LBND(LLBND(11)+3) = 164
      LBND(LLBND(11)+4) = 171
      LBND(LLBND(11)+5) = 177
      LBND(LLBND(11)+6) = 216
      LBND(LLBND(11)+7) = 225
      LBND(LLBND(11)+8) = 234
C Down and up boundary plane pointers
      DO 260 I = 0, 2
         DO 270 J = 0, 2
            LBND(LLBND(10)+I*3+J) = 236 - (I*(NX+1)+J)
            LBND(LLBND(12)+I*3+J) = 115 + I*(NX-1) + J
  270    CONTINUE
  260 CONTINUE
C Front and back boundary plane pointers
      LBND(LLBND(13)  ) = 129
      LBND(LLBND(13)+1) = 130
      LBND(LLBND(13)+2) = 131
      LBND(LLBND(13)+3) = 177
      LBND(LLBND(13)+4) = 178
      LBND(LLBND(13)+5) = 179
      LBND(LLBND(13)+6) = 234
      LBND(LLBND(13)+7) = 235
      LBND(LLBND(13)+8) = 236
      LBND(LLBND(14)  ) = 115
      LBND(LLBND(14)+1) = 116
      LBND(LLBND(14)+2) = 117
      LBND(LLBND(14)+3) = 164
      LBND(LLBND(14)+4) = 165
      LBND(LLBND(14)+5) = 166
      LBND(LLBND(14)+6) = 216
      LBND(LLBND(14)+7) = 217
      LBND(LLBND(14)+8) = 218
C
      LLBND(NBNDS+2) = LLBND(NBNDS+1)
      PRINT *, 'Input domain:'
      CALL PRDOM (LPLN, IPLN, LROW, IROW, ICOL, LLBND, ILBND, LBND,
     +   IDOM, NX, NY, NZ)

      RETURN
      END
      SUBROUTINE CHSPCM (T, LEVEL, NPTS, X, Y, Z, NPDE, U, SPCMON, TOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LEVEL, NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE), SPCMON(NPTS), TOL
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
      INTEGER I
C
      IF (LEVEL .GE. 3) RETURN
      DO 10 I = 1, NPTS
         IF (ABS(X(I)-1.0) .LT. 0.0001 .AND.
     +       ABS(Y(I)-0.5) .LT. 0.0001 .AND.
     +       ABS(Z(I)) .LT. 0.0001) THEN
            SPCMON(I) = 2*TOL
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE MONITR (T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB,
     +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
      REAL T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB, SOL(*)
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
C             LPLN(0) = NPLNS: Actual # planes in LROW
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
C
C Local arrays:
      INTEGER MAXPTS, NPDE
      PARAMETER (MAXPTS=250000, NPDE=3)
      REAL X(MAXPTS), Y(MAXPTS), Z(MAXPTS), UEX(MAXPTS*NPDE)
C
C-----------------------------------------------------------------------
C
      INTEGER MAXLEV, LEVEL, LLPLN, LIPLN, LLROW, LIROW, LICOL,
     +   NPLNS, NROWS, NPTS
      REAL DX, DY, DZ
C
C Loop over the grid levels from coarse to fine.
C Get physical coordinates of grid points
C Compute ||err||_max
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
         IF (NPTS .GT. MAXPTS) stop 'MONITR'
         LIROW  = LLROW+NROWS+1
         LICOL  = LIROW+NROWS
         CALL SETXYZ (XL, YF, ZD, DX, DY, DZ,
     +      ISTRUC(LLPLN), ISTRUC(LIPLN),
     +      ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL), X, Y, Z)
         DX = DX/2
         DY = DY/2
         DZ = DZ/2
         CALL PRERR (LEVEL, T, NPTS, NPDE, X, Y, Z, SOL(LSOL(LEVEL)+1),
     +      UEX)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PRERR (LEVEL, T, NPTS, NPDE, X, Y, Z, U, UEX)
      INTEGER LEVEL, NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE), UEX(NPTS,NPDE)

      INTEGER I,J
      REAL RMAX(3)

      CALL PDEIV (T, X, Y, Z, UEX, NPTS, NPDE)
      DO 1 J = 1,NPDE
      RMAX(J) = 0.0
      DO 10 I = 1, NPTS
         RMAX(J) = MAX(RMAX(J),ABS(UEX(I,J)-U(I,J)))
   10 CONTINUE
    1 CONTINUE
      WRITE(28,'(''Error at T='',E9.3,'', level='',I1,'' :'',
     +           I10,3E12.3)')
     +   T, LEVEL, NPTS, (RMAX(J), J=1, NPDE)

      RETURN
      END
      SUBROUTINE PDEIV (T, X, Y, Z, U, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE)
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
C
      INTEGER I
      REAL EPS
      PARAMETER (EPS = 5.0E-3)

      DO 10 I = 1, NPTS
         U(I,1) = 1-0.5/(1+EXP((-X(I)+Y(I)+Z(I)-0.75*T)/(4*EPS)))
         U(I,2) = 1.5-U(I,1)
         U(I,3) = 1.5-U(I,1)
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PDEF (T, X, Y, Z, U,
     +   UT, UX, UY, UZ, UXX, UYY, UZZ, UXY, UXZ, UYZ, RES, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
     +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE), UZ(NPTS,NPDE),
     +     UXX(NPTS,NPDE), UYY(NPTS,NPDE), UZZ(NPTS,NPDE),
     +     UXY(NPTS,NPDE), UXZ(NPTS,NPDE), UYZ(NPTS,NPDE),
     +     RES(NPTS,NPDE)
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
      INTEGER I
      REAL EPS
      PARAMETER (EPS = 5.0E-3)

      DO 10 I = 1, NPTS
         RES(I,1) = UT(I,1) + U(I,1)*UX(I,1) +
     +              U(I,2)*UY(I,1) + U(I,3)*UZ(I,1) -
     +              EPS*(UXX(I,1)+UYY(I,1)+UZZ(I,1))
         RES(I,2) = UT(I,2) + U(I,1)*UX(I,2) +
     +              U(I,2)*UY(I,2) + U(I,3)*UZ(I,2) -
     +              EPS*(UXX(I,2)+UYY(I,2)+UZZ(I,2))
         RES(I,3) = UT(I,3) + U(I,1)*UX(I,3) +
     +              U(I,2)*UY(I,3) + U(I,3)*UZ(I,3) -
     +              EPS*(UXX(I,3)+UYY(I,3)+UZZ(I,3))
   10 CONTINUE

      RETURN
      END
      SUBROUTINE PDEBC (T, X, Y, Z, U, UT, UX, UY, UZ, RES, NPTS, NPDE,
     +   LLBND, ILBND, LBND)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE, LLBND(0:*), ILBND(*), LBND(*)
      REAL T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
     +     UT(NPTS,NPDE), UX(NPTS,NPDE), UY(NPTS,NPDE), UZ(NPTS,NPDE),
     +     RES(NPTS,NPDE)
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
C-----------------------------------------------------------------------
C
      INTEGER I, J, K, NBNDS
      REAL EPS, UI
      PARAMETER (EPS = 5.0E-3)

      NBNDS = LLBND(0)
      DO 10 J = 1, NBNDS
         DO 20 K = LLBND(J), LLBND(J+1)-1
            I = LBND(K)
            UI = 1-0.5/(1+EXP((-X(I)+Y(I)+Z(I)-0.75*T)/(4*EPS)))
            RES(I,1) = U(I,1) - UI
            RES(I,2) = U(I,2) - (1.5-UI)
            RES(I,3) = U(I,3) - (1.5-UI)
   20    CONTINUE
   10 CONTINUE

      RETURN
      END
