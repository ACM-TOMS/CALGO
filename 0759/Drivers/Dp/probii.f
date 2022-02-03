      PROGRAM EXII
C
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK, LENLWK
      PARAMETER (MXLEV=2, NPD=2, NPTS=70000)
C     PARAMETER (LENIWK=NPTS*(7*MXLEV+43),
C    +           LENRWK=NPTS*NPD*(5*MXLEV+38*NPD+13),
      PARAMETER (LENIWK=NPTS*(7*MXLEV+24),
     +           LENRWK=NPTS*NPD*(5*MXLEV+37+21),
     +           LENLWK=2*NPTS)
C
C-----------------------------------------------------------------------
C
      INTEGER LUNDMP
      PARAMETER (LUNDMP = 89)
C
      CHARACTER*2 NR
      INTEGER I, NPDE, INFO(7), IWK(LENIWK), MNTR
      LOGICAL LWK(LENLWK)
      DOUBLE PRECISION T, TOUT, TE, DT, XL, YF, ZD, XR, YB, ZU, DX, DY,
     +   DZ,
     +   TOLS, TOLT, RINFO(2+3*NPD), RWK(LENRWK)

C First call of VLUGR3
      MNTR = 0
      NPDE = 2
      T    = 0.0
      TOUT = 1.0
      DT   = 1D-5
      XL   = 0.0
      YF   = 0.0
      ZD   = 0.0
      XR   = 1.0
      YB   = 1.0
      ZU   = 1.0
      DX   = 0.1
      DY   = 0.1
      DZ   = 0.1
      TOLS = 0.1
      TOLT = 0.1
      INFO(1) = 1
C MAXLEV
      INFO(2) = 4
C Domain a rectangular prism
      INFO(3) = 0
C Linear system solver
      print *, 'Lin.sys.solver BiCGStab or GCRO (0 / 10,11,12,13) ?'
      read *, INFO(4)
      OPEN (UNIT=61,FILE='RunInfo')
C Write integration history to unit # 61
      INFO(5) = 61
C Write Newton info to unit # 61
      INFO(6) = 61
C Write linear solver info to unit # 61
      INFO(7) = 61
C DTMIN = 1D-9
      RINFO(1) = 1.0D-9
C DTMAX = TE - TS
      RINFO(2) = 0.0
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
      DO 10 I = 1, 10
         TE = T + TOUT/10
      CALL  VLUGR3 (NPDE, T, TE, DT, XL,YF,ZD, XR,YB,ZU, DX, DY, DZ,
     +   TOLS, TOLT, INFO, RINFO, RWK, LENRWK, IWK, LENIWK, LWK, LENLWK,
     +   MNTR)
      PRINT *, 'VLUGR3 returned with MNTR=', MNTR
C
C Save info on file
      WRITE(NR,'(I2.2)') I
      OPEN(UNIT=LUNDMP,FILE='DUMP'//NR,FORM='UNFORMATTED')
      CALL DUMP (LUNDMP, RWK, IWK)
      CLOSE(LUNDMP)
   10 CONTINUE
      END
      SUBROUTINE MONITR (T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB,
     +   LGRID, ISTRUC, LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LGRID(0:*), ISTRUC(*), LSOL(*)
      DOUBLE PRECISION T, DT, DTNEW, XL, YF, ZD, DXB, DYB, DZB, SOL(*)
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
      PARAMETER (MAXPTS=50000, NPDE=2)
      DOUBLE PRECISION X(MAXPTS), Y(MAXPTS), Z(MAXPTS), UEX(MAXPTS*NPDE)
C
C-----------------------------------------------------------------------
C
      INTEGER MAXLEV, LEVEL, LLPLN, LIPLN, LLROW, LIROW, LICOL,
     +   NPLNS, NROWS, NPTS
      DOUBLE PRECISION DX, DY, DZ
C
C Loop over the grid levels from coarse to fine.
C Get physical coordinates of grid points
C Compute ||err||_max and clip negative solution values
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
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
     +   UEX(NPTS,NPDE)

      INTEGER I
      DOUBLE PRECISION RMAX1, RMAX2

      CALL PDEIV (T, X, Y, Z, UEX, NPTS, NPDE)
      RMAX1 = 0.0
      RMAX2 = 0.0
      DO 10 I = 1, NPTS
         RMAX1 = MAX(RMAX1,ABS(UEX(I,1)-U(I,1)))
         RMAX2 = MAX(RMAX2,ABS(UEX(I,2)-U(I,2)))
         IF (U(I,1) .LT. 0) U(I,1) = 0.0
         IF (U(I,2) .LT. 0) U(I,2) = 0.0
   10 CONTINUE
      WRITE(28,'("Error at T=",E9.3,", level=",I1," :",2E10.3,I10)')
     +   T, LEVEL, RMAX1, RMAX2, NPTS

      RETURN
      END
      SUBROUTINE PDEIV (T, X, Y, Z, U, NPTS, NPDE)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPDE
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE)
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
      DOUBLE PRECISION PI, S2, R, S, D, K1, K2, C2, EDKT
      PARAMETER (PI = 3.141592653589793, S2 = 1.414213562373095)
      PARAMETER (K1 = 1000.0, K2 = 1.0)
      PARAMETER (C2 = -K1/K2)
      INTEGER I

      R = (2+SIN(2*PI*T))/4
      S = (2+S2/2*COS(2*PI*T))/4
      DO 10 I = 1, NPTS
         D = EXP(-80*((X(I)-R)**2+(Y(I)-S)**2+(Z(I)-S)**2))
         EDKT = EXP(-D*K2*T)
         U(I,1) = D/(K1+K2) * (K1*(1-C2)+(K1+K2)*C2*EDKT)/(1-C2+C2*EDKT)
         U(I,2) = D - U(I,1)
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
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
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
      DOUBLE PRECISION PI, S2, K1, K2
      PARAMETER (PI = 3.141592653589793, S2 = 1.414213562373095)
      PARAMETER (K1 = 1000.0, K2 = 1.0)
      INTEGER I

      DO 10 I = 1, NPTS
         RES(I,1) = UT(I,1) +
     +              2*PI*S2*((Y(I)+Z(I))/2-0.5)*UX(I,1) -
     +              2*PI/S2*(X(I)-0.5)*(UY(I,1)+UZ(I,1))
     +      - (-K2*U(I,1)*U(I,2) + K1*U(I,2)*U(I,2))
         RES(I,2) = UT(I,2) +
     +              2*PI*S2*((Y(I)+Z(I))/2-0.5)*UX(I,2) -
     +              2*PI/S2*(X(I)-0.5)*(UY(I,2)+UZ(I,2))
     +      - (-K1*U(I,2)*U(I,2) + K2*U(I,1)*U(I,2))
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
      DOUBLE PRECISION T, X(NPTS), Y(NPTS), Z(NPTS), U(NPTS,NPDE),
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
      DOUBLE PRECISION PI, S2, R, S, D, K1, K2, C2, EDKT, UI1
      PARAMETER (PI = 3.141592653589793, S2 = 1.414213562373095)
      PARAMETER (K1 = 1000.0, K2 = 1.0)
      PARAMETER (C2 = -K1/K2)
      INTEGER I, J, K, NBNDS

      R = (2+SIN(2*PI*T))/4
      S = (2+S2/2*COS(2*PI*T))/4
C
      NBNDS = LLBND(0)
      DO 10 J = 1, NBNDS
C
C Change the inflow boundaries
         IF (ILBND(J) .EQ. 1) THEN
            DO 20 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF (Y(I)+Z(I) .GE. 1.0) THEN
                  D = EXP(-80*((X(I)-R)**2+(Y(I)-S)**2+(Z(I)-S)**2))
                  EDKT = EXP(-D*K2*T)
                  UI1 = D/(K1+K2) *
     +                  (K1*(1-C2)+(K1+K2)*C2*EDKT)/(1-C2+C2*EDKT)
                  RES(I,1) = U(I,1) - UI1
                  RES(I,2) = U(I,2) - (D - UI1)
               ENDIF
   20       CONTINUE
         ELSE IF (ILBND(J) .EQ. 2 .OR. ILBND(J) .EQ. 5) THEN
            DO 30 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF (X(I) .LE. 0.5) THEN
                  D = EXP(-80*((X(I)-R)**2+(Y(I)-S)**2+(Z(I)-S)**2))
                  EDKT = EXP(-D*K2*T)
                  UI1 = D/(K1+K2) *
     +                  (K1*(1-C2)+(K1+K2)*C2*EDKT)/(1-C2+C2*EDKT)
                  RES(I,1) = U(I,1) - UI1
                  RES(I,2) = U(I,2) - (D - UI1)
               ENDIF
   30       CONTINUE
         ELSE IF (ILBND(J) .EQ. 3)THEN
            DO 40 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF (Y(I)+Z(I) .LE. 1.0) THEN
                  D = EXP(-80*((X(I)-R)**2+(Y(I)-S)**2+(Z(I)-S)**2))
                  EDKT = EXP(-D*K2*T)
                  UI1 = D/(K1+K2) *
     +                  (K1*(1-C2)+(K1+K2)*C2*EDKT)/(1-C2+C2*EDKT)
                  RES(I,1) = U(I,1) - UI1
                  RES(I,2) = U(I,2) - (D - UI1)
               ENDIF
   40       CONTINUE
         ELSE IF (ILBND(J) .EQ. 4 .OR. ILBND(J) .EQ. 6) THEN
            DO 50 K = LLBND(J), LLBND(J+1)-1
               I = LBND(K)
               IF (X(I) .GE. 0.5) THEN
                  D = EXP(-80*((X(I)-R)**2+(Y(I)-S)**2+(Z(I)-S)**2))
                  EDKT = EXP(-D*K2*T)
                  UI1 = D/(K1+K2) *
     +                  (K1*(1-C2)+(K1+K2)*C2*EDKT)/(1-C2+C2*EDKT)
                  RES(I,1) = U(I,1) - UI1
                  RES(I,2) = U(I,2) - (D - UI1)
               ENDIF
   50       CONTINUE
         ENDIF
   10 CONTINUE

      RETURN
      END
