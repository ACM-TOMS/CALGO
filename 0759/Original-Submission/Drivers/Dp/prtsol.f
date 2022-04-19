      PROGRAM PRTSOL
C
C-----------------------------------------------------------------------
C
Ccc This program reads a file made by subroutine DUMP and prints the
C solution on an output file. Both filenames are read from standard
C input.
C
Ccc EXTERNALS USED:
      EXTERNAL PRSOL, RDDUMP
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
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK
      PARAMETER (MXLEV=3, NPD=1, NPTS=100000)
      PARAMETER (LENIWK=NPTS*(7*MXLEV+23),
     +           LENRWK=5*NPTS*NPD*MXLEV)
C
      CHARACTER FILE*128
      INTEGER IWK(LENIWK),
     +   LSGNM1, LSGN, LSGNP1, LSUNM1, LSSN, LSUN
      DOUBLE PRECISION RWK(LENRWK)

      PRINT *, 'DUMP file?'
      READ '(A)', FILE
C 
      OPEN(UNIT=62,FILE=FILE,FORM='UNFORMATTED')
      CALL RDDUMP (62, RWK, LENRWK, IWK, LENIWK)
      CLOSE(62)
C 
C Setup work storage
      LSGNM1 = 1
      LSGN   = LSGNM1 + MAXLVW+1
      LSGNP1 = LSGN   + MAXLVW+1
      LSUNM1 = LSGNP1 + MAXLVW+1
      LSSN   = LSUNM1 + MAXLVW
      LSUN   = LSSN   + MAXLVW
C
C call print routine
      PRINT *, 'output file?'
      READ '(A)', FILE
C 
      OPEN(UNIT=61,FILE=FILE)
      CALL PRSOL (61, TW, NPDEW, XLW, YFW, ZDW, DXB, DYB, DZB,
     +   IWK(LSGN), IWK(LIWKPS), IWK(LSUN), RWK(LRWKPS))
      CLOSE(61)
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
