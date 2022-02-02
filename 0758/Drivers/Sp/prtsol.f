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
      REAL T0, TW, TEW, DTW, XLW, YLW, XRW, YUW, DXB, DYB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0, TW, TEW, DTW, XLW,YLW, XRW,YUW, DXB, DYB, DTO
      SAVE /WRITIF/, /WRITLF/, /WRITRF/
C
C end INCLUDE 'CMNWRITEF'
C
C
C-----------------------------------------------------------------------
C
      INTEGER MXLEV, NPD, NPTS, LENIWK, LENRWK
      PARAMETER (MXLEV=5, NPD=3, NPTS=10000)
      PARAMETER (LENIWK=NPTS*(7*MXLEV+20),
     +           LENRWK=5*NPTS*NPD*MXLEV)
C
      CHARACTER FILE*128
      INTEGER IWK(LENIWK),
     +   LSGNM1, LSGN, LSGNP1, LSUNM1, LSSN, LSUN
      REAL RWK(LENRWK)

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
      CALL PRSOL (61, TW, NPDEW, XLW, YLW, DXB, DYB,
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
      REAL RWK(LENRWK)
C
Ccc PURPOSE:
C Read all information necessary for a restart of VLUGR2 from file
C
Ccc PARAMETER DESCRIPTION:
C LUNDMP : IN.  Logical unit number of dumpfile. Should be opened as an
C          unformatted file.
C RWK    : OUT. Real workstorage intended to pass to VLUGR2
C LENRWK : IN.  Dimension of RWK.
C IWK    : OUT. Integer workstorage intended to pass to VLUGR2
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
      REAL T0, TW, TEW, DTW, XLW, YLW, XRW, YUW, DXB, DYB, DTO
      COMMON /WRITIF/ MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      COMMON /WRITLF/ FIRST, SECOND
      COMMON /WRITRF/ T0, TW, TEW, DTW, XLW,YLW, XRW,YUW, DXB, DYB, DTO
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
     +   T0, TW, TEW, DTW, XLW, YLW, XRW, YUW, DXB, DYB, DTO
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
      SUBROUTINE PRSOL (LUN, T, NPDE, XL, YL, DXB, DYB, LGRID, ISTRUC,
     +   LSOL, SOL)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUN, NPDE, LGRID(0:*), ISTRUC(*), LSOL(*)
      REAL T, XL, YL, DXB, DYB, SOL(*)
C
Ccc PURPOSE:
C Print solution and coordinate values at all grid levels.
C
Ccc PARAMETER DESCRIPTION:
C LUN    : IN.  Logical unit number of print file
C T      : IN.  Current value of time variable
C NPDE   : IN.  # PDE components
C XL     : IN.  X-coordinate of lowerleft corner of (virtual) domain
C YL     : IN.  Y-coordinate of lowerleft corner of (virtual) domain
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
Ccc EXTERNALS USED:
      EXTERNAL PRSOLL
C
C-----------------------------------------------------------------------
C
      INTEGER MAXLEV, LEVEL, LLROW, NROWS, NPTS, LIROW, LICOL
      REAL DX, DY

      MAXLEV = LGRID(0)
      DX = DXB
      DY = DYB
      DO 10 LEVEL = 1, MAXLEV
         LLROW  = LGRID(LEVEL)
         NROWS  = ISTRUC(LLROW)
         NPTS   = ISTRUC(LLROW+NROWS+1)-1
         LIROW  = LLROW+NROWS+2
         LICOL  = LIROW+NROWS
         CALL PRSOLL (LUN, LEVEL, T, NPTS, NPDE, XL, YL, DX, DY,
     +      ISTRUC(LLROW), ISTRUC(LIROW), ISTRUC(LICOL),
     +      SOL(LSOL(LEVEL)+1))
         DX = DX/2
         DY = DY/2
   10 CONTINUE
      RETURN
      END
      SUBROUTINE PRSOLL (LUN, LEVEL, T, NPTS, NPDE, XL, YL, DX, DY,
     +   LROW, IROW, ICOL, U)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUN, LEVEL, NPTS, NPDE, LROW(0:*), IROW(*), ICOL(*)
      REAL T, XL, YL, DX, DY, U(NPTS,NPDE)
C
Ccc PURPOSE:
C Print solution and  X- and Y-coordinates of gridlevel LEVEL.
C
Ccc PARAMETER DESCRIPTION:
C LUN    : IN.  Logical unit number of print file
C LEVEL  : IN.  Grid level corresponding with solution U.
C T      : IN.  Current value of time variable
C NPTS   : IN.  # grid points at this level
C NPDE   : IN.  # PDE components
C XL     : IN.  X-coordinate of lower-left point of virtual rectangle
C YL     : IN.  Y-coordinate of lower-left point of virtual rectangle
C DX     : IN.  Grid width in X-direction
C DY     : IN.  Grid width in Y-direction
C LROW   : IN.  (0:LROW(0)+1)
C          LROW(0) = NROWS: Actual # rows in grid
C          LROW(1:NROWS): pointers to the start of a row in the grid
C          LROW(NROWS+1) = NPTS+1: Actual # nodes in grid + 1
C IROW   : IN.  (NROWS)
C          IROW(IR): row number of row IR in virtual rectangle
C ICOL   : IN.  (NPTS)
C          ICOL(IPT): column number of grid point IPT in virtual
C                     rectangle
C U      : IN.  Solution on this grid level
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, IPT, IR, NROWS
      REAL X, Y
C
      NROWS = LROW(0)

      WRITE(LUN,'(//// T10,A,T30,A,T46,A,T62,A,T71,A //)')
     +   'Level', 't', 'Y', 'X', 'Solution'
      IR = 1
         Y = YL + IROW(IR)*DY
         IPT = LROW(IR)
            X = XL + ICOL(IPT)*DX
            WRITE(LUN,
     +         '(T13,I2,T21,E12.5,T37,E12.5,T53,E12.5,T69,E12.5)')
     +         LEVEL, T, Y, X, U(IPT,1)
            DO 10 IC = 2, NPDE
               WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   10       CONTINUE
         DO 20 IPT = LROW(IR)+1, LROW(IR+1)-1
            X = XL + ICOL(IPT)*DX
            WRITE(LUN,'(T53,E12.5,T69,E12.5)') X, U(IPT,1)
            DO 30 IC = 2, NPDE
               WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   30       CONTINUE
   20    CONTINUE
      DO 40 IR = 2, NROWS
         Y = YL + IROW(IR)*DY
         IPT = LROW(IR)
            X = XL + ICOL(IPT)*DX
            WRITE(LUN,
     +         '(T21,E12.5,T37,E12.5,T53,E12.5,T69,E12.5)')
     +         T, Y, X, U(IPT,1)
            DO 50 IC = 2, NPDE
               WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   50       CONTINUE
         DO 60 IPT = LROW(IR)+1, LROW(IR+1)-1
            X = XL + ICOL(IPT)*DX
            WRITE(LUN,'(T53,E12.5,T69,E12.5)') X, U(IPT,1)
            DO 70 IC = 2, NPDE
               WRITE(LUN,'(T69,E12.5)') U(IPT,IC)
   70       CONTINUE
   60    CONTINUE
   40 CONTINUE

      RETURN
      END
