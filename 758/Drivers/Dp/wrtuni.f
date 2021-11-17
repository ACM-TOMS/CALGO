      PROGRAM WRTUNI
C
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!                                                                !!!
C !!! In subroutine WRUNI the constant NONVAL should be adjusted to  !!!
C !!! the data (NONVAL = impossible value for the first componenent) !!!
C !!!                                                                !!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C-----------------------------------------------------------------------
C
Ccc This program reads a file made by subroutine DUMP and writes the
C (interpolated) solution on a uniform grid of a specified grid level
C to the output file sol.dat. The maximum grid level used in each point
C is written to the file grid.dat.
C NB. This program is not correct for a domain with holes in it with
C a size of the width of the base grid, e.g. it will ignore some holes
C in the domain of the example problem.
C
Ccc EXTERNALS USED:
      EXTERNAL WRUNI, RDDUMP
C
C
Ccc   INCLUDE 'CMNWRITEF'
C
C CMNWRITEF
C
C COMMON needed for continuation calls
      INTEGER MAXLVW, NPDEW, LRWKPS, LIWKPS, LRWKB, LIWKB
      LOGICAL FIRST, SECOND
      DOUBLE PRECISION T0, TW, TEW, DTW, XLW, YLW, XRW, YUW, DXB, DYB,
     +   DTO
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
     +   LSGNM1, LSGN, LSGNP1, LSUNM1, LSSN, LSUN,
     +   LUNI, MAXLEV, NX, NXB, NY, NYB, UNILEV
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
C Check workspace
      MAXLEV = IWK(LSGN)
      PRINT *, 'Max. grid level?'
      READ *, UNILEV
      UNILEV = MIN(UNILEV,MAXLEV)
      NXB = NINT((XRW - XLW)/DXB)
      NYB = NINT((YUW - YLW)/DYB)
      NX = NXB * 2**(UNILEV-1)
      NY = NYB * 2**(UNILEV-1)
      LUNI   = LENRWK - (NX+1)*(NY+1)*NPDEW
      IF (LUNI .LT. IWK(LSUN+MAXLVW)) STOP 'workspace'
C
C Write problem info to standard output and write the interpolated
C solution and grid levels to the files
      PRINT *, 'T, NPDE, XL, YL, DXB, DYB, NXB, NYB'
      PRINT *, TW, NPDEW, XLW, YLW, DXB, DYB, NXB, NYB
      FILE = 'sol.dat'
      OPEN(UNIT=61,FILE=FILE)
      FILE = 'grid.dat'
      OPEN(UNIT=63,FILE=FILE)
      CALL WRUNI (61, 63, UNILEV,
     +   TW, NPDEW, XLW, YLW, DXB, DYB, NXB, NYB,
     +   IWK(LSGN), IWK(LIWKPS), IWK(LSUN), RWK(LRWKPS),
     +   RWK(LUNI), NX, NY)
      CLOSE(61)
      CLOSE(63)
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
      DOUBLE PRECISION T0, TW, TEW, DTW, XLW, YLW, XRW, YUW, DXB, DYB,
     +   DTO
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
      SUBROUTINE WRUNI (LUNS, LUNG, UNILEV,
     +   T, NPDE, XL, YL, DXB, DYB, NXB, NYB,
     +   LGRID, ISTRUC, LSOL, SOL, UNIFRM, NX, NY)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUNS, LUNG, UNILEV,
     +   NPDE, NXB, NYB, LGRID(0:*), ISTRUC(*), LSOL(*), NX, NY
      DOUBLE PRECISION T, XL, YL, DXB, DYB, SOL(*),
     +   UNIFRM(0:NX,0:NY,NPDE)
C
Ccc PURPOSE:
C Write (interpolated) solution values at grid level UNILEV to file
C LUNS.
C Write maximum gridlevel used in each point to file LUNG.
C NB. The data will not be correct for a domain with holes in it with
C a size of the width of the base grid, e.g. it will ignore some holes
C in the domain of the example problem.
C
Ccc PARAMETER DESCRIPTION:
C LUNS   : IN.  Logical unit number of solution file
C LUNG   : IN.  Logical unit number of grid level file
C UNILEV : IN.  Maximum grid level to be used to generate the data
C T      : IN.  Value of time variable
C NPDE   : IN.  # PDE components
C XL     : IN.  X-coordinate of lower left corner of (virtual) domain
C YL     : IN.  Y-coordinate of lower left corner of (virtual) domain
C DXB    : IN.  Cell width in X-direction of base grid
C DYB    : IN.  Cell width in Y-direction of base grid
C NXB,NYB: IN. # gridcells in X- and Y-direction, resp., on base grid
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
C UNIFRM : WORK. (Interpolated) solution on level UNILEV / max. grid
C          level used.
C NX, NY : IN. # gridcells in X- and Y-direction, resp., on grid of
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

      INTEGER I, IC, ICOL, IMUL, IPT, IR, IROW, J,
     +   LEVEL, LLROW, LIROW, LICOL, MAXLEV, NROWS, NPTS

      DO 1 IC = 1, NPDE
      DO 1 IROW = 0, NY
      DO 1 ICOL = 0, NX
         UNIFRM(ICOL,IROW,IC) = NONVAL
    1 CONTINUE

      MAXLEV = LGRID(0)
      DO 10 LEVEL = 1, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         LLROW  = LGRID(LEVEL)
         NROWS  = ISTRUC(LLROW)
         NPTS   = ISTRUC(LLROW+NROWS+1)-1
         LIROW  = LLROW+NROWS+2
         LICOL  = LIROW+NROWS
         DO 20 IR= 1, NROWS
            IROW = ISTRUC(LIROW-1+IR)*IMUL
            DO 30 IPT = ISTRUC(LLROW+IR), ISTRUC(LLROW+IR+1)-1
               ICOL = ISTRUC(LICOL-1+IPT)*IMUL
               DO 40 IC = 1, NPDE
                  UNIFRM(ICOL,IROW,IC) =
     +               SOL(LSOL(LEVEL)+(IC-1)*NPTS+IPT)
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      DO 100 LEVEL = 2, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         DO 110 J = IMUL, NY, IMUL*2
         DO 110 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,1) .EQ. NONVAL) THEN
               DO 120 IC = 1, NPDE
                  UNIFRM(I,J,IC) = 
     +               (UNIFRM(I,J-IMUL,IC)+UNIFRM(I,J+IMUL,IC))/2
  120          CONTINUE
            ENDIF
  110    CONTINUE
         DO 130 J = 0, NY, IMUL
         DO 130 I = IMUL, NX, IMUL*2
            IF (UNIFRM(I,J,1) .EQ. NONVAL) THEN
               DO 140 IC = 1, NPDE
                  UNIFRM(I,J,IC) = 
     +               (UNIFRM(I-IMUL,J,IC)+UNIFRM(I+IMUL,J,IC))/2
  140          CONTINUE
            ENDIF
  130    CONTINUE
  100 CONTINUE

      DO 150 J = 0, NY
      DO 150 I = 0, NX
         WRITE(LUNS,'(100E13.3)') (UNIFRM(I,J,IC), IC = 1, NPDE)
  150 CONTINUE
C
C Grids
      DO 201 IROW = 0, NY
      DO 201 ICOL = 0, NX
         UNIFRM(ICOL,IROW,1) = 0
  201 CONTINUE
      DO 210 LEVEL = 1, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         LLROW  = LGRID(LEVEL)
         NROWS  = ISTRUC(LLROW)
         NPTS   = ISTRUC(LLROW+NROWS+1)-1
         LIROW  = LLROW+NROWS+2
         LICOL  = LIROW+NROWS
         DO 220 IR= 1, NROWS
            IROW = ISTRUC(LIROW-1+IR)*IMUL
            DO 230 IPT = ISTRUC(LLROW+IR), ISTRUC(LLROW+IR+1)-1
               ICOL = ISTRUC(LICOL-1+IPT)*IMUL
               UNIFRM(ICOL,IROW,1) = LEVEL
  230       CONTINUE
  220    CONTINUE
  210 CONTINUE
      DO 300 LEVEL = 2, UNILEV
         IMUL = 2**(UNILEV-LEVEL)
         DO 310 J = IMUL, NY, IMUL*2
         DO 310 I = 0, NX, IMUL*2
            IF (UNIFRM(I,J,1) .LT. LEVEL) THEN
               UNIFRM(I,J,1) =
     +            MIN(UNIFRM(I,J-IMUL,1),UNIFRM(I,J+IMUL,1))
            ENDIF
  310    CONTINUE
         DO 330 J = 0, NY, IMUL
         DO 330 I = IMUL, NX, IMUL*2
            IF (UNIFRM(I,J,1) .LT. LEVEL) THEN
               UNIFRM(I,J,1) =
     +            MIN(UNIFRM(I-IMUL,J,1),UNIFRM(I+IMUL,J,1))
            ENDIF
  330    CONTINUE
  300 CONTINUE

      DO 350 J = 0, NY
      DO 350 I = 0, NX
         WRITE(LUNG,'(I2)') NINT(UNIFRM(I,J,1))
  350 CONTINUE

      RETURN
      END
