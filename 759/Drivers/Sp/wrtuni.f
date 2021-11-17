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
C a size of the width of the base grid.
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
      REAL T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,DZB, DTO
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
     +   LSGNM1, LSGN, LSGNP1, LSUNM1, LSSN, LSUN,
     +   LUNI, MAXLEV, NX, NXB, NY, NYB, NZ, NZB, UNILEV
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
C Check workspace
      MAXLEV = IWK(LSGN)
      PRINT *, 'Max. grid level?'
      READ *, UNILEV
      UNILEV = MIN(UNILEV,MAXLEV)
      NXB = NINT((XRW - XLW)/DXB)
      NYB = NINT((YBW - YFW)/DYB)
      NZB = NINT((ZUW - ZDW)/DZB)
      NX = NXB * 2**(UNILEV-1)
      NY = NYB * 2**(UNILEV-1)
      NZ = NZB * 2**(UNILEV-1)
      LUNI   = LENRWK - (NX+1)*(NY+1)*(NZ+1)*NPDEW
      IF (LUNI .LT. IWK(LSUN+MAXLVW)) STOP 'workspace'
C
C Write problem info to standard output and write the interpolated
C solution and grid levels to the files
      PRINT *, 'T, NPDE, XL, YF, ZD, DXB, DYB, DZB, NXB, NYB, NZB'
      PRINT *, TW, NPDEW, XLW, YFW, ZDW, DXB, DYB, DZB, NXB, NYB, NZB
      FILE = 'sol.dat'
      OPEN(UNIT=61,FILE=FILE)
      FILE = 'grid.dat'
      OPEN(UNIT=63,FILE=FILE)
      CALL WRUNI (61, 63, UNILEV,
     +   TW, NPDEW, XLW, YFW, ZDW, DXB, DYB, DZB, NXB, NYB, NZB,
     +   IWK(LSGN), IWK(LIWKPS), IWK(LSUN), RWK(LRWKPS),
     +   RWK(LUNI), NX, NY, NZ)
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
      REAL RWK(LENRWK)
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
      REAL T0,TW,TEW,DTW, XLW,YFW,ZDW, XRW,YBW,ZUW, DXB,DYB,DZB, DTO
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
      SUBROUTINE WRUNI (LUNS, LUNG, UNILEV,
     +   T, NPDE, XL, YF, ZD, DXB, DYB, DZB, NXB, NYB, NZB,
     +   LGRID, ISTRUC, LSOL, SOL, UNIFRM, NX, NY, NZ)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER LUNS, LUNG, UNILEV,
     +   NPDE, NXB, NYB, NZB, LGRID(0:*), ISTRUC(*), LSOL(*), NX, NY, NZ
      REAL T, XL, YF, ZD, DXB, DYB, DZB, SOL(*),
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
      REAL NONVAL
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
