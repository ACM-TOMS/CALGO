C
C
C                        CSRFTEST
C                        06/25/02
C
C   This program is a test driver for the software package
C CSRFPACK which constructs a convex once-continuously
C differentiable bivariate function that interpolates a set
C of data values at arbitrarily distributed nodes in the
C plane.
C
C   The user must provide an input data set consisting of
C the number of data points N (Format I5) and a sequence of
C N ordered triples (X,Y,Z) with nodal coordinates (X,Y) and
C data values Z:  Format 3E23.15.
C
C   Provided the data defines a convex function, the output
C consists of the following files:
C
C csrftest.prt:  Print file with parameter values and error
C                messages.
C csrftest.eps:  Encapsulated Postscript file containing a
C                contour plot of the interpolatory surface.
C csrftest.out:  Data set suitable for a surface plotting
C                package consisting of interpolated values
C                on a triangulation of a uniform grid in the
C                bounding box (the smallest rectangle with
C                sides parallel to the axes that contains
C                the nodes).
C tplot.eps:  Encapsulated Postscript file containing the
C             convexity-preserving triangulation of the
C             nodes.
C gplot.eps:  Encapsulated Postscript file containing the
C             gradient feasibility diagram and nodal
C             gradients
C tgplot.eps:  Encapsulated Postscript file containing the
C              gradient triangulation.
C rplot.eps:  Encapsulated Postscript file containing the
C             cell diagram and nodes.
C
C Previously existing files with the above names are over-
C written.
C
      INTEGER LC, LWK, N2, N6, NA, NMAX, NR, NX, NY
C
C Maximum number of nodes NMAX:
C
      PARAMETER (NMAX=300, N2=2*NMAX, N6=6*NMAX)
C
C Number of grid points NX (horizontally) and NY
C   (vertically) defining the uniform grid.
C
      PARAMETER (NX=40, NY=40, LWK=NX*NY+(NX*NY)/2,
     .           LC=NX*NY+LWK)
C
C Number of quadrature points in the radial and angular
C   directions:
C
      PARAMETER (NR=8, NA=24)
C
      CHARACTER*60 FNAME
      DOUBLE PRECISION C(NMAX), D(NMAX), DMIN, DX, DXL(N2),
     .                 DY, DYL(N2), EPS, F(NX,NY), GX(NMAX),
     .                 GY(NMAX), PLTSIZ, PX(NX), PY(NY),
     .                 QX(NR,NA), QY(NR,NA), TOLBE, W(NR),
     .                 X(NMAX), XC(LC), XMAX, XMIN, Y(NMAX),
     .                 YC(LC), YMAX, YMIN, Z(NMAX)
      INTEGER I, IER, IPLOT, IWK(LWK), J, K, LCON,
     .        LEND(NMAX), LIN, LIST(N6), LISTV(N6), LOUT,
     .        LPRT, LPTR(N6), LNEW, N, NCON, ND, NEAR(NMAX),
     .        NEXT(NMAX), NV
      LOGICAL STRICT
C
C Data:
C
C Plot option in the range 0 to 15:  IPLOT = (b3,b2,b1,b0),
C   with bits specifying which Postscript files are created:
C
C   b0:  tplot.eps
C   b1:  gplot.eps
C   b2:  tgplot.eps
C   b3:  rplot.eps
C
      DATA IPLOT/15/
C
C Number of contour values in the contour plot.
C
      DATA NCON/10/
C
C PostScript plot size:  1.0 <= PLTSIZ <= 7.5.
C
      DATA PLTSIZ/3.0/
C
C Strict convexity option (quadratic correction term
C   requiring an O(N**2) algorithm).
C
      DATA STRICT/.TRUE./
C
C Positive tolerance for removing boundary triangles.
C
      DATA TOLBE/1.D-2/
C
C LIN,LOUT = Logical unit numbers for input and output.
C
      DATA LCON/1/, LIN/5/, LOUT/6/, LPRT/4/
C
C Input/output formats:
C
  300 FORMAT (I5)
  310 FORMAT (3E23.15)
  320 FORMAT (A60)
  330 FORMAT (3I6)
C
C Get an input file name and open it.
C
C   1 WRITE (*,100)
C 100 FORMAT (///13X,'CSRFTEST:  CSRFPACK Test Program'//
C    .        5X,'Specify a data set file name (at most 60',
C    .           ' characters):'/)
C     READ (*,320,ERR=1) FNAME
C     OPEN (LIN,FILE=FNAME,STATUS='OLD',ERR=1)
C
C Read N and the nodal coordinates and data values from
C   unit LIN.
C
      READ (LIN,300) N
      IF (N .LT. 3  .OR.  N .GT. NMAX) THEN
        WRITE (*,110) NMAX
  110   FORMAT (//5X,'Invalid data.  N must be in the ',
     .               'range 3 to ',I5,'.')
        STOP
      ENDIF
      READ (LIN,310) (X(I),Y(I),Z(I), I = 1,N)
C
C Compute the bounding box corner coordinates (XMIN,YMIN)
C   and (XMAX,YMAX) and the uniform grid points
C   (((PX(i),PY(j)), i = 1,NX), j = 1,NY).
C
      XMIN = X(1)
      XMAX = X(1)
      YMIN = Y(1)
      YMAX = Y(1)
      DO 2 I = 2,N
        IF (X(I) .LT. XMIN) XMIN = X(I)
        IF (X(I) .GT. XMAX) XMAX = X(I)
        IF (Y(I) .LT. YMIN) YMIN = Y(I)
        IF (Y(I) .GT. YMAX) YMAX = Y(I)
    2   CONTINUE
      DX = (XMAX-XMIN)/DBLE(NX-1)
      DY = (YMAX-YMIN)/DBLE(NY-1)
      DO 3 I = 1,NX
        PX(I) = XMIN + DBLE(I-1)*DX
    3   CONTINUE
      DO 4 J = 1,NY
        PY(J) = YMIN + DBLE(J-1)*DY
    4   CONTINUE
C
C Create a print file with parameter values.
C
      OPEN (LPRT,FILE='csrftest.prt')
      WRITE (LPRT,120) N, XMIN, XMAX, YMIN, YMAX, TOLBE,
     .                 STRICT, NR, NA
  120 FORMAT (///27X,'CSRFTEST.PRT'//
     .        10X,'N = ',I5/
     .        10X,'XMIN = ',D10.3,3X,'XMAX = ',D10.3/
     .        10X,'YMIN = ',D10.3,3X,'YMAX = ',D10.3/
     .        10X,'TOLBE = ',D10.3/
     .        10X,'STRICT = ',L1/
     .        10X,'NR = ',I2/
     .        10X,'NA = ',I2)
C
C Construct the interpolatory surface F.
C
      CALL CSURF (N,X,Y,TOLBE,IPLOT,PLTSIZ,STRICT,NR,NA, ND,
     .            Z,C,LIST,LPTR,LEND,LNEW,NEAR,NEXT,NV,
     .            LISTV,DXL,DYL,GX,GY,EPS,D,DMIN,W,QX,QY,
     .            IER)
      IF (IER .NE. 0) THEN
        IF (IER .NE. -11) GO TO 10
C
C DMIN = 0.
C
        WRITE (*,130)
        WRITE (LPRT,130)
      ENDIF
  130 FORMAT (//10X,'*** The surface is not strictly ',
     .              'convex:  DMIN = 0. ***'//)
C
C Output computed parameter values to the print file.
C
      WRITE (LPRT,140) ND, EPS, DMIN
      WRITE (*,140) ND, EPS, DMIN
  140 FORMAT (//10X,'Number of boundary edge deletions:  ',
     .              'ND = ',I3/
     .          10X,'Scale factor for quadratic term:  ',
     .              'EPS = ',D10.3/
     .          10X,'Min. distance node to cell boundary: ',
     .              ' DMIN = ',D10.3)
C
C Compute interpolated values on the uniform grid.
C
      CALL FGRID (NX,NY,PX,PY,EPS,GX,GY,C,LIST,LPTR,
     .            LEND,DMIN,NR,NA,W,QX,QY, F,IER)
C
C Create a contour plot of the interpolatory surface on the
C   NX by NY uniform grid.
C
      OPEN (LCON,FILE='csrftest.eps')
      CALL PLTCNT (LCON,PLTSIZ,NX,NY,PX,PY,F,NCON,IWK,XC,
     .             YC, IER)
C
C Create the output data set using the natural ordering of
C   the rectangular grid points:  left-to-right within
C   bottom-to-top.
C
C Write the number of vertices followed by the sequence of
C   vertex coordinate triples.
C
C     OPEN (LOUT,FILE='csrftest.out')
      WRITE (LOUT,300) NX*NY
      DO 6 J = 1,NY
        DO 5 I = 1,NX
          WRITE (LOUT,310) PX(I), PY(J), F(I,J)
    5     CONTINUE
    6   CONTINUE
C
C Write the number of triangles followed by the sequence of
C   triangles (CCW-ordered vertex indices).  For each of the
C   (NX-1)*(NY-1) cells, k is the vertex index of the upper
C   right corner, and the cell is partitioned by the
C   diagonal with slope 1.
C
      WRITE (LOUT,300) 2*(NX-1)*(NY-1)
      DO 8 J = 2,NY
        DO 7 I = 2,NX
          K = NX*(J-1) + I
          WRITE (LOUT,330) K-NX-1, K-NX, K
          WRITE (LOUT,330) K-NX-1, K, K-1
    7     CONTINUE
    8   CONTINUE
      STOP
C
C Error in CSURF.
C
   10 WRITE (LPRT,400) IER
      WRITE (*,400) IER
  400 FORMAT (///10X,'Error flag ',I4,' returned by CSURF.')
      IF (IER .EQ. -2) THEN
        WRITE (LPRT,410)
        WRITE (*,410)
      ENDIF
  410 FORMAT (/10X,'The first three nodes are collinear.')
      IF (IER .EQ. -3) THEN
        WRITE (LPRT,420)
        WRITE (*,420)
      ENDIF
  420 FORMAT (/10X,'A convex triangulation does not exist ',
     .             '(ADDNDC error -3).')
      IF (IER .EQ. -4) THEN
        WRITE (LPRT,430)
        WRITE (*,430)
      ENDIF
  430 FORMAT (/10X,'A convex triangulation does not exist ',
     .             '(ADDNDC error -4).')
      IF (IER .GT. 0) THEN
        WRITE (LPRT,440)
        WRITE (*,440)
      ENDIF
  440 FORMAT (/10X,'Duplicate nodes.')
      STOP
      END
