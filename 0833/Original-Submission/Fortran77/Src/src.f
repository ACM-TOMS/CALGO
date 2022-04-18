      SUBROUTINE CSURF (N,X,Y,TOLBE,IPLOT,PLTSIZ,STRICT,NR,
     .                  NA, ND,Z,C,LIST,LPTR,LEND,LNEW,NEAR,
     .                  NEXT,NV,LISTV,DXL,DYL,GX,GY,EPS,D,
     .                  DMIN,W,QX,QY,IER)
      INTEGER N, IPLOT, NR, NA, ND, LIST(*), LPTR(*),
     .        LEND(N), LNEW, NEAR(N), NEXT(N), NV, LISTV(*),
     .        IER
      LOGICAL STRICT
      DOUBLE PRECISION X(N), Y(N), TOLBE, PLTSIZ, Z(N),
     .                 C(N), DXL(*), DYL(*), GX(N), GY(N),
     .                 EPS, D(N), DMIN, W(NR), QX(NR,NA),
     .                 QY(NR,NA)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/98
C
C   This subroutine constructs a once-continuously differ-
C entiable convex bivariate function F that interpolates a
C set of N strictly convex data points (data points for
C which the piecewise linear interpolant on some triangula-
C tion is convex, and no four points are coplanar).  Refer
C to Subroutines FVAL and FGRID for a means of evaluating F.
C
C   The interpolant F is obtained by applying convolution
C smoothing to a piecewise linear Hermite interpolant H of
C the data values and a set of convexity-preserving nodal
C gradients.  The procedure consists of the following ten
C steps:
C
C 1)  Construct a convexity-preserving triangulation T, if
C       it exists, of the data points (nodes and data
C       values):  Subroutine TRMSHC.  Remove extraneous
C       boundary edges (nearly null triangles):  Subroutine
C       DELBE.
C 2)  Construct a gradient feasibility diagram (the straight
C       line dual of T):  a set of N convex polygons that
C       partition the plane and such that, for any choice of
C       nodal gradients, one from each polygon, the piece-
C       wise linear Hermite interpolant H of the nodal
C       values and gradients is convex:  Subroutine VLIST.
C 3)  Select a set of nodal gradients by taking the cen-
C       troids of the gradient feasibility regions or (in
C       the case of boundary nodes) truncated regions:  Sub-
C       routine GLIST.
C 4)  Optionally, create level-2 Encapsulated PostScript
C       files containing plots of T and/or the dual of T
C       (along with the nodal gradients):  Subroutines PLTTR
C       and PLTGR.
C 5)  Optionally, compute a scale factor EPS defining a
C       quadratic function q(p) = EPS*<p,p>, and adjust the
C       data values and gradients by subtracting nodal
C       values and gradients of q, where EPS is defined so
C       that a strictly convex Hermite interpolant of the
C       original data can be constructed from a convex Her-
C       mite interpolant of the adjusted data (Subroutine
C       ADDQT).
C 6)  Construct a convexity-preserving triangulation Tg of
C       the (possibly adjusted) nodal gradients and negative
C       z-intercepts of the affine nodal functions defined
C       by the (adjusted) data values and nodal gradients:
C       Subroutine TRMSHC.
C 7)  Construct a cell diagram (the straight line dual of
C       Tg):  a set of N convex polygons R_i that partition
C       the plane and such that the restriction of H to R_i
C       is the i-th affine nodal function f_i, where f_i(p)
C       = <g_i,p-p_i> + z_i for nodes p_i, (adjusted) data
C       values z_i and (adjusted) gradients g_i for i = 1 to
C       N and H(p) = Max(i=1,N){f_i(p)}:  Subroutine VLIST.
C 8)  Optionally, create level-2 Encapsulated PostScript
C       files containing plots of Tg and/or the dual of Tg
C       (along with the nodes):  Subroutines PLTTR and
C       PLTGR.
C 9)  Compute the distance D_i from each node p_i to the
C       nearest boundary point of cell R_i:  Subroutine
C       DELTAI.
C 10) Compute weights W_i and quadrature points q_ij de-
C       fining the quadrature rule
C         Q(p) = Sum(i=1,NR) [ W_i*Sum(j=1,NA) H(p+q_ij) ]
C       for approximating F(p) = Integral[H(p+q)*Phi(q)]dq,
C       where the integral is over the disk D0 of radius
C       Dmin = Min(i=1,N){D_i} centered at the origin, and
C       Phi(q) = phi(Norm(q)/Dmin), normalized to have
C       integral 1, for phi(t) = 1 - 3*t**2 + 2*t**3 (t in
C       [0,1]):  Subroutine GETQW.
C
C   The following is a list of related subprograms which a
C user may wish to call directly:
C
C  ADDQT  - Given a set of strictly convex Hermite data
C             (nodes, data values, and nodal gradients such
C             that there exists a convex Hermite interpolant
C             of the data values and gradients, and no four
C             data points are coplanar), adjusts the data
C             values and gradients such that a strictly con-
C             vex Hermite interpolant of the data can be
C             constructed from a convex Hermite interpolant
C             of the adjusted data.
C
C  PLTCNT - Given a set of function values Z = f(x,y) at the
C             vertices of a rectangular grid, creates a
C             level-2 Encapsulated PostScript plot contain-
C             ing a contour plot of the piecewise bilinear
C             interpolant of the function values.  Refer to
C             Subroutine FGRID.
C
C  DELTAI - Given a set of nodes p_i, data values z_i, and
C             nodal gradients g_i for which there exists a
C             convex Hermite interpolatory surface, let R_i
C             denote the set of points p such that f_i(p) >=
C             f_j(p) for all j, where f_i is the affine
C             nodal function that interpolates the Hermite
C             data:  f_i(p) = <g_i,p-p_i> + z_i.  DELTAI
C             returns the perpendicular distance from each
C             node p_i to the boundary of R_i.
C
C  FGRID  - Evaluates F at the vertices of a rectangular
C             grid.
C
C  FVAL   - Given a point p and the output parameters from
C             Subroutine CSURF, returns the value and gradi-
C             ent of F at p.
C
C  GETQW  - Computes weights and quadrature points defining
C             a rule for approximating values of F.
C
C  GLIST  - Given a convexity-preserving triangulation and
C             its straight-line dual computed by VLIST, re-
C             turns a set of nodal gradients for which there
C             exists a convex Hermite interpolant of the
C             data values and gradients.  The nodal gradi-
C             ents are taken to be the centroids of the
C             gradient feasibility regions (or truncated
C             regions).
C
C  LGRAD  - Computes the gradient of the linear interpolant
C             of data values at the vertices of a user-
C             specified triangle.
C
C  PLTGR - Creates a level-2 Encapsulated PostScript plot
C            of the straight-line dual of a triangulation
C            along with a set of nodal gradients.
C
C  PLTTR  - Creates a level-2 Encapsulated PostScript plot
C             of a triangulation.
C
C  TRMSHC - Constructs a convexity-preserving triangulation,
C             if it exists, of a set of N arbitrarily dis-
C             tributed points in the plane (referred to as
C             nodes) with associated data values.
C
C  VLIST  - Given a convexity-preserving triangulation, com-
C             putes its straight-line dual:  the set of all
C             nodal gradients for which there exists a con-
C             vex Hermite interpolant of the data values and
C             gradients.
C
C   Refer to the header comments in Subroutine TRMSHC for a
C list of additional user-callable subprograms related to
C the triangulations.  The remaining subprograms in this
C package appear in alphabetical order.
C
C
C On input:
C
C       N = Number of nodes and data values.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.  (X(K),Y(K)) is re-
C             ferred to as node K, and K is referred to as
C             a nodal index.  The first three nodes must not
C             be collinear.
C
C       TOLBE = Positive tolerance used (in Subroutine
C               DELBE) as a lower bound on the acceptable
C               aspect ratio (ratio of the radius of the
C               inscribed circle to that of the circum-
C               circle) of boundary triangles.  This should
C               be at least 1.E-6 in order to avoid failure
C               in LGRAD (called by VLIST) and should be
C               larger (1.E-2) to avoid large errors caused
C               by steep gradients on the boundary.  A large
C               value, on the other hand, may lead to fail-
C               ure with IER = -7 or -8 if the data is
C               convex but not strictly convex.
C
C       IPLOT = Plot option in the range 0 to 15.  Denoting
C               the sequence of bits (high order to low
C               order) by (b3,b2,b1,b0), each bit specifies
C               an option as follows:
C                 b0 = 1 iff the triangulation T is to be
C                        plotted (tplot.eps),
C                 b1 = 1 iff the gradient feasibility dia-
C                        gram with the nodal gradients is to
C                        be plotted (gplot.eps),
C                 b2 = 1 iff the gradient triangulation Tg
C                        is to be plotted (tgplot.eps),
C                 b3 = 1 iff the cell diagram with the nodes
C                        is to be plotted (rplot.eps).
C               The file names are shown in parentheses for
C               each option.  If a requested plot is not
C               created (and no other error occurred) it is
C               because an error occurred in opening or
C               writing to the plot file.  No error flag is
C               returned in this case.
C
C       PLTSIZ = Plot size in inches (or dummy parameter if
C                IPLOT = 0).  A window containing the data
C                is mapped, with aspect ratio preserved, to
C                a rectangular viewport with maximum side-
C                length equal to .88*PLTSIZ (leaving room
C                for labels outside the viewport).  The
C                viewport is centered on the 8.5 by 11 inch
C                page, and its boundary is drawn.  1.0 .LE.
C                PLTSIZ .LE. 8.5.
C
C       STRICT = Logical variable with value TRUE iff the
C                data is to be adjusted so that a strictly
C                convex interpolant can be computed.  This
C                requires an O(N**2) algorithm to compute
C                EPS and dominates the cost of constructing
C                the interpolant.
C
C       NR = Number of intervals in the radial direction
C            defining a partitioning of the disk D0, and
C            number of quadrature weights.  NR = 8 should
C            be sufficient.  NR > 0.
C
C       NA = Number of intervals in the angular direction
C            defining the partitioning of D0.  NA = 3*NR is
C            recommended.  NA > 1.
C
C The above parameters are not altered by this routine.
C
C       Z = Array of length N containing the data values.
C
C       C = Array of length at least N.
C
C       LIST,LPTR = Integer arrays of length at least 6N-12.
C
C       LEND,NEAR,NEXT = Integer arrays of length at least
C                        N.  NEAR and NEXT are work space
C                        arrays used by Subroutine TRMSHC.
C
C       LISTV = Integer array of length at least 6*N-12.
C
C       DXL,DYL = Arrays of length at least 2*N-2.
C
C       GX,GY,D = Arrays of length at least N.
C
C       W = Array of length at least NR.
C
C       QX,QY = Arrays of length at least NR*NA.
C
C On output:
C
C       ND = Number of boundary edges deleted by DELBE.
C
C       Z = Array of data values or, if STRICT = TRUE,
C           adjusted data values.
C
C       C = Array containing negative z-intercepts (constant
C           terms) of the affine nodal functions:  C(i) =
C           <g_i,p_i> - z_i for i = 1 to N.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             gradient triangulation Tg.
C                             Refer to Subroutine TRMSHC.
C
C       NEAR,NEXT = Garbage.
C
C       NV = Number of cell diagram vertices, including a
C            pseudo-vertex on each boundary edge (used to
C            truncate the unbounded regions):  NV = NT+NB =
C            2*N-2, where NB is the number of boundary nodes
C            and NT = 2*N-NB-2 is the number of triangles.
C
C       LISTV = Array containing vertex indexes (indexes to
C               DXL and DYL) stored in one-to-one corres-
C               pondence with LIST/LPTR entries.  Refer to
C               Subroutine VLIST.
C
C       DXL,DYL = Arrays of length NV containing the ver-
C                 tices of the cell diagram in the first NT
C                 positions and pseudo-vertices in the last
C                 NB positions.
C
C       GX,GY = Arrays containing the components of the
C               nodal gradients or, if STRICT = TRUE, the
C               adjusted nodal gradients.
C
C       EPS = Scale factor defining the quadratic term (0 if
C             STRICT = FALSE).
C
C       D = Array containing the distance from each node to
C           its cell boundary.
C
C       DMIN = Smallest element of D.
C
C       W = Array containing quadrature weights W_i.
C
C       QX,QY = Arrays dimensioned NR by NA containing the
C               components of the quadrature points:  q_ij =
C               (QX(i,j),QY(i,j)) for i = 1 to NR and j = 1
C               to NA.
C
C       IER = Error indicator:
C             IER = 0 if no error was encountered.
C             IER = -1 if N, IPLOT, PLTSIZ, NR, or NA is
C                      outside its valid range on input.
C             IER = -2 if the first three nodes are
C                      collinear.
C             IER = -3 if a convex triangulation does not
C                      exist (and Subroutine ADDNDC returned
C                      IER = -3).
C             IER = -4 if a convex triangulation does not
C                      exist (and Subroutine ADDNDC returned
C                      IER = -4).
C             IER =  L if nodes L and M coincide for some
C                      M > L.  The linked list represents
C                      a triangulation of nodes 1 to M-1
C                      in this case.
C             IER = -5 if a nonpositive triangle area was
C                      returned by Subroutine LGRAD (called
C                      by VLIST) implying an invalid triang-
C                      ulation T.
C             IER = -6 if the first three nodal gradients
C                      are collinear.
C             IER = -7 if a convex triangulation of the
C                      nodal gradients does not exist (and
C                      Subroutine ADDNDC returned IER = -3).
C             IER = -8 if a convex triangulation of the
C                      nodal gradients does not exist (and
C                      Subroutine ADDNDC returned IER = -4).
C             IER = -9 if a nonpositive triangle area was
C                      returned by Subroutine LGRAD (called
C                      by VLIST) implying an invalid triang-
C                      ulation Tg.
C             IER = -10 if the Hermite data is not convex
C                       (Subroutine DELTAI returned IER =
C                       3).
C             IER = -11 if the Hermite data is not strictly
C                       convex (D includes a zero entry and
C                       the domain D0 therefore has radius
C                       0).
C
C   The error conditions are tested in the order specified.
C Error flags -5, -6, ..., -10 are unlikely but could occur
C as a result of floating-point inaccuracy.  Error flag -11
C may occur if the data is convex but not strictly convex.
C W, QX, and QY are not altered in this case, but F is well-
C defined (it coincides with H) and may be evaluated by FVAL
C or FGRID.  Also, note the comments regarding error flags
C -7 and -8 in the description of TOLBE above.
C
C Modules required by CSURF:  ADDNDC, ADDQT, BDYADD, CIRCUM,
C                               DELARC, DELBE, DELNB,
C                               DELTAI, DSTORE, GETQW,
C                               GLIST, INSERT, INTADD,
C                               JRAND, LEFT, LGRAD, LSTPTR,
C                               PLTGR, PLTTR, SWAP, SWPTC,
C                               TRFIND, TRMSHC, VLIST
C
C***********************************************************
C
      CHARACTER*80 TITLE
      DOUBLE PRECISION DX, DY, WX1, WX2, WY1, WY2
      INTEGER I, IERR, IPLT, LU0, LU1, LU2, LU3
      LOGICAL NUMBR
C
C Local parameters:
C
C I =          Nodal index
C IERR =       Error flag for calls to VLIST
C IPLT =       Local copy of IPLOT (right-shifted 0 to 3
C                times)
C LU0,..,LU3 = Logical unit numbers for PostScript plots
C NUMBR =      Logical variable with value TRUE iff nodal
C                indexes are to be displayed on the plots
C WX1,WY1 =    Lower left corner of a window against which
C                a plot is clipped
C WX2,WY2 =    Upper rigjt corner of a window against which
C                a plot is clipped
C TITLE =      Title for the PostScript plots -- must be
C                enclosed in parentheses
C
      DATA LU0/70/,  LU1/71/,  LU2/72/,  LU3/73/
      NUMBR = N .LE. 100
C
C Test for an invalid input parameter.
C
      IER = -1
      IF (N .LT. 3  .OR.  IPLOT .LT. 0  .OR.  IPLOT .GT. 15
     .    .OR.  PLTSIZ .LT. 1.D0  .OR.  PLTSIZ .GT. 8.5D0
     .    .OR.  NR .LT. 1  .OR.  NA .LT. 2) RETURN
C
C Construct the convexity-preserving triangulation T.
C
      CALL TRMSHC (N,X,Y,Z, LIST,LPTR,LEND,LNEW,NEAR,NEXT,
     .             IER)
      IF (IER .NE. 0) RETURN
C
C Remove nearly null triangles from the boundary.  These
C   could result in very innaccurate vertices in the
C   gradient feasibility diagram.
C
      CALL DELBE (N,X,Y,TOLBE, LIST,LPTR,LEND,LNEW, ND)
C
C Construct the gradient feasibility diagram (the straight
C   line dual of T).
C
      CALL VLIST (N,X,Y,Z,LIST,LPTR,LEND, NV,LISTV,DXL,DYL,
     .            IERR)
      IF (IERR .NE. 0) THEN
        IER = -5
        RETURN
      ENDIF
C
C Select a set of nodal gradients by taking the centroids of
C   the gradient feasibility regions or (in the case of
C   boundary nodes) truncated regions.
C
      CALL GLIST (N,LIST,LPTR,LEND,LISTV,DXL,DYL, GX,GY,
     .            IERR)
C
C Plot T and/or its dual (along with the nodal gradients).
C
      IPLT = IPLOT
      IF (IPLT .NE. 2*(IPLT/2)) THEN
        OPEN (LU0,FILE='tplot.eps')
        WX1 = X(1)
        WX2 = WX1
        WY1 = Y(1)
        WY2 = WY1
        DO 1 I = 2,N
          IF (X(I) .LT. WX1) WX1 = X(I)
          IF (X(I) .GT. WX2) WX2 = X(I)
          IF (Y(I) .LT. WY1) WY1 = Y(I)
          IF (Y(I) .GT. WY2) WY2 = Y(I)
    1     CONTINUE
        DX = WX2-WX1
        DY = WY2-WY1
        WX1 = WX1 - 0.1D0*DX
        WX2 = WX2 + 0.1D0*DX
        WY1 = WY1 - 0.1D0*DY
        WY2 = WY2 + 0.1D0*DY
        TITLE = '()'
        CALL PLTTR (LU0,PLTSIZ,WX1,WX2,WY1,WY2,N,X,Y,LIST,
     .              LPTR,LEND,TITLE,NUMBR, IERR)
      ENDIF
      IPLT = IPLT/2
      IF (IPLT .NE. 2*(IPLT/2)) THEN
        OPEN (LU1,FILE='gplot.eps')
        WX1 = DXL(1)
        WX2 = WX1
        WY1 = DYL(1)
        WY2 = WY1
        DO 2 I = 2,NV
          IF (DXL(I) .LT. WX1) WX1 = DXL(I)
          IF (DXL(I) .GT. WX2) WX2 = DXL(I)
          IF (DYL(I) .LT. WY1) WY1 = DYL(I)
          IF (DYL(I) .GT. WY2) WY2 = DYL(I)
    2     CONTINUE
        DX = WX2-WX1
        DY = WY2-WY1
        WX1 = WX1 - 0.1D0*DX
        WX2 = WX2 + 0.1D0*DX
        WY1 = WY1 - 0.1D0*DY
        WY2 = WY2 + 0.1D0*DY
        TITLE = '()'
        CALL PLTGR (LU1,PLTSIZ,WX1,WX2,WY1,WY2,N,LIST,LPTR,
     .              LEND,NV,LISTV,DXL,DYL,GX,GY,TITLE,
     .              NUMBR, IERR)
      ENDIF
C
C Compute the scale factor EPS for the quadratic function
C   q(p) = EPS*<p,p>, and adjust the data values and
C   gradients if STRICT = TRUE.
C
      IF (STRICT) THEN
        CALL ADDQT (N,X,Y, Z,GX,GY, EPS,IERR)
      ELSE
        EPS = 0.
      ENDIF
C
C Compute the negative z-intercepts (constant terms)
C   of the affine nodal functions:  C(i) = <g_i,p_i> - z_i
C   for i = 1 to N.
C
      DO 3 I = 1,N
        C(I) = GX(I)*X(I) + GY(I)*Y(I) - Z(I)
    3   CONTINUE
C
C Construct the convexity-preserving triangulation Tg of
C   the nodal gradients and negative z-intercepts.  The
C   nodal triangulation data structure is overwritten.
C
      CALL TRMSHC (N,GX,GY,C, LIST,LPTR,LEND,LNEW,NEAR,
     .             NEXT,IERR)
      IF (IERR .NE. 0) THEN
        IER = -12
        IF (IERR .EQ. -2) IER = -6
        IF (IERR .EQ. -3) IER = -7
        IF (IERR .EQ. -4) IER = -8
        RETURN
      ENDIF
C
C Construct the cell diagram (the straight line dual of Tg).
C   The data structure defining the dual of T is over-
C   written.
C
      CALL VLIST (N,GX,GY,C,LIST,LPTR,LEND, NV,LISTV,DXL,
     .            DYL,IERR)
      IF (IERR .NE. 0) THEN
        IER = -9
        RETURN
      ENDIF
C
C Plot Tg and/or its dual (along with the nodes).
C
      IPLT = IPLT/2
      IF (IPLT .NE. 2*(IPLT/2)) THEN
        OPEN (LU2,FILE='tgplot.eps')
        WX1 = GX(1)
        WX2 = WX1
        WY1 = GY(1)
        WY2 = WY1
        DO 4 I = 2,N
          IF (GX(I) .LT. WX1) WX1 = GX(I)
          IF (GX(I) .GT. WX2) WX2 = GX(I)
          IF (GY(I) .LT. WY1) WY1 = GY(I)
          IF (GY(I) .GT. WY2) WY2 = GY(I)
    4     CONTINUE
        DX = WX2-WX1
        DY = WY2-WY1
        WX1 = WX1 - 0.1D0*DX
        WX2 = WX2 + 0.1D0*DX
        WY1 = WY1 - 0.1D0*DY
        WY2 = WY2 + 0.1D0*DY
        TITLE = '()'
        CALL PLTTR (LU2,PLTSIZ,WX1,WX2,WY1,WY2,N,GX,GY,
     .              LIST,LPTR,LEND,TITLE,NUMBR, IERR)
      ENDIF
      IPLT = IPLT/2
      IF (IPLT .NE. 2*(IPLT/2)) THEN
        OPEN (LU3,FILE='rplot.eps')
        WX1 = DXL(1)
        WX2 = WX1
        WY1 = DYL(1)
        WY2 = WY1
        DO 5 I = 2,NV
          IF (DXL(I) .LT. WX1) WX1 = DXL(I)
          IF (DXL(I) .GT. WX2) WX2 = DXL(I)
          IF (DYL(I) .LT. WY1) WY1 = DYL(I)
          IF (DYL(I) .GT. WY2) WY2 = DYL(I)
    5     CONTINUE
        DX = WX2-WX1
        DY = WY2-WY1
        WX1 = WX1 - 0.1D0*DX
        WX2 = WX2 + 0.1D0*DX
        WY1 = WY1 - 0.1D0*DY
        WY2 = WY2 + 0.1D0*DY
        TITLE = '()'
        CALL PLTGR (LU3,PLTSIZ,WX1,WX2,WY1,WY2,N,LIST,LPTR,
     .              LEND,NV,LISTV,DXL,DYL,X,Y,TITLE,
     .              NUMBR, IERR)
      ENDIF
C
C Compute the distance from each node to its cell boundary,
C   and compute DMIN.
C
      CALL DELTAI (N,X,Y,Z,GX,GY,LIST,LPTR,LEND, D,IERR)
      IF (IERR .NE. 0) THEN
        IER = -10
        RETURN
      ENDIF
      DMIN = 1.D20
      DO 6 I = 1,N
        IF (D(I) .LT. DMIN) DMIN = D(I)
    6   CONTINUE
      IF (DMIN .LE. 0.) THEN
        IER = -11
        RETURN
      ENDIF
C
C Compute weights and quadrature points defining the rule
C   for approximating values of F.
C
      CALL GETQW (NR,NA,DMIN, W,QX,QY,IERR)
      RETURN
      END
      SUBROUTINE ADDNDC (NST,K,X,Y,Z, LIST,LPTR,LEND,
     .                   LNEW, IER)
      INTEGER NST, K, LIST(*), LPTR(*), LEND(K), LNEW, IER
      DOUBLE PRECISION X(K), Y(K), Z(K)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   This subroutine adds node K to a convex (data-dependent)
C triangulation of nodes 1,...,K-1, producing a convex tri-
C angulation, if it exists, of nodes 1,...,K.  Refer to
C Subroutine TRMSHC.
C
C   The algorithm consists of the following steps:  node K
C is located relative to the triangulation (TRFIND), its
C index is added to the data structure (INTADD or BDYADD),
C and a sequence of swaps (SWPTC and SWAP) are applied to
C the edges opposite K so that all edges incident on node K
C and opposite node K are locally optimal (satisfy the con-
C vexity test).
C
C
C On input:
C
C       NST = Index of a node at which TRFIND begins its
C             search.  Search time depends on the proximity
C             of this node to K.  If NST < 1, the search is
C             begun at node K-1.
C
C       K = Nodal index (index for X, Y, Z, and LEND) of the
C           new node to be added.  K .GE. 4.
C
C       X,Y,Z = Arrays of length .GE. K containing Car-
C               tesian coordinates of the nodes and data
C               values.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure associated with
C                             the triangulation of nodes 1
C                             to K-1.  The array lengths are
C                             assumed to be large enough to
C                             add node K.  Refer to Subrou-
C                             tine TRMSHC.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node K as the
C                             last entry unless IER .NE. 0,
C                             in which case the arrays are
C                             not altered (except in the
C                             case of IER = -4).
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if K is outside its valid range
C                      on input.
C             IER = -2 if all nodes (including K) are col-
C                      linear.
C             IER = -3 if node K is contained in a triangle
C                      and data point K lies above the cor-
C                      responding surface triangle.  A
C                      convex triangulation does not exist
C                      in this case.
C             IER = -4 if, following the addition of node K
C                      to the data structure, an error con-
C                      dition was flagged by Function SWPTC
C                      indicating that a convex triangula-
C                      tion does not exist.
C             IER =  L if nodes L and K coincide for some
C                      L < K.
C
C Modules required by ADDNDC:  BDYADD, DSTORE, INSERT,
C                                INTADD, JRAND, LEFT,
C                                LSTPTR, SWAP, SWPTC, TRFIND
C
C Intrinsic function called by ADDNDC:  ABS
C
C***********************************************************
C
      DOUBLE PRECISION DX1, DX2, DX3, DY1, DY2, DY3, DZ1,
     .                 DZ2, DZ3, XK, YK, ZK
      INTEGER LSTPTR
      INTEGER I1, I2, I3, IERR, IN1, IO1, IO2, IST, KK, KM1,
     .        L, LP, LPF, LPO1
      LOGICAL SWPTC
C
C Local parameters:
C
C DX1,..,DZ3 = Components of the vectors from surface
C                point K to the vertices of a surface
C                triangle
C I1,I2,I3 =   Vertex indexes of a triangle containing K
C IERR =       Error flag for calls to SWPTC
C IN1 =        Vertex opposite K:  first neighbor of IO2
C                that precedes IO1.  IN1,IO1,IO2 are in
C                counterclockwise order.
C IO1,IO2 =    Adjacent neighbors of K defining an edge to
C                be tested for a swap
C IST =        Index of node at which TRFIND begins its
C                search
C KK =         Local copy of K
C KM1 =        K-1
C L =          Vertex index (I1, I2, or I3) returned in IER
C                if node K coincides with a vertex
C LP =         LIST index (pointer)
C LPF =        LIST pointer to the first neighbor of K
C LPO1 =       LIST pointer to IO1
C XK,YK =      Cartesian coordinates of node K
C
      KK = K
      IF (KK .LT. 4) GO TO 3
C
C Initialization:
C
      KM1 = KK - 1
      IST = NST
      IF (IST .LT. 1) IST = KM1
      XK = X(KK)
      YK = Y(KK)
      ZK = Z(KK)
C
C Find a triangle (I1,I2,I3) containing K or the rightmost
C   (I1) and leftmost (I2) visible boundary nodes as viewed
C   from node K.
C
      CALL TRFIND (IST,XK,YK,KM1,X,Y,LIST,LPTR,LEND, I1,I2,
     .             I3)
C
C Test for collinear nodes.
C
      IF (I1 .EQ. 0) GO TO 4
      IF (I3 .NE. 0) THEN
C
C Node K is contained in triangle (I1,I2,I3).
C
C   Test for duplicate nodes.
C
        L = I1
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 5
        L = I2
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 5
        L = I3
        IF (XK .EQ. X(L)  .AND.  YK .EQ. Y(L)) GO TO 5
C
C   Test for a nonconvex surface.
C
        DX1 = X(I1) - XK
        DX2 = X(I2) - XK
        DX3 = X(I3) - XK
C
        DY1 = Y(I1) - YK
        DY2 = Y(I2) - YK
        DY3 = Y(I3) - YK
C
        DZ1 = Z(I1) - ZK
        DZ2 = Z(I2) - ZK
        DZ3 = Z(I3) - ZK
C
        IF (DX1*(DY2*DZ3-DY3*DZ2) - DY1*(DX2*DZ3-DX3*DZ2) +
     .      DZ1*(DX2*DY3-DX3*DY2) .LT. 0.) GO TO 6
C
C   Add node K to the data structure.
C
        CALL INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      ELSE
C
C Node K is exterior to the triangulation.
C
        CALL BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
      ENDIF
      IER = 0
C
C Initialize variables for optimization of the
C   triangulation.
C
      LP = LEND(KK)
      LPF = LPTR(LP)
      IO2 = LIST(LPF)
      LPO1 = LPTR(LPF)
      IO1 = ABS(LIST(LPO1))
C
C Begin loop:  find the node opposite K.
C
    1 LP = LSTPTR(LEND(IO1),IO2,LIST,LPTR)
        IF (LIST(LP) .LT. 0) GO TO 2
        LP = LPTR(LP)
        IN1 = ABS(LIST(LP))
C
C Swap test:  if a swap occurs, two new edges are
C             opposite K and must be tested.
C
        IF ( SWPTC(IN1,KK,IO1,IO2,X,Y,Z, IERR) ) THEN
          CALL SWAP (IN1,KK,IO1,IO2, LIST,LPTR,LEND, LPO1)
          IO1 = IN1
          GO TO 1
        ELSE
C
C   Test for a nonconvex triangulation.
C
          IF (IERR .NE. 0) GO TO 7
        ENDIF
C
C No swap occurred.  Test for termination and reset
C   IO2 and IO1.
C
    2   IF (LPO1 .EQ. LPF  .OR.  LIST(LPO1) .LT. 0) RETURN
        IO2 = IO1
        LPO1 = LPTR(LPO1)
        IO1 = ABS(LIST(LPO1))
        GO TO 1
C
C KK < 4.
C
    3 IER = -1
      RETURN
C
C All nodes are collinear.
C
    4 IER = -2
      RETURN
C
C Nodes L and K coincide.
C
    5 IER = L
      RETURN
C
C Nonconvex triangulation.
C
    6 IER = -3
      RETURN
C
C Error flag returned by SWPTC.
C
    7 IER = -4
      RETURN
      END
      SUBROUTINE ADDQT (N,X,Y, Z,GX,GY, EPS,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), Z(N), GX(N), GY(N), EPS
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/27/96
C
C   Given a set of strictly convex Hermite data (nodes p_i,
C data values z_i, and nodal gradients g_i, i = 1 to N, such
C that there exists a convex Hermite interpolant of the data
C values and gradients, and no four data points are co-
C planar), this subroutine adjusts the data values and
C gradients such that a strictly convex Hermite interpolant
C of the data can be constructed from a convex Hermite
C interpolant of the adjusted data.  The adjusted data val-
C ues and gradients are defined as follows:
C
C   Z_i = z_i - eps*<p_i,p_i>  and  G_i = g_i - 2*eps*p_i
C
C for eps < min {(z_j-z_i-<g_i,p_j-p_i>)/<p_j-p_i,p_j-p_i>},
C
C where the minimum is over all distinct values of i and j
C in the range 1 to N.  By convexity of the data, eps > 0,
C and it follows from the definition of eps that Z_j - Z_i -
C <G_i,p_j-p_i> > 0, and the adjusted data is therefore also
C strictly convex.  Letting h denote a convex interpolant of
C the adjusted data, the function F(p) = h(p) + eps*<p,p> is
C a strictly convex interpolant of the original data since
C eps*<p,p> is strictly convex.
C
C   The value of eps is arbitrarily chosen to be 0.95 times
C its upper bound.  Computation of the bound requires O(N*N)
C operations.
C
C
C On input:
C
C       N = Number of nodes, data values, and gradients.
C           N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of a set of distinct nodes.
C
C The above parameters are not altered by this routine.
C
C       Z = Array of length N containing the data values.
C
C       GX,GY = Arrays of length N containing the components
C               of the nodal gradients.
C
C On output:
C
C       Z = Adjusted data values unless IER > 0.
C
C       GX,GY = Adjusted nodal gradients unless IER > 0.
C
C       EPS = Scale factor for the quadratic term unless
C             IER > 0.  EPS = 0 if the data is not strictly
C             convex.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C             IER = 2 if a pair of nodes coincide.
C
C Modules required by ADDQT:  None
C
C Intrinsic function called by ADDQT:  MIN
C
C***********************************************************
C
      DOUBLE PRECISION DS, DX, DY, DZ, S, T1, T2, TMIN
      INTEGER I, J
      DATA S/0.95D0/
C
C Local parameters:
C
C DS =    <p_j-p_i,p_j-p_i>
C DX,DY = Components of p_j-p_i
C DZ =    z_j-z_i
C I,J =   Nodal indexes and DO-loop indexes
C S =     Scale factor for the upper bound on EPS
C T1 =    Term associated with i,j in the bound on EPS
C T2 =    Term associated with j,i in the bound on EPS
C TMIN =  Bound on EPS
C
      IF (N .LT. 3) GO TO 11
      TMIN = 1.D20
C
C Loop on pairs of nodal indexes I,J in the upper triangle.
C
      DO 2 I = 1,N-1
        DO 1 J = I+1,N
          DX = X(J)-X(I)
          DY = Y(J)-Y(I)
          DS = DX*DX + DY*DY
          IF (DS .EQ. 0.) GO TO 12
          DZ = Z(J)-Z(I)
          T1 = ( DZ - GX(I)*DX-GY(I)*DY )/DS
          T2 = -( DZ - GX(J)*DX-GY(J)*DY )/DS
          IF (T1 .LE. 0.  .OR.  T2 .LE. 0.) GO TO 4
          TMIN = MIN(TMIN,T1,T2)
    1     CONTINUE
    2   CONTINUE
      EPS = S*TMIN
C
C Adjust the data values and nodal gradients.
C
      DO 3 I = 1,N
        Z(I) = Z(I) - EPS*(X(I)**2 + Y(I)**2)
        GX(I) = GX(I) - 2.D0*EPS*X(I)
        GY(I) = GY(I) - 2.D0*EPS*Y(I)
    3   CONTINUE
      GO TO 5
C
C The data is not strictly convex.
C
    4 EPS = 0.
C
C No error encountered.
C
    5 IER = 0
      RETURN
C
C N < 3.
C
   11 IER = 1
      RETURN
C
C Nodes I and J coincide.
C
   12 IER = 2
      RETURN
      END
      DOUBLE PRECISION FUNCTION AREAP (X,Y,NB,NODES)
      INTEGER NB, NODES(NB)
      DOUBLE PRECISION X(*), Y(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/21/90
C
C   Given a sequence of NB points in the plane, this func-
C tion computes the signed area bounded by the closed poly-
C gonal curve which passes through the points in the
C specified order.  Each simple closed curve is positively
C oriented (bounds positive area) if and only if the points
C are specified in counterclockwise order.  The last point
C of the curve is taken to be the first point specified, and
C this point should therefore not be specified twice.
C
C   The area of a triangulation may be computed by calling
C AREAP with values of NB and NODES determined by Subroutine
C BNODES.
C
C
C On input:
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of a set of points in the plane
C             for some N .GE. NB.
C
C       NB = Length of NODES.
C
C       NODES = Array of length NB containing the ordered
C               sequence of nodal indexes (in the range
C               1 to N) which define the polygonal curve.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       AREAP = Signed area bounded by the polygonal curve,
C              or zero if NB < 3.
C
C Modules required by AREAP:  None
C
C***********************************************************
C
      DOUBLE PRECISION A
      INTEGER I, ND1, ND2, NNB
C
C Local parameters:
C
C A =       Partial sum of signed (and doubled) trapezoid
C             areas
C I =       DO-loop and NODES index
C ND1,ND2 = Elements of NODES
C NNB =     Local copy of NB
C
      NNB = NB
      A = 0.
      IF (NNB .LT. 3) GO TO 2
      ND2 = NODES(NNB)
C
C Loop on line segments NODES(I-1) -> NODES(I), where
C   NODES(0) = NODES(NB), adding twice the signed trapezoid
C   areas (integrals of the linear interpolants) to A.
C
      DO 1 I = 1,NNB
        ND1 = ND2
        ND2 = NODES(I)
        A = A + (X(ND2)-X(ND1))*(Y(ND1)+Y(ND2))
    1   CONTINUE
C
C A contains twice the negative signed area of the region.
C
    2 AREAP = -A/2.D0
      RETURN
      END
      SUBROUTINE BDYADD (KK,I1,I2, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, LIST(*), LPTR(*), LEND(*), LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/91
C
C   This subroutine adds a boundary node to a triangulation
C of a set of points in the plane.  The data structure is
C updated with the insertion of node KK, but no optimization
C is performed.
C
C
C On input:
C
C       KK = Index of a node to be connected to the sequence
C            of all visible boundary nodes.  KK .GE. 1 and
C            KK must not be equal to I1 or I2.
C
C       I1 = First (rightmost as viewed from KK) boundary
C            node in the triangulation which is visible from
C            node KK (the line segment KK-I1 intersects no
C            edges.
C
C       I2 = Last (leftmost) boundary node which is visible
C            from node KK.  I1 and I2 may be determined by
C            Subroutine TRFIND.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure.
C                             Nodes I1 and I2 must be in-
C                             cluded in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  Node
C                             KK is connected to I1, I2, and
C                             all boundary nodes in between.
C
C Module required by BDYADD:  INSERT
C
C***********************************************************
C
      INTEGER K, LP, LSAV, N1, N2, NEXT, NSAV
C
C Local parameters:
C
C K =     Local copy of KK
C LP =    LIST pointer
C LSAV =  LIST pointer
C N1,N2 = Local copies of I1 and I2, respectively
C NEXT =  Boundary node visible from K
C NSAV =  Boundary node visible from K
C
      K = KK
      N1 = I1
      N2 = I2
C
C Add K as the last neighbor of N1.
C
      LP = LEND(N1)
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = -K
      LPTR(LNEW) = LSAV
      LEND(N1) = LNEW
      LNEW = LNEW + 1
      NEXT = -LIST(LP)
      LIST(LP) = NEXT
      NSAV = NEXT
C
C Loop on the remaining boundary nodes between N1 and N2,
C   adding K as the first neighbor.
C
    1 LP = LEND(NEXT)
        CALL INSERT (K,LP,LIST,LPTR,LNEW)
        IF (NEXT .EQ. N2) GO TO 2
        NEXT = -LIST(LP)
        LIST(LP) = NEXT
        GO TO 1
C
C Add the boundary nodes between N1 and N2 as neighbors
C   of node K.
C
    2 LSAV = LNEW
      LIST(LNEW) = N1
      LPTR(LNEW) = LNEW + 1
      LNEW = LNEW + 1
      NEXT = NSAV
C
    3 IF (NEXT .EQ. N2) GO TO 4
        LIST(LNEW) = NEXT
        LPTR(LNEW) = LNEW + 1
        LNEW = LNEW + 1
        LP = LEND(NEXT)
        NEXT = LIST(LP)
        GO TO 3
C
    4 LIST(LNEW) = -N2
      LPTR(LNEW) = LSAV
      LEND(K) = LNEW
      LNEW = LNEW + 1
      RETURN
      END
      SUBROUTINE BNODES (N,LIST,LPTR,LEND, NODES,NB,NE,NT)
      INTEGER N, LIST(*), LPTR(*), LEND(N), NODES(*), NB,
     .        NE, NT
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   Given a triangulation of N points in the plane, this
C subroutine returns an array containing the indexes, in
C counterclockwise order, of the nodes on the boundary of
C the convex hull of the set of points.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.
C
C The above parameters are not altered by this routine.
C
C       NODES = Integer array of length at least NB
C               (NB .LE. N).
C
C On output:
C
C       NODES = Ordered sequence of boundary node indexes
C               in the range 1 to N.
C
C       NB = Number of boundary nodes.
C
C       NE,NT = Number of edges and triangles, respectively,
C               in the triangulation.
C
C Modules required by BNODES:  None
C
C***********************************************************
C
      INTEGER K, LP, N0, NST
C
C Local parameters:
C
C K =   NODES index
C LP =  LIST pointer
C N0 =  Boundary node to be added to NODES
C NST = First element of nodes (arbitrarily chosen to be
C         the one with smallest index)
C
      NST = 1
    1 LP = LEND(NST)
        IF (LIST(LP) .LT. 0) GO TO 2
        NST = NST + 1
        GO TO 1
C
C Initialization.
C
    2 NODES(1) = NST
      K = 1
      N0 = NST
C
C Traverse the boundary in counterclockwise order.
C
    3 LP = LEND(N0)
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (N0 .EQ. NST) GO TO 4
        K = K + 1
        NODES(K) = N0
        GO TO 3
C
C Termination.
C
    4 NB = K
      NT = 2*N - NB - 2
      NE = NT + N - 1
      RETURN
      END
      SUBROUTINE CIRCUM (X1,Y1,X2,Y2,X3,Y3,RATIO, XC,YC,CR,
     .                   SA,AR)
      LOGICAL RATIO
      DOUBLE PRECISION X1, Y1, X2, Y2, X3, Y3, XC, YC, CR,
     .                 SA, AR
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   12/10/96
C
C   Given three vertices defining a triangle, this subrou-
C tine returns the circumcenter, circumradius, signed
C triangle area, and, optionally, the aspect ratio of the
C triangle.
C
C
C On input:
C
C       X1,...,Y3 = Cartesian coordinates of the vertices.
C
C       RATIO = Logical variable with value TRUE if and only
C               if the aspect ratio is to be computed.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       XC,YC = Cartesian coordinates of the circumcenter
C               (center of the circle defined by the three
C               points) unless SA = 0, in which XC and YC
C               are not altered.
C
C       CR = Circumradius (radius of the circle defined by
C            the three points) unless SA = 0 (infinite
C            radius), in which case CR is not altered.
C
C       SA = Signed triangle area with positive value if
C            and only if the vertices are specified in
C            counterclockwise order:  (X3,Y3) is strictly
C            to the left of the directed line from (X1,Y1)
C            toward (X2,Y2).
C
C       AR = Aspect ratio r/CR, where r is the radius of the
C            inscribed circle, unless RATIO = FALSE, in
C            which case AR is not altered.  AR is in the
C            range 0 to .5, with value 0 iff SA = 0 and
C            value .5 iff the vertices define an equilateral
C            triangle.
C
C Modules required by CIRCUM:  None
C
C Intrinsic functions called by CIRCUM:  ABS, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION DS(3), FX, FY, U(3), V(3)
      INTEGER I
C
C Set U(K) and V(K) to the x and y components, respectively,
C   of the directed edge opposite vertex K.
C
      U(1) = X3 - X2
      U(2) = X1 - X3
      U(3) = X2 - X1
      V(1) = Y3 - Y2
      V(2) = Y1 - Y3
      V(3) = Y2 - Y1
C
C Set SA to the signed triangle area.
C
      SA = (U(1)*V(2) - U(2)*V(1))/2.D0
      IF (SA .EQ. 0.) THEN
        IF (RATIO) AR = 0.
        RETURN
      ENDIF
C
C Set DS(K) to the squared distance from the origin to
C   vertex K.
C
      DS(1) = X1*X1 + Y1*Y1
      DS(2) = X2*X2 + Y2*Y2
      DS(3) = X3*X3 + Y3*Y3
C
C Compute factors of XC and YC.
C
      FX = 0.
      FY = 0.
      DO 1 I = 1,3
        FX = FX - DS(I)*V(I)
        FY = FY + DS(I)*U(I)
    1   CONTINUE
      XC = FX/(4.D0*SA)
      YC = FY/(4.D0*SA)
      CR = SQRT( (XC-X1)**2 + (YC-Y1)**2 )
      IF (.NOT. RATIO) RETURN
C
C Compute the squared edge lengths and aspect ratio.
C
      DO 2 I = 1,3
        DS(I) = U(I)*U(I) + V(I)*V(I)
    2   CONTINUE
      AR = 2.D0*ABS(SA)/
     .     ( (SQRT(DS(1)) + SQRT(DS(2)) + SQRT(DS(3)))*CR )
      RETURN
      END
      SUBROUTINE CNTOUR (NX,NY,X,Y,Z,CVAL,LC,NCMAX,IWK, XC,
     .			 YC,ILC,NC,IER)
      INTEGER NX, NY, LC, NCMAX, IWK(NX,*), ILC(NCMAX), NC,
     .	      IER
      DOUBLE PRECISION X(NX), Y(NY), Z(NX,NY), CVAL, XC(LC),
     .                 YC(LC)
C
C***********************************************************
C
C                                              From CSRFPACK
C					     Robert J. Renka
C				   Dept. of Computer Science
C					Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   04/28/90
C
C   Given a set of function values Z = F(X,Y) at the verti-
C ces of an NX by NY rectangular grid, this subroutine de-
C termines a set of contour lines associated with F = CVAL.
C A contour line is specified by an ordered sequence of
C points (XC,YC), each lying on a grid edge and computed
C from the linear interpolant of the function values at the
C endpoints of the edge.  The accuracy of the contour lines
C is thus directly related to the number of grid points.  If
C a contour line forms a closed curve, the first point coin-
C cides with the last point.  Otherwise, the first and last
C points lie on the grid boundary.
C
C   Note that the problem is ill-conditioned in the vicinity
C of a double zero of F-CVAL.  Thus, if a grid cell is
C crossed by two contour lines (all four sides intersected),
C three different configurations are possible, corresponding
C to a local minimum, a local maximum, or a saddle point.
C It is arbitrarily assumed in this case that the contour
C lines intersect, representing a saddle point.  Also, in
C order to treat the case of F = CVAL at a vertex in a con-
C sistent manner, this case is always treated as F > CVAL.
C Hence, if F takes on the same value at both ends of an
C edge, it is assumed that no contour line intersects that
C edge.  In particular, a constant function, including
C F = CVAL, results in no contour lines.
C
C On input:
C
C	NX = Number of grid points in the x direction.
C	     NX .GE. 2.
C
C	NY = Number of grid points in the y direction.
C	     NY .GE. 2.
C
C	X = Array of length NX containing a strictly in-
C	    creasing sequence of values.
C
C	Y = Array of length NY containing a strictly in-
C	    creasing sequence of values.
C
C	Z = Array of function values at the vertices of the
C	    rectangular grid.  Z(I,J) = F(X(I),Y(J)) for
C	    I = 1,...,NX and J = 1,...,NY.
C
C	CVAL = Constant function value defining a contour
C	       line as the set of points (X,Y) such that
C	       F(X,Y) = CVAL.
C
C	LC = Length of arrays XC and YC, and maximum allow-
C	     able number of points defining contour lines.
C	     LC = 2(NX-1)(NY-1) + (NX*NY+1)/2 is (probably
C	     more than) sufficient.  LC .GE. 2.
C
C	NCMAX = Length of array ILC, and maximum allowable
C		number of contour lines.  NCMAX = (NX*NY+1)/
C		2 is sufficient.  NCMAX .GE. 1.
C
C The above parameters are not altered by this routine.
C
C	IWK = Integer array of length .GE. NX*(NY-1) to be
C	      used as work space.
C
C	XC,YC = Arrays of length LC.
C
C	ILC = Integer array of length NCMAX.
C
C On output:
C
C	XC,YC = Arrays containing the coordinates of NC con-
C		tour lines.  For K = 1,...,NC, contour line
C		K is defined by the sequence of points with
C               indexes ILC(K-1)+1,...,ILC(K) where ILC(0) =
C		0.
C
C       ILC = Array containing the indexes (to XC and YC)
C	      associated with the terminal point of contour
C	      line K in position K for K = 1,...,NC (if NC
C	      .GT. 0).
C
C	NC = Number of contour lines whose points are stored
C	     in XC and YC.
C
C       IER = Error indicator:
C	      IER = 0 if no errors were encountered and all
C		      contour lines were found.
C             IER = 1 if NX, NY, LC, or NCMAX is outside its
C                     valid range.  NC = 0 and XC, YC, and
C                     ILC are not altered in this case.
C	      IER = 2 if X or Y is not strictly increasing.
C		      NC = 0 and XC, YC, and ILC are not
C		      altered in this case.
C             IER = K for K > LC, where K is the required
C		      length of XC and YC, if more storage
C		      space is required to complete the
C		      specification of contour line NC and/
C		      or additional contour lines up to a
C		      total of NCMAX.  NC .GE. 1 and ILC(NC)
C		      = LC in this case.
C	     IER = -1 if more than NCMAX contour lines are
C		      present (more space is required in
C		      ILC).  NC = NCMAX, and LC may or may
C		      not be sufficient for the additional
C		      contour lines in this case.  (This is
C		      not determined.)
C
C   In the unlikely event of an internal failure, a message
C is printed on logical unit LUN (specified in the DATA
C statement below).  IER may be 0 in this case.
C
C Modules required by CNTOUR:  None
C
C***********************************************************
C
      DOUBLE PRECISION CV, W, XF, XN, XP, YF, YN, YP, Z1, Z2
      INTEGER I, I1, I2, IB, IN, IND, ISID, ISIDB, ISIDN,
     .        J, J1, J2, JB, JN, K, LCON, LMX, LUN, NCMX,
     .        NCON, NI, NIM1, NJ, NJM1
      LOGICAL BDRY
      DATA    LUN/0/
C
C Store parameters in local variables.
C
      NI = NX
      NJ = NY
      NIM1 = NI - 1
      NJM1 = NJ - 1
      CV = CVAL
      LMX = LC
      NCMX = NCMAX
      NC = 0
C
C Test for invalid input parameters.
C
      IER = 1
      IF (NI .LT. 2  .OR.  NJ .LT. 2  .OR.  LMX .LT. 2	.OR.
     .	  NCMX .LT. 1) RETURN
C
C Test for nonincreasing values of X or Y.
C
      IER = 2
      DO 1 I = 2,NI
        IF (X(I) .LE. X(I-1)) RETURN
    1   CONTINUE
      DO 2 J = 2,NJ
        IF (Y(J) .LE. Y(J-1)) RETURN
    2   CONTINUE
C
C Loop on grid cells, initializing edge indicators (stored
C   in IWK) to zeros.  For each cell, the indicator IND is a
C   4-bit integer with each bit corresponding to an edge of
C   the cell, and having value 1 iff the edge has been pro-
C   cessed.  Note that two IND values must be adjusted when
C   an interior edge is processed.  The cell sides (edges)
C   are numbered (1,2,4,8) in counterclockwise order start-
C   ing from the bottom.  This corresponds to an ordering of
C   the weighted IND bits from low order to high order.
C   Grid cells are identified with their lower left corners.
C
      DO 4 J = 1,NJM1
	DO 3 I = 1,NIM1
          IWK(I,J) = 0
    3     CONTINUE
    4   CONTINUE
C
C First determine open contours by looping on boundary edges
C   in counterclockwise order starting from the lower left.
C   For each unprocessed boundary edge intersected by a con-
C   tour line, the contour line is determined and IWK is up-
C   dated to reflect the edges intersected.  The boundary
C   cell (lower left corner) is indexed by (IB,JB) and the
C   boundary edge is specified by ISIDB.  NCON and LCON are
C   local variables containing the number of contour lines
C   encountered and the current length of XC and YC.
C
      NCON = 0
      LCON = 0
      ISIDB = 1
      IB = 1
      JB = 1
C
C Top of loop on boundary edges.  The edge has been
C   processed iff IND/ISIDB is odd.
C
    5 IND = IWK(IB,JB)
      IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 9
C
C Update the edge indicator and store the vertex indexes of
C   the endpoints of the edge.
C
      IWK(IB,JB) = IND + ISIDB
      IF (ISIDB .EQ. 1) THEN
	I1 = IB
	J1 = JB
	I2 = IB + 1
	J2 = JB
      ELSEIF (ISIDB .EQ. 2) THEN
	I1 = IB + 1
	J1 = JB
	I2 = IB + 1
	J2 = JB + 1
      ELSEIF (ISIDB .EQ. 4) THEN
	I1 = IB + 1
	J1 = JB + 1
	I2 = IB
	J2 = JB + 1
      ELSE
	I1 = IB
	J1 = JB + 1
	I2 = IB
	J2 = JB
      ENDIF
C
C Proceed to the next edge if there is no intersection.
C
      Z1 = Z(I1,J1)
      Z2 = Z(I2,J2)
      IF ((Z1 .LT. CV  .AND.  Z2 .LT. CV)  .OR.
     .    (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 9
C
C Store the zero of the linear interpolant of Z1-CV and
C   Z2-CV as the first point of an open contour unless
C   NCMAX contour lines have been found or there is in-
C   sufficient space reserved for XC and YC.
C
      IF (NCON .EQ. NCMX) THEN
	IER = -1
        GO TO 16
      ENDIF
      NCON = NCON + 1
      LCON = LCON + 1
      W = (CV-Z1)/(Z2-Z1)
      XP = X(I1) + W*(X(I2)-X(I1))
      YP = Y(J1) + W*(Y(J2)-Y(J1))
      IF (LCON .LE. LMX) THEN
	XC(LCON) = XP
	YC(LCON) = YP
      ENDIF
C
C Initialize for loop on cells intersected by the open
C   contour line.
C
      I = IB
      J = JB
      ISID = ISIDB
C
C Traverse the contour line.  Cell (I,J) was entered on side
C   ISID = (I1,J1)->(I2,J2).  Find an exit edge E (unproces-
C   sed edge intersected by the contour) by looping on the
C   remaining three sides, starting with the side opposite
C   ISID.
C
    6 IND = IWK(I,J)
      DO 7 K = 1,3
	ISID = 2*ISID
	IF (K .NE. 2) ISID = 2*ISID
	IF (ISID .GT. 15) ISID = ISID/16
	IF (ISID .EQ. 1) THEN
	  I1 = I
	  J1 = J
	  I2 = I + 1
	  J2 = J
	ELSEIF (ISID .EQ. 2) THEN
	  I1 = I + 1
	  J1 = J
	  I2 = I + 1
	  J2 = J + 1
	ELSEIF (ISID .EQ. 4) THEN
	  I1 = I + 1
	  J1 = J + 1
	  I2 = I
	  J2 = J + 1
	ELSE
	  I1 = I
	  J1 = J + 1
	  I2 = I
	  J2 = J
	ENDIF
C
C Test for a 1 in bit position ISID of cell (I,J) and bypass
C   the edge if it has been previously encountered.
C
        IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 7
C
C Update IWK for edge E = (I1,J1)->(I2,J2).  (IN,JN) indexes
C   the cell which shares E with cell (I,J), and ISIDN is
C   the side number of E in (IN,JN).  BDRY is true iff E is
C   a boundary edge (with no neighboring cell).
C
	IWK(I,J) = IWK(I,J) + ISID
	IF (ISID .LE. 2) THEN
	  IN = I1
	  JN = J2 - 1
	  ISIDN = 4*ISID
	ELSE
	  IN = I1 - 1
	  JN = J2
	  ISIDN = ISID/4
	ENDIF
	BDRY = IN .EQ. 0  .OR.	IN .EQ. NI  .OR.
     .	       JN .EQ. 0  .OR.	JN .EQ. NJ
	IF (.NOT. BDRY) IWK(IN,JN) = IWK(IN,JN) + ISIDN
C
C Exit the loop on sides if E is intersected by the contour.
C
	Z1 = Z(I1,J1)
	Z2 = Z(I2,J2)
	IF ((Z1 .LT. CV  .AND.	Z2 .GE. CV)  .OR.
     .      (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 8
    7   CONTINUE
C*
C Error -- No exit point found.  Print a message and exit
C	   the contour traversal loop.
C
      WRITE (LUN,100) NCON
  100 FORMAT (///5X,'Error in CNTOUR:  Contour line L ',
     .	      'begins on the boundary'/5X,'and terminates ',
     .        'in the interior for L =',I4/)
      ILC(NCON) = LCON
      GO TO 9
C*
C Add the intersection point (XN,YN) to the list unless it
C   coincides with the previous point (XP,YP) or there is
C   not enough space in XC and YC.
C
    8 W = (CV-Z1)/(Z2-Z1)
      XN = X(I1) + W*(X(I2)-X(I1))
      YN = Y(J1) + W*(Y(J2)-Y(J1))
      IF (XN .NE. XP  .OR.  YN .NE. YP) THEN
	LCON = LCON + 1
	XP = XN
	YP = YN
	IF (LCON .LE. LMX) THEN
	  XC(LCON) = XN
	  YC(LCON) = YN
	ENDIF
      ENDIF
C
C Bottom of contour traversal loop.  If E is not a boundary
C   edge, reverse the edge direction (endpoint indexes) and
C   update the cell index and side number.
C
      IF (.NOT. BDRY) THEN
	I = I1
	J = J1
	I1 = I2
	J1 = J2
	I2 = I
	J2 = J
	I = IN
	J = JN
	ISID = ISIDN
        GO TO 6
      ENDIF
C
C Update ILC with a pointer to the end of the contour line.
C
      ILC(NCON) = LCON
C
C Bottom of loop on boundary edges.  Update the boundary
C   cell index and side number, and test for termination.
C
    9 IF (ISIDB .EQ. 1) THEN
	IF (IB .LT. NIM1) THEN
	  IB = IB + 1
	ELSE
	  ISIDB = 2
	ENDIF
      ELSEIF (ISIDB .EQ. 2) THEN
	IF (JB .LT. NJM1) THEN
	  JB = JB + 1
	ELSE
	  ISIDB = 4
	ENDIF
      ELSEIF (ISIDB .EQ. 4) THEN
	IF (IB .GT. 1) THEN
	  IB = IB - 1
	ELSE
	  ISIDB = 8
	ENDIF
      ELSE
	IF (JB .GT. 1) THEN
	  JB = JB - 1
	ELSE
	  ISIDB = 16
	ENDIF
      ENDIF
      IF (ISIDB .LT. 16) GO TO 5
C
C Determine closed contours by looping on interior edges --
C   the first two sides (bottom and right) of each cell,
C   excluding boundary edges.  The beginning cell is indexed
C   by (IB,JB), and the beginning side number is ISIDB.
C
      DO 15 JB = 1,NJM1
      DO 14 IB = 1,NIM1
      DO 13 ISIDB = 1,2
        IF (JB .EQ. 1  .AND.  ISIDB .EQ. 1) GO TO 13
        IF (IB .EQ. NIM1  .AND.  ISIDB .EQ. 2) GO TO 13
C
C Bypass the edge if it was previously encountered
C   (IND/ISIDB odd).
C
	IND = IWK(IB,JB)
        IF (IND/ISIDB .NE. 2*((IND/ISIDB)/2)) GO TO 13
C
C Determine the endpoint indexes of the beginning edge E =
C   (I1,J1)->(I2,J2), find the index (I,J) and side number
C   ISID of the cell which shares E with (IB,JB), and up-
C   date IWK.
C
	IF (ISIDB .EQ. 1) THEN
	  I1 = IB
	  J1 = JB
	  I2 = IB + 1
	  J2 = JB
	  I = IB
	  J = JB - 1
	  ISID = 4
	ELSE
	  I1 = IB + 1
	  J1 = JB
	  I2 = IB + 1
	  J2 = JB + 1
	  I = I1
	  J = J1
	  ISID = 8
	ENDIF
	IWK(IB,JB) = IND + ISIDB
	IWK(I,J) = IWK(I,J) + ISID
C
C Proceed to the next interior edge if there is no
C   intersection.
C
	Z1 = Z(I1,J1)
	Z2 = Z(I2,J2)
	IF ((Z1 .LT. CV  .AND.	Z2 .LT. CV)  .OR.
     .      (Z1 .GE. CV  .AND.  Z2 .GE. CV)) GO TO 13
C
C Store the intersection point as the first point of a
C   closed contour unless NCMAX contour lines have been
C   found or there is insufficient space in XC and YC.
C
	IF (NCON .EQ. NCMX) THEN
	  IER = -1
          GO TO 16
	ENDIF
	NCON = NCON + 1
	LCON = LCON + 1
	W = (CV-Z1)/(Z2-Z1)
	XP = X(I1) + W*(X(I2)-X(I1))
	YP = Y(J1) + W*(Y(J2)-Y(J1))
	IF (LCON .LE. LMX) THEN
	  XC(LCON) = XP
	  YC(LCON) = YP
	ENDIF
	XF = XP
	YF = YP
C
C Traverse the contour line.  Cell (I,J) was entered on side
C   ISID = edge (I2,J2)->(I1,J1).  Reverse the edge direc-
C   tion.
C
   10   IN = I1
	JN = J1
	I1 = I2
	J1 = J2
	I2 = IN
	J2 = JN
	IND = IWK(I,J)
C
C Find an exit edge E by looping on the remaining three
C   sides, starting with the side opposite ISID.
C
        DO 11 K = 1,3
	  ISID = 2*ISID
	  IF (K .NE. 2) ISID = 2*ISID
	  IF (ISID .GT. 15) ISID = ISID/16
	  IF (ISID .EQ. 1) THEN
	    I1 = I
	    J1 = J
	    I2 = I + 1
	    J2 = J
	  ELSEIF (ISID .EQ. 2) THEN
	    I1 = I + 1
	    J1 = J
	    I2 = I + 1
	    J2 = J + 1
	  ELSEIF (ISID .EQ. 4) THEN
	    I1 = I + 1
	    J1 = J + 1
	    I2 = I
	    J2 = J + 1
	  ELSE
	    I1 = I
	    J1 = J + 1
	    I2 = I
	    J2 = J
	  ENDIF
C
C Bypass the edge if it has been previously encountered.
C
          IF (IND/ISID .NE. 2*((IND/ISID)/2)) GO TO 11
C
C Determine the index (IN,JN) and side number ISIDN of the
C   cell which shares edge E = (I1,J1)->(I2,J2) with cell
C   (I,J), and update IWK.
C
	  IF (ISID .LE. 2) THEN
	    IN = I1
	    JN = J2 - 1
	    ISIDN = 4*ISID
	  ELSE
	    IN = I1 - 1
	    JN = J2
	    ISIDN = ISID/4
	  ENDIF
	  IWK(I,J) = IWK(I,J) + ISID
	  IWK(IN,JN) = IWK(IN,JN) + ISIDN
C
C Exit the loop on sides if E is intersected.
C
	  Z1 = Z(I1,J1)
	  Z2 = Z(I2,J2)
	  IF ((Z1 .LT. CV  .AND.  Z2 .GE. CV)  .OR.
     .        (Z1 .GE. CV  .AND.  Z2 .LT. CV)) GO TO 12
   11     CONTINUE
C*
C Error -- No exit point found.  Print a message and exit
C	   the contour traversal loop.
C
        WRITE (LUN,110) NCON
  110   FORMAT (///5X,'Error in CNTOUR:  Contour line L ',
     .		'is open but'/5X,'does not intersect the ',
     .          'boundary for L =',I4/)
	ILC(NCON) = LCON
        GO TO 13
C*
C Add the intersection point to the list unless it coincides
C   with the previous point or there is not enough space in
C   XC and YC.
C
   12   W = (CV-Z1)/(Z2-Z1)
	XN = X(I1) + W*(X(I2)-X(I1))
	YN = Y(J1) + W*(Y(J2)-Y(J1))
	IF (XN .NE. XP	.OR.  YN .NE. YP) THEN
	  LCON = LCON + 1
	  XP = XN
	  YP = YN
	  IF (LCON .LE. LMX) THEN
	    XC(LCON) = XN
	    YC(LCON) = YN
	  ENDIF
	ENDIF
C
C Bottom of contour traversal loop.  If the next cell is not
C   the beginning cell, update the cell index and side num-
C   ber.
C
	IF (IN .NE. IB	.OR.  JN .NE. JB) THEN
	  I = IN
	  J = JN
	  ISID = ISIDN
          GO TO 10
	ENDIF
C
C Add the first point as the last point (unless the first
C   and last points already coincide), and update ILC.
C
	IF (XP .NE. XF	.OR.  YP .NE. YF) THEN
	  LCON = LCON + 1
	  IF (LCON .LE. LMX) THEN
	    XC(LCON) = XF
	    YC(LCON) = YF
	  ENDIF
	ENDIF
	ILC(NCON) = LCON
C
C Bottom of loop on interior edges.
C
   13   CONTINUE
   14   CONTINUE
   15   CONTINUE
      IER = 0
C
C Test for insufficient storage reserved for XC and YC.
C
   16 IF (LCON .GT. LMX) IER = LCON
      NC = NCON
      RETURN
      END
      SUBROUTINE DELARC (N,IO1,IO2, LIST,LPTR,LEND,
     .                   LNEW, IER)
      INTEGER N, IO1, IO2, LIST(*), LPTR(*), LEND(N), LNEW,
     .        IER
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/12/94
C
C   This subroutine deletes a boundary arc from a triangula-
C tion.  It may be used to remove a null triangle from the
C convex hull boundary.  Note, however, that if the union of
C triangles is rendered nonconvex, Subroutines DELNDC and
C TRFIND may fail.  Thus, Subroutines ADDNDC and DELNDC
C should not be called following an arc deletion.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 4.
C
C       IO1,IO2 = Indexes (in the range 1 to N) of a pair of
C                 adjacent boundary nodes defining the arc
C                 to be removed.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Triangulation data structure
C                             created by Subroutine TRMSHC.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of arc IO1-IO2
C                             unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, IO1, or IO2 is outside its valid
C                     range, or IO1 = IO2.
C             IER = 2 if IO1-IO2 is not a boundary arc.
C             IER = 3 if the node opposite IO1-IO2 is al-
C                     ready a boundary node, and thus IO1
C                     or IO2 has only two neighbors or a
C                     deletion would result in two triangu-
C                     lations sharing a single node.
C             IER = 4 if one of the nodes is a neighbor of
C                     the other, but not vice versa, imply-
C                     ing an invalid triangulation data
C                     structure.
C
C Modules required by DELARC:  DELNB, LSTPTR
C
C Intrinsic function called by DELARC:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPL, N1, N2, N3
C
C Local parameters:
C
C LP =    LIST index (pointer)
C LPH =   Pointer (or flag) returned by DELNB
C LPL =   Pointer to the last neighbor of N1, N2, or N3
C N1,N2 = Local copies of IO1 and IO2 ordered so that N1-N2
C           is a (directed) edge
C N3 =    Node opposite N1->N2
C
      N1 = IO1
      N2 = IO2
C
C Test for errors, and set N1->N2 to the directed boundary
C   edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
C   for some N3.
C
      IF (N .LT. 4  .OR.  N1 .LT. 1  .OR.  N1 .GT. N  .OR.
     .    N2 .LT. 1  .OR.  N2 .GT. N  .OR.  N1 .EQ. N2) THEN
        IER = 1
        RETURN
      ENDIF
C
      LPL = LEND(N2)
      IF (-LIST(LPL) .NE. N1) THEN
        N1 = N2
        N2 = IO1
        LPL = LEND(N2)
        IF (-LIST(LPL) .NE. N1) THEN
          IER = 2
          RETURN
        ENDIF
      ENDIF
C
C Set N3 to the node opposite N1->N2 (the second neighbor
C   of N1), and test for error 3 (N3 already a boundary
C   node).
C
      LPL = LEND(N1)
      LP = LPTR(LPL)
      LP = LPTR(LP)
      N3 = ABS(LIST(LP))
      LPL = LEND(N3)
      IF (LIST(LPL) .LE. 0) THEN
        IER = 3
        RETURN
      ENDIF
C
C Delete N2 as a neighbor of N1, making N3 the first
C   neighbor, and test for error 4 (N2 not a neighbor
C   of N1).  Note that previously computed pointers may
C   no longer be valid following the call to DELNB.
C
      CALL DELNB (N1,N2,N, LIST,LPTR,LEND,LNEW, LPH)
      IF (LPH .LT. 0) THEN
        IER = 4
        RETURN
      ENDIF
C
C Delete N1 as a neighbor of N2, making N3 the new last
C   neighbor.
C
      CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNEW, LPH)
C
C Make N3 a boundary node with first neighbor N2 and last
C   neighbor N1.
C
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      LEND(N3) = LP
      LIST(LP) = -N1
C
C No errors encountered.
C
      IER = 0
      RETURN
      END
      SUBROUTINE DELBE (N,X,Y,TOL, LIST,LPTR,LEND,LNEW, ND)
      INTEGER N, LIST(*), LPTR(*), LEND(N), LNEW, ND
      DOUBLE PRECISION X(N), Y(N), TOL
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/27/96
C
C   This subroutine deletes extraneous boundary edges from a
C triangulation.  For each triangle (N1,N2,N3) such that
C N1->N2 is a boundary edge and N3 is an interior node, if
C the aspect ratio (ratio of the radius of the inscribed
C circle to that of the circumcircle) falls below tolerance
C TOL, edge N1->N2 (and triangle (N1,N2,N3)) is removed.
C (Following the removal, N1->N3 and N3->N2 become candi-
C dates for removal.)
C   Triangles with aspect ratio zero (signed area zero) lead
C to failure in Subroutine LGRAD (and hence VLIST), and tri-
C angles with small aspect ratio lead to inaccuracy in
C triangle gradients (vertices of the gradient feasibility
C diagram computed by VLIST).  On the other hand, if the
C triangulation is rendered nonconvex by removing a boundary
C edge, Subroutines DELNDC and TRFIND (and hence ADDNDC) may
C fail.  Thus, if it is necessary to add or delete a node
C following a call to this routine, the tolerance should be
C very small.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  No comp-
C           utation takes place if N < 4 or TOL <= 0.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       TOL = Positive tolerance used as a lower bound on
C             the acceptable aspect ratio of boundary tri-
C             angles.  A reasonable value is 1.E-5.
c
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with ND
C                             boundary edge removals.
C
C       ND = Number of boundary edges removed, or -1 if an
C            error flag was returned by DELARC indicating
C            an invalid triangulation data structure.
C
C Modules required by DELBE:  CIRCUM, DELARC, DELNB, LSTPTR
C
C***********************************************************
C
      DOUBLE PRECISION AR, CR, SA, XC, YC
      INTEGER IER, LP, LPL, N0, N1, N2, N3
C
C Local parameters:
C
C AR =    Aspect ratio of triangle (N1,N2,N3)
C CR =    Circumradius of triangle (N1,N2,N3)
C IER =   Error flag for calls to DELARC
C LP =    LIST index of (pointer to) the first or second
C           neighbor of N1
C LPL =   Pointer to the last neighbor of N0 or N1
C N0 =    Nodal index of the first boundary node (the one
C           with smallest index)
C N1,N2 = Indexes of the endpoints of a boundary edge
C N3 =    Index of the node opposite N1->N2
C SA =    Signed area of triangle (N1,N2,N3)
C XC,YC = Components of the circumcenter of triangle
C           (N1,N2,N3)
C
      ND = 0
C
C Test for invalid input.
C
      IF (N .LT. 4  .OR.  TOL .LE. 0.) RETURN
C
C Set N0 and N2 to the first boundary node encountered.
C
      N0 = 0
    1 N0 = N0 + 1
        LPL = LEND(N0)
        IF (LIST(LPL) .GT. 0) GO TO 1
      N2 = N0
C
C CCW loop on boundary edges N1->N2.
C
    2 N1 = N2
        LPL = LEND(N1)
        LP = LPTR(LPL)
        N2 = LIST(LP)
C
C Set N3 to the node opposite N1->N2 -- the second neighbor
C   of N1.
C
        LP = LPTR(LP)
        N3 = LIST(LP)
C
C Bypass triangles for which N3 is a boundary node.
C
        IF (N3 .LT. 0) GO TO 3
        IF (LIST(LEND(N3)) .LT. 0) GO TO 3
C
C Compute the aspect ratio AR of triangle (N1,N2,N3).
C
        CALL CIRCUM (X(N1),Y(N1),X(N2),Y(N2),X(N3),Y(N3),
     .               .TRUE., XC,YC,CR,SA,AR)
        IF (AR .LT. TOL) THEN
C
C Remove edge N1-N2 and test for an invalid data structure.
C
          CALL DELARC (N,N1,N2, LIST,LPTR,LEND,LNEW, IER)
          IF (IER .NE. 0) THEN
            ND = -1
            RETURN
          ENDIF
C
C Update ND and process edge N1->N3 next.
C
          ND = ND + 1
          N2 = N1
          GO TO 2
        ENDIF
C
C Bottom of loop on boundary edges.
C
    3   IF (N2 .NE. N0) GO TO 2
C
C No error encountered.
C
      RETURN
      END
      SUBROUTINE DELNB (N0,NB,N, LIST,LPTR,LEND,LNEW, LPH)
      INTEGER N0, NB, N, LIST(*), LPTR(*), LEND(N), LNEW,
     .        LPH
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/30/98
C
C   This subroutine deletes a neighbor NB from the adjacency
C list of node N0 (but N0 is not deleted from the adjacency
C list of NB) and, if NB is a boundary node, makes N0 a
C boundary node.  For pointer (LIST index) LPH to NB as a
C neighbor of N0, the empty LIST,LPTR location LPH is filled
C in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
C possibly in LEND) is changed to LPH, and LNEW is decremen-
C ted.  This requires a search of LEND and LPTR entailing an
C expected operation count of O(N).
C
C
C On input:
C
C       N0,NB = Indexes, in the range 1 to N, of a pair of
C               nodes such that NB is a neighbor of N0.
C               (N0 need not be a neighbor of NB.)
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the removal of NB from the ad-
C                             jacency list of N0 unless
C                             LPH < 0.
C
C       LPH = List pointer to the hole (NB as a neighbor of
C             N0) filled in by the values at LNEW-1 or error
C             indicator:
C             LPH > 0 if no errors were encountered.
C             LPH = -1 if N0, NB, or N is outside its valid
C                      range.
C             LPH = -2 if NB is not a neighbor of N0.
C
C Modules required by DELNB:  None
C
C Intrinsic function called by DELNB:  ABS
C
C***********************************************************
C
      INTEGER I, LNW, LP, LPB, LPL, LPP, NN
C
C Local parameters:
C
C I =   DO-loop index
C LNW = LNEW-1 (output value of LNEW)
C LP =  LIST pointer of the last neighbor of NB
C LPB = Pointer to NB as a neighbor of N0
C LPL = Pointer to the last neighbor of N0
C LPP = Pointer to the neighbor of N0 that precedes NB
C NN =  Local copy of N
C
      NN = N
C
C Test for error 1.
C
      IF (N0 .LT. 1  .OR.  N0 .GT. NN  .OR.  NB .LT. 1  .OR.
     .    NB .GT. NN  .OR.  NN .LT. 3) THEN
        LPH = -1
        RETURN
      ENDIF
C
C   Find pointers to neighbors of N0:
C
C     LPL points to the last neighbor,
C     LPP points to the neighbor NP preceding NB, and
C     LPB points to NB.
C
      LPL = LEND(N0)
      LPP = LPL
      LPB = LPTR(LPP)
    1 IF (LIST(LPB) .EQ. NB) GO TO 2
        LPP = LPB
        LPB = LPTR(LPP)
        IF (LPB .NE. LPL) GO TO 1
C
C   Test for error 2 (NB not found).
C
      IF (ABS(LIST(LPB)) .NE. NB) THEN
        LPH = -2
        RETURN
      ENDIF
C
C   NB is the last neighbor of N0.  Make NP the new last
C     neighbor and, if NB is a boundary node, then make N0
C     a boundary node.
C
      LEND(N0) = LPP
      LP = LEND(NB)
      IF (LIST(LP) .LT. 0) LIST(LPP) = -LIST(LPP)
      GO TO 3
C
C   NB is not the last neighbor of N0.  If NB is a boundary
C     node and N0 is not, then make N0 a boundary node with
C     last neighbor NP.
C
    2 LP = LEND(NB)
      IF (LIST(LP) .LT. 0  .AND.  LIST(LPL) .GT. 0) THEN
        LEND(N0) = LPP
        LIST(LPP) = -LIST(LPP)
      ENDIF
C
C   Update LPTR so that the neighbor following NB now fol-
C     lows NP, and fill in the hole at location LPB.
C
    3 LPTR(LPP) = LPTR(LPB)
      LNW = LNEW-1
      LIST(LPB) = LIST(LNW)
      LPTR(LPB) = LPTR(LNW)
      DO 4 I = NN,1,-1
        IF (LEND(I) .EQ. LNW) THEN
          LEND(I) = LPB
          GO TO 5
        ENDIF
    4   CONTINUE
C
    5 DO 6 I = 1,LNW-1
        IF (LPTR(I) .EQ. LNW) THEN
          LPTR(I) = LPB
        ENDIF
    6   CONTINUE
C
C No errors encountered.
C
      LNEW = LNW
      LPH = LPB
      RETURN
      END
      SUBROUTINE DELNDC (K, N,X,Y,Z,LIST,LPTR,LEND,LNEW,LWK,
     .                   IWK, IER)
      INTEGER K, N, LIST(*), LPTR(*), LEND(*), LNEW, LWK,
     .        IWK(2,*), IER
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/29/98
C
C   This subroutine deletes node K (along with all edges
C incident on node K) from a convex triangulation of N
C nodes, and inserts edges as necessary to produce a tri-
C angulation of the remaining N-1 nodes.  If a convex
C triangulation is input, a convex triangulation will re-
C sult, and thus, DELNDC reverses the effect of a call to
C Subroutine ADDNDC.
C
C
C On input:
C
C       K = Index (for X, Y, and Z) of the node to be
C           deleted.  1 .LE. K .LE. N.
C
C K is not altered by this routine.
C
C       N = Number of nodes in the triangulation on input.
C           N .GE. 4.  Note that N will be decremented
C           following the deletion.
C
C       X,Y,Z = Arrays of length N containing the coordi-
C               nates of the nodes and data values.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Refer to Sub-
C                             routine TRMSHC.
C
C       LWK = Number of columns reserved for IWK.  LWK must
C             be at least NNB-3, where NNB is the number of
C             neighbors of node K, including an extra
C             pseudo-node if K is a boundary node.
C
C       IWK = Integer work array dimensioned 2 by LWK (or
C             array of length .GE. 2*LWK).
C
C On output:
C
C       N = Number of nodes in the triangulation on output.
C           The input value is decremented unless 1 .LE. IER
C           .LE. 4.
C
C       X,Y,Z = Updated arrays (with elements K+1,...,N+1
C               shifted up one position, thus overwriting
C               element K) unless 1 .LE. IER .LE. 4.
C
C       LIST,LPTR,LEND,LNEW = Updated triangulation data
C                             structure reflecting the dele-
C                             tion unless IER > 0.  Note
C                             that the data structure may
C                             have been altered if IER
C                             .GE. 3.
C
C       LWK = Number of IWK columns required unless IER = 1
C             or IER = 3.
C
C       IWK = Indexes of the endpoints of the new edges
C             added unless LWK = 0 or 1 .LE. IER .LE. 4.
C             (Edges are associated with columns, or pairs
C             of adjacent elements if IWK is declared as a
C             singly-subscripted array.)
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if K or N is outside its valid range
C                     or LWK < 0 on input.
C             IER = 2 if more space is required in IWK.
C                     Refer to LWK.
C             IER = 3 if the triangulation data structure is
C                     invalid on input.
C             IER = 4 if K is an interior node with four or
C                     more neighbors, and the number of
C                     neighbors could not be reduced to
C                     three by swaps.  This could be caused
C                     by floating point errors with collin-
C                     ear nodes or by an invalid data
C                     structure.
C             IER = 5 if an error flag was returned by
C                     OPTIM.  An error message is written to
C                     the standard output unit in this event.
C
C   Note that the deletion may result in all remaining nodes
C being collinear.  This situation is not flagged.
C
C Modules required by DELNDC:  DELNB, LEFT, LSTPTR, NBCNT,
C                                OPTIM, SWAP, SWPTC
C
C Intrinsic function called by DELNDC:  ABS
C
C***********************************************************
C
      DOUBLE PRECISION X1, X2, XL, XR, Y1, Y2, YL, YR
      INTEGER LSTPTR, NBCNT
      LOGICAL LEFT
      INTEGER I, IERR, IWL, J, LNW, LP, LP21, LPF, LPH, LPL,
     .        LPL2, LPN, LWKL, N1, N2, NFRST, NIT, NL, NN,
     .        NNB, NR
      LOGICAL BDRY
C
C Local parameters:
C
C BDRY =    Logical variable with value TRUE iff N1 is a
C             boundary node
C I,J =     DO-loop indexes
C IERR =    Error flag returned by OPTIM
C IWL =     Number of IWK columns containing edges
C LNW =     Local copy of LNEW
C LP =      LIST pointer
C LP21 =    LIST pointer returned by SWAP
C LPF,LPL = Pointers to the first and last neighbors of N1
C LPH =     Pointer (or flag) returned by DELNB
C LPL2 =    Pointer to the last neighbor of N2
C LPN =     Pointer to a neighbor of N1
C LWKL =    Input value of LWK
C N1 =      Local copy of K
C N2 =      Neighbor of N1
C NFRST =   First neighbor of N1:  LIST(LPF)
C NIT =     Number of iterations in OPTIM
C NR,NL =   Neighbors of N1 preceding (to the right of) and
C             following (to the left of) N2, respectively
C NN =      Number of nodes in the triangulation
C NNB =     Number of neighbors of N1 (including a pseudo-
C             node representing the boundary if N1 is a
C             boundary node)
C X1,Y1 =   Coordinates of N1
C X2,Y2 =   Coordinates of N2
C XL,YL =   Coordinates of NL
C XR,YR =   Coordinates of NR
C
C
C Set N1 to K and NNB to the number of neighbors of N1 (plus
C   one if N1 is a boundary node), and test for errors.  LPF
C   and LPL are LIST indexes of the first and last neighbors
C   of N1, IWL is the number of IWK columns containing
C   edges, and BDRY is TRUE iff N1 is a boundary node.
C
      N1 = K
      NN = N
      IF (N1 .LT. 1  .OR.  N1 .GT. NN  .OR.  NN .LT. 4  .OR.
     .    LWK .LT. 0) GO TO 21
      LPL = LEND(N1)
      LPF = LPTR(LPL)
      NNB = NBCNT(LPL,LPTR)
      BDRY = LIST(LPL) .LT. 0
      IF (BDRY) NNB = NNB + 1
      IF (NNB .LT. 3) GO TO 23
      LWKL = LWK
      LWK = NNB - 3
      IF (LWKL .LT. LWK) GO TO 22
      IWL = 0
      IF (NNB .EQ. 3) GO TO 5
C
C Initialize for loop on edges N1-N2 for neighbors N2 of N1,
C   beginning with the second neighbor.  NR and NL are the
C   neighbors preceding and following N2, respectively, and
C   LP indexes NL.  The loop is exited when all possible
C   swaps have been applied to edges incident on N1.  If N1
C   is interior, the number of neighbors will be reduced
C   to 3.
C
      X1 = X(N1)
      Y1 = Y(N1)
      NFRST = LIST(LPF)
      NR = NFRST
      XR = X(NR)
      YR = Y(NR)
      LP = LPTR(LPF)
      N2 = LIST(LP)
      X2 = X(N2)
      Y2 = Y(N2)
      LP = LPTR(LP)
C
C Top of loop:  set NL to the neighbor following N2.
C
    2 NL = ABS(LIST(LP))
      IF (NL .EQ. NFRST  .AND.  BDRY) GO TO 5
      XL = X(NL)
      YL = Y(NL)
C
C   Test for a convex quadrilateral.  To avoid an incorrect
C     test caused by collinearity, use the fact that if N1
C     is a boundary node, then N1 LEFT NR->NL and if N2 is
C     a boundary node, then N2 LEFT NL->NR.
C
      LPL2 = LEND(N2)
      IF ( (BDRY  .OR.  LEFT(XR,YR,XL,YL,X1,Y1))  .AND.
     .     (LIST(LPL2) .LT. 0  .OR.
     .      LEFT(XL,YL,XR,YR,X2,Y2)) ) GO TO 3
C
C   Nonconvex quadrilateral -- no swap is possible.
C
      NR = N2
      XR = X2
      YR = Y2
      GO TO 4
C
C   The quadrilateral defined by adjacent triangles
C     (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
C     NL-NR and store it in IWK.  Indexes larger than N1
C     must be decremented since N1 will be deleted from
C     X and Y.
C
    3 CALL SWAP (NL,NR,N1,N2, LIST,LPTR,LEND, LP21)
      IWL = IWL + 1
      IF (NL .LE. N1) THEN
        IWK(1,IWL) = NL
      ELSE
        IWK(1,IWL) = NL - 1
      ENDIF
      IF (NR .LE. N1) THEN
        IWK(2,IWL) = NR
      ELSE
        IWK(2,IWL) = NR - 1
      ENDIF
C
C   Recompute the LIST indexes LPL,LP and decrement NNB.
C
      LPL = LEND(N1)
      NNB = NNB - 1
      IF (NNB .EQ. 3) GO TO 5
      LP = LSTPTR(LPL,NL,LIST,LPTR)
      IF (NR .EQ. NFRST) GO TO 4
C
C   NR is not the first neighbor of N1.
C     Back up and test N1-NR for a swap again:  Set N2 to
C     NR and NR to the previous neighbor of N1 -- the
C     neighbor of NR which follows N1.  LP21 points to NL
C     as a neighbor of NR.
C
      N2 = NR
      X2 = XR
      Y2 = YR
      LP21 = LPTR(LP21)
      LP21 = LPTR(LP21)
      NR = ABS(LIST(LP21))
      XR = X(NR)
      YR = Y(NR)
      GO TO 2
C
C   Bottom of loop -- test for invalid termination.
C
    4 IF (N2 .EQ. NFRST) GO TO 24
      N2 = NL
      X2 = XL
      Y2 = YL
      LP = LPTR(LP)
      GO TO 2
C
C Delete N1 from the adjacency list of N2 for all neighbors
C   N2 of N1.  LPL points to the last neighbor of N1.
C   LNEW is stored in local variable LNW.
C
    5 LP = LPL
      LNW = LNEW
C
C Loop on neighbors N2 of N1, beginning with the first.
C
    6 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        CALL DELNB (N2,N1,N, LIST,LPTR,LEND,LNW, LPH)
        IF (LPH .LT. 0) GO TO 23
C
C   LP and LPL may require alteration.
C
        IF (LPL .EQ. LNW) LPL = LPH
        IF (LP .EQ. LNW) LP = LPH
        IF (LP .NE. LPL) GO TO 6
C
C Delete N1 from X, Y, Z, and LEND, and remove its adjacency
C   list from LIST and LPTR.  LIST entries (nodal indexes)
C   which are larger than N1 must be decremented.
C
      NN = NN - 1
      IF (N1 .GT. NN) GO TO 9
      DO 7 I = N1,NN
        X(I) = X(I+1)
        Y(I) = Y(I+1)
        Z(I) = Z(I+1)
        LEND(I) = LEND(I+1)
    7   CONTINUE
C
      DO 8 I = 1,LNW-1
        IF (LIST(I) .GT. N1) LIST(I) = LIST(I) - 1
        IF (LIST(I) .LT. -N1) LIST(I) = LIST(I) + 1
    8   CONTINUE
C
C   For LPN = first to last neighbors of N1, delete the
C     preceding neighbor (indexed by LP).
C
C   Each empty LIST,LPTR location LP is filled in with the
C     values at LNW-1, and LNW is decremented.  All pointers
C     (including those in LPTR and LEND) with value LNW-1
C     must be changed to LP.
C
C  LPL points to the last neighbor of N1.
C
    9 IF (BDRY) NNB = NNB - 1
      LPN = LPL
      DO 13 J = 1,NNB
        LNW = LNW - 1
        LP = LPN
        LPN = LPTR(LP)
        LIST(LP) = LIST(LNW)
        LPTR(LP) = LPTR(LNW)
        IF (LPTR(LPN) .EQ. LNW) LPTR(LPN) = LP
        IF (LPN .EQ. LNW) LPN = LP
        DO 10 I = NN,1,-1
          IF (LEND(I) .EQ. LNW) THEN
            LEND(I) = LP
            GO TO 11
          ENDIF
   10     CONTINUE
C
   11   DO 12 I = LNW-1,1,-1
          IF (LPTR(I) .EQ. LNW) LPTR(I) = LP
   12     CONTINUE
   13   CONTINUE
C
C Update N and LNEW, and optimize the patch of triangles
C   containing K (on input) by applying swaps to the edges
C   in IWK.
C
      N = NN
      LNEW = LNW
      IF (IWL .GT. 0) THEN
        NIT = 4*IWL
        CALL OPTIM (X,Y,Z,IWL, LIST,LPTR,LEND,NIT,IWK, IERR)
        IF (IERR .NE. 0) GO TO 25
      ENDIF
C
C Successful termination.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   21 IER = 1
      RETURN
C
C Insufficient space reserved for IWK.
C
   22 IER = 2
      RETURN
C
C Invalid triangulation data structure.  NNB < 3 on input or
C   N2 is a neighbor of N1 but N1 is not a neighbor of N2.
C
   23 IER = 3
      RETURN
C
C K is an interior node with 4 or more neighbors, but the
C   number of neighbors could not be reduced.
C
   24 IER = 4
      RETURN
C
C Error flag returned by OPTIM.
C
   25 IER = 5
      WRITE (*,100) NIT, IERR
      RETURN
  100 FORMAT (//5X,'*** DELNDC:  Error flag returned by ',
     .        'OPTIM:  NIT = ',I4,', IER = ',I1,' ***'/)
      END
      SUBROUTINE DELTAI (N,X,Y,Z,GX,GY,LIST,LPTR,LEND, D,
     .                   IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), IER
      DOUBLE PRECISION X(N), Y(N), Z(N), GX(N), GY(N), D(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/27/96
C
C   Given a set of nodes p_i, data values z_i, and nodal
C gradients g_i, i = 1 to N, for which there exists a convex
C Hermite interpolatory surface, the affine nodal function
C f_i(p) = <g_i,p-p_i> + z_i defines a tangent plane at
C (p_i,z_i), and H(p) = max (over i = 1 to N) {f_i(p)} is a
C convex piecewise linear interpolant of the Hermite data.
C Let R_i denote the set of points p such that f_i(p) >=
C f_j(p) for j = 1 to N.  Then f_i is the restriction of H
C to R_i.  For strictly convex data, R_i is a convex region
C containing a neighborhood of p_i.  The set of such regions
C partitions the plane into a cell diagram.
C
C   Remarkably, it can be shown that the cell diagram is the
C straight-line dual of a convexity-preserving triangulation
C of the nodal gradients and negative z-intercepts (constant
C terms) of the affine nodal functions.
C
C   This subroutine returns a set of radii d_i such that
C H(p) = f_i(p) for all p in a disk of radius d_i centered
C at p_i; i.e., d_i is the perpendicular distance from p_i
C to the boundary of R_i.  The algorithm has complexity
C O(N).
C
C
C On input:
C
C       N = Number of nodes, data values, and gradients.
C           N .GE. 3.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       Z = Array of length N containing the data values.
C
C       GX,GY = Arrays of length N containing the components
C               of the nodal gradients.  These must be dis-
C               tinct.
C
C       LIST,LPTR,LEND = Data structure defining the gradi-
C                        ent triangulation.  Refer to
C                        Subroutine TRMSHC.
C
C The above parameters are not altered by this routine.
C
C       D = Array of length at least N.
C
C On output:
C
C       D = Array containing the distance from each node to
C           its cell boundary unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C             IER = 2 if a pair of nodal gradients coincide.
C             IER = 3 if the Hermite data is not convex.
C
C Modules required by DELTAI:  None
C
C Intrinsic functions called by DELTAI:  ABS, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION DGS, DGX, DGY, DI, DIJ, TOL
      INTEGER I, J, LP, LPL
C
C Local parameters:
C
C DGS =      Euclidean norm squared of g_j - g_i
C DGX, DGY = Components of g_j - g_i
C DI =       Local copy of D(I) -- minimum over J of the
C              DIJ values
C DIJ =      Perpendicular distance from node p_i to cell
C              edge e_ij
C I,J =      Indexes of a pair of adjacent nodal gradients
C              corresponding to cell edge e_ij
C LP =       LIST index of J as a neighbor of I
C LPL =      Pointer to the last neighbor of vertex I
C TOL =      Tolerance for distinguishing between a zero
C              distance and a negative value (error 3)
C
      DATA TOL/1.D-6/
C
C Test for error 1.
C
      IF (N .LT. 3) GO TO 11
C
C Loop on nodal indexes I.
C
      DO 2 I = 1,N
        DI = 1.D20
C
C Loop on indexes J of neighbors of I in the gradient trian-
C   gulation.  The triangulation edge corresponds to the
C   cell edge e_ij shared by cells I and J.  LPL points to
C   the last neighbor of I.
C
        LPL = LEND(I)
        LP = LPL
    1   LP = LPTR(LP)
          J = ABS(LIST(LP))
C
C Compute the perpendicular distance DIJ from p_i to edge
C   e_ij:
C           DIJ = |L(p_i)|/|grad(L)| =
C                 (z_i - f_j(p_i))/|g_j - g_i| =
C                 (z_i - z_j - <g_j,p_i-p_j>)/|g_j - g_i|,
C
C where L(p) = f_i(p) - f_j(p) is a linear function with
C zeros on e_ij.
C
          DGX = GX(J)-GX(I)
          DGY = GY(J)-GY(I)
          DGS = DGX*DGX + DGY*DGY
          IF (DGS .EQ. 0.) GO TO 12
          DIJ = ( Z(I)-Z(J) - GX(J)*(X(I)-X(J))
     .                      - GY(J)*(Y(I)-Y(J)) )/SQRT(DGS)
          IF (DIJ .LT. 0.) THEN
C
C Treat DIJ as 0 if it is in [-TOL,0].
C
            IF (DIJ .LT. -TOL) GO TO 13
            DIJ = 0.
          ENDIF
          IF (DIJ .LT. DI) DI = DIJ
          IF (LP .NE. LPL) GO TO 1
C
C Store DI.
C
        D(I) = DI
    2   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C N < 3.
C
   11 IER = 1
      RETURN
C
C Nodal gradients I and J coincide.
C
   12 IER = 2
      RETURN
C
C f_j(p_i) > z_i, and the Hermite data is therefore not
C   convex.
C
   13 IER = 3
      RETURN
      END
      DOUBLE PRECISION FUNCTION DSTORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/25/96
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C On input:
C
C       X = Double precision value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       DSTORE = Value of X after it has been stored and
C                possibly truncated or rounded to the double
C                precision word length.
C
C Modules required by DSTORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      DSTORE = Y
      RETURN
      END
      SUBROUTINE FGRID (NX,NY,PX,PY,EPS,GX,GY,C,LIST,LPTR,
     .                  LEND,DMIN,NR,NA,W,QX,QY, F,IER)
      INTEGER NX, NY, LIST(*), LPTR(*), LEND(*), NR, NA, IER
      DOUBLE PRECISION PX(NX), PY(NY), EPS, GX(*), GY(*),
     .                 C(*), DMIN, W(NR), QX(NR,NA),
     .                 QY(NR,NA), F(NX,NY)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/27/96
C
C   This subroutine returns the values of the convex bivar-
C iate function F (defined in Subroutine CSURF) at the
C vertices of an NX by NY rectangular grid.
C
C
C On input:
C
C       NX = Number of grid points in the x direction.
C	     NX .GE. 2.
C
C	NY = Number of grid points in the y direction.
C	     NY .GE. 2.
C
C       PX = Array of length NX containing a strictly in-
C            creasing sequence of values.
C
C       PY = Array of length NY containing a strictly in-
C            creasing sequence of values.
C
C       EPS,...,QY = Parameters defining the interpolant F.
C                    These should be input unaltered from a
C                    call to Subroutine CSURF (with IER = 0
C                    or IER = -11 on output).
C
C The above parameters are not altered by this routine.
C
C       F = Array of length at least NX*NY.
C
C On output:
C
C       F = Array of function values at the vertices of the
C           rectangular grid.  F(I,J) = F(X(I),Y(J)) for
C	    I = 1,...,NX and J = 1,...,NY.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if NX < 2, NY < 2, EPS < 0, DMIN < 0,
C                     NR < 1, or NA < 2.
C             IER = 2 if PX or PY is not strictly
C                     increasing.
C
C Module required by FGRID:  FVAL
C
C***********************************************************
C
      DOUBLE PRECISION FX, FY
      INTEGER I, J, KST
C
C Local parameters:
C
C FX,FY = Components of the gradient of F -- not saved
C I,J =   Indexes for PX and PY, respectively
C KST =   Nodal index used as the starting point for the
C           search in Subroutine FVAL
C
      KST = 1
C
C Test for invalid input parameters.
C
      IF (NX .LT. 2  .OR.  NY .LT. 2  .OR.  EPS .LT. 0.
     .    .OR.  DMIN .LT. 0.  .OR.  NR .LT. 1  .OR.
     .    NA .LT. 2) GO TO 11
C
C Test for nonincreasing values of PX or PY.
C
      DO 1 I = 2,NX
        IF (PX(I) .LE. PX(I-1)) GO TO 12
    1   CONTINUE
      DO 2 J = 2,NY
        IF (PY(J) .LE. PY(J-1)) GO TO 12
    2   CONTINUE
C
C Loop on grid points.
C
      DO 4 J = 1,NY
        DO 3 I = 1,NX
          CALL FVAL (PX(I),PY(J),EPS,GX,GY,C,LIST,LPTR,LEND,
     .               DMIN,NR,NA,W,QX,QY, KST, F(I,J),FX,FY)
    3     CONTINUE
    4   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C PX or PY is not strictly increasing.
C
   12 IER = 2
      RETURN
      END
      SUBROUTINE FVAL (PX,PY,EPS,GX,GY,C,LIST,LPTR,LEND,
     .                 DMIN,NR,NA,W,QX,QY, KST, F,FX,FY)
      INTEGER LIST(*), LPTR(*), LEND(*), NR, NA, KST
      DOUBLE PRECISION PX, PY, EPS, GX(*), GY(*), C(*),
     .                 DMIN, W(NR), QX(NR,NA), QY(NR,NA), F,
     .                 FX, FY
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/20/96
C
C   This subroutine returns the value and gradient of F at a
C point p, where F is defined in Subroutine CSURF.  Evalua-
C tion of the gradient adds relatively little to the cost of
C evaluating F.
C
C
C On input:
C
C       PX,PY = Components of the evaluation point p.
C
C       EPS,...,QY = Parameters defining the interpolant F.
C                    These should be input unaltered from a
C                    call to Subroutine CSURF (with IER = 0
C                    on output).
C
C The above parameters are not altered by this routine.
C
C       KST = Nodal index in the range 1 to N defining a
C             cell R_k, k = KST, at which a search for p
C             is started.  If F is evaluated at the elements
C             of an ordered sequence of points, the value of
C             KST returned by a previous call to FVAL is
C             a good choice.
C
C   It is assumed without a test that all input parameters
C are within their valid ranges.
C
C On output:
C
C       KST = Index of the cell R_k that contains p.
C
C       F = Value of the interpolant at p.
C
C       FX,FY = Components of the gradient of F at p.
C
C Modules required by FVAL:  None
C
C Intrinsic functions called by FVAL:  ABS, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION FK, FL, FS, FXS, FYS, QXIJ, QYIJ
      INTEGER I, J, K, L, LP, LPL
C
C Local parameters:
C
C FK =        Value of the affine nodal function f_k at p or
C               at p+q_ij for quadrature point q_ij
C FL =        Value of the affine nodal function f_l at p or
C               at p+q_ij for quadrature point q_ij
C FS =        Sum(j=1,NA){H(p+q_ij)} for some i (1 to NR)
C FXS,FYS =   Components of the gradient of FS
C I =         Index for W and row index for QX and QY
C J =         Column index for QX and QY
C K,L =       Nodal (cell) indexes, where L is a neighbor of
C               K in the gradient triangulation
C LP =        Pointer (LIST index) to L as a neighbor of K
C LPL =       Pointer to the last neighbor of K in the
C               gradient triangulation
C QXIJ,QYIJ = Components of p+q_ij
C
      K = KST
C
C Find k such that p is in R_k:  f_k(p) GE f_l(p) for all
C   neighbors l of k, where f_k(p) = <g_k,p> - c_k.
C
      FK = GX(K)*PX + GY(K)*PY - C(K)
    1 LPL = LEND(K)
      LP = LPL
    2 LP = LPTR(LP)
        L = ABS(LIST(LP))
        FL = GX(L)*PX + GY(L)*PY - C(L)
        IF (FK .LT. FL) THEN
          K = L
          FK = FL
          GO TO 1
        ENDIF
        IF (LP .NE. LPL) GO TO 2
C
C Update KST, and initialize F and (FX,FY) to the quadratic
C   function q(p) = EPS*<p,p> and its gradient 2*EPS*p,
C   respectively.  (Refer to Subroutine ADDQT.)  The inter-
C   polated value and gradient will be added.
C
      KST = K
      F = EPS*(PX*PX + PY*PY)
      FX = 2.D0*EPS*PX
      FY = 2.D0*EPS*PY
C
C F(p) = H(p) = f_k(p) if the disk of radius DMIN centered
C   at p is contained in R_k; i.e., the distance from p to
C   each edge e_kl of R_k is at least DMIN; i.e., for each
C   neighbor l of k, L(p)/abs(grad(L)) >= DMIN, where L(p) =
C   f_k(p)-f_l(p) and grad(L) denotes the gradient of L (a
C   linear function with zeros on e_kl).
C
    3 LP = LPTR(LP)
        L = ABS(LIST(LP))
        FL = GX(L)*PX + GY(L)*PY - C(L)
        IF (FK-FL .LT. DMIN*SQRT((GX(K)-GX(L))**2 +
     .                           (GY(K)-GY(L))**2)) GO TO 4
        IF (LP .NE. LPL) GO TO 3
C
C F(p) = H(p).
C
      F = F + FK
      FX = FX + GX(K)
      FY = FY + GY(K)
      RETURN
C
C F(p) must be approximated by the quadrature rule
C
C   Sum(i=1,NR){ W_i*Sum(j=1,NA){H(p+q_ij)} }
C
C   (Refer to Subroutine GETQW.)
C
    4 DO 7 I = 1,NR
C
C Accumulate Sum(j=1,NA){H(p+q_ij)} and its gradient in
C   FS, FXS, and FYS.
C
        FS = 0.
        FXS = 0.
        FYS = 0.
        DO 6 J = 1,NA
C
C Store the translated quadrature point p+q_ij.
C
          QXIJ = PX + QX(I,J)
          QYIJ = PY + QY(I,J)
C
C Find k such that p+q_ij is in R_k.  The search begins with
C   the index of the cell containing the previous quadrature
C   point (or p if I = J = 1).
C
          FK = GX(K)*QXIJ + GY(K)*QYIJ - C(K)
          LP = LPL
    5     LP = LPTR(LP)
            L = ABS(LIST(LP))
            FL = GX(L)*QXIJ + GY(L)*QYIJ - C(L)
            IF (FK .LT. FL) THEN
              K = L
              FK = FL
              LPL = LEND(K)
              LP = LPL
              GO TO 5
            ENDIF
            IF (LP .NE. LPL) GO TO 5
C
C Update FS and (FXS,FYS) with the value and gradient of f_k
C   at (p+q_ij).
C
          FS = FS + FK
          FXS = FXS + GX(K)
          FYS = FYS + GY(K)
    6     CONTINUE
C
C Update F and (FX,FY).
C
        F = F + W(I)*FS
        FX = FX + W(I)*FXS
        FY = FY + W(I)*FYS
    7   CONTINUE
      RETURN
      END
      SUBROUTINE GETQW (NR,NA,D, W,QX,QY,IER)
      INTEGER NR, NA, IER
      DOUBLE PRECISION D, W(NR), QX(NR,NA), QY(NR,NA)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/29/96
C
C   This subroutine computes weights and abscissae defining
C a quadrature rule for approximating
C
C     F(p) = Integral [H(p+q)*Phi(q)]dq,
C
C where the integral is over the disk D0 of radius D center-
C ed at the origin, and Phi(q) = phi(Norm(q)/D), normalized
C to have integral 1, for phi(t) = 1 - 3*t**2 + 2*t**3 (t in
C [0,1]).  The disk D0 is partitioned uniformly by radius
C and angle resulting in Nr*Na subregions, where region
C (i,j) is defined by polar coordinates r in [(i-1)*D/Nr,
C i*D/Nr] and angle a in [(j-1)*(2*pi)/Na,j*(2*pi)/Na].  The
C quadrature points (abscissae) are the centroids of the
C subregions, and the weights are the centroid values of Phi
C scaled by the subregion areas.  The rule is thus
C
C     Sum(i = 1 to Nr) [W_i*Sum(j = 1 to Na) H(p+q_ij) ],
C
C where W_i = (10/3)*[(2*i-1)/(Nr**2*Na)]*phi((2*i-1)/(2*Nr))
C and q_ij = [(2*i-1)/(2*Nr)]*D*(cos(a_j),sin(a_j)) for
C a_j = [(2*j-1)/(2*Na)]*(2*pi).
C
C   If the radius D of D0 depends on the evaluation point p,
C it is still only necessary to call this routine once with
C D = 1, and then to scale the quadrature points appropri-
C ately for each evaluation of F; i.e., the approximation to
C F(p) is
C
C     Sum(i = 1 to Nr) [W_i*Sum(j = 1 to Na) H(p+D*q_ij) ].
C
C
C On input:
C
C       NR = Number of intervals Nr in the radial direction
C            defining the partitioning of D0.  NR = 8 should
C            be sufficient.  NR > 0.
C
C       NA = Number of intervals Na in the angular direction
C            defining the partitioning of D0.  NA = 3*NR is
C            recommended.  NA > 1.
C
C       D = Radius of the disk D0 -- the domain of integra-
C           tion and support of Phi.  D > 0.
C
C The above parameters are not altered by this routine.
C
C       W = Array of length at least NR.
C
C       QX,QY = Arrays of length at least NR*NA.
C
C On output:
C
C       W = Array containing quadrature weights W_i unless
C           IER > 0.
C
C       QX,QY = Arrays dimensioned NR by NA containing the
C               components of the quadrature points:  q_ij =
C               (QX(i,j),QY(i,j)) for i = 1 to NR and j = 1
C               to NA, unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if NR .LE. 0, NA .LE. 1, or D .LE. 0.
C
C Modules required by GETQW:  None
C
C Intrinsic functions called by GETQW:  ATAN, COS, DBLE, SIN
C
C***********************************************************
C
      DOUBLE PRECISION AJ, PHIT, PIDNA, RI, S, T, TWONR
      INTEGER I, J
C
C Local parameters:
C
C AJ =    Angular coordinate of a quadrature point:
C           ((2*j-1)/NA)*pi
C I,J =   Indexes for W, QX, and QY
C PHIT =  phi(T)
C PIDNA = pi/NA
C RI =    Radial coordinate of a quadrature point:  T*D
C S =     Scale factor for weight W(I):  (10/3)/(NR**2*NA)
C T =     Evaluation point for phi:  (2*i-1)/(2*NR)
C TWONR = 2*NR
C
      IF (NR .LE. 0  .OR.  NA .LE. 1  .OR.  D .LE. 0.)
     .  GO TO 11
C
C Compute constants.
C
      PIDNA = 4.D0*ATAN(1.D0)/DBLE(NA)
      S = 10.D0/(3.D0*DBLE(NR*NR*NA))
      TWONR = 2*NR
C
C Outer loop on radial index I.
C
      DO 2 I = 1,NR
        T = DBLE(2*I-1)/TWONR
        RI = T*D
        PHIT = (2.D0*T-3.D0)*T*T + 1.D0
        W(I) = S*DBLE(2*I-1)*PHIT
C
C Inner loop on angular coordinate J.
C
        DO 1 J = 1,NA
          AJ = DBLE(2*J-1)*PIDNA
          QX(I,J) = RI*COS(AJ)
          QY(I,J) = RI*SIN(AJ)
    1     CONTINUE
    2   CONTINUE
C
C No errors encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
      END
      SUBROUTINE GLIST (N,LIST,LPTR,LEND,LISTV,DXL,DYL, GX,
     .                  GY,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), LISTV(*), IER
      DOUBLE PRECISION DXL(*), DYL(*), GX(N), GY(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/13/96
C
C   Given a convexity-preserving triangulation of a set of
C nodes and data values, along with its dual, the gradient
C feasibility diagram, this subroutine returns a set of
C nodal gradients for which there exists a convex Hermite
C interpolant of the data values and gradients.
C
C   The nodal gradients are taken to be the centroids of
C the feasibility regions or truncated feasibility regions.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C           Note that, if N = 3, the three truncated feasi-
C           bility regions degenerate to the same point (the
C           gradient of the linear interpolant) which is
C           returned as the gradient of all three nodes.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMSHC.
C
C   The above parameters should be unaltered from the values
C input to Subroutine VLIST (which should be called to ob-
C tain LISTV, DXL, and DYL).
C
C       LISTV = Integer array containing feasibility region
C               vertex indexes (indexes to DXL and DYL)
C               stored in one-to-one correspondence with
C               LIST/LPTR entries.  Refer to Subroutine
C               VLIST.
C
C       DXL,DYL = Arrays of length NV = 2*N-2 containing
C                 the vertices (components of the gradients
C                 associated with triangles) and pseudo-
C                 vertices.  Refer to Subroutine VLIST.
C
C The above parameters are not altered by this routine.
C
C       GX,GY = Arrays of length at least N.
C
C On output:
C
C       GX,GY = Arrays of length N containing the components
C               of the nodal gradients (centroids of the
C               gradient feasibility regions or truncated
C               regions) unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C
C Modules required by GLIST:  None
C
C Intrinsic function called by GLIST:  DBLE
C
C***********************************************************
C
      DOUBLE PRECISION SX, SY
      INTEGER KV, LP, LP1L, LP2L, N1, N2, NN, NNB
C
C Local parameters:
C
C KV =    Index of a vertex of N1's gradient feasibility
C           region
C LP =    LIST index (pointer)
C LP2L =  Pointer to the last neighbor of N2
C LP1L =  Pointer to the last neighbor of N1
C N1,N2 = Nodal indexes
C NN =    Local copy of N
C NNB =   Number of vertices in the feasibility region
C           associated with node N1
C SX,SY = Components of the sum of vertices in the
C           feasibility region associated with N1
C
      NN = N
      IF (NN .LT. 3) GO TO 11
C
C Compute the centroids of the gradient feasibility regions.
C
C   Loop on nodes N1, accumulating the sums SX and SY of the
C     components of the NNB vertices of the associated
C     feasibility region.
C
      DO 2 N1 = 1,NN
        NNB = 0
        SX = 0.
        SY = 0.
        LP1L = LEND(N1)
        LP = LP1L
C
C Loop on neighbors of N1.
C
    1   LP = LPTR(LP)
          KV = LISTV(LP)
          NNB = NNB + 1
          SX = SX + DXL(KV)
          SY = SY + DYL(KV)
          IF (LP .NE. LP1L) GO TO 1
        IF (LIST(LP1L) .LT. 0) THEN
C
C   N1 is a boundary node.  Add in the contributions from
C     the pseudo-vertex associated with its first neighbor
C     N2.
C
          LP = LPTR(LP1L)
          N2 = LIST(LP)
          LP2L = LEND(N2)
          KV = LISTV(LP2L)
          NNB = NNB + 1
          SX = SX + DXL(KV)
          SY = SY + DYL(KV)
        ENDIF
C
C   Compute and store the centroid.
C
        GX(N1) = SX/DBLE(NNB)
        GY(N1) = SY/DBLE(NNB)
    2   CONTINUE
C
C No error encountered.
C
      IER = 0
      RETURN
C
C N < 3.
C
   11 IER = 1
      RETURN
      END
      SUBROUTINE INSERT (K,LP, LIST,LPTR,LNEW )
      INTEGER K, LP, LIST(*), LPTR(*), LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine inserts K as a neighbor of N1 following
C N2, where LP is the LIST pointer of N2 as a neighbor of
C N1.  Note that, if N2 is the last neighbor of N1, K will
C become the first neighbor (even if N1 is a boundary node).
C
C
C On input:
C
C       K = Index of the node to be inserted.
C
C       LP = LIST pointer of N2 as a neighbor of N1.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LNEW = Data structure defining the trian-
C                        gulation.
C
C On output:
C
C       LIST,LPTR,LNEW = Data structure updated with the
C                        addition of node K.
C
C Modules required by INSERT:  None
C
C***********************************************************
C
      INTEGER LSAV
C
      LSAV = LPTR(LP)
      LPTR(LP) = LNEW
      LIST(LNEW) = K
      LPTR(LNEW) = LSAV
      LNEW = LNEW + 1
      RETURN
      END
      SUBROUTINE INTADD (KK,I1,I2,I3, LIST,LPTR,LEND,LNEW )
      INTEGER KK, I1, I2, I3, LIST(*), LPTR(*), LEND(*),
     .        LNEW
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/22/91
C
C   This subroutine adds an interior node to a triangulation
C of a set of points in the plane.  The data structure is
C updated with the insertion of node KK into the triangle
C whose vertices are I1, I2, and I3.  No optimization of the
C triangulation is performed.
C
C
C On input:
C
C       KK = Index of the node to be inserted.  KK .GE. 1
C            and KK must not be equal to I1, I2, or I3.
C
C       I1,I2,I3 = Indexes of the counterclockwise-ordered
C                  sequence of vertices of a triangle which
C                  contains node KK.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND,LNEW = Data structure defining the
C                             triangulation.  Triangle
C                             (I1,I2,I3) must be included
C                             in the triangulation.
C
C On output:
C
C       LIST,LPTR,LEND,LNEW = Data structure updated with
C                             the addition of node KK.  KK
C                             will be connected to nodes I1,
C                             I2, and I3.
C
C Modules required by INTADD:  INSERT, LSTPTR
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER K, LP, N1, N2, N3
C
C Local parameters:
C
C K =        Local copy of KK
C LP =       LIST pointer
C N1,N2,N3 = Local copies of I1, I2, and I3
C
      K = KK
C
C Initialization.
C
      N1 = I1
      N2 = I2
      N3 = I3
C
C Add K as a neighbor of I1, I2, and I3.
C
      LP = LSTPTR(LEND(N1),N2,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
      LP = LSTPTR(LEND(N2),N3,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
      LP = LSTPTR(LEND(N3),N1,LIST,LPTR)
      CALL INSERT (K,LP, LIST,LPTR,LNEW )
C
C Add I1, I2, and I3 as neighbors of K.
C
      LIST(LNEW) = N1
      LIST(LNEW+1) = N2
      LIST(LNEW+2) = N3
      LPTR(LNEW) = LNEW + 1
      LPTR(LNEW+1) = LNEW + 2
      LPTR(LNEW+2) = LNEW
      LEND(K) = LNEW + 2
      LNEW = LNEW + 3
      RETURN
      END
      INTEGER FUNCTION JRAND (N, IX,IY,IZ )
      INTEGER N, IX, IY, IZ
C
C***********************************************************
C
C                                              From STRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   This function returns a uniformly distributed pseudo-
C random integer in the range 1 to N.
C
C
C On input:
C
C       N = Maximum value to be returned.
C
C N is not altered by this function.
C
C       IX,IY,IZ = Integer seeds initialized to values in
C                  the range 1 to 30,000 before the first
C                  call to JRAND, and not altered between
C                  subsequent calls (unless a sequence of
C                  random numbers is to be repeated by
C                  reinitializing the seeds).
C
C On output:
C
C       IX,IY,IZ = Updated integer seeds.
C
C       JRAND = Random integer in the range 1 to N.
C
C Reference:  B. A. Wichmann and I. D. Hill, "An Efficient
C             and Portable Pseudo-random Number Generator",
C             Applied Statistics, Vol. 31, No. 2, 1982,
C             pp. 188-190.
C
C Modules required by JRAND:  None
C
C Intrinsic functions called by JRAND:  INT, MOD, REAL
C
C***********************************************************
C
      REAL U, X
C
C Local parameters:
C
C U = Pseudo-random number uniformly distributed in the
C     interval (0,1).
C X = Pseudo-random number in the range 0 to 3 whose frac-
C       tional part is U.
C
      IX = MOD(171*IX,30269)
      IY = MOD(172*IY,30307)
      IZ = MOD(170*IZ,30323)
      X = (REAL(IX)/30269.) + (REAL(IY)/30307.) +
     .    (REAL(IZ)/30323.)
      U = X - INT(X)
      JRAND = REAL(N)*U + 1.
      RETURN
      END
      LOGICAL FUNCTION LEFT (X1,Y1,X2,Y2,X0,Y0)
      DOUBLE PRECISION X1, Y1, X2, Y2, X0, Y0
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function determines whether node N0 is to the left
C or to the right of the line through N1-N2 as viewed by an
C observer at N1 facing N2.
C
C
C On input:
C
C       X1,Y1 = Coordinates of N1.
C
C       X2,Y2 = Coordinates of N2.
C
C       X0,Y0 = Coordinates of N0.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LEFT = TRUE if and only if (X0,Y0) is on or to the
C              left of the directed line N1->N2.
C
C Modules required by LEFT:  None
C
C***********************************************************
C
      DOUBLE PRECISION DX1, DY1, DX2, DY2
C
C Local parameters:
C
C DX1,DY1 = X,Y components of the vector N1->N2
C DX2,DY2 = X,Y components of the vector N1->N0
C
      DX1 = X2-X1
      DY1 = Y2-Y1
      DX2 = X0-X1
      DY2 = Y0-Y1
C
C If the sign of the vector cross product of N1->N2 and
C   N1->N0 is positive, then sin(A) > 0, where A is the
C   angle between the vectors, and thus A is in the range
C   (0,180) degrees.
C
      LEFT = DX1*DY2 .GE. DX2*DY1
      RETURN
      END
      SUBROUTINE LGRAD (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3, FX,FY,
     .                  SA)
      DOUBLE PRECISION X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3,
     .                 FX, FY, SA
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/13/96
C
C   Given three points defining a non-vertical plane, and
C hence a bivariate linear function F, this subroutine re-
C turns the gradient of F and the signed area of the
C triangle obtained by projecting the points onto the x-y
C plane.
C
C
C On input:
C
C       X1,Y1,Z1 = Cartesian coordinates of the first vertex
C                  (X1,Y1) and associated function value
C                  Z1 = F(X1,Y1).
C
C       X2,Y2,Z2 = Cartesian coordinates of the second ver-
C                  tex (X2,Y2) and associated function value
C                  Z2 = F(X2,Y2).
C
C       X3,Y3,Z3 = Cartesian coordinates of the third vertex
C                  (X3,Y3) and associated function value
C                  Z3 = F(X3,Y3).
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       FX,FY = Components of the gradient of F -- partial
C               derivatives with respect to x and y, respec-
C               tively, unless SA = 0, in which case the
C               gradient is not defined.
C
C       SA = Signed area of the triangle with vertices
C            (X1,Y1), (X2,Y2), and (X3,Y3).  SA > 0 if and
C            only if the vertices are in counterclockwise
C            order:  (X3,Y3) is strictly to the left of the
C            directed line from (X1,Y1) to (X2,Y2).
C
C Modules required by LGRAD:  None
C
C***********************************************************
C
      DOUBLE PRECISION A, DX1, DX2, DY1, DY2, DZ1, DZ2
C
C Local parameters:
C
C A       = Twice the signed triangle area:  z-component of
C             (V1-V3) X (V2-V3) for vertices V1, V2, and V3
C DX1,DY1 = Components of V1-V3
C DX2,DY2 = Components of V2-V3
C DZ1,DZ2 = Z1-Z3 and Z2-Z3, respectively
C
      DX1 = X1-X3
      DY1 = Y1-Y3
      DX2 = X2-X3
      DY2 = Y2-Y3
C
      A = DX1*DY2 - DX2*DY1
      SA = A/2.D0
      IF (A .EQ. 0.) RETURN
C
      DZ1 = Z1-Z3
      DZ2 = Z2-Z3
      FX = (DZ1*DY2 - DZ2*DY1)/A
      FY = (DX1*DZ2 - DX2*DZ1)/A
      RETURN
      END
      INTEGER FUNCTION LSTPTR (LPL,NB,LIST,LPTR)
      INTEGER LPL, NB, LIST(*), LPTR(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function returns the index (LIST pointer) of NB in
C the adjacency list for N0, where LPL = LEND(N0).
C
C
C On input:
C
C       LPL = LEND(N0)
C
C       NB = Index of the node whose pointer is to be re-
C            turned.  NB must be connected to N0.
C
C       LIST,LPTR = Data structure defining the triangula-
C                   tion.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       LSTPTR = Pointer such that LIST(LSTPTR) = NB or
C                LIST(LSTPTR) = -NB, unless NB is not a
C                neighbor of N0, in which case LSTPTR = LPL.
C
C Modules required by LSTPTR:  None
C
C***********************************************************
C
      INTEGER LP, ND
C
C Local parameters:
C
C LP = LIST pointer
C ND = Nodal index
C
      LP = LPTR(LPL)
    1 ND = LIST(LP)
        IF (ND .EQ. NB) GO TO 2
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 1
C
    2 LSTPTR = LP
      RETURN
      END
      INTEGER FUNCTION NBCNT (LPL,LPTR)
      INTEGER LPL, LPTR(*)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function returns the number of neighbors of a node
C N0 in a triangulation.
C
C
C On input:
C
C       LPL = LIST pointer to the last neighbor of N0 --
C             LPL = LEND(N0).
C
C       LPTR = Array of pointers associated with LIST.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       NBCNT = Number of neighbors of N0.
C
C Modules required by NBCNT:  None
C
C***********************************************************
C
      INTEGER K, LP
C
C Local parameters:
C
C K =  Counter for computing the number of neighbors
C LP = LIST pointer
C
      LP = LPL
      K = 1
C
    1 LP = LPTR(LP)
        IF (LP .EQ. LPL) GO TO 2
        K = K + 1
        GO TO 1
C
    2 NBCNT = K
      RETURN
      END
      SUBROUTINE OPTIM (X,Y,Z,NE, LIST,LPTR,LEND,NIT,
     .                  IWK, IER)
      INTEGER NE, LIST(*), LPTR(*), LEND(*), NIT, IWK(2,NE),
     .        IER
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/12/96
C
C   Given a set of NE triangulation edges, this subroutine
C optimizes the portion of the triangulation consisting of
C the quadrilaterals (pairs of adjacent triangles) which
C have the edges as diagonals by applying the swap test
C and appropriate swaps to the edges.
C
C   An iteration consists of applying the swap test and
C swaps to all NE edges in the order in which they are
C stored.  The iteration is repeated until no swap occurs
C or NIT iterations have been performed.  The bound on the
C number of iterations may be necessary to prevent an
C infinite loop caused by cycling (reversing the effect of a
C previous swap) due to floating point inaccuracy when four
C or more data points are nearly coplanar.
C
C
C On input:
C
C       X,Y,Z = Arrays containing the nodal coordinates and
C               data values.
C
C       NE = Number of edges in the set.  NE .GE. 0.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMSHC.
C
C       NIT = Maximum number of iterations to be performed.
C             A reasonable value is 3*NE.  NIT .GE. 1.
C
C       IWK = Integer array dimensioned 2 by NE containing
C             the nodal indexes of the edge endpoints (pairs
C             of endpoints are stored in columns).
C
C On output:
C
C       LIST,LPTR,LEND = Updated triangulation data struc-
C                        ture reflecting the swaps.
C
C       NIT = Number of iterations performed.
C
C       IWK = Endpoint indexes of the new set of edges
C             reflecting the swaps.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if a swap occurred on the last of
C                     MAXIT iterations, where MAXIT is the
C                     value of NIT on input.  The new set
C                     of edges in not necessarily optimal
C                     in this case.
C             IER = 2 if NE < 0 or NIT < 1 on input.
C             IER = 3 if IWK(2,I) is not a neighbor of
C                     IWK(1,I) for some I in the range 1
C                     to NE.  A swap may have occurred in
C                     this case.
C             IER = 4 if an error flag was returned by
C                     Function SWPTC:  the triangulation
C                     is not convex.
C
C Modules required by OPTIM:  LSTPTR, SWAP, SWPTC
C
C Intrinsic function called by OPTIM:  ABS
C
C***********************************************************
C
      INTEGER I, IERR, IO1, IO2, ITER, LP, LP21, LPL, LPP,
     .        MAXIT, N1, N2, NNE
      LOGICAL SWPTC
      LOGICAL SWP
C
C Local parameters:
C
C I =       Column index for IWK
C IERR =    Error flag for calls to SWPTC
C IO1,IO2 = Nodal indexes of the endpoints of an edge in IWK
C ITER =    Iteration count
C LP =      LIST pointer
C LP21 =    Parameter returned by SWAP (not used)
C LPL =     Pointer to the last neighbor of IO1
C LPP =     Pointer to the node preceding IO2 as a neighbor
C             of IO1
C MAXIT =   Input value of NIT
C N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1,
C             respectively
C NNE =     Local copy of NE
C SWP =     Flag set to TRUE iff a swap occurs in the
C             optimization loop
C
      NNE = NE
      MAXIT = NIT
      IF (NNE .LT. 0  .OR.  MAXIT .LT. 1) GO TO 7
C
C Initialize iteration count ITER and test for NE = 0.
C
      ITER = 0
      IF (NNE .EQ. 0) GO TO 5
C
C Top of loop --
C   SWP = TRUE iff a swap occurred in the current iteration.
C
    1 IF (ITER .EQ. MAXIT) GO TO 6
      ITER = ITER + 1
      SWP = .FALSE.
C
C   Inner loop on edges IO1-IO2 --
C
      DO 4 I = 1,NNE
        IO1 = IWK(1,I)
        IO2 = IWK(2,I)
C
C   Set N1 and N2 to the nodes opposite IO1->IO2 and
C     IO2->IO1, respectively.  Determine the following:
C
C     LPL = pointer to the last neighbor of IO1,
C     LP = pointer to IO2 as a neighbor of IO1, and
C     LPP = pointer to the node N2 preceding IO2.
C
        LPL = LEND(IO1)
        LPP = LPL
        LP = LPTR(LPP)
    2   IF (LIST(LP) .EQ. IO2) GO TO 3
          LPP = LP
          LP = LPTR(LPP)
          IF (LP .NE. LPL) GO TO 2
C
C   IO2 should be the last neighbor of IO1.  Test for no
C     edge and bypass the swap test if IO1 is a boundary
C     node.
C
        IF (ABS(LIST(LP)) .NE. IO2) GO TO 8
        IF (LIST(LP) .LT. 0) GO TO 4
C
C   Store N1 and N2, or bypass the swap test if IO1 is a
C     boundary node and IO2 is its first neighbor.
C
    3   N2 = LIST(LPP)
        IF (N2 .LT. 0) GO TO 4
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
C
C   Test IO1-IO2 for a swap, and update IWK if necessary.
C
        IF ( SWPTC(N1,N2,IO1,IO2,X,Y,Z, IERR) ) THEN
          SWP = .TRUE.
          CALL SWAP (N1,N2,IO1,IO2, LIST,LPTR,LEND, LP21)
          IWK(1,I) = N1
          IWK(2,I) = N2
        ENDIF
        IF (IERR .NE. 0) GO TO 9
    4   CONTINUE
      IF (SWP) GO TO 1
C
C Successful termination.
C
    5 NIT = ITER
      IER = 0
      RETURN
C
C MAXIT iterations performed without convergence.
C
    6 NIT = MAXIT
      IER = 1
      RETURN
C
C Invalid input parameter.
C
    7 NIT = 0
      IER = 2
      RETURN
C
C IO2 is not a neighbor of IO1.
C
    8 NIT = ITER
      IER = 3
      RETURN
C
C Non-convex triangulation.
C
    9 NIT = ITER
      IER = 4
      RETURN
      END
      SUBROUTINE PLTCNT (LUN,PLTSIZ,NX,NY,PX,PY,PZ,NCON,IWK,
     .                   XC,YC, IER)
      INTEGER LUN, NX, NY, NCON, IWK(*), IER
      DOUBLE PRECISION PLTSIZ, PX(NX), PY(NY), PZ(NX,NY),
     .                 XC(*), YC(*)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   04/12/97
C
C   Given a set of function values PZ = F(X,Y) at the ver-
C tices of an NX by NY rectangular grid, this subroutine
C creates a level-2 Encapsulated PostScript (EPS) file
C containing a contour plot of the piecewise bilinear inter-
C polant of the function values.
C
C   The accuracy of the contour lines increases with the
C number of grid points.  Refer to Subroutine CNTOUR for
C further details.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  A window containing
C                the plot is mapped, with aspect ratio
C                preserved, to a rectangular viewport with
C                maximum side-length PLTSIZ.  The viewport
C                is centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  1.0 .LE. PLTSIZ
C                .LE. 7.5.
C
C       NX = Number of grid points in the x direction.
C	     NX .GE. 2.
C
C	NY = Number of grid points in the y direction.
C	     NY .GE. 2.
C
C       PX = Array of length NX containing a strictly in-
C            creasing sequence of values.
C
C       PY = Array of length NY containing a strictly in-
C            creasing sequence of values.
C
C       PZ = Array of function values at the vertices of the
C            rectangular grid.  PZ(I,J) = F(PX(I),PY(J)) for
C            I = 1,...,NX and J = 1,...,NY.
C
C       NCON = Number of contour values.  The contour values
C              are uniformly distributed over the range of
C              PZ values.  NCON .GE. 1.
C
C The above parameters are not altered by this routine.
C
C       IWK = Integer array of length at least 1.5*NX*NY to
C             be used as work space.
C
C       XC,YC = Arrays of length at least 2.5*NX*NY to
C               be used as work space.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, NX, NY, or NCON is
C                     outside its valid range.
C             IER = 2 if PX or PY is not strictly
C                     increasing.
C             IER = 3 if the range of PZ values has zero
C                     width (F is constant).
C             IER = 4 if an error was encountered in writing
C                     to unit LUN.
C             IER = 5 if an unexpected error flag was re-
C                     turned by Subroutine CNTOUR.  This
C                     should not occur.
C
C   In the unlikely event of an internal failure, a message
C is printed on the standard output device.  IER may be 0
C in this case.
C
C Module required by PLTCNT:  CNTOUR
C
C Intrinsic functions called by PLTCNT:  CHAR, DBLE
C
C***********************************************************
C
      DOUBLE PRECISION CVAL, DX, DY, DZ, PZIJ, R, SFX, SFY,
     .                 T, TX, TY, ZMAX, ZMIN
      INTEGER I, IC, IERR, IH, IPX1, IPX2, IPY1, IPY2, IW,
     .        J, K, KV, LC, NC, NCMAX
C
C Local parameters:
C
C CVAL =      Contour value between ZMIN and ZMAX
C DX =        Window width PX(NX)-PX(1)N
C DY =        Window height PY(NY)-PY(1)
C DZ =        Interval between contour values:
C               (ZMAX-ZMIN)/(NCON+1)
C I,J =       Row and column indexes for PZ
C IC =        Index (for IWK) of a contour line associated
C               with contour value CVAL:  1 to NC
C IERR =      Error flag for calls to CNTOUR
C IH =        Height of the bounding box (viewport) in
C               points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box
C IW =        Width of the bounding box in points
C K =         Index (for XC and YC) of a point on a contour
C               line
C KV =        DO-loop index for loop on contour values
C LC =        Length of arrays XC and YC
C NC =        Number of contour lines associated with
C               contour value CVAL
C NCMAX =     Maximum allowable value of NC
C PZIJ =      PZ(I,J)
C R =         Aspect ratio DX/DY
C SFX,SFY =   Scale factors for mapping window coordinates
C               to viewport coordinates
C T =         Temporary variable
C TX,TY =     Translation vector for mapping window coordi-
C               nates to viewport coordinates
C ZMIN,ZMAX = Minimum and maximum of the PZ values
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.D0  .OR.  PLTSIZ .GT. 7.5D0  .OR.
     .    NX .LT. 2  .OR.  NY .LT. 2  .OR.  NCON .LT. 1)
     .  GO TO 11
C
C Compute the aspect ratio of the window.
C
      DX = PX(NX) - PX(1)
      DY = PY(NY) - PY(1)
      IF (DX .EQ. 0.  .OR.  DY .EQ. 0.) GO TO 12
      R = DX/DY
C
C Compute the range of PZ values and the interval between
C   contour values.
C
      ZMIN = PZ(1,1)
      ZMAX = ZMIN
      DO 2 J = 1,NY
        DO 1 I = 1,NX
          PZIJ = PZ(I,J)
          IF (PZIJ .LT. ZMIN) ZMIN = PZIJ
          IF (PZIJ .GT. ZMAX) ZMAX = PZIJ
    1     CONTINUE
    2   CONTINUE
      DZ = (ZMAX-ZMIN)/DBLE(NCON+1)
      IF (DZ .LE. 0.) GO TO 13
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box
C   (the viewport).  The coordinates, specified in default
C   user space units (points, at 72 points/inch with origin
C   at the lower left corner of the page), are chosen to
C   preserve the aspect ratio R, and to center the plot on
C   the 8.5 by 11 inch page.  The center of the page is
C   (306,396), and T = PLTSIZ/2 in points.
C
      T = 36.D0*PLTSIZ
      IF (R .GE. 1.D0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.
C
      WRITE (LUN,100,ERR=14) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Contour Plot'/
     .        '%%Creator:  CSRFPACK'/
     .        '%%EndComments')
C
C Draw the bounding box.
C
      WRITE (LUN,110,ERR=14) IPX1, IPY1
      WRITE (LUN,120,ERR=14) IPX1, IPY2
      WRITE (LUN,120,ERR=14) IPX2, IPY2
      WRITE (LUN,120,ERR=14) IPX2, IPY1
      WRITE (LUN,130,ERR=14)
      WRITE (LUN,140,ERR=14)
  110 FORMAT (2I4,' moveto')
  120 FORMAT (2I4,' lineto')
  130 FORMAT ('closepath')
  140 FORMAT ('stroke')
C
C Set up a mapping from the window to the viewport.
C
      IW = IPX2 - IPX1
      IH = IPY2 - IPY1
      SFX = DBLE(IW)/DX
      SFY = DBLE(IH)/DY
      TX = IPX1 - SFX*PX(1)
      TY = IPY1 - SFY*PY(1)
      WRITE (LUN,150,ERR=14) TX, TY, SFX, SFY
  150 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C Set the line thickness to 2 points.  (Since the scale
C   factors are applied to everything, the width must be
C   specified in world coordinates.)
C
      T = 4.D0/(SFX+SFY)
      WRITE (LUN,160,ERR=14) T
  160 FORMAT (F12.6,' setlinewidth')
C
C Compute parameters for CNTOUR:
C
C   NCMAX = Maximum allowable number of contour lines
C           associated with each contour value.
C   LC = Length of arrays XC and YC and maximum allowable
C        number of points defining all the contour lines
C        associated with a contour value.
C
      NCMAX = (NX*NY+1)/2
      LC = 2*(NX-1)*(NY-1) + NCMAX
C
C Loop on contour values CVAL uniformly spaced in the open
C   interval (ZMIN,ZMAX).
C
      CVAL = ZMIN
      DO 5 KV = 1,NCON
        CVAL = CVAL + DZ
C
C Compute a sequence of NC contour lines associated with
C   F = CVAL.  For IC = 1 to NC, IWK(IC) is the index (for
C   XC and YC) of the last point of contour IC.
C
        CALL CNTOUR (NX,NY,PX,PY,PZ,CVAL,LC,NCMAX,
     .               IWK(NCMAX+1), XC,YC,IWK,NC,IERR)
        IF (IERR .EQ. 2) GO TO 12
        IF (IERR .NE. 0) GO TO 15
C
C Draw the NC contours.
C
        IC = 0
        K = 0
    3   IC = IC + 1
          K = K + 1
C
C   Create a path consisting of contour IC.
C
          WRITE (LUN,170,ERR=14) XC(K), YC(K)
  170     FORMAT (2F12.6,' moveto')
    4     K = K + 1
            WRITE (LUN,180,ERR=14) XC(K), YC(K)
  180       FORMAT (2F12.6,' lineto')
            IF (K .NE. IWK(IC)) GO TO 4
C
C   Paint the path.
C
          WRITE (LUN,140,ERR=14)
          IF (IC .NE. NC) GO TO 3
    5   CONTINUE
C
C Output the showpage command and end-of-file indicator.
C
      WRITE (LUN,200,ERR=14)
  200 FORMAT ('showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,210,ERR=14) CHAR(4)
  210 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C PX or PY is not strictly increasing.
C
   12 IER = 2
      RETURN
C
C DZ = 0.
C
   13 IER = 3
      RETURN
C
C Error writing to unit LUN.
C
   14 IER = 4
      RETURN
C
C Error flag returned by CNTOUR.
C
   15 IER = 5
      RETURN
      END
      SUBROUTINE PLTGR (LUN,PLTSIZ,WX1,WX2,WY1,WY2,N,LIST,
     .                  LPTR,LEND,NV,LISTV,DXL,DYL,GX,GY,
     .                  TITLE,NUMBR, IER)
      CHARACTER*(*) TITLE
      INTEGER LUN, N, LIST(*), LPTR(*), LEND(N), NV,
     .        LISTV(*), IER
      LOGICAL NUMBR
      DOUBLE PRECISION PLTSIZ, WX1, WX2, WY1, WY2, DXL(NV),
     .                 DYL(NV), GX(N), GY(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/98
C
C   Given a convexity-preserving triangulation of a set of
C nodes and data values, along with its straight-line dual,
C the gradient feasibility diagram, and a set of nodal grad-
C ients, this subroutine creates a level-2 Encapsulated
C PostScript (EPS) file containing a plot of the feasibility
C diagram and nodal gradients.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  The window is mapped,
C                with aspect ratio preserved, to a rectangu-
C                lar viewport with maximum side-length equal
C                to .88*PLTSIZ (leaving room for labels out-
C                side the viewport).  The viewport is
C                centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  1.0 .LE. PLTSIZ
C                .LE. 8.5.
C
C       WX1,WX2,WY1,WY2 = Parameters defining a rectangular
C                         window against which the feasibil-
C                         ity diagram is clipped.  Only the
C                         portion of the diagram that lies
C                         in the window is drawn.  (WX1,WY1)
C                         and (WX2,WY2) are the lower left
C                         and upper right corners, respec-
C                         tively.  WX1 < WX2 and WY1 < WY2.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMSHC.
C
C       NV = Number of feasibility region vertices, includ-
C            ing pseudo-vertices associated with boundary
C            edges:  NV = NT+NB = 2*N-2, where NT = 2*N-NB-2
C            is the number of triangles.
C
C       LISTV = Integer array containing feasibility region
C               vertex indexes (indexes to DXL and DYL)
C               stored in one-to-one correspondence with
C               LIST/LPTR entries.  Refer to Subroutine
C               VLIST.
C
C       DXL,DYL = Arrays of length NV containing the
C                 vertices (components of the gradients
C                 associated with triangles) and pseudo-
C                 vertices.  Refer to Subroutine VLIST.
C
C       GX,GY = Arrays of length N containing the components
C               of the nodal gradients (centroids of the
C               gradient feasibility regions or truncated
C               regions).  Refer to Subroutine GLIST.
C
C       TITLE = Type CHARACTER variable or constant contain-
C               ing a string to be centered above the plot.
C               The string must be enclosed in parentheses;
C               i.e., the first and last characters must be
C               '(' and ')', respectively, but these are not
C               displayed.  TITLE may have at most 80 char-
C               acters including the parentheses.
C
C       NUMBR = Option indicator:  If NUMBR = TRUE, the
C               nodal indexes are plotted next to the nodal
C               gradients.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, N, or NV is outside
C                     its valid range.
C             IER = 2 if WX1 >= WX2 or WY1 >= WY2.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C   Various plotting options can be controlled by altering
C the data statement below.
C
C Modules required by PLTGR:  None
C
C Intrinsic functions called by PLTGR:  CHAR, DBLE, NINT
C
C***********************************************************
C
      DOUBLE PRECISION DASHL, DX, DY, FSIZN, FSIZT, R, SFX,
     .                 SFY, T, TX, TY, X1, X2, XB1, XB2, XP,
     .                 Y1, Y2, YB1, YB2, YP
      INTEGER IH, IPX1, IPX2, IPY1, IPY2, IW, KV0, KV1, KV2,
     .        LP, LP1L, LPL, N0, N1
      LOGICAL ANNOT
C
      DATA    ANNOT/.TRUE./,  DASHL/4.D0/,  FSIZN/10.D0/,
     .        FSIZT/16.D0/
C
C Local parameters:
C
C ANNOT =     Logical variable with value TRUE iff the plot
C               is to be annotated with the values of WX1,
C               WX2, WY1, and WY2
C DASHL =     Length (in points, at 72 points per inch) of
C               dashes and spaces in a dashed line pattern
C               used for connecting pseudo-vertices
C DX =        Window width WX2-WX1
C DY =        Window height WY2-WY1
C FSIZN =     Font size in points for labeling nodal gradi-
C               ents with their indexes if NUMBR = TRUE
C FSIZT =     Font size in points for the title (and
C               annotation if ANNOT = TRUE)
C IH =        Height of the viewport in points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box or viewport
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box or viewport
C IW =        Width of the viewport in points
C KV0 =       Index of the first vertex of an unbounded
C               feasibility region
C KV1,KV2 =   Endpoint indexes of a feasibility region edge
C LP =        LIST index (pointer) of a neighbor of N0
C LP1L =      Pointer to the last neighbor of N1
C LPL =       Pointer to the last neighbor of N0
C N0 =        Nodal index and DO-loop index
C N1 =        Index of the first neighbor of N0
C R =         Aspect ratio DX/DY
C SFX,SFY =   Scale factors for mapping world coordinates
C               (window coordinates in [WX1,WX2] X [WY1,WY2])
C               to viewport coordinates in [IPX1,IPX2] X
C               [IPY1,IPY2]
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C X1,Y1 =     Components of vertex KV1:  DXL(KV1), DYL(KV1)
C X2,Y2 =     Components of vertex KV2:  DXL(KV2), DYL(KV2)
C XB1,YB1 =   Components of the lower left corner of the
C               viewport in world coordinates
C XB2,YB2 =   Components of the upper right corner of the
C               viewport in world coordinates
C XP,YP     = Components of the point at which an edge
C               of an unbounded region intersects the
C               viewport, or components of a nodal gradient
C               and/or label location
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.D0  .OR.  PLTSIZ .GT. 7.5D0  .OR.
     .    N .LT. 3  .OR.  NV .NE. 2*N-2) GO TO 11
C
C Compute the aspect ratio of the window.
C
      DX = WX2 - WX1
      DY = WY2 - WY1
      IF (DX .LE. 0.  .OR.  DY .LE. 0.) GO TO 12
      R = DX/DY
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box.
C   The coordinates, specified in default user space units
C   (points, at 72 points/inch with origin at the lower
C   left corner of the page), are chosen to preserve the
C   aspect ratio R, and to center the plot on the 8.5 by 11
C   inch page.  The center of the page is (306,396), and
C   T = PLTSIZ/2 in points.
C
      T = 36.D0*PLTSIZ
      IF (R .GE. 1.D0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.
C
      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Triangulation Dual'/
     .        '%%Creator:  CSRFPACK'/
     .        '%%EndComments')
C
C Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
C   of a viewport obtained by shrinking the bounding box by
C   12% in each dimension.
C
      IW = NINT(0.88D0*DBLE(IPX2-IPX1))
      IH = NINT(0.88D0*DBLE(IPY2-IPY1))
      IPX1 = 306 - IW/2
      IPX2 = 306 + IW/2
      IPY1 = 396 - IH/2
      IPY2 = 396 + IH/2
C
C Set the line thickness to 2 points, and draw the
C   viewport boundary.
C
      T = 2.D0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IPX1, IPY1
      WRITE (LUN,130,ERR=13) IPX1, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY1
      WRITE (LUN,140,ERR=13)
      WRITE (LUN,150,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT (2I4,' moveto')
  130 FORMAT (2I4,' lineto')
  140 FORMAT ('closepath')
  150 FORMAT ('stroke')
C
C Set up a mapping from the window to the viewport.
C
      SFX = DBLE(IW)/DX
      SFY = DBLE(IH)/DY
      TX = IPX1 - SFX*WX1
      TY = IPY1 - SFY*WY1
      WRITE (LUN,160,ERR=13) TX, TY, SFX, SFY
  160 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C Compute the viewport corners (XB1,YB1) and (XB2,YB2)
C   in world coordinates (by applying the inverse window-
C   to-viewport mapping to (IPX1,IPY1) and (IPX2,IPY2)).
C
      XB1 = (IPX1-TX)/SFX
      XB2 = (IPX2-TX)/SFX
      YB1 = (IPY1-TY)/SFY
      YB2 = (IPY2-TY)/SFY
C
C Set the line thickness to 1 point.  (Since the scale
C   factors are applied to everything, the width must be
C   specified in world coordinates.)
C
      T = 2.D0/(SFX+SFY)
      WRITE (LUN,110,ERR=13) T
C
C Save the current graphics state, and set the clip path to
C   the boundary of the window.
C
      WRITE (LUN,170,ERR=13)
      WRITE (LUN,180,ERR=13) WX1, WY1
      WRITE (LUN,190,ERR=13) WX2, WY1
      WRITE (LUN,190,ERR=13) WX2, WY2
      WRITE (LUN,190,ERR=13) WX1, WY2
      WRITE (LUN,200,ERR=13)
  170 FORMAT ('gsave')
  180 FORMAT (2F12.6,' moveto')
  190 FORMAT (2F12.6,' lineto')
  200 FORMAT ('closepath clip newpath')
C
C Set T to the world-coordinate equivalent of DASHL
C   (for drawing the boundary of the truncated feasibility
C   diagram).
C
      T = 2.D0*DASHL/(SFX+SFY)
C
C Loop on nodes N0 corresponding to feasibility regions.
C   LPL indexes the last neighbor of N0.
C
      DO 3 N0 = 1,N
        LPL = LEND(N0)
C
C Set KV2 to the first vertex index.
C
        IF (LIST(LPL) .GT. 0) THEN
C
C   N0 is an interior node.  The first vertex coincides with
C     the last.
C
          KV2 = LISTV(LPL)
        ELSE
C
C   N0 is a boundary node.  The first vertex is the last
C     vertex of the first neighbor N1 of N0.  A copy of KV2
C     is saved in KV0.
C
          LP = LPTR(LPL)
          N1 = LIST(LP)
          LP1L = LEND(N1)
          KV2 = LISTV(LP1L)
          KV0 = KV2
        ENDIF
C
C Loop on neighbors NB of N0.  For each triangulation edge
C   N0-NB, set KV1-KV2 to the corresponding feasibility
C   region edge.
C
        LP = LPL
    1   LP = LPTR(LP)
          KV1 = KV2
          KV2 = LISTV(LP)
C
C Add edge KV1-KV2 to the path iff KV2 > KV1.
C
          IF (KV2 .GT. KV1) WRITE (LUN,210,ERR=13) DXL(KV1),
     .                        DYL(KV1), DXL(KV2), DYL(KV2)
  210     FORMAT (2F12.6,' moveto',2F12.6,' lineto')
          IF (LP .NE. LPL) GO TO 1
        IF (LIST(LPL) .LT. 0) THEN
C
C N0 is a boundary node.  Paint the current path, draw a
C   dashed line between KV2 and KV0, and restore the solid
C   line attribute.
C
          WRITE (LUN,150,ERR=13)
          WRITE (LUN,220,ERR=13) T, T, DXL(KV2), DYL(KV2),
     .                           DXL(KV0), DYL(KV0)
  220     FORMAT ('[',2F12.6,'] 0 setdash'/
     .            2F12.6,' moveto',2F12.6,' lineto'/
     .            '  stroke'/
     .            '  [] 0 setdash')
C
C Extend KV1->KV2 to the boundary of the viewport.
C
          X1 = DXL(KV1)
          Y1 = DYL(KV1)
          X2 = DXL(KV2)
          Y2 = DYL(KV2)
          IF (X1 .NE. X2) THEN
C
C   Find the intersection (XP,YP) with the left or right
C     boundary.
C
            IF (X1 .LT. X2) THEN
              XP = XB2
            ELSE
              XP = XB1
            ENDIF
            YP = Y1 + ((XP-X1)/(X2-X1))*(Y2-Y1)
C
C   Test for the point of intersection inside the viewport.
C
            IF (YP .GE. YB1  .AND.  YP .LE. YB2) GO TO 2
          ENDIF
C
C   Find the intersection with the top or bottom boundary.
C
          IF (Y1 .LT. Y2) THEN
            YP = YB2
          ELSE
            YP = YB1
          ENDIF
          XP = X1 + ((YP-Y1)/(Y2-Y1))*(X2-X1)
C
C   Draw the line segment from KV2 to (XP,YP).
C
    2     WRITE (LUN,210,ERR=13) X2, Y2, XP, YP
          WRITE (LUN,150,ERR=13)
        ENDIF
    3   CONTINUE
C
C Paint the path and restore the saved graphics state (with
C   no clip path).
C
      WRITE (LUN,150,ERR=13)
      WRITE (LUN,230,ERR=13)
  230 FORMAT ('grestore')
      IF (NUMBR) THEN
C
C Nodal gradients in the window are to be labeled with their
C   indexes.  Convert FSIZN from points to world coordinates
C   and output the commands to select a font and scale it.
C
        T = FSIZN*2.D0/(SFX+SFY)
        WRITE (LUN,240,ERR=13) T
  240   FORMAT ('/Helvetica findfont'/
     .          F12.6,' scalefont setfont')
      ENDIF
C
C Draw dots (degenerate paths with a rounded line cap)
C   corresponding to the nodal gradients (XP,YP).
C
      T = 8.D0/(SFX+SFY)
      WRITE (LUN,250,ERR=13) T
  250 FORMAT (F12.6,' setlinewidth'/
     .        '1 setlinecap')
      DO 4 N0 = 1,N
        XP = GX(N0)
        YP = GY(N0)
        IF (XP .LT. WX1  .OR.  XP .GT. WX2  .OR.
     .      YP .LT. WY1  .OR.  YP .GT. WY2) GO TO 4
        WRITE (LUN,210,ERR=13) XP, YP, XP, YP
        IF (NUMBR) THEN
C
C   Draw the label N0.  The first character will have its
C     lower left corner about one character width to the
C     right of the nodal position.
C
          WRITE (LUN,260,ERR=13) N0
  260     FORMAT ('(',I3,') show')
        ENDIF
    4   CONTINUE
C
C Paint the path.
C
      WRITE (LUN,150,ERR=13)
C
C Convert FSIZT from points to world coordinates, and output
C   the commands to select a font and scale it.
C
      T = FSIZT*2.D0/(SFX+SFY)
      WRITE (LUN,240,ERR=13) T
C
C Display TITLE centered above the plot:
C
      YP = WY2 + 3.D0*T
      WRITE (LUN,270,ERR=13) TITLE, (WX1+WX2)/2.D0, YP
  270 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6,
     .        ' add ',F12.6,' moveto')
      WRITE (LUN,280,ERR=13) TITLE
  280 FORMAT (A80/'  show')
      IF (ANNOT) THEN
C
C Display the window extrema below the plot.
C
        XP = WX1
        YP = WY1 - 100.D0/(SFX+SFY)
        WRITE (LUN,180,ERR=13) XP, YP
        WRITE (LUN,290,ERR=13) WX1, WX2
        YP = YP - 2.D0*T
        WRITE (LUN,300,ERR=13) XP, YP, WY1, WY2
  290   FORMAT ('(Window:   WX1 = ',E9.3,',   WX2 = ',E9.3,
     .          ') show')
  300   FORMAT ('(Window:  ) stringwidth pop ',F12.6,' add',
     .          F12.6,' moveto'/
     .          '( WY1 = ',E9.3,',   WY2 = ',E9.3,') show')
      ENDIF
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,310,ERR=13)
  310 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,320,ERR=13) CHAR(4)
  320 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C DX or DY is not positive.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE PLTTR (LUN,PLTSIZ,WX1,WX2,WY1,WY2,N,X,Y,
     .                  LIST,LPTR,LEND,TITLE,NUMBR, IER)
      CHARACTER*(*) TITLE
      INTEGER LUN, N, LIST(*), LPTR(*), LEND(N), IER
      LOGICAL NUMBR
      DOUBLE PRECISION PLTSIZ, WX1, WX2, WY1, WY2, X(N),
     .                 Y(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/15/98
C
C   This subroutine creates a level-2 Encapsulated Post-
C Script (EPS) file containing a triangulation plot.
C
C
C On input:
C
C       LUN = Logical unit number in the range 0 to 99.
C             The unit should be opened with an appropriate
C             file name before the call to this routine.
C
C       PLTSIZ = Plot size in inches.  The window is mapped,
C                with aspect ratio preserved, to a rectangu-
C                lar viewport with maximum side-length equal
C                to .88*PLTSIZ (leaving room for labels out-
C                side the viewport).  The viewport is
C                centered on the 8.5 by 11 inch page, and
C                its boundary is drawn.  1.0 .LE. PLTSIZ
C                .LE. 8.5.
C
C       WX1,WX2,WY1,WY2 = Parameters defining a rectangular
C                         window against which the triangu-
C                         lation is clipped.  (Only the
C                         portion of the triangulation that
C                         lies in the window is drawn.)
C                         (WX1,WY1) and (WX2,WY2) are the
C                         lower left and upper right cor-
C                         ners, respectively.  WX1 < WX2 and
C                         WY1 < WY2.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the nodal coor-
C             dinates.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMSHC.
C
C       TITLE = Type CHARACTER variable or constant contain-
C               ing a string to be centered above the plot.
C               The string must be enclosed in parentheses;
C               i.e., the first and last characters must be
C               '(' and ')', respectively, but these are not
C               displayed.  TITLE may have at most 80 char-
C               acters including the parentheses.
C
C       NUMBR = Option indicator:  If NUMBR = TRUE, the
C               nodal indexes are plotted next to the nodes.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if LUN, PLTSIZ, or N is outside its
C                     valid range.
C             IER = 2 if WX1 >= WX2 or WY1 >= WY2.
C             IER = 3 if an error was encountered in writing
C                     to unit LUN.
C
C   Various plotting options can be controlled by altering
C the data statement below.
C
C Modules required by PLTTR:  None
C
C Intrinsic functions called by PLTTR:  ABS, CHAR, DBLE,
C                                         NINT
C
C***********************************************************
C
      DOUBLE PRECISION DX, DY, FSIZN, FSIZT, R, SFX, SFY, T,
     .                 TX, TY, X0, Y0
      INTEGER IH, IPX1, IPX2, IPY1, IPY2, IW, LP, LPL, N0,
     .        N1
      LOGICAL ANNOT
C
      DATA    ANNOT/.TRUE./,  FSIZN/10.D0/,  FSIZT/16.D0/
C
C Local parameters:
C
C ANNOT =     Logical variable with value TRUE iff the plot
C               is to be annotated with the values of WX1,
C               WX2, WY1, and WY2
C DX =        Window width XMAX-XMIN
C DY =        Window height YMAX-YMIN
C FSIZN =     Font size in points for labeling nodes with
C               their indexes if NUMBR = TRUE
C FSIZT =     Font size in points for the title (and
C               annotation if ANNOT = TRUE)
C IH =        Height of the viewport in points
C IPX1,IPY1 = X and y coordinates (in points) of the lower
C               left corner of the bounding box or viewport
C IPX2,IPY2 = X and y coordinates (in points) of the upper
C               right corner of the bounding box or viewport
C IW =        Width of the viewport in points
C LP =        LIST index (pointer)
C LPL =       Pointer to the last neighbor of N0
C N0 =        Nodal index and DO-loop index
C N1 =        Index of a neighbor of N0
C R =         Aspect ratio DX/DY
C SFX,SFY =   Scale factors for mapping world coordinates
C               (window coordinates in [WX1,WX2] X [WY1,WY2])
C               to viewport coordinates in [IPX1,IPX2] X
C               [IPY1,IPY2]
C T =         Temporary variable
C TX,TY =     Translation vector for mapping world coordi-
C               nates to viewport coordinates
C X0,Y0 =     X(N0),Y(N0) or label location
C
C
C Test for error 1.
C
      IF (LUN .LT. 0  .OR.  LUN .GT. 99  .OR.
     .    PLTSIZ .LT. 1.D0  .OR.  PLTSIZ .GT. 7.5D0  .OR.
     .    N .LT. 3) GO TO 11
C
C Compute the aspect ratio of the window.
C
      DX = WX2 - WX1
      DY = WY2 - WY1
      IF (DX .LE. 0.  .OR.  DY .LE. 0.) GO TO 12
      R = DX/DY
C
C Compute the lower left (IPX1,IPY1) and upper right
C   (IPX2,IPY2) corner coordinates of the bounding box.
C   The coordinates, specified in default user space units
C   (points, at 72 points/inch with origin at the lower
C   left corner of the page), are chosen to preserve the
C   aspect ratio R, and to center the plot on the 8.5 by 11
C   inch page.  The center of the page is (306,396), and
C   T = PLTSIZ/2 in points.
C
      T = 36.D0*PLTSIZ
      IF (R .GE. 1.D0) THEN
        IPX1 = 306 - NINT(T)
        IPX2 = 306 + NINT(T)
        IPY1 = 396 - NINT(T/R)
        IPY2 = 396 + NINT(T/R)
      ELSE
        IPX1 = 306 - NINT(T*R)
        IPX2 = 306 + NINT(T*R)
        IPY1 = 396 - NINT(T)
        IPY2 = 396 + NINT(T)
      ENDIF
C
C Output header comments.
C
      WRITE (LUN,100,ERR=13) IPX1, IPY1, IPX2, IPY2
  100 FORMAT ('%!PS-Adobe-3.0 EPSF-3.0'/
     .        '%%BoundingBox:',4I4/
     .        '%%Title:  Triangulation'/
     .        '%%Creator:  CSRFPACK'/
     .        '%%EndComments')
C
C Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
C   of a viewport obtained by shrinking the bounding box by
C   12% in each dimension.
C
      IW = NINT(0.88D0*DBLE(IPX2-IPX1))
      IH = NINT(0.88D0*DBLE(IPY2-IPY1))
      IPX1 = 306 - IW/2
      IPX2 = 306 + IW/2
      IPY1 = 396 - IH/2
      IPY2 = 396 + IH/2
C
C Set the line thickness to 2 points, and draw the
C   viewport boundary.
C
      T = 2.D0
      WRITE (LUN,110,ERR=13) T
      WRITE (LUN,120,ERR=13) IPX1, IPY1
      WRITE (LUN,130,ERR=13) IPX1, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY2
      WRITE (LUN,130,ERR=13) IPX2, IPY1
      WRITE (LUN,140,ERR=13)
      WRITE (LUN,150,ERR=13)
  110 FORMAT (F12.6,' setlinewidth')
  120 FORMAT (2I4,' moveto')
  130 FORMAT (2I4,' lineto')
  140 FORMAT ('closepath')
  150 FORMAT ('stroke')
C
C Set up a mapping from the window to the viewport.
C
      SFX = DBLE(IW)/DX
      SFY = DBLE(IH)/DY
      TX = IPX1 - SFX*WX1
      TY = IPY1 - SFY*WY1
      WRITE (LUN,160,ERR=13) TX, TY, SFX, SFY
  160 FORMAT (2F12.6,' translate'/
     .        2F12.6,' scale')
C
C Set the line thickness to 1 point.  (Since the scale
C   factors are applied to everything, the width must be
C   specified in world coordinates.)
C
      T = 2.D0/(SFX+SFY)
      WRITE (LUN,110,ERR=13) T
C
C Save the current graphics state, and set the clip path to
C   the boundary of the window.
C
      WRITE (LUN,170,ERR=13)
      WRITE (LUN,180,ERR=13) WX1, WY1
      WRITE (LUN,190,ERR=13) WX2, WY1
      WRITE (LUN,190,ERR=13) WX2, WY2
      WRITE (LUN,190,ERR=13) WX1, WY2
      WRITE (LUN,200,ERR=13)
  170 FORMAT ('gsave')
  180 FORMAT (2F12.6,' moveto')
  190 FORMAT (2F12.6,' lineto')
  200 FORMAT ('closepath clip newpath')
C
C Draw the edges N0->N1, where N1 > N0.  LPL points to the
C   last neighbor of N0.
C
      DO 2 N0 = 1,N
        X0 = X(N0)
        Y0 = Y(N0)
	LPL = LEND(N0)
	LP = LPL
C
C   Loop on neighbors N1 of N0.
C
    1   LP = LPTR(LP)
          N1 = ABS(LIST(LP))
	  IF (N1 .GT. N0) THEN
C
C   Add the edge to the path.
C
            WRITE (LUN,210,ERR=13) X0, Y0, X(N1), Y(N1)
  210       FORMAT (2F12.6,' moveto',2F12.6,' lineto')
          ENDIF
          IF (LP .NE. LPL) GO TO 1
    2   CONTINUE
C
C Paint the path and restore the saved graphics state (with
C   no clip path).
C
      WRITE (LUN,150,ERR=13)
      WRITE (LUN,230,ERR=13)
  230 FORMAT ('grestore')
      IF (NUMBR) THEN
C
C Nodes in the window are to be labeled with their indexes.
C   Convert FSIZN from points to world coordinates, and
C   output the commands to select a font and scale it.
C
        T = FSIZN*2.D0/(SFX+SFY)
        WRITE (LUN,240,ERR=13) T
  240   FORMAT ('/Helvetica findfont'/
     .          F12.6,' scalefont setfont')
C
C   Loop on nodes N0 with coordinates (X0,Y0).
C
        DO 3 N0 = 1,N
          X0 = X(N0)
          Y0 = Y(N0)
          IF (X0 .LT. WX1  .OR.  X0 .GT. WX2  .OR.
     .        Y0 .LT. WY1  .OR.  Y0 .GT. WY2) GO TO 3
C
C   Move to (X0,Y0), and draw the label N0.  The first char-
C     acter will have its lower left corner about one
C     character width to the right of the nodal position.
C
          WRITE (LUN,180,ERR=13) X0, Y0
          WRITE (LUN,250,ERR=13) N0
  250     FORMAT ('(',I3,') show')
    3     CONTINUE
      ENDIF
C
C Convert FSIZT from points to world coordinates, and output
C   the commands to select a font and scale it.
C
      T = FSIZT*2.D0/(SFX+SFY)
      WRITE (LUN,240,ERR=13) T
C
C Display TITLE centered above the plot:
C
      Y0 = WY2 + 3.D0*T
      WRITE (LUN,260,ERR=13) TITLE, (WX1+WX2)/2.D0, Y0
  260 FORMAT (A80/'  stringwidth pop 2 div neg ',F12.6,
     .        ' add ',F12.6,' moveto')
      WRITE (LUN,270,ERR=13) TITLE
  270 FORMAT (A80/'  show')
      IF (ANNOT) THEN
C
C Display the window extrema below the plot.
C
        X0 = WX1
        Y0 = WY1 - 100.D0/(SFX+SFY)
        WRITE (LUN,180,ERR=13) X0, Y0
        WRITE (LUN,280,ERR=13) WX1, WX2
        Y0 = Y0 - 2.D0*T
        WRITE (LUN,290,ERR=13) X0, Y0, WY1, WY2
  280   FORMAT ('(Window:   WX1 = ',E9.3,',   WX2 = ',E9.3,
     .          ') show')
  290   FORMAT ('(Window:  ) stringwidth pop ',F12.6,' add',
     .          F12.6,' moveto'/
     .          '( WY1 = ',E9.3,',   WY2 = ',E9.3,') show')
      ENDIF
C
C Paint the path and output the showpage command and
C   end-of-file indicator.
C
      WRITE (LUN,300,ERR=13)
  300 FORMAT ('stroke'/
     .        'showpage'/
     .        '%%EOF')
C
C HP's interpreters require a one-byte End-of-PostScript-Job
C   indicator (to eliminate a timeout error message):
C   ASCII 4.
C
      WRITE (LUN,310,ERR=13) CHAR(4)
  310 FORMAT (A1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
   11 IER = 1
      RETURN
C
C DX or DY is not positive.
C
   12 IER = 2
      RETURN
C
C Error writing to unit LUN.
C
   13 IER = 3
      RETURN
      END
      SUBROUTINE SWAP (IN1,IN2,IO1,IO2, LIST,LPTR,
     .                 LEND, LP21)
      INTEGER IN1, IN2, IO1, IO2, LIST(*), LPTR(*), LEND(*),
     .        LP21
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/22/98
C
C   Given a triangulation of a set of points on the unit
C sphere, this subroutine replaces a diagonal arc in a
C strictly convex quadrilateral (defined by a pair of adja-
C cent triangles) with the other diagonal.  Equivalently, a
C pair of adjacent triangles is replaced by another pair
C having the same union.
C
C
C On input:
C
C       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
C                         the quadrilateral.  IO1-IO2 is re-
C                         placed by IN1-IN2.  (IO1,IO2,IN1)
C                         and (IO2,IO1,IN2) must be trian-
C                         gles on input.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C On output:
C
C       LIST,LPTR,LEND = Data structure updated with the
C                        swap -- triangles (IO1,IO2,IN1) and
C                        (IO2,IO1,IN2) are replaced by
C                        (IN1,IN2,IO2) and (IN2,IN1,IO1)
C                        unless LP21 = 0.
C
C       LP21 = Index of IN1 as a neighbor of IN2 after the
C              swap is performed unless IN1 and IN2 are
C              adjacent on input, in which case LP21 = 0.
C
C Module required by SWAP:  LSTPTR
C
C Intrinsic function called by SWAP:  ABS
C
C***********************************************************
C
      INTEGER LSTPTR
      INTEGER LP, LPH, LPSAV
C
C Local parameters:
C
C LP,LPH,LPSAV = LIST pointers
C
C
C Test for IN1 and IN2 adjacent.
C
      LP = LSTPTR(LEND(IN1),IN2,LIST,LPTR)
      IF (ABS(LIST(LP)) .EQ. IN2) THEN
        LP21 = 0
        RETURN
      ENDIF
C
C Delete IO2 as a neighbor of IO1.
C
      LP = LSTPTR(LEND(IO1),IN2,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO2 is the last neighbor of IO1, make IN2 the
C   last neighbor.
C
      IF (LEND(IO1) .EQ. LPH) LEND(IO1) = LP
C
C Insert IN2 as a neighbor of IN1 following IO1
C   using the hole created above.
C
      LP = LSTPTR(LEND(IN1),IO1,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN2
      LPTR(LPH) = LPSAV
C
C Delete IO1 as a neighbor of IO2.
C
      LP = LSTPTR(LEND(IO2),IN1,LIST,LPTR)
      LPH = LPTR(LP)
      LPTR(LP) = LPTR(LPH)
C
C If IO1 is the last neighbor of IO2, make IN1 the
C   last neighbor.
C
      IF (LEND(IO2) .EQ. LPH) LEND(IO2) = LP
C
C Insert IN1 as a neighbor of IN2 following IO2.
C
      LP = LSTPTR(LEND(IN2),IO2,LIST,LPTR)
      LPSAV = LPTR(LP)
      LPTR(LP) = LPH
      LIST(LPH) = IN1
      LPTR(LPH) = LPSAV
      LP21 = LPH
      RETURN
      END
      LOGICAL FUNCTION SWPTC (IN1,IN2,IO1,IO2,X,Y,Z, IER)
      INTEGER IN1, IN2, IO1, IO2, IER
      DOUBLE PRECISION X(*), Y(*), Z(*)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/12/96
C
C   This function applies a swap test to a pair of adjacent
C surface triangles.  The diagonal edge (shared triangle
C side in the quadrilateral obtained by projecting onto the
C x-y plane) should be swapped for the other diagonal if and
C only if the the pair of triangles is nonconvex (the fourth
C vertex lies below the plane defined by one of the tri-
C angles).
C
C
C On input:
C
C       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
C                         the quadrilateral.  IO1-IO2 is the
C                         triangulation edge (shared trian-
C                         gle side) to be replaced by
C                         IN1-IN2 if the decision is to
C                         swap.  The triples (IO1,IO2,IN1)
C                         and (IO2,IO1,IN2) must define tri-
C                         angles (be in counterclockwise
C                         order) on input.
C
C       X,Y,Z = Arrays containing the nodal coordinates and
C               data values (coordinates of the surface
C               points).
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       SWPTC = TRUE if and only if triangles (IO1,IO2,IN1)
C               and (IO2,IO1,IN2) should be replaced by
C               (IN1,IN2,IO2) and (IN2,IN1,IO1).
C
C       IER = Error indicator:
C             IER = 0 if no error was encountered.
C             IER = 1 if the pair of triangles is nonconvex
C                     but a swap cannot be applied because
C                     the quadrilateral is not convex.  A
C                     convex triangulation does not exists
C                     in this case.
C
C Modules required by SWPTC:  None
C
C***********************************************************
C
      DOUBLE PRECISION DX1, DX2, DX3, DY1, DY2, DY3, DZ1,
     .                 DZ2, DZ3
C
      DX1 = X(IO1) - X(IN1)
      DX2 = X(IN2) - X(IN1)
      DX3 = X(IO2) - X(IN1)
C
      DY1 = Y(IO1) - Y(IN1)
      DY2 = Y(IN2) - Y(IN1)
      DY3 = Y(IO2) - Y(IN1)
C
      DZ1 = Z(IO1) - Z(IN1)
      DZ2 = Z(IN2) - Z(IN1)
      DZ3 = Z(IO2) - Z(IN1)
C
C Test for a convex surface.
C
      IF (DX1*(DY2*DZ3-DY3*DZ2) - DY1*(DX2*DZ3-DX3*DZ2) +
     .    DZ1*(DX2*DY3-DX3*DY2) .LE. 0.) THEN
        SWPTC = .FALSE.
        IER = 0
      ELSE
C
C Test for a nonconvex quadrilateral.
C
        IF (DX1*DY2-DX2*DY1 .LT. 0.  .OR.
     .      DX2*DY3-DX3*DY2 .LT. 0.) THEN
C
C Nonconvex surface and nonconvex quadrilateral.
C
          SWPTC = .FALSE.
          IER = 1
        ELSE
          SWPTC = .TRUE.
          IER = 0
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE TRFIND (NST,PX,PY,N,X,Y,LIST,LPTR,LEND, I1,
     .                   I2,I3)
      INTEGER NST, N, LIST(*), LPTR(*), LEND(N), I1, I2, I3
      DOUBLE PRECISION PX, PY, X(N), Y(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/28/98
C
C   This subroutine locates a point P relative to a triangu-
C lation created by Subroutine TRMESH or TRMSHR.  If P is
C contained in a triangle, the three vertex indexes are
C returned.  Otherwise, the indexes of the rightmost and
C leftmost visible boundary nodes are returned.
C
C
C On input:
C
C       NST = Index of a node at which TRFIND begins the
C             search.  Search time depends on the proximity
C             of this node to P.
C
C       PX,PY = X and y coordinates of the point P to be
C               located.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the nodes in the triangulation.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMESH.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       I1,I2,I3 = Nodal indexes, in counterclockwise order,
C                  of the vertices of a triangle containing
C                  P if P is contained in a triangle.  If P
C                  is not in the convex hull of the nodes,
C                  I1 indexes the rightmost visible boundary
C                  node, I2 indexes the leftmost visible
C                  boundary node, and I3 = 0.  Rightmost and
C                  leftmost are defined from the perspective
C                  of P, and a pair of points are visible
C                  from each other if and only if the line
C                  segment joining them intersects no trian-
C                  gulation arc.  If P and all of the nodes
C                  lie on a common line, then I1 = I2 = I3 =
C                  0 on output.
C
C Modules required by TRFIND:  DSTORE, JRAND, LEFT, LSTPTR
C
C Intrinsic function called by TRFIND:  ABS
C
C***********************************************************
C
      DOUBLE PRECISION DSTORE
      INTEGER JRAND, LSTPTR
      LOGICAL LEFT
      DOUBLE PRECISION B1, B2, XA, XB, XC, XP, YA, YB, YC,
     .                 YP
      INTEGER IX, IY, IZ, LP, N0, N1, N1S, N2, N2S, N3, N4,
     .        NB, NF, NL, NP, NPP
      LOGICAL FRWRD
C
      SAVE    IX, IY, IZ
      DATA    IX/1/, IY/2/, IZ/3/
C
C Local parameters:
C
C B1,B2 =    Unnormalized barycentric coordinates of P with
C              respect to (N1,N2,N3)
C IX,IY,IZ = Integer seeds for JRAND
C LP =       LIST pointer
C N0,N1,N2 = Nodes in counterclockwise order defining a
C              cone (with vertex N0) containing P
C N1S,N2S =  Saved values of N1 and N2
C N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
C NB =       Index of a boundary node -- first neighbor of
C              NF or last neighbor of NL in the boundary
C              traversal loops
C NF,NL =    First and last neighbors of N0, or first
C              (rightmost) and last (leftmost) nodes
C              visible from P when P is exterior to the
C              triangulation
C NP,NPP =   Indexes of boundary nodes used in the boundary
C              traversal loops
C XA,XB,XC = Dummy arguments for FRWRD
C YA,YB,YC = Dummy arguments for FRWRD
C XP,YP =    Local variables containing the components of P
C
C Statement function:
C
C FRWRD = TRUE iff C is forward of A->B
C              iff <A->B,A->C> .GE. 0.
C
      FRWRD(XA,YA,XB,YB,XC,YC) = (XB-XA)*(XC-XA) +
     .                           (YB-YA)*(YC-YA) .GE. 0.
C
C Initialize variables.
C
      XP = PX
      YP = PY
      N0 = NST
      IF (N0 .LT. 1  .OR.  N0 .GT. N)
     .  N0 = JRAND(N, IX,IY,IZ )
C
C Set NF and NL to the first and last neighbors of N0, and
C   initialize N1 = NF.
C
    1 LP = LEND(N0)
      NL = LIST(LP)
      LP = LPTR(LP)
      NF = LIST(LP)
      N1 = NF
C
C Find a pair of adjacent neighbors N1,N2 of N0 that define
C   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
C
      IF (NL .GT. 0) GO TO 2
C
C   N0 is a boundary node.  Test for P exterior.
C
      NL = -NL
      IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) ) THEN
        NL = N0
        GO TO 9
      ENDIF
      IF ( .NOT. LEFT(X(NL),Y(NL),X(N0),Y(N0),XP,YP) ) THEN
        NB = NF
        NF = N0
        NP = NL
        NPP = N0
        GO TO 11
      ENDIF
      GO TO 3
C
C   N0 is an interior node.  Find N1.
C
    2 IF ( LEFT(X(N0),Y(N0),X(N1),Y(N1),XP,YP) ) GO TO 3
        LP = LPTR(LP)
        N1 = LIST(LP)
        IF (N1 .EQ. NL) GO TO 6
        GO TO 2
C
C   P is to the left of edge N0->N1.  Initialize N2 to the
C     next neighbor of N0.
C
    3 LP = LPTR(LP)
        N2 = ABS(LIST(LP))
        IF ( .NOT. LEFT(X(N0),Y(N0),X(N2),Y(N2),XP,YP) )
     .    GO TO 7
        N1 = N2
        IF (N1 .NE. NL) GO TO 3
      IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) )
     .  GO TO 6
      IF (XP .EQ. X(N0) .AND. YP .EQ. Y(N0)) GO TO 5
C
C   P is left of or on edges N0->NB for all neighbors NB
C     of N0.
C   All points are collinear iff P is left of NB->N0 for
C     all neighbors NB of N0.  Search the neighbors of N0.
C     NOTE -- N1 = NL and LP points to NL.
C
    4 IF ( .NOT. LEFT(X(N1),Y(N1),X(N0),Y(N0),XP,YP) )
     .  GO TO 5
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
        IF (N1 .EQ. NL) GO TO 17
        GO TO 4
C
C   P is to the right of N1->N0, or P=N0.  Set N0 to N1 and
C     start over.
C
    5 N0 = N1
      GO TO 1
C
C   P is between edges N0->N1 and N0->NF.
C
    6 N2 = NF
C
C P is contained in the wedge defined by line segments
C   N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
C   N3 to the node opposite N1->N2, and save N1 and N2 to
C   test for cycling.
C
    7 N3 = N0
      N1S = N1
      N2S = N2
C
C Top of edge hopping loop.  Test for termination.
C
    8 IF ( LEFT(X(N1),Y(N1),X(N2),Y(N2),XP,YP) ) THEN
C
C   P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
C     error resulted from floating point inaccuracy and
C     collinearity.  Compute the unnormalized barycentric
C     coordinates of P with respect to (N1,N2,N3).
C
        B1 = (X(N3)-X(N2))*(YP-Y(N2)) -
     .       (XP-X(N2))*(Y(N3)-Y(N2))
        B2 = (X(N1)-X(N3))*(YP-Y(N3)) -
     .       (XP-X(N3))*(Y(N1)-Y(N3))
        IF (DSTORE(B1+1.D0) .GE. 1.D0  .AND.
     .      DSTORE(B2+1.D0) .GE. 1.D0) GO TO 16
C
C   Restart with N0 randomly selected.
C
        N0 = JRAND(N, IX,IY,IZ )
        GO TO 1
      ENDIF
C
C   Set N4 to the neighbor of N2 which follows N1 (node
C     opposite N2->N1) unless N1->N2 is a boundary edge.
C
      LP = LSTPTR(LEND(N2),N1,LIST,LPTR)
      IF (LIST(LP) .LT. 0) THEN
        NF = N2
        NL = N1
        GO TO 9
      ENDIF
      LP = LPTR(LP)
      N4 = ABS(LIST(LP))
C
C   Select the new edge N1->N2 which intersects the line
C     segment N0-P, and set N3 to the node opposite N1->N2.
C
      IF ( LEFT(X(N0),Y(N0),X(N4),Y(N4),XP,YP) ) THEN
        N3 = N1
        N1 = N4
        N2S = N2
        IF (N1 .NE. N1S  .AND.  N1 .NE. N0) GO TO 8
      ELSE
        N3 = N2
        N2 = N4
        N1S = N1
        IF (N2 .NE. N2S  .AND.  N2 .NE. N0) GO TO 8
      ENDIF
C
C   The starting node N0 or edge N1-N2 was encountered
C     again, implying a cycle (infinite loop).  Restart
C     with N0 randomly selected.
C
      N0 = JRAND(N, IX,IY,IZ )
      GO TO 1
C
C Boundary traversal loops.  NL->NF is a boundary edge and
C   P RIGHT NL->NF.  Save NL and NF.

    9 NP = NL
      NPP = NF
C
C Find the first (rightmost) visible boundary node NF.  NB
C   is set to the first neighbor of NF, and NP is the last
C   neighbor.
C
   10 LP = LEND(NF)
      LP = LPTR(LP)
      NB = LIST(LP)
      IF ( .NOT. LEFT(X(NF),Y(NF),X(NB),Y(NB),XP,YP) )
     .  GO TO 12
C
C   P LEFT NF->NB and thus NB is not visible unless an error
C     resulted from floating point inaccuracy and collinear-
C     ity of the 4 points NP, NF, NB, and P.
C
   11 IF ( FRWRD(X(NF),Y(NF),X(NP),Y(NP),XP,YP)  .OR.
     .     FRWRD(X(NF),Y(NF),X(NP),Y(NP),X(NB),Y(NB)) ) THEN
        I1 = NF
        GO TO 13
      ENDIF
C
C   Bottom of loop.
C
   12 NP = NF
      NF = NB
      GO TO 10
C
C Find the last (leftmost) visible boundary node NL.  NB
C   is set to the last neighbor of NL, and NPP is the first
C   neighbor.
C
   13 LP = LEND(NL)
      NB = -LIST(LP)
      IF ( .NOT. LEFT(X(NB),Y(NB),X(NL),Y(NL),XP,YP) )
     .  GO TO 14
C
C   P LEFT NB->NL and thus NB is not visible unless an error
C     resulted from floating point inaccuracy and collinear-
C     ity of the 4 points P, NB, NL, and NPP.
C
      IF ( FRWRD(X(NL),Y(NL),X(NPP),Y(NPP),XP,YP)  .OR.
     .     FRWRD(X(NL),Y(NL),X(NPP),Y(NPP),X(NB),Y(NB)) )
     .  GO TO 15
C
C   Bottom of loop.
C
   14 NPP = NL
      NL = NB
      GO TO 13
C
C NL is the leftmost visible boundary node.
C
   15 I2 = NL
      I3 = 0
      RETURN
C
C P is in the triangle (N1,N2,N3).
C
   16 I1 = N1
      I2 = N2
      I3 = N3
      RETURN
C
C All points are collinear.
C
   17 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
      END
      SUBROUTINE TRMSHC (N,X,Y,Z, LIST,LPTR,LEND,LNEW,NEAR,
     .                   NEXT,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), LNEW, NEAR(N),
     .        NEXT(N), IER
      DOUBLE PRECISION X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/27/98
C
C   This subroutine creates a convexity-preserving triangu-
C lation, if it exists, of a set of N arbitrarily
C distributed points in the plane (referred to as nodes)
C with associated data values.  The triangulation is defined
C as a set of triangles with the following four properties:
C
C  1)  The triangle vertices are nodes.
C  2)  No triangle contains a node other than its vertices.
C  3)  The interiors of the triangles are pairwise disjoint.
C  4)  The union of triangles is the convex hull of the set
C        of nodes (the smallest convex set which contains
C        the nodes).
C
C The triangulation is said to be convex if the correspond-
C ing triangulated surface (the triangle-based piecewise
C linear interpolant of the data values) is convex.
C
C   Since the algorithm proceeds by adding nodes incremen-
C tally, the triangulation may be updated with the addition
C or deletion of a node very efficiently.  The adjacency
C information representing the triangulation is stored as a
C linked list requiring approximately 13N storage locations.
C
C The algorithm has expected time complexity O(N*log(N)).
C
C
C   The following is a list of related subprograms which a
C user may wish to call directly:
C
C  ADDNDC - Updates the triangulation by adding a new node.
C
C  AREAP  - Computes the area bounded by a closed polygonal
C             curve such as the boundary of the triangula-
C             tion.
C
C  BNODES - Returns an array containing the indexes of the
C             boundary nodes in counterclockwise order.
C             Counts of boundary nodes, triangles, and edges
C             are also returned.
C
C  CIRCUM - Computes the area, circumcenter, circumradius,
C             and, optionally, the aspect ratio of a trian-
C             gle defined by user-specified vertices.
C
C  DELBE  - Deletes extraneous boundary edges (associated
C             with nearly null triangles) from the triangu-
C             lation.
C
C  DELNDC - Updates the triangulation with the deletion of a
C             node.
C
C  JRAND  - Generates a uniformly distributed pseudo-random
C             integer.
C
C  LEFT   - Locates a point relative to a line.
C
C  TRPRNT - Prints the triangulation data structure and,
C             optionally, the nodal coordinates and/or
C             data values.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of the nodes and data values.
C               (X(K),Y(K)) is referred to as node K, and K
C               is referred to as a nodal index.  The first
C               three nodes must not be collinear.
C
C The above parameters are not altered by this routine.
C
C       LIST,LPTR = Arrays of length at least 6N-12.
C
C       LEND = Array of length at least N.
C
C       NEAR,NEXT = Work space arrays of length at least N.
C                   The space is used to efficiently find
C                   a triangle containing a node.
C
C On output:
C
C       LIST = Set of nodal indexes which, along with LPTR,
C              LEND, and LNEW, define the triangulation as a
C              set of N adjacency lists -- counterclockwise-
C              ordered sequences of neighboring nodes such
C              that the first and last neighbors of a bound-
C              ary node are boundary nodes (the first neigh-
C              bor of an interior node is arbitrary).  In
C              order to distinguish between interior and
C              boundary nodes, the last neighbor of each
C              boundary node is represented by the negative
C              of its index.
C
C       LPTR = Set of pointers (LIST indexes) in one-to-one
C              correspondence with the elements of LIST.
C              LIST(LPTR(I)) indexes the node which follows
C              LIST(I) in cyclical counterclockwise order
C              (the first neighbor follows the last neigh-
C              bor).
C
C       LEND = Set of pointers to adjacency lists.  LEND(K)
C              points to the last neighbor of node K for
C              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and
C              only if K is a boundary node.
C
C       LNEW = Pointer to the first empty location in LIST
C              and LPTR (list length plus one).  LIST, LPTR,
C              LEND, and LNEW are not altered if IER = -1 or
C              IER = -2, and are incomplete if IER .NE. 0.
C
C       NEAR,NEXT = Garbage.
C
C       IER = Error indicator:
C             IER =  0 if no errors were encountered.
C             IER = -1 if N < 3 on input.
C             IER = -2 if the first three nodes are
C                      collinear.
C             IER = -3 if a convex triangulation does not
C                      exist (and Subroutine ADDNDC returned
C                      IER = -3).
C             IER = -4 if a convex triangulation does not
C                      exist (and Subroutine ADDNDC returned
C                      IER = -4).
C             IER =  L if nodes L and M coincide for some
C                      M > L.  The linked list represents
C                      a triangulation of nodes 1 to M-1
C                      in this case.
C
C Modules required by TRMSHC:  ADDNDC, BDYADD, DSTORE,
C                                INSERT, INTADD, JRAND,
C                                LEFT, LSTPTR, SWAP, SWPTC,
C                                TRFIND
C
C Intrinsic function called by TRMSHC:  ABS
C
C***********************************************************
C
      DOUBLE PRECISION D1, D2, D3
      INTEGER I, I0, J, JP, JS, K, LP, LPL, LPS, NEXTI, NN
      LOGICAL LEFT
C
C Local parameters:
C
C D1,D2,D3 = Squared distances from node K to nodes 1, 2,
C              and 3, respectively
C I,J =      Nodal indexes
C JP,JS =    Neighbors of node K such that (JP,J,JS) is a
C              sequence of adjacent neighbors of K
C I0 =       Index of the node preceding I in a sequence of
C              unprocessed nodes:  I = NEXT(I0)
C K =        Index of node to be added and DO-loop index:
C              K > 3
C LP =       LIST index (pointer) of a neighbor of K
C LPL =      Pointer to the last neighbor of K
C LPS =      Pointer to JS as a neighbor of K
C NEXTI =    NEXT(I)
C NN =       Local copy of N
C
      NN = N
      IF (NN .LT. 3) THEN
        IER = -1
        RETURN
      ENDIF
C
C Store the first triangle in the linked list.
C
      IF ( .NOT. LEFT(X(1),Y(1),X(2),Y(2),X(3),Y(3)) ) THEN
C
C   The initial triangle is (1,3,2).
C
        LIST(1) = 3
        LPTR(1) = 2
        LIST(2) = -2
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 1
        LPTR(3) = 4
        LIST(4) = -3
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 2
        LPTR(5) = 6
        LIST(6) = -1
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSEIF ( .NOT. LEFT(X(2),Y(2),X(1),Y(1),X(3),Y(3)) )
     .       THEN
C
C   The initial triangle is (1,2,3).
C
        LIST(1) = 2
        LPTR(1) = 2
        LIST(2) = -3
        LPTR(2) = 1
        LEND(1) = 2
C
        LIST(3) = 3
        LPTR(3) = 4
        LIST(4) = -1
        LPTR(4) = 3
        LEND(2) = 4
C
        LIST(5) = 1
        LPTR(5) = 6
        LIST(6) = -2
        LPTR(6) = 5
        LEND(3) = 6
C
      ELSE
C
C   The first three nodes are collinear.
C
        IER = -2
        RETURN
      ENDIF
C
C Initialize LNEW and test for N = 3.
C
      LNEW = 7
      IF (NN .EQ. 3) THEN
        IER = 0
        RETURN
      ENDIF
C
C A near-node data structure (NEAR and NEXT) is used to
C   obtain an expected-time (N*log(N)) incremental algorithm
C   by enabling constant search time for locating each new
C   node in the triangulation.
C
C For each unprocessed node I, NEAR(I) is the index of a
C   vertex of a triangle, if any, that contains I, or a
C   boundary node that is visible from I otherwise.
C   (NEAR(I) is used as the starting point for the search
C   in Subroutine TRFIND.)  Thus, the unprocessed nodes are
C   partitioned into subsets, each associated with a
C   triangulation node.
C
C Since it is necessary to efficiently find the subset of
C   unprocessed nodes associated with each triangulation
C   node J (those that have J as their NEAR entries), the
C   subsets are stored in NEAR and NEXT as follows:  for
C   each node J in the triangulation, I = NEAR(J) is the
C   first unprocessed node in J's set (with I = 0 if the
C   set is empty), L = NEXT(I) (if I > 0) is the second,
C   NEXT(L) (if L > 0) is the third, etc.  The nodes in each
C   set are initially ordered by increasing indexes (which
C   maximizes efficiency) but that ordering is not main-
C   tained as the data structure is updated.
C
C Initialize the data structure for the single triangle.
C   For each unprocessed node K, NEAR(K) is set to the near-
C   est vertex.  This is not necessarily visible but it's
C   close enough.
C
      NEAR(1) = 0
      NEAR(2) = 0
      NEAR(3) = 0
      DO 1 K = NN,4,-1
        D1 = (X(K)-X(1))**2 + (Y(K)-Y(1))**2
        D2 = (X(K)-X(2))**2 + (Y(K)-Y(2))**2
        D3 = (X(K)-X(3))**2 + (Y(K)-Y(3))**2
        IF (D1 .LE. D2  .AND.  D1 .LE. D3) THEN
          NEAR(K) = 1
          NEXT(K) = NEAR(1)
          NEAR(1) = K
        ELSEIF (D2 .LE. D1  .AND.  D2 .LE. D3) THEN
          NEAR(K) = 2
          NEXT(K) = NEAR(2)
          NEAR(2) = K
        ELSE
          NEAR(K) = 3
          NEXT(K) = NEAR(3)
          NEAR(3) = K
        ENDIF
    1   CONTINUE
C
C Add the remaining nodes
C
      DO 7 K = 4,NN
        CALL ADDNDC (NEAR(K),K,X,Y,Z, LIST,LPTR,LEND,
     .               LNEW, IER)
        IF (IER .NE. 0) RETURN
C
C Remove K from the set of unprocessed nodes associated
C   with NEAR(K).
C
        I = NEAR(K)
        IF (NEAR(I) .EQ. K) THEN
          NEAR(I) = NEXT(K)
        ELSE
          I = NEAR(I)
    2     I0 = I
            I = NEXT(I0)
            IF (I .NE. K) GO TO 2
          NEXT(I0) = NEXT(K)
        ENDIF
        NEAR(K) = 0
C
C Loop on neighbors J of node K.  JP and JS are the
C   predecessor and successor of J, respectively.
C
        LPL = LEND(K)
        JP = ABS(LIST(LPL))
        LP = LPL
    3   LP = LPTR(LP)
          J = ABS(LIST(LP))
          LPS = LPTR(LP)
          JS = ABS(LIST(LPS))
C
C Loop on elements I in the sequence of unprocessed nodes
C   associated with J:  K is a candidate for replacing J
C   as the value of NEAR(I).  The next value of I in the
C   sequence, NEXT(I), must be saved before I is moved
C   because it is altered by adding I to K's set.
C
          I = NEAR(J)
    4     IF (I .EQ. 0) GO TO 6
          NEXTI = NEXT(I)
C
C Test for J no longer visible from I:  J separated from I
C   by edge K-JP or K-JS.
C
          IF ( LEFT(X(K),Y(K),X(J),Y(J),X(I),Y(I)) ) THEN
            IF ( LEFT(X(JS),Y(JS),X(K),Y(K),X(I),Y(I)) )
     .        GO TO 5
            IF ( LEFT(X(JS),Y(JS),X(J),Y(J),X(I),Y(I)) )
     .        GO TO 5
          ELSE
            IF ( LEFT(X(K),Y(K),X(JP),Y(JP),X(I),Y(I)) )
     .        GO TO 5
            IF ( LEFT(X(J),Y(J),X(JP),Y(JP),X(I),Y(I)) )
     .        GO TO 5
          ENDIF
C
C Replace J by K as the NEAR entry for I:  update NEAR(I),
C   and remove I from J's set of unprocessed nodes and add
C   it to K's set.
C
          NEAR(I) = K
          IF (I .EQ. NEAR(J)) THEN
            NEAR(J) = NEXTI
          ELSE
            NEXT(I0) = NEXTI
          ENDIF
          NEXT(I) = NEAR(K)
          NEAR(K) = I
          I = I0
C
C Bottom of loop on I.
C
    5     I0 = I
          I = NEXTI
          GO TO 4
C
C Bottom of loop on neighbors J.
C
    6     IF (LP .NE. LPL) THEN
            JP = J
            GO TO 3
          ENDIF
    7   CONTINUE
      RETURN
      END
      SUBROUTINE TRPRNT (N,X,Y,Z,IFLAG,LIST,LPTR,LEND,LOUT)
      INTEGER N, IFLAG, LIST(*), LPTR(*), LEND(N), LOUT
      DOUBLE PRECISION X(N), Y(N), Z(N)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/25/98
C
C   This subroutine prints triangulation adjacency lists
C and, optionally, nodal coordinates and data values on log-
C ical unit LOUT.  The list of neighbors of a boundary node
C is followed by index 0.  The numbers of boundary nodes,
C triangles, and edges are also printed.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3
C           and N .LE. 9999.
C
C       X,Y,Z = Arrays of length N containing the nodal
C               coordinates and data values if IFLAG = 0,
C               or (X and Y only) arrays of length N con-
C               taining the nodal coordinates if IFLAG > 0,
C               or unused dummy parameters if IFLAG < 0.
C
C       IFLAG = Coordinate option indicator:
C               IFLAG = 0 if X, Y, and Z are to be printed
C                         (to 6 decimal places).
C               IFLAG > 0 if only X and Y are to be printed
C                         (to 6 decimal places).
C               IFLAG < 0 if only the adjacency lists are to
C                         be printed.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.
C
C       LOUT = Logical unit for output.  If LOUT is not in
C              the range 0 to 99, output is written to
C              logical unit 6.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C   The adjacency lists and coordinates (as specified by
C IFLAG) are written to unit LOUT.
C
C Modules required by TRPRNT:  None
C
C***********************************************************
C
      INTEGER I, INC, K, LP, LPL, LUN, NABOR(100), NB, ND,
     .        NE, NL, NLMAX, NMAX, NODE, NN, NT
      DATA  NMAX/9999/,  NLMAX/58/
C
C Local parameters:
C
C I =     NABOR index (1 to K)
C INC =   Increment for NL associated with an adjacency list
C K =     Counter and number of neighbors of NODE
C LP =    LIST pointer of a neighbor of NODE
C LPL =   Pointer to the last neighbor of NODE
C LUN =   Logical unit for output (copy of LOUT)
C NABOR = Array containing the adjacency list associated
C           with NODE, with zero appended if NODE is a
C           boundary node
C NB =    Number of boundary nodes encountered
C ND =    Index of a neighbor of NODE (or negative index)
C NE =    Number of edges in the triangulation
C NL =    Number of lines that have been printed on the
C           current page
C NLMAX = Maximum number of print lines per page (except
C           for the last page which may have two addi-
C           tional lines)
C NMAX =  Upper bound on N (allows 4-digit indexes)
C NODE =  Index of a node and DO-loop index (1 to N)
C NN =    Local copy of N
C NT =    Number of triangles in the triangulation
C
      NN = N
      LUN = LOUT
      IF (LUN .LT. 0  .OR.  LUN .GT. 99) LUN = 6
C
C Print a heading and test the range of N.
C
      WRITE (LUN,100) NN
      IF (NN .LT. 3  .OR.  NN .GT. NMAX) THEN
C
C N is outside its valid range.
C
        WRITE (LUN,110)
        RETURN
      ENDIF
C
C Initialize NL (the number of lines printed on the current
C   page) and NB (the number of boundary nodes encountered).
C
      NL = 6
      NB = 0
      IF (IFLAG .LT. 0) THEN
C
C Print LIST only.  K is the number of neighbors of NODE
C   that have been stored in NABOR.
C
        WRITE (LUN,101)
        DO 2 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    1     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 1
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.  Correct the sign of the last
C     neighbor, add 0 to the end of the list, and increment
C     NB.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print the list of neighbors.
C
          INC = (K-1)/14 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,104) NODE, (NABOR(I), I = 1,K)
          IF (K .NE. 14) WRITE (LUN,107)
    2     CONTINUE
      ELSEIF (IFLAG .GT. 0) THEN
C
C Print X, Y, and LIST.
C
        WRITE (LUN,102)
        DO 4 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    3     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 3
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print X, Y, and NABOR.
C
          INC = (K-1)/8 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,105) NODE, X(NODE), Y(NODE),
     .                    (NABOR(I), I = 1,K)
          IF (K .NE. 8) WRITE (LUN,107)
    4     CONTINUE
      ELSE
C
C Print X, Y, Z, and LIST.
C
        WRITE (LUN,103)
        DO 6 NODE = 1,NN
          LPL = LEND(NODE)
          LP = LPL
          K = 0
C
    5     K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 5
          IF (ND .LE. 0) THEN
C
C   NODE is a boundary node.
C
            NABOR(K) = -ND
            K = K + 1
            NABOR(K) = 0
            NB = NB + 1
          ENDIF
C
C   Increment NL and print X, Y, Z, and NABOR.
C
          INC = (K-1)/5 + 2
          NL = NL + INC
          IF (NL .GT. NLMAX) THEN
            WRITE (LUN,108)
            NL = INC
          ENDIF
          WRITE (LUN,106) NODE, X(NODE), Y(NODE),
     .                    Z(NODE), (NABOR(I), I = 1,K)
          IF (K .NE. 5) WRITE (LUN,107)
    6     CONTINUE
      ENDIF
C
C Print NB, NE, and NT (boundary nodes, edges, and
C   triangles).
C
      IF (NB .NE. 0) THEN
        NE = 3*NN - NB - 3
        NT = 2*NN - NB - 2
      ELSE
        NE = 3*NN - 6
        NT = 2*NN - 4
      ENDIF
      WRITE (LUN,109) NB, NE, NT
      RETURN
C
C Print formats:
C
  100 FORMAT (///15X,'CSRFPACK Triangulation Data ',
     .        'Structure,  N = ',I5//)
  101 FORMAT (1X,'Node',31X,'Neighbors of Node'//)
  102 FORMAT (1X,'Node',5X,'X(Node)',8X,'Y(Node)',
     .        18X,'Neighbors of Node'//)
  103 FORMAT (1X,'Node',5X,'X(Node)',8X,'Y(Node)',8X,
     .        'Z(Node)',11X,'Neighbors of Node'//)
  104 FORMAT (1X,I4,4X,14I5/(1X,8X,14I5))
  105 FORMAT (1X,I4,2E15.6,4X,8I5/(1X,38X,8I5))
  106 FORMAT (1X,I4,3E15.6,4X,5I5/(1X,53X,5I5))
  107 FORMAT (1X)
  108 FORMAT (///)
  109 FORMAT (/1X,'NB = ',I4,' Boundary Nodes',5X,
     .        'NE = ',I5,' Edges',5X,'NT = ',I5,
     .        ' Triangles')
  110 FORMAT (1X,10X,'*** N is outside its valid',
     .        ' range ***')
      END
      SUBROUTINE VLIST (N,X,Y,Z,LIST,LPTR,LEND, NV,LISTV,
     .                  DXL,DYL,IER)
      INTEGER N, LIST(*), LPTR(*), LEND(N), NV, LISTV(*),
     .        IER
      DOUBLE PRECISION X(N), Y(N), Z(N), DXL(*), DYL(*)
C
C***********************************************************
C
C                                              From CSRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/21/96
C
C   Given a convexity-preserving triangulation of a set of
C nodes and data values, this subroutine constructs its
C dual, a gradient feasibility diagram.  The feasible nodal
C gradients are those for which there exists a convex Her-
C mite interpolant of the data values and gradients.
C
C   The feasibility region associated with node Pi is the
C set of gradients G for which <G,Pj-Pi> <= Zj-Zi for all
C j such that Pj is a neighbor of Pi.  The feasibility re-
C gions are convex and, if the data is strictly convex (no
C four data points are coplanar), have nonzero (finite or
C infinite) area.  Furthermore, the set of regions parti-
C tions the plane, and the diagram is dual to the
C triangulation:  there are one-to-one correspondences be-
C tween feasibility regions and nodes, edges of feasibility
C regions and triangulation edges (the two being perpendicu-
C lar to each other), and vertices of feasibility regions
C and triangles.  The feasibility region vertices are the
C gradients of the linear polynomials that comprise the
C triangle-based piecewise linear interpolant of the data
C values.
C
C   While the feasibility regions associated with interior
C nodes have finite area, those associated with boundary
C nodes are unbounded and have infinite-length edges.  In
C order to define these edges and to facilitate the selec-
C tion of nodal gradients as centroids, each feasibility
C region edge corresponding to a boundary edge is terminated
C at a pseudo-vertex chosen so that its length is the aver-
C age of the lengths of the other two edges that share the
C same vertex.
C
C
C On input:
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C           Note that, if N = 3, the three truncated feasi-
C           bility regions degenerate to the same point (the
C           gradient of the linear interpolant).
C
C       X,Y,Z = Arrays of length N containing the nodal
C               coordinates (X,Y) and data values Z.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to Subroutine
C                        TRMSHC.
C
C The above parameters are not altered by this routine.
C
C       LISTV = Integer array of length at least 6*N-2*NB-6,
C               where NB is the number of boundary nodes.
C
C       DXL,DYL = Arrays of length at least 2*N-2.
C
C On output:
C
C       NV = Number of feasibility region vertices, includ-
C            ing pseudo-vertices associated with boundary
C            edges:  NV = NT+NB = 2*N-2, where NT = 2*N-NB-2
C            is the number of triangles.
C
C       LISTV = Array containing vertex indexes (indexes to
C               DXL and DYL) stored in one-to-one corres-
C               pondence with LIST/LPTR entries.  The index
C               associated with triangle (N1,N2,N3) is
C               stored in LISTV(I), LISTV(J), and LISTV(K),
C               where LIST(I), LIST(J), and LIST(K) are the
C               indexes of N2 as a neighbor of N1, N3 as a
C               neighbor of N2, and N1 as a neighbor of N3.
C               Pseudo-vertices have indexes in the range
C               NT+1 to NV and each of these is stored only
C               once in the position associated with a boun-
C               dary edge N1->N2:  LISTV(LEND(N2)), where
C               LIST(LEND(N2)) = -N1.  Thus, the feasibility
C               region associated with a node is defined by
C               the CCW-ordered sequence of vertices in
C               one-to-one correspondence with its adjacency
C               list.  The truncated feasibility region
C               associated with a boundary node has two
C               pseudo-vertices (including the one stored
C               with its first neighbor).
C
C       DXL,DYL = Arrays of length NV containing the ver-
C                 tices (components of the gradients
C                 associated with triangles) in the first
C                 NT positions and pseudo-vertices in the
C                 last NB positions.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.  No other parameters are
C                     altered in this case.
C             IER = 2 if a nonpositive triangle area was
C                     returned by Subroutine LGRAD implying
C                     an invalid triangulation.  Output
C                     parameters are not valid in this case.
C
C Modules required by VLIST:  LGRAD, LSTPTR
C
C Intrinsic functions called by VLIST:  ABS, SQRT
C
C***********************************************************
C
      INTEGER LSTPTR
      DOUBLE PRECISION D, D1, D2, RN, RX, RY, SA, XN1, YN1,
     .                 ZN1
      INTEGER KV, KV1, KV2, L1, L2, L3, LP, LP1F, LP2, LP2L,
     .        LP3, LPL, N0, N1, N2, N3, NM2, NT
C
C Local parameters:
C
C D =           Distance between a vertex and a pseudo-
C                 vertex -- average of D1 and D2
C D1,D2 =       Distances from vertex KV to KV1 and to KV2,
C                 respectively
C KV =          Index of a vertex
C KV1,KV2 =     Indexes of the vertices associated with the
C                 triangles opposite boundary nodes N1 and
C                 N2, respectively
C L1,L2,L3 =    Pointers (LIST/LISTV indexes) to neighbors
C                 of N2
C LP =          Pointer
C LP1F =        Pointer to the first neighbor of N1
C LP2 =         Pointer to N2 as a neighbor of N1
C LP2L =        Pointer to the last neighbor of N2
C LP3 =         Pointer to N3 as a neighbor of N1
C LPL =         Pointer to the last neighbor of N1
C N0 =          Index of a boundary node
C N1,N2,N3 =    Nodal indexes of the vertices of a triangle
C                 (in CCW order)
C NM2 =         N-2
C NT =          Number of triangles in the triangulation (or
C                 number of triangles whose feasible region
C                 vertices have been computed
C RN =          Length of boundary edge N1->N2
C RX,RY =       X and y components of boundary edge N1->N2
C SA =          Signed triangle area returned by LGRAD
C XN1,YN1,ZN1 = Nodal coordinates and data value for node N1
C
      IF (N .LT. 3) GO TO 11
      NT = 0
C
C Loop on triangles (N1,N2,N3) such that N2 and N3 have
C   larger indexes than N1.
C
      NM2 = N - 2
      DO 3 N1 = 1,NM2
        XN1 = X(N1)
        YN1 = Y(N1)
        ZN1 = Z(N1)
C
C Loop on neighbors of N1.
C
        LPL = LEND(N1)
        LP2 = LPL
    1   LP2 = LPTR(LP2)
          N2 = LIST(LP2)
          LP3 = LPTR(LP2)
          N3 = ABS(LIST(LP3))
          IF (N2 .LT. N1  .OR.  N3 .LT. N1) GO TO 2
          NT = NT + 1
C
C   Store the triangle index NT in the appropriate LISTV
C     locations.
C
          LISTV(LP2) = NT
          LP = LSTPTR (LEND(N2),N3,LIST,LPTR)
          LISTV(LP) = NT
          LP = LSTPTR (LEND(N3),N1,LIST,LPTR)
          LISTV(LP) = NT
C
C   Compute and store vertex NT.
C
          CALL LGRAD (XN1,YN1,ZN1,X(N2),Y(N2),Z(N2),X(N3),
     .                Y(N3),Z(N3), DXL(NT),DYL(NT),SA)
          IF (SA .LE. 0.) GO TO 12
C
C   Bottom of loop on neighbors.
C
    2     IF (LP2 .NE. LPL) GO TO 1
    3   CONTINUE
      NV = NT
C
C Add the NB pseudo-vertices.
C
C   Set N0 and N2 to the first boundary node encountered.
C
      N0 = 0
    4 N0 = N0 + 1
        LP2L = LEND(N0)
        IF (LIST(LP2L) .GT. 0) GO TO 4
      N2 = N0
C
C   CCW loop on boundary edges N1->N2:  LP1F and LP2L point
C     to the first neighbor of N1 and last neighbor of N2,
C     respectively.
C
    5 N1 = N2
        LP1F = LPTR(LP2L)
        N2 = LIST(LP1F)
        LP2L = LEND(N2)
C
C   Set KV to the index of the vertex (gradient) associated
C     with the triangle containing N1->N2, set KV1 and KV2
C     to the indexes of the vertices associated with the
C     triangles opposite N1 and N2, respectively (if they
C     exist), and compute the distances D1 and D2 from ver-
C     tex KV to vertices KV1 and KV2 (with zero values in
C     the case of no triangle).  Note that D1 = D2 = 0 if
C     N = 3.
C
        KV = LISTV(LP1F)
C
C   Compute D1.  Find the pointer L1 to the third from last
C     neighbor of N2 if it exists.
C
        L1 = LPTR(LP2L)
        L2 = LPTR(L1)
        IF (L2 .EQ. LP2L) THEN
C
C   N2 has only two neighbors:  no triangle opposite N1.
C
          D1 = 0.
        ELSE
C
C   Loop on triples of adjacent neighbors of N2 until L3
C     points to the last neighbor.
C
    6     L3 = LPTR(L2)
            IF (L3 .EQ. LP2L) GO TO 7
            L1 = L2
            L2 = L3
            GO TO 6
C
    7     KV1 = LISTV(L1)
          D1 = SQRT( (DXL(KV)-DXL(KV1))**2 +
     .               (DYL(KV)-DYL(KV1))**2 )
        ENDIF
C
C   Compute D2.  Set LP to the pointer to the second
C     neighbor of N1.
C
        LP = LPTR(LP1F)
        IF (LP .EQ. LEND(N1)) THEN
C
C   N1 has only two neighbors:  no triangle opposite N2.
C
          D2 = 0.
        ELSE
          KV2 = LISTV(LP)
          D2 = SQRT( (DXL(KV)-DXL(KV2))**2 +
     .               (DYL(KV)-DYL(KV2))**2 )
        ENDIF
C
C   Compute the average D of D1 and D2.
C
        D = (D1+D2)/2.D0
C
C   Compute and store the pseudo-vertex V3 associated with
C     boundary edge N1->N2.  V3 = V + D*R, where V is vertex
C     KV and R is a unit vector orthogonal to the edge.
C
        RX = X(N2)-X(N1)
        RY = Y(N2)-Y(N1)
        RN = SQRT(RX*RX + RY*RY)
        NV = NV + 1
        DXL(NV) = DXL(KV) + D*RY/RN
        DYL(NV) = DYL(KV) - D*RX/RN
        LISTV(LP2L) = NV
C
C   Bottom of loop on boundary edges.
C
        IF (N2 .NE. N0) GO TO 5
C
C No error encountered.
C
      IER = 0
      RETURN
C
C N < 3.
C
   11 IER = 1
      RETURN
C
C Triangle (N1,N2,N3) has nonpositive area.
C
   12 IER = 2
      RETURN
      END
