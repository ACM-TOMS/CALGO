      SUBROUTINE TSHEP2 (N,X,Y,F,NC,NW,NR, LCELL,LNEXT,XMIN,
     .                   YMIN,DX,DY,SX,SY,RMAX,RW,A,IER)
      INTEGER N, NC, NW, NR, LCELL(NR,NR), LNEXT(N), IER
      DOUBLE PRECISION  X(N), Y(N), F(N), XMIN, YMIN, DX,
     .                  DY, SX, SY, RMAX, RW(N), A(10,N)
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/23/97
C
C   This subroutine computes a set of parameters defining a
C C2 (twice continuously differentiable) bivariate function
C C(X,Y) which interpolates data values F at a set of N
C arbitrarily distributed points (X,Y) in the plane (nodes).
C The interpolant C may be evaluated at an arbitrary point
C by function TS2VAL, and its first partial derivatives are
C computed by Subroutine TS2GRD.
C
C   The interpolation scheme is a modified Shepard method
C based on a cosine series:
C
C C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]
C
C for bivariate functions W(k) and C(k).  The nodal func-
C tions are given by
C
C  C(k)(x,y) = A(1,k) + A(2,k)*cos(px) + A(3,k)*cos(py) +
C              A(4,k)*cos(2*px) + A(5,k)*cos(px)*cos(py) +
C              A(6,k)*cos(2*py) + A(7,k)*cos(3*px) +
C              A(8,k)*cos(2*px)*cos(py) +
C              A(9,k)*cos(px)*cos(2*py) + A(10,k)*cos(3*py),
C
C for px = ((x-XMIN)/(XMAX-XMIN))*Pi and
C     py = ((y-YMIN)/(YMAX-YMIN))*Pi ,
C
C where [XMIN,XMAX] X [YMIN,YMAX] is the smallest rectangle
C that contains the nodes.  The method exactly reproduces
C true cubic functions of cos(px) and cos(py).
C
C The coefficients A(,k) of C(k) are obtained by a weighted
C least squares fit to the closest NC data points (along
C with point k) with weights related to inverse distance.
C Note that the radius of influence for the least squares
C fit is fixed for each k, but varies with k.
C
C The weights are taken to be
C
C   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,
C
C where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is
C the Euclidean distance between (x,y) and (X(k),Y(k)).  The
C radius of influence R(k) varies with k and is chosen so
C that NW nodes are within the radius.  Note that W(k) is
C not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
C as (x,y) approaches (X(k),Y(k)).
C
C On input:
C
C       N = Number of nodes and data values.  N .GE. 10.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       F = Array of length N containing the data values
C           in one-to-one correspondence with the nodes.
C
C       NC = Number of data points (in addition to point K)
C            to be used in the least squares fit for coeffi-
C            cients defining the nodal functions C(k).
C            Values found to be optimal for test data sets
C            ranged from 14 to 19.  The recommended value is
C            NC = 18, and NC must be in the range 9 to
C            Min(40,N-1).
C
C       NW = Number of nodes within (and defining) the radii
C            of influence R(k) which enter into the weights
C            W(k).  For N sufficiently large, a recommended
C            value is NW = 32.  In general, NW should be at
C            least 1.5*NC.  1 .LE. NW .LE. Min(40,N-1).
C
C       NR = Number of rows and columns in the cell grid de-
C            fined in Subroutine STORE2.  A rectangle con-
C            taining the nodes is partitioned into cells in
C            order to increase search efficiency.  NR =
C            Sqrt(N/3) is recommended.  NR .GE. 1.
C
C The above parameters are not altered by this routine.
C
C       LCELL = Array of length .GE. NR**2.
C
C       LNEXT = Array of length .GE. N.
C
C       RW = Array of length .GE. N.
C
C       A = Array of length .GE. 10N.
C
C On output:
C
C       LCELL = NR by NR array of nodal indexes associated
C               with cells.  Refer to Subroutine STORE2.
C
C       LNEXT = Array of length N containing next-node
C               indexes.  Refer to Subroutine STORE2.
C
C       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
C                         dimensions.  Refer to Subroutine
C                         STORE2.
C
C       SX,SY = Scale factors for mapping [XMIN,XMAX] and
C               [YMIN,YMAX] to [0,Pi]:
C               SX = Pi/(XMAX-XMIN) and SY = Pi/(YMAX-YMIN).
C
C       RMAX = Largest element in RW -- maximum radius R(k).
C
C       RW = Array containing the the radii R(k) which enter
C            into the weights W(k).
C
C       A = 10 by N array containing the coefficients for
C           nodal function C(k) in column k.
C
C   Note that the output parameters described above are not
C defined unless IER = 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, NC, NW, or NR is outside its
C                     valid range.
C             IER = 2 if duplicate nodes were encountered.
C             IER = 3 if all nodes are collinear.
C
C Modules required by TSHEP2:  GETNP2, GIVENS, ROTATE,
C                                SETUPC, STORE2
C
C Intrinsic functions called by TSHEP2:  ABS, ACOS, DBLE,
C                                          MAX, MIN, SQRT
C
C***********************************************************
C
      INTEGER LMX
      PARAMETER (LMX=40)
      INTEGER I, IERR, IP1, IRM1, IROW, J, JP1, K, LMAX,
     .        LNP, NEQ, NN, NNC, NNW, NP, NPTS(LMX),
     .        NCWMAX
      DOUBLE PRECISION B(11,11), C, D, DMIN, DTOL, FK, PI,
     .                 RC, RS(0:LMX), RSMX, RTOL, RWS, S,
     .                 STF, T, W, WMAX, WSF, XK, XMAX, XNP,
     .                 YK, YMAX, YNP
C
      DATA    DTOL/.01/,  RTOL/1.D-5/,  WSF/1.D10/
C
C Local parameters:
C
C B =          Transpose of the augmented regression matrix
C C =          First component of the plane rotation used to
C                zero the lower triangle of B**T -- computed
C                by Subroutine GIVENS
C D =          Distance between nodes K and NP
C DMIN =       Minimum of the magnitudes of the diagonal
C                elements of the regression matrix after
C                zeros are introduced below the diagonal
C DTOL =       Tolerance for detecting an ill-conditioned
C                system.  The system is accepted when
C                DMIN .GE. DTOL.
C FK =         Data value at mode K -- F(K)
C I =          Index for A, B, and NPTS
C IERR =       Error flag for the call to Subroutine STORE2
C IP1 =        I+1
C IRM1 =       IROW-1
C IROW =       Row index for B**T
C J =          Index for A and B
C JP1 =        J+1
C K =          Nodal function index and column index for A
C LMAX =       Maximum number of NPTS elements
C LMX =        Maximum value of LMAX (value of LMAX for N
C                sufficiently large
C LNP =        Current length of NPTS
C NEQ =        Number of equations in the least squares fit
C NN,NNC,NNW = Local copies of N, NC, and NW
C NP =         NPTS element
C NPTS =       Array containing the indexes of a sequence of
C                nodes to be used (along with node K) in the
C                least squares fit associated with C(k), or
C                to compute RW.  The nodes are ordered by
C                distance from K, and the last element
C                (usually indexed by LNP) is used only to
C                determine RC, or RW(K) if NW > NC.
C NCWMAX =     Max(NC,NW)
C PI =         Pi
C RC =         Radius of influence which enters into the
C                weights W for C(K)
C RS =         Array containing squared distances between
C                nodes K and NPTS(LNP) -- used to compute
C                RC and RW(K)
C RSMX =       Maximum squared RW element encountered
C RTOL =       Tolerance for detecting a sufficiently large
C                relative change in consecutive RS elements.
C                If the change is not greater than RTOL, the
C                nodes are treated as being the same
C                distance from K
C RWS =        Current squared value of RW(K)
C STF =        Marquardt stabilization factor used to damp
C                out the last four solution components when
C                the system is ill-conditioned.
C T =          Temporary variable for accumulating a scalar
C                product in the back solve
C W =          Weight associated with node NP in the least
C                squares system for nodal function K:  W =
C                (RC-D)/(RC*D) = 1/D - 1/RC
C WMAX =       Weight associated with node K in the least
C                squares system for nodal function K.  A
C                large weight is used to force interpolation
C                of the data value F(K).  WMAX is the
C                largest W value scaled by WSF.
C WSF =        Scale factor for WMAX
C XK,YK =      Coordinates of node K -- X(K), Y(K)
C XMAX,YMAX =  Maximum nodal coordinates
C XNP,YNP =    Coordinates of node NP transformed by the
C                affine transformation that maps [XMIN,XMAX]
C                X [YMIN,YMAX] to [0,Pi] X [0,Pi]
C
      NN = N
      NNC = NC
      NNW = NW
      NCWMAX = MAX(NNC,NNW)
      LMAX = MIN(LMX,NN-1)
      IF (NNC .LT. 9  .OR.  NNW .LT. 1  .OR.  NCWMAX .GT.
     .    LMAX  .OR.  NR .LT. 1) GO TO 21
C
C Create the cell data structure, compute XMAX, YMAX, SX,
C   and SY, and initialize RSMX.
C
      CALL STORE2 (NN,X,Y,NR, LCELL,LNEXT,XMIN,YMIN,DX,DY,
     .             IERR)
      XMAX = XMIN + DBLE(NR)*DX
      YMAX = YMIN + DBLE(NR)*DY
      PI = ACOS(-1.D0)
      SX = PI/(XMAX-XMIN)
      SY = PI/(YMAX-YMIN)
      IF (IERR .NE. 0) GO TO 23
      RSMX = 0.
C
C Outer loop on node K:
C
      DO 15 K = 1,NN
        XK = X(K)
        YK = Y(K)
        FK = F(K)
C
C Mark node K to exclude it from the search for nearest
C   neighbors.
C
        LNEXT(K) = -LNEXT(K)
C
C Initialize for loop on NPTS elements.
C
        RWS = 0.
        RC = 0.
        LNP = 0
        RS(0) = 0.
C
C Compute NPTS, LNP, RWS, NEQ, and RC.
C
    1     IF (LNP .EQ. LMAX) GO TO 2
          LNP = LNP + 1
          CALL GETNP2 (XK,YK,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,
     .                 DX,DY, NP,RS(LNP))
          IF (RS(LNP) .EQ. 0.) GO TO 22
          NPTS(LNP) = NP
          IF ((RS(LNP)-RS(LNP-1))/RS(LNP) .LT. RTOL) GO TO 1
          IF (RWS .EQ. 0.  .AND.  LNP .GT. NNW)
     .      RWS = RS(LNP)
          IF (RC .EQ. 0.  .AND.  LNP .GT. NNC) THEN
C
C   RC = 0 (not yet computed) and LNP > NC.  RC =
C     Sqrt(RS(LNP)) is sufficiently large to (strictly)
C     include NC nodes.  The least squares fit will
C     include NEQ = LNP equations for 9 .LE. NC .LT. NEQ
C     .LE. LMAX .LE. N-1.
C
            NEQ = LNP
            RC = SQRT(RS(LNP))
          ENDIF
C
C   Bottom of loop -- test for termination.
C
          IF (LNP .GT. NCWMAX) GO TO 3
          GO TO 1
C
C All LMAX nodes are included in NPTS.  RWS and/or RC**2 is
C   (arbitrarily) taken to be 10 percent larger than the
C   distance to the last node included.
C
    2   IF (RWS .EQ. 0.) RWS = 1.1*RS(LMAX)
        IF (RC .EQ. 0.) THEN
          NEQ = LMAX + 1
          RC = SQRT(1.1*RS(LMAX))
        ENDIF
C
C Store RW(K) and update RSMX if necessary.
C
    3   RW(K) = SQRT(RWS)
        IF (RWS .GT. RSMX) RSMX = RWS
C
C A Q-R decomposition is used to solve the least squares
C   system.  The transpose of the augmented regression
C   matrix is stored in B.
C
C Store the equation associated with node K in the first
C   row.
C
        WMAX = WSF*(1./SQRT(RS(1)) - 1./RC)
        XNP = SX*(XK-XMIN)
        YNP = SY*(YK-YMIN)
        CALL SETUPC (XNP,YNP,FK,WMAX, B(1,1))
C
C Set up the remaining equations and zero out the lower
C   triangle with Givens rotations.
C
        I = 0
    4     I = I + 1
          NP = NPTS(I)
          IROW = MIN(I+1,11)
          D = SQRT(RS(I))
          W = 1./D - 1./RC
          XNP = SX*(X(NP)-XMIN)
          YNP = SY*(Y(NP)-YMIN)
          CALL SETUPC (XNP,YNP,F(NP),W, B(1,IROW))
          IRM1 = IROW-1
          DO 5 J = 1,IRM1
            JP1 = J + 1
            CALL GIVENS (B(J,J),B(J,IROW),C,S)
            CALL ROTATE (11-J,C,S,B(JP1,J),B(JP1,IROW))
    5       CONTINUE
          IF (I+1 .LT. NEQ) GO TO 4
C
C Test the system for ill-conditioning.
C
        DMIN = MIN( ABS(B(1,1)),ABS(B(2,2)),ABS(B(3,3)),
     .              ABS(B(4,4)),ABS(B(5,5)),ABS(B(6,6)),
     .              ABS(B(7,7)),ABS(B(8,8)),ABS(B(9,9)),
     .              ABS(B(10,10)) )
        IF (DMIN .GE. DTOL) GO TO 11
        IF (NEQ .GT. LMAX) GO TO 7
C
C Increase RC and add another equation to the system to
C   improve the conditioning.  The number of NPTS elements
C   is also increased if necessary.
C
    6   NEQ = NEQ + 1
        IF (NEQ .EQ. LMAX+1) THEN
          RC = SQRT(1.1*RS(LMAX))
          GO TO 4
        ENDIF
        IF (NEQ .LE. LNP) THEN
C
C   NEQ .LE. LNP.
C
          IF ((RS(NEQ)-RS(NEQ-1))/RS(NEQ) .LT. RTOL) GO TO 6
          RC = SQRT(RS(NEQ))
          GO TO 4
        ENDIF
C
C   NEQ = LNP+1.  Add an element to NPTS.
C
        LNP = LNP + 1
        CALL GETNP2 (XK,YK,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,
     .               DX,DY, NP,RS(LNP))
        IF (NP .EQ. 0) GO TO 22
        NPTS(LNP) = NP
        IF ( (RS(LNP)-RS(LNP-1))/RS(LNP) .LT. RTOL ) GO TO 6
        RC = SQRT(RS(LNP))
        GO TO 4
C
C Stabilize the system by damping the last four basis
C   functions -- add multiples of the last four unit
C   vectors to the last four equations.
C
    7   STF = 1.0/RC
        DO 10 I = 7,10
          B(I,11) = STF
          IP1 = I + 1
          DO 8 J = IP1,11
            B(J,11) = 0.
    8       CONTINUE
          DO 9 J = I,10
            JP1 = J + 1
            CALL GIVENS (B(J,J),B(J,11),C,S)
            CALL ROTATE (11-J,C,S,B(JP1,J),B(JP1,11))
    9       CONTINUE
   10     CONTINUE
C
C Test the damped system for ill-conditioning.
C
        DMIN = MIN( ABS(B(5,5)),ABS(B(6,6)),ABS(B(7,7)),
     .              ABS(B(8,8)),ABS(B(9,9)),ABS(B(10,10)) )
        IF (DMIN .LT. DTOL) GO TO 23
C
C Solve the order-10 triangular system for the coefficients.
C
   11   DO 13 I = 10,1,-1
          T = 0.
          IF (I .NE. 10) THEN
            IP1 = I + 1
            DO 12 J = IP1,10
              T = T + B(J,I)*A(J,K)
   12         CONTINUE
          ENDIF
          A(I,K) = (B(11,I)-T)/B(I,I)
   13     CONTINUE
C
C Unmark K and the elements of NPTS.
C
        LNEXT(K) = -LNEXT(K)
        DO 14 I = 1,LNP
          NP = NPTS(I)
          LNEXT(NP) = -LNEXT(NP)
   14     CONTINUE
   15   CONTINUE
C
C No errors encountered.
C
      RMAX = SQRT(RSMX)
      IER = 0
      RETURN
C
C N, NC, NW, or NR is outside its valid range.
C
   21 IER = 1
      RETURN
C
C Duplicate nodes were encountered by GETNP2.
C
   22 IER = 2
      RETURN
C
C No unique solution due to collinear nodes.
C
   23 IER = 3
      RETURN
      END
      DOUBLE PRECISION FUNCTION TS2VAL (PX,PY,N,X,Y,F,NR,
     .                 LCELL,LNEXT,XMIN,YMIN,DX,DY,SX,SY,
     .                 RMAX,RW,A)
      INTEGER N, NR, LCELL(NR,NR), LNEXT(N)
      DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     .                 DX, DY, SX, SY, RMAX, RW(N), A(10,N)
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/19/97
C
C   This function returns the value C(PX,PY), where C is the
C weighted sum of nodal functions defined in Subroutine
C TSHEP2.  TS2GRD may be called to compute a gradient
C of C along with the value, and/or to test for errors.
C TS2HES may be called to compute a value, first partial
C derivatives, and second partial derivatives at a point.
C
C On input:
C
C       PX,PY = Cartesian coordinates of the point P at
C               which C is to be evaluated.
C
C       N = Number of nodes and data values defining C.
C           N .GE. 10.
C
C       X,Y,F = Arrays of length N containing the nodes and
C               data values interpolated by C.
C
C       NR = Number of rows and columns in the cell grid.
C            Refer to Subroutine STORE2.  NR .GE. 1.
C
C       LCELL = NR by NR array of nodal indexes associated
C               with cells.  Refer to Subroutine STORE2.
C
C       LNEXT = Array of length N containing next-node
C               indexes.  Refer to Subroutine STORE2.
C
C       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
C                         dimensions.  DX and DY must be
C                         positive.  Refer to Subroutine
C                         STORE2.
C
C       SX,SY = Scale factors for mapping [XMIN,XMAX] and
C               [YMIN,YMAX] to [0,Pi]:
C               SX = Pi/(XMAX-XMIN) and SY = Pi/(YMAX-YMIN).
C
C       RMAX = Largest element in RW -- maximum radius R(k).
C
C       RW = Array containing the the radii R(k) which enter
C            into the weights W(k) defining C.
C
C       A = 10 by N array containing the coefficients for
C           nodal function C(k) in column k.
C
C   Input parameters are not altered by this function.  The
C parameters other than PX and PY should be input unaltered
C from their values on output from TSHEP2.  This function
C should not be called if a nonzero error flag was returned
C by TSHEP2.
C
C On output:
C
C       TS2VAL = Function value C(PX,PY) unless N, NR, DX,
C                DY, or RMAX is invalid, in which case no
C                value is returned.
C
C Modules required by TS2VAL:  NONE
C
C Intrinsic functions called by TS2VAL:  COS, INT, SQRT
C
C***********************************************************
C
      INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      DOUBLE PRECISION CX1, CX2, CX3, CY1, CY2, CY3, D, R,
     .                 SW, SWC, W, XP, YP
C
C Local parameters:
C
C CX1 =       Cos(XP)
C CX2 =       Cos(2*XP)
C CX3 =       Cos(3*XP)
C CY1 =       Cos(YP)
C CY2 =       Cos(2*YP)
C CY3 =       Cos(3*YP)
C D =         Distance between P and node K
C I =         Cell row index in the range IMIN to IMAX
C IMIN,IMAX = Range of cell row indexes of the cells
C               intersected by a disk of radius RMAX
C               centered at P
C J =         Cell column index in the range JMIN to JMAX
C JMIN,JMAX = Range of cell column indexes of the cells
C               intersected by a disk of radius RMAX
C               centered at P
C K =         Index of a node in cell (I,J)
C KP =        Previous value of K in the sequence of nodes
C               in cell (I,J)
C R =         Radius of influence for node K
C SW =        Sum of weights W(K)
C SWC =       Sum of weighted nodal function values at P
C W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
C               where (R-D)+ = 0 if R < D
C XP,YP =    Coordinates of node NP transformed by the
C                affine transformation that maps [XMIN,XMAX]
C                X [YMIN,YMAX] to [0,Pi] X [0,Pi]
C
      XP = SX*(PX-XMIN)
      YP = SY*(PY-YMIN)
      IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR.
     .    DY .LE. 0.  .OR.  RMAX .LT. 0.) RETURN
C
C Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
C   the range of the search for nodes whose radii include
C   P.  The cells which must be searched are those inter-
C   sected by (or contained in) a circle of radius RMAX
C   centered at P.
C
      IMIN = INT((PX-XMIN-RMAX)/DX) + 1
      IMAX = INT((PX-XMIN+RMAX)/DX) + 1
      IF (IMIN .LT. 1) IMIN = 1
      IF (IMAX .GT. NR) IMAX = NR
      JMIN = INT((PY-YMIN-RMAX)/DY) + 1
      JMAX = INT((PY-YMIN+RMAX)/DY) + 1
      IF (JMIN .LT. 1) JMIN = 1
      IF (JMAX .GT. NR) JMAX = NR
C
C The following is a test for no cells within the circle
C   of radius RMAX.
C
      IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 6
C
C Compute trig function values at (XP,YP).
C
      CX1 = COS(XP)
      CY1 = COS(YP)
      CX2 = COS(2.0*XP)
      CY2 = COS(2.0*YP)
      CX3 = COS(3.0*XP)
      CY3 = COS(3.0*YP)
C
C Accumulate weight values in SW and weighted nodal function
C   values in SWC.  The weights are W(K) = ((R-D)+/(R*D))**3
C   for R = RW(K) and D = distance between P and node K.
C
      SW = 0.
      SWC = 0.
C
C Outer loop on cells (I,J).
C
      DO 4 J = JMIN,JMAX
        DO 3 I = IMIN,IMAX
          K = LCELL(I,J)
          IF (K .EQ. 0) GO TO 3
C
C Inner loop on nodes K.
C
    1     D = SQRT((PX-X(K))**2 + (PY-Y(K))**2)
          R = RW(K)
          IF (D .GE. R) GO TO 2
          IF (D .EQ. 0.) GO TO 5
          W = (1.0/D - 1.0/R)**3
          SW = SW + W
          SWC = SWC + W*(A(1,K) + A(2,K)*CX1 + A(3,K)*CY1 +
     .                   A(4,K)*CX2 + A(5,K)*CX1*CY1 +
     .                   A(6,K)*CY2 + A(7,K)*CX3 +
     .                   A(8,K)*CX2*CY1 + A(9,K)*CX1*CY2 +
     .                   A(10,K)*CY3)
C
C Bottom of loop on nodes in cell (I,J).
C
    2     KP = K
          K = LNEXT(KP)
          IF (K .NE. KP) GO TO 1
    3     CONTINUE
    4   CONTINUE
C
C SW = 0 iff P is not within the radius R(K) for any node K.
C
      IF (SW .EQ. 0.) GO TO 6
      TS2VAL = SWC/SW
      RETURN
C
C (PX,PY) = (X(K),Y(K)).
C
    5 TS2VAL = F(K)
      RETURN
C
C All weights are 0 at P.
C
    6 TS2VAL = 0.
      RETURN
      END
      SUBROUTINE TS2GRD (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     .                   YMIN,DX,DY,SX,SY,RMAX,RW,A, C,CX,
     .                   CY,IER)
      INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
      DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     .                 DX, DY, SX, SY, RMAX, RW(N), A(10,N),
     .                 C, CX, CY
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/19/97
C
C   This subroutine computes the value and gradient at P =
C (PX,PY) of the interpolatory function C defined in Sub-
C routine TSHEP2.  C is a weighted sum of cosine series
C nodal functions.
C
C On input:
C
C       PX,PY = Cartesian coordinates of the point P at
C               which C and its partial derivatives are
C               to be evaluated.
C
C       N = Number of nodes and data values defining C.
C           N .GE. 10.
C
C       X,Y,F = Arrays of length N containing the nodes and
C               data values interpolated by C.
C
C       NR = Number of rows and columns in the cell grid.
C            Refer to Subroutine STORE2.  NR .GE. 1.
C
C       LCELL = NR by NR array of nodal indexes associated
C               with cells.  Refer to Subroutine STORE2.
C
C       LNEXT = Array of length N containing next-node
C               indexes.  Refer to Subroutine STORE2.
C
C       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
C                         dimensions.  DX and DY must be
C                         positive.  Refer to Subroutine
C                         STORE2.
C
C       SX,SY = Scale factors for mapping [XMIN,XMAX] and
C               [YMIN,YMAX] to [0,Pi]:
C               SX = Pi/(XMAX-XMIN) and SY = Pi/(YMAX-YMIN).
C
C       RMAX = Largest element in RW -- maximum radius R(k).
C
C       RW = Array of length N containing the the radii R(k)
C            which enter into the weights W(k) defining C.
C
C       A = 10 by N array containing the coefficients for
C           nodal function C(k) in column k.
C
C   Input parameters are not altered by this subroutine.
C The parameters other than PX and PY should be input
C unaltered from their values on output from TSHEP2.  This
C subroutine should not be called if a nonzero error flag
C was returned by TSHEP2.
C
C On output:
C
C       C = Value of C at (PX,PY) unless IER .EQ. 1, in
C           which case no values are returned.
C
C       CX,CY = First partial derivatives of C at (PX,PY)
C               unless IER .EQ. 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, NR, DX, DY or RMAX is invalid.
C             IER = 2 if no errors were encountered but
C                     (PX,PY) is not within the radius R(k)
C                     for any node k (and thus C=CX=CY=0).
C
C Modules required by TS2GRD:  None
C
C Intrinsic functions called by TS2GRD:  COS, INT, SIN, SQRT
C
C***********************************************************
C
      INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      DOUBLE PRECISION CK, CKX, CKY, CX1, CX2, CX3, CY1,
     .                 CY2, CY3, D, DELX, DELY, R, SW, SWC,
     .                 SWCX, SWCY, SWS, SWX, SWY, SX1, SX2,
     .                 SX3, SY1, SY2, SY3, T, W, WX, WY, XP,
     .                 YP
C
C Local parameters:
C
C CK =          Value of nodal function C(K) at P
C CKX,CKY =     Partial derivatives of C(K) with respect to
C                 X and Y, respectively
C CX1,CX2,CX3 = Cos(XP), Cos(2*XP), and Cos(3*XP),
C                 respectively
C CY1,CY2,CY3 = Cos(YP), Cos(2*YP), and Cos(3*YP),
C                 respectively
C D =           Distance between P and node K
C DELX =        PX - X(K)
C DELY =        PY - Y(K)
C I =           Cell row index in the range IMIN to IMAX
C IMIN,IMAX =   Range of cell row indexes of the cells
C                 intersected by a disk of radius RMAX
C                 centered at P
C J =           Cell column index in the range JMIN to JMAX
C JMIN,JMAX =   Range of cell column indexes of the cells
C                 intersected by a disk of radius RMAX
C                 centered at P
C K =           Index of a node in cell (I,J)
C KP =          Previous value of K in the sequence of nodes
C                 in cell (I,J)
C R =           Radius of influence for node K
C SW =          Sum of weights W(K)
C SWC =         Sum of weighted nodal function values at P
C SWCX,SWCY =   Partial derivatives of SWC with respect to X
C                 and Y, respectively
C SWS =         SW**2
C SWX,SWY =     Partial derivatives of SW with respect to X
C                 and Y, respectively
C SX1,SX2,SX3 = Sin(XP), Sin(2*XP), and Sin(3*XP),
C                 respectively
C SY1,SY2,SY3 = Sin(YP), Sin(2*YP), and Sin(3*YP),
C                 respectively
C T =           Temporary variable
C W =           Weight W(K) value at P:  ((R-D)+/(R*D))**3,
C                 where (R-D)+ = 0 if R < D
C WX,WY =       Partial derivatives of W with respect to X
C                 and Y, respectively
C XP,YP =    Coordinates of node NP transformed by the
C                affine transformation that maps [XMIN,XMAX]
C                X [YMIN,YMAX] to [0,Pi] X [0,Pi]
C
      XP = SX*(PX-XMIN)
      YP = SY*(PY-YMIN)
      IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR.
     .    DY .LE. 0.  .OR.  RMAX .LT. 0.) GO TO 6
C
C Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
C   the range of the search for nodes whose radii include
C   P.  The cells which must be searched are those inter-
C   sected by (or contained in) a circle of radius RMAX
C   centered at P.
C
      IMIN = INT((PX-XMIN-RMAX)/DX) + 1
      IMAX = INT((PX-XMIN+RMAX)/DX) + 1
      IF (IMIN .LT. 1) IMIN = 1
      IF (IMAX .GT. NR) IMAX = NR
      JMIN = INT((PY-YMIN-RMAX)/DY) + 1
      JMAX = INT((PY-YMIN+RMAX)/DY) + 1
      IF (JMIN .LT. 1) JMIN = 1
      IF (JMAX .GT. NR) JMAX = NR
C
C The following is a test for no cells within the circle
C   of radius RMAX.
C
      IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 7
C
C Compute trig function values at (XP,YP).
C
      CX1 = COS(XP)
      CY1 = COS(YP)
      CX2 = COS(2.0*XP)
      CY2 = COS(2.0*YP)
      CX3 = COS(3.0*XP)
      CY3 = COS(3.0*YP)
      SX1 = SIN(XP)
      SY1 = SIN(YP)
      SX2 = SIN(2.0*XP)
      SY2 = SIN(2.0*YP)
      SX3 = SIN(3.0*XP)
      SY3 = SIN(3.0*YP)
C
C C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
C   from K = 1 to N, C(K) is the nodal function value,
C   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
C   ance D(K).  Thus
C
C        CX = (SWCX*SW - SWC*SWX)/SW**2  and
C        CY = (SWCY*SW - SWC*SWY)/SW**2
C
C   where SWCX and SWX are partial derivatives with respect
C   to X of SWC and SW, respectively.  SWCY and SWY are de-
C   fined similarly.
C
      SW = 0.
      SWX = 0.
      SWY = 0.
      SWC = 0.
      SWCX = 0.
      SWCY = 0.
C
C Outer loop on cells (I,J).
C
      DO 4 J = JMIN,JMAX
        DO 3 I = IMIN,IMAX
          K = LCELL(I,J)
          IF (K .EQ. 0) GO TO 3
C
C Inner loop on nodes K.
C
    1     DELX = PX - X(K)
          DELY = PY - Y(K)
          D = SQRT(DELX*DELX + DELY*DELY)
          R = RW(K)
          IF (D .GE. R) GO TO 2
C
          CK = A(1,K) + A(2,K)*CX1 + A(3,K)*CY1 +
     .         A(4,K)*CX2 + A(5,K)*CX1*CY1 + A(6,K)*CY2 +
     .         A(7,K)*CX3 + A(8,K)*CX2*CY1 +
     .         A(9,K)*CX1*CY2 + A(10,K)*CY3
          CKX = -SX*(A(2,K)*SX1 + 2.0*A(4,K)*SX2 +
     .               A(5,K)*SX1*CY1 + 3.0*A(7,K)*SX3 +
     .               2.0*A(8,K)*SX2*CY1 + A(9,K)*SX1*CY2)
          CKY = -SY*(A(3,K)*SY1 + A(5,K)*CX1*SY1 +
     .               2.0*A(6,K)*SY2 + A(8,K)*CX2*SY1 +
     .               2.0*A(9,K)*CX1*SY2 + 3.0*A(10,K)*SY3)
C
          IF (D .EQ. 0.) GO TO 5
          T = (1.0/D - 1.0/R)
          W = T**3
          T = -3.0*T*T/(D**3)
          WX = DELX*T
          WY = DELY*T
C
          SW = SW + W
          SWX = SWX + WX
          SWY = SWY + WY
          SWC = SWC + W*CK
          SWCX = SWCX + WX*CK + W*CKX
          SWCY = SWCY + WY*CK + W*CKY
C
C Bottom of loop on nodes in cell (I,J).
C
    2     KP = K
          K = LNEXT(KP)
          IF (K .NE. KP) GO TO 1
    3     CONTINUE
    4   CONTINUE
C
C SW = 0 iff P is not within the radius R(K) for any node K.
C
      IF (SW .EQ. 0.) GO TO 7
      C = SWC/SW
      SWS = SW*SW
      CX = (SWCX*SW - SWC*SWX)/SWS
      CY = (SWCY*SW - SWC*SWY)/SWS
      IER = 0
      RETURN
C
C (PX,PY) = (X(K),Y(K)).
C
    5 C = F(K)
      CX = CKX
      CY = CKY
      IER = 0
      RETURN
C
C Invalid input parameter.
C
    6 IER = 1
      RETURN
C
C No cells contain a point within RMAX of P, or
C   SW = 0 and thus D .GE. RW(K) for all K.
C
    7 C = 0.
      CX = 0.
      CY = 0.
      IER = 2
      RETURN
      END
      SUBROUTINE TS2HES (PX,PY,N,X,Y,F,NR,LCELL,LNEXT,XMIN,
     .                   YMIN,DX,DY,SX,SY,RMAX,RW,A, C,CX,
     .                   CY,CXX,CXY,CYY,IER)
      INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
      DOUBLE PRECISION PX, PY, X(N), Y(N), F(N), XMIN, YMIN,
     .                 DX, DY, SX, SY, RMAX, RW(N), A(10,N),
     .                 C, CX, CY, CXX, CXY, CYY
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/19/97
C
C   This subroutine computes the value, gradient, and
C Hessian at P = (PX,PY) of the interpolatory function C
C defined in Subroutine TSHEP2.  C is a weighted sum of
C cosine series nodal functions.
C
C On input:
C
C       PX,PY = Cartesian coordinates of the point P at
C               which C and its partial derivatives are
C               to be evaluated.
C
C       N = Number of nodes and data values defining C.
C           N .GE. 10.
C
C       X,Y,F = Arrays of length N containing the nodes and
C               data values interpolated by C.
C
C       NR = Number of rows and columns in the cell grid.
C            Refer to Subroutine STORE2.  NR .GE. 1.
C
C       LCELL = NR by NR array of nodal indexes associated
C               with cells.  Refer to Subroutine STORE2.
C
C       LNEXT = Array of length N containing next-node
C               indexes.  Refer to Subroutine STORE2.
C
C       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
C                         dimensions.  DX and DY must be
C                         positive.  Refer to Subroutine
C                         STORE2.
C
C       SX,SY = Scale factors for mapping [XMIN,XMAX] and
C               [YMIN,YMAX] to [0,Pi]:
C               SX = Pi/(XMAX-XMIN) and SY = Pi/(YMAX-YMIN).
C
C       RMAX = Largest element in RW -- maximum radius R(k).
C
C       RW = Array of length N containing the the radii R(k)
C            which enter into the weights W(k) defining C.
C
C       A = 10 by N array containing the coefficients for
C           nodal function C(k) in column k.
C
C   Input parameters are not altered by this subroutine.
C The parameters other than PX and PY should be input
C unaltered from their values on output from TSHEP2.  This
C subroutine should not be called if a nonzero error flag
C was returned by TSHEP2.
C
C On output:
C
C       C = Value of C at (PX,PY) unless IER .EQ. 1, in
C           which case no values are returned.
C
C       CX,CY = First partial derivatives of C at (PX,PY)
C               unless IER .EQ. 1.
C
C       CXX,CXY,CYY = Second partial derivatives of C at
C                     (PX,PY) unless IER .EQ. 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, NR, DX, DY or RMAX is invalid.
C             IER = 2 if no errors were encountered but
C                     (PX,PY) is not within the radius R(k)
C                     for any node k (and thus C = 0).
C
C Modules required by TS2HES:  None
C
C Intrinsic functions called by TS2HES:  COS, INT, SIN, SQRT
C
C***********************************************************
C
      INTEGER I, IMAX, IMIN, J, JMAX, JMIN, K, KP
      DOUBLE PRECISION CK, CKX, CKXX, CKXY, CKY, CKYY, CX1,
     .                 CX2, CX3, CY1, CY2, CY3, D, DELX,
     .                 DELY, DXSQ, DYSQ, R, SW, SWC, SWCX,
     .                 SWCXX, SWCXY, SWCY, SWCYY, SWS,
     .                 SWX, SWXX, SWXY, SWY, SWYY, SX1, SX2,
     .                 SX3, SY1, SY2, SY3, T1, T2, W, WX,
     .                 WXX, WXY, WY, WYY, XP, YP
C
C Local parameters:
C
C CK =             Value of nodal function C(K) at P
C CKX,CKY =        Partial derivatives of C(K) with respect
C                    to X and Y, respectively
C CKXX,CKXY,CKYY = Second partial derivatives of CK
C CX1,CX2,CX3 =    Cos(XP), Cos(2*XP), and Cos(3*XP),
C                    respectively
C CY1,CY2,CY3 =    Cos(YP), Cos(2*YP), and Cos(3*YP),
C                    respectively
C D =              Distance between P and node K
C DELX =           PX - X(K)
C DELY =           PY - Y(K)
C DXSQ,DYSQ =      DELX**2, DELY**2
C I =              Cell row index in the range IMIN to IMAX
C IMIN,IMAX =      Range of cell row indexes of the cells
C                    intersected by a disk of radius RMAX
C                    centered at P
C J =              Cell column index in the range JMIN to
C                    JMAX
C JMIN,JMAX =      Range of cell column indexes of the cells
C                    intersected by a disk of radius RMAX
C                    centered at P
C K =              Index of a node in cell (I,J)
C KP =             Previous value of K in the sequence of
C                    nodes in cell (I,J)
C R =              Radius of influence for node K
C SW =             Sum of weights W(K)
C SWC =            Sum of weighted nodal function values
C                    at P
C SWCX,SWCY =      Partial derivatives of SWC with respect
C                    to X and Y, respectively
C SWCXX,SWCXY,SWCYY = Second partial derivatives of SWC
C SWS =            SW**2
C SWX,SWY =        Partial derivatives of SW with respect to
C                    X and Y, respectively
C SWXX,SWXY,SWYY = Second partial derivatives of SW
C SX1,SX2,SX3 =    Sin(XP), Sin(2*XP), and Sin(3*XP),
C                    respectively
C SY1,SY2,SY3 =    Sin(YP), Sin(2*YP), and Sin(3*YP),
C                    respectively
C T1,T2 =          Temporary variables
C W =              Weight W(K) value at P:
C                    ((R-D)+/(R*D))**3, where (R-D)+ = 0
C                    if R < D
C WX,WY =          Partial derivatives of W with respect to
C                    X and Y, respectively
C WXX,WXY,WYY =    Second partial derivatives of W
C XP,YP =          Coordinates of node NP transformed by the
C                    affine transformation that maps
C                    [XMIN,XMAX] X [YMIN,YMAX] to [0,Pi] X
C                    [0,Pi]
C
      XP = SX*(PX-XMIN)
      YP = SY*(PY-YMIN)
      IF (N .LT. 10  .OR.  NR .LT. 1  .OR.  DX .LE. 0.  .OR.
     .    DY .LE. 0.  .OR.  RMAX .LT. 0.) GO TO 6
C
C Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining
C   the range of the search for nodes whose radii include
C   P.  The cells which must be searched are those inter-
C   sected by (or contained in) a circle of radius RMAX
C   centered at P.
C
      IMIN = INT((PX-XMIN-RMAX)/DX) + 1
      IMAX = INT((PX-XMIN+RMAX)/DX) + 1
      IF (IMIN .LT. 1) IMIN = 1
      IF (IMAX .GT. NR) IMAX = NR
      JMIN = INT((PY-YMIN-RMAX)/DY) + 1
      JMAX = INT((PY-YMIN+RMAX)/DY) + 1
      IF (JMIN .LT. 1) JMIN = 1
      IF (JMAX .GT. NR) JMAX = NR
C
C The following is a test for no cells within the circle
C   of radius RMAX.
C
      IF (IMIN .GT. IMAX  .OR.  JMIN .GT. JMAX) GO TO 7
C
C Compute trig function values at (XP,YP).
C
      CX1 = COS(XP)
      CY1 = COS(YP)
      CX2 = COS(2.0*XP)
      CY2 = COS(2.0*YP)
      CX3 = COS(3.0*XP)
      CY3 = COS(3.0*YP)
      SX1 = SIN(XP)
      SY1 = SIN(YP)
      SX2 = SIN(2.0*XP)
      SY2 = SIN(2.0*YP)
      SX3 = SIN(3.0*XP)
      SY3 = SIN(3.0*YP)
C
C C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is
C   from K = 1 to N, C(K) is the nodal function value,
C   and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and dist-
C   ance D(K).  Thus
C
C        CX = (SWCX*SW - SWC*SWX)/SW**2  and
C        CY = (SWCY*SW - SWC*SWY)/SW**2
C
C   where SWCX and SWX are partial derivatives with respect
C   to x of SWC and SW, respectively.  SWCY and SWY are de-
C   fined similarly.  The second partials are
C
C        CXX = ( SW*(SWCXX -    2*SWX*CX) - SWC*SWXX )/SW**2
C        CXY = ( SW*(SWCXY-SWX*CY-SWY*CX) - SWC*SWXY )/SW**2
C        CYY = ( SW*(SWCYY -    2*SWY*CY) - SWC*SWYY )/SW**2
C
C   where SWCXX and SWXX are second partials with respect
C   to x, SWCXY and SWXY are mixed partials, and SWCYY and
C   SWYY are second partials with respect to y.
C
      SW = 0.
      SWX = 0.
      SWY = 0.
      SWXX = 0.
      SWXY = 0.
      SWYY = 0.
      SWC = 0.
      SWCX = 0.
      SWCY = 0.
      SWCXX = 0.
      SWCXY = 0.
      SWCYY = 0.
C
C Outer loop on cells (I,J).
C
      DO 4 J = JMIN,JMAX
        DO 3 I = IMIN,IMAX
          K = LCELL(I,J)
          IF (K .EQ. 0) GO TO 3
C
C Inner loop on nodes K.
C
    1     DELX = PX - X(K)
          DELY = PY - Y(K)
          DXSQ = DELX*DELX
          DYSQ = DELY*DELY
          D = SQRT(DXSQ + DYSQ)
          R = RW(K)
          IF (D .GE. R) GO TO 2
C
          CK = A(1,K) + A(2,K)*CX1 + A(3,K)*CY1 +
     .         A(4,K)*CX2 + A(5,K)*CX1*CY1 + A(6,K)*CY2 +
     .         A(7,K)*CX3 + A(8,K)*CX2*CY1 +
     .         A(9,K)*CX1*CY2 + A(10,K)*CY3
          CKX = -SX*(A(2,K)*SX1 + 2.0*A(4,K)*SX2 +
     .               A(5,K)*SX1*CY1 + 3.0*A(7,K)*SX3 +
     .               2.0*A(8,K)*SX2*CY1 + A(9,K)*SX1*CY2)
          CKY = -SY*(A(3,K)*SY1 + A(5,K)*CX1*SY1 +
     .               2.0*A(6,K)*SY2 + A(8,K)*CX2*SY1 +
     .               2.0*A(9,K)*CX1*SY2 + 3.0*A(10,K)*SY3)
          CKXX = -SX*SX*(A(2,K)*CX1 + 4.0*A(4,K)*CX2 +
     .                   A(5,K)*CX1*CY1 + 9.0*A(7,K)*CX3 +
     .                   4.0*A(8,K)*CX2*CY1 + A(9,K)*CX1*CY2)
          CKXY = SX*SY*(A(5,K)*SX1*SY1 + 2.0*A(8,K)*SX2*SY1 +
     .                  2.0*A(9,K)*SX1*SY2)
          CKYY = -SY*SY*(A(3,K)*CY1 + A(5,K)*CX1*CY1 +
     .                   4.0*A(6,K)*CY2 + A(8,K)*CX2*CY1 +
     .                   4.0*A(9,K)*CX1*CY2 + 9.0*A(10,K)*CY3)
C
          IF (D .EQ. 0.) GO TO 5
          T1 = (1.0/D - 1.0/R)
          W = T1**3
          T2 = -3.0*T1*T1/(D**3)
          WX = DELX*T2
          WY = DELY*T2
          T1 = 3.0*T1*(2.0+3.0*D*T1)/(D**6)
          WXX = T1*DXSQ + T2
          WXY = T1*DELX*DELY
          WYY = T1*DYSQ + T2
C
          SW = SW + W
          SWX = SWX + WX
          SWY = SWY + WY
          SWXX = SWXX + WXX
          SWXY = SWXY + WXY
          SWYY = SWYY + WYY
          SWC = SWC + W*CK
          SWCX = SWCX + WX*CK + W*CKX
          SWCY = SWCY + WY*CK + W*CKY
          SWCXX = SWCXX + W*CKXX + 2.0*WX*CKX + CK*WXX
          SWCXY = SWCXY + W*CKXY + WX*CKY + WY*CKX + CK*WXY
          SWCYY = SWCYY + W*CKYY + 2.0*WY*CKY + CK*WYY
C
C Bottom of loop on nodes in cell (I,J).
C
    2     KP = K
          K = LNEXT(KP)
          IF (K .NE. KP) GO TO 1
    3     CONTINUE
    4   CONTINUE
C
C SW = 0 iff P is not within the radius R(K) for any node K.
C
      IF (SW .EQ. 0.) GO TO 7
      C = SWC/SW
      SWS = SW*SW
      CX = (SWCX*SW - SWC*SWX)/SWS
      CY = (SWCY*SW - SWC*SWY)/SWS
      CXX = (SW*(SWCXX-2.0*SWX*CX) - SWC*SWXX)/SWS
      CXY = (SW*(SWCXY-SWY*CX-SWX*CY) - SWC*SWXY)/SWS
      CYY = (SW*(SWCYY-2.0*SWY*CY) - SWC*SWYY)/SWS
      IER = 0
      RETURN
C
C (PX,PY) = (X(K),Y(K)).
C
    5 C = F(K)
      CX = CKX
      CY = CKY
      CXX = CKXX
      CXY = CKXY
      CYY = CKYY
      IER = 0
      RETURN
C
C Invalid input parameter.
C
    6 IER = 1
      RETURN
C
C No cells contain a point within RMAX of P, or
C   SW = 0 and thus D .GE. RW(K) for all K.
C
    7 C = 0.
      CX = 0.
      CY = 0.
      CXX = 0.
      CXY = 0.
      CYY = 0.
      IER = 2
      RETURN
      END
      SUBROUTINE GETNP2 (PX,PY,X,Y,NR,LCELL,LNEXT,XMIN,YMIN,
     .                   DX,DY, NP,DSQ)
      INTEGER NR, LCELL(NR,NR), LNEXT(*), NP
      DOUBLE PRECISION PX, PY, X(*), Y(*), XMIN, YMIN, DX,
     .                 DY, DSQ
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/03/97
C
C   Given a set of N nodes and the data structure defined in
C Subroutine STORE2, this subroutine uses the cell method to
C find the closest unmarked node NP to a specified point P.
C NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
C node is marked if and only if the corresponding LNEXT ele-
C ment is negative.  The absolute values of LNEXT elements,
C however, must be preserved.)  Thus, the closest M nodes to
C P may be determined by a sequence of M calls to this rou-
C tine.  Note that if the nearest neighbor to node K is to
C be determined (PX = X(K) and PY = Y(K)), then K should be
C marked before the call to this routine.
C
C   The search is begun in the cell containing (or closest
C to) P and proceeds outward in rectangular layers until all
C cells which contain points within distance R of P have
C been searched, where R is the distance from P to the first
C unmarked node encountered (infinite if no unmarked nodes
C are present).
C
C   This code is essentially unaltered from the subroutine
C of the same name in QSHEP2D.
C
C On input:
C
C       PX,PY = Cartesian coordinates of the point P whose
C               nearest unmarked neighbor is to be found.
C
C       X,Y = Arrays of length N, for N .GE. 2, containing
C             the Cartesian coordinates of the nodes.
C
C       NR = Number of rows and columns in the cell grid.
C            Refer to Subroutine STORE2.  NR .GE. 1.
C
C       LCELL = NR by NR array of nodal indexes associated
C               with cells.  Refer to Subroutine STORE2.
C
C       LNEXT = Array of length N containing next-node
C               indexes (or their negatives).  Refer to
C               Subroutine STORE2.
C
C       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
C                         dimensions.  DX and DY must be
C                         positive.  Refer to Subroutine
C                         STORE2.
C
C   Input parameters other than LNEXT are not altered by
C this routine.  With the exception of (PX,PY) and the signs
C of LNEXT elements, these parameters should be unaltered
C from their values on output from Subroutine STORE2.
C
C On output:
C
C       NP = Index (for X and Y) of the nearest unmarked
C            node to P, or 0 if all nodes are marked or NR
C            .LT. 1 or DX .LE. 0 or DY .LE. 0.  LNEXT(NP)
C            .LT. 0 IF NP .NE. 0.
C
C       DSQ = Squared Euclidean distance between P and node
C             NP, or 0 if NP = 0.
C
C Modules required by GETNP2:  None
C
C Intrinsic functions called by GETNP2:  ABS, INT, SQRT
C
C***********************************************************
C
      INTEGER I, I0, I1, I2, IMAX, IMIN, J, J0, J1, J2,
     .        JMAX, JMIN, L, LMIN, LN
      LOGICAL FIRST
      DOUBLE PRECISION DELX, DELY, R, RSMIN, RSQ, XP, YP
C
C Local parameters:
C
C DELX,DELY =   PX-XMIN, PY-YMIN
C FIRST =       Logical variable with value TRUE iff the
C                 first unmarked node has yet to be
C                 encountered
C I,J =         Cell indexes in the range [I1,I2] X [J1,J2]
C I0,J0 =       Indexes of the cell containing or closest
C                 to P
C I1,I2,J1,J2 = Range of cell indexes defining the layer
C                 whose intersection with the range
C                 [IMIN,IMAX] X [JMIN,JMAX] is currently
C                 being searched
C IMIN,IMAX =   Cell row indexes defining the range of the
C                 search
C JMIN,JMAX =   Cell column indexes defining the range of
C                 the search
C L,LN =        Indexes of nodes in cell (I,J)
C LMIN =        Current candidate for NP
C R =           Distance from P to node LMIN
C RSMIN =       Squared distance from P to node LMIN
C RSQ =         Squared distance from P to node L
C XP,YP =       Local copy of PX,PY -- coordinates of P
C
      XP = PX
      YP = PY
C
C Test for invalid input parameters.
C
      IF (NR .LT. 1  .OR.  DX .LE. 0.  .OR.  DY .LE. 0.)
     .  GO TO 9
C
C Initialize parameters.
C
      FIRST = .TRUE.
      IMIN = 1
      IMAX = NR
      JMIN = 1
      JMAX = NR
      DELX = XP - XMIN
      DELY = YP - YMIN
      I0 = INT(DELX/DX) + 1
      IF (I0 .LT. 1) I0 = 1
      IF (I0 .GT. NR) I0 = NR
      J0 = INT(DELY/DY) + 1
      IF (J0 .LT. 1) J0 = 1
      IF (J0 .GT. NR) J0 = NR
      I1 = I0
      I2 = I0
      J1 = J0
      J2 = J0
C
C Outer loop on layers, inner loop on layer cells, excluding
C   those outside the range [IMIN,IMAX] X [JMIN,JMAX].
C
    1 DO 6 J = J1,J2
        IF (J .GT. JMAX) GO TO 7
        IF (J .LT. JMIN) GO TO 6
        DO 5 I = I1,I2
          IF (I .GT. IMAX) GO TO 6
          IF (I .LT. IMIN) GO TO 5
          IF (J .NE. J1  .AND.  J .NE. J2  .AND.  I .NE. I1
     .        .AND.  I .NE. I2) GO TO 5
C
C Search cell (I,J) for unmarked nodes L.
C
          L = LCELL(I,J)
          IF (L .EQ. 0) GO TO 5
C
C   Loop on nodes in cell (I,J).
C
    2     LN = LNEXT(L)
          IF (LN .LT. 0) GO TO 4
C
C   Node L is not marked.
C
          RSQ = (X(L)-XP)**2 + (Y(L)-YP)**2
          IF (.NOT. FIRST) GO TO 3
C
C   Node L is the first unmarked neighbor of P encountered.
C     Initialize LMIN to the current candidate for NP, and
C     RSMIN to the squared distance from P to LMIN.  IMIN,
C     IMAX, JMIN, and JMAX are updated to define the smal-
C     lest rectangle containing a circle of radius R =
C     Sqrt(RSMIN) centered at P, and contained in [1,NR] X
C     [1,NR] (except that, if P is outside the rectangle
C     defined by the nodes, it is possible that IMIN > NR,
C     IMAX < 1, JMIN > NR, or JMAX < 1).  FIRST is reset to
C     FALSE.
C
          LMIN = L
          RSMIN = RSQ
          R = SQRT(RSMIN)
          IMIN = INT((DELX-R)/DX) + 1
          IF (IMIN .LT. 1) IMIN = 1
          IMAX = INT((DELX+R)/DX) + 1
          IF (IMAX .GT. NR) IMAX = NR
          JMIN = INT((DELY-R)/DY) + 1
          IF (JMIN .LT. 1) JMIN = 1
          JMAX = INT((DELY+R)/DY) + 1
          IF (JMAX .GT. NR) JMAX = NR
          FIRST = .FALSE.
          GO TO 4
C
C   Test for node L closer than LMIN to P.
C
    3     IF (RSQ .GE. RSMIN) GO TO 4
C
C   Update LMIN and RSMIN.
C
          LMIN = L
          RSMIN = RSQ
C
C   Test for termination of loop on nodes in cell (I,J).
C
    4     IF (ABS(LN) .EQ. L) GO TO 5
          L = ABS(LN)
          GO TO 2
    5     CONTINUE
    6   CONTINUE
C
C Test for termination of loop on cell layers.
C
    7 IF (I1 .LE. IMIN  .AND.  I2 .GE. IMAX  .AND.
     .    J1 .LE. JMIN  .AND.  J2 .GE. JMAX) GO TO 8
      I1 = I1 - 1
      I2 = I2 + 1
      J1 = J1 - 1
      J2 = J2 + 1
      GO TO 1
C
C Unless no unmarked nodes were encountered, LMIN is the
C   closest unmarked node to P.
C
    8 IF (FIRST) GO TO 9
      NP = LMIN
      DSQ = RSMIN
      LNEXT(LMIN) = -LNEXT(LMIN)
      RETURN
C
C Error:  NR, DX, or DY is invalid or all nodes are marked.
C
    9 NP = 0
      DSQ = 0.
      RETURN
      END
      SUBROUTINE GIVENS ( A,B, C,S)
      DOUBLE PRECISION A, B, C, S
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This subroutine constructs the Givens plane rotation,
C
C           ( C  S)
C       G = (     ) , where C*C + S*S = 1,
C           (-S  C)
C
C which zeros the second component of the vector (A,B)**T
C (transposed).  Subroutine ROTATE may be called to apply
C the transformation to a 2 by N matrix.
C
C   This routine is identical to subroutine SROTG from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       A,B = Components of the vector defining the rota-
C             tion.  These are overwritten by values R
C             and Z (described below) which define C and S.
C
C On output:
C
C       A = Signed Euclidean norm R of the input vector:
C           R = +/-SQRT(A*A + B*B)
C
C       B = Value Z such that:
C             C = SQRT(1-Z*Z) and S=Z if ABS(Z) .LE. 1, and
C             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.
C
C       C = +/-(A/R) or 1 if R = 0.
C
C       S = +/-(B/R) or 0 if R = 0.
C
C Modules required by GIVENS:  None
C
C Intrinsic functions called by GIVENS:  ABS, SQRT
C
C***********************************************************
C
      DOUBLE PRECISION AA, BB, R, U, V
C
C Local parameters:
C
C AA,BB = Local copies of A and B
C R =     C*A + S*B = +/-SQRT(A*A+B*B)
C U,V =   Variables used to scale A and B for computing R
C
      AA = A
      BB = B
      IF (ABS(AA) .LE. ABS(BB)) GO TO 1
C
C ABS(A) > ABS(B).
C
      U = AA + AA
      V = BB/U
      R = SQRT(.25 + V*V) * U
      C = AA/R
      S = V * (C + C)
C
C Note that R has the sign of A, C > 0, and S has
C   SIGN(A)*SIGN(B).
C
      B = S
      A = R
      RETURN
C
C ABS(A) .LE. ABS(B).
C
    1 IF (BB .EQ. 0.) GO TO 2
      U = BB + BB
      V = AA/U
C
C Store R in A.
C
      A = SQRT(.25 + V*V) * U
      S = BB/A
      C = V * (S + S)
C
C Note that R has the sign of B, S > 0, and C has
C   SIGN(A)*SIGN(B).
C
      B = 1.
      IF (C .NE. 0.) B = 1./C
      RETURN
C
C A = B = 0.
C
    2 C = 1.
      S = 0.
      RETURN
      END
      SUBROUTINE ROTATE (N,C,S, X,Y )
      INTEGER N
      DOUBLE PRECISION C, S, X(N), Y(N)
C
C***********************************************************
C
C                                               From SRFPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C                                                ( C  S)
C   This subroutine applies the Givens rotation  (     )  to
C                                                (-S  C)
C                    (X(1) ... X(N))
C the 2 by N matrix  (             ) .
C                    (Y(1) ... Y(N))
C
C   This routine is identical to subroutine SROT from the
C LINPACK BLAS (Basic Linear Algebra Subroutines).
C
C On input:
C
C       N = Number of columns to be rotated.
C
C       C,S = Elements of the Givens rotation.  Refer to
C             subroutine GIVENS.
C
C The above parameters are not altered by this routine.
C
C       X,Y = Arrays of length .GE. N containing the compo-
C             nents of the vectors to be rotated.
C
C On output:
C
C       X,Y = Arrays containing the rotated vectors (not
C             altered if N < 1).
C
C Modules required by ROTATE:  None
C
C***********************************************************
C
      INTEGER I
      DOUBLE PRECISION XI, YI
C
      DO 1 I = 1,N
        XI = X(I)
        YI = Y(I)
        X(I) = C*XI + S*YI
        Y(I) = -S*XI + C*YI
    1   CONTINUE
      RETURN
      END
      SUBROUTINE SETUPC (XI,YI,ZI,W, ROW)
      DOUBLE PRECISION XI, YI, ZI, W, ROW(11)
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   02/19/97
C
C   This subroutine sets up the I-th row of an augmented re-
C gression matrix for a weighted least squares fit of a
C cosine series function f(x,y) to a set of data values z.
C
C On input:
C
C       XI,YI = Normalized nodal coordinates (in the range
C               0 to Pi).
C
C       ZI = Data value at node I.
C
C       W = Weight associated with node I.
C
C The above parameters are not altered by this routine.
C
C       ROW = Array of length 11.
C
C On output:
C
C       ROW = Array containing a row of the augmented re-
C             gression matrix.
C
C Modules required by SETUPC:  None
C
C Intrinsic function called by SETUPC:  COS
C
C***********************************************************
C
      DOUBLE PRECISION CX1, CX2, CX3, CY1, CY2, CY3
C
C Local parameters:
C
C CX1 = Cos(XI)
C CX2 = Cos(2*XI)
C CX3 = Cos(3*XI)
C CY1 = Cos(YI)
C CY2 = Cos(2*YI)
C CY3 = Cos(3*YI)
C
      CX1 = COS(XI)
      CY1 = COS(YI)
      CX2 = COS(2.0*XI)
      CY2 = COS(2.0*YI)
      CX3 = COS(3.0*XI)
      CY3 = COS(3.0*YI)
C
      ROW(1) = W
      ROW(2) = W*CX1
      ROW(3) = W*CY1
      ROW(4) = W*CX2
      ROW(5) = W*CX1*CY1
      ROW(6) = W*CY2
      ROW(7) = W*CX3
      ROW(8) = W*CX2*CY1
      ROW(9) = W*CX1*CY2
      ROW(10) = W*CY3
      ROW(11) = W*ZI
      RETURN
      END
      SUBROUTINE STORE2 (N,X,Y,NR, LCELL,LNEXT,XMIN,YMIN,DX,
     .                   DY,IER)
      INTEGER N, NR, LCELL(NR,NR), LNEXT(N), IER
      DOUBLE PRECISION X(N), Y(N), XMIN, YMIN, DX, DY
C
C***********************************************************
C
C                                               From TSHEP2D
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   03/28/97
C
C   Given a set of N arbitrarily distributed nodes in the
C plane, this subroutine creates a data structure for a
C cell-based method of solving closest-point problems.  The
C smallest rectangle containing the nodes is partitioned
C into an NR by NR uniform grid of cells, and nodes are as-
C sociated with cells.  In particular, the data structure
C stores the indexes of the nodes contained in each cell.
C For a uniform random distribution of nodes, the nearest
C node to an arbitrary point can be determined in constant
C expected time.
C
C   This code is essentially unaltered from the subroutine
C of the same name in QSHEP2D.
C
C On input:
C
C       N = Number of nodes.  N .GE. 2.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the nodes.
C
C       NR = Number of rows and columns in the grid.  The
C            cell density (average number of nodes per cell)
C            is D = N/(NR**2).  A recommended value, based
C            on empirical evidence, is D = 3 -- NR =
C            Sqrt(N/3).  NR .GE. 1.
C
C The above parameters are not altered by this routine.
C
C       LCELL = Array of length .GE. NR**2.
C
C       LNEXT = Array of length .GE. N.
C
C On output:
C
C       LCELL = NR by NR cell array such that LCELL(I,J)
C               contains the index (for X and Y) of the
C               first node (node with smallest index) in
C               cell (I,J), or LCELL(I,J) = 0 if no nodes
C               are contained in the cell.  The upper right
C               corner of cell (I,J) has coordinates (XMIN+
C               I*DX,YMIN+J*DY).  LCELL is not defined if
C               IER .NE. 0.
C
C       LNEXT = Array of next-node indexes such that
C               LNEXT(K) contains the index of the next node
C               in the cell which contains node K, or
C               LNEXT(K) = K if K is the last node in the
C               cell for K = 1,...,N.  (The nodes contained
C               in a cell are ordered by their indexes.)
C               If, for example, cell (I,J) contains nodes
C               2, 3, and 5 (and no others), then LCELL(I,J)
C               = 2, LNEXT(2) = 3, LNEXT(3) = 5, and
C               LNEXT(5) = 5.  LNEXT is not defined if
C               IER .NE. 0.
C
C       XMIN,YMIN = Cartesian coordinates of the lower left
C                   corner of the rectangle defined by the
C                   nodes (smallest nodal coordinates) un-
C                   less IER = 1.  The upper right corner is
C                   (XMAX,YMAX) for XMAX = XMIN + NR*DX and
C                   YMAX = YMIN + NR*DY.
C
C       DX,DY = Dimensions of the cells unless IER = 1.  DX
C               = (XMAX-XMIN)/NR and DY = (YMAX-YMIN)/NR,
C               where XMIN, XMAX, YMIN, and YMAX are the
C               extrema of X and Y.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 2 or NR < 1.
C             IER = 2 if DX = 0 or DY = 0.
C
C Modules required by STORE2:  None
C
C Intrinsic functions called by STORE2:  DBLE, INT
C
C***********************************************************
C
      INTEGER I, J, K, L, NN, NNR
      DOUBLE PRECISION DELX, DELY, XMN, XMX, YMN, YMX
C
C Local parameters:
C
C DELX,DELY = Components of the cell dimensions -- local
C               copies of DX,DY
C I,J =       Cell indexes
C K =         Nodal index
C L =         Index of a node in cell (I,J)
C NN =        Local copy of N
C NNR =       Local copy of NR
C XMN,XMX =   Range of nodal X coordinates
C YMN,YMX =   Range of nodal Y coordinates
C
      NN = N
      NNR = NR
      IF (NN .LT. 2  .OR.  NNR .LT. 1) GO TO 5
C
C Compute the dimensions of the rectangle containing the
C   nodes.
C
      XMN = X(1)
      XMX = XMN
      YMN = Y(1)
      YMX = YMN
      DO 1 K = 2,NN
        IF (X(K) .LT. XMN) XMN = X(K)
        IF (X(K) .GT. XMX) XMX = X(K)
        IF (Y(K) .LT. YMN) YMN = Y(K)
        IF (Y(K) .GT. YMX) YMX = Y(K)
    1   CONTINUE
      XMIN = XMN
      YMIN = YMN
C
C Compute cell dimensions and test for zero area.
C
      DELX = (XMX-XMN)/DBLE(NNR)
      DELY = (YMX-YMN)/DBLE(NNR)
      DX = DELX
      DY = DELY
      IF (DELX .EQ. 0.  .OR.  DELY .EQ. 0.) GO TO 6
C
C Initialize LCELL.
C
      DO 3 J = 1,NNR
        DO 2 I = 1,NNR
          LCELL(I,J) = 0
    2     CONTINUE
    3   CONTINUE
C
C Loop on nodes, storing indexes in LCELL and LNEXT.
C
      DO 4 K = NN,1,-1
        I = INT((X(K)-XMN)/DELX) + 1
        IF (I .GT. NNR) I = NNR
        J = INT((Y(K)-YMN)/DELY) + 1
        IF (J .GT. NNR) J = NNR
        L = LCELL(I,J)
        LNEXT(K) = L
        IF (L .EQ. 0) LNEXT(K) = K
        LCELL(I,J) = K
    4   CONTINUE
C
C No errors encountered.
C
      IER = 0
      RETURN
C
C Invalid input parameter.
C
    5 IER = 1
      RETURN
C
C DX = 0 or DY = 0.
C
    6 IER = 2
      RETURN
      END
