C----------------------------------------------------------------
      SUBROUTINE DMTOMS(DM,N,MST,IMST,CST)
C
C     THIS SUBROUTINE FINDS A SET OF EDGES OF A LINEAR GRAPH
C     COMPRISING A TREE WITH MINIMAL TOTAL EDGE LENGTH. THE
C     GRAPH IS SPECIFIED AS AN ARRAY OF INTER-NODE EDGE LENGTHS.
C     THE EDGES OF THE MINIMAL SPANNING TREE OF THE GRAPH ARE
C     PLACED IN ARRAY MST.  EXECUTION TIME IS PROPORTIONAL TO
C     THE SQUARE OF THE NUMBER OF NODES.
C
C
C ADDITIONAL NOTES:
C -----------------
C
C This code is taken from "Algorithm 422 Minimal Spanning Tree [H]"
C by V. Kevin M. Whitney, Communications of the ACM, April 1972,
C Volume 15, Number 4, pages 273-274.  It has been digitized on
C September 16, 1999 so it can be added to the CALGO collection of
C algorithms made available by the ACM.  The code has been modified
C slightly so it could be compiled and tested using Fortran 77.
C
C The algorithm uses a technique suggested by E.W. Dijkstra in
C "A note on two problems in connection with graphs",
C Numer. Math. 1, 5 (Oct. 1959), 269-271.
C
C from Whitney:
C "The Dijkstra algorithm grows a minimal spanning tree by successively
C adjoining the nearest remaining node to a partially formed tree until
C all nodes of the graph are included in the tree.  At each iterative
C step the nodes not yet included in the tree are stored in array NIT.
C The node of the partially completed tree nearest to node NIT(I) is
C stored in JI(I), and the length of edge from NIT(I) to JI(I) is stored
C in UI(I). Hence the node not yet in the tree which is nearest to a
C node of the tree may be found by searching for the minimal element of
C array UI. That node, KP, is added to the tree and removed from array
C NIT. For each node remaining in array NIT, the distance from the
C nearest node of the tree (stored in array UI) is compared to the
C distance from KP, the new node of the tree, and arrays UI and JI are
C updated if the new distance is shorter.  The nearest node selection
C and list updating are performed N - 1 times until all nodes are in
C the tree."
C
C Diagonal elements of array DM are not used.  The edges of the output
C minimal spanning tree are specified by pairs of nodes in array MST.
C If the graph represented by the inter-node edge length array DM is
C not connected, the procedure will generate a minimal spanning forest
C containing the minimal spanning trees of the disjoint components
C joined together by edges of length 10.**10.  A disconnected graph
C is indicated by a value 10.**10 for variable UK at step 500 during
C execution of the algorithm.  This algorithm can also be used to find
C a maximal spanning tree by changing the loop between statements 300
C and 400 to search for the most distant rather than for the nearest
C remaining node to be adjoined to the partially completed tree.
C
C
C     CALLING SEQUENCE VARIABLES ARE:
C
C     DM    ARRAY OF INTER-NODE EDGE LENGTHS.
C           DM(I,J) (1 .LE. I,J .LE. IN) IS THE LENGTH OF
C           AN EDGE FROM NODE I TO NODE J. IF THERE IS NO
C           EDGE FROM NODE I TO NODE J, SET DM(I,J)=10.**10
C     N     NODES ARE NUMBERED  1, 2, ..., N.
C
C     MST   ARRAY IN WHICH EDGE LIST OF MST IS PLACED. MST(1,I)
C           IS THE ORIGINAL NODE AND MST(2,I) IS THE TERMINAL
C           NODE OF EDGE I FOR 1 .LE. I .LE. IMST.
C     IMST  NUMBER OF EDGES IN ARRAY MST.
C     CST   SUM OF EDGE LENGTHS OF EDGES OF TREE.
C
C     PROGRAM VARIABLES    :
C
C     NIT   ARRAY OF NODES NOT YET IN TREE.
C     NITP  NUMBER OF NODES IN ARRAY NIT.
C     JI(I) NODE OF PARTIAL MST CLOSEST TO NODE NIT(I).
C     UI(I) LENGTH OF EDGE FROM NIT(I) TO JI(I).
C     KP    NEXT NODE TO BE ADDED TO ARRAY MST.
C
C
C     INITIALIZE NODE LABEL ARRAYS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION CST
      INTEGER IMST,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DM(NMAX,NMAX)
      INTEGER MST(2,NMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION D,UK
      INTEGER I,K,KP,NI,NITP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION UI(NMAX)
      INTEGER JI(NMAX),NIT(NMAX)
C     ..
C     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=50)
C     ..
      CST = 0.0D0
      NITP = N - 1
      KP = N
      IMST = 0
      DO 10 I = 1,NITP
          NIT(I) = I
          UI(I) = DM(I,KP)
          JI(I) = KP
   10 CONTINUE
C
C     UPDATE LABELS OF NODES NOT YET IN TREE.
C
   20 DO 30 I = 1,NITP
          NI = NIT(I)
          D = DM(NI,KP)
          IF (UI(I).LE.D) GO TO 30
          UI(I) = D
          JI(I) = KP
   30 CONTINUE
C
C     FIND NODE OUTSIDE TREE NEAREST TO TREE.
C
      UK = UI(1)
      DO 40 I = 1,NITP
          IF (UI(I).GT.UK) GO TO 40
          UK = UI(I)
          K = I
   40 CONTINUE
C
C     PUT NODES OF APPROPRIATE EDGE INTO ARRAY MST.
C
      IMST = IMST + 1
      MST(1,IMST) = NIT(K)
      MST(2,IMST) = JI(K)
      CST = CST + UK
      KP = NIT(K)
C
C     DELETE NEW TREE NODE FROM ARRAY NIT.
C
      UI(K) = UI(NITP)
      NIT(K) = NIT(NITP)
      JI(K) = JI(NITP)
      NITP = NITP - 1
   50 IF (NITP.NE.0) GO TO 20
C
C     WHEN ALL NODES ARE IN TREE, QUIT.
C
      RETURN

      END
