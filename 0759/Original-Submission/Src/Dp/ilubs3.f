      SUBROUTINE ILU (NPTS, NPD, A, LLDG, LSL, LLSL)
      INTEGER NPDE
      PARAMETER (NPDE = 3)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPD, LLDG(NPTS,-9:-2), LSL(*), LLSL(0:*)
      DOUBLE PRECISION A(NPTS,NPDE,NPDE,-9:9)
C
Ccc PURPOSE:
C Incomplete LU decomposition of A
C    A stems from a 19-point stencil on a grid with irregular row and
C    plane sizes
C A = ILU, vectorized by `hyperplane' ordening, where the `hyperplanes'
C    are sets of points that can be treated simultaneously.
C A((i,j,k),1:NPDE,1:NPDE,.) contains a block of NPDE.NPDE elements
C    corresponding with node (i,j,k)
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : INOUT.
C          IN: A(.,1:NPDE,1:NPDE,-9:-1)  : lower block diagonals
C              A(.,1:NPDE,1:NPDE,0)      :   main block diagonal
C              A(.,1:NPDE,1:NPDE, 1: 9)  : upper block diagonals
C          OUT: A(.,1:NPDE,1:NPDE,-9:-1)  : lower block diagonals of L
C               A(.,ic,jc,0):      jc < ic: main block diagonal of L
C                                           main diagonal L == I
C                                  jc >=ic: main block diagonal of U
C                                           main diagonal U inverted
C               A(.,1:NPDE,1:NPDE, 1: 9)  : upper block diagonals of U
C LLDG   : IN. Block-column index of lower 8 block-diagonals
C              If block ld does not exist the LLDG(N,ld) = N
C LSL    : IN. (NPTS)
C          LSL(LLSL(m-1)+1:LLSL(m)): pointers to set of points S_m
C          in actual grid that can be treated at the m-th iteration
C LLSL   : IN. (0:LLSL(0))
C          LLSL(0) = # iterations needed
C          LLSL(1:LLSL(0)): pointers to the start of a list in LSL
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, JC, LC, N, L, M
C
      IF (NPDE .NE. NPD) STOP 'Wrong ILUBS loaded.'
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # N = S_m(LLSL(l))
C
C S_1 = {(i,j,k)| (i,j,k) not dependent on (i-ii,j-jj,k-kk),
C                                   ii,jj,kk >= 0, not ii=jj=kk=0;
C                 i.e., Dirichlet points and left/down/front corners}
      M = 1
C
C   Compute main block diagonal
         DO 550 IC = 1, NPDE
            DO 554 LC = 1, IC-1
            DO 555 JC = IC, NPDE
CDIR$ IVDEP
            DO 551 L = 1, LLSL(M)
               N = LSL(L)
               A(N,IC,JC,0) = A(N,IC,JC,0)
     +                        - A(N,IC,LC,0)*A(N,LC,JC,0)
  551       CONTINUE
  555       CONTINUE
            DO 556 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 552 L = 1, LLSL(M)
               N = LSL(L)
               A(N,JC,IC,0) = A(N,JC,IC,0)
     +                        - A(N,JC,LC,0)*A(N,LC,IC,0)
  552       CONTINUE
  556       CONTINUE
  554       CONTINUE
CDIR$ IVDEP
            DO 553 L = 1, LLSL(M)
               N = LSL(L)
               A(N,IC,IC,0) = 1.0 / A(N,IC,IC,0)
  553       CONTINUE
            DO 557 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 559 L = 1, LLSL(M)
               N = LSL(L)
               A(N,JC,IC,0) = A(N,JC,IC,0) * A(N,IC,IC,0)
  559       CONTINUE
  557       CONTINUE
  550    CONTINUE
C
C   Compute upper block diagonals
         DO 560 IC = 1, NPDE
            DO 563 LC = 1, IC-1
            DO 564 JC = 1, NPDE
CDIR$ IVDEP
            DO 561 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,IC,JC,1) = A(N,IC,JC,1)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,1)
               A(N,IC,JC,2) = A(N,IC,JC,2)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,2)
               A(N,IC,JC,3) = A(N,IC,JC,3)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,3)
               A(N,IC,JC,4) = A(N,IC,JC,4)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,4)
               A(N,IC,JC,5) = A(N,IC,JC,5)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,5)
               A(N,IC,JC,6) = A(N,IC,JC,6)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,6)
               A(N,IC,JC,7) = A(N,IC,JC,7)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,7)
               A(N,IC,JC,8) = A(N,IC,JC,8)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,8)
               A(N,IC,JC,9) = A(N,IC,JC,9)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,9)
  561       CONTINUE
  564       CONTINUE
  563       CONTINUE
  560    CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the first row and the first
C     point of the second row, since N > 1 in the loop and for those
C     points LLDG(N,.) = N, (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 10 M = 2, LLSL(0)
C
C   Compute lower diagonals
         DO 120 JC = 1, NPDE
            DO 121 LC = 1, JC-1
CDIR$ IVDEP
            DO 20 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 122 IC = 1, NPDE
               A(N,IC,JC,-9) = A(N,IC,JC,-9)
     +                       - A(N,IC,LC,-9) * A(LLDG(N,-9),LC,JC,0)
  122       CONTINUE
   20       CONTINUE
  121       CONTINUE
CDIR$ IVDEP
            DO 21 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 123 IC = 1, NPDE
               A(N,IC,JC,-9) = A(N,IC,JC,-9) * A(LLDG(N,-9),JC,JC,0)
  123       CONTINUE
   21       CONTINUE
  120    CONTINUE
         DO 130 JC = 1, NPDE
CDIR$ IVDEP
            DO 30 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 131 LC = 1, NPDE
CFPP$ UNROLL
            DO 132 IC = 1, NPDE
               A(N,IC,JC,-8) = A(N,IC,JC,-8)
     +                       - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,2)
  132       CONTINUE
  131       CONTINUE
   30       CONTINUE
            DO 133 LC = 1, JC-1
CDIR$ IVDEP
            DO 31 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 134 IC = 1, NPDE
               A(N,IC,JC,-8) = A(N,IC,JC,-8)
     +                       - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,0)
  134       CONTINUE
   31       CONTINUE
  133       CONTINUE
CDIR$ IVDEP
            DO 32 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 135 IC = 1, NPDE
               A(N,IC,JC,-8) = A(N,IC,JC,-8) * A(LLDG(N,-8),JC,JC,0)
  135       CONTINUE
   32       CONTINUE
  130    CONTINUE
         DO 140 JC = 1, NPDE
CDIR$ IVDEP
            DO 40 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 141 LC = 1, NPDE
CFPP$ UNROLL
            DO 142 IC = 1, NPDE
               A(N,IC,JC,-7) = A(N,IC,JC,-7)
     +                       - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,3)
     +                       - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,1)
  142       CONTINUE
  141       CONTINUE
   40       CONTINUE
            DO 143 LC = 1, JC-1
CDIR$ IVDEP
            DO 41 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 144 IC = 1, NPDE
               A(N,IC,JC,-7) = A(N,IC,JC,-7)
     +                       - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,0)
  144       CONTINUE
   41       CONTINUE
  143       CONTINUE
CDIR$ IVDEP
            DO 42 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 145 IC = 1, NPDE
               A(N,IC,JC,-7) = A(N,IC,JC,-7) * A(LLDG(N,-7),JC,JC,0)
  145       CONTINUE
   42       CONTINUE
  140    CONTINUE
         DO 150 JC = 1, NPDE
CDIR$ IVDEP
            DO 50 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 151 LC = 1, NPDE
CFPP$ UNROLL
            DO 152 IC = 1, NPDE
               A(N,IC,JC,-6) = A(N,IC,JC,-6)
     +                       - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,4)
     +                       - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,1)
  152       CONTINUE
  151       CONTINUE
   50       CONTINUE
            DO 153 LC = 1, JC-1
CDIR$ IVDEP
            DO 51 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 154 IC = 1, NPDE
               A(N,IC,JC,-6) = A(N,IC,JC,-6)
     +                       - A(N,IC,LC,-6)*A(LLDG(N,-6),LC,JC,0)
  154       CONTINUE
   51       CONTINUE
  153       CONTINUE
CDIR$ IVDEP
            DO 52 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 155 IC = 1, NPDE
               A(N,IC,JC,-6) = A(N,IC,JC,-6) * A(LLDG(N,-6),JC,JC,0)
  155       CONTINUE
   52       CONTINUE
  150    CONTINUE
         DO 160 JC = 1, NPDE
CDIR$ IVDEP
            DO 60 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 161 LC = 1, NPDE
CFPP$ UNROLL
            DO 162 IC = 1, NPDE
               A(N,IC,JC,-5) = A(N,IC,JC,-5)
     +                 - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,4)
     +                 - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,3)
     +                 - A(N,IC,LC,-6)*A(LLDG(N,-6),LC,JC,2)
               A(N,IC,JC,-4) = A(N,IC,JC,-4)
     +                 - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,6)
     +                 - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,5)
  162       CONTINUE
  161       CONTINUE
   60       CONTINUE
            DO 163 LC = 1, JC-1
CDIR$ IVDEP
            DO 61 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 164 IC = 1, NPDE
               A(N,IC,JC,-5) = A(N,IC,JC,-5)
     +                       - A(N,IC,LC,-5)*A(LLDG(N,-5),LC,JC,0)
               A(N,IC,JC,-4) = A(N,IC,JC,-4)
     +                       - A(N,IC,LC,-4)*A(LLDG(N,-4),LC,JC,0)
  164       CONTINUE
   61       CONTINUE
  163       CONTINUE
CDIR$ IVDEP
            DO 62 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 165 IC = 1, NPDE
               A(N,IC,JC,-5) = A(N,IC,JC,-5) * A(LLDG(N,-5),JC,JC,0)
               A(N,IC,JC,-4) = A(N,IC,JC,-4) * A(LLDG(N,-4),JC,JC,0)
  165       CONTINUE
   62       CONTINUE
  160    CONTINUE
         DO 170 JC = 1, NPDE
CDIR$ IVDEP
            DO 70 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 171 LC = 1, NPDE
CFPP$ UNROLL
            DO 172 IC = 1, NPDE
               A(N,IC,JC,-3) = A(N,IC,JC,-3)
     +                 - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,8)
     +                 - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,5)
     +                 - A(N,IC,LC,-4)*A(LLDG(N,-4),LC,JC,1)
  172       CONTINUE
  171       CONTINUE
   70       CONTINUE
            DO 173 LC = 1, JC-1
CDIR$ IVDEP
            DO 71 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 174 IC = 1, NPDE
               A(N,IC,JC,-3) = A(N,IC,JC,-3)
     +                       - A(N,IC,LC,-3)*A(LLDG(N,-3),LC,JC,0)
  174       CONTINUE
   71       CONTINUE
  173       CONTINUE
CDIR$ IVDEP
            DO 72 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 175 IC = 1, NPDE
               A(N,IC,JC,-3) = A(N,IC,JC,-3) * A(LLDG(N,-3),JC,JC,0)
  175       CONTINUE
   72       CONTINUE
  170    CONTINUE
         DO 180 JC = 1, NPDE
CDIR$ IVDEP
            DO 80 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 181 LC = 1, NPDE
CFPP$ UNROLL
            DO 182 IC = 1, NPDE
               A(N,IC,JC,-2) = A(N,IC,JC,-2)
     +                 - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,8)
     +                 - A(N,IC,LC,-6)*A(LLDG(N,-6),LC,JC,5)
     +                 - A(N,IC,LC,-3)*A(LLDG(N,-3),LC,JC,1)
               A(N,IC,JC,-1) = A(N,IC,JC,-1)
     +                 - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,7)
     +                 - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,6)
     +                 - A(N,IC,LC,-4)*A(LLDG(N,-4),LC,JC,3)
     +                 - A(N,IC,LC,-3)*A(LLDG(N,-3),LC,JC,2)
  182       CONTINUE
  181       CONTINUE
   80       CONTINUE
            DO 183 LC = 1, JC-1
CDIR$ IVDEP
            DO 81 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 184 IC = 1, NPDE
               A(N,IC,JC,-2) = A(N,IC,JC,-2)
     +                       - A(N,IC,LC,-2)*A(LLDG(N,-2),LC,JC,0)
               A(N,IC,JC,-1) = A(N,IC,JC,-1)
     +                       - A(N,IC,LC,-1)*A(N-1       ,LC,JC,0)
  184       CONTINUE
   81       CONTINUE
  183       CONTINUE
CDIR$ IVDEP
            DO 82 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 185 IC = 1, NPDE
               A(N,IC,JC,-2) = A(N,IC,JC,-2) * A(LLDG(N,-2),JC,JC,0)
               A(N,IC,JC,-1) = A(N,IC,JC,-1) * A(N-1       ,JC,JC,0)
  185       CONTINUE
   82       CONTINUE
  180    CONTINUE
C
C   Compute main diagonal
         DO 300 IC = 1, NPDE
            DO 302 JC = IC, NPDE
CDIR$ IVDEP
            DO 200 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 301 LC = 1, NPDE
               A(N,IC,JC,0) = A(N,IC,JC, 0)
     +                      - A(N,IC,LC,-9)*A(LLDG(N,-9),LC,JC,9)
     +                      - A(N,IC,LC,-8)*A(LLDG(N,-8),LC,JC,8)
     +                      - A(N,IC,LC,-7)*A(LLDG(N,-7),LC,JC,7)
     +                      - A(N,IC,LC,-6)*A(LLDG(N,-6),LC,JC,6)
     +                      - A(N,IC,LC,-5)*A(LLDG(N,-5),LC,JC,5)
     +                      - A(N,IC,LC,-4)*A(LLDG(N,-4),LC,JC,4)
     +                      - A(N,IC,LC,-3)*A(LLDG(N,-3),LC,JC,3)
     +                      - A(N,IC,LC,-2)*A(LLDG(N,-2),LC,JC,2)
     +                      - A(N,IC,LC,-1)*A(N-1       ,LC,JC,1)
  301       CONTINUE
  200       CONTINUE
  302       CONTINUE
            DO 303 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 201 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 1301 LC = 1, NPDE
               A(N,JC,IC,0) = A(N,JC,IC, 0)
     +                      - A(N,JC,LC,-9)*A(LLDG(N,-9),LC,IC,9)
     +                      - A(N,JC,LC,-8)*A(LLDG(N,-8),LC,IC,8)
     +                      - A(N,JC,LC,-7)*A(LLDG(N,-7),LC,IC,7)
     +                      - A(N,JC,LC,-6)*A(LLDG(N,-6),LC,IC,6)
     +                      - A(N,JC,LC,-5)*A(LLDG(N,-5),LC,IC,5)
     +                      - A(N,JC,LC,-4)*A(LLDG(N,-4),LC,IC,4)
     +                      - A(N,JC,LC,-3)*A(LLDG(N,-3),LC,IC,3)
     +                      - A(N,JC,LC,-2)*A(LLDG(N,-2),LC,IC,2)
     +                      - A(N,JC,LC,-1)*A(N-1       ,LC,IC,1)
 1301       CONTINUE
  201       CONTINUE
  303       CONTINUE
            DO 304 LC = 1, IC-1
            DO 305 JC = IC, NPDE
CDIR$ IVDEP
            DO 202 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,IC,JC,0) = A(N,IC,JC,0)
     +                      - A(N,IC,LC,0)*A(N        ,LC,JC,0)
  202       CONTINUE
  305       CONTINUE
            DO 306 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 203 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,JC,IC,0) = A(N,JC,IC,0)
     +                      - A(N,JC,LC,0)*A(N        ,LC,IC,0)
  203       CONTINUE
  306       CONTINUE
  304       CONTINUE
CDIR$ IVDEP
            DO 204 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,IC,IC,0) = 1.0 / A(N,IC,IC,0)
  204       CONTINUE
            DO 307 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 205 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,JC,IC,0) = A(N,JC,IC,0) * A(N,IC,IC,0)
  205       CONTINUE
  307       CONTINUE
  300    CONTINUE
C
C   Compute upper diagonals
         DO 500 IC = 1, NPDE
CDIR$ IVDEP
            DO 400 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
      A(N,IC,1,1) = A(N,IC,1,1) - A(N,IC,1,-7)*A(LLDG(N,-7),1,1,8) - A
     1   (N,IC,1,-6)*A(LLDG(N,-6),1,1,7) - A(N,IC,1,-3)*A(LLDG(N,-3),
     2   1,1,4) - A(N,IC,1,-2)*A(LLDG(N,-2),1,1,3)
      A(N,IC,1,2) = A(N,IC,1,2) - A(N,IC,1,-8)*A(LLDG(N,-8),1,1,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,1,6) - A(N,IC,1,-1)*A(N-1,1,1,3)
      A(N,IC,1,3) = A(N,IC,1,3) - A(N,IC,1,-7)*A(LLDG(N,-7),1,1,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,1,7) - A(N,IC,1,-1)*A(N-1,1,1,4)
      A(N,IC,1,4) = A(N,IC,1,4) - A(N,IC,1,-6)*A(LLDG(N,-6),1,1,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,1,8)
      A(N,IC,1,5) = A(N,IC,1,5) - A(N,IC,1,-9)*A(LLDG(N,-9),1,1,8) - A
     1   (N,IC,1,-3)*A(LLDG(N,-3),1,1,7) - A(N,IC,1,-2)*A(LLDG(N,-2),
     2   1,1,6)
      A(N,IC,1,6) = A(N,IC,1,6) - A(N,IC,1,-4)*A(LLDG(N,-4),1,1,9) - A
     1   (N,IC,1,-1)*A(N-1,1,1,7)
      A(N,IC,1,7) = A(N,IC,1,7) - A(N,IC,1,-3)*A(LLDG(N,-3),1,1,9) - A
     1   (N,IC,1,-1)*A(N-1,1,1,8)
      A(N,IC,1,8) = A(N,IC,1,8) - A(N,IC,1,-2)*A(LLDG(N,-2),1,1,9)
      A(N,IC,2,1) = A(N,IC,2,1) - A(N,IC,1,-7)*A(LLDG(N,-7),1,2,8) - A
     1   (N,IC,1,-6)*A(LLDG(N,-6),1,2,7) - A(N,IC,1,-3)*A(LLDG(N,-3),
     2   1,2,4) - A(N,IC,1,-2)*A(LLDG(N,-2),1,2,3)
      A(N,IC,2,2) = A(N,IC,2,2) - A(N,IC,1,-8)*A(LLDG(N,-8),1,2,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,2,6) - A(N,IC,1,-1)*A(N-1,1,2,3)
      A(N,IC,2,3) = A(N,IC,2,3) - A(N,IC,1,-7)*A(LLDG(N,-7),1,2,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,2,7) - A(N,IC,1,-1)*A(N-1,1,2,4)
      A(N,IC,2,4) = A(N,IC,2,4) - A(N,IC,1,-6)*A(LLDG(N,-6),1,2,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,2,8)
      A(N,IC,2,5) = A(N,IC,2,5) - A(N,IC,1,-9)*A(LLDG(N,-9),1,2,8) - A
     1   (N,IC,1,-3)*A(LLDG(N,-3),1,2,7) - A(N,IC,1,-2)*A(LLDG(N,-2),
     2   1,2,6)
      A(N,IC,2,6) = A(N,IC,2,6) - A(N,IC,1,-4)*A(LLDG(N,-4),1,2,9) - A
     1   (N,IC,1,-1)*A(N-1,1,2,7)
      A(N,IC,2,7) = A(N,IC,2,7) - A(N,IC,1,-3)*A(LLDG(N,-3),1,2,9) - A
     1   (N,IC,1,-1)*A(N-1,1,2,8)
      A(N,IC,2,8) = A(N,IC,2,8) - A(N,IC,1,-2)*A(LLDG(N,-2),1,2,9)
      A(N,IC,3,1) = A(N,IC,3,1) - A(N,IC,1,-7)*A(LLDG(N,-7),1,3,8) - A
     1   (N,IC,1,-6)*A(LLDG(N,-6),1,3,7) - A(N,IC,1,-3)*A(LLDG(N,-3),
     2   1,3,4) - A(N,IC,1,-2)*A(LLDG(N,-2),1,3,3)
      A(N,IC,3,2) = A(N,IC,3,2) - A(N,IC,1,-8)*A(LLDG(N,-8),1,3,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,3,6) - A(N,IC,1,-1)*A(N-1,1,3,3)
      A(N,IC,3,3) = A(N,IC,3,3) - A(N,IC,1,-7)*A(LLDG(N,-7),1,3,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,3,7) - A(N,IC,1,-1)*A(N-1,1,3,4)
      A(N,IC,3,4) = A(N,IC,3,4) - A(N,IC,1,-6)*A(LLDG(N,-6),1,3,9) - A
     1   (N,IC,1,-5)*A(LLDG(N,-5),1,3,8)
      A(N,IC,3,5) = A(N,IC,3,5) - A(N,IC,1,-9)*A(LLDG(N,-9),1,3,8) - A
     1   (N,IC,1,-3)*A(LLDG(N,-3),1,3,7) - A(N,IC,1,-2)*A(LLDG(N,-2),
     2   1,3,6)
      A(N,IC,3,6) = A(N,IC,3,6) - A(N,IC,1,-4)*A(LLDG(N,-4),1,3,9) - A
     1   (N,IC,1,-1)*A(N-1,1,3,7)
      A(N,IC,3,7) = A(N,IC,3,7) - A(N,IC,1,-3)*A(LLDG(N,-3),1,3,9) - A
     1   (N,IC,1,-1)*A(N-1,1,3,8)
      A(N,IC,3,8) = A(N,IC,3,8) - A(N,IC,1,-2)*A(LLDG(N,-2),1,3,9)

      A(N,IC,1,1) = A(N,IC,1,1) - A(N,IC,2,-7)*A(LLDG(N,-7),2,1,8) - A
     1   (N,IC,2,-6)*A(LLDG(N,-6),2,1,7) - A(N,IC,2,-3)*A(LLDG(N,-3),
     2   2,1,4) - A(N,IC,2,-2)*A(LLDG(N,-2),2,1,3)
      A(N,IC,1,2) = A(N,IC,1,2) - A(N,IC,2,-8)*A(LLDG(N,-8),2,1,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,1,6) - A(N,IC,2,-1)*A(N-1,2,1,3)
      A(N,IC,1,3) = A(N,IC,1,3) - A(N,IC,2,-7)*A(LLDG(N,-7),2,1,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,1,7) - A(N,IC,2,-1)*A(N-1,2,1,4)
      A(N,IC,1,4) = A(N,IC,1,4) - A(N,IC,2,-6)*A(LLDG(N,-6),2,1,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,1,8)
      A(N,IC,1,5) = A(N,IC,1,5) - A(N,IC,2,-9)*A(LLDG(N,-9),2,1,8) - A
     1   (N,IC,2,-3)*A(LLDG(N,-3),2,1,7) - A(N,IC,2,-2)*A(LLDG(N,-2),
     2   2,1,6)
      A(N,IC,1,6) = A(N,IC,1,6) - A(N,IC,2,-4)*A(LLDG(N,-4),2,1,9) - A
     1   (N,IC,2,-1)*A(N-1,2,1,7)
      A(N,IC,1,7) = A(N,IC,1,7) - A(N,IC,2,-3)*A(LLDG(N,-3),2,1,9) - A
     1   (N,IC,2,-1)*A(N-1,2,1,8)
      A(N,IC,1,8) = A(N,IC,1,8) - A(N,IC,2,-2)*A(LLDG(N,-2),2,1,9)
      A(N,IC,2,1) = A(N,IC,2,1) - A(N,IC,2,-7)*A(LLDG(N,-7),2,2,8) - A
     1   (N,IC,2,-6)*A(LLDG(N,-6),2,2,7) - A(N,IC,2,-3)*A(LLDG(N,-3),
     2   2,2,4) - A(N,IC,2,-2)*A(LLDG(N,-2),2,2,3)
      A(N,IC,2,2) = A(N,IC,2,2) - A(N,IC,2,-8)*A(LLDG(N,-8),2,2,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,2,6) - A(N,IC,2,-1)*A(N-1,2,2,3)
      A(N,IC,2,3) = A(N,IC,2,3) - A(N,IC,2,-7)*A(LLDG(N,-7),2,2,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,2,7) - A(N,IC,2,-1)*A(N-1,2,2,4)
      A(N,IC,2,4) = A(N,IC,2,4) - A(N,IC,2,-6)*A(LLDG(N,-6),2,2,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,2,8)
      A(N,IC,2,5) = A(N,IC,2,5) - A(N,IC,2,-9)*A(LLDG(N,-9),2,2,8) - A
     1   (N,IC,2,-3)*A(LLDG(N,-3),2,2,7) - A(N,IC,2,-2)*A(LLDG(N,-2),
     2   2,2,6)
      A(N,IC,2,6) = A(N,IC,2,6) - A(N,IC,2,-4)*A(LLDG(N,-4),2,2,9) - A
     1   (N,IC,2,-1)*A(N-1,2,2,7)
      A(N,IC,2,7) = A(N,IC,2,7) - A(N,IC,2,-3)*A(LLDG(N,-3),2,2,9) - A
     1   (N,IC,2,-1)*A(N-1,2,2,8)
      A(N,IC,2,8) = A(N,IC,2,8) - A(N,IC,2,-2)*A(LLDG(N,-2),2,2,9)
      A(N,IC,3,1) = A(N,IC,3,1) - A(N,IC,2,-7)*A(LLDG(N,-7),2,3,8) - A
     1   (N,IC,2,-6)*A(LLDG(N,-6),2,3,7) - A(N,IC,2,-3)*A(LLDG(N,-3),
     2   2,3,4) - A(N,IC,2,-2)*A(LLDG(N,-2),2,3,3)
      A(N,IC,3,2) = A(N,IC,3,2) - A(N,IC,2,-8)*A(LLDG(N,-8),2,3,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,3,6) - A(N,IC,2,-1)*A(N-1,2,3,3)
      A(N,IC,3,3) = A(N,IC,3,3) - A(N,IC,2,-7)*A(LLDG(N,-7),2,3,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,3,7) - A(N,IC,2,-1)*A(N-1,2,3,4)
      A(N,IC,3,4) = A(N,IC,3,4) - A(N,IC,2,-6)*A(LLDG(N,-6),2,3,9) - A
     1   (N,IC,2,-5)*A(LLDG(N,-5),2,3,8)
      A(N,IC,3,5) = A(N,IC,3,5) - A(N,IC,2,-9)*A(LLDG(N,-9),2,3,8) - A
     1   (N,IC,2,-3)*A(LLDG(N,-3),2,3,7) - A(N,IC,2,-2)*A(LLDG(N,-2),
     2   2,3,6)
      A(N,IC,3,6) = A(N,IC,3,6) - A(N,IC,2,-4)*A(LLDG(N,-4),2,3,9) - A
     1   (N,IC,2,-1)*A(N-1,2,3,7)
      A(N,IC,3,7) = A(N,IC,3,7) - A(N,IC,2,-3)*A(LLDG(N,-3),2,3,9) - A
     1   (N,IC,2,-1)*A(N-1,2,3,8)
      A(N,IC,3,8) = A(N,IC,3,8) - A(N,IC,2,-2)*A(LLDG(N,-2),2,3,9)

      A(N,IC,1,1) = A(N,IC,1,1) - A(N,IC,3,-7)*A(LLDG(N,-7),3,1,8) - A
     1   (N,IC,3,-6)*A(LLDG(N,-6),3,1,7) - A(N,IC,3,-3)*A(LLDG(N,-3),
     2   3,1,4) - A(N,IC,3,-2)*A(LLDG(N,-2),3,1,3)
      A(N,IC,1,2) = A(N,IC,1,2) - A(N,IC,3,-8)*A(LLDG(N,-8),3,1,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,1,6) - A(N,IC,3,-1)*A(N-1,3,1,3)
      A(N,IC,1,3) = A(N,IC,1,3) - A(N,IC,3,-7)*A(LLDG(N,-7),3,1,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,1,7) - A(N,IC,3,-1)*A(N-1,3,1,4)
      A(N,IC,1,4) = A(N,IC,1,4) - A(N,IC,3,-6)*A(LLDG(N,-6),3,1,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,1,8)
      A(N,IC,1,5) = A(N,IC,1,5) - A(N,IC,3,-9)*A(LLDG(N,-9),3,1,8) - A
     1   (N,IC,3,-3)*A(LLDG(N,-3),3,1,7) - A(N,IC,3,-2)*A(LLDG(N,-2),
     2   3,1,6)
      A(N,IC,1,6) = A(N,IC,1,6) - A(N,IC,3,-4)*A(LLDG(N,-4),3,1,9) - A
     1   (N,IC,3,-1)*A(N-1,3,1,7)
      A(N,IC,1,7) = A(N,IC,1,7) - A(N,IC,3,-3)*A(LLDG(N,-3),3,1,9) - A
     1   (N,IC,3,-1)*A(N-1,3,1,8)
      A(N,IC,1,8) = A(N,IC,1,8) - A(N,IC,3,-2)*A(LLDG(N,-2),3,1,9)
      A(N,IC,2,1) = A(N,IC,2,1) - A(N,IC,3,-7)*A(LLDG(N,-7),3,2,8) - A
     1   (N,IC,3,-6)*A(LLDG(N,-6),3,2,7) - A(N,IC,3,-3)*A(LLDG(N,-3),
     2   3,2,4) - A(N,IC,3,-2)*A(LLDG(N,-2),3,2,3)
      A(N,IC,2,2) = A(N,IC,2,2) - A(N,IC,3,-8)*A(LLDG(N,-8),3,2,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,2,6) - A(N,IC,3,-1)*A(N-1,3,2,3)
      A(N,IC,2,3) = A(N,IC,2,3) - A(N,IC,3,-7)*A(LLDG(N,-7),3,2,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,2,7) - A(N,IC,3,-1)*A(N-1,3,2,4)
      A(N,IC,2,4) = A(N,IC,2,4) - A(N,IC,3,-6)*A(LLDG(N,-6),3,2,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,2,8)
      A(N,IC,2,5) = A(N,IC,2,5) - A(N,IC,3,-9)*A(LLDG(N,-9),3,2,8) - A
     1   (N,IC,3,-3)*A(LLDG(N,-3),3,2,7) - A(N,IC,3,-2)*A(LLDG(N,-2),
     2   3,2,6)
      A(N,IC,2,6) = A(N,IC,2,6) - A(N,IC,3,-4)*A(LLDG(N,-4),3,2,9) - A
     1   (N,IC,3,-1)*A(N-1,3,2,7)
      A(N,IC,2,7) = A(N,IC,2,7) - A(N,IC,3,-3)*A(LLDG(N,-3),3,2,9) - A
     1   (N,IC,3,-1)*A(N-1,3,2,8)
      A(N,IC,2,8) = A(N,IC,2,8) - A(N,IC,3,-2)*A(LLDG(N,-2),3,2,9)
      A(N,IC,3,1) = A(N,IC,3,1) - A(N,IC,3,-7)*A(LLDG(N,-7),3,3,8) - A
     1   (N,IC,3,-6)*A(LLDG(N,-6),3,3,7) - A(N,IC,3,-3)*A(LLDG(N,-3),
     2   3,3,4) - A(N,IC,3,-2)*A(LLDG(N,-2),3,3,3)
      A(N,IC,3,2) = A(N,IC,3,2) - A(N,IC,3,-8)*A(LLDG(N,-8),3,3,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,3,6) - A(N,IC,3,-1)*A(N-1,3,3,3)
      A(N,IC,3,3) = A(N,IC,3,3) - A(N,IC,3,-7)*A(LLDG(N,-7),3,3,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,3,7) - A(N,IC,3,-1)*A(N-1,3,3,4)
      A(N,IC,3,4) = A(N,IC,3,4) - A(N,IC,3,-6)*A(LLDG(N,-6),3,3,9) - A
     1   (N,IC,3,-5)*A(LLDG(N,-5),3,3,8)
      A(N,IC,3,5) = A(N,IC,3,5) - A(N,IC,3,-9)*A(LLDG(N,-9),3,3,8) - A
     1   (N,IC,3,-3)*A(LLDG(N,-3),3,3,7) - A(N,IC,3,-2)*A(LLDG(N,-2),
     2   3,3,6)
      A(N,IC,3,6) = A(N,IC,3,6) - A(N,IC,3,-4)*A(LLDG(N,-4),3,3,9) - A
     1   (N,IC,3,-1)*A(N-1,3,3,7)
      A(N,IC,3,7) = A(N,IC,3,7) - A(N,IC,3,-3)*A(LLDG(N,-3),3,3,9) - A
     1   (N,IC,3,-1)*A(N-1,3,3,8)
      A(N,IC,3,8) = A(N,IC,3,8) - A(N,IC,3,-2)*A(LLDG(N,-2),3,3,9)
  400       CONTINUE
            DO 503 LC = 1, IC-1
CDIR$ IVDEP
            DO 401 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
CFPP$ UNROLL
            DO 504 JC = 1, NPDE
               A(N,IC,JC,1) = A(N,IC,JC,1)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,1)
               A(N,IC,JC,2) = A(N,IC,JC,2)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,2)
               A(N,IC,JC,3) = A(N,IC,JC,3)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,3)
               A(N,IC,JC,4) = A(N,IC,JC,4)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,4)
               A(N,IC,JC,5) = A(N,IC,JC,5)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,5)
               A(N,IC,JC,6) = A(N,IC,JC,6)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,6)
               A(N,IC,JC,7) = A(N,IC,JC,7)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,7)
               A(N,IC,JC,8) = A(N,IC,JC,8)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,8)
               A(N,IC,JC,9) = A(N,IC,JC,9)
     +                      - A(N,IC,LC,0)*A(N,LC,JC,9)
  504       CONTINUE
  401       CONTINUE
  503       CONTINUE
  500    CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKSLV (NPTS, NPD, A, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   B)
      INTEGER NPDE
      PARAMETER (NPDE = 3)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPD, LLDG(NPTS,-9:-2), LUDG(NPTS,2:9),
     +   LSL(*), LLSL(0:*), LSU(*), LLSU(0:*)
      DOUBLE PRECISION A(NPTS,NPDE,NPDE,-9:9), B(NPTS,NPDE)
C
Ccc PURPOSE:
C Solve LUx = b
C    A stems from a 19-point stencil on a grid with irregular row and
C    plane sizes
C A = ILU, vectorized by `hyperplane' ordening, where the `hyperplanes'
C    are sets of points that can be treated simultaneously.
C A((i,j,k),1:NPDE,1:NPDE,.) contains a block of NPDE.NPDE elements
C    corresponding with node (i,j,k)
C
Ccc PARAMETER DESCRIPTION:
C NPTS   : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : IN.  A(.,1:NPDE,1:NPDE,-9:-1)  : lower block diagonals of L
C               A(.,ic,jc,0):      jc < ic: main block diagonal of L
C                                           main diagonal L == I
C                                  jc >=ic: main block diagonal of U
C                                           main diagonal U inverted
C               A(.,1:NPDE,1:NPDE, 1: 9)  : upper block diagonals of U
C LLDG   : IN. Block-column index of lower 8 block-diagonals
C              If block ld does not exist the LLDG(N,ld) = N
C LUDG   : IN. Block-column index of upper 8 block-diagonals
C              If block ud does not exist the LUDG(N,lu) = N
C LSL    : IN. (NPTS)
C          LSL(LLSL(m-1)+1:LLSL(m)): pointers to set of points S_m
C          in actual grid that can be treated at the m-th iteration
C          of Ly = b
C LLSL   : IN. (0:LLSL(0))
C          LLSL(0) = # iterations needed
C          LLSL(1:LLSL(0)): pointers to the start of a list in LSL
C LSU    : IN. (NPTS)
C          LSU(LLSU(m-1)+1:LLSU(m)): pointers to set of points S_m
C          in actual grid that can be treated at the m-th iteration
C          of Ux = y
C LLSU   : IN. (0:LLSU(0))
C          LLSU(0) = # iterations needed
C          LLSU(1:LLSU(0)): pointers to the start of a list in LSU
C B      : INOUT.
C          IN: right-hand side vector b
C          OUT: solution vector x
C
Ccc EXTERNALS USED: NONE
C
C-----------------------------------------------------------------------
C
      INTEGER IC, JC, N, L, M
C
CCC Ly = b
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # N = LSL_m(LLSL(l))
C
C LSL_1 = {(i,j,k)| (i,j,k) not dependent on (i-ii,j-jj,k-kk),
C                                   ii,jj,kk >= 0, not ii=jj=kk=0;
C                   i.e., Dirichlet points and left/down/front corners}
      M = 1
      DO 100 IC = 2, NPDE
      DO 101 JC = 1, IC-1
CDIR$ IVDEP
      DO 1 L = 1, LLSL(M)
         N = LSL(L)
         B(N,IC) = B(N,IC) - A(N,IC,JC,0)*B(N,JC)
    1 CONTINUE
  101 CONTINUE
  100 CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the first row and the first
C     point of the second row, since N > 1 in the loop and for those
C     points LLDG(N,.) = N, (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 10 M = 2, LLSL(0)
C
C   Compute y elements in this set
CDIR$ IVDEP
         DO 20 L = LLSL(M-1)+1, LLSL(M)
            N = LSL(L)
CFPP$ UNROLL
         DO 120 IC = 1, NPDE
CFPP$ UNROLL
         DO 121 JC = 1, NPDE
               B(N,IC) = B(N,IC) - A(N,IC,JC,-1)*B(N- 1,JC)
     +                           - A(N,IC,JC,-2)*B(LLDG(N,-2),JC)
     +                           - A(N,IC,JC,-3)*B(LLDG(N,-3),JC)
     +                           - A(N,IC,JC,-4)*B(LLDG(N,-4),JC)
     +                           - A(N,IC,JC,-5)*B(LLDG(N,-5),JC)
     +                           - A(N,IC,JC,-6)*B(LLDG(N,-6),JC)
     +                           - A(N,IC,JC,-7)*B(LLDG(N,-7),JC)
     +                           - A(N,IC,JC,-8)*B(LLDG(N,-8),JC)
     +                           - A(N,IC,JC,-9)*B(LLDG(N,-9),JC)
  121    CONTINUE
  120    CONTINUE
   20    CONTINUE
         DO 123 IC = 2, NPDE
         DO 122 JC = 1, IC-1
CDIR$ IVDEP
            DO 21 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               B(N,IC) = B(N,IC) - A(N,IC,JC,0)*B(N,JC)
   21       CONTINUE
  122    CONTINUE
  123    CONTINUE
C
   10 CONTINUE
C
CCC Ux = y
C
C Loop over `hyperplanes' LSU_m, m = 1, LLSU(0)
C    Node # N = LSU_m(LLSU(l))
C
C LSU_1 = {(i,j,k)| (i,j,k) not dependent on (i+ii,j+jj,k+kk),
C                                     ii,jj,kk >= 0, not ii=jj=kk=0;
C                   e.g., Dirichlet points and right/up/back corners}
C
      M = 1
      DO 130 IC = NPDE, 1, -1
      DO 131 JC = NPDE, IC+1, -1
CDIR$ IVDEP
      DO 132 L = 1, LLSU(M)
         N = LSU(L)
         B(N,IC) = B(N,IC) - A(N,IC,JC,0)*B(N,JC)
  132 CONTINUE
  131 CONTINUE
CDIR$ IVDEP
      DO 133 L = 1, LLSU(M)
         N = LSU(L)
         B(N,IC) = B(N,IC) * A(N,IC,IC,0)
  133 CONTINUE
  130 CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the last row and the first
C     point of the second last row, since N < NPTS in the loop and for
C     those points LUDG(N,.) = N (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 30 M = 2, LLSU(0)
C
C   Compute x elements in this set
CDIR$ IVDEP
         DO 40 L = LLSU(M-1)+1, LLSU(M)
            N = LSU(L)
CFPP$ UNROLL
         DO 150 IC = NPDE, 1, -1
CFPP$ UNROLL
         DO 151 JC = NPDE, 1, -1
               B(N,IC) = B(N,IC) - A(N,IC,JC,1)*B(N+1      ,JC)
     +                           - A(N,IC,JC,2)*B(LUDG(N,2),JC)
     +                           - A(N,IC,JC,3)*B(LUDG(N,3),JC)
     +                           - A(N,IC,JC,4)*B(LUDG(N,4),JC)
     +                           - A(N,IC,JC,5)*B(LUDG(N,5),JC)
     +                           - A(N,IC,JC,6)*B(LUDG(N,6),JC)
     +                           - A(N,IC,JC,7)*B(LUDG(N,7),JC)
     +                           - A(N,IC,JC,8)*B(LUDG(N,8),JC)
     +                           - A(N,IC,JC,9)*B(LUDG(N,9),JC)
  151    CONTINUE
  150    CONTINUE
   40    CONTINUE
         DO 1150 IC = NPDE, 1, -1
         DO 152 JC = NPDE, IC+1, -1
CDIR$ IVDEP
         DO 51 L = LLSU(M-1)+1, LLSU(M)
            N = LSU(L)
               B(N,IC) = B(N,IC) - A(N,IC,JC,0)*B(N,JC)
   51    CONTINUE
  152    CONTINUE
CDIR$ IVDEP
         DO 52 L = LLSU(M-1)+1, LLSU(M)
            N = LSU(L)
               B(N,IC) = B(N,IC) * A(N,IC,IC,0)
   52    CONTINUE
 1150    CONTINUE
C
   30 CONTINUE
C
      RETURN
      END
