      SUBROUTINE ILU (N, NPD, A, LLDG, LSL, LLSL)
      INTEGER NPDE
      PARAMETER (NPDE = 2)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, NPD, LLDG(N,-4:-2), LSL(*), LLSL(0:*)
      REAL A(N,NPDE,NPDE,-4:4)
C
Ccc PURPOSE:
C Incomplete LU decomposition of A
C    A stems from a 9-point stencil on a grid with irregular row lengths
C A = ILU, vectorized by `hyperplane' ordening, where the `hyperplanes'
C    are sets of points that can be treated simultaneously.
C A((i,j),1:NPDE,1:NPDE,.) contains a block of NPDE.NPDE elements
C    corresponding with node (i,j)
C
Ccc PARAMETER DESCRIPTION:
C N      : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : INOUT.
C          IN: A(.,1:NPDE,1:NPDE,-4:-1)  : lower block diagonals
C              A(.,1:NPDE,1:NPDE,0)      :   main block diagonal
C              A(.,1:NPDE,1:NPDE, 1: 4)  : upper block diagonals
C          OUT: A(.,1:NPDE,1:NPDE,-4:-1)  : lower block diagonals of L
C               A(.,ic,jc,0):      jc < ic: main block diagonal of L
C                                           main diagonal L == I
C                                  jc >=ic: main block diagonal of U
C                                           main diagonal U inverted
C               A(.,1:NPDE,1:NPDE, 1: 4)  : upper block diagonals of U
C LLDG   : IN. Block-column index of lower 3 block-diagonals
C              LLDG(K,.) = K, K = 1, NX
C LSL    : IN. (NPTS)
C          LSL(LLSL(m-1)+1:LLSL(m)): pointers to set of points S_m
C          in actual grid that can be treated at the m-th iteration
C LLSL   : IN. (0:LLSL(0))
C          LLSL(0) = # iterations needed
C          LLSL(1:LLSL(0)): pointers to the start of a list in LSL
C
Ccc EXTERNALS USED: NONE
C
C ----------------------------------------------------------------------
C
      INTEGER IC, JC, LC, K, L, M
C
      IF (NPDE .NE. NPD) STOP 'Wrong ILUBS loaded.'
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # K = S_m(LLSL(l))
C
C S_1 = {(i,j)| (i,j) not dependent on (i-ii,j-jj),
C                                   ii,jj >= 0, not ii=jj=0;
C                i.e., Dirichlet points and lower left corner (K=1)}
      M = 1
C
C   Compute main block diagonal
         DO 550 IC = 1, NPDE
            DO 554 LC = 1, IC-1
            DO 555 JC = IC, NPDE
CDIR$ IVDEP
            DO 551 L = 1, LLSL(M)
               K = LSL(L)
               A(K,IC,JC,0) = A(K,IC,JC,0)
     +                        - A(K,IC,LC,0)*A(K,LC,JC,0)
  551       CONTINUE
  555       CONTINUE
            DO 556 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 552 L = 1, LLSL(M)
               K = LSL(L)
               A(K,JC,IC,0) = A(K,JC,IC,0)
     +                        - A(K,JC,LC,0)*A(K,LC,IC,0)
  552       CONTINUE
  556       CONTINUE
  554       CONTINUE
CDIR$ IVDEP
            DO 553 L = 1, LLSL(M)
               K = LSL(L)
               A(K,IC,IC,0) = 1.0 / A(K,IC,IC,0)
  553       CONTINUE
            DO 557 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 559 L = 1, LLSL(M)
               K = LSL(L)
               A(K,JC,IC,0) = A(K,JC,IC,0) * A(K,IC,IC,0)
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
               K = LSL(L)
               A(K,IC,JC,1) = A(K,IC,JC,1)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,1)
               A(K,IC,JC,2) = A(K,IC,JC,2)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,2)
               A(K,IC,JC,3) = A(K,IC,JC,3)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,3)
               A(K,IC,JC,4) = A(K,IC,JC,4)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,4)
  561       CONTINUE
  564       CONTINUE
  563       CONTINUE
  560    CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the first row and the first
C     point of the second row, since K > 1 in the loop and
C     LLDG(K,.) = K, K = 1, NX (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 10 M = 2, LLSL(0)
C
C   Compute lower block diagonals
         DO 120 JC = 1, NPDE
            DO 121 LC = 1, JC-1
CDIR$ IVDEP
            DO 20 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 122 IC = 1, NPDE
               A(K,IC,JC,-4) = A(K,IC,JC,-4)
     +                       - A(K,IC,LC,-4)*A(LLDG(K,-4),LC,JC,0)
  122       CONTINUE
   20       CONTINUE
  121       CONTINUE
CDIR$ IVDEP
            DO 21 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 123 IC = 1, NPDE
               A(K,IC,JC,-4) = A(K,IC,JC,-4) * A(LLDG(K,-4),JC,JC,0)
  123       CONTINUE
   21       CONTINUE
  120    CONTINUE
         DO 130 JC = 1, NPDE
CDIR$ IVDEP
            DO 30 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 131 LC = 1, NPDE
CFPP$ UNROLL
            DO 132 IC = 1, NPDE
               A(K,IC,JC,-3) = A(K,IC,JC,-3)
     +                       - A(K,IC,LC,-4)*A(LLDG(K,-4),LC,JC,1)
  132       CONTINUE
  131       CONTINUE
   30       CONTINUE
            DO 133 LC = 1, JC-1
CDIR$ IVDEP
            DO 31 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 134 IC = 1, NPDE
               A(K,IC,JC,-3) = A(K,IC,JC,-3)
     +                       - A(K,IC,LC,-3)*A(LLDG(K,-3),LC,JC,0)
  134       CONTINUE
   31       CONTINUE
  133       CONTINUE
CDIR$ IVDEP
            DO 32 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 135 IC = 1, NPDE
               A(K,IC,JC,-3) = A(K,IC,JC,-3) * A(LLDG(K,-3),JC,JC,0)
  135       CONTINUE
   32       CONTINUE
  130    CONTINUE
         DO 140 JC = 1, NPDE
CDIR$ IVDEP
            DO 40 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 141 LC = 1, NPDE
CFPP$ UNROLL
            DO 142 IC = 1, NPDE
               A(K,IC,JC,-2) = A(K,IC,JC,-2)
     +                       - A(K,IC,LC,-3)*A(LLDG(K,-3),LC,JC,1)
               A(K,IC,JC,-1) = A(K,IC,JC,-1)
     +                       - A(K,IC,LC,-4)*A(LLDG(K,-4),LC,JC,3)
     +                       - A(K,IC,LC,-3)*A(LLDG(K,-3),LC,JC,2)
  142       CONTINUE
  141       CONTINUE
   40       CONTINUE
            DO 143 LC = 1, JC-1
CDIR$ IVDEP
            DO 41 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 144 IC = 1, NPDE
               A(K,IC,JC,-2) = A(K,IC,JC,-2)
     +                       - A(K,IC,LC,-2)*A(LLDG(K,-2),LC,JC,0)
               A(K,IC,JC,-1) = A(K,IC,JC,-1)
     +                       - A(K,IC,LC,-1)*A(K-1      ,LC,JC,0)
  144       CONTINUE
   41       CONTINUE
  143       CONTINUE
CDIR$ IVDEP
            DO 42 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 145 IC = 1, NPDE
               A(K,IC,JC,-2) = A(K,IC,JC,-2) * A(LLDG(K,-2),JC,JC,0)
               A(K,IC,JC,-1) = A(K,IC,JC,-1) * A(K-1      ,JC,JC,0)
  145       CONTINUE
   42       CONTINUE
  140    CONTINUE
C
C   Compute main block diagonal
         DO 150 IC = 1, NPDE
            DO 152 JC = IC, NPDE
CDIR$ IVDEP
            DO 50 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 151 LC = 1, NPDE
               A(K,IC,JC,0) = A(K,IC,JC, 0)
     +                      - A(K,IC,LC,-4)*A(LLDG(K,-4),LC,JC,4)
     +                      - A(K,IC,LC,-3)*A(LLDG(K,-3),LC,JC,3)
     +                      - A(K,IC,LC,-2)*A(LLDG(K,-2),LC,JC,2)
     +                      - A(K,IC,LC,-1)*A(K-1      ,LC,JC,1)
  151       CONTINUE
   50       CONTINUE
  152       CONTINUE
            DO 153 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 51 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 1151 LC = 1, NPDE
               A(K,JC,IC,0) = A(K,JC,IC, 0)
     +                      - A(K,JC,LC,-4)*A(LLDG(K,-4),LC,IC,4)
     +                      - A(K,JC,LC,-3)*A(LLDG(K,-3),LC,IC,3)
     +                      - A(K,JC,LC,-2)*A(LLDG(K,-2),LC,IC,2)
     +                      - A(K,JC,LC,-1)*A(K-1      ,LC,IC,1)
 1151       CONTINUE
   51       CONTINUE
  153       CONTINUE
            DO 154 LC = 1, IC-1
            DO 155 JC = IC, NPDE
CDIR$ IVDEP
            DO 52 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               A(K,IC,JC,0) = A(K,IC,JC,0)
     +                      - A(K,IC,LC,0)*A(K        ,LC,JC,0)
   52       CONTINUE
  155       CONTINUE
            DO 156 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 53 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               A(K,JC,IC,0) = A(K,JC,IC,0)
     +                      - A(K,JC,LC,0)*A(K        ,LC,IC,0)
   53       CONTINUE
  156       CONTINUE
  154       CONTINUE
CDIR$ IVDEP
            DO 54 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               A(K,IC,IC,0) = 1.0 / A(K,IC,IC,0)
   54       CONTINUE
            DO 157 JC = IC+1, NPDE
CDIR$ IVDEP
            DO 55 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               A(K,JC,IC,0) = A(K,JC,IC,0) * A(K,IC,IC,0)
   55       CONTINUE
  157       CONTINUE
  150    CONTINUE
C
C   Compute upper block diagonals
         DO 160 IC = 1, NPDE
CDIR$ IVDEP
            DO 60 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 161 LC = 1, NPDE
CFPP$ UNROLL
            DO 162 JC = 1, NPDE
               A(K,IC,JC,1) = A(K,IC,JC, 1)
     +                      - A(K,IC,LC,-3)*A(LLDG(K,-3),LC,JC,4)
     +                      - A(K,IC,LC,-2)*A(LLDG(K,-2),LC,JC,3)
               A(K,IC,JC,2) = A(K,IC,JC, 2)
     +                      - A(K,IC,LC,-1)*A(K-1      ,LC,JC,3)
               A(K,IC,JC,3) = A(K,IC,JC, 3)
     +                      - A(K,IC,LC,-1)*A(K-1      ,LC,JC,4)
  162       CONTINUE
  161       CONTINUE
   60       CONTINUE
            DO 163 LC = 1, IC-1
CDIR$ IVDEP
            DO 61 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
CFPP$ UNROLL
            DO 164 JC = 1, NPDE
               A(K,IC,JC,1) = A(K,IC,JC,1)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,1)
               A(K,IC,JC,2) = A(K,IC,JC,2)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,2)
               A(K,IC,JC,3) = A(K,IC,JC,3)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,3)
               A(K,IC,JC,4) = A(K,IC,JC,4)
     +                      - A(K,IC,LC,0)*A(K,LC,JC,4)
  164       CONTINUE
   61       CONTINUE
  163       CONTINUE
  160    CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKSLV (N, NPD, A, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   B)
      INTEGER NPDE
      PARAMETER (NPDE = 2)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, NPD, LLDG(N,-4:-2), LUDG(N,2:4),
     +   LSL(*), LLSL(0:*), LSU(*), LLSU(0:*)
      REAL A(N,NPDE,NPDE,-4:4), B(N,NPDE)
C
Ccc PURPOSE:
C Solve LUx = b
C    A stems from a 9-point stencil on a grid with irregular row lengths
C A = ILU, vectorized by `hyperplane' ordening, where the `hyperplanes'
C    are sets of points that can be treated simultaneously.
C A((i,j),1:NPDE,1:NPDE,.) contains a block of NPDE.NPDE elements
C    corresponding with node (i,j)
C
Ccc PARAMETER DESCRIPTION:
C N      : IN. # gridpoints
C NPDE   : IN. # components of the PDE
C A      : IN. A(.,1:NPDE,1:NPDE,-4:-1)  : lower block diagonals of L
C              A(.,ic,jc,0):      jc < ic: main block diagonal of L
C                                          main diagonal L == I
C                                 jc >=ic: main block diagonal of U
C                                          main diagonal U inverted
C              A(.,1:NPDE,1:NPDE, 1: 4)  : upper block diagonals of U
C LLDG   : IN. Block-column index of lower 3 block-diagonals
C              LLDG(K,.) = K, K = 1, NX
C LUDG   : IN. Block-column index of upper 3 block-diagonals
C              LUDG(K,.) = K, K = N-NX+1, N
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
C ----------------------------------------------------------------------
C
      INTEGER IC, JC, K, L, M
C
CCC Ly = b
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # K = LSL_m(LLSL(l))
C
C LSL_1 = {(i,j)| (i,j) not dependent on (i-ii,j-jj),
C                                     ii,jj >= 0, not ii=jj=0;
C                 i.e., Dirichlet points and lower left corner (K=1)}
      M = 1
      DO 100 IC = 2, NPDE
      DO 101 JC = 1, IC-1
CDIR$ IVDEP
      DO 1 L = 1, LLSL(M)
         K = LSL(L)
         B(K,IC) = B(K,IC) - A(K,IC,JC,0)*B(K,JC)
    1 CONTINUE
  101 CONTINUE
  100 CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the first row and the first
C     point of the second row, since K > 1 in the loop and
C     LLDG(K,.) = K, K = 1, NX (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 10 M = 2, LLSL(0)
C
C   Compute y elements in this set
CDIR$ IVDEP
         DO 20 L = LLSL(M-1)+1, LLSL(M)
            K = LSL(L)
CFPP$ UNROLL
         DO 120 IC = 1, NPDE
CFPP$ UNROLL
         DO 121 JC = 1, NPDE
               B(K,IC) = B(K,IC) - A(K,IC,JC,-1)*B(K- 1     ,JC)
     +                           - A(K,IC,JC,-2)*B(LLDG(K,-2),JC)
     +                           - A(K,IC,JC,-3)*B(LLDG(K,-3),JC)
     +                           - A(K,IC,JC,-4)*B(LLDG(K,-4),JC)
  121    CONTINUE
  120    CONTINUE
   20    CONTINUE
         DO 123 IC = 2, NPDE
         DO 122 JC = 1, IC-1
CDIR$ IVDEP
            DO 21 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               B(K,IC) = B(K,IC) - A(K,IC,JC,0)*B(K,JC)
   21       CONTINUE
  122    CONTINUE
  123    CONTINUE
C
   10 CONTINUE
C
CCC Ux = y
C
C Loop over `hyperplanes' LSU_m, m = 1, LLSU(0)
C    Node # K = LSU_m(LLSU(l))
C
C LSU_1 = {(i,j)| (i,j) not dependent on (i+ii,j+jj),
C                                     ii,jj >= 0, not ii=jj=0;
C                 e.g., Dirichlet points and upper right corner (K=N)}
C
      M = 1
      DO 130 IC = NPDE, 1, -1
      DO 131 JC = NPDE, IC+1, -1
CDIR$ IVDEP
      DO 132 L = 1, LLSU(M)
         K = LSU(L)
         B(K,IC) = B(K,IC) - A(K,IC,JC,0)*B(K,JC)
  132 CONTINUE
  131 CONTINUE
CDIR$ IVDEP
      DO 133 L = 1, LLSU(M)
         K = LSU(L)
         B(K,IC) = B(K,IC) * A(K,IC,IC,0)
  133 CONTINUE
  130 CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the last row and the first
C     point of the second last row, since K < N in the loop and
C     LUDG(K,.) = K, K = N-NX+1, N (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 30 M = 2, LLSU(0)
C
C   Compute x elements in this set
CDIR$ IVDEP
         DO 40 L = LLSU(M-1)+1, LLSU(M)
            K = LSU(L)
CFPP$ UNROLL
         DO 150 IC = NPDE, 1, -1
CFPP$ UNROLL
         DO 151 JC = NPDE, 1, -1
               B(K,IC) = B(K,IC) - A(K,IC,JC,1)*B(K+1     ,JC)
     +                           - A(K,IC,JC,2)*B(LUDG(K,2),JC)
     +                           - A(K,IC,JC,3)*B(LUDG(K,3),JC)
     +                           - A(K,IC,JC,4)*B(LUDG(K,4),JC)
  151    CONTINUE
  150    CONTINUE
   40    CONTINUE
         DO 1150 IC = NPDE, 1, -1
         DO 152 JC = NPDE, IC+1, -1
CDIR$ IVDEP
         DO 51 L = LLSU(M-1)+1, LLSU(M)
            K = LSU(L)
               B(K,IC) = B(K,IC) - A(K,IC,JC,0)*B(K,JC)
   51    CONTINUE
  152    CONTINUE
CDIR$ IVDEP
         DO 52 L = LLSU(M-1)+1, LLSU(M)
            K = LSU(L)
               B(K,IC) = B(K,IC) * A(K,IC,IC,0)
   52    CONTINUE
 1150    CONTINUE
C
   30 CONTINUE
C
      RETURN
      END
