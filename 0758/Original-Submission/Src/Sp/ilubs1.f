      SUBROUTINE ILU (N, NPD, A, LLDG, LSL, LLSL)
      INTEGER NPDE
      PARAMETER (NPDE = 1)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, NPD, LLDG(N,-4:-2), LSL(*), LLSL(0:*)
      REAL A(N,-4:4)
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
C-----------------------------------------------------------------------
C
      INTEGER K, L, M
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
CDIR$ IVDEP
            DO 553 L = 1, LLSL(M)
               K = LSL(L)
               A(K,0) = 1.0 / A(K,0)
  553       CONTINUE
C
C Loop over rest of `hyperplane' sets
C NB. No exception handling is necessary for the first row and the first
C     point of the second row, since K > 1 in the loop and
C     LLDG(K,.) = K, K = 1, NX (=> no array index problems),
C     and the multiplicator of `non existing' array elements is zero.
C
      DO 10 M = 2, LLSL(0)
C
C   Compute lower diagonals
CDIR$ IVDEP
         DO 20 L = LLSL(M-1)+1, LLSL(M)
               K = LSL(L)
               A(K,-4) = A(K,-4) * A(LLDG(K,-4),0)
               A(K,-3) = (A(K,-3) - A(K,-4)*A(LLDG(K,-4),1))
     +                   * A(LLDG(K,-3),0)
               A(K,-2) = (A(K,-2) - A(K,-3)*A(LLDG(K,-3),1))
     +                   * A(LLDG(K,-2),0)
               A(K,-1) = (A(K,-1) - A(K,-4)*A(LLDG(K,-4),3)
     +                            - A(K,-3)*A(LLDG(K,-3),2)) * A(K-1,0)
C
C   Compute main diagonal
               A(K,0) = 1.0 / (A(K, 0) - A(K,-4)*A(LLDG(K,-4),4)
     +                                 - A(K,-3)*A(LLDG(K,-3),3)
     +                                 - A(K,-2)*A(LLDG(K,-2),2)
     +                                 - A(K,-1)*A(K-1       ,1))
C
C   Compute upper diagonals
               A(K,1) = A(K, 1) - A(K,-3)*A(LLDG(K,-3),4)
     +                          - A(K,-2)*A(LLDG(K,-2),3)
               A(K,2) = A(K, 2) - A(K,-1)*A(K-1       ,3)
               A(K,3) = A(K, 3) - A(K,-1)*A(K-1       ,4)
   20    CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKSLV (N, NPD, A, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   B)
      INTEGER NPDE
      PARAMETER (NPDE = 1)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER N, NPD, LLDG(N,-4:-2), LUDG(N,2:4),
     +   LSL(*), LLSL(0:*), LSU(*), LLSU(0:*)
      REAL A(N,-4:4), B(N)
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
C-----------------------------------------------------------------------
C
      INTEGER K, L, M
C
CCC Ly = b
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # K = LSL_m(LLSL(l))
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
               B(K) = B(K) - A(K,-1)*B(K- 1)
     +                     - A(K,-2)*B(LLDG(K,-2))
     +                     - A(K,-3)*B(LLDG(K,-3))
     +                     - A(K,-4)*B(LLDG(K,-4))
   20    CONTINUE
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
CDIR$ IVDEP
      DO 133 L = 1, LLSU(M)
         K = LSU(L)
         B(K) = B(K) * A(K,0)
  133 CONTINUE
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
               B(K) = (B(K) - A(K,1)*B(K+1     )
     +                      - A(K,2)*B(LUDG(K,2))
     +                      - A(K,3)*B(LUDG(K,3))
     +                      - A(K,4)*B(LUDG(K,4))) * A(K,0)
   40    CONTINUE
C
   30 CONTINUE
C
      RETURN
      END
