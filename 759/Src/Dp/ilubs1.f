      SUBROUTINE ILU (NPTS, NPD, A, LLDG, LSL, LLSL)
      INTEGER NPDE
      PARAMETER (NPDE = 1)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPD, LLDG(NPTS,-9:-2), LSL(*), LLSL(0:*)
      DOUBLE PRECISION A(NPTS,-9:9)
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
C              If block ld does not exist the LLDG(K,ld) = K
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
      INTEGER N, L, M
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
CDIR$ IVDEP
            DO 553 L = 1, LLSL(M)
               N = LSL(L)
               A(N,0) = 1.0 / A(N,0)
  553       CONTINUE
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
CDIR$ IVDEP
         DO 20 L = LLSL(M-1)+1, LLSL(M)
               N = LSL(L)
               A(N,-9) = (A(N,-9)) * A(LLDG(N,-9),0)
               A(N,-8) = (A(N,-8) - A(N,-9)*A(LLDG(N,-9),2))
     +                   * A(LLDG(N,-8),0)
               A(N,-7) = (A(N,-7) - A(N,-9)*A(LLDG(N,-9),3)
     +                            - A(N,-8)*A(LLDG(N,-8),1))
     +                   * A(LLDG(N,-7),0)
               A(N,-6) = (A(N,-6) - A(N,-9)*A(LLDG(N,-9),4)
     +                            - A(N,-7)*A(LLDG(N,-7),1))
     +                   * A(LLDG(N,-6),0)
               A(N,-5) = (A(N,-5) - A(N,-8)*A(LLDG(N,-8),4)
     +                            - A(N,-7)*A(LLDG(N,-7),3)
     +                            - A(N,-6)*A(LLDG(N,-6),2))
     +                   * A(LLDG(N,-5),0)
               A(N,-4) = (A(N,-4) - A(N,-9)*A(LLDG(N,-9),6)
     +                            - A(N,-8)*A(LLDG(N,-8),5))
     +                   * A(LLDG(N,-4),0)
               A(N,-3) = (A(N,-3) - A(N,-9)*A(LLDG(N,-9),8)
     +                            - A(N,-7)*A(LLDG(N,-7),5)
     +                            - A(N,-4)*A(LLDG(N,-4),1))
     +                   * A(LLDG(N,-3),0)
               A(N,-2) = (A(N,-2) - A(N,-9)*A(LLDG(N,-9),8)
     +                            - A(N,-6)*A(LLDG(N,-6),5)
     +                            - A(N,-3)*A(LLDG(N,-3),1))
     +                   * A(LLDG(N,-2),0)
               A(N,-1) = (A(N,-1) - A(N,-8)*A(LLDG(N,-8),7)
     +                            - A(N,-7)*A(LLDG(N,-7),6)
     +                            - A(N,-4)*A(LLDG(N,-4),3)
     +                            - A(N,-3)*A(LLDG(N,-3),2))
     +                   * A(N-1,0)
C
C   Compute main diagonal
               A(N,0) = 1.0 / (A(N, 0) - A(N,-9)*A(LLDG(N,-9),9)
     +                                 - A(N,-8)*A(LLDG(N,-8),8)
     +                                 - A(N,-7)*A(LLDG(N,-7),7)
     +                                 - A(N,-6)*A(LLDG(N,-6),6)
     +                                 - A(N,-5)*A(LLDG(N,-5),5)
     +                                 - A(N,-4)*A(LLDG(N,-4),4)
     +                                 - A(N,-3)*A(LLDG(N,-3),3)
     +                                 - A(N,-2)*A(LLDG(N,-2),2)
     +                                 - A(N,-1)*A(N-1       ,1))
C
C   Compute upper diagonals
               A(N,1) = A(N, 1) - A(N,-7)*A(LLDG(N,-7),8)
     +                          - A(N,-6)*A(LLDG(N,-6),7)
     +                          - A(N,-3)*A(LLDG(N,-3),4)
     +                          - A(N,-2)*A(LLDG(N,-2),3)
               A(N,2) = A(N, 2) - A(N,-8)*A(LLDG(N,-8),9)
     +                          - A(N,-5)*A(LLDG(N,-5),6)
     +                          - A(N,-1)*A(N-1       ,3)
               A(N,3) = A(N, 3) - A(N,-7)*A(LLDG(N,-7),9)
     +                          - A(N,-5)*A(LLDG(N,-5),7)
     +                          - A(N,-1)*A(N-1       ,4)
               A(N,4) = A(N, 4) - A(N,-6)*A(LLDG(N,-6),9)
     +                          - A(N,-5)*A(LLDG(N,-5),8)
               A(N,5) = A(N, 5) - A(N,-9)*A(LLDG(N,-9),8)
     +                          - A(N,-3)*A(LLDG(N,-3),7)
     +                          - A(N,-2)*A(LLDG(N,-2),6)
               A(N,6) = A(N, 6) - A(N,-4)*A(LLDG(N,-4),9)
     +                          - A(N,-1)*A(N-1       ,7)
               A(N,7) = A(N, 7) - A(N,-3)*A(LLDG(N,-3),9)
     +                          - A(N,-1)*A(N-1       ,8)
               A(N,8) = A(N, 8) - A(N,-2)*A(LLDG(N,-2),9)
   20    CONTINUE
C
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKSLV (NPTS, NPD, A, LLDG, LUDG, LSL, LLSL, LSU, LLSU,
     +   B)
      INTEGER NPDE
      PARAMETER (NPDE = 1)
C
C-----------------------------------------------------------------------
C
Ccc PARAMETER SPECIFICATION:
      INTEGER NPTS, NPD, LLDG(NPTS,-9:-2), LUDG(NPTS,2:9),
     +   LSL(*), LLSL(0:*), LSU(*), LLSU(0:*)
      DOUBLE PRECISION A(NPTS,-9:9), B(NPTS)
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
C              If block ld does not exist the LLDG(K,ld) = K
C LUDG   : IN. Block-column index of upper 8 block-diagonals
C              If block ud does not exist the LUDG(K,lu) = K
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
      INTEGER N, L, M
C
CCC Ly = b
C
C Loop over `hyperplanes' S_m, m = 1, LLSL(0)
C    Node # N = LSL_m(LLSL(l))
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
               B(N) = B(N) - A(N,-1)*B(N- 1)
     +                     - A(N,-2)*B(LLDG(N,-2))
     +                     - A(N,-3)*B(LLDG(N,-3))
     +                     - A(N,-4)*B(LLDG(N,-4))
     +                     - A(N,-5)*B(LLDG(N,-5))
     +                     - A(N,-6)*B(LLDG(N,-6))
     +                     - A(N,-7)*B(LLDG(N,-7))
     +                     - A(N,-8)*B(LLDG(N,-8))
     +                     - A(N,-9)*B(LLDG(N,-9))
   20    CONTINUE
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
CDIR$ IVDEP
      DO 133 L = 1, LLSU(M)
         N = LSU(L)
         B(N) = B(N) * A(N,0)
  133 CONTINUE
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
               B(N) = (B(N) - A(N,1)*B(N+1     )
     +                      - A(N,2)*B(LUDG(N,2))
     +                      - A(N,3)*B(LUDG(N,3))
     +                      - A(N,4)*B(LUDG(N,4))
     +                      - A(N,5)*B(LUDG(N,5))
     +                      - A(N,6)*B(LUDG(N,6))
     +                      - A(N,7)*B(LUDG(N,7))
     +                      - A(N,8)*B(LUDG(N,8))
     +                      - A(N,9)*B(LUDG(N,9))) * A(N,0)
   40    CONTINUE
C
   30 CONTINUE
C
      RETURN
      END
