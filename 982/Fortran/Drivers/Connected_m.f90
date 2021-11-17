module Connected_m

  ! Find the connected components of a graph represented by a matrix

  implicit NONE
  public :: Connected

  interface Connected
    module procedure Connected_d, Connected_s
  end interface

  interface BFS
    module procedure BFS_d, BFS_s
  end interface

contains

  subroutine Connected_d ( A, P, Blocks, NB )
    ! Find the connected components of A.  P is a permutation that
    ! brings the components together.
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: A(:,:)     ! Square Array that represents the graph.
                                       ! If A(I,J) is nonzero there is an edge
                                       ! from I to J (and A(I,J) is its label).
    integer, intent(out) :: P(:)       ! Permutation of A (Same size as first
                                       ! dimension of A).
    integer, intent(out) :: Blocks(0:) ! Block N of P is bounded by Block(n-1)+1
                                       ! and Block(n) (UBOUND(.) same as for P)
    integer, intent(out) :: NB         ! Number of blocks found.

    integer :: I, N
    logical :: Visited(size(p))

    visited = .false.
    blocks(0) = 0
    n = 0 ! Last used element of P
    nb = 0
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( count(a(i,:) /= 0) == 1 .and. a(i,i) /= 0 ) then
        call bfs ( a, i, visited, p, n )
        nb = nb + 1
        blocks(nb) = n ! Last used element of P
      end if
    end do
    ! Now fill in P for the ones that weren't visited
    if ( blocks(nb) < size(a,1) ) then
      do i = 1, size(a,1)
        if ( .not. visited(i) ) then
          n = n + 1
          p(n) = i
        end if
      end do
      nb = nb + 1
      blocks(nb) = size(a,1)
    end if
  end subroutine Connected_d

  subroutine Connected_s ( A, P, Blocks, NB )
    ! Find the connected components of A.  P is a permutation that
    ! brings the components together.
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: A(:,:)     ! Square Array that represents the graph.
                                       ! If A(I,J) is nonzero there is an edge
                                       ! from I to J (and A(I,J) is its label).
    integer, intent(out) :: P(:)       ! Permutation of A (Same size as first
                                       ! dimension of A).
    integer, intent(out) :: Blocks(0:) ! Block N of P is bounded by Block(n-1)+1
                                       ! and Block(n) (UBOUND(.) same as for P)
    integer, intent(out) :: NB         ! Number of blocks found.

    integer :: I, N
    logical :: Visited(size(p))

    visited = .false.
    blocks(0) = 0
    n = 0 ! Last used element of P
    nb = 0
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( count(a(i,:) /= 0) == 1 .and. a(i,i) /= 0 ) then
        call bfs ( a, i, visited, p, n )
        nb = nb + 1
        blocks(nb) = n ! Last used element of P
      end if
    end do
    ! Now fill in P for the ones that weren't visited
    if ( blocks(nb) < size(a,1) ) then
      do i = 1, size(a,1)
        if ( .not. visited(i) ) then
          n = n + 1
          p(n) = i
        end if
      end do
      nb = nb + 1
      blocks(nb) = size(a,1)
    end if
  end subroutine Connected_s

  recursive subroutine BFS_d ( A, Start, Visited, P, N )
    ! Do a breadth-first traversal of the adjacency graph represented by A,
    ! starting at Start
    integer, parameter :: RK = kind(0.0d0)
    real(rk),  intent(in) :: A(:,:)    ! Array to be permuted (Square)
    integer, intent(in) :: Start
    logical, intent(inout) :: Visited(:)
    integer, intent(inout) :: P(:)     ! Permutation of A
    integer, intent(inout) :: N        ! Last element of P used so far

    integer :: I, J

    visited(start) = .true.
    n = n + 1
    p(n) = start
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( a(i,start) /= 0 ) call bfs ( a, i, visited, p, n )
    end do
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( a(start,i) /= 0 ) call bfs ( a, i, visited, p, n )
    end do
  end subroutine BFS_d

  recursive subroutine BFS_s ( A, Start, Visited, P, N )
    ! Do a breadth-first traversal of the adjacency graph represented by A,
    ! starting at Start
    integer, parameter :: RK = kind(0.0e0)
    real(rk),  intent(in) :: A(:,:)    ! Array to be permuted (Square)
    integer, intent(in) :: Start
    logical, intent(inout) :: Visited(:)
    integer, intent(inout) :: P(:)     ! Permutation of A
    integer, intent(inout) :: N        ! Last element of P used so far

    integer :: I, J

    visited(start) = .true.
    n = n + 1
    p(n) = start
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( a(i,start) /= 0 ) call bfs ( a, i, visited, p, n )
    end do
    do i = 1, size(a,1)
      if ( visited(i) ) cycle
      if ( a(start,i) /= 0 ) call bfs ( a, i, visited, p, n )
    end do
  end subroutine BFS_s

end module Connected_m
