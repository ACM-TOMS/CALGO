module Triangularize_m

  ! Compute a permutation that puts a matrix that represents a DAG
  ! into lower-triangular form.

  implicit NONE
  public :: Triangularize

  interface Triangularize
    module procedure Triangularize_d, Triangularize_s
  end interface

contains

  subroutine Triangularize_d ( A, P, Blocks )
    ! Permute blocks of A to triangular form
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: A(:,:)    ! Array to be permuted (Square)
    integer, intent(inout) :: P(:)    ! Permutation of A (Same size as first
                                      ! dimension of A)
    integer, intent(in) :: Blocks(0:) ! Block N of P is bounded by Block(n-1)+1
                                      ! and Block(n) (UBOUND(.) same as for P)

    integer :: I, J, K, L, M, NB
    integer :: N                      ! Column of negative unvisited element
    integer :: NNZ                    ! Number of nonzero unvisited elements
    logical :: Visited(size(p))
    logical :: X                      ! A permutation was exchanged

    visited = .false.
    nb = ubound(blocks,1)
    do i = 1, nb
      ! Look for an unvisited row of A(P) that has only one unvisited nonzero
      ! element, and that element is negative.
      m = blocks(i-1)
      x = .true.
      do while ( x )
        x = .false.
        do j = blocks(i-1)+1, blocks(i)
          if ( visited(p(j)) ) cycle
          nnz = 0
          n = 0
          do k = blocks(i-1)+1, blocks(i)
            if ( visited(p(k)) ) cycle
            if ( a(p(j),p(k)) /= 0 ) nnz = nnz + 1
            if ( a(p(j),p(k)) < 0 ) n = k
          end do
          if ( nnz == 1 .and. n /= 0 ) then
            visited(p(n)) = .true.
            m = m + 1
            ! Exchange P(m) and P(n) and mark P(n) visited
            l = p(m)
            p(m) = p(n)
            p(n) = l
            x = .true.
          end if
        end do
      end do
    end do
  end subroutine Triangularize_d

  subroutine Triangularize_s ( A, P, Blocks )
    ! Permute blocks of A to triangular form
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: A(:,:)    ! Array to be permuted (Square)
    integer, intent(inout) :: P(:)    ! Permutation of A (Same size as first
                                      ! dimension of A)
    integer, intent(in) :: Blocks(0:) ! Block N of P is bounded by Block(n-1)+1
                                      ! and Block(n) (UBOUND(.) same as for P)

    integer :: I, J, K, L, M, NB
    integer :: N                      ! Column of negative unvisited element
    integer :: NNZ                    ! Number of nonzero unvisited elements
    logical :: Visited(size(p))
    logical :: X                      ! A permutation was exchanged

    visited = .false.
    nb = ubound(blocks,1)
    do i = 1, nb
      ! Look for an unvisited row of A(P) that has only one unvisited nonzero
      ! element, and that element is negative.
      m = blocks(i-1)
      x = .true.
      do while ( x )
        x = .false.
        do j = blocks(i-1)+1, blocks(i)
          if ( visited(p(j)) ) cycle
          nnz = 0
          n = 0
          do k = blocks(i-1)+1, blocks(i)
            if ( visited(p(k)) ) cycle
            if ( a(p(j),p(k)) /= 0 ) nnz = nnz + 1
            if ( a(p(j),p(k)) < 0 ) n = k
          end do
          if ( nnz == 1 .and. n /= 0 ) then
            visited(p(n)) = .true.
            m = m + 1
            ! Exchange P(m) and P(n) and mark P(n) visited
            l = p(m)
            p(m) = p(n)
            p(n) = l
            x = .true.
          end if
        end do
      end do
    end do
  end subroutine Triangularize_s

end module Triangularize_m
