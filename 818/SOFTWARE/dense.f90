      module mod_dense_mat_algos
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : DENSE MATRIX ALGORITHMS FOR BLOCK SPARSE MATRICES
! **********************************************************************
      use properties
      implicit none
      interface block_mult_vec
        module procedure iblock_mult_vec
        module procedure sblock_mult_vec
        module procedure dblock_mult_vec
        module procedure cblock_mult_vec
        module procedure zblock_mult_vec
      end interface
      interface block_Z_mult_vec
        module procedure iblock_Z_mult_vec
        module procedure sblock_Z_mult_vec
        module procedure dblock_Z_mult_vec
        module procedure cblock_Z_mult_vec
        module procedure zblock_Z_mult_vec
      end interface
      interface block_T_mult_vec
        module procedure iblock_T_mult_vec
        module procedure sblock_T_mult_vec
        module procedure dblock_T_mult_vec
        module procedure cblock_T_mult_vec
        module procedure zblock_T_mult_vec
      end interface
      interface block_H_mult_vec
        module procedure iblock_H_mult_vec
        module procedure sblock_H_mult_vec
        module procedure dblock_H_mult_vec
        module procedure cblock_H_mult_vec
        module procedure zblock_H_mult_vec
      end interface
      interface invert_left_lower
        module procedure iinvert_left_lower
        module procedure sinvert_left_lower
        module procedure dinvert_left_lower
        module procedure cinvert_left_lower
        module procedure zinvert_left_lower
      end interface
      interface invert_T_left_lower
        module procedure iinvert_T_left_lower
        module procedure sinvert_T_left_lower
        module procedure dinvert_T_left_lower
        module procedure cinvert_T_left_lower
        module procedure zinvert_T_left_lower
      end interface
      interface invert_right_upper
        module procedure iinvert_right_upper
        module procedure sinvert_right_upper
        module procedure dinvert_right_upper
        module procedure cinvert_right_upper
        module procedure zinvert_right_upper
      end interface
      interface invert_T_right_upper
        module procedure iinvert_T_right_upper
        module procedure sinvert_T_right_upper
        module procedure dinvert_T_right_upper
        module procedure cinvert_T_right_upper
        module procedure zinvert_T_right_upper
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iblock_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iblock_r_mult_vec (A,x,n,y,m,ierr)
      else
         call iblock_l_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine iblock_mult_vec 
! ***
      subroutine iblock_Z_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iblock_r_mult_vec (A, (x),n,y,m,ierr)
      else
         call iblock_l_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine iblock_Z_mult_vec 
! ***
      subroutine iblock_T_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iblock_l_mult_vec (A,x,n,y,m,ierr)
      else
         call iblock_r_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine iblock_T_mult_vec 
! ***
      subroutine iblock_H_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iblock_l_mult_vec (A, (x),n,y,m,ierr)
      else
         call iblock_r_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine iblock_H_mult_vec 
! ***
      subroutine iinvert_left_lower (A,x,n,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iinvert_r_left_lower (A,x,n,ierr)
      else
         call iinvert_l_right_upper (A,x,n,ierr)
      end if
      end subroutine iinvert_left_lower 
! ***
      subroutine iinvert_T_left_lower (A,x,n,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iinvert_l_left_lower (A,x,n,ierr)
      else
         call iinvert_r_right_upper (A,x,n,ierr)
      end if
      end subroutine iinvert_T_left_lower 
! ***
      subroutine iinvert_right_upper (A,x,n,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iinvert_r_right_upper (A,x,n,ierr)
      else
         call iinvert_l_left_lower (A,x,n,ierr)
      end if
      end subroutine iinvert_right_upper 
! ***
      subroutine iinvert_T_right_upper (A,x,n,store,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call iinvert_l_right_upper (A,x,n,ierr)
      else
         call iinvert_r_left_lower (A,x,n,ierr)
      end if
      end subroutine iinvert_T_right_upper 
! ***
! ***
! ***
      subroutine iblock_r_mult_vec (A,x,n,y,m,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(size(y).ne.m).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y = y + A((j-1)*m+1:j*m) * x(j)
      end do
      ierr = 0
      end subroutine iblock_r_mult_vec 
! ***
      subroutine iblock_l_mult_vec (A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      integer , dimension(:), intent(in) :: A,x
      integer , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.m).or.(size(y).ne.n).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y(j) = y(j) + dot_product(A((j-1)*m+1:j*m),x)
      end do
      ierr = 0
      end subroutine iblock_l_mult_vec 
! ***
      subroutine iinvert_r_left_lower (A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(j+1:n) = x(j+1:n) - x(j) * a((j-1)*n+j+1:j*n)
      end do
      ierr = 0
      end subroutine iinvert_r_left_lower 
! ***
      subroutine iinvert_l_left_lower (A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         x(j) = x(j) - dot_product(a((j-1)*n+j+1:j*n),x(j+1:n))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine iinvert_l_left_lower 
! ***
      subroutine iinvert_r_right_upper (A,x,n,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(1:j-1) = x(1:j-1) - x(j) * a((j-1)*n+1:(j-1)*n+j-1)
      end do
      ierr = 0
      end subroutine iinvert_r_right_upper 
! ***
      subroutine iinvert_l_right_upper (A,x,n,ierr)
      implicit none
      integer , dimension(:), intent(in) :: A
      integer , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         x(j) = x(j) - dot_product(a((j-1)*n+1:(j-1)*n+j-1),x(1:j-1))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine iinvert_l_right_upper 
! **********************************************************************
! **********************************************************************
      subroutine sblock_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sblock_r_mult_vec (A,x,n,y,m,ierr)
      else
         call sblock_l_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine sblock_mult_vec 
! ***
      subroutine sblock_Z_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sblock_r_mult_vec (A, (x),n,y,m,ierr)
      else
         call sblock_l_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine sblock_Z_mult_vec 
! ***
      subroutine sblock_T_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sblock_l_mult_vec (A,x,n,y,m,ierr)
      else
         call sblock_r_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine sblock_T_mult_vec 
! ***
      subroutine sblock_H_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sblock_l_mult_vec (A, (x),n,y,m,ierr)
      else
         call sblock_r_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine sblock_H_mult_vec 
! ***
      subroutine sinvert_left_lower (A,x,n,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sinvert_r_left_lower (A,x,n,ierr)
      else
         call sinvert_l_right_upper (A,x,n,ierr)
      end if
      end subroutine sinvert_left_lower 
! ***
      subroutine sinvert_T_left_lower (A,x,n,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sinvert_l_left_lower (A,x,n,ierr)
      else
         call sinvert_r_right_upper (A,x,n,ierr)
      end if
      end subroutine sinvert_T_left_lower 
! ***
      subroutine sinvert_right_upper (A,x,n,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sinvert_r_right_upper (A,x,n,ierr)
      else
         call sinvert_l_left_lower (A,x,n,ierr)
      end if
      end subroutine sinvert_right_upper 
! ***
      subroutine sinvert_T_right_upper (A,x,n,store,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call sinvert_l_right_upper (A,x,n,ierr)
      else
         call sinvert_r_left_lower (A,x,n,ierr)
      end if
      end subroutine sinvert_T_right_upper 
! ***
! ***
! ***
      subroutine sblock_r_mult_vec (A,x,n,y,m,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(size(y).ne.m).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y = y + A((j-1)*m+1:j*m) * x(j)
      end do
      ierr = 0
      end subroutine sblock_r_mult_vec 
! ***
      subroutine sblock_l_mult_vec (A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      real(KIND=sp) , dimension(:), intent(in) :: A,x
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.m).or.(size(y).ne.n).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y(j) = y(j) + dot_product(A((j-1)*m+1:j*m),x)
      end do
      ierr = 0
      end subroutine sblock_l_mult_vec 
! ***
      subroutine sinvert_r_left_lower (A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(j+1:n) = x(j+1:n) - x(j) * a((j-1)*n+j+1:j*n)
      end do
      ierr = 0
      end subroutine sinvert_r_left_lower 
! ***
      subroutine sinvert_l_left_lower (A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         x(j) = x(j) - dot_product(a((j-1)*n+j+1:j*n),x(j+1:n))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine sinvert_l_left_lower 
! ***
      subroutine sinvert_r_right_upper (A,x,n,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(1:j-1) = x(1:j-1) - x(j) * a((j-1)*n+1:(j-1)*n+j-1)
      end do
      ierr = 0
      end subroutine sinvert_r_right_upper 
! ***
      subroutine sinvert_l_right_upper (A,x,n,ierr)
      implicit none
      real(KIND=sp) , dimension(:), intent(in) :: A
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         x(j) = x(j) - dot_product(a((j-1)*n+1:(j-1)*n+j-1),x(1:j-1))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine sinvert_l_right_upper 
! **********************************************************************
! **********************************************************************
      subroutine dblock_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dblock_r_mult_vec (A,x,n,y,m,ierr)
      else
         call dblock_l_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine dblock_mult_vec 
! ***
      subroutine dblock_Z_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dblock_r_mult_vec (A, (x),n,y,m,ierr)
      else
         call dblock_l_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine dblock_Z_mult_vec 
! ***
      subroutine dblock_T_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dblock_l_mult_vec (A,x,n,y,m,ierr)
      else
         call dblock_r_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine dblock_T_mult_vec 
! ***
      subroutine dblock_H_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dblock_l_mult_vec (A, (x),n,y,m,ierr)
      else
         call dblock_r_mult_vec (A, (x),n,y,m,ierr)
      end if
      y=  (y)
      end subroutine dblock_H_mult_vec 
! ***
      subroutine dinvert_left_lower (A,x,n,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dinvert_r_left_lower (A,x,n,ierr)
      else
         call dinvert_l_right_upper (A,x,n,ierr)
      end if
      end subroutine dinvert_left_lower 
! ***
      subroutine dinvert_T_left_lower (A,x,n,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dinvert_l_left_lower (A,x,n,ierr)
      else
         call dinvert_r_right_upper (A,x,n,ierr)
      end if
      end subroutine dinvert_T_left_lower 
! ***
      subroutine dinvert_right_upper (A,x,n,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dinvert_r_right_upper (A,x,n,ierr)
      else
         call dinvert_l_left_lower (A,x,n,ierr)
      end if
      end subroutine dinvert_right_upper 
! ***
      subroutine dinvert_T_right_upper (A,x,n,store,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call dinvert_l_right_upper (A,x,n,ierr)
      else
         call dinvert_r_left_lower (A,x,n,ierr)
      end if
      end subroutine dinvert_T_right_upper 
! ***
! ***
! ***
      subroutine dblock_r_mult_vec (A,x,n,y,m,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(size(y).ne.m).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y = y + A((j-1)*m+1:j*m) * x(j)
      end do
      ierr = 0
      end subroutine dblock_r_mult_vec 
! ***
      subroutine dblock_l_mult_vec (A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      real(KIND=dp) , dimension(:), intent(in) :: A,x
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.m).or.(size(y).ne.n).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y(j) = y(j) + dot_product(A((j-1)*m+1:j*m),x)
      end do
      ierr = 0
      end subroutine dblock_l_mult_vec 
! ***
      subroutine dinvert_r_left_lower (A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(j+1:n) = x(j+1:n) - x(j) * a((j-1)*n+j+1:j*n)
      end do
      ierr = 0
      end subroutine dinvert_r_left_lower 
! ***
      subroutine dinvert_l_left_lower (A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         x(j) = x(j) - dot_product(a((j-1)*n+j+1:j*n),x(j+1:n))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine dinvert_l_left_lower 
! ***
      subroutine dinvert_r_right_upper (A,x,n,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(1:j-1) = x(1:j-1) - x(j) * a((j-1)*n+1:(j-1)*n+j-1)
      end do
      ierr = 0
      end subroutine dinvert_r_right_upper 
! ***
      subroutine dinvert_l_right_upper (A,x,n,ierr)
      implicit none
      real(KIND=dp) , dimension(:), intent(in) :: A
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         x(j) = x(j) - dot_product(a((j-1)*n+1:(j-1)*n+j-1),x(1:j-1))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine dinvert_l_right_upper 
! **********************************************************************
! **********************************************************************
      subroutine cblock_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cblock_r_mult_vec (A,x,n,y,m,ierr)
      else
         call cblock_l_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine cblock_mult_vec 
! ***
      subroutine cblock_Z_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cblock_r_mult_vec (A,conjg (x),n,y,m,ierr)
      else
         call cblock_l_mult_vec (A,conjg (x),n,y,m,ierr)
      end if
      y= conjg (y)
      end subroutine cblock_Z_mult_vec 
! ***
      subroutine cblock_T_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cblock_l_mult_vec (A,x,n,y,m,ierr)
      else
         call cblock_r_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine cblock_T_mult_vec 
! ***
      subroutine cblock_H_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cblock_l_mult_vec (A,conjg (x),n,y,m,ierr)
      else
         call cblock_r_mult_vec (A,conjg (x),n,y,m,ierr)
      end if
      y= conjg (y)
      end subroutine cblock_H_mult_vec 
! ***
      subroutine cinvert_left_lower (A,x,n,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cinvert_r_left_lower (A,x,n,ierr)
      else
         call cinvert_l_right_upper (A,x,n,ierr)
      end if
      end subroutine cinvert_left_lower 
! ***
      subroutine cinvert_T_left_lower (A,x,n,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cinvert_l_left_lower (A,x,n,ierr)
      else
         call cinvert_r_right_upper (A,x,n,ierr)
      end if
      end subroutine cinvert_T_left_lower 
! ***
      subroutine cinvert_right_upper (A,x,n,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cinvert_r_right_upper (A,x,n,ierr)
      else
         call cinvert_l_left_lower (A,x,n,ierr)
      end if
      end subroutine cinvert_right_upper 
! ***
      subroutine cinvert_T_right_upper (A,x,n,store,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call cinvert_l_right_upper (A,x,n,ierr)
      else
         call cinvert_r_left_lower (A,x,n,ierr)
      end if
      end subroutine cinvert_T_right_upper 
! ***
! ***
! ***
      subroutine cblock_r_mult_vec (A,x,n,y,m,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(size(y).ne.m).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y = y + A((j-1)*m+1:j*m) * x(j)
      end do
      ierr = 0
      end subroutine cblock_r_mult_vec 
! ***
      subroutine cblock_l_mult_vec (A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      complex(KIND=sp) , dimension(:), intent(in) :: A,x
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.m).or.(size(y).ne.n).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y(j) = y(j) + dot_product(A((j-1)*m+1:j*m),x)
      end do
      ierr = 0
      end subroutine cblock_l_mult_vec 
! ***
      subroutine cinvert_r_left_lower (A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(j+1:n) = x(j+1:n) - x(j) * a((j-1)*n+j+1:j*n)
      end do
      ierr = 0
      end subroutine cinvert_r_left_lower 
! ***
      subroutine cinvert_l_left_lower (A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         x(j) = x(j) - dot_product(a((j-1)*n+j+1:j*n),x(j+1:n))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine cinvert_l_left_lower 
! ***
      subroutine cinvert_r_right_upper (A,x,n,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(1:j-1) = x(1:j-1) - x(j) * a((j-1)*n+1:(j-1)*n+j-1)
      end do
      ierr = 0
      end subroutine cinvert_r_right_upper 
! ***
      subroutine cinvert_l_right_upper (A,x,n,ierr)
      implicit none
      complex(KIND=sp) , dimension(:), intent(in) :: A
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         x(j) = x(j) - dot_product(a((j-1)*n+1:(j-1)*n+j-1),x(1:j-1))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine cinvert_l_right_upper 
! **********************************************************************
! **********************************************************************
      subroutine zblock_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zblock_r_mult_vec (A,x,n,y,m,ierr)
      else
         call zblock_l_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine zblock_mult_vec 
! ***
      subroutine zblock_Z_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zblock_r_mult_vec (A,conjg (x),n,y,m,ierr)
      else
         call zblock_l_mult_vec (A,conjg (x),n,y,m,ierr)
      end if
      y= conjg (y)
      end subroutine zblock_Z_mult_vec 
! ***
      subroutine zblock_T_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zblock_l_mult_vec (A,x,n,y,m,ierr)
      else
         call zblock_r_mult_vec (A,x,n,y,m,ierr)
      end if
      end subroutine zblock_T_mult_vec 
! ***
      subroutine zblock_H_mult_vec (A,x,n,y,m,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zblock_l_mult_vec (A,conjg (x),n,y,m,ierr)
      else
         call zblock_r_mult_vec (A,conjg (x),n,y,m,ierr)
      end if
      y= conjg (y)
      end subroutine zblock_H_mult_vec 
! ***
      subroutine zinvert_left_lower (A,x,n,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zinvert_r_left_lower (A,x,n,ierr)
      else
         call zinvert_l_right_upper (A,x,n,ierr)
      end if
      end subroutine zinvert_left_lower 
! ***
      subroutine zinvert_T_left_lower (A,x,n,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zinvert_l_left_lower (A,x,n,ierr)
      else
         call zinvert_r_right_upper (A,x,n,ierr)
      end if
      end subroutine zinvert_T_left_lower 
! ***
      subroutine zinvert_right_upper (A,x,n,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zinvert_r_right_upper (A,x,n,ierr)
      else
         call zinvert_l_left_lower (A,x,n,ierr)
      end if
      end subroutine zinvert_right_upper 
! ***
      subroutine zinvert_T_right_upper (A,x,n,store,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call zinvert_l_right_upper (A,x,n,ierr)
      else
         call zinvert_r_left_lower (A,x,n,ierr)
      end if
      end subroutine zinvert_T_right_upper 
! ***
! ***
! ***
      subroutine zblock_r_mult_vec (A,x,n,y,m,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(size(y).ne.m).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y = y + A((j-1)*m+1:j*m) * x(j)
      end do
      ierr = 0
      end subroutine zblock_r_mult_vec 
! ***
      subroutine zblock_l_mult_vec (A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      complex(KIND=dp) , dimension(:), intent(in) :: A,x
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.m).or.(size(y).ne.n).or.(n*m.ne.size(A))) then
         return
      end if
      do j=1,n
         y(j) = y(j) + dot_product(A((j-1)*m+1:j*m),x)
      end do
      ierr = 0
      end subroutine zblock_l_mult_vec 
! ***
      subroutine zinvert_r_left_lower (A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(j+1:n) = x(j+1:n) - x(j) * a((j-1)*n+j+1:j*n)
      end do
      ierr = 0
      end subroutine zinvert_r_left_lower 
! ***
      subroutine zinvert_l_left_lower (A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         x(j) = x(j) - dot_product(a((j-1)*n+j+1:j*n),x(j+1:n))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine zinvert_l_left_lower 
! ***
      subroutine zinvert_r_right_upper (A,x,n,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = n,1,-1
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
         x(1:j-1) = x(1:j-1) - x(j) * a((j-1)*n+1:(j-1)*n+j-1)
      end do
      ierr = 0
      end subroutine zinvert_r_right_upper 
! ***
      subroutine zinvert_l_right_upper (A,x,n,ierr)
      implicit none
      complex(KIND=dp) , dimension(:), intent(in) :: A
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: j
      ierr = -1
      if((size(x).ne.n).or.(n*n.ne.size(A))) then
         return
      end if
      do j = 1,n
         x(j) = x(j) - dot_product(a((j-1)*n+1:(j-1)*n+j-1),x(1:j-1))
         if (a((j-1)*n+j).ne.0) then
            x(j) = x(j)/a((j-1)*n+j)
         else
            ierr=blas_error_singtria 
            return
         end if
      end do
      ierr = 0
      end subroutine zinvert_l_right_upper 
! **********************************************************************
! **********************************************************************
      end module mod_dense_mat_algos
