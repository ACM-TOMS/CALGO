      subroutine XBLOCK_MULT(A,x,n,y,m,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XBLOCK_R_MULT(A,x,n,y,m,ierr)
      else
         call XBLOCK_L_MULT(A,x,n,y,m,ierr)
      end if
      end subroutine XBLOCK_MULT
! ***
      subroutine XBLOCK_Z_MULT(A,x,n,y,m,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XBLOCK_R_MULT(A,CONJUG(x),n,y,m,ierr)
      else
         call XBLOCK_L_MULT(A,CONJUG(x),n,y,m,ierr)
      end if
      y=CONJUG(y)
      end subroutine XBLOCK_Z_MULT
! ***
      subroutine XBLOCK_T_MULT(A,x,n,y,m,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XBLOCK_L_MULT(A,x,n,y,m,ierr)
      else
         call XBLOCK_R_MULT(A,x,n,y,m,ierr)
      end if
      end subroutine XBLOCK_T_MULT
! ***
      subroutine XBLOCK_H_MULT(A,x,n,y,m,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
      integer, intent(in) :: m,n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XBLOCK_L_MULT(A,CONJUG(x),n,y,m,ierr)
      else
         call XBLOCK_R_MULT(A,CONJUG(x),n,y,m,ierr)
      end if
      y=CONJUG(y)
      end subroutine XBLOCK_H_MULT
! ***
      subroutine XINV_LL(A,x,n,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XINV_R_LL(A,x,n,ierr)
      else
         call XINV_L_RU(A,x,n,ierr)
      end if
      end subroutine XINV_LL
! ***
      subroutine XINV_T_LL(A,x,n,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XINV_L_LL(A,x,n,ierr)
      else
         call XINV_R_RU(A,x,n,ierr)
      end if
      end subroutine XINV_T_LL
! ***
      subroutine XINV_RU(A,x,n,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XINV_R_RU(A,x,n,ierr)
      else
         call XINV_L_LL(A,x,n,ierr)
      end if
      end subroutine XINV_RU
! ***
      subroutine XINV_T_RU(A,x,n,store,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(in) :: n
      character, intent(in) :: store
      integer, intent(out) :: ierr
      if (store.eq.'C') then
         call XINV_L_RU(A,x,n,ierr)
      else
         call XINV_R_LL(A,x,n,ierr)
      end if
      end subroutine XINV_T_RU
! ***
! ***
! ***
      subroutine XBLOCK_R_MULT(A,x,n,y,m,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
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
      end subroutine XBLOCK_R_MULT
! ***
      subroutine XBLOCK_L_MULT(A,x,m,y,n,ierr)
      implicit none
      intrinsic dot_product
      DCOMPLEX, dimension(:), intent(in) :: A,x
      DCOMPLEX, dimension(:), intent(inout) :: y
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
      end subroutine XBLOCK_L_MULT
! ***
      subroutine XINV_R_LL(A,x,n,ierr)
      !left_lower, stored column-wise
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
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
      end subroutine XINV_R_LL
! ***
      subroutine XINV_L_LL(A,x,n,ierr)
      !left_lower, stored row-wise
      implicit none
      intrinsic dot_product
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
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
      end subroutine XINV_L_LL
! ***
      subroutine XINV_R_RU(A,x,n,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
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
      end subroutine XINV_R_RU
! ***
      subroutine XINV_L_RU(A,x,n,ierr)
      implicit none
      DCOMPLEX, dimension(:), intent(in) :: A
      DCOMPLEX, dimension(:), intent(inout) :: x
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
      end subroutine XINV_L_RU
