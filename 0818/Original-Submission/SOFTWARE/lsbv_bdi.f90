      module mod_lsbv_bdi
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BDI'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface lsbv_bdi
        module procedure ilsbv_bdi
        module procedure slsbv_bdi
        module procedure dlsbv_bdi
        module procedure clsbv_bdi
        module procedure zlsbv_bdi
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilsbv_bdi (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,mm,nn,blda,nbdiag,dd,nn_sq
      character :: diag,part,store
      integer , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn * nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0 
      if (part.eq.'L') then
         do i=blda,1,-1
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((-mat%IA1(j).lt.blda-i+1).and.&
                     (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((-mat%IA1(j).lt.blda-i+1).and.&
                          (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_left_lower(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i=1,blda
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((mat%IA1(j).lt.i).and.&
                    (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((mat%IA1(j).lt.i).and.&
                          (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_right_upper(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine ilsbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine slsbv_bdi (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,mm,nn,blda,nbdiag,dd,nn_sq
      character :: diag,part,store
      real(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn * nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0e0 
      if (part.eq.'L') then
         do i=blda,1,-1
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((-mat%IA1(j).lt.blda-i+1).and.&
                     (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((-mat%IA1(j).lt.blda-i+1).and.&
                          (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_left_lower(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i=1,blda
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((mat%IA1(j).lt.i).and.&
                    (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((mat%IA1(j).lt.i).and.&
                          (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_right_upper(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine slsbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine dlsbv_bdi (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,mm,nn,blda,nbdiag,dd,nn_sq
      character :: diag,part,store
      real(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn * nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0d0 
      if (part.eq.'L') then
         do i=blda,1,-1
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((-mat%IA1(j).lt.blda-i+1).and.&
                     (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((-mat%IA1(j).lt.blda-i+1).and.&
                          (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_left_lower(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i=1,blda
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((mat%IA1(j).lt.i).and.&
                    (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((mat%IA1(j).lt.i).and.&
                          (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_right_upper(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine dlsbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine clsbv_bdi (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,mm,nn,blda,nbdiag,dd,nn_sq
      character :: diag,part,store
      complex(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn * nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0e0, 0.0e0) 
      if (part.eq.'L') then
         do i=blda,1,-1
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((-mat%IA1(j).lt.blda-i+1).and.&
                     (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((-mat%IA1(j).lt.blda-i+1).and.&
                          (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_left_lower(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i=1,blda
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((mat%IA1(j).lt.i).and.&
                    (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((mat%IA1(j).lt.i).and.&
                          (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_right_upper(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine clsbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine zlsbv_bdi (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,mm,nn,blda,nbdiag,dd,nn_sq
      character :: diag,part,store
      complex(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn * nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0d0, 0.0d0) 
      if (part.eq.'L') then
         do i=blda,1,-1
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((-mat%IA1(j).lt.blda-i+1).and.&
                     (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((-mat%IA1(j).lt.blda-i+1).and.&
                          (mat%IA1(j).lt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_left_lower(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i=1,blda
            if (diag.eq.'U') then
               do j = 1,nbdiag
                  if((mat%IA1(j).lt.i).and.&
                    (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                  x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
            else
               dd = -1
               do j = 1,nbdiag
                  if (mat%IA1(j).eq.0) then
                     dd = j
                  else if((mat%IA1(j).lt.i).and.&
                          (mat%IA1(j).gt.0)) then
                     call block_T_mult_vec(&
      mat%A((blda*(j-1)-mat%IA1(j)+i-1)*nn_sq+1:&
            (blda*(j-1)-mat%IA1(j)+i)*nn_sq),&
      x((i-mat%IA1(j)-1)*nn+1:(i-mat%IA1(j))*nn),nn,y,nn,store,ierr)
                     x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
                  end if
               end do
               if (dd.ne.-1) then 
                  call invert_T_right_upper(&
      mat%A((blda*(dd-1)+i-1)*nn_sq+1:(blda*(dd-1)+i)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine zlsbv_bdi 
! **********************************************************************
! **********************************************************************
      end module mod_lsbv_bdi
