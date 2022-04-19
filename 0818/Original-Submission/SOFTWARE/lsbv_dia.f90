      module mod_lsbv_dia
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'DIA'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface lsbv_dia
        module procedure ilsbv_dia
        module procedure slsbv_dia
        module procedure dlsbv_dia
        module procedure clsbv_dia
        module procedure zlsbv_dia
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilsbv_dia (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      integer  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
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
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = 0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = 0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine ilsbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine slsbv_dia (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      real(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
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
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = 0.0e0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0.0e0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = 0.0e0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0.0e0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine slsbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine dlsbv_dia (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      real(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
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
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = 0.0d0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0.0d0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = 0.0d0 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. 0.0d0 ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine dlsbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine clsbv_dia (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      complex(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
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
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = (0.0e0, 0.0e0) 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. (0.0e0, 0.0e0) ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = (0.0e0, 0.0e0) 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. (0.0e0, 0.0e0) ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine clsbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine zlsbv_dia (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      complex(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
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
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ierr = -1
      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = (0.0d0, 0.0d0) 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. (0.0d0, 0.0d0) ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = (0.0d0, 0.0d0) 
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne. (0.0d0, 0.0d0) ) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine zlsbv_dia 
! **********************************************************************
! **********************************************************************
      end module mod_lsbv_dia
