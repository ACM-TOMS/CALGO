      module mod_lsbv_csc
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'CSC'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************      
      use representation_of_data
      use properties
      implicit none
      interface lsbv_csc
        module procedure ilsbv_csc
        module procedure slsbv_csc
        module procedure dlsbv_csc
        module procedure clsbv_csc
        module procedure zlsbv_csc
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilsbv_csc (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      integer  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      ierr = -1
      if (part.eq.'L') then
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) & 
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine ilsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine slsbv_csc (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      real(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      ierr = -1
      if (part.eq.'L') then
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0.0e0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) & 
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0.0e0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0.0e0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0.0e0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine slsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine dlsbv_csc (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      real(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      ierr = -1
      if (part.eq.'L') then
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0.0d0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) & 
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0.0d0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = 0.0d0 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. 0.0d0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine dlsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine clsbv_csc (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      complex(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      ierr = -1
      if (part.eq.'L') then
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = (0.0e0, 0.0e0) 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) & 
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. (0.0e0, 0.0e0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = (0.0e0, 0.0e0) 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. (0.0e0, 0.0e0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine clsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine zlsbv_csc (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      complex(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      ierr = -1
      if (part.eq.'L') then
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = (0.0d0, 0.0d0) 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) & 
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. (0.0d0, 0.0d0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  pntr = pntr + 1
               end do
            else
               de = (0.0d0, 0.0d0) 
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.eq.i) then
                     de = mat%A(pntr + ofs)
                  else
                     x(i) = x(i) &
                    - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq. (0.0d0, 0.0d0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if
      end subroutine zlsbv_csc 
! **********************************************************************
! **********************************************************************
      end module mod_lsbv_csc
