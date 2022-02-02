      module mod_rsbv_csr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface rsbv_csr
        module procedure irsbv_csr
        module procedure srsbv_csr
        module procedure drsbv_csr
        module procedure crsbv_csr
        module procedure zrsbv_csr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irsbv_csr (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      integer  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0  
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
      end subroutine irsbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine srsbv_csr (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      real(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0.0e0  
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0.0e0 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
      end subroutine srsbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine drsbv_csr (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      real(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0.0d0  
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = 0.0d0 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
      end subroutine drsbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine crsbv_csr (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      complex(KIND=sp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = (0.0e0, 0.0e0)  
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = (0.0e0, 0.0e0) 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
      end subroutine crsbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine zrsbv_csr (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      complex(KIND=dp)  :: de
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = (0.0d0, 0.0d0)  
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = (0.0d0, 0.0d0) 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
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
      end subroutine zrsbv_csr 
! **********************************************************************
! **********************************************************************
      end module mod_rsbv_csr
