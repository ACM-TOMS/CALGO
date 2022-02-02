      module mod_rsbv_csc
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSC'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface rsbv_csc
        module procedure irsbv_csc
        module procedure srsbv_csc
        module procedure drsbv_csc
        module procedure crsbv_csc
        module procedure zrsbv_csc
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irsbv_csc (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if
      end subroutine irsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine srsbv_csc (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0.0e0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0.0e0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if
      end subroutine srsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine drsbv_csc (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0.0d0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. 0.0d0 ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if
      end subroutine drsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine crsbv_csc (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. (0.0e0, 0.0e0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. (0.0e0, 0.0e0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if
      end subroutine crsbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine zrsbv_csc (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. (0.0d0, 0.0d0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq. (0.0d0, 0.0d0) ) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if
      end subroutine zrsbv_csc 
! **********************************************************************
! **********************************************************************
      end module mod_rsbv_csc
