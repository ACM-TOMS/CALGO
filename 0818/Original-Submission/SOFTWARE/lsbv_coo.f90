      module mod_lsbv_coo
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'COO'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties
      implicit none
      interface lsbv_coo
        module procedure ilsbv_coo
        module procedure slsbv_coo
        module procedure dlsbv_coo
        module procedure clsbv_coo
        module procedure zlsbv_coo
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilsbv_coo (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,nnz
      character :: diag,part
      type(capsule), pointer :: dummy
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
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
      call setup_hash(n,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,i,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = n,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1, n
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      end if
      call remove_hash(ierr)
      end subroutine ilsbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine slsbv_coo (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,nnz
      character :: diag,part
      type(capsule), pointer :: dummy
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
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
      call setup_hash(n,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,i,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = n,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0.0e0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1, n
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0.0e0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      end if
      call remove_hash(ierr)
      end subroutine slsbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine dlsbv_coo (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,nnz
      character :: diag,part
      type(capsule), pointer :: dummy
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
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
      call setup_hash(n,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,i,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = n,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0.0d0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1, n
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. 0.0d0 ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      end if
      call remove_hash(ierr)
      end subroutine dlsbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine clsbv_coo (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,nnz
      character :: diag,part
      type(capsule), pointer :: dummy
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
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
      call setup_hash(n,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,i,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = n,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. (0.0e0, 0.0e0) ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1, n
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. (0.0e0, 0.0e0) ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      end if
      call remove_hash(ierr)
      end subroutine clsbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine zlsbv_coo (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,nnz
      character :: diag,part
      type(capsule), pointer :: dummy
      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
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
      call setup_hash(n,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,i,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = n,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. (0.0d0, 0.0d0) ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      else 
         do i = 1, n
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               x(i) = x(i) - x(dummy%jndx) * mat%A(dummy%val_pos)
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  if (mat%A(hash(i)%val_pos).ne. (0.0d0, 0.0d0) ) then
                     x(i) = x(i)/mat%A(hash(i)%val_pos)
                  else
                     ierr = blas_error_singtria
                     return
                  end if
               end if
            end if
         end do
         ierr = 0
      end if
      call remove_hash(ierr)
      end subroutine zlsbv_coo 
! **********************************************************************
! **********************************************************************
      end module mod_lsbv_coo
