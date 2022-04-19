      subroutine XLSBV_COO(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
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
                  if (mat%A(hash(i)%val_pos).ne.ZZERO) then
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
                  if (mat%A(hash(i)%val_pos).ne.ZZERO) then
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

      end subroutine XLSBV_COO
