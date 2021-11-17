      subroutine XLSBV_BCO(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      DCOMPLEX, allocatable, dimension(:) :: y !extra stor.!

      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
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

      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if

      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = ZZERO

      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if

      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do

      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
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
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
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
      call remove_hash(ierr)

      end subroutine XLSBV_BCO
