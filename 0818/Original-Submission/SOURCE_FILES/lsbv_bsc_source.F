      subroutine XLSBV_BSC(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      DCOMPLEX, allocatable, dimension(:) :: y

      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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

      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((j-1)*nn+1:j*nn) = x((j-1)*nn+1:j*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while(pntr.lt.mat%pe(j))
                  if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                     dd = pntr
                  else 
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((j-1)*nn+1:j*nn) = x((j-1)*nn+1:j*nn) - y
                  end if
                  pntr = pntr + 1
               end do
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,store,ierr)
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
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((j-1)*nn+1:j*nn) = x((j-1)*nn+1:j*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while(pntr.lt.mat%pe(j))
                  if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                     dd = pntr
                  else 
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((j-1)*nn+1:j*nn) = x((j-1)*nn+1:j*nn) - y
                  end if
                  pntr = pntr + 1
               end do
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,store,ierr)
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

      end subroutine XLSBV_BSC
