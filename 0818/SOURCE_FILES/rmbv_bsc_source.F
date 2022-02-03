      subroutine XRMBV_BSC(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,bofs,j,pntr,mm,nn,mb,nb,nn_sq
      character :: diag,type,part,store
      
      ierr = -1
      m = size(y)
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
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if

      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if

      y = ZZERO
      nn_sq = nn*nn

      if (diag.eq.'U') then !process unstored diagonal
         if (m.eq.n) then
            y = x
         else
            ierr = blas_error_param
            return
         end if
      end if

      if ((type.eq.'S').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do j = 1, nb
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine XRMBV_BSC
