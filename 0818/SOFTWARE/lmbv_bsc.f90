      module mod_lmbv_bsc
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSC'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface lmbv_bsc
        module procedure ilmbv_bsc
        module procedure slmbv_bsc
        module procedure dlmbv_bsc
        module procedure clmbv_bsc
        module procedure zlmbv_bsc
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilmbv_bsc (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
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
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      y = 0 
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
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
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
               call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine ilmbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine slmbv_bsc (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
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
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      y = 0.0e0 
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
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
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
               call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine slmbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine dlmbv_bsc (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
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
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      y = 0.0d0 
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
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
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
               call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine dlmbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine clmbv_bsc (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
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
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      y = (0.0e0, 0.0e0) 
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
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
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
               call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine clmbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine zlmbv_bsc (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
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
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
         (mm.ne.nn)) then
         ierr = blas_error_param
         return
      end if
      y = (0.0d0, 0.0d0) 
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
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do j = 1, nb
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(j.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                  else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                     call block_Z_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),&
      nn,y((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+&
      bofs+1)*nn),nn,store,ierr) 
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
               call block_T_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y((j-1)*nn+1:j*nn),nn,store,ierr) 
                pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine zlmbv_bsc 
! **********************************************************************
! **********************************************************************
      end module mod_lmbv_bsc
