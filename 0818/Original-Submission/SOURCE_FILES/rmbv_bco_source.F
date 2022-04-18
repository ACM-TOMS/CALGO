      subroutine XRMBV_BCO(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      
      ierr = -1
      m = size(y)
      n = size(x)

      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
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

      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_H_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_H_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine XRMBV_BCO
