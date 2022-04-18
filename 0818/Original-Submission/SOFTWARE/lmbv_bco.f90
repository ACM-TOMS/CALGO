      module mod_lmbv_bco
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BCO'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface lmbv_bco
        module procedure ilmbv_bco
        module procedure slmbv_bco
        module procedure dlmbv_bco
        module procedure clmbv_bco
        module procedure zlmbv_bco
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilmbv_bco (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      bofs = -base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
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
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine ilmbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine slmbv_bco (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      bofs = -base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
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
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine slmbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine dlmbv_bco (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      bofs = -base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
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
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine dlmbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine clmbv_bco (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      bofs = -base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
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
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine clmbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine zlmbv_bco (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,bofs,i,mm,nn,nnz,nn_sq
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      bofs = -base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'t',type,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr=blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.m).or.&
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
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
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
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).lt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
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
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
             else if (mat%IA1(i).gt.mat%IA2(i)) then
                call block_Z_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,&
      y((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,store,ierr) 
                call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn, &
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            call block_T_mult_vec(mat%A((i-1)*nn_sq+1:i*nn_sq),&
      x((mat%IA1(i)+bofs)*nn+1:(mat%IA1(i)+bofs+1)*nn),nn,&
      y((mat%IA2(i)+bofs)*nn+1:(mat%IA2(i)+bofs+1)*nn),nn,store,ierr) 
         end do
         ierr = 0
      end if
      end subroutine zlmbv_bco 
! **********************************************************************
! **********************************************************************
      end module mod_lmbv_bco
