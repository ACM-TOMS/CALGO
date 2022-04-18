      module mod_lmbv_coo
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH TRANSPOSE IN 'COO'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T = x^TA
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface lmbv_coo
        module procedure ilmbv_coo
        module procedure slmbv_coo
        module procedure dlmbv_coo
        module procedure clmbv_coo
        module procedure zlmbv_coo
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilmbv_coo (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,nnz,base,i,ofs
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
      y = 0 
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                               + mat%A(i) * x(mat%IA1(i)+ofs)
         end do
         ierr = 0
      end if 
      end subroutine ilmbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine slmbv_coo (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,nnz,base,i,ofs
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
      y = 0.0e0 
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                               + mat%A(i) * x(mat%IA1(i)+ofs)
         end do
         ierr = 0
      end if 
      end subroutine slmbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine dlmbv_coo (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,nnz,base,i,ofs
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
      y = 0.0d0 
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           +  (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                               + mat%A(i) * x(mat%IA1(i)+ofs)
         end do
         ierr = 0
      end if 
      end subroutine dlmbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine clmbv_coo (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,nnz,base,i,ofs
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
      y = (0.0e0, 0.0e0) 
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           + conjg (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           + conjg (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                               + mat%A(i) * x(mat%IA1(i)+ofs)
         end do
         ierr = 0
      end if 
      end subroutine clmbv_coo 
! **********************************************************************
! **********************************************************************
      subroutine zlmbv_coo (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,nnz,base,i,ofs
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'COO').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
      y = (0.0d0, 0.0d0) 
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
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                                   + mat%A(i) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).gt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           + conjg (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         else
           do i = 1, nnz
             if (mat%IA1(i).eq.mat%IA2(i)) then
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             else if (mat%IA1(i).lt.mat%IA2(i)) then 
                y(mat%IA1(i)+ofs) = y(mat%IA1(i)+ofs) &
                           + conjg (mat%A(i)) * x(mat%IA2(i)+ofs)
                y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                                   + mat%A(i) * x(mat%IA1(i)+ofs)
             end if
           end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, nnz
            y(mat%IA2(i)+ofs) = y(mat%IA2(i)+ofs) &
                               + mat%A(i) * x(mat%IA1(i)+ofs)
         end do
         ierr = 0
      end if 
      end subroutine zlmbv_coo 
! **********************************************************************
! **********************************************************************
      end module mod_lmbv_coo
