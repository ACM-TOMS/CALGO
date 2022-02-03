      module mod_lmbv_csc
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSC'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface lmbv_csc
        module procedure ilmbv_csc
        module procedure slmbv_csc
        module procedure dlmbv_csc
        module procedure clmbv_csc
        module procedure zlmbv_csc
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilmbv_csc (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                            +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                      +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else 
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine ilmbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine slmbv_csc (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                            +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                      +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else 
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine slmbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine dlmbv_csc (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                            +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                      +  (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else 
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine dlmbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine clmbv_csc (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                            + conjg (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                      + conjg (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else 
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine clmbv_csc 
! **********************************************************************
! **********************************************************************
      subroutine zlmbv_csc (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr + ofs) + ofs) = &
               y(mat%IA1(pntr + ofs) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then 
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                            + conjg (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then  
               y(mat%IA1(pntr+ofs)+ofs)= y(mat%IA1(pntr+ofs)+ofs) &
                      + conjg (mat%A(pntr + ofs)) * x(j) 
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else 
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               y(j) = y(j) + &
                     mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine zlmbv_csc 
! **********************************************************************
! **********************************************************************
      end module mod_lmbv_csc
