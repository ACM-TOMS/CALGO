      module mod_rmbv_csr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface rmbv_csr
        module procedure irmbv_csr
        module procedure srmbv_csr
        module procedure drmbv_csr
        module procedure crmbv_csr
        module procedure zrmbv_csr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irmbv_csr (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) =  &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                    y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else
         do i = 1, mat%M
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
              y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine irmbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine srmbv_csr (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) =  &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                    y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else
         do i = 1, mat%M
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
              y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine srmbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine drmbv_csr (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) =  &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                    y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           +  (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else
         do i = 1, mat%M
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
              y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine drmbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine crmbv_csr (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) =  &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                    y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           + conjg (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           + conjg (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else
         do i = 1, mat%M
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
              y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine crmbv_csr 
! **********************************************************************
! **********************************************************************
      subroutine zrmbv_csr (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) =  &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                    y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           + conjg (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                   y(mat%IA1(pntr+ofs)+ofs)=y(mat%IA1(pntr+ofs)+ofs) &
                           + conjg (mat%A(pntr + ofs)) * x(i) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else
         do i = 1, mat%M
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
              y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine zrmbv_csr 
! **********************************************************************
! **********************************************************************
      end module mod_rmbv_csr
