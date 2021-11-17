      module mod_rmbv_dia
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'DIA'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface rmbv_dia
        module procedure irmbv_dia
        module procedure srmbv_dia
        module procedure drmbv_dia
        module procedure crmbv_dia
        module procedure zrmbv_dia
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irmbv_dia (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j
      integer :: lda,ndiag,start_a,end_a,start_x,start_y
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
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
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + mat%A(start_a+j) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         +  (mat%A(start_a+j)) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
                  end_a = i*lda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     y(start_x +j) = y(start_x+j) &
                         +  (mat%A(start_a+j)) * x(start_y+j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,ndiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.mat%K-lda) then
               start_a = (i-1)*lda
               end_a = i*lda -mat%IA1(i)+mat%K-lda
            else if (mat%IA1(i).lt.-mat%M+lda) then
               start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
               end_a = i*lda               
            else
               start_a = (i-1)*lda
               end_a = i*lda
            end if
            j = 1
            do while((start_a + j).le.end_a)
               y(start_y+j) = y(start_y+j) &
                           + mat%A(start_a+j) * x(start_x+j)
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine irmbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine srmbv_dia (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j
      integer :: lda,ndiag,start_a,end_a,start_x,start_y
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
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
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + mat%A(start_a+j) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         +  (mat%A(start_a+j)) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
                  end_a = i*lda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     y(start_x +j) = y(start_x+j) &
                         +  (mat%A(start_a+j)) * x(start_y+j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,ndiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.mat%K-lda) then
               start_a = (i-1)*lda
               end_a = i*lda -mat%IA1(i)+mat%K-lda
            else if (mat%IA1(i).lt.-mat%M+lda) then
               start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
               end_a = i*lda               
            else
               start_a = (i-1)*lda
               end_a = i*lda
            end if
            j = 1
            do while((start_a + j).le.end_a)
               y(start_y+j) = y(start_y+j) &
                           + mat%A(start_a+j) * x(start_x+j)
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine srmbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine drmbv_dia (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j
      integer :: lda,ndiag,start_a,end_a,start_x,start_y
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
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
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + mat%A(start_a+j) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         +  (mat%A(start_a+j)) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
                  end_a = i*lda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     y(start_x +j) = y(start_x+j) &
                         +  (mat%A(start_a+j)) * x(start_y+j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,ndiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.mat%K-lda) then
               start_a = (i-1)*lda
               end_a = i*lda -mat%IA1(i)+mat%K-lda
            else if (mat%IA1(i).lt.-mat%M+lda) then
               start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
               end_a = i*lda               
            else
               start_a = (i-1)*lda
               end_a = i*lda
            end if
            j = 1
            do while((start_a + j).le.end_a)
               y(start_y+j) = y(start_y+j) &
                           + mat%A(start_a+j) * x(start_x+j)
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine drmbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine crmbv_dia (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j
      integer :: lda,ndiag,start_a,end_a,start_x,start_y
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
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
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + mat%A(start_a+j) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + conjg (mat%A(start_a+j)) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
                  end_a = i*lda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     y(start_x +j) = y(start_x+j) &
                         + conjg (mat%A(start_a+j)) * x(start_y+j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,ndiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.mat%K-lda) then
               start_a = (i-1)*lda
               end_a = i*lda -mat%IA1(i)+mat%K-lda
            else if (mat%IA1(i).lt.-mat%M+lda) then
               start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
               end_a = i*lda               
            else
               start_a = (i-1)*lda
               end_a = i*lda
            end if
            j = 1
            do while((start_a + j).le.end_a)
               y(start_y+j) = y(start_y+j) &
                           + mat%A(start_a+j) * x(start_x+j)
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine crmbv_dia 
! **********************************************************************
! **********************************************************************
      subroutine zrmbv_dia (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j
      integer :: lda,ndiag,start_a,end_a,start_x,start_y
      character :: diag,type,part
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'d',lda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',ndiag,ierr)
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
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + mat%A(start_a+j) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'U') then
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda -mat%IA1(i)+mat%K-lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y +j) = y(start_y +j) &
                         + mat%A(start_a+j) * x(start_x +j)
                     y(start_x +j) = y(start_x +j) &
                         + conjg (mat%A(start_a+j)) * x(start_y +j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,ndiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
                  end_a = i*lda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     y(start_x +j) = y(start_x+j) &
                         + conjg (mat%A(start_a+j)) * x(start_y+j)
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*lda
                  end_a = i*lda
                  j = 1
                  do while((start_a + j).le.end_a)
                     y(start_y+j) = y(start_y+j) &
                         + mat%A(start_a+j) * x(start_x+j)
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,ndiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.mat%K-lda) then
               start_a = (i-1)*lda
               end_a = i*lda -mat%IA1(i)+mat%K-lda
            else if (mat%IA1(i).lt.-mat%M+lda) then
               start_a = (i-1)*lda -mat%IA1(i)-mat%M+lda
               end_a = i*lda               
            else
               start_a = (i-1)*lda
               end_a = i*lda
            end if
            j = 1
            do while((start_a + j).le.end_a)
               y(start_y+j) = y(start_y+j) &
                           + mat%A(start_a+j) * x(start_x+j)
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine zrmbv_dia 
! **********************************************************************
! **********************************************************************
      end module mod_rmbv_dia
