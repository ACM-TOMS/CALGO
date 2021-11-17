      module mod_rmbv_bdi
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'BDI'-STORAGE
!                   rmbv = Right Multiplication By Vector: y=Ax
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface rmbv_bdi
        module procedure irmbv_bdi
        module procedure srmbv_bdi
        module procedure drmbv_bdi
        module procedure crmbv_bdi
        module procedure zrmbv_bdi
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irmbv_bdi (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j,mm,nn,nn_sq
      integer :: blda,nbdiag,start_a,end_a,start_x,start_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
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
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
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
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,nbdiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.(mat%K/nn)-blda) then
               start_a = (i-1)*blda
               end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
            else if (mat%IA1(i).lt.-(mat%M/nn)+blda) then
               start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
               end_a = i*blda               
            else
               start_a = (i-1)*blda
               end_a = i*blda
            end if
            j = 1
            do while((start_a+j).le.end_a)
               call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine irmbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine srmbv_bdi (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j,mm,nn,nn_sq
      integer :: blda,nbdiag,start_a,end_a,start_x,start_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
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
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
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
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,nbdiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.(mat%K/nn)-blda) then
               start_a = (i-1)*blda
               end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
            else if (mat%IA1(i).lt.-(mat%M/nn)+blda) then
               start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
               end_a = i*blda               
            else
               start_a = (i-1)*blda
               end_a = i*blda
            end if
            j = 1
            do while((start_a+j).le.end_a)
               call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine srmbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine drmbv_bdi (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j,mm,nn,nn_sq
      integer :: blda,nbdiag,start_a,end_a,start_x,start_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
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
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
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
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,nbdiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.(mat%K/nn)-blda) then
               start_a = (i-1)*blda
               end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
            else if (mat%IA1(i).lt.-(mat%M/nn)+blda) then
               start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
               end_a = i*blda               
            else
               start_a = (i-1)*blda
               end_a = i*blda
            end if
            j = 1
            do while((start_a+j).le.end_a)
               call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine drmbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine crmbv_bdi (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j,mm,nn,nn_sq
      integer :: blda,nbdiag,start_a,end_a,start_x,start_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
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
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
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
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,nbdiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.(mat%K/nn)-blda) then
               start_a = (i-1)*blda
               end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
            else if (mat%IA1(i).lt.-(mat%M/nn)+blda) then
               start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
               end_a = i*blda               
            else
               start_a = (i-1)*blda
               end_a = i*blda
            end if
            j = 1
            do while((start_a+j).le.end_a)
               call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine crmbv_bdi 
! **********************************************************************
! **********************************************************************
      subroutine zrmbv_bdi (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,i,j,mm,nn,nn_sq
      integer :: blda,nbdiag,start_a,end_a,start_x,start_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
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
      call get_infoa(mat%INFOA,'f',blda,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nbdiag,ierr)
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
      if ((mat%FIDA.ne.'BDI').or.(mat%M.ne.m).or.(mat%K.ne.n).or.&
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
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_T_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).gt.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         else !(part.eq.'L')
            do i=1,nbdiag
               start_x = max(0,mat%IA1(i))
               start_y = max(0,-mat%IA1(i))
               if (mat%IA1(i).lt.0) then
                  start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
                  end_a = i*blda               
                  j = 1
                  do while((start_a + j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     call block_H_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_y+j-1)*nn+1:(start_y+j)*nn),&
      nn,y((start_x+j-1)*nn+1:(start_x+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else if (mat%IA1(i).eq.0) then
                  start_a = (i-1)*blda
                  end_a = i*blda
                  j = 1
                  do while((start_a+j).le.end_a)
                     call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
                     j = j+1
                  end do
               else
                  cycle
               end if
            end do
         end if
         ierr = 0
      else !no symmetry
         do i=1,nbdiag
            start_x = max(0,mat%IA1(i))
            start_y = max(0,-mat%IA1(i))
            if (mat%IA1(i).gt.(mat%K/nn)-blda) then
               start_a = (i-1)*blda
               end_a = i*blda -mat%IA1(i)+(mat%K/nn)-blda
            else if (mat%IA1(i).lt.-(mat%M/nn)+blda) then
               start_a = (i-1)*blda -mat%IA1(i)-(mat%M/nn)+blda
               end_a = i*blda               
            else
               start_a = (i-1)*blda
               end_a = i*blda
            end if
            j = 1
            do while((start_a+j).le.end_a)
               call block_mult_vec(&
      mat%A((start_a+j-1)*nn_sq+1:(start_a+j)*nn_sq),&
      x((start_x+j-1)*nn+1:(start_x+j)*nn),&
      nn,y((start_y+j-1)*nn+1:(start_y+j)*nn),nn,store,ierr) 
               j = j+1
            end do
         end do
         ierr = 0
      end if
      end subroutine zrmbv_bdi 
! **********************************************************************
! **********************************************************************
      end module mod_rmbv_bdi
