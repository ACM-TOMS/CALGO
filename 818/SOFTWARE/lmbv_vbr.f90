      module mod_lmbv_vbr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS MV MULT. WITH MATRIX IN 'VBR'-STORAGE
!                   lmbv = Left Multiplication By Vector: y^T=x^TA
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface lmbv_vbr
        module procedure ilmbv_vbr
        module procedure slmbv_vbr
        module procedure dlmbv_vbr
        module procedure clmbv_vbr
        module procedure zlmbv_vbr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilmbv_vbr (mat,x,y,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(in) :: x
      integer , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr,mb,nb
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
         y = 0 
         start_a = -1
         end_a = -1
         start_x = -1
         end_x = -1
         start_y = -1
         end_y = -1
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
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, mb
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
               start_x = mat%bp1(i) + ofs
               end_x = mat%bp1(i+1) + ofs -1
               len_x = end_x - start_x + 1
               start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
               end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
               len_y = end_y - start_y + 1
               start_a = mat%IA2(pntr+ofs) + ofs
               end_a = mat%IA2(pntr+ofs+1) + ofs - 1
               call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine ilmbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine slmbv_vbr (mat,x,y,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(in) :: x
      real(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr,mb,nb
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
         y = 0.0e0 
         start_a = -1
         end_a = -1
         start_x = -1
         end_x = -1
         start_y = -1
         end_y = -1
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
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, mb
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
               start_x = mat%bp1(i) + ofs
               end_x = mat%bp1(i+1) + ofs -1
               len_x = end_x - start_x + 1
               start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
               end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
               len_y = end_y - start_y + 1
               start_a = mat%IA2(pntr+ofs) + ofs
               end_a = mat%IA2(pntr+ofs+1) + ofs - 1
               call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine slmbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine dlmbv_vbr (mat,x,y,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(in) :: x
      real(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr,mb,nb
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
         y = 0.0d0 
         start_a = -1
         end_a = -1
         start_x = -1
         end_x = -1
         start_y = -1
         end_y = -1
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
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, mb
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
               start_x = mat%bp1(i) + ofs
               end_x = mat%bp1(i+1) + ofs -1
               len_x = end_x - start_x + 1
               start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
               end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
               len_y = end_y - start_y + 1
               start_a = mat%IA2(pntr+ofs) + ofs
               end_a = mat%IA2(pntr+ofs+1) + ofs - 1
               call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine dlmbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine clmbv_vbr (mat,x,y,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(in) :: x
      complex(KIND=sp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr,mb,nb
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
         y = (0.0e0, 0.0e0) 
         start_a = -1
         end_a = -1
         start_x = -1
         end_x = -1
         start_y = -1
         end_y = -1
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
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, mb
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
               start_x = mat%bp1(i) + ofs
               end_x = mat%bp1(i+1) + ofs -1
               len_x = end_x - start_x + 1
               start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
               end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
               len_y = end_y - start_y + 1
               start_a = mat%IA2(pntr+ofs) + ofs
               end_a = mat%IA2(pntr+ofs+1) + ofs - 1
               call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine clmbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine zlmbv_vbr (mat,x,y,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(in) :: x
      complex(KIND=dp) , dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr,mb,nb
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,type,part,store
      ierr = -1
      m = size(y)
      n = size(x)
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
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
         y = (0.0d0, 0.0d0) 
         start_a = -1
         end_a = -1
         start_x = -1
         end_x = -1
         start_y = -1
         end_y = -1
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
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'L') then
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
            x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
            x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         else
            do i = 1, mb
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then
                     call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
                     call block_Z_mult_vec(mat%A(start_a:end_a),&
             x(start_y:end_y),len_y,y(start_x:end_x),len_x,store,ierr) 
                  end if
                  pntr = pntr + 1
               end do
            end do
         end if
         ierr = 0
      else !no symmetry
         do i = 1, mb
            pntr = mat%pb(i)
            do while(pntr.lt.mat%pe(i))
               start_x = mat%bp1(i) + ofs
               end_x = mat%bp1(i+1) + ofs -1
               len_x = end_x - start_x + 1
               start_y = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
               end_y = mat%bp2(mat%IA1(pntr+ofs)+ofs+1) + ofs - 1
               len_y = end_y - start_y + 1
               start_a = mat%IA2(pntr+ofs) + ofs
               end_a = mat%IA2(pntr+ofs+1) + ofs - 1
               call block_T_mult_vec(mat%A(start_a:end_a),&
             x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine zlmbv_vbr 
! **********************************************************************
! **********************************************************************
      end module mod_lmbv_vbr
