      subroutine XLMBV_COO(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
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

      y = ZZERO

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
                           + CONJUG(mat%A(i)) * x(mat%IA2(i)+ofs)
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
                           + CONJUG(mat%A(i)) * x(mat%IA2(i)+ofs)
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
      end subroutine XLMBV_COO
