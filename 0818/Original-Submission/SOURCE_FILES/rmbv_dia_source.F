      subroutine XRMBV_DIA(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
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
                         + CONJUG(mat%A(start_a+j)) * x(start_y +j)
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
                         + CONJUG(mat%A(start_a+j)) * x(start_y+j)
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
      end subroutine XRMBV_DIA
