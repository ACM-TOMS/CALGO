      subroutine XLSBV_DIA(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,j,n,lda,ndiag
      character :: diag,part
      DCOMPLEX :: de

      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'DIA').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
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
      ierr = -1

      if (part.eq.'L') then
         do i=n,1,-1
            if (diag.eq.'U') then
               do j = 1,ndiag
                 if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                 end if
               end do
            else
               de = ZZERO
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
               else if((-mat%IA1(j).lt.n-i+1).and.(mat%IA1(j).lt.0)) &
                  then
             x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j))*x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne.ZZERO) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      else
         do i=1,n
            if (diag.eq.'U') then
               do j = 1,ndiag
                  if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
            else
               de = ZZERO
               do j = 1,ndiag
                  if (mat%IA1(j).eq.0) then
                     de = mat%A(lda*(j-1) + i)
                  else if((mat%IA1(j).lt.i).and.(mat%IA1(j).gt.0)) &
                  then
            x(i) = x(i)-mat%A(lda*(j-1)+i-mat%IA1(j)) *x(i-mat%IA1(j))
                  end if
               end do
               if (de.ne.ZZERO) then 
                  x(i) = x(i)/de
               else
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         ierr = 0
      end if

      end subroutine XLSBV_DIA
