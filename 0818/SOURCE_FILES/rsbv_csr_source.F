      subroutine XRSBV_CSR(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,ofs,pntr
      character :: diag,part
      DCOMPLEX :: de

      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      ierr = -1

      if (part.eq.'L') then
         do i = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = ZZERO 
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq.ZZERO) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      else 
         do i = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(i)
               de = ZZERO
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     x(i) = x(i) &
                  - mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
                  else
                     de = mat%A(pntr + ofs)
                  end if
                  pntr = pntr + 1
               end do
               if(de.eq.ZZERO) then
                  ierr = blas_error_singtria
                  return
               else
                  x(i) = x(i)/de
               end if
            end if
         end do
         ierr = 0
      end if

      end subroutine XRSBV_CSR
