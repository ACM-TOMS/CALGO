      subroutine XRSBV_CSC(mat,x,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,ofs,pntr
      character :: diag,part
      DCOMPLEX :: de

      ierr = -1
      n = size(x)
      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.n).or.(mat%K.ne.n)) then
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
         do j = 1,n
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq.ZZERO) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0

      else
         do j = n,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
               end do
            else
               pntr = mat%pb(j)
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr + ofs) + ofs.ne.j))
                  pntr = pntr + 1
               end do
               if(mat%IA1(pntr + ofs) + ofs.eq.j) then
                  de = mat%A(pntr + ofs)
               else
                  ierr = blas_error_singtria
                  return
               end if
               if(de.eq.ZZERO) then
                  ierr = blas_error_singtria
                  return
               else
                  x(j) = x(j)/de
                  de = x(j)
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  x(mat%IA1(pntr + ofs) + ofs) = &
              x(mat%IA1(pntr + ofs ) + ofs) - mat%A(pntr + ofs) * de
                  pntr = pntr + 1
               end do
               x(j) = de
            end if
         end do
         ierr = 0
      end if

      end subroutine XRSBV_CSC
