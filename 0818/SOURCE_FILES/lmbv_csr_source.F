      subroutine XLMBV_CSR(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,i,pntr
      character :: diag,type,part
      
      ierr = -1
      m = size(y)
      n = size(x)

      if ((mat%FIDA.ne.'CSR').or.(mat%M.ne.n).or.(mat%K.ne.m)) then
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
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then 
                     y(i) = y(i) &
                  + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
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
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then 
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
      else if ((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) &
         then 
         if (part.eq.'L') then
            do i = 1, mat%M
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  if(i.eq.mat%IA1(pntr + ofs) + ofs) then
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  else if (i.gt.mat%IA1(pntr + ofs) + ofs) then 
                     y(i) = y(i) + CONJUG(mat%A(pntr + ofs)) &
                                   * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
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
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
                  else if (i.lt.mat%IA1(pntr + ofs) + ofs) then 
                     y(i) = y(i) + CONJUG(mat%A(pntr + ofs))&
                                   * x(mat%IA1(pntr + ofs ) + ofs) 
                     y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
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
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(i) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine XLMBV_CSR
