      subroutine XRMBV_CSC(mat,x,y,ierr)
      implicit none
      type(XSPMAT), pointer :: mat
      DCOMPLEX, dimension(:), intent(in) :: x
      DCOMPLEX, dimension(:), intent(out) :: y
      integer, intent(out) :: ierr
      integer :: m,n,base,ofs,j,pntr
      character :: diag,type,part
      
      ierr = -1
      m = size(y)
      n = size(x)

      if ((mat%FIDA.ne.'CSC').or.(mat%M.ne.m).or.(mat%K.ne.n)) then
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
         if (part.eq.'U') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) &
                   + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) &
                   + mat%A(pntr + ofs) * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else if((type.eq.'H').and.(.not.(part.eq.'B')).and.(m.eq.n)) then 
         if (part.eq.'U') then
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               else if (j.gt.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               y(j) = y(j) &
              +CONJUG(mat%A(pntr + ofs)) * x(mat%IA1(pntr + ofs)+ ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
               if(j.eq.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               else if (j.lt.mat%IA1(pntr + ofs) + ofs) then
               y(mat%IA1(pntr + ofs)+ofs) =y(mat%IA1(pntr+ofs)+ofs) &
                                          +mat%A(pntr + ofs) * x(j) 
                    y(j) = y(j) + CONJUG(mat%A(pntr+ofs)) &
                                   * x(mat%IA1(pntr + ofs) + ofs) 
               end if
               pntr = pntr + 1
            end do
         end do
         end if
         ierr = 0
      else
         do j = 1, mat%K
            pntr = mat%pb(j)
            do while(pntr.lt.mat%pe(j))
              y(mat%IA1(pntr + ofs) + ofs) = &
              y(mat%IA1(pntr + ofs ) + ofs) + mat%A(pntr + ofs) * x(j) 
               pntr = pntr + 1
            end do
         end do
         ierr = 0
      end if
      end subroutine XRMBV_CSC
