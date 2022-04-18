      module mod_rsbv_vbr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'VBR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface rsbv_vbr
        module procedure irsbv_vbr
        module procedure srsbv_vbr
        module procedure drsbv_vbr
        module procedure crsbv_vbr
        module procedure zrsbv_vbr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irsbv_vbr (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,mb,nb,dd
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,part,store
      integer , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      allocate(y(size(x)),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0 
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
          (size(mat%bp1).ne.size(mat%bp2)).or.&
          (maxval(abs(mat%bp1-mat%bp2)).ne.0)) then
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
      start_a = -1
      end_a = -1
      start_x = -1
      end_x = -1
      start_y = -1
      end_y = -1
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_left_lower(mat%A(start_a:end_a),&
                               x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = mb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_right_upper(mat%A(start_a:end_a),&
                                  x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine irsbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine srsbv_vbr (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,mb,nb,dd
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,part,store
      real(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      allocate(y(size(x)),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0e0 
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
          (size(mat%bp1).ne.size(mat%bp2)).or.&
          (maxval(abs(mat%bp1-mat%bp2)).ne.0)) then
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
      start_a = -1
      end_a = -1
      start_x = -1
      end_x = -1
      start_y = -1
      end_y = -1
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_left_lower(mat%A(start_a:end_a),&
                               x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = mb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_right_upper(mat%A(start_a:end_a),&
                                  x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine srsbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine drsbv_vbr (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,mb,nb,dd
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,part,store
      real(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      allocate(y(size(x)),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0d0 
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
          (size(mat%bp1).ne.size(mat%bp2)).or.&
          (maxval(abs(mat%bp1-mat%bp2)).ne.0)) then
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
      start_a = -1
      end_a = -1
      start_x = -1
      end_x = -1
      start_y = -1
      end_y = -1
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_left_lower(mat%A(start_a:end_a),&
                               x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = mb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_right_upper(mat%A(start_a:end_a),&
                                  x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine drsbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine crsbv_vbr (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,mb,nb,dd
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,part,store
      complex(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      allocate(y(size(x)),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0e0, 0.0e0) 
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
          (size(mat%bp1).ne.size(mat%bp2)).or.&
          (maxval(abs(mat%bp1-mat%bp2)).ne.0)) then
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
      start_a = -1
      end_a = -1
      start_x = -1
      end_x = -1
      start_y = -1
      end_y = -1
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_left_lower(mat%A(start_a:end_a),&
                               x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = mb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_right_upper(mat%A(start_a:end_a),&
                                  x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine crsbv_vbr 
! **********************************************************************
! **********************************************************************
      subroutine zrsbv_vbr (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,mb,nb,dd
      integer :: start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y
      character :: diag,part,store
      complex(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      allocate(y(size(x)),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0d0, 0.0d0) 
      if ((mat%FIDA.ne.'VBR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
          (size(mat%bp1).ne.size(mat%bp2)).or.&
          (maxval(abs(mat%bp1-mat%bp2)).ne.0)) then
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
      start_a = -1
      end_a = -1
      start_x = -1
      end_x = -1
      start_y = -1
      end_y = -1
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
          x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_left_lower(mat%A(start_a:end_a),&
                               x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = mb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                  end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                  len_x = end_x - start_x + 1
                  start_y = mat%bp1(i) + ofs
                  end_y = mat%bp1(i+1) + ofs -1
                  len_y = end_y - start_y + 1
                  start_a = mat%IA2(pntr+ofs) + ofs
                  end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                  call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     start_x = mat%bp2(mat%IA1(pntr+ofs)+ofs) + ofs
                     end_x = mat%bp2(mat%IA1(pntr+ofs)+ofs+1)+ofs-1
                     len_x = end_x - start_x + 1
                     start_y = mat%bp1(i) + ofs
                     end_y = mat%bp1(i+1) + ofs -1
                     len_y = end_y - start_y + 1
                     start_a = mat%IA2(pntr+ofs) + ofs
                     end_a = mat%IA2(pntr+ofs+1) + ofs - 1
                     call block_mult_vec(mat%A(start_a:end_a),&
           x(start_x:end_x),len_x,y(start_y:end_y),len_y,store,ierr) 
             x(start_y:end_y) = x(start_y:end_y) - y(start_y:end_y)
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  start_x = mat%bp1(i) + ofs
                  end_x = mat%bp1(i+1) + ofs -1
                  len_x = end_x - start_x + 1
                  start_a = mat%IA2(dd+ofs) + ofs
                  end_a = mat%IA2(dd+ofs+1) + ofs - 1
                  call invert_right_upper(mat%A(start_a:end_a),&
                                  x(start_x:end_x),len_x,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine zrsbv_vbr 
! **********************************************************************
! **********************************************************************
      end module mod_rsbv_vbr
