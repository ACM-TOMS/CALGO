      module mod_rsbv_bsr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSR'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface rsbv_bsr
        module procedure irsbv_bsr
        module procedure srsbv_bsr
        module procedure drsbv_bsr
        module procedure crsbv_bsr
        module procedure zrsbv_bsr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irsbv_bsr (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      integer , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      if ((mat%FIDA.ne.'BSR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0 
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
      end subroutine irsbv_bsr 
! **********************************************************************
! **********************************************************************
      subroutine srsbv_bsr (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      real(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      if ((mat%FIDA.ne.'BSR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0e0 
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
      end subroutine srsbv_bsr 
! **********************************************************************
! **********************************************************************
      subroutine drsbv_bsr (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      real(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      if ((mat%FIDA.ne.'BSR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0d0 
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
      end subroutine drsbv_bsr 
! **********************************************************************
! **********************************************************************
      subroutine crsbv_bsr (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      complex(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      if ((mat%FIDA.ne.'BSR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0e0, 0.0e0) 
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
      end subroutine crsbv_bsr 
! **********************************************************************
! **********************************************************************
      subroutine zrsbv_bsr (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      complex(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      if ((mat%FIDA.ne.'BSR').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0d0, 0.0d0) 
      if (part.eq.'L') then
         do i = 1,mb
            if (diag.eq.'U') then
               pntr = mat%pb(i)
               do while(pntr.lt.mat%pe(i))
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
                  call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                  x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(i)
               dd = -1
               do while(pntr.lt.mat%pe(i))
                  if(mat%IA1(pntr + ofs) + ofs.ne.i) then
                     call block_mult_vec(&
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn),&
      nn,y,nn,store,ierr) 
                     x((i-1)*nn+1:i*nn)= x((i-1)*nn+1:i*nn) - y
                  else
                     dd = pntr
                  end if 
                  pntr = pntr + 1
               end do 
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper(&
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
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
      end subroutine zrsbv_bsr 
! **********************************************************************
! **********************************************************************
      end module mod_rsbv_bsr
