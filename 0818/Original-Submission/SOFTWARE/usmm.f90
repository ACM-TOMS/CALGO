      module mod_usmm
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : MM MULTIPLICATION, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************
      use representation_of_data
      use properties
      use mod_mbv
      implicit none
      interface usmm
        module procedure iusmm
        module procedure susmm
        module procedure dusmm
        module procedure cusmm
        module procedure zusmm
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iusmm (a,b,c,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      integer , dimension(:,:), intent(in) :: b      
      integer , dimension(:,:), intent(inout) :: c
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      integer , intent(in), optional :: alpha
      integer , dimension(:), allocatable :: z
      type(ispmat ), pointer :: dspmtx
      integer transa_work,i
      integer  :: alpha_work
      ierr = -1
      if (present(transa)) then
         transa_work = transa
      else
         transa_work = ORIGIN_MATRIX
      end if
      if (present(alpha)) then
         alpha_work = alpha
      else
         alpha_work = 1.
      end if
      if (alpha_work.eq. 0 ) then
         !no matrix multiplication necessary
      else
         call accessdata_isp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(c,1)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rmbv_coo(dspmtx,b(:,i),z,ierr)
               case('CSC')
                  call rmbv_csc(dspmtx,b(:,i),z,ierr)
               case('CSR')
                  call rmbv_csr(dspmtx,b(:,i),z,ierr)
               case('DIA')
                  call rmbv_dia(dspmtx,b(:,i),z,ierr)
               case('BCO')
                  call rmbv_bco(dspmtx,b(:,i),z,ierr)
               case('BSC')
                  call rmbv_bsc(dspmtx,b(:,i),z,ierr)
               case('BSR')
                  call rmbv_bsr(dspmtx,b(:,i),z,ierr)
               case('BDI')
                  call rmbv_bdi(dspmtx,b(:,i),z,ierr)
               case('VBR')
                  call rmbv_vbr(dspmtx,b(:,i),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lmbv_coo(dspmtx, (b(:,i)),z,ierr)
               case('CSC')
                  call lmbv_csc(dspmtx, (b(:,i)),z,ierr)
               case('CSR')
                  call lmbv_csr(dspmtx, (b(:,i)),z,ierr)
               case('DIA')
                  call lmbv_dia(dspmtx, (b(:,i)),z,ierr)
               case('BCO')
                  call lmbv_bco(dspmtx, (b(:,i)),z,ierr)
               case('BSC')
                  call lmbv_bsc(dspmtx, (b(:,i)),z,ierr)
               case('BSR')
                  call lmbv_bsr(dspmtx, (b(:,i)),z,ierr)
               case('BDI')
                  call lmbv_bdi(dspmtx, (b(:,i)),z,ierr)
               case('VBR')
                  call lmbv_vbr(dspmtx, (b(:,i)),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               deallocate(z,STAT=ierr)
               return
            else
               if(transa_work.eq.ORIGIN_MATRIX) then
                  c(:,i) = alpha_work * z + c(:,i)
               else
                  c(:,i) = alpha_work * ( (z)) + c(:,i)
               end if
            end if      
         end do
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine iusmm 
! **********************************************************************
! **********************************************************************
      subroutine susmm (a,b,c,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      real(KIND=sp) , dimension(:,:), intent(in) :: b      
      real(KIND=sp) , dimension(:,:), intent(inout) :: c
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=sp) , intent(in), optional :: alpha
      real(KIND=sp) , dimension(:), allocatable :: z
      type(sspmat ), pointer :: dspmtx
      integer transa_work,i
      real(KIND=sp)  :: alpha_work
      ierr = -1
      if (present(transa)) then
         transa_work = transa
      else
         transa_work = ORIGIN_MATRIX
      end if
      if (present(alpha)) then
         alpha_work = alpha
      else
         alpha_work = 1.
      end if
      if (alpha_work.eq. 0.0e0 ) then
         !no matrix multiplication necessary
      else
         call accessdata_ssp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(c,1)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rmbv_coo(dspmtx,b(:,i),z,ierr)
               case('CSC')
                  call rmbv_csc(dspmtx,b(:,i),z,ierr)
               case('CSR')
                  call rmbv_csr(dspmtx,b(:,i),z,ierr)
               case('DIA')
                  call rmbv_dia(dspmtx,b(:,i),z,ierr)
               case('BCO')
                  call rmbv_bco(dspmtx,b(:,i),z,ierr)
               case('BSC')
                  call rmbv_bsc(dspmtx,b(:,i),z,ierr)
               case('BSR')
                  call rmbv_bsr(dspmtx,b(:,i),z,ierr)
               case('BDI')
                  call rmbv_bdi(dspmtx,b(:,i),z,ierr)
               case('VBR')
                  call rmbv_vbr(dspmtx,b(:,i),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lmbv_coo(dspmtx, (b(:,i)),z,ierr)
               case('CSC')
                  call lmbv_csc(dspmtx, (b(:,i)),z,ierr)
               case('CSR')
                  call lmbv_csr(dspmtx, (b(:,i)),z,ierr)
               case('DIA')
                  call lmbv_dia(dspmtx, (b(:,i)),z,ierr)
               case('BCO')
                  call lmbv_bco(dspmtx, (b(:,i)),z,ierr)
               case('BSC')
                  call lmbv_bsc(dspmtx, (b(:,i)),z,ierr)
               case('BSR')
                  call lmbv_bsr(dspmtx, (b(:,i)),z,ierr)
               case('BDI')
                  call lmbv_bdi(dspmtx, (b(:,i)),z,ierr)
               case('VBR')
                  call lmbv_vbr(dspmtx, (b(:,i)),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               deallocate(z,STAT=ierr)
               return
            else
               if(transa_work.eq.ORIGIN_MATRIX) then
                  c(:,i) = alpha_work * z + c(:,i)
               else
                  c(:,i) = alpha_work * ( (z)) + c(:,i)
               end if
            end if      
         end do
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine susmm 
! **********************************************************************
! **********************************************************************
      subroutine dusmm (a,b,c,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      real(KIND=dp) , dimension(:,:), intent(in) :: b      
      real(KIND=dp) , dimension(:,:), intent(inout) :: c
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=dp) , intent(in), optional :: alpha
      real(KIND=dp) , dimension(:), allocatable :: z
      type(dspmat ), pointer :: dspmtx
      integer transa_work,i
      real(KIND=dp)  :: alpha_work
      ierr = -1
      if (present(transa)) then
         transa_work = transa
      else
         transa_work = ORIGIN_MATRIX
      end if
      if (present(alpha)) then
         alpha_work = alpha
      else
         alpha_work = 1.
      end if
      if (alpha_work.eq. 0.0d0 ) then
         !no matrix multiplication necessary
      else
         call accessdata_dsp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(c,1)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rmbv_coo(dspmtx,b(:,i),z,ierr)
               case('CSC')
                  call rmbv_csc(dspmtx,b(:,i),z,ierr)
               case('CSR')
                  call rmbv_csr(dspmtx,b(:,i),z,ierr)
               case('DIA')
                  call rmbv_dia(dspmtx,b(:,i),z,ierr)
               case('BCO')
                  call rmbv_bco(dspmtx,b(:,i),z,ierr)
               case('BSC')
                  call rmbv_bsc(dspmtx,b(:,i),z,ierr)
               case('BSR')
                  call rmbv_bsr(dspmtx,b(:,i),z,ierr)
               case('BDI')
                  call rmbv_bdi(dspmtx,b(:,i),z,ierr)
               case('VBR')
                  call rmbv_vbr(dspmtx,b(:,i),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lmbv_coo(dspmtx, (b(:,i)),z,ierr)
               case('CSC')
                  call lmbv_csc(dspmtx, (b(:,i)),z,ierr)
               case('CSR')
                  call lmbv_csr(dspmtx, (b(:,i)),z,ierr)
               case('DIA')
                  call lmbv_dia(dspmtx, (b(:,i)),z,ierr)
               case('BCO')
                  call lmbv_bco(dspmtx, (b(:,i)),z,ierr)
               case('BSC')
                  call lmbv_bsc(dspmtx, (b(:,i)),z,ierr)
               case('BSR')
                  call lmbv_bsr(dspmtx, (b(:,i)),z,ierr)
               case('BDI')
                  call lmbv_bdi(dspmtx, (b(:,i)),z,ierr)
               case('VBR')
                  call lmbv_vbr(dspmtx, (b(:,i)),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               deallocate(z,STAT=ierr)
               return
            else
               if(transa_work.eq.ORIGIN_MATRIX) then
                  c(:,i) = alpha_work * z + c(:,i)
               else
                  c(:,i) = alpha_work * ( (z)) + c(:,i)
               end if
            end if      
         end do
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine dusmm 
! **********************************************************************
! **********************************************************************
      subroutine cusmm (a,b,c,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      complex(KIND=sp) , dimension(:,:), intent(in) :: b      
      complex(KIND=sp) , dimension(:,:), intent(inout) :: c
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=sp) , intent(in), optional :: alpha
      complex(KIND=sp) , dimension(:), allocatable :: z
      type(cspmat ), pointer :: dspmtx
      integer transa_work,i
      complex(KIND=sp)  :: alpha_work
      ierr = -1
      if (present(transa)) then
         transa_work = transa
      else
         transa_work = ORIGIN_MATRIX
      end if
      if (present(alpha)) then
         alpha_work = alpha
      else
         alpha_work = 1.
      end if
      if (alpha_work.eq. (0.0e0, 0.0e0) ) then
         !no matrix multiplication necessary
      else
         call accessdata_csp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(c,1)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rmbv_coo(dspmtx,b(:,i),z,ierr)
               case('CSC')
                  call rmbv_csc(dspmtx,b(:,i),z,ierr)
               case('CSR')
                  call rmbv_csr(dspmtx,b(:,i),z,ierr)
               case('DIA')
                  call rmbv_dia(dspmtx,b(:,i),z,ierr)
               case('BCO')
                  call rmbv_bco(dspmtx,b(:,i),z,ierr)
               case('BSC')
                  call rmbv_bsc(dspmtx,b(:,i),z,ierr)
               case('BSR')
                  call rmbv_bsr(dspmtx,b(:,i),z,ierr)
               case('BDI')
                  call rmbv_bdi(dspmtx,b(:,i),z,ierr)
               case('VBR')
                  call rmbv_vbr(dspmtx,b(:,i),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lmbv_coo(dspmtx,conjg (b(:,i)),z,ierr)
               case('CSC')
                  call lmbv_csc(dspmtx,conjg (b(:,i)),z,ierr)
               case('CSR')
                  call lmbv_csr(dspmtx,conjg (b(:,i)),z,ierr)
               case('DIA')
                  call lmbv_dia(dspmtx,conjg (b(:,i)),z,ierr)
               case('BCO')
                  call lmbv_bco(dspmtx,conjg (b(:,i)),z,ierr)
               case('BSC')
                  call lmbv_bsc(dspmtx,conjg (b(:,i)),z,ierr)
               case('BSR')
                  call lmbv_bsr(dspmtx,conjg (b(:,i)),z,ierr)
               case('BDI')
                  call lmbv_bdi(dspmtx,conjg (b(:,i)),z,ierr)
               case('VBR')
                  call lmbv_vbr(dspmtx,conjg (b(:,i)),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               deallocate(z,STAT=ierr)
               return
            else
               if(transa_work.eq.ORIGIN_MATRIX) then
                  c(:,i) = alpha_work * z + c(:,i)
               else
                  c(:,i) = alpha_work * (conjg (z)) + c(:,i)
               end if
            end if      
         end do
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine cusmm 
! **********************************************************************
! **********************************************************************
      subroutine zusmm (a,b,c,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      complex(KIND=dp) , dimension(:,:), intent(in) :: b      
      complex(KIND=dp) , dimension(:,:), intent(inout) :: c
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=dp) , intent(in), optional :: alpha
      complex(KIND=dp) , dimension(:), allocatable :: z
      type(zspmat ), pointer :: dspmtx
      integer transa_work,i
      complex(KIND=dp)  :: alpha_work
      ierr = -1
      if (present(transa)) then
         transa_work = transa
      else
         transa_work = ORIGIN_MATRIX
      end if
      if (present(alpha)) then
         alpha_work = alpha
      else
         alpha_work = 1.
      end if
      if (alpha_work.eq. (0.0d0, 0.0d0) ) then
         !no matrix multiplication necessary
      else
         call accessdata_zsp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(c,1)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rmbv_coo(dspmtx,b(:,i),z,ierr)
               case('CSC')
                  call rmbv_csc(dspmtx,b(:,i),z,ierr)
               case('CSR')
                  call rmbv_csr(dspmtx,b(:,i),z,ierr)
               case('DIA')
                  call rmbv_dia(dspmtx,b(:,i),z,ierr)
               case('BCO')
                  call rmbv_bco(dspmtx,b(:,i),z,ierr)
               case('BSC')
                  call rmbv_bsc(dspmtx,b(:,i),z,ierr)
               case('BSR')
                  call rmbv_bsr(dspmtx,b(:,i),z,ierr)
               case('BDI')
                  call rmbv_bdi(dspmtx,b(:,i),z,ierr)
               case('VBR')
                  call rmbv_vbr(dspmtx,b(:,i),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lmbv_coo(dspmtx,conjg (b(:,i)),z,ierr)
               case('CSC')
                  call lmbv_csc(dspmtx,conjg (b(:,i)),z,ierr)
               case('CSR')
                  call lmbv_csr(dspmtx,conjg (b(:,i)),z,ierr)
               case('DIA')
                  call lmbv_dia(dspmtx,conjg (b(:,i)),z,ierr)
               case('BCO')
                  call lmbv_bco(dspmtx,conjg (b(:,i)),z,ierr)
               case('BSC')
                  call lmbv_bsc(dspmtx,conjg (b(:,i)),z,ierr)
               case('BSR')
                  call lmbv_bsr(dspmtx,conjg (b(:,i)),z,ierr)
               case('BDI')
                  call lmbv_bdi(dspmtx,conjg (b(:,i)),z,ierr)
               case('VBR')
                  call lmbv_vbr(dspmtx,conjg (b(:,i)),z,ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               deallocate(z,STAT=ierr)
               return
            else
               if(transa_work.eq.ORIGIN_MATRIX) then
                  c(:,i) = alpha_work * z + c(:,i)
               else
                  c(:,i) = alpha_work * (conjg (z)) + c(:,i)
               end if
            end if      
         end do
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine zusmm 
! **********************************************************************
! **********************************************************************
      end module mod_usmm
