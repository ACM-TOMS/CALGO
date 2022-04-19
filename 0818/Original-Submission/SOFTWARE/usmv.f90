      module mod_usmv
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : MV MULTIPLICATION, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************
      use representation_of_data
      use properties
      use mod_mbv
      implicit none
      interface usmv
        module procedure iusmv
        module procedure susmv
        module procedure dusmv
        module procedure cusmv
        module procedure zusmv
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iusmv (a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      integer , dimension(:), intent(in) :: x      
      integer , dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      integer , intent(in), optional :: alpha
      integer , dimension(:), allocatable :: z
      type(ispmat ), pointer :: dspmtx
      integer :: transa_work
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
         allocate(z(size(y)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rmbv_coo(dspmtx,x,z,ierr)
            case('CSC')
               call rmbv_csc(dspmtx,x,z,ierr)
            case('CSR')
               call rmbv_csr(dspmtx,x,z,ierr)
            case('DIA')
               call rmbv_dia(dspmtx,x,z,ierr)
            case('BCO')
               call rmbv_bco(dspmtx,x,z,ierr)
            case('BSC')
               call rmbv_bsc(dspmtx,x,z,ierr)
            case('BSR')
               call rmbv_bsr(dspmtx,x,z,ierr)
            case('BDI')
               call rmbv_bdi(dspmtx,x,z,ierr)
            case('VBR')
               call rmbv_vbr(dspmtx,x,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lmbv_coo(dspmtx, (x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx, (x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx, (x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx, (x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx, (x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx, (x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx, (x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx, (x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx, (x),z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
               ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            deallocate(z,STAT=ierr)
            return
         end if            
         if(transa_work.eq.ORIGIN_MATRIX) then
            y = alpha_work * z + y
         else
            y = alpha_work * ( (z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine iusmv 
! **********************************************************************
! **********************************************************************
      subroutine susmv (a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      real(KIND=sp) , dimension(:), intent(in) :: x      
      real(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=sp) , intent(in), optional :: alpha
      real(KIND=sp) , dimension(:), allocatable :: z
      type(sspmat ), pointer :: dspmtx
      integer :: transa_work
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
         allocate(z(size(y)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rmbv_coo(dspmtx,x,z,ierr)
            case('CSC')
               call rmbv_csc(dspmtx,x,z,ierr)
            case('CSR')
               call rmbv_csr(dspmtx,x,z,ierr)
            case('DIA')
               call rmbv_dia(dspmtx,x,z,ierr)
            case('BCO')
               call rmbv_bco(dspmtx,x,z,ierr)
            case('BSC')
               call rmbv_bsc(dspmtx,x,z,ierr)
            case('BSR')
               call rmbv_bsr(dspmtx,x,z,ierr)
            case('BDI')
               call rmbv_bdi(dspmtx,x,z,ierr)
            case('VBR')
               call rmbv_vbr(dspmtx,x,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lmbv_coo(dspmtx, (x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx, (x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx, (x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx, (x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx, (x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx, (x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx, (x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx, (x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx, (x),z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
               ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            deallocate(z,STAT=ierr)
            return
         end if            
         if(transa_work.eq.ORIGIN_MATRIX) then
            y = alpha_work * z + y
         else
            y = alpha_work * ( (z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine susmv 
! **********************************************************************
! **********************************************************************
      subroutine dusmv (a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      real(KIND=dp) , dimension(:), intent(in) :: x      
      real(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=dp) , intent(in), optional :: alpha
      real(KIND=dp) , dimension(:), allocatable :: z
      type(dspmat ), pointer :: dspmtx
      integer :: transa_work
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
         allocate(z(size(y)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rmbv_coo(dspmtx,x,z,ierr)
            case('CSC')
               call rmbv_csc(dspmtx,x,z,ierr)
            case('CSR')
               call rmbv_csr(dspmtx,x,z,ierr)
            case('DIA')
               call rmbv_dia(dspmtx,x,z,ierr)
            case('BCO')
               call rmbv_bco(dspmtx,x,z,ierr)
            case('BSC')
               call rmbv_bsc(dspmtx,x,z,ierr)
            case('BSR')
               call rmbv_bsr(dspmtx,x,z,ierr)
            case('BDI')
               call rmbv_bdi(dspmtx,x,z,ierr)
            case('VBR')
               call rmbv_vbr(dspmtx,x,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lmbv_coo(dspmtx, (x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx, (x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx, (x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx, (x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx, (x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx, (x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx, (x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx, (x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx, (x),z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
               ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            deallocate(z,STAT=ierr)
            return
         end if            
         if(transa_work.eq.ORIGIN_MATRIX) then
            y = alpha_work * z + y
         else
            y = alpha_work * ( (z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine dusmv 
! **********************************************************************
! **********************************************************************
      subroutine cusmv (a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      complex(KIND=sp) , dimension(:), intent(in) :: x      
      complex(KIND=sp) , dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=sp) , intent(in), optional :: alpha
      complex(KIND=sp) , dimension(:), allocatable :: z
      type(cspmat ), pointer :: dspmtx
      integer :: transa_work
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
         allocate(z(size(y)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rmbv_coo(dspmtx,x,z,ierr)
            case('CSC')
               call rmbv_csc(dspmtx,x,z,ierr)
            case('CSR')
               call rmbv_csr(dspmtx,x,z,ierr)
            case('DIA')
               call rmbv_dia(dspmtx,x,z,ierr)
            case('BCO')
               call rmbv_bco(dspmtx,x,z,ierr)
            case('BSC')
               call rmbv_bsc(dspmtx,x,z,ierr)
            case('BSR')
               call rmbv_bsr(dspmtx,x,z,ierr)
            case('BDI')
               call rmbv_bdi(dspmtx,x,z,ierr)
            case('VBR')
               call rmbv_vbr(dspmtx,x,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lmbv_coo(dspmtx,conjg (x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx,conjg (x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx,conjg (x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx,conjg (x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx,conjg (x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx,conjg (x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx,conjg (x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx,conjg (x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx,conjg (x),z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
               ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            deallocate(z,STAT=ierr)
            return
         end if            
         if(transa_work.eq.ORIGIN_MATRIX) then
            y = alpha_work * z + y
         else
            y = alpha_work * (conjg (z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine cusmv 
! **********************************************************************
! **********************************************************************
      subroutine zusmv (a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      complex(KIND=dp) , dimension(:), intent(in) :: x      
      complex(KIND=dp) , dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=dp) , intent(in), optional :: alpha
      complex(KIND=dp) , dimension(:), allocatable :: z
      type(zspmat ), pointer :: dspmtx
      integer :: transa_work
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
         allocate(z(size(y)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rmbv_coo(dspmtx,x,z,ierr)
            case('CSC')
               call rmbv_csc(dspmtx,x,z,ierr)
            case('CSR')
               call rmbv_csr(dspmtx,x,z,ierr)
            case('DIA')
               call rmbv_dia(dspmtx,x,z,ierr)
            case('BCO')
               call rmbv_bco(dspmtx,x,z,ierr)
            case('BSC')
               call rmbv_bsc(dspmtx,x,z,ierr)
            case('BSR')
               call rmbv_bsr(dspmtx,x,z,ierr)
            case('BDI')
               call rmbv_bdi(dspmtx,x,z,ierr)
            case('VBR')
               call rmbv_vbr(dspmtx,x,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lmbv_coo(dspmtx,conjg (x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx,conjg (x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx,conjg (x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx,conjg (x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx,conjg (x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx,conjg (x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx,conjg (x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx,conjg (x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx,conjg (x),z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
               ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            deallocate(z,STAT=ierr)
            return
         end if            
         if(transa_work.eq.ORIGIN_MATRIX) then
            y = alpha_work * z + y
         else
            y = alpha_work * (conjg (z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine zusmv 
! **********************************************************************
! **********************************************************************
      end module mod_usmv
