      module mod_ussv
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : TRI. SOLVE, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************
      use representation_of_data
      use properties
      use mod_sbv
      implicit none
      interface ussv
        module procedure iussv
        module procedure sussv
        module procedure dussv
        module procedure cussv
        module procedure zussv
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iussv (a,x,ierr,transa,alpha)
      integer, intent(in) :: a
      integer , intent(inout) :: x(:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      integer , intent(in), optional :: alpha
      integer :: transa_work
      integer  :: alpha_work
      integer , dimension(:), allocatable :: z
      type(ispmat ), pointer :: dspmtx
      character :: type
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
      if (alpha_work.ne. 0 ) then
         call accessdata_isp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         call get_descra(dspmtx%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         if (type.ne.'T') then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(x)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         z=  (x)
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rsbv_coo(dspmtx,x,ierr)
            case('CSC')
               call rsbv_csc(dspmtx,x,ierr)
            case('CSR')
               call rsbv_csr(dspmtx,x,ierr)
            case('DIA')
               call rsbv_dia(dspmtx,x,ierr)
            case('BCO')
               call rsbv_bco(dspmtx,x,ierr)
            case('BSC')
               call rsbv_bsc(dspmtx,x,ierr)
            case('BSR')
               call rsbv_bsr(dspmtx,x,ierr)
            case('BDI')
               call rsbv_bdi(dspmtx,x,ierr)
            case('VBR')
               call rsbv_vbr(dspmtx,x,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lsbv_coo(dspmtx,z,ierr)
            case('CSC')
               call lsbv_csc(dspmtx,z,ierr)
            case('CSR')
               call lsbv_csr(dspmtx,z,ierr)
            case('DIA')
               call lsbv_dia(dspmtx,z,ierr)
            case('BCO')
               call lsbv_bco(dspmtx,z,ierr)
            case('BSC')
               call lsbv_bsc(dspmtx,z,ierr)
            case('BSR')
               call lsbv_bsr(dspmtx,z,ierr)
            case('BDI')
               call lsbv_bdi(dspmtx,z,ierr)
            case('VBR')
               call lsbv_vbr(dspmtx,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
            ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            return
         end if            
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         x = alpha_work * x
      else
         x = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine iussv 
! **********************************************************************
! **********************************************************************
      subroutine sussv (a,x,ierr,transa,alpha)
      integer, intent(in) :: a
      real(KIND=sp) , intent(inout) :: x(:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=sp) , intent(in), optional :: alpha
      integer :: transa_work
      real(KIND=sp)  :: alpha_work
      real(KIND=sp) , dimension(:), allocatable :: z
      type(sspmat ), pointer :: dspmtx
      character :: type
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
      if (alpha_work.ne. 0.0e0 ) then
         call accessdata_ssp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         call get_descra(dspmtx%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         if (type.ne.'T') then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(x)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         z=  (x)
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rsbv_coo(dspmtx,x,ierr)
            case('CSC')
               call rsbv_csc(dspmtx,x,ierr)
            case('CSR')
               call rsbv_csr(dspmtx,x,ierr)
            case('DIA')
               call rsbv_dia(dspmtx,x,ierr)
            case('BCO')
               call rsbv_bco(dspmtx,x,ierr)
            case('BSC')
               call rsbv_bsc(dspmtx,x,ierr)
            case('BSR')
               call rsbv_bsr(dspmtx,x,ierr)
            case('BDI')
               call rsbv_bdi(dspmtx,x,ierr)
            case('VBR')
               call rsbv_vbr(dspmtx,x,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lsbv_coo(dspmtx,z,ierr)
            case('CSC')
               call lsbv_csc(dspmtx,z,ierr)
            case('CSR')
               call lsbv_csr(dspmtx,z,ierr)
            case('DIA')
               call lsbv_dia(dspmtx,z,ierr)
            case('BCO')
               call lsbv_bco(dspmtx,z,ierr)
            case('BSC')
               call lsbv_bsc(dspmtx,z,ierr)
            case('BSR')
               call lsbv_bsr(dspmtx,z,ierr)
            case('BDI')
               call lsbv_bdi(dspmtx,z,ierr)
            case('VBR')
               call lsbv_vbr(dspmtx,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
            ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            return
         end if            
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         x = alpha_work * x
      else
         x = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine sussv 
! **********************************************************************
! **********************************************************************
      subroutine dussv (a,x,ierr,transa,alpha)
      integer, intent(in) :: a
      real(KIND=dp) , intent(inout) :: x(:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=dp) , intent(in), optional :: alpha
      integer :: transa_work
      real(KIND=dp)  :: alpha_work
      real(KIND=dp) , dimension(:), allocatable :: z
      type(dspmat ), pointer :: dspmtx
      character :: type
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
      if (alpha_work.ne. 0.0d0 ) then
         call accessdata_dsp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         call get_descra(dspmtx%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         if (type.ne.'T') then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(x)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         z=  (x)
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rsbv_coo(dspmtx,x,ierr)
            case('CSC')
               call rsbv_csc(dspmtx,x,ierr)
            case('CSR')
               call rsbv_csr(dspmtx,x,ierr)
            case('DIA')
               call rsbv_dia(dspmtx,x,ierr)
            case('BCO')
               call rsbv_bco(dspmtx,x,ierr)
            case('BSC')
               call rsbv_bsc(dspmtx,x,ierr)
            case('BSR')
               call rsbv_bsr(dspmtx,x,ierr)
            case('BDI')
               call rsbv_bdi(dspmtx,x,ierr)
            case('VBR')
               call rsbv_vbr(dspmtx,x,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lsbv_coo(dspmtx,z,ierr)
            case('CSC')
               call lsbv_csc(dspmtx,z,ierr)
            case('CSR')
               call lsbv_csr(dspmtx,z,ierr)
            case('DIA')
               call lsbv_dia(dspmtx,z,ierr)
            case('BCO')
               call lsbv_bco(dspmtx,z,ierr)
            case('BSC')
               call lsbv_bsc(dspmtx,z,ierr)
            case('BSR')
               call lsbv_bsr(dspmtx,z,ierr)
            case('BDI')
               call lsbv_bdi(dspmtx,z,ierr)
            case('VBR')
               call lsbv_vbr(dspmtx,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
            ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            return
         end if            
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         x = alpha_work * x
      else
         x = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine dussv 
! **********************************************************************
! **********************************************************************
      subroutine cussv (a,x,ierr,transa,alpha)
      integer, intent(in) :: a
      complex(KIND=sp) , intent(inout) :: x(:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=sp) , intent(in), optional :: alpha
      integer :: transa_work
      complex(KIND=sp)  :: alpha_work
      complex(KIND=sp) , dimension(:), allocatable :: z
      type(cspmat ), pointer :: dspmtx
      character :: type
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
      if (alpha_work.ne. (0.0e0, 0.0e0) ) then
         call accessdata_csp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         call get_descra(dspmtx%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         if (type.ne.'T') then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(x)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         z= conjg (x)
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rsbv_coo(dspmtx,x,ierr)
            case('CSC')
               call rsbv_csc(dspmtx,x,ierr)
            case('CSR')
               call rsbv_csr(dspmtx,x,ierr)
            case('DIA')
               call rsbv_dia(dspmtx,x,ierr)
            case('BCO')
               call rsbv_bco(dspmtx,x,ierr)
            case('BSC')
               call rsbv_bsc(dspmtx,x,ierr)
            case('BSR')
               call rsbv_bsr(dspmtx,x,ierr)
            case('BDI')
               call rsbv_bdi(dspmtx,x,ierr)
            case('VBR')
               call rsbv_vbr(dspmtx,x,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lsbv_coo(dspmtx,z,ierr)
            case('CSC')
               call lsbv_csc(dspmtx,z,ierr)
            case('CSR')
               call lsbv_csr(dspmtx,z,ierr)
            case('DIA')
               call lsbv_dia(dspmtx,z,ierr)
            case('BCO')
               call lsbv_bco(dspmtx,z,ierr)
            case('BSC')
               call lsbv_bsc(dspmtx,z,ierr)
            case('BSR')
               call lsbv_bsr(dspmtx,z,ierr)
            case('BDI')
               call lsbv_bdi(dspmtx,z,ierr)
            case('VBR')
               call lsbv_vbr(dspmtx,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
            ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            return
         end if            
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         x = alpha_work * x
      else
         x = alpha_work * (conjg (z))
      end if
      ierr = 0
      end subroutine cussv 
! **********************************************************************
! **********************************************************************
      subroutine zussv (a,x,ierr,transa,alpha)
      integer, intent(in) :: a
      complex(KIND=dp) , intent(inout) :: x(:)
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=dp) , intent(in), optional :: alpha
      integer :: transa_work
      complex(KIND=dp)  :: alpha_work
      complex(KIND=dp) , dimension(:), allocatable :: z
      type(zspmat ), pointer :: dspmtx
      character :: type
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
      if (alpha_work.ne. (0.0d0, 0.0d0) ) then
         call accessdata_zsp (dspmtx,a,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param 
            return
         end if
         call get_descra(dspmtx%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         if (type.ne.'T') then
            ierr = blas_error_param
            return
         end if
         allocate(z(size(x)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc 
            return
         end if
         z= conjg (x)
         select case(transa_work)
         case(ORIGIN_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call rsbv_coo(dspmtx,x,ierr)
            case('CSC')
               call rsbv_csc(dspmtx,x,ierr)
            case('CSR')
               call rsbv_csr(dspmtx,x,ierr)
            case('DIA')
               call rsbv_dia(dspmtx,x,ierr)
            case('BCO')
               call rsbv_bco(dspmtx,x,ierr)
            case('BSC')
               call rsbv_bsc(dspmtx,x,ierr)
            case('BSR')
               call rsbv_bsr(dspmtx,x,ierr)
            case('BDI')
               call rsbv_bdi(dspmtx,x,ierr)
            case('VBR')
               call rsbv_vbr(dspmtx,x,ierr)
            case default
               ierr = blas_error_param
            end select
         case(TRANSP_MATRIX)
            select case(dspmtx%FIDA)
            case('COO')
               call lsbv_coo(dspmtx,z,ierr)
            case('CSC')
               call lsbv_csc(dspmtx,z,ierr)
            case('CSR')
               call lsbv_csr(dspmtx,z,ierr)
            case('DIA')
               call lsbv_dia(dspmtx,z,ierr)
            case('BCO')
               call lsbv_bco(dspmtx,z,ierr)
            case('BSC')
               call lsbv_bsc(dspmtx,z,ierr)
            case('BSR')
               call lsbv_bsr(dspmtx,z,ierr)
            case('BDI')
               call lsbv_bdi(dspmtx,z,ierr)
            case('VBR')
               call lsbv_vbr(dspmtx,z,ierr)
            case default
               ierr = blas_error_param
            end select
         case default
            ierr = blas_error_param
         end select
         if (ierr.ne.0) then
            return
         end if            
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         x = alpha_work * x
      else
         x = alpha_work * (conjg (z))
      end if
      ierr = 0
      end subroutine zussv 
! **********************************************************************
! **********************************************************************
      end module  mod_ussv
