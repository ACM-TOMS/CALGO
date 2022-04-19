      module mod_ussm
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : TRI. SOLVE, CHOOSES APPROPRIATE SUBROUTINES
! **********************************************************************
      use representation_of_data
      use properties
      use mod_sbv
      implicit none
      interface ussm
        module procedure iussm
        module procedure sussm
        module procedure dussm
        module procedure cussm
        module procedure zussm
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iussm (a,b,ierr,transa,alpha)
      integer, intent(in) :: a
      integer , dimension(:,:), intent(inout) :: b
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      integer , intent(in), optional :: alpha
      integer :: transa_work,i
      integer  :: alpha_work
      integer , dimension(:,:), allocatable :: z
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
      if (alpha_work.ne.0.) then
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
         allocate(z(size(b,1),size(b,2)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         z=  (b)
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rsbv_coo(dspmtx,b(:,i),ierr)
               case('CSC')
                  call rsbv_csc(dspmtx,b(:,i),ierr)
               case('CSR')
                  call rsbv_csr(dspmtx,b(:,i),ierr)
               case('DIA')
                  call rsbv_dia(dspmtx,b(:,i),ierr)
               case('BCO')
                  call rsbv_bco(dspmtx,b(:,i),ierr)
               case('BSC')
                  call rsbv_bsc(dspmtx,b(:,i),ierr)
               case('BSR')
                  call rsbv_bsr(dspmtx,b(:,i),ierr)
               case('BDI')
                  call rsbv_bdi(dspmtx,b(:,i),ierr)
               case('VBR')
                  call rsbv_vbr(dspmtx,b(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lsbv_coo(dspmtx,z(:,i),ierr)
               case('CSC')
                  call lsbv_csc(dspmtx,z(:,i),ierr)
               case('CSR')
                  call lsbv_csr(dspmtx,z(:,i),ierr)
               case('DIA')
                  call lsbv_dia(dspmtx,z(:,i),ierr)
               case('BCO')
                  call lsbv_bco(dspmtx,z(:,i),ierr)
               case('BSC')
                  call lsbv_bsc(dspmtx,z(:,i),ierr)
               case('BSR')
                  call lsbv_bsr(dspmtx,z(:,i),ierr)
               case('BDI')
                  call lsbv_bdi(dspmtx,z(:,i),ierr)
               case('VBR')
                  call lsbv_vbr(dspmtx,z(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               return
            end if            
         end do
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         b = alpha_work * b
      else
         b = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine iussm 
! **********************************************************************
! **********************************************************************
      subroutine sussm (a,b,ierr,transa,alpha)
      integer, intent(in) :: a
      real(KIND=sp) , dimension(:,:), intent(inout) :: b
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=sp) , intent(in), optional :: alpha
      integer :: transa_work,i
      real(KIND=sp)  :: alpha_work
      real(KIND=sp) , dimension(:,:), allocatable :: z
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
      if (alpha_work.ne.0.) then
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
         allocate(z(size(b,1),size(b,2)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         z=  (b)
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rsbv_coo(dspmtx,b(:,i),ierr)
               case('CSC')
                  call rsbv_csc(dspmtx,b(:,i),ierr)
               case('CSR')
                  call rsbv_csr(dspmtx,b(:,i),ierr)
               case('DIA')
                  call rsbv_dia(dspmtx,b(:,i),ierr)
               case('BCO')
                  call rsbv_bco(dspmtx,b(:,i),ierr)
               case('BSC')
                  call rsbv_bsc(dspmtx,b(:,i),ierr)
               case('BSR')
                  call rsbv_bsr(dspmtx,b(:,i),ierr)
               case('BDI')
                  call rsbv_bdi(dspmtx,b(:,i),ierr)
               case('VBR')
                  call rsbv_vbr(dspmtx,b(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lsbv_coo(dspmtx,z(:,i),ierr)
               case('CSC')
                  call lsbv_csc(dspmtx,z(:,i),ierr)
               case('CSR')
                  call lsbv_csr(dspmtx,z(:,i),ierr)
               case('DIA')
                  call lsbv_dia(dspmtx,z(:,i),ierr)
               case('BCO')
                  call lsbv_bco(dspmtx,z(:,i),ierr)
               case('BSC')
                  call lsbv_bsc(dspmtx,z(:,i),ierr)
               case('BSR')
                  call lsbv_bsr(dspmtx,z(:,i),ierr)
               case('BDI')
                  call lsbv_bdi(dspmtx,z(:,i),ierr)
               case('VBR')
                  call lsbv_vbr(dspmtx,z(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               return
            end if            
         end do
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         b = alpha_work * b
      else
         b = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine sussm 
! **********************************************************************
! **********************************************************************
      subroutine dussm (a,b,ierr,transa,alpha)
      integer, intent(in) :: a
      real(KIND=dp) , dimension(:,:), intent(inout) :: b
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      real(KIND=dp) , intent(in), optional :: alpha
      integer :: transa_work,i
      real(KIND=dp)  :: alpha_work
      real(KIND=dp) , dimension(:,:), allocatable :: z
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
      if (alpha_work.ne.0.) then
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
         allocate(z(size(b,1),size(b,2)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         z=  (b)
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rsbv_coo(dspmtx,b(:,i),ierr)
               case('CSC')
                  call rsbv_csc(dspmtx,b(:,i),ierr)
               case('CSR')
                  call rsbv_csr(dspmtx,b(:,i),ierr)
               case('DIA')
                  call rsbv_dia(dspmtx,b(:,i),ierr)
               case('BCO')
                  call rsbv_bco(dspmtx,b(:,i),ierr)
               case('BSC')
                  call rsbv_bsc(dspmtx,b(:,i),ierr)
               case('BSR')
                  call rsbv_bsr(dspmtx,b(:,i),ierr)
               case('BDI')
                  call rsbv_bdi(dspmtx,b(:,i),ierr)
               case('VBR')
                  call rsbv_vbr(dspmtx,b(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lsbv_coo(dspmtx,z(:,i),ierr)
               case('CSC')
                  call lsbv_csc(dspmtx,z(:,i),ierr)
               case('CSR')
                  call lsbv_csr(dspmtx,z(:,i),ierr)
               case('DIA')
                  call lsbv_dia(dspmtx,z(:,i),ierr)
               case('BCO')
                  call lsbv_bco(dspmtx,z(:,i),ierr)
               case('BSC')
                  call lsbv_bsc(dspmtx,z(:,i),ierr)
               case('BSR')
                  call lsbv_bsr(dspmtx,z(:,i),ierr)
               case('BDI')
                  call lsbv_bdi(dspmtx,z(:,i),ierr)
               case('VBR')
                  call lsbv_vbr(dspmtx,z(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               return
            end if            
         end do
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         b = alpha_work * b
      else
         b = alpha_work * ( (z))
      end if
      ierr = 0
      end subroutine dussm 
! **********************************************************************
! **********************************************************************
      subroutine cussm (a,b,ierr,transa,alpha)
      integer, intent(in) :: a
      complex(KIND=sp) , dimension(:,:), intent(inout) :: b
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=sp) , intent(in), optional :: alpha
      integer :: transa_work,i
      complex(KIND=sp)  :: alpha_work
      complex(KIND=sp) , dimension(:,:), allocatable :: z
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
      if (alpha_work.ne.0.) then
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
         allocate(z(size(b,1),size(b,2)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         z= conjg (b)
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rsbv_coo(dspmtx,b(:,i),ierr)
               case('CSC')
                  call rsbv_csc(dspmtx,b(:,i),ierr)
               case('CSR')
                  call rsbv_csr(dspmtx,b(:,i),ierr)
               case('DIA')
                  call rsbv_dia(dspmtx,b(:,i),ierr)
               case('BCO')
                  call rsbv_bco(dspmtx,b(:,i),ierr)
               case('BSC')
                  call rsbv_bsc(dspmtx,b(:,i),ierr)
               case('BSR')
                  call rsbv_bsr(dspmtx,b(:,i),ierr)
               case('BDI')
                  call rsbv_bdi(dspmtx,b(:,i),ierr)
               case('VBR')
                  call rsbv_vbr(dspmtx,b(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lsbv_coo(dspmtx,z(:,i),ierr)
               case('CSC')
                  call lsbv_csc(dspmtx,z(:,i),ierr)
               case('CSR')
                  call lsbv_csr(dspmtx,z(:,i),ierr)
               case('DIA')
                  call lsbv_dia(dspmtx,z(:,i),ierr)
               case('BCO')
                  call lsbv_bco(dspmtx,z(:,i),ierr)
               case('BSC')
                  call lsbv_bsc(dspmtx,z(:,i),ierr)
               case('BSR')
                  call lsbv_bsr(dspmtx,z(:,i),ierr)
               case('BDI')
                  call lsbv_bdi(dspmtx,z(:,i),ierr)
               case('VBR')
                  call lsbv_vbr(dspmtx,z(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               return
            end if            
         end do
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         b = alpha_work * b
      else
         b = alpha_work * (conjg (z))
      end if
      ierr = 0
      end subroutine cussm 
! **********************************************************************
! **********************************************************************
      subroutine zussm (a,b,ierr,transa,alpha)
      integer, intent(in) :: a
      complex(KIND=dp) , dimension(:,:), intent(inout) :: b
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      complex(KIND=dp) , intent(in), optional :: alpha
      integer :: transa_work,i
      complex(KIND=dp)  :: alpha_work
      complex(KIND=dp) , dimension(:,:), allocatable :: z
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
      if (alpha_work.ne.0.) then
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
         allocate(z(size(b,1),size(b,2)),STAT=ierr)
         if (ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         z= conjg (b)
         do i = 1,size(b,2)
            select case(transa_work)
            case(ORIGIN_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call rsbv_coo(dspmtx,b(:,i),ierr)
               case('CSC')
                  call rsbv_csc(dspmtx,b(:,i),ierr)
               case('CSR')
                  call rsbv_csr(dspmtx,b(:,i),ierr)
               case('DIA')
                  call rsbv_dia(dspmtx,b(:,i),ierr)
               case('BCO')
                  call rsbv_bco(dspmtx,b(:,i),ierr)
               case('BSC')
                  call rsbv_bsc(dspmtx,b(:,i),ierr)
               case('BSR')
                  call rsbv_bsr(dspmtx,b(:,i),ierr)
               case('BDI')
                  call rsbv_bdi(dspmtx,b(:,i),ierr)
               case('VBR')
                  call rsbv_vbr(dspmtx,b(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case(TRANSP_MATRIX)
               select case(dspmtx%FIDA)
               case('COO')
                  call lsbv_coo(dspmtx,z(:,i),ierr)
               case('CSC')
                  call lsbv_csc(dspmtx,z(:,i),ierr)
               case('CSR')
                  call lsbv_csr(dspmtx,z(:,i),ierr)
               case('DIA')
                  call lsbv_dia(dspmtx,z(:,i),ierr)
               case('BCO')
                  call lsbv_bco(dspmtx,z(:,i),ierr)
               case('BSC')
                  call lsbv_bsc(dspmtx,z(:,i),ierr)
               case('BSR')
                  call lsbv_bsr(dspmtx,z(:,i),ierr)
               case('BDI')
                  call lsbv_bdi(dspmtx,z(:,i),ierr)
               case('VBR')
                  call lsbv_vbr(dspmtx,z(:,i),ierr)
               case default
                  ierr = blas_error_param
               end select
            case default
               ierr = blas_error_param
            end select
            if (ierr.ne.0) then
               return
            end if            
         end do
      end if
      if(transa_work.eq.ORIGIN_MATRIX) then
         b = alpha_work * b
      else
         b = alpha_work * (conjg (z))
      end if
      ierr = 0
      end subroutine zussm 
! **********************************************************************
! **********************************************************************
      end module  mod_ussm
