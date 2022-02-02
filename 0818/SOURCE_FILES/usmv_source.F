      subroutine XUSMV(a,x,y,ierr,transa,alpha)
      implicit none
      integer, intent(in) :: a
      DCOMPLEX, dimension(:), intent(in) :: x      
      DCOMPLEX, dimension(:), intent(inout) :: y
      integer, intent(out) :: ierr
      integer, intent(in), optional :: transa
      DCOMPLEX, intent(in), optional :: alpha
      DCOMPLEX, dimension(:), allocatable :: z
      type(XSPMAT), pointer :: dspmtx
      integer :: transa_work
      DCOMPLEX :: alpha_work

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
      if (alpha_work.eq.ZZERO) then
         !no matrix multiplication necessary
      else
         call ACCESSDATA_XSP(dspmtx,a,ierr)
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
               call lmbv_coo(dspmtx,CONJUG(x),z,ierr)
            case('CSC')
               call lmbv_csc(dspmtx,CONJUG(x),z,ierr)
            case('CSR')
               call lmbv_csr(dspmtx,CONJUG(x),z,ierr)
            case('DIA')
               call lmbv_dia(dspmtx,CONJUG(x),z,ierr)
            case('BCO')
               call lmbv_bco(dspmtx,CONJUG(x),z,ierr)
            case('BSC')
               call lmbv_bsc(dspmtx,CONJUG(x),z,ierr)
            case('BSR')
               call lmbv_bsr(dspmtx,CONJUG(x),z,ierr)
            case('BDI')
               call lmbv_bdi(dspmtx,CONJUG(x),z,ierr)
            case('VBR')
               call lmbv_vbr(dspmtx,CONJUG(x),z,ierr)
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
            y = alpha_work * (CONJUG(z)) + y
         end if
         deallocate(z,STAT=ierr)
      end if
      ierr = 0
      end subroutine XUSMV
