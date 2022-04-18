      module mod_uscr_dia
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'DIA'-FORMAT
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      interface uscr_dia
        module procedure iuscr_dia
        module procedure suscr_dia
        module procedure duscr_dia
        module procedure cuscr_dia
        module procedure zuscr_dia
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_dia (m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      integer , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(ispmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if
      call new_isp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'DIA'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine iuscr_dia 
! **********************************************************************
! **********************************************************************
      subroutine suscr_dia (m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      real(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(sspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if
      call new_ssp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'DIA'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine suscr_dia 
! **********************************************************************
! **********************************************************************
      subroutine duscr_dia (m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      real(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(dspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if
      call new_dsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'DIA'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine duscr_dia 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_dia (m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      complex(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(cspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if
      call new_csp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'DIA'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine cuscr_dia 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_dia (m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      complex(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(zspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if
      call new_zsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'DIA'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine zuscr_dia 
! **********************************************************************
! **********************************************************************
      end module mod_uscr_dia
