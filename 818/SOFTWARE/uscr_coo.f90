      module mod_uscr_coo
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'COO'-FORMAT
! **********************************************************************
      use representation_of_data
      use properties
      use mod_usds
      implicit none
      interface uscr_coo
        module procedure iuscr_coo
        module procedure suscr_coo
        module procedure duscr_coo
        module procedure cuscr_coo
        module procedure zuscr_coo
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,nnz,prpty
      integer, dimension(:), intent(inout),target :: indx,jndx
      integer , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options,base
      logical :: COPY
      character :: message
      type(ispmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
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
      dsp_data%FIDA = 'COO'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
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
       if((nnz.ne.size(indx)).or.(nnz.ne.size(jndx)).or.&
         (nnz.ne.size(val)).or.(minval(indx).lt.base).or.&
         (minval(jndx).lt.base).or.(maxval(indx).gt.m-1+base).or.&
         (maxval(jndx).gt.n-1+base)) then
         call usds(nmb,ierr)
         ierr = blas_error_param
         return
       end if
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
         dsp_data%IA1 => indx
         dsp_data%IA2 => jndx
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(indx)),&
                             dsp_data%IA2(size(jndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = indx
         dsp_data%IA2 = jndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine iuscr_coo 
! **********************************************************************
! **********************************************************************
      subroutine suscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,nnz,prpty
      integer, dimension(:), intent(inout),target :: indx,jndx
      real(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options,base
      logical :: COPY
      character :: message
      type(sspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
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
      dsp_data%FIDA = 'COO'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
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
       if((nnz.ne.size(indx)).or.(nnz.ne.size(jndx)).or.&
         (nnz.ne.size(val)).or.(minval(indx).lt.base).or.&
         (minval(jndx).lt.base).or.(maxval(indx).gt.m-1+base).or.&
         (maxval(jndx).gt.n-1+base)) then
         call usds(nmb,ierr)
         ierr = blas_error_param
         return
       end if
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
         dsp_data%IA1 => indx
         dsp_data%IA2 => jndx
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(indx)),&
                             dsp_data%IA2(size(jndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = indx
         dsp_data%IA2 = jndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine suscr_coo 
! **********************************************************************
! **********************************************************************
      subroutine duscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,nnz,prpty
      integer, dimension(:), intent(inout),target :: indx,jndx
      real(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options,base
      logical :: COPY
      character :: message
      type(dspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
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
      dsp_data%FIDA = 'COO'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
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
       if((nnz.ne.size(indx)).or.(nnz.ne.size(jndx)).or.&
         (nnz.ne.size(val)).or.(minval(indx).lt.base).or.&
         (minval(jndx).lt.base).or.(maxval(indx).gt.m-1+base).or.&
         (maxval(jndx).gt.n-1+base)) then
         call usds(nmb,ierr)
         ierr = blas_error_param
         return
       end if
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
         dsp_data%IA1 => indx
         dsp_data%IA2 => jndx
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(indx)),&
                             dsp_data%IA2(size(jndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = indx
         dsp_data%IA2 = jndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine duscr_coo 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,nnz,prpty
      integer, dimension(:), intent(inout),target :: indx,jndx
      complex(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options,base
      logical :: COPY
      character :: message
      type(cspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
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
      dsp_data%FIDA = 'COO'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
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
       if((nnz.ne.size(indx)).or.(nnz.ne.size(jndx)).or.&
         (nnz.ne.size(val)).or.(minval(indx).lt.base).or.&
         (minval(jndx).lt.base).or.(maxval(indx).gt.m-1+base).or.&
         (maxval(jndx).gt.n-1+base)) then
         call usds(nmb,ierr)
         ierr = blas_error_param
         return
       end if
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
         dsp_data%IA1 => indx
         dsp_data%IA2 => jndx
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(indx)),&
                             dsp_data%IA2(size(jndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = indx
         dsp_data%IA2 = jndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine cuscr_coo 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,nnz,prpty
      integer, dimension(:), intent(inout),target :: indx,jndx
      complex(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options,base
      logical :: COPY
      character :: message
      type(zspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
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
      dsp_data%FIDA = 'COO'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
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
       if((nnz.ne.size(indx)).or.(nnz.ne.size(jndx)).or.&
         (nnz.ne.size(val)).or.(minval(indx).lt.base).or.&
         (minval(jndx).lt.base).or.(maxval(indx).gt.m-1+base).or.&
         (maxval(jndx).gt.n-1+base)) then
         call usds(nmb,ierr)
         ierr = blas_error_param
         return
       end if
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
         dsp_data%IA1 => indx
         dsp_data%IA2 => jndx
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(indx)),&
                             dsp_data%IA2(size(jndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = indx
         dsp_data%IA2 = jndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine zuscr_coo 
! **********************************************************************
! **********************************************************************
      end module mod_uscr_coo
