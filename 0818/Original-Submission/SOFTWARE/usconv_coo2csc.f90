      module mod_usconv_coo2csc
      use properties
      use mod_conv_tools
      use representation_of_data
      contains
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine iusconv_coo2csc (a,ierr) 
      integer,intent(inout) :: a 
      type( ispmat ), pointer :: dspmtx
      integer ,intent(inout)::ierr
      integer :: res
      ierr=-1
      call accessdata(dspmtx,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dspmtx%FIDA=='COO') then
         call get_infoa(dspmtx%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dspmtx%FIDA='CSC'
            allocate(dspmtx%PB(dspmtx%K))
            allocate(dspmtx%PE(dspmtx%K))
            call ipre_usconv_coo2csc ( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine iusconv_coo2csc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine susconv_coo2csc (a,ierr) 
      integer,intent(inout) :: a 
      type( sspmat ), pointer :: dspmtx
      integer ,intent(inout)::ierr
      integer :: res
      ierr=-1
      call accessdata(dspmtx,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dspmtx%FIDA=='COO') then
         call get_infoa(dspmtx%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dspmtx%FIDA='CSC'
            allocate(dspmtx%PB(dspmtx%K))
            allocate(dspmtx%PE(dspmtx%K))
            call spre_usconv_coo2csc ( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine susconv_coo2csc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine dusconv_coo2csc (a,ierr) 
      integer,intent(inout) :: a 
      type( dspmat ), pointer :: dspmtx
      integer ,intent(inout)::ierr
      integer :: res
      ierr=-1
      call accessdata(dspmtx,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dspmtx%FIDA=='COO') then
         call get_infoa(dspmtx%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dspmtx%FIDA='CSC'
            allocate(dspmtx%PB(dspmtx%K))
            allocate(dspmtx%PE(dspmtx%K))
            call dpre_usconv_coo2csc ( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine dusconv_coo2csc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine cusconv_coo2csc (a,ierr) 
      integer,intent(inout) :: a 
      type( cspmat ), pointer :: dspmtx
      integer ,intent(inout)::ierr
      integer :: res
      ierr=-1
      call accessdata(dspmtx,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dspmtx%FIDA=='COO') then
         call get_infoa(dspmtx%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dspmtx%FIDA='CSC'
            allocate(dspmtx%PB(dspmtx%K))
            allocate(dspmtx%PE(dspmtx%K))
            call cpre_usconv_coo2csc ( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine cusconv_coo2csc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine zusconv_coo2csc (a,ierr) 
      integer,intent(inout) :: a 
      type( zspmat ), pointer :: dspmtx
      integer ,intent(inout)::ierr
      integer :: res
      ierr=-1
      call accessdata(dspmtx,a,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if(dspmtx%FIDA=='COO') then
         call get_infoa(dspmtx%INFOA ,'c',res,ierr)
         if(res.eq.COP_OF_SOURCE) then
            dspmtx%FIDA='CSC'
            allocate(dspmtx%PB(dspmtx%K))
            allocate(dspmtx%PE(dspmtx%K))
            call zpre_usconv_coo2csc ( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine zusconv_coo2csc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      end module mod_usconv_coo2csc
