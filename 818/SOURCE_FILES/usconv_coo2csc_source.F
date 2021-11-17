      subroutine XCOO2CSC(a,ierr) 
      integer,intent(inout) :: a 
      type( XSPMAT), pointer :: dspmtx
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
            call XPRE_COO2CSC( dspmtx%A, dspmtx%IA1, dspmtx%IA2,& 
                dspmtx%K, dspmtx%PB, dspmtx%PE)
            nullify(dspmtx%IA2)
         end if
      else 
         ierr = blas_error_param
         return
      end if
      end subroutine XCOO2CSC
 
