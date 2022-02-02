! *** Allocate new memory
      subroutine NEW_XSP(nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(XSP_LINKNODE), pointer :: help

      if(.not.XSP_INIT) then
         nullify(XSP_FIRST)
         XSP_INIT = .TRUE.
      endif

      if (.not.associated(XSP_FIRST)) then
         allocate(XSP_FIRST,STAT=ierr)
         XSP_FIRST%number = XSP_MATRIX
         nullify(XSP_FIRST%pntr)
         XSP_LAST => XSP_FIRST
      else
         allocate(help,STAT=ierr)
         XSP_LAST%pntr => help
         help%number = XSP_LAST%number + no_of_types
         nullify(help%pntr)
         XSP_LAST => help
      end if
      nullify(XSP_LAST%contents%A,XSP_LAST%contents%IA1,&
              XSP_LAST%contents%IA2,XSP_LAST%contents%PB,&
              XSP_LAST%contents%PE,XSP_LAST%contents%BP1,&
              XSP_LAST%contents%BP2)
      XSP_LAST%contents%FIDA =''
      XSP_LAST%contents%DESCRA =''
      XSP_LAST%contents%INFOA = 0
      nmb = XSP_LAST%number
      end subroutine NEW_XSP

! *** Deallocate unused memory
      subroutine DEL_XSP(nmb,ierr)
      type(XSP_LINKNODE), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (XSP_FIRST%number.eq.nmb) then
         ierr=0
         if (associated(XSP_FIRST,XSP_LAST)) then
            deallocate(XSP_FIRST)
            nullify(XSP_FIRST,XSP_LAST)
         else
            help2 => XSP_FIRST%pntr
            deallocate(XSP_FIRST)
            XSP_FIRST => help2
         end if
      else
         help => XSP_FIRST
         do while((ierr.eq.-1).and.(associated(help%pntr%pntr)))
            if (help%pntr%number.eq.nmb) then
               help2 => help%pntr
               help%pntr => help%pntr%pntr
               deallocate(help2)
               ierr = 0 
            else
               help => help%pntr
            end if
         end do
         if((ierr.eq.-1).and.(help%pntr%number.eq.nmb)) then
            ierr = 0
            help2 => help%pntr
            XSP_LAST => help
            nullify(XSP_LAST%pntr)
            deallocate(help2)
         end if
      end if
      end subroutine DEL_XSP

! *** access contents for given number nmb
      subroutine ACCESSDATA_XSP(dspmtx,nmb,ierr)
      type(XSPMAT), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(XSP_LINKNODE), pointer :: XSP_HANDLE
      ierr = -1
      XSP_HANDLE => XSP_FIRST
      do while((XSP_HANDLE%number.ne.nmb).and.&
               (associated(XSP_HANDLE%pntr)))
         XSP_HANDLE => XSP_HANDLE%pntr
      end do
      if (XSP_HANDLE%number.eq.nmb) then
         ierr = 0 
         dspmtx => XSP_HANDLE%contents
      else
         nullify(dspmtx)
      end if
      end subroutine ACCESSDATA_XSP
