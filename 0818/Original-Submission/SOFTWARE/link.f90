      module  representation_of_data
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : THE PRINCIPAL DATA STRUCTURE
!                   Matrix data is stored in nodes of a linked list
!                   Node number is the handle number
!                   new: creates new node WITHOUT initialization
!                   del: frees unused memory, does NOT care if there
!                        is other memory that should be freed first
!                   accessdata: given a handle number, it returns a
!                        pointer to the matrix inside the relevant node
! **********************************************************************
      use types
      use properties
      implicit none
      interface accessdata 
         module procedure accessdata_isp
         module procedure accessdata_ssp
         module procedure accessdata_dsp
         module procedure accessdata_csp
         module procedure accessdata_zsp
      end interface
      type isp_linknode
        type(ispmat) :: contents
        integer :: number
        type(isp_linknode), pointer :: pntr
      end type isp_linknode
      type ssp_linknode
        type(sspmat) :: contents
        integer :: number
        type(ssp_linknode), pointer :: pntr
      end type ssp_linknode
      type dsp_linknode
        type(dspmat) :: contents
        integer :: number
        type(dsp_linknode), pointer :: pntr
      end type dsp_linknode
      type csp_linknode
        type(cspmat) :: contents
        integer :: number
        type(csp_linknode), pointer :: pntr
      end type csp_linknode
      type zsp_linknode
        type(zspmat) :: contents
        integer :: number
        type(zsp_linknode), pointer :: pntr
      end type zsp_linknode
      type(isp_linknode), pointer,SAVE,PRIVATE :: isp_first, isp_last
      type(ssp_linknode), pointer,SAVE,PRIVATE :: ssp_first, ssp_last
      type(dsp_linknode), pointer,SAVE,PRIVATE :: dsp_first, dsp_last
      type(csp_linknode), pointer,SAVE,PRIVATE :: csp_first, csp_last
      type(zsp_linknode), pointer,SAVE,PRIVATE :: zsp_first, zsp_last
      logical,SAVE,PRIVATE :: isp_init = .FALSE.
      logical,SAVE,PRIVATE :: ssp_init = .FALSE.
      logical,SAVE,PRIVATE :: dsp_init = .FALSE.
      logical,SAVE,PRIVATE :: csp_init = .FALSE.
      logical,SAVE,PRIVATE :: zsp_init = .FALSE.
      contains 
! **********************************************************************
! **********************************************************************
! *** Allocate new memory
      subroutine new_isp (nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(isp_linknode ), pointer :: help
      if(.not. isp_init ) then
         nullify(isp_first )
         isp_init  = .TRUE.
      endif
      if (.not.associated(isp_first )) then
         allocate(isp_first ,STAT=ierr)
         isp_first %number = ISP_MATRIX 
         nullify(isp_first %pntr)
         isp_last  => isp_first 
      else
         allocate(help,STAT=ierr)
         isp_last %pntr => help
         help%number = isp_last %number + no_of_types
         nullify(help%pntr)
         isp_last  => help
      end if
      nullify(isp_last %contents%A,isp_last %contents%IA1,&
              isp_last %contents%IA2,isp_last %contents%PB,&
              isp_last %contents%PE,isp_last %contents%BP1,&
              isp_last %contents%BP2)
      isp_last %contents%FIDA =''
      isp_last %contents%DESCRA =''
      isp_last %contents%INFOA = 0
      nmb = isp_last %number
      end subroutine new_isp 
! *** Deallocate unused memory
      subroutine del_isp (nmb,ierr)
      type(isp_linknode ), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (isp_first %number.eq.nmb) then
         ierr=0
         if (associated(isp_first ,isp_last )) then
            deallocate(isp_first )
            nullify(isp_first ,isp_last )
         else
            help2 => isp_first %pntr
            deallocate(isp_first )
            isp_first  => help2
         end if
      else
         help => isp_first 
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
            isp_last  => help
            nullify(isp_last %pntr)
            deallocate(help2)
         end if
      end if
      end subroutine del_isp 
! *** access contents for given number nmb
      subroutine accessdata_isp (dspmtx,nmb,ierr)
      type(ispmat ), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(isp_linknode ), pointer :: isp_handle 
      ierr = -1
      isp_handle  => isp_first 
      do while((isp_handle %number.ne.nmb).and.&
               (associated(isp_handle %pntr)))
         isp_handle  => isp_handle %pntr
      end do
      if (isp_handle %number.eq.nmb) then
         ierr = 0 
         dspmtx => isp_handle %contents
      else
         nullify(dspmtx)
      end if
      end subroutine accessdata_isp 
! **********************************************************************
! **********************************************************************
! *** Allocate new memory
      subroutine new_ssp (nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(ssp_linknode ), pointer :: help
      if(.not. ssp_init ) then
         nullify(ssp_first )
         ssp_init  = .TRUE.
      endif
      if (.not.associated(ssp_first )) then
         allocate(ssp_first ,STAT=ierr)
         ssp_first %number = SSP_MATRIX 
         nullify(ssp_first %pntr)
         ssp_last  => ssp_first 
      else
         allocate(help,STAT=ierr)
         ssp_last %pntr => help
         help%number = ssp_last %number + no_of_types
         nullify(help%pntr)
         ssp_last  => help
      end if
      nullify(ssp_last %contents%A,ssp_last %contents%IA1,&
              ssp_last %contents%IA2,ssp_last %contents%PB,&
              ssp_last %contents%PE,ssp_last %contents%BP1,&
              ssp_last %contents%BP2)
      ssp_last %contents%FIDA =''
      ssp_last %contents%DESCRA =''
      ssp_last %contents%INFOA = 0
      nmb = ssp_last %number
      end subroutine new_ssp 
! *** Deallocate unused memory
      subroutine del_ssp (nmb,ierr)
      type(ssp_linknode ), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (ssp_first %number.eq.nmb) then
         ierr=0
         if (associated(ssp_first ,ssp_last )) then
            deallocate(ssp_first )
            nullify(ssp_first ,ssp_last )
         else
            help2 => ssp_first %pntr
            deallocate(ssp_first )
            ssp_first  => help2
         end if
      else
         help => ssp_first 
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
            ssp_last  => help
            nullify(ssp_last %pntr)
            deallocate(help2)
         end if
      end if
      end subroutine del_ssp 
! *** access contents for given number nmb
      subroutine accessdata_ssp (dspmtx,nmb,ierr)
      type(sspmat ), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(ssp_linknode ), pointer :: ssp_handle 
      ierr = -1
      ssp_handle  => ssp_first 
      do while((ssp_handle %number.ne.nmb).and.&
               (associated(ssp_handle %pntr)))
         ssp_handle  => ssp_handle %pntr
      end do
      if (ssp_handle %number.eq.nmb) then
         ierr = 0 
         dspmtx => ssp_handle %contents
      else
         nullify(dspmtx)
      end if
      end subroutine accessdata_ssp 
! **********************************************************************
! **********************************************************************
! *** Allocate new memory
      subroutine new_dsp (nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(dsp_linknode ), pointer :: help
      if(.not. dsp_init ) then
         nullify(dsp_first )
         dsp_init  = .TRUE.
      endif
      if (.not.associated(dsp_first )) then
         allocate(dsp_first ,STAT=ierr)
         dsp_first %number = DSP_MATRIX 
         nullify(dsp_first %pntr)
         dsp_last  => dsp_first 
      else
         allocate(help,STAT=ierr)
         dsp_last %pntr => help
         help%number = dsp_last %number + no_of_types
         nullify(help%pntr)
         dsp_last  => help
      end if
      nullify(dsp_last %contents%A,dsp_last %contents%IA1,&
              dsp_last %contents%IA2,dsp_last %contents%PB,&
              dsp_last %contents%PE,dsp_last %contents%BP1,&
              dsp_last %contents%BP2)
      dsp_last %contents%FIDA =''
      dsp_last %contents%DESCRA =''
      dsp_last %contents%INFOA = 0
      nmb = dsp_last %number
      end subroutine new_dsp 
! *** Deallocate unused memory
      subroutine del_dsp (nmb,ierr)
      type(dsp_linknode ), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (dsp_first %number.eq.nmb) then
         ierr=0
         if (associated(dsp_first ,dsp_last )) then
            deallocate(dsp_first )
            nullify(dsp_first ,dsp_last )
         else
            help2 => dsp_first %pntr
            deallocate(dsp_first )
            dsp_first  => help2
         end if
      else
         help => dsp_first 
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
            dsp_last  => help
            nullify(dsp_last %pntr)
            deallocate(help2)
         end if
      end if
      end subroutine del_dsp 
! *** access contents for given number nmb
      subroutine accessdata_dsp (dspmtx,nmb,ierr)
      type(dspmat ), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(dsp_linknode ), pointer :: dsp_handle 
      ierr = -1
      dsp_handle  => dsp_first 
      do while((dsp_handle %number.ne.nmb).and.&
               (associated(dsp_handle %pntr)))
         dsp_handle  => dsp_handle %pntr
      end do
      if (dsp_handle %number.eq.nmb) then
         ierr = 0 
         dspmtx => dsp_handle %contents
      else
         nullify(dspmtx)
      end if
      end subroutine accessdata_dsp 
! **********************************************************************
! **********************************************************************
! *** Allocate new memory
      subroutine new_csp (nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(csp_linknode ), pointer :: help
      if(.not. csp_init ) then
         nullify(csp_first )
         csp_init  = .TRUE.
      endif
      if (.not.associated(csp_first )) then
         allocate(csp_first ,STAT=ierr)
         csp_first %number = CSP_MATRIX 
         nullify(csp_first %pntr)
         csp_last  => csp_first 
      else
         allocate(help,STAT=ierr)
         csp_last %pntr => help
         help%number = csp_last %number + no_of_types
         nullify(help%pntr)
         csp_last  => help
      end if
      nullify(csp_last %contents%A,csp_last %contents%IA1,&
              csp_last %contents%IA2,csp_last %contents%PB,&
              csp_last %contents%PE,csp_last %contents%BP1,&
              csp_last %contents%BP2)
      csp_last %contents%FIDA =''
      csp_last %contents%DESCRA =''
      csp_last %contents%INFOA = 0
      nmb = csp_last %number
      end subroutine new_csp 
! *** Deallocate unused memory
      subroutine del_csp (nmb,ierr)
      type(csp_linknode ), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (csp_first %number.eq.nmb) then
         ierr=0
         if (associated(csp_first ,csp_last )) then
            deallocate(csp_first )
            nullify(csp_first ,csp_last )
         else
            help2 => csp_first %pntr
            deallocate(csp_first )
            csp_first  => help2
         end if
      else
         help => csp_first 
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
            csp_last  => help
            nullify(csp_last %pntr)
            deallocate(help2)
         end if
      end if
      end subroutine del_csp 
! *** access contents for given number nmb
      subroutine accessdata_csp (dspmtx,nmb,ierr)
      type(cspmat ), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(csp_linknode ), pointer :: csp_handle 
      ierr = -1
      csp_handle  => csp_first 
      do while((csp_handle %number.ne.nmb).and.&
               (associated(csp_handle %pntr)))
         csp_handle  => csp_handle %pntr
      end do
      if (csp_handle %number.eq.nmb) then
         ierr = 0 
         dspmtx => csp_handle %contents
      else
         nullify(dspmtx)
      end if
      end subroutine accessdata_csp 
! **********************************************************************
! **********************************************************************
! *** Allocate new memory
      subroutine new_zsp (nmb,ierr)
      integer, intent(out) :: nmb,ierr
      type(zsp_linknode ), pointer :: help
      if(.not. zsp_init ) then
         nullify(zsp_first )
         zsp_init  = .TRUE.
      endif
      if (.not.associated(zsp_first )) then
         allocate(zsp_first ,STAT=ierr)
         zsp_first %number = ZSP_MATRIX 
         nullify(zsp_first %pntr)
         zsp_last  => zsp_first 
      else
         allocate(help,STAT=ierr)
         zsp_last %pntr => help
         help%number = zsp_last %number + no_of_types
         nullify(help%pntr)
         zsp_last  => help
      end if
      nullify(zsp_last %contents%A,zsp_last %contents%IA1,&
              zsp_last %contents%IA2,zsp_last %contents%PB,&
              zsp_last %contents%PE,zsp_last %contents%BP1,&
              zsp_last %contents%BP2)
      zsp_last %contents%FIDA =''
      zsp_last %contents%DESCRA =''
      zsp_last %contents%INFOA = 0
      nmb = zsp_last %number
      end subroutine new_zsp 
! *** Deallocate unused memory
      subroutine del_zsp (nmb,ierr)
      type(zsp_linknode ), pointer :: help,help2
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      ierr = -1
      if (zsp_first %number.eq.nmb) then
         ierr=0
         if (associated(zsp_first ,zsp_last )) then
            deallocate(zsp_first )
            nullify(zsp_first ,zsp_last )
         else
            help2 => zsp_first %pntr
            deallocate(zsp_first )
            zsp_first  => help2
         end if
      else
         help => zsp_first 
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
            zsp_last  => help
            nullify(zsp_last %pntr)
            deallocate(help2)
         end if
      end if
      end subroutine del_zsp 
! *** access contents for given number nmb
      subroutine accessdata_zsp (dspmtx,nmb,ierr)
      type(zspmat ), pointer :: dspmtx
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(zsp_linknode ), pointer :: zsp_handle 
      ierr = -1
      zsp_handle  => zsp_first 
      do while((zsp_handle %number.ne.nmb).and.&
               (associated(zsp_handle %pntr)))
         zsp_handle  => zsp_handle %pntr
      end do
      if (zsp_handle %number.eq.nmb) then
         ierr = 0 
         dspmtx => zsp_handle %contents
      else
         nullify(dspmtx)
      end if
      end subroutine accessdata_zsp 
! **********************************************************************
! **********************************************************************
      end module representation_of_data
