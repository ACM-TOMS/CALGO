      module mod_usgz
      use mod_usga
      use blas_sparse_namedconstants
      interface usgz
      module procedure iusgz
      module procedure susgz
      module procedure dusgz
      module procedure cusgz
      module procedure zusgz
      end interface
      contains
! **********************************************************************
! **********************************************************************
subroutine iusgz (y,x,indx)
     integer  ,dimension(:),intent(out) ::x
     integer  ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine iusgz 
! **********************************************************************
! **********************************************************************
subroutine susgz (y,x,indx)
     real(KIND=sp)  ,dimension(:),intent(out) ::x
     real(KIND=sp)  ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine susgz 
! **********************************************************************
! **********************************************************************
subroutine dusgz (y,x,indx)
     real(KIND=dp)  ,dimension(:),intent(out) ::x
     real(KIND=dp)  ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine dusgz 
! **********************************************************************
! **********************************************************************
subroutine cusgz (y,x,indx)
     complex(KIND=sp)  ,dimension(:),intent(out) ::x
     complex(KIND=sp)  ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine cusgz 
! **********************************************************************
! **********************************************************************
subroutine zusgz (y,x,indx)
     complex(KIND=dp)  ,dimension(:),intent(out) ::x
     complex(KIND=dp)  ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine zusgz 
! **********************************************************************
! **********************************************************************
      end module mod_usgz
