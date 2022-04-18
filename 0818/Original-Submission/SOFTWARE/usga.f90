      module mod_usga
      use blas_sparse_namedconstants
      interface usga
      module procedure iusga
      module procedure susga
      module procedure dusga
      module procedure cusga
      module procedure zusga
      end interface
      contains
! **********************************************************************
! **********************************************************************
 subroutine iusga (y,x,indx)
    integer  ,dimension(:),intent(inout) ::x
    integer  ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine iusga 
! **********************************************************************
! **********************************************************************
 subroutine susga (y,x,indx)
    real(KIND=sp)  ,dimension(:),intent(inout) ::x
    real(KIND=sp)  ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine susga 
! **********************************************************************
! **********************************************************************
 subroutine dusga (y,x,indx)
    real(KIND=dp)  ,dimension(:),intent(inout) ::x
    real(KIND=dp)  ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine dusga 
! **********************************************************************
! **********************************************************************
 subroutine cusga (y,x,indx)
    complex(KIND=sp)  ,dimension(:),intent(inout) ::x
    complex(KIND=sp)  ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine cusga 
! **********************************************************************
! **********************************************************************
 subroutine zusga (y,x,indx)
    complex(KIND=dp)  ,dimension(:),intent(inout) ::x
    complex(KIND=dp)  ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine zusga 
! **********************************************************************
! **********************************************************************
      end module mod_usga
