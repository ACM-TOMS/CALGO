      module mod_ussc
      use blas_sparse_namedconstants
      interface ussc
      module procedure iussc
      module procedure sussc
      module procedure dussc
      module procedure cussc
      module procedure zussc
      end interface
      contains
! **********************************************************************
! **********************************************************************
  subroutine iussc (x,y,indx)
 integer ,dimension(:),intent(in) ::x
 integer  ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  end subroutine iussc 
! **********************************************************************
! **********************************************************************
  subroutine sussc (x,y,indx)
 real(KIND=sp) ,dimension(:),intent(in) ::x
 real(KIND=sp)  ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  end subroutine sussc 
! **********************************************************************
! **********************************************************************
  subroutine dussc (x,y,indx)
 real(KIND=dp) ,dimension(:),intent(in) ::x
 real(KIND=dp)  ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  end subroutine dussc 
! **********************************************************************
! **********************************************************************
  subroutine cussc (x,y,indx)
 complex(KIND=sp) ,dimension(:),intent(in) ::x
 complex(KIND=sp)  ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  end subroutine cussc 
! **********************************************************************
! **********************************************************************
  subroutine zussc (x,y,indx)
 complex(KIND=dp) ,dimension(:),intent(in) ::x
 complex(KIND=dp)  ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  end subroutine zussc 
! **********************************************************************
! **********************************************************************
      end module mod_ussc
