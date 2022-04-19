      module mod_usdot
      use blas_sparse_namedconstants
      interface usdot
      module procedure iusdot
      module procedure susdot
      module procedure dusdot
      module procedure cusdot
      module procedure zusdot
      end interface
      contains
! **********************************************************************
! **********************************************************************
integer  function iusdot (x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    integer  ,dimension(:),intent(in) ::x,y
    integer  ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
    t=size(indx)
    if(t.le.0) then
       iusdot =0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          iusdot =dot_product(x,zy)
       else
          iusdot =dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function iusdot 
! **********************************************************************
! **********************************************************************
real(KIND=sp)  function susdot (x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    real(KIND=sp)  ,dimension(:),intent(in) ::x,y
    real(KIND=sp)  ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
    t=size(indx)
    if(t.le.0) then
       susdot =0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          susdot =dot_product(x,zy)
       else
          susdot =dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function susdot 
! **********************************************************************
! **********************************************************************
real(KIND=dp)  function dusdot (x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    real(KIND=dp)  ,dimension(:),intent(in) ::x,y
    real(KIND=dp)  ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
    t=size(indx)
    if(t.le.0) then
       dusdot =0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          dusdot =dot_product(x,zy)
       else
          dusdot =dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function dusdot 
! **********************************************************************
! **********************************************************************
complex(KIND=sp)  function cusdot (x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    complex(KIND=sp)  ,dimension(:),intent(in) ::x,y
    complex(KIND=sp)  ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
    t=size(indx)
    if(t.le.0) then
       cusdot =0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          cusdot =dot_product(x,zy)
       else
          cusdot =dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function cusdot 
! **********************************************************************
! **********************************************************************
complex(KIND=dp)  function zusdot (x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    complex(KIND=dp)  ,dimension(:),intent(in) ::x,y
    complex(KIND=dp)  ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
    t=size(indx)
    if(t.le.0) then
       zusdot =0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          zusdot =dot_product(x,zy)
       else
          zusdot =dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function zusdot 
! **********************************************************************
! **********************************************************************
      end module mod_usdot
