DCOMPLEX function XUSDOT(x,indx,y,conj)
    implicit none 
    integer,dimension(:),intent(in) :: indx
    DCOMPLEX ,dimension(:),intent(in) ::x,y
    DCOMPLEX ,dimension(:),allocatable :: zy
    integer,optional ::conj
    integer                        ::t
    intrinsic dot_product,conjg,cmplx
   
    t=size(indx)
    if(t.le.0) then
       XUSDOT=0.
    else
       allocate(zy(t))
       zy= y(indx)
       if(present(conj)) then
          XUSDOT=dot_product(x,zy)
       else
          XUSDOT=dot_product(conjg(cmplx(x)),zy)
       end if
       deallocate(zy)
    end if
  end function XUSDOT
   
