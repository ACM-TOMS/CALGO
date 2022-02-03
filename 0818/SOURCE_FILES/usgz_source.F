subroutine XUSGZ(y,x,indx)
     DCOMPLEX ,dimension(:),intent(out) ::x
     DCOMPLEX ,dimension(:),intent(inout) ::y
     integer,dimension(:),intent(in)::indx   
     integer  ::i,t
     
     t=size(indx)
     if(t.gt.0) then               
        call  usga(y,x,indx)
        do i=1,t
           y(indx(i))=0
        end do
     end if
   end subroutine XUSGZ
