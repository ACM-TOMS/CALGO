Program Driver

!use numerical_libraries
use constants
use Blanch_Complex



implicit none
integer:: N_Max,k_max,parity,der,N,choice, radial_kind, sign
real(kind=double):: u,x
complex (kind=double):: q,Se, norm, Re
complex(kind=double), dimension(0:MAX_) :: D_m
character(80):: filename


!Sample computation of Mathieu angular function

parity=0
der=0
n=1
q=-0.9D0+j*0.9D0
x=1.0D0
choice=1

Se=MathieuAngular(parity,der,n,Q,x,k_max, choice, D_m, norm)

print *, "The value of the Mathieu angular function of even parity of"
print *, "order n=",n
print *, "for q=",q
print *, "and argument x=",x
print *, "is ", Se


!Sample computation of Mathieu radial function

radial_kind=4
sign=1

Re=MathieuRadial(parity,radial_kind,der,n,Q,x,k_max, choice, sign)

print *, "The value of the Mathieu radial function of even parity of"
print *, "order n=",n
print *, "kind=", radial_kind
print *, "for q=",q
print *, "and argument x=",x
print *, "is ", Re


!Sample of call to the subroutine that computes the wronskian check
N_max=20
u=2.0D0
filename="Wronskian_Table_I_Linux_Cluster.dat"
call Test_Wronskian(q, N_max, u,filename)



end program Driver


!************************************************************************************
subroutine Test_Wronskian(q, N_max, u,filename)

! Evaluates the Wronskian of the radial Mathieu functions computed by Mathieu_Radial_2


use constants
use Blanch_Complex
implicit none

complex(kind=double), intent(in) ::q
integer, intent(in):: N_max

real(kind=double), intent(in):: u
character(80):: filename

integer:: k_max, sign, choice, n
complex(kind=double):: Re1, Re1p, Re2, Re2p, Ro1, Ro1p, Ro2, Ro2p,We, Wo



sign=1
choice=1


open(unit=12,file=filename,status="replace",action="write")

! Treat separately n=0 because odd functions only exist for n>0
n=0
                Re1=  MathieuRadial(0,1,0,n,Q,u,k_max, choice, sign)
                Re1p= MathieuRadial(0,1,1,n,Q,u,k_max, choice, sign)
        
                Re2=  MathieuRadial(0,2,0,n,Q,u,k_max, choice, sign)
                Re2p= MathieuRadial(0,2,1,n,Q,u,k_max, choice, sign)

                We= Re1*Re2p-Re2*Re1p-one
                write(12,'(I5, '' & '', e13.6e2, " & ", e13.6e2, "\\")') n, real(We), aimag(We)

                We= Re1*Re2p-Re2*Re1p-one


Do n=1,N_max
        
!               EVEN FUNCTIONS

                Re1=  MathieuRadial(0,1,0,n,Q,u,k_max, choice, sign)
                Re1p= MathieuRadial(0,1,1,n,Q,u,k_max, choice, sign)
        
                Re2=  MathieuRadial(0,2,0,n,Q,u,k_max, choice, sign)
                Re2p= MathieuRadial(0,2,1,n,Q,u,k_max, choice, sign)
        
!               ODD FUNCTIONS

                Ro1=  MathieuRadial(1,1,0,n,Q,u,k_max, choice, sign)
                Ro1p= MathieuRadial(1,1,1,n,Q,u,k_max, choice, sign)
        
                Ro2=  MathieuRadial(1,2,0,n,Q,u,k_max, choice, sign)
                Ro2p= MathieuRadial(1,2,1,n,Q,u,k_max, choice, sign)


!               Test wronskian
                We= Re1*Re2p-Re2*Re1p-one
                Wo= Ro1*Ro2p-Ro2*Ro1p-one

write(12,'(I5, '' & '', e13.6e2, " & ", e13.6e2, " & ", e13.6e2, " & ", e13.6e2, " \\" )') n,real(We),aimag(We),real(Wo),aimag(Wo)

end do

close(12)

end subroutine







