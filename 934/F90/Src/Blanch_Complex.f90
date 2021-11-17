Module Blanch_Complex


use constants
use FEIGEN0 
implicit none
private

public:: MathieuRadial, MathieuAngular
   

contains


!   ==========================================================            

   function MathieuAngular(parity,der,n,Q,x,k_max, choice, D_m, norm)
!   
!     copyright by Danilo Erricolo
!     University of Illinois at Chicago
!     Oct 1, 2002
!   Modified on Dec. 19, 2008
!   ==========================================================
!   Purpose
!   Compute the Mathieu angular functions Se_m, So_m
!   parity:    0 Even, Re, ..
!          1 Odd, Ro, ...
!   der:     0 No derivatives required
!          1 Derivatives are required
!   n:      order of the Mathieu function 
!   Q:      parameter of the Mathieu function, according to the notation of NcLachlan
!   x       argument of the Mathieu function
!     k_max:    index of the last element of D_m, i.e. D_m contains k_max+1 elements
!   choice:      type of Normalization 
!                    1 Stratton
!                    2 Ince
!                    3 Neutral
!     D_m      expansion coefficients
!   norm     Normalization 
!   =============================================================

    
    use constants
   IMPLICIT NONE
   integer,intent(in) :: parity, der, n, choice
     integer, intent(out):: K_max
     complex(kind=double), intent(in):: Q
     real(kind=double), intent(in):: x
   complex(kind=double), intent(out)::norm
   complex(kind=double) :: MathieuAngular
     complex(kind=double), dimension(0:MAX_), intent(out) :: D_m
   
   integer :: steps,p,case_type, comp_coeff
   real(kind=double):: Q0,A0,iDeltaQ, X1, PI2,PI
     complex(kind=double):: Q1,QT, A, S1, S2, DS1, DS2, Snorm
     
   pi=dacos(-one)

     if (parity==0) then
        if (n==2*int(n/2)) then
            case_type=0
        else  
            case_type=1
        end if
     else
        if (n==2*int(n/2)) then
            case_type=2
        else
            case_type=3
        end if
    end if
    
    if (real(q)>= zero) then
        Q0=real(Q)
        call init_eigenval_approx(parity, n, Q0, A0)

        steps=int(abs(aimag(Q))/0.5D0)+1
        iDeltaQ=aimag(Q)/steps

        A=A0
        comp_coeff=0
        do p=1, steps-1
            QT=Q0+cmplx(zero,p*iDeltaQ,double)
          call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
        end do
        comp_coeff=1
        call Mathieu_Eigenval_Coeff(A,Q,n,case_type,D_m,k_max,norm,choice,comp_coeff)
        call Compute_Angular_functions(case_type,X,D_m,K_max,S1, DS1)  
    
      if (der.eq.0) then
            MathieuAngular=S1
        else
            MathieuAngular=DS1
    end if
    else
        Pi2=PI/2.0D0
        X1=PI2-x        
        Q1=-Q
    Q0=real(Q1)           
    
        if (parity==0) then
            if (N==2*int(n/2)) then
        
                call init_eigenval_approx(parity, n, Q0, A0)

                steps=int(abs(aimag(Q1))/0.5D0)+1
                iDeltaQ=aimag(Q1)/steps

                A=A0
                comp_coeff=0
                do p=1, steps-1
                    QT=Q0+cmplx(zero,p*iDeltaQ,double)
                    call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
                end do

                comp_coeff=1
                call Mathieu_Eigenval_Coeff(A,Q1,n,case_type,D_m,k_max,norm,choice,comp_coeff)

                call Compute_Angular_functions(case_type,X1,D_m,K_max,S1, DS1)
                
                call Compute_Angular_functions(case_type,PI2,D_m,K_max,Snorm, DS2)

            else
                call init_eigenval_approx(1, n, Q0, A0)

        steps=int(abs(aimag(Q1))/0.5D0)+1
                iDeltaQ=aimag(Q1)/steps

                case_type=3
                A=A0
                comp_coeff=0
                do p=1, steps-1
                    QT=Q0+cmplx(zero,p*iDeltaQ,double)
                    call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
                end do

                comp_coeff=1
                call Mathieu_Eigenval_Coeff(A,Q1,n,case_type,D_m,k_max,norm,choice,comp_coeff)

                call Compute_Angular_functions(case_type,X1,D_m,K_max,S1, DS1)
                    
                call Compute_Angular_functions(case_type,PI2,D_m,K_max,Snorm,DS2)
            end if
                    
            if (der==0) then
                MathieuAngular=S1/Snorm
            else
                MathieuAngular=-DS1/Snorm
            end if
        else
            if (N==2*int(n/2)) then
                call init_eigenval_approx(parity, n, Q0, A0)

        steps=int(abs(aimag(Q1))/0.5D0)+1
                iDeltaQ=aimag(Q1)/steps

        A=A0
                comp_coeff=0
                do p=1, steps-1
                    QT=Q0+cmplx(zero,p*iDeltaQ,double)
                    call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
                end do
                
                comp_coeff=1
                call Mathieu_Eigenval_Coeff(A,Q1,n,case_type,D_m,k_max,norm,choice,comp_coeff)

                call Compute_Angular_functions(case_type,X1,D_m,K_max,S1, DS1)

        call Compute_Angular_functions(case_type,PI2,D_m,K_max,S2, Snorm)
            else
                call init_eigenval_approx(0, n, Q0, A0)

        steps=int(abs(aimag(Q1))/0.5D0)+1
                iDeltaQ=aimag(Q1)/steps

                case_type=1
                A=A0
                comp_coeff=0
                do p=1, steps-1
                    QT=Q0+cmplx(zero,p*iDeltaQ,double)
                    call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
                end do

                comp_coeff=1
                call Mathieu_Eigenval_Coeff(A,Q1,n,case_type,D_m,k_max,norm,choice,comp_coeff)

                call Compute_Angular_functions(case_type,X1,D_m,K_max,S1, DS1)
                call Compute_Angular_functions(case_type,PI2,D_m,K_max,S2, Snorm)
            end if

            if (der==0) then
                MathieuAngular=-S1/Snorm
            else
                MathieuAngular= DS1/Snorm
            end if
        end if
    end if
    
     
end function MathieuAngular
!***************************************************************************







!   ==========================================================            

   function MathieuRadial(parity,radial_kind,der,n,Q,x,k_max, choice, sign)

!   
!     copyright by Danilo Erricolo
!     University of Illinois at Chicago
!     Oct 1, 2002
!   Modified on Feb. 7, 2009

!   ==========================================================
!   Purpose
!   Compute the Mathieu radial functions Re1,Re2,Re3,Re4,Ro1,Ro2,Ro3,Ro4
!   parity:      0 Even, Re, ..
!                    1 Odd, Ro, ...
!   radial_kind  1,2,3,4 
!   der            0 No derivatives required
!                    1 Derivatives are required
!   n             order of the Mathieu function 
!   Q             parameter of the Mathieu function, according to the notation of Blanch
!   x             argument of the Mathieu function
!     k_max         index of the last element of D_m, i.e. D_m contains k_max+1 elements
!   choice:      type of Normalization 
!                    1 Stratton
!                    2 Ince
!                    3 Neutral
!     sign         1 or -1 depending on the sign of sqrt(Q) that appears in the argument of the radial functions
!   =============================================================

    
     use constants
   IMPLICIT NONE

   integer,intent(in) :: parity, radial_kind, der, n, choice, sign
     integer, intent(out):: K_max
     complex(kind=double), intent(in):: Q
     real(kind=double), intent(in):: x
   complex(kind=double) :: MathieuRadial
   
   integer :: KC, steps,p,case_type, comp_coeff, case_type_coeff, parity_coeff
   real(kind=double):: Q0,A0,iDeltaQ
     complex(kind=double):: QT, Q1, A, R1,R1p,R2,R2p, norm
     complex(kind=double), dimension(0:MAX_) :: D_m
    

    if (parity==0) then
        if (n==2*int(n/2)) then
            case_type=0
        else  
            case_type=1
        end if
    else
        if (n==2*int(n/2)) then
            case_type=2
        else
            case_type=3
        end if
    end if
    
    KC=3

    if (real(q)>0) then
        Q1=Q
        Q0=real(Q1)

        call init_eigenval_approx(parity, n, Q0, A0)

    steps=int(abs(aimag(Q1))/0.5D0)+1
        iDeltaQ=aimag(Q1)/steps
        A=A0
        comp_coeff=0
        do p=1, steps-1
            QT=Q0+cmplx(zero,p*iDeltaQ,double)
          call Mathieu_Eigenval_Coeff(A,QT,n,case_type,D_m,k_max,norm,choice,comp_coeff)
        end do

        comp_coeff=1
        call Mathieu_Eigenval_Coeff(A,Q1,n,case_type,D_m,k_max,norm,choice,comp_coeff)
    call Compute_Radial_functions(KC,case_type,N,Q1,X,D_M,K_max,R1,R1p,R2,R2p,sign)     
    else
        Q1=-Q
        q0=real(Q1)
            select case (case_type)
                case(0)     
                            case_type_coeff=0
                            parity_coeff=0
                case(1)     
                            case_type_coeff=3    
                            parity_coeff=1
                case(2)     
                            case_type_coeff=2
                            parity_coeff=1 
                case(3)     
                            case_type_coeff=1
                            parity_coeff=0 
            end select
        
        steps=int(abs(aimag(Q1))/0.5D0)+1
        iDeltaQ=aimag(Q1)/steps

        call init_eigenval_approx(parity_coeff, n, Q0, A0)
        
        A=A0
        comp_coeff=0
        do p=1, steps-1
            QT=Q0+cmplx(zero,p*iDeltaQ,double)
          call Mathieu_Eigenval_Coeff(A,QT,n,case_type_coeff,D_m,k_max,norm,choice,comp_coeff)
        end do
        
        comp_coeff=1
        call Mathieu_Eigenval_Coeff(A,Q1,n,case_type_coeff,D_m,k_max,norm,choice,comp_coeff)
        call Comp_Rad_func_neg_par_cmplx(KC,case_type, N,Q1,X,D_m,K_max,R1,R1p,R2,R2p)
    end if
    
      if (der.eq.0) then
    select case (radial_kind)
      case (1)
      MathieuRadial=R1
      case (2)
      MathieuRadial=R2
      case (3)
      MathieuRadial=R1+j*R2
      case (4)
      MathieuRadial=R1-j*R2
     end select
        else
    select case (radial_kind)
      case (1)
      
          MathieuRadial=R1p
      case (2)
      MathieuRadial=R2p
      case (3)
      MathieuRadial=R1p+j*R2p
      case (4)
      MathieuRadial=R1p-j*R2p
    end select
    
    end if
end function MathieuRadial
!***************************************************************************




!*******************************************************************
subroutine sort_eigenvalues_real(Eig, length)
! Returns eigenvalues ordered in increasing order

use constants
implicit none

integer, intent(in):: length
real(kind=double), dimension(max_size),intent(inout)::EIG


integer,dimension(max_size):: ind
real(kind=double), dimension(max_size)::EIG_temp

integer:: k

call qsortd(Eig,ind, length)

eig_temp=zero
do k=1,length
    eig_temp(k)=eig(ind(K))
end do

eig=eig_temp

end subroutine sort_eigenvalues_real
!***************************************************************************************************



SUBROUTINE qsortd(x,ind,n)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18 Time: 11:55:47

use constants
IMPLICIT NONE


INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), dimension(max_size), INTENT(IN) :: x
INTEGER, dimension(max_size),INTENT(OUT)  :: ind
INTEGER, INTENT(IN)  :: n

!***************************************************************************

!                             ROBERT RENKA
!                         OAK RIDGE NATL. LAB.

!  THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
! ARRAY X INTO INCREASING ORDER. THE ALGORITHM IS AS FOLLOWS. IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND. X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T. THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER. INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION. THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION. WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -  N - LENGTH OF THE ARRAY X.

!           X - VECTOR OF LENGTH N TO BE SORTED.

!          IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!             FASHION AS X WOULD BE. THUS, THE ORDERING ON
!             X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

INTEGER  :: iu(21), il(21)
INTEGER  :: m, i, j_, k, l, ij, it, itt, indx
REAL   :: r
REAL (dp) :: t

! LOCAL PARAMETERS -

! IU,IL = TEMPORARY STORAGE FOR THE UPPER AND LOWER
!      INDICES OF PORTIONS OF THE ARRAY X
! M =   INDEX FOR IU AND IL
! I,J =  LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =  INDICES IN THE RANGE I,...,J
! IJ =   RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =  TEMPORARY INDEX FOR X
! R =   PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =   CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO i = 1, n
 ind(i) = i
END DO
m = 1
i = 1
j_ = n
r = .375

! TOP OF LOOP

20 IF (i >= j_) GO TO 70
IF (r <= .5898437) THEN
 r = r + .0390625
ELSE
 r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + r*(j_-i)
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!  INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
 ind(ij) = indx
 ind(i) = it
 it = indx
 t = x(it)
END IF

! INITIALIZE L

l = j_

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!  INTERCHANGE IT WITH T

indx = ind(j_)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j_) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!  INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!  NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!  ARRAY YET TO BE SORTED

IF (l-i > j_-k) THEN
 il(m) = i
 iu(m) = l
 i = k
 m = m + 1
 GO TO 80
END IF

il(m) = k
iu(m) = j_
j_ = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j_ = iu(m)

80 IF (j_-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J. NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j_) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90
END SUBROUTINE qsortd








!=========================================
function cm(m)
!=========================================
!     copyright by Danilo Erricolo
!     University of Illinois at Chicago
!     July 19, 2002

implicit none

integer m,cm
 
 if (m==0) then
  cm=2
 else
  cm=1
 end if
end function cm
!=========================================


!=========================================
function V(m,a,q)
!=========================================
use constants
implicit none

complex(kind=double):: V
complex(kind=double),INTENT(IN):: A,Q
integer, intent(IN):: m

V=(a-m*m)/Q

end function V
!=========================================


!===========================================
function tail(M2S,alpha,q)
!===========================================

use constants
implicit none
complex(kind=double):: tail
integer, intent(IN) ::M2S
complex(kind=double),intent(IN):: alpha,q

integer:: k,mf
complex(kind=double):: A,Ap,App,B,Bp,Bpp,bk,C,Q_k,Q_p

    App= one
    Ap= zero
    Bpp= zero
    Bp= one

    k=1
    bk=V(M2S,alpha,Q)

    A=bk*Ap+App
    B=bk*Bp+Bpp
    
    Q_k=A/B
    Q_p=zero
    
    do while (abs(Q_k-Q_p)>Tol_Tail)
        App=Ap
        Ap=A
        Bpp=Bp
        Bp=B
    
        k=k+1
        mf=M2S+2*(k-1)
        bk=(alpha-mf**2)/Q
 
        A=bk*Ap-App
        B=bk*Bp-Bpp

        Q_p=Q_k
        Q_k=A/B
    end do


    C=zero
    k=mf-2
    do 
        C=1/(V(k,alpha,Q)-C)
        if (k==M2S) exit
        k=k-2
    end do

  tail=0.5D0*(Q_k+C)

end function tail



!===========================================
function tail_2(M2S,alpha,q)
!===========================================

use constants
implicit none
complex(kind=double):: tail_2
integer, intent(IN) ::M2S
complex(kind=double),intent(IN):: alpha,q

integer:: k,mf
complex(kind=double):: A,Ap,App,B,Bp,Bpp,b_k,b_p,C,W_p,W_k,d_k,d_p,nu_k,S_k,S_p

    k=1
    b_k=(alpha-(M2S)**2)/Q
    
    App= one
    Ap= zero
    Bpp= zero
    Bp= one

    A=b_k*Ap+App
    B=b_k*Bp+Bpp

    W_p=A/B
    W_k=W_p

    d_p=W_p

    do while (abs(W_k-W_p)>Tol_Tail)
        k=k+1
        mf=M2S+2*(k-1)

        b_p=b_k
        b_k=(alpha-mf**2)/Q
        
        nu_k=-1/(b_k*b_p)

        d_k=-nu_k*(1.0D0+d_p)/(1.0D0+nu_k*(1.0D0+d_p))

        S_k=d_k*S_p
        W_k=W_p+S_k

        S_p=S_k
        W_p=W_k

    end do

    C=zero
    k=mf-2
    do 
        C=1/(V(k,alpha,Q)-C)
        if (k==M2S) exit
        k=k-2
    end do

  tail_2=0.5D0*(W_k+C)

end function tail_2



!=========================================
function Hm2(m,a,q,G_2,m_start)
!
!     copyright by Danilo Erricolo
!     University of Illinois at Chicago
!     July 19, 2002
!=========================================
use constants
implicit none


integer,intent(in):: m,m_start
complex(kind=double),intent(in):: A,Q
complex(kind=double),dimension(1:Max_), intent(in):: G_2

integer:: k
complex(kind=double):: Hm2,Vm


k=((m+2)-m_start)/2+1
Vm=(a-m*m)/Q

Hm2=(Vm-G_2(k))/cm(m-2)

end function Hm2
!=========================================

!=========================================
function Fm(m,a,q,H_m,m_start)
!=========================================
use constants
implicit none


integer,INTENT(in):: m,m_start
complex(kind=double),INTENT(in):: A,q
complex(kind=double),DIMENSION(1:Max_),intent(IN)::H_M

integer :: k
complex(kind=double):: Fm,Vm

k=((m+2)-m_start)/2+1
if (m>=m_start+2) then
    Vm=(a-m*m)/Q
    Fm=Vm*H_m(k)-one
ELSE
    Fm=(Vm*H_m(k)-one)/2.0D0
end if

end function
!=========================================





!************************************************************************************************
subroutine normalization(case_type,choice,A_m,m1,m_max,B_m,norm)
!     copyright by Danilo Erricolo
!     University of Illinois at Chicago
!     July 19, 2002

! choice==1 Stratton
! choice==2 Ince
! choice==3 Neutral
! norm returns a normalization coefficient only for Stratton's normalization

use constants
implicit none



integer, intent(in):: case_type,m1,m_max,choice
COMPLEX(kind=double), dimension(0:MAX_),intent(in) :: A_m
COMPLEX(kind=double), dimension(0:MAX_),intent(out):: B_m
complex(kind=double),intent(out):: norm

integer:: k,k1,k2,k3,m_start,p
REAL(KIND=DOUBLE):: PI
COMPLEX(kind=double):: A,sum,sum1,sum2,sum_sq,sum_sq1,sum_sq2

pi=acos(-one)

select case(case_type)
    case(0) 
        m_start=2
    case(1)
        m_start=3
    case(2)
        m_start=4
        p=2
    case(3)
        m_start=3
        p=1
end select


if (choice==1) then
    sum1=zero
    sum2=zero
    sum_sq1=zero
    sum_sq2=zero
    k1=((m1-2)-m_start)/2+1
    k2=(m1-m_start)/2+1
    k3=(m_max-m_start)/2+1

    if (case_type<=1) then
        do k=0,k1
            sum1=sum1+A_m(k)
        end do
        do k=k3,k2,-1
            sum2=sum2+A_m(k)
        end do
    
    else
        do k=0,k1
            sum1=sum1+(2*k+p)*A_m(k)
        end do
        do k=k3,k2,-1
            sum2=sum2+(2*k+p)*A_m(k)
        end do
    end if
    
    do k=0,k1        
        sum_sq1=sum_sq1+A_m(k)**2
    end do
    do k=k3,k2,-1
        sum_sq2=sum_sq2+A_m(k)**2
  end do

    sum=sum1+sum2
    do k=0,k3
        B_m(k)=A_m(k)/sum
    end do

    if (case_type==0) then
        Norm=pi*(sum_sq1+sum_sq2+A_m(0)**2)/sum**2
    else
        Norm=pi*(sum_sq1+sum_sq2)/sum**2
    end if
else if (choice==2) then
    sum_sq1=zero
    k1=((m1-2)-m_start)/2+1
    do k=0,k1
        sum_sq1=sum_sq1+A_m(k)**2
    end do

    sum_sq2=zero
  k3=(m_max-m_start)/2+1
    k2=(m1-m_start)/2+1
    do k=k3,k2,-1
        sum_sq2=sum_sq2+A_m(k)**2
    end do
    
    if (case_type==0) then
        sum_sq=sqrt(sum_sq1+sum_sq2+A_m(0)**2)*A_m(0)/abs(A_m(0))
    else
        sum_sq=sqrt(sum_sq1+sum_sq2)*A_m(0)/abs(A_m(0))
    end if

    do k=0,k3
        B_m(k)=A_m(k)/sum_sq
    end do
    Norm=one
else
  A=A_m(0)
    k1=(m_max-m_start)/2+1
    do k=1,k1
        if (abs(A_m(k))>abs(A)) then
         A=A_m(k)
        end if
    end do

    do k=0,k1
        B_m(k)=A_m(k)/A
    end do
    Norm=-one
end if
end subroutine normalization
!************************************************************************************************




!************************************************************************************************

subroutine init_eigenval_approx(parity, order, Q, A)
! parity : 0 for even and 1 for odd
! order : the integer order of the Mathieu function
! Q   : the complex parameter of Mathieu equation
! A   : the approximate value of the eigenvalue of Mathieu function for the parameter q


use constants
implicit none
integer, intent(in) :: parity, order
real(kind=double), intent(in):: Q
real(kind=double), intent(out)::A

integer:: position, size, length 
real(kind=double), dimension(max_size) ::EIG

size=10+2*order

call Mathieu_eigenvalue_real_par(parity, size, q,EIG,length) 
call sort_eigenvalues_real(Eig, length)         

        if (parity==0) then
      position=order+1  
        else
            position=order  
        end if
        A=eig(position)

end subroutine




!=======================================================
subroutine Mathieu_Eigenval_Coeff(A,Q,order,case_type,D_m,k_max,norm,choice,comp_coeff)
!=======================================================
! copyright by Danilo Erricolo
! University of Illinois at Chicago
! Oct 9, 2008
!
! INPUT
! A: eigenvalue of the Mathieu equation that corresponds to a solution of order "m" and type given by "kind"
! Q: parameter 
! order: order of the Mathieu function 
! case_type: assumes four possible values
!    0 - even solution of period pi     (Type 0)
!        1 - even solution with period 2*pi (Type 1)
!        2 - odd solution of period pi     (Type 2)
!    3 - odd solution of period 2*pi    (Type 3)
! 
! choice:parameter that selects the type of normalization according to
!        1 - Morse-Stratton
!    2 - Ince
!    3 - Neutral
! comp_coeff: 0 does not compute the expansion coefficients, but only the Mathieu eigenvalue
!             1     computes the expansion coefficients and the Mathieu eigenvalue          

! OUTPUT
! D_m:  array of expansion coefficients    
! k_max: Index of last element of D_m
! norm: normalization coefficient, only provided when "choice==2"

! This subroutine uses the integer parameter "exp_coeff_factor" defined in the module constants.f90
! Its purpose is to force the computation of a number of expansion coeffiencients that is at least exp_coeff_factor*order/2
! When this subroutine is run at double precision, a good number for exp_coeff_factor is 4. In case of quadruple precision,
! a good number is 7.


use constants
implicit none

complex(kind=double), intent(inout) :: A
complex(kind=double), intent(in) :: Q
integer, intent(in):: order,case_type,choice, comp_coeff
complex(kind=double), dimension(0:MAX_),intent(out) :: D_m
integer, intent(out):: k_max 
complex(kind=double),intent(out):: norm


integer :: k,m,m_max,m_add,m_start,M2S,M1,counter,k1
real(kind=double):: F_m,pi,U,FL
complex(kind=double):: P0, P1, P2, Temp, A1, A2
complex(kind=double),dimension(1:MAX_) :: G_1,G_1p,G_1pp,G_2,G_2p,G_2pp,H_m,Tm,A_array
COMPLEX(kind=double), dimension(0:MAX_):: A_m
logical:: loop, additional_terms


pi=acos(-one)

G_2=cmplx(zero,zero,double)
G_1=cmplx(zero,zero,double )

FL=2.0D0**126


if (abs(Q)>zero) then

counter=0
    do
    counter=counter+1
    A_array(counter)=A

    select case(case_type)
        case(0) 
        m_start=2
        G_1(1)=V(0,A,Q)
        case(1)
        m_start=3
        G_1(1)=V(1,A,Q)-one
    case(2)
        m_start=4
        G_1(1)=V(2,A,Q)
    case(3)
        m_start=3
        G_1(1)=V(1,A,Q)+one
    end select
    m=m_start
    U=0.5*sqrt(3.0D0*Q+A) 
    M2s=2*floor(U)+m+4

  k=1
    do while ((abs(G_1(k))>=one) .and. (m<M2S))
        H_m(k)=1/G_1(k)         
        M=M+2
        k=k+1
        G_1(k)=V(M-2,A,Q)-cm(m-4)*H_m(k-1)
    END DO

    if (abs(G_1(k))>=one) THEN
        H_m(k)=1/G_1(k)     
    END IF

    M1=M 

    IF ((M==M2S) .and. (abs(G_1(k))>one)) then
        G_2(k)=tail(M2S,A,Q)
        Tm(K)=G_2(k)-G_1(k) 
                
        do while (m-m_start>0) 
            m=m-2
            k=k-1
            G_2(k)=cm(m-2)/(V(m,a,Q)-G_2(k+1)) 
            Tm(k)=G_2(k)-G_1(k)
            if (abs(Tm(k))>abs(Tm(k+1))) then    
                m1=m+2
                exit
            else
                m1=m
            end if
        end do     
    end if 

    if (m1<M2S) then !2.41
        k=(m2s-m_start)/2+1
        G_2(k)=tail(M2S,a,q) 
        m=M2S-2
        k=k-1
        do while (m>=m1) 
            H_m(k)=Hm2(m,a,q,G_2,m_start) 
            if (abs(H_m(k))>=one) then
                G_2(k)=1/H_m(k)         
                H_m(k)=FL
            else
                loop=.true.
                do while (loop)
                    if (m==m1) then 
                        if (abs(H_m(k))==zero) then   
                            print *,"Fatal error"
                        else                      
                            G_2(k)=1/H_m(k)                     
                            H_m(k)=FL
                        end if
                        loop=.false.
                    else                                  
                        G_2(k)=FL
                        m=m-2
                        k=k-1
                        F_m=Fm(M,A,Q,H_m,m_start)    
                        if (abs(F_m)>=abs(H_m(k+1))) then 
                            G_2(k)=H_m(k+1)/F_m      
                            loop=.false.
                            H_m(k)=FL
                        else
                            H_m(k)=F_m/H_m(k+1)       
                            G_2(k)=FL
                        end if
                    end if
                end do
            end if
            m=m-2      
            k=k-1
        end do
  else 
        m=m2s
        k=(m2s-m_start)/2+1
        Tm(k)=G_2(k)-G_1(k)    
    end if


k1=(m1-m_start)/2+1
P0=G_2(k1)-G_1(k1) 

if (abs(P0)<1.0D-6) then
        exit
    else
        G_1p(1)=1/q
        do m=m_start+2, m2s, 2
            k=(m-m_start)/2+1
            G_1p(k)=1/q+cm(m-4)*H_m(k-1)**2*G_1p(k-1)            
        end do

    m=m2s
        k=(m-m_start)/2+1
        G_2p(k)=zero
        do while (m>m_start)
            m=m-2
            k=k-1
            G_2p(k)= -G_2(k)**2*(1/q-G_2p(k+1))/cm(m-2)
        end do

        P1=G_2p(k1) -G_1p(k1) 
        if (abs(P1)>0.1) then
            a=a-P0/P1

      if (abs(A_array(counter)-a)<1.0D-6) then
             exit
            end if
   
        else
            G_1pp(1)=zero

            do m=m_start+2, m2s, 2
                k=(m-m_start)/2+1
                G_1pp(k)=cm(m-4)*(-2.0D0*H_m(k-1)**3*G_1p(k-1)**2+H_m(k-1)**2*G_1pp(k-1))
            end do

            m=m2s
            k=(m-m_start)/2+1            
            G_2pp(k)=zero
            do while (m>m_start)
                m=m-2
                k=k-1
                G_2pp(k)=-(2.0D0*G_2(k)*G_2p(k)*(1/q-G_2p(k+1))-G_2(k)**2*G_2pp(k+1))/cm(m-2)
            end do
    
            P2=G_2pp(k1)-G_1pp(k1) 
            TEMP=SQRT((P1/P2)**2-2.0D0*P0/P2)
            A1=-P1/P2+TEMP
            A2=-P1/P2-TEMP
            if (abs(A1)>abs(A2)) then
                a=a+a2 
            else
                a=a+a1
            end if

        end if
    end if

    end do

    if (comp_coeff==1) then
        m=m1-2
        k=(m-m_start)/2+1
        A_m(k)=one
        k=k-1
        do while (k>=0)
            A_m(k)=A_m(k+1)/G_1(k+1)
            k=k-1
        end do
    
        m=m1
        k=(m-m_start)/2+1
        do while (m<=M2S)
            if (ABS(G_2(k))>=FL) then
                A_m(k)=A_m(k-2)/Fm(M,A,Q,H_m,m_start)
            else
                A_m(k)=G_2(k)*A_m(k-1)
            end if 
            m=m+2 
            k=k+1
        end do
    
      k_max=k-1

        m=m2s
        k=(m-m_start)/2+1
        m_add=ceiling(log(abs(A_m(k)/TOL))/log(1/abs(G_2(k))))+1

      additional_terms=.false.
        if (m_add>0) then
            m_max=m2s+m_add+2
            additional_terms=.true.
        end if

        if (m_max<exp_coeff_factor*order) then
            m_max=exp_coeff_factor*order
            additional_terms=.true.
        end if

        if (additional_terms) then

            m=m_max
            k=(m-m_start)/2+1
            G_2(k)=tail(m,A,Q)
            m=m-2
            k=k-1
            do while (m>M2S)    
                G_2(k)=cm(m-2)/(V(m,A,Q)-G_2(k+1))
                m=m-2
                k=k-1
            end do

            do while (m<=m_max)
                A_m(k)=G_2(k)*A_m(k-1)
                m=m+2
                k=k+1
            end do

            k_max=k-1
        end if


        call normalization(case_type,choice,A_m,m1,m_max,D_m,norm)

    end if 

else
    A=order**2

    if (comp_coeff==1) then
        k=(order-m_start)/2+1
        m_max=order


        if (choice==1) then
            if (case_type<=1) then
                D_m(k)=one
            else
                D_m(k)=one/order
            end if
            if ((case_type==1) .and. (order==0)) then
            Norm=2*pi*D_m(0)**2
            else
                Norm=pi*D_m(k)**2
            end if
        else if (choice==2) then
            if ((case_type==0) .and. (order==0)) then
                D_m(k)=1/sqrt(2.0D0)
            else
                D_m(k)=one
            end if
        else
            D_m(k)=one
        end if
    
        k_max=k         
    end if

end if


end subroutine
!******************************************************************************************





!******************************************************************************************
SUBROUTINE Compute_Angular_functions(case_type,X,D_m,K_max,S, DS)
!
!     by Danilo Erricolo
!     University of Illinois at Chicago
!     Dec. 24, 2008
!    ==============================================================
!    Purpose: Compute Mathieu Angular functions and their derivatives
!
!    Input:  
!                case_type --- 0,1,2,3
!        x --- Argument of Mathieu functions
!                D_m --- Array of expansion coefficients
!                k_max -- index of last element of D_m

!    Output: S --- Se_m(q,x) or So_m(q,x)
!        DS --- d/dx(Se_m(q,x)) or d/dx(So_m(q,x))
!
!    ==============================================================

  use constants
    IMPLICIT NONE
    
    integer:: K
    
    integer, intent(in):: case_type, K_max
    real(kind=double), intent(in):: X
    complex(kind=double), intent(out):: S, DS
    
    complex(kind=double),dimension(0:Max_)::D_m
    
  S=cmplx(zero,zero,double)


    
IF (case_type==0) THEN
  DO K=K_Max,0,-1
       S=S+D_m(K)*COS((2*K)*X)
  END DO
ELSE IF (case_type==1) THEN
    DO K=K_Max,0,-1
       S=S+D_m(K)*COS((2*K+1)*X)
    END DO
ELSE IF (case_type==2) THEN
  DO K=K_Max,0,-1
    S=S+D_m(K)*SIN((2*K+2)*X)
    END DO
ELSE 
    DO K=K_Max,0,-1
        S=S+D_m(K)*SIN((2*K+1)*X)
    END DO
ENDIF

DS=zero

IF (case_type==0) THEN
    DO K=K_Max,0,-1
      DS=DS-(2*K)*D_m(K)*SIN((2*K)*X)
    END DO
ELSE IF (case_type==1) THEN
  DO K=K_Max,0,-1
      DS=DS-(2*K+1)*D_m(K)*SIN((2*K+1)*X)
  END DO
ELSE IF (case_type==2) THEN
  DO K=K_Max,0,-1
    DS=DS+(2.0D0*K+2)*D_m(K)*COS((2*K+2)*X)
  END DO
ELSE 
  DO K=K_Max,0,-1
      DS=DS+(2*K+1)*D_m(K)*COS((2*K+1)*X)
  END DO

ENDIF
    
END subroutine Compute_Angular_functions
!===============================================================




!******************************************************************************************
SUBROUTINE Compute_Radial_functions(KC,case_type, N,Q,X,D_m,KM,R1, R1p, R2, R2p,sign)
!
!     by Danilo Erricolo
!     University of Illinois at Chicago
!     Oct 1, 2002
!
!   this is a modification with permission from the authors of the subroutine
!   MTU12
!     originally developed and copyrighted by 
!     Shanjie Zhang and Jianming Jin and described in
!   COMPUTATION OF SPECIAL FUNCTIONS
!     JOHN WILEY & SONS, 1996
!     ISBN 0-471-11963-6 
!   
!   This version was downloaded from the URL
!     http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!   where it is written that "..All the programs and subroutines contained
!   in this archive are copyrighted. However, we give permission to the user
!   who downloads these routines to incorporate any of these routines into his
!   or her programs provided that the copyright is acknowledged. "
!
!
!     The modifications consist of:
!     1) introduced IMPLICIT NONE 
!   2) explicitly declared all variables
!     3) introduced function KIND to allow for change of precision
!     4) Removed GOTO statements
!   5) Replaced call to the subroutine FCOEF that computes the Mathieu coefficients 
!    with a call to the subroutine Blanch_Coefficients and made
!    related changes
!
!    ==============================================================
!    Purpose: Compute modified Mathieu functions of the first and
!        second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
!        and their derivatives
!    Input:  
!        KC --- Function Code
!            KC=1 for computing the first kind
!            KC=2 for computing the second kind
!              or Msm(2)(x,q) and Msm(2)'(x,q)
!            KC=3 for computing both the first
!              and second kinds
!                case_type --- 0,1,2,3
!        n --- Order of Mathieu functions
!        q --- Parameter of Mathieu functions ( q ò 0 ), according to the notation of McLachlan
!        x --- Argument of Mathieu functions
!                D_m --- Array of expansion coefficients
!                km   --- index of last element of D_m
!                sign -- 1 or -1 depending on the sign of sqrt(S)

!    Output: R1 --- Mcm(1)(x,q) or Msm(1)(x,q)
!        R1p --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
!        R2 --- Mcm(2)(x,q) or Msm(2)(x,q)
!        R2p --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
!                
!
!    ==============================================================

  use constants
    IMPLICIT NONE
    
    integer:: R,K
    
    integer, intent(in)::KC, case_type, N, KM, sign
    
    real(kind=double), intent(in):: X
    complex(kind=double), intent(in) :: Q
    complex(kind=double),intent(out):: R1, R1p, R2, R2p
    complex(kind=double)::U1,U2,S 

    real(kind=double):: PI
    complex(kind=double),dimension(0:Max_+2)::BJ1,BJ2,BY1,DJ1,DJ2,DY1
    complex(kind=double),dimension(0:Max_)::D_m

    integer:: NZ, IERR
    real(kind=double):: ZR1,ZR2,ZI1,ZI2
    real(kind=double), dimension(0:Max_+2)::CWRKR, CWRKI,BY1R,BY1I,BY2R,BY2I,BJ1R,BJ1I,BJ2R,BJ2I
    

    R=INT(N/2)
    S=4.0D0*Q

    PI=Dacos(-1.0D0)

    if ( (real(s)<0) .and. (aimag(s)<0) ) then
        U1=-sign*SQRT(S)*DEXP( X)/2.0D0 
        U2=-sign*SQRT(S)*DEXP(-X)/2.0D0
    else
        U1=sign*SQRT(S)*DEXP( X)/2.0D0 
        U2=sign*SQRT(S)*DEXP(-X)/2.0D0
    end if

    ZR1=real(U1)
    ZI1=AIMAG(U1)
    ZR2=real(U2)
    ZI2=AIMAG(U2) 
    BJ1R=zero
    BJ1I=zero
    BJ2R=zero
    BJ2I=zero
    BY1R=zero
    BY1I=zero
    BY2R=zero
    BY2I=zero
   
    call ZBESJ(ZR1, ZI1, zero, 1, KM+3, BJ1R, BJ1I, NZ, IERR)
  BJ1=cmplx(BJ1R,BJ1I,double)
  call ZBESJ(ZR2, ZI2, zero, 1, KM+3, BJ2R, BJ2I, NZ, IERR)
  BJ2=cmplx(BJ2R,BJ2I,double)
  
    call ZBESY(ZR1, ZI1, zero, 1, KM+3, BY1R, BY1I, NZ, CWRKR, CWRKI, IERR)
    BY1=cmplx(BY1R,BY1I,double)

    CWRKR=cmplx(zero,zero,double)
    CWRKI=cmplx(zero,zero,double)
    call ZBESY(ZR2, ZI2, zero, 1, KM+3, BY2R, BY2I, NZ, CWRKR, CWRKI, IERR)

    DO k=0, KM+2
        DJ1(K)=-BJ1(K+1)+K/U1*BJ1(K)
        DJ2(K)=-BJ2(K+1)+K/U2*BJ2(K)
        DY1(K)=-BY1(K+1)+K/U1*BY1(K)
    END DO

  IF ( (KC==1) .OR. (KC==3)) THEN
        R1=zero
        R1p=zero
        IF (case_type==0) THEN
            DO K=KM,0,-1
                R1=R1+(-1)**K*D_m(K)*BJ1(K)*BJ2(K)
                R1p=R1p+(-1)**K*D_m(K)*(U1*DJ1(K)*BJ2(K)-U2*BJ1(K)*DJ2(K))
            END DO
        ELSE IF (case_type==1) THEN
            DO K=KM,0,-1
                R1=R1+(-1)**K*D_m(K)*(BJ1(K+1)*BJ2(K)+BJ1(K)*BJ2(K+1)) !(3.06)
                R1p=R1p+(-1)**K*D_m(K)*(U1*DJ1(K+1)*BJ2(K)-U2*BJ1(K+1)*DJ2(K)+U1*DJ1(K)*BJ2(K+1)-U2*BJ1(K)*DJ2(K+1))
            END DO
        ELSE IF (case_type==2) THEN
            DO K=KM,0,-1
                R1=R1+ (-1)**(K+1)*D_m(K)*(BJ1(K+2)*BJ2(K)-BJ2(K+2)*BJ1(K))
                R1p=R1p+(-1)**(K+1)*D_m(K)*(U1*DJ1(K+2)*BJ2(K)-U2*BJ1(K+2)*DJ2(K)+U2*DJ2(K+2)*BJ1(K)-U1*BJ2(K+2)*DJ1(K))
            END DO
        ELSE
            DO K=KM,0,-1
                R1=R1+(-1)**K*D_m(K)*(BJ1(K+1)*BJ2(K)-BJ1(K)*BJ2(K+1))
                R1p=R1p+(-1)**K*D_M(K)*(U1*DJ1(K+1)*BJ2(K)-U2*BJ1(K+1)*DJ2(K)-U1*DJ1(K)*BJ2(K+1)+U2*BJ1(K)*DJ2(K+1))
            END DO
        END IF
        R1= (-1)**R*dsqrt(pi/2.0D0)*R1 /D_m(0)
        R1p=(-1)**R*dsqrt(pi/2.0D0)*R1p/D_m(0)
    END IF     

    IF ( (KC==2) .OR. (KC==3)) THEN
        R2=zero
        R2p=zero

        IF (case_type==0) THEN
            DO K=KM,0,-1
                R2= R2+(-1)**K*D_m(K)*BY1(K)*BJ2(K)
                R2p=R2p+(-1)**K*D_m(K)*(U1*DY1(K)*BJ2(K)-U2*BY1(K)*DJ2(K))
            END DO
        ELSE IF (case_type==1) THEN
            DO K=KM,0,-1
                R2= R2+(-1)**K*D_m(K)*(BY1(K+1)*BJ2(K)+BY1(K)*BJ2(K+1))
                R2p=R2p+(-1)**K*D_m(K)*(U1*DY1(K+1)*BJ2(K)-U2*BY1(K+1)*DJ2(K)+U1*DY1(K)*BJ2(K+1)-U2*BY1(K)*DJ2(K+1))
            END DO
        ELSE if (case_type==2) then
            DO K=KM,0,-1
                R2= R2+(-1)**(K+1)*D_m(K)*(BY1(K+2)*BJ2(K)-BY1(K)*BJ2(K+2))
                R2p=R2p+(-1)**(K+1)*D_m(K)*(U1*DY1(K+2)*BJ2(K)-U2*BY1(K+2)*DJ2(K)-U1*DY1(K)*BJ2(K+2)+U2*BY1(K)*DJ2(K+2))
            END DO
        ELSE
            DO K=KM,0,-1
                R2= R2+(-1)**K*D_m(K)*(BY1(K+1)*BJ2(K)-BY1(K)*BJ2(K+1))
                R2p=R2p+(-1)**K*D_m(K)*(U1*DY1(K+1)*BJ2(K)-U2*BY1(K+1)*DJ2(K)-U1*DY1(K)*BJ2(K+1)+U2*BY1(K)*DJ2(K+1))
            END DO
        END IF
        R2= (-1)**R*dsqrt(pi/2.0D0)*R2 /D_m(0)
        R2p=(-1)**R*dsqrt(pi/2.0D0)*R2p/D_m(0)
    END IF
END subroutine Compute_Radial_functions
!===============================================================





SUBROUTINE Comp_Rad_func_neg_par_cmplx(KC,case_type, N,Q,X,D_m,k_max,R1,R1p,R2,R2p)
!
!     by Danilo Erricolo
!     University of Illinois at Chicago
!     Dec. 20, 2008
!
!   A small portion of this code was obtained with permission from the authors of the subroutine
!   MTU12
!     originally developed and copyrighted by 
!     Shanjie Zhang and Jianming Jin and described in
!   COMPUTATION OF SPECIAL FUNCTIONS
!     JOHN WILEY & SONS, 1996
!     ISBN 0-471-11963-6 
!   
!   This version was downloaded from the URL
!     http://iris-lee3.ece.uiuc.edu/~jjin/jin_home.html
!   where it is written that "..All the programs and subroutines contained
!   in this archive are copyrighted. However, we give permission to the user
!   who downloads these routines to incorporate any of these routines into his
!   or her programs provided that the copyright is acknowledged. "
!
!    Input:  
!        KC --- Function Code
!            KC=1 for computing radial functions of the first kind (Re_n^{(1)}(-Q,x) or Ro_n^{(1)}(-Q,x))
!            KC=2 for computing radial functions of the second kind (Re_n^{(2)}(-Q,x) or Ro_n^{(2)}(-Q,x))
!            KC=3 for computing radial functions of both the first and second kinds
!                case_type: 0,1,2 or 3
!        n --- Order of Mathieu functions
!        q --- Parameter of Mathieu functions according to the notation of McLachlan.
!                        Because this subroutine computes radial functions for Q<0, we assume that when it is invoked, 
!            one passes only the absolute value of Q.
!        x --- Argument of Mathieu functions
!                D_m ---Array of expansion coefficients
!                k_max index of the last element of D_m, hence D_m contains K_max+1 elements

!    Output: R1 --- Value of the Mathieu radial function of the first kind: Re_n^{(1)}(-Q,x) or Ro_n^{(1)}(-Q,x)
!        R1p --- Value of the derivative of the Mathieu radial function of the first kind: d/dx(Re_n^{(1)}(-Q,x))
!            or d/dx(Ro_n^{(1)}(-Q,x))
!        R2 --- Value of the Mathieu radial function of the second kind: Re_n^{(2)}(-Q,x) or Ro_n^{(2)}(-Q,x)
!        R2p --- Derivative of the Mathieu radial function of the second kind: d/dx(Re_n^{(2)}(-Q,x)) 
!            or d/dx(Ro_n^{(2)}(-Q,x))
!                
!
!    ==============================================================

  use constants
    IMPLICIT NONE
    
    integer:: K
    
    integer, intent(in)::KC, case_type, N, K_max
    real(kind=double), intent(in):: X
    complex(kind=double), intent(in):: Q
    complex(kind=double), intent(out):: R1,R1p,R2,R2p

    integer:: R, NZ, IERR
    real(kind=double),dimension(0:Max_+3)::IKR1,IKI1,IKR2,IKI2,KJR1,KJI1,KJR2,KJI2
    complex(kind=double),dimension(0:Max_+3)::IK1, IK2, KJ1,DIK1,DIK2,DKJ1

    complex(kind=double),dimension(0:Max_)::D_m
    real(kind=double):: pi,ZR1,ZI1,ZR2,ZI2
    complex(kind=double):: s,u1,u2

    S=4.0D0*Q
        
    pi=dacos(-1.0D0)

    R=INT(N/2)

    U1=SQRT(S)*exp(X)/2.0D0
    U2=SQRT(S)*exp(-x)/2.0D0

    ZR1= real(U1)
    ZI1=AIMAG(U1)
    ZR2= real(U2)
    ZI2=AIMAG(U2)
    IKR1=zero
    IKI1=zero
    IKR2=zero
    IKI2=zero
    CALL ZBESI(ZR1, ZI1, zero, 1, K_Max+3, IKR1, IKI1, NZ, IERR)
    CALL ZBESI(ZR2, ZI2, zero, 1, K_Max+3, IKR2, IKI2, NZ, IERR)
    IK1=cmplx(IKR1,IKI1,double)
    IK2=cmplx(IKR2,IKI2,double)

    KJR1=zero
    KJI1=zero
    KJR2=zero
    KJI2=zero
    CALL ZBESK(ZR1, ZI1, zero, 1, K_Max+3, KJR1, KJI1, NZ, IERR)
    CALL ZBESK(ZR2, ZI2, zero, 1, K_Max+3, KJR2, KJI2, NZ, IERR)
    KJ1=cmplx(KJR1,KJI1,double)

  DO k=0, K_Max+2
        DIK1(K)= IK1(K+1)+K/U1*IK1(K)
        DIK2(K)= IK2(K+1)+K/U2*IK2(K)
        DKJ1(K)=-KJ1(K+1)+K/U1*KJ1(K)
    END DO

  R1=cmplx(zero,zero,double) 
    R1p=cmplx(zero,zero,double)

IF (case_type==0) THEN
            DO K=K_Max,0, -1
                    R1= R1 +(-1)**K*D_m(K)*IK1(K)*IK2(K) !(7.11)
                    R1p=R1p+(-1)**k*D_m(K)*(u1*DIK1(K)*IK2(K)-u2*IK1(K)*DIK2(K)) 
            END DO
            R1= (-1)**R*dsqrt(pi/2.0D0)* R1/D_m(0) 
            R1p=(-1)**R*dsqrt(pi/2.0D0)*R1p/D_m(0)
            

ELSE IF (case_type==1) THEN

            DO K=K_Max,0,-1
                    R1= R1+(-1)**k*D_m(K)*( IK1(K+1)*IK2(K) + IK2(K+1)*IK1(K)) 
                    R1p=R1p+(-1)**k*D_m(K)*( u1*DIK1(K+1)*IK2(K) - u2*IK1(K+1)*DIK2(K) - u2*DIK2(K+1)*IK1(K)+ U1*IK2(K+1)*DIK1(K))
            END DO
            R1= j*(-1)**R*dsqrt(pi/2.0D0)* R1/D_m(0)    
            R1p=j*(-1)**R*dsqrt(pi/2.0D0)*R1p/D_m(0)

            
ELSE IF (case_type==2) THEN
            DO K=K_Max,0,-1

                    R1= R1+(-1)**(k+1)*D_m(K)*(IK1(K)*IK2(K+2)-IK1(K+2)*IK2(K)) 
                    R1p=R1p+(-1)**(K+1)*D_m(K)*(U1*DIK1(K)*IK2(K+2)-u2*IK1(K)*DIK2(K+2)-U1*DIK1(K+2)*IK2(K)+U2*IK1(K+2)*DIK2(K))
            END DO         
            R1= (-1)**R*dsqrt(pi/2.0D0)* R1/D_m(0) 
            R1p=(-1)**R*dsqrt(pi/2.0D0)*R1p/D_m(0)     
ELSE !case_type=3
            DO K=K_Max,0,-1
                    R1= R1+(-1)**K*D_m(K)*(IK1(K+1)*IK2(K)-IK1(K)*IK2(K+1)) 
                    R1p=R1p+(-1)**K*D_m(K)*(u1*DIK1(K+1)*IK2(K)-U2*IK1(K+1)*DIK2(K)-U1*DIK1(K)*IK2(K+1)+U2*IK1(K)*DIK2(K+1))
            END DO
            R1= j*(-1)**R*dsqrt(pi/2.0D0)* R1/D_m(0)    
            R1p=j*(-1)**R*dsqrt(pi/2.0D0)*R1p/D_m(0)        
            
END IF



IF ( (KC==2) .OR. (KC==3)) THEN
        R2 =cmplx(zero,zero,double)
        R2p=cmplx(zero, zero,double)
        
        IF (case_type==0) THEN
            DO K=K_max,0,-1
                    R2=R2+D_m(K)*IK2(K)*KJ1(K) 
                    R2p=R2p+D_m(K)*(-U2*DIK2(K)*KJ1(K)+U1*IK2(K)*DKJ1(K))
            END DO
            R2=R2/D_m(0)                         
            R2=j*R1-(-1)**R*dsqrt(2.0D0/pi)*R2    
            R2p=R2p/D_m(0)                        
            R2p=j*R1p-(-1)**R*dsqrt(2.0D0/pi)*R2p  
                        
        ELSE IF (case_type==1) THEN

            DO K=K_Max,0,-1
                    R2=R2+ D_m(K)*(IK2(K)*KJ1(K+1)-IK2(K+1)*KJ1(K))
                    R2p=R2p+D_m(K)*( -u2*DIK2(K)*KJ1(K+1)+U1*IK2(K)*DKJ1(K+1)+U2*DIK2(K+1)*KJ1(K)-U1*IK2(K+1)*DKJ1(K))
            END DO
            R2= R2/D_m(0) 
            R2p=R2p/D_m(0)  
            R2= j*R1 +j*(-1)**R*dsqrt(2.0D0/pi)*R2 
      R2p=j*R1p+j*(-1)**R*dsqrt(2.0D0/pi)*R2p 
        ELSE IF (case_type==2) THEN
            DO K=K_Max,0,-1
                    R2=R2+ D_m(K)*(IK2(K)*KJ1(K+2)-IK2(K+2)*KJ1(K))
                    R2p=R2p+D_m(K)*(-U2*DIK2(K)*KJ1(K+2)+U1*IK2(K)*DKJ1(K+2)+u2*DIK2(K+2)*KJ1(K)-U1*IK2(K+2)*DKJ1(K))
            END DO
            R2= R2/D_m(0)      
            R2p=R2p/D_m(0)   
            R2= j*R1 -(-1)**R*dsqrt(2.0D0/pi)*R2 
      R2p=j*R1p-(-1)**R*dsqrt(2.0D0/pi)*R2p 
        ELSE 
            DO K=K_Max,0,-1

                    R2=R2+ D_m(K)*( IK2(K+1)*KJ1(K)+IK2(K)*KJ1(K+1))        
                    R2p=R2p+D_m(K)*( -U2*DIK2(K+1)*KJ1(K)+U1*IK2(K+1)*DKJ1(K)-u2*DIK2(K)*KJ1(K+1)+U1*IK2(K)*DKJ1(K+1))
            END DO
            R2= R2/D_m(0)      
            R2p=R2p/D_m(0)    
            R2= j*R1+ j*(-1)**R*DSQRT(2.0D0/pi)*R2 
            R2p=j*R1p+j*(-1)**R*DSQRT(2.0D0/pi)*R2p 
        ENDIF      
        


END IF
END subroutine Comp_Rad_func_neg_par_cmplx
!===============================================================


!*******************************************************************
subroutine Mathieu_eigenvalue_real_par(parity, m_max, q,EIG,length)
!*******************************************************************
!
! Author: Danilo Erricolo
!     University of Illinois at Chicago
!         September 25, 2008

use constants
implicit none

integer, intent(in):: parity, m_max
real(kind=double), intent(in):: Q
real(kind=double), dimension(max_size), intent(out)::EIG
integer, intent(out) :: length

integer :: length_A,length_B
integer:: m,k,k_max,k_max_A,k_max_B
complex(kind=double), dimension(max_size)::EIG_A,EIG_B


real(kind=double), dimension(max_size,max_size):: A,B


A=zero
B=zero

if (parity.eq.0) then
    if (m_max.eq.2*int(m_max/2)) then
        k_max_A=m_max/2
        k_max_B=k_max_A-1
  else
      k_max_A=(m_max-1)/2
        k_max_B=k_max_A
    end if

    k_max=k_max_A
  A(1,2)=Q

    A(2,1)=2*Q
    A(2,2)=4
    A(2,3)=Q
            
    DO K=2,k_max-1 
        m=(2*k)
        A(k+1,K)=Q
        A(k+1,k+1)=m**2
        A(k+1,k+2)=Q              
    END DO

    A(k_max+1,k_max)=Q
  m=2*k_max
    A(k_max+1,k_max+1)=m**2
        
    k_max=k_max_B
  B(1,1)=Q+1.0D0
    B(1,2)=Q
    
    B(2,1)=Q
    B(2,2)=9.0D0
    B(2,3)=Q

    DO k=2,k_max-1
        m=2*k+1
        B(K+1,K)=Q
        B(K+1,K+1)=m**2
        B(K+1,K+2)=Q              
    END DO
    
        B(k_max+1,k_max)=Q
    m=2*k_max+1 
      B(k_max+1,k_max+1)=m**2
else
    
    if (m_max.eq.2*int(m_max/2)) then
        k_max_A=(m_max-2)/2
        k_max_B=k_max_A
  else
      k_max_B=(m_max-1)/2
        k_max_A=k_max_B-1
    end if 
    
    k_max=k_max_A
            
    A(1,1)=4.0D0
    A(1,2)=Q

    A(2,1)=Q
    A(2,2)=16.0D0
    A(2,3)=Q
            
    DO k=2,k_max-1
        m=2*K+2
        A(k+1,k)=Q
        A(k+1,K+1)=m**2
        A(k+1,k+2)=Q              
    END DO
    A(k_max+1,k_max)=Q
  m=2*k_max+2
    A(k_max+1,k_max+1)=m**2
    
    k_max=k_max_B
    B(1,1)=1.0D0-Q
    B(1,2)=Q
    
    B(2,1)=Q
    B(2,2)=9.0D0
    B(2,3)=Q

    DO k=2,k_max-1 
        m=2*k+1
        B(k+1,k)=Q
        B(k+1,k+1)=m**2
        B(k+1,k+2)=Q              
    END DO
    
    B(k_max+1,k_max)=Q
    m=2*k_max+1 
    B(k_max+1,k_max+1)=m**2
end if

length_A=k_max_A+1 
length_B=k_max_B+1

call DEVLRG(length_A,A,max_size,EIG_A)
call DEVLRG(length_B,B,max_size,EIG_B)

do k=1,length_A
    eig(k)=real(eig_A(k)) 
end do
do k=1,length_B
    eig(length_a+k)=real(eig_B(k))
end do
length=length_A+length_B


end subroutine
!*******************************************************************



end module
