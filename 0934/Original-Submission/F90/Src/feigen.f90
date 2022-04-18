! Modified on 2011/12/30 by Danilo Erricolo
! Modifications consisted of
!
! introducing "implicit none" statements
! replacing declarations real*8 with real(kind=double)
! introducing "intent(in)" and "intent(out) declarations
! adding type delcarations


!* ---------------------- MODULE feigen0.f90 ------------------------ *
!* Reference: "Numerical Algorithms with C by G. Engeln-Mueller and   *
!*             F. Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
!*                                                                    *
!*                                 F90 version by J-P Moreau, Paris.  *
!* ------------------------------------------------------------------ *
MODULE FEIGEN0         !double precision version with index starting from
                       !zero and dynamic allocations in calling program


implicit none

public::devlrg

!global constants
private
! The following line is to be used with double precision
integer, parameter:: double=kind(0.0D0)
!The following line is to be used with quadruple precision
!integer, parameter:: double=kind(0.0Q0)
real(kind=double), parameter::zero=0.0_double,one=1.0_double,TWO=2.d0,XMACH_EPS = 2.22d-16
integer,parameter:: MAXIT=50


CONTAINS

!*--------------------------------------------------------------------*
!* Auxiliary subroutines for subroutine eigen()                       *
!*--------------------------------------------------------------------*

Subroutine RSWAP(a,b)

use constants
implicit none
real(kind=double):: a,b, t

 t=a; a=b; b=t
End Subroutine

!*************************************************************8

subroutine devlrg(max_M, MatEi, maxsize, EigVal)

!subroutine eigcalc computes the complex eigenvalues EigVal of the real matrix MatEi
!
!           max_M is an integer representing the order of the matrix
!                       MatEi is a real(kind=double) matrix having
!                       size (1:maxsize, 1:maxsize)
!                       maxsize is the size of the matrix MatEi
!                       EigVal is the complex(kind=double) eigenvalues vector 
!                       (1:maxsize) containing the eigenvalues of the matrix MatEi


use constants

implicit none

integer:: max_M, maxsize
real(kind=double), intent(in):: MatEi(1:maxsize, 1:maxsize)
complex(kind=double) EigVal(1:maxsize), I
real(kind=double) BMat(0:(maxsize),0:(maxsize))
real(kind=double) EigVR(0:(maxsize)), EigVI(0:(maxsize))
integer cnt(0:(maxsize))
real(kind=double) TempR, TempI
real(kind=double), dimension(0:(maxsize),0:(maxsize)) :: EiVL

integer l,r, N

I=(0, 1.0D0)

do r=0, maxsize
        do l=0, maxsize
                BMat(r,l) = 0
        end do
end do

do r=0, max_M-1
        do l=0, max_M-1
                BMat(r,l) = MatEi(r+1, l+1)
        end do
end do


  call eigen( 0, 0, maxsize, BMat, EiVL, EigVR, EigVI, cnt, N)

do r=1, max_M
        TempR = EigVR(r-1)
        TempI = EigVI(r-1)
        EigVal(r) = TempR+I*TempI
end do

do r=max_M+1, maxsize
        EigVal(r)=0
end do



end subroutine

!****************************************

Subroutine balance       &  !balance a matrix
                 (n,     &  !size of matrix         
                  mat,   &  !input matrix
                  scal,  &  !Scaling data
                  low,   &  !first relevant row index
                  high   &  !last relevant row index
                 )
implicit none
integer,intent(in):: n
real(kind=double),intent(inout):: mat(0:n,0:n)

integer,intent(out):: low,high

real(kind=double),intent(out):: scal(0:n)

!*====================================================================*
!*                                                                    *
!*  balance balances the matrix mat so that the rows with zero entries*
!*  off the diagonal are isolated and the remaining columns and rows  *
!*  are resized to have one norm close to 1.                          *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      n        integer;  ( n > 0 )                                  *
!*               Dimension of mat                                     *
!*      mat      n x n input matrix                                   *
!*                                                                    *
!*   Output parameters:                                               *
!*   ==================                                               *
!*      mat      n x n scaled matrix                                  *
!*      low      integer;                                             *
!*      high     integer;                                             *
!*               the rows 0 to low-1 and those from high to n-1       *
!*               contain isolated eigenvalues (only nonzero entry on  *
!*               the diagonal)                                        *
!*      scal     vector of size n                                     *
!*               the vector scal contains the isolated eigenvalues in *
!*               the positions 0 to low-1 and high to n-1, its other  *
!*               components contain the scaling factors for           *
!*               transforming mat.                                    *
!*                                                                    *
!*====================================================================*
integer, parameter:: basis = 2       !Floating point basis of used CPU
integer:: i,iter, j, k,m
real(kind=double):: b2, r, c, f, g, s

  b2 = basis * basis
  m = 0
  k = n - 1
  iter=1
  do while(iter==1)
    iter = 0
    do j = k, 0, -1
      r = ZERO
      do i = 0, k
        if (i.ne.j)  r = r + DABS(mat(j,i))
      end do
      if (r == ZERO) then
        scal(k) = j
        if (j.ne.k) then
          do i = 0, k 
                    call RSWAP(mat(i,j), mat(i,k))
          end do
          do i = m, n-1 
                    call RSWAP(mat(j,i), mat(k,i))
          end do
        end if
        k=k-1
        iter = 1
      end if
    end do !j loop
  end do !while iter=1

  iter=1
  do while (iter==1)
    iter = 0
    do j = m, k
      c = ZERO
      do i = m, k
        if (i.ne.j)  c = c + DABS(mat(i,j))
      end do
      if (c == ZERO) then
        scal(m) = j
        if (j.ne.m) then
          do i = 0, k 
                    call RSWAP(mat(i,j), mat(i,m))
          end do
          do i = m, n-1 
                    call RSWAP(mat(j,i), mat(m,i))
          end do
        end if
        m = m + 1
        iter = 1
      end if
    end do !j loop
  end do !while iter=1

  low = m
  high = k
  do i = m, k 
    scal(i) = ONE
  end do         

  iter=1
  do while (iter==1)
    iter = 0
    do i = m, k
      c=ZERO; r=ZERO
      do j = m, k
        if (j.ne.i) then
          c = c + DABS(mat(j,i))
          r = r + DABS(mat(i,j))
        end if
      end do
      g = r / basis
      f = ONE
      s = c + r

      do while (c < g)
        f = f * basis
        c = c * b2
      end do

      g = r * basis
      do while (c >= g)
        f = f / basis
        c = c / b2
      end do

      if ((c + r) / f < 0.95 * s) then
        g = ONE / f
        scal(i) = scal(i) * f
        iter = 1
        do j = m, n-1 
                  mat(i,j) = mat(i,j) * g
                end do
        do j = 0, k  
                  mat(j,i) = mat(j,i) * f
        end do
      end if
    end do !i loop
  end do !while iter=1
  return
End Subroutine !balance()


!*************************************************************
Subroutine balback        &  !reverse balancing ...........
                  (n,     &  !Dimension of matrix .........
                   low,   &  !first nonzero row ...........
                   high,  &  !last nonzero row ............
                   scal,  &  !Scaling data ................
                   eivec  &  !Eigenvectors ................
                  )
implicit none

integer, intent(in):: n, low, high
real(kind=double),intent(in)::  scal(0:n)
real(kind=double),intent(out):: eivec(0:n,0:n)
!*====================================================================*
!*                                                                    *
!*  balback reverses the balancing of balance for the eigenvactors.   *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ----------------                                                 *
!*      n        integer;  ( n > 0 )                                  *
!*               Dimension of mat                                     *
!*      low      integer;                                             *
!*      high     integer;   see balance                               *
!*      eivec    n x n matrix of eigenvectors, as computed in  qr2    *
!*      scal     vector of size n;                                    *
!*               Scaling data from  balance                           *
!*                                                                    *
!*   Output parameter:                                                *
!*   ----------------                                                 *
!*      eivec    n x n matrix;                                        *
!*               Non-normalized eigenvectors of the original matrix   *
!*                                                                    *
!*   Subroutine used:   RSWAP                                         *
!*   ---------------                                                  *
!*                                                                    *
!*====================================================================*
  integer:: i,j,k
  real(kind=double):: s

  do i = low, high
    s = scal(i)
    do j = 0, n-1  
          eivec(i,j) = eivec(i,j) * s
    end do
  end do

  do i = low-1, 0, -1
    k = Int(scal(i))
    if (k.ne.i) then
      do j = 0, n-1
            call RSWAP(eivec(i,j), eivec(k,j))
      end do
    end if
  end do

  do i = high + 1, n-1
    k = Int(scal(i))
    if (k.ne.i) then
      do j = 0, n-1
            call RSWAP(eivec(i,j), eivec(k,j))
      end do
    end if
  end do
  return
End Subroutine

!******************************************************************

Subroutine elmhes      &  !reduce matrix to upper Hessenberg form
                (n,    &  !Dimension of matrix .........
                 low,  &  !first nonzero row ...........
                 high, &  !last nonzero row ............
                 mat,  &  !input/output matrix .........
                 perm  &  !Permutation vector ..........
                )

implicit none

integer,intent(in)::n,low,high
integer,dimension(0:n), intent(out):: perm
real(kind=double), intent(inout):: mat(0:n,0:n)
!*====================================================================*
!*                                                                    *
!*  elmhes transforms the matrix mat to upper Hessenberg form.        *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      n        integer;  ( n > 0 )                                  *
!*               Dimension of mat                                     *
!*      low      integer;                                             *
!*      high     integer; see  balance                                *
!*      mat      n x n matrix                                         *
!*                                                                    *
!*   Output parameter:                                                *
!*   =================                                                *
!*      mat      n x n matrix;                                        *
!*               upper Hessenberg matrix; additional information on   *
!*               the transformation is stored in the lower triangle   *
!*      perm     integer vector of size n;                            *
!*               Permutation vector for elmtrans                      *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Subroutine used:   RSWAP                                         *
!*   ===============                                                  *
!*                                                                    *
!*====================================================================*
  integer:: i,j,m
  real(kind=double)::  x, y
  do m = low + 1, high-1
    i = m
    x = ZERO
    do j = m, high
      if (DABS(mat(j,m-1)) > DABS (x)) then
        x = mat(j,m-1)
        i = j
      end if
    end do

    perm(m) = i
    if (i.ne.m) then
      do j = m - 1, n-1 
            call RSWAP(mat(i,j), mat(m,j))
      end do 
      do j = 0, high 
            call RSWAP(mat(j,i), mat(j,m))
      end do
    end if

    if (x.ne.ZERO) then
      do i = m + 1, high
        y = mat(i,m-1)
        if (y.ne.ZERO) then
          y = y / x
          mat(i,m-1) = y
          do j = m, n-1 
                    mat(i,j) = mat(i,j) - y * mat(m,j)
          end do
          do j = 0, high 
                    mat(j,m) = mat(j,m) + y * mat(j,i)
          end do
        end if
      end do !i loop
    end if !x <> ZERO
  end do !m loop
  return
End Subroutine

!****************************************************************
Subroutine elmtrans      &  !copy to Hessenberg form......
                  (n,    &  !Dimension of matrix .........
                   low,  &  !first nonzero row ...........
                   high, &  !last nonzero row ............
                   mat,  &  !input matrix ................
                   perm, &  !row permutations ............
                   h     &  !Hessenberg matrix ...........
                   )

implicit none

integer, intent(in)::n,low,high
integer, intent(in),dimension(0:n) ::perm
real(kind=double),intent(in),dimension(0:n,0:n) :: mat
real(kind=double),intent(out),dimension(0:n,0:n) ::h

!*====================================================================*
!*                                                                    *
!*  elmtrans copies the Hessenberg matrix stored in mat to h.         *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      n        integer;  ( n > 0 )                                  *
!*               Dimension of  mat and eivec                          *
!*      low      integer;                                             *
!*      high     integer; see  balance                                *
!*      mat      n x n input matrix                                   *
!*      perm     Integer vector of size n;                            *
!*               Permutation data from  elmhes                        *
!*                                                                    *
!*   Output parameter:                                                *
!*   ================                                                 *
!*      h        n x n matrix;                                        *
!*               Hessenberg matrix                                    *
!*                                                                    *
!*====================================================================*
integer:: i,j,k
  do i = 0, n-1
    do k = 0, n-1 
          h(i,k) = ZERO
    end do
    h(i,i) = ONE
  end do

  do i = high - 1, low+1, -1
    j = perm(i)
    do k = i + 1, high 
          h(k,i) = mat(k,i-1)
    end do
    if (i.ne.j) then
      do k = i, high
        h(i,k) = h(j,k)
        h(j,k) = ZERO
      end do
      h(j,i) = ONE
    end if
  end do
  return
End Subroutine

!****************************************************************
Subroutine Comdiv   &       !Complex division ................
           (        &
            ar,     &       !Real part of numerator ..........
            ai,     &       !Imaginary part of numerator .....
            br,     &       !Real part of denominator ........
            bi,     &       !Imaginary part of denominator ...
            cr,     &       !Real part of quotient ...........
            ci,     &       !Imaginary part of quotient ......
            rc      &       !return code .....................
           )

implicit none

real(kind=double)::  ar,ai,br,bi,cr,ci
integer, intent(out):: rc
!*====================================================================*
!*                                                                    *
!*  Complex division  c = a / b                                       *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      ar,ai    real*8;                                              *
!*               real, imaginary parts of numerator                   *
!*      br,bi    real*8;                                              *
!*               real, imaginary parts of denominator                 *
!*                                                                    *
!*   Output parameters:                                               *
!*   ==================                                               *
!*      cr,ci    real*8;                                              *
!*               real , imaginary parts of the quotient               *
!*                                                                    *
!*   Return code value:                                               *
!*   =================                                                *
!*      = 0      ok                                                   *
!*      = 1      division by 0                                        *
!*                                                                    *
!*====================================================================*
  real(kind=double):: tmp
  if (br == ZERO.AND.bi == ZERO) then
    rc = 1
    return
  end if
  if (dabs(br) > dabs(bi)) then
    tmp = bi / br
    br  = tmp * bi + br
    cr  = (ar + tmp * ai) / br
    ci  = (ai - tmp * ar) / br
  else
    tmp = br / bi
    bi  = tmp * br + bi
    cr  = (tmp * ar + ai) / bi
    ci  = (tmp * ai - ar) / bi
  end if
 rc = 0
 return
End Subroutine !Comdiv

!*====================================================================*

FUNCTION Comabs        &  !Complex absolute value ..........
              (ar,            &  !Real part .......................
               ai             &  !Imaginary part ..................
              )
implicit none

real(kind=double):: Comabs
real(kind=double):: ar,ai
!*====================================================================*
!*                                                                    *
!*  Complex absolute value of   a                                     *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      ar,ai    real*8;                                              *
!*               Real, imaginary parts of  a                          *
!*                                                                    *
!*   Return value :                                                   *
!*   =============                                                    *
!*      Absolute value of a (real*8)                                  *
!*                                                                    *
!*   Subroutine used:    RSWAP                                        *
!*   ===============                                                  *
!*                                                                    *
!*====================================================================*
  if (ar == ZERO.and.ai == ZERO) then
    Comabs = ZERO
    return
  end if

  ar = DABS(ar)
  ai = DABS(ai)

  if (ai > ar) then                                  !Switch  ai and ar
    call RSWAP(ai, ar)
  end if

  if (ai == ZERO) then
    Comabs = ar
  else
    Comabs = ar * DSQRT(ONE + ai / ar * ai / ar)
  end if
  return
End Function

!*====================================================================*
Subroutine  hqrvec        & !compute eigenvectors ......
                  (n,     & !Dimension of matrix .......
                   low,   & !first nonzero row .........
                   high,  & !last nonzero row ..........
                   h,     & !upper Hessenberg matrix ...
                   wr,    & !Real parts of evalues .....
                   wi,    & !Imaginary parts of evalues 
                   eivec, & !Eigenvectors ..............
                   rc     & !return code ...............

                  )

implicit none

integer, intent(in):: n, low, high

real(kind=double),dimension(0:n),intent(in):: wr,wi
real(kind=double),dimension(0:n,0:n):: h

integer, intent(out):: rc
real(kind=double),dimension(0:n,0:n),intent(out)::eivec

!*====================================================================*
!*                                                                    *
!*  hqrvec computes the eigenvectors for the eigenvalues found in hqr2*
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      n        int n;  ( n > 0 )                                    *
!*               Dimension of  mat and eivec, number of eigenvalues.  *
!*      low      int low;                                             *
!*      high     int high; see  balance                               *
!*      h        n x n upper Hessenberg matrix                        *
!*      wr       vector of size n;                                    *
!*               Real parts of the n eigenvalues.                     *
!*      wi       vector of size n;                                    *
!*               Imaginary parts of the n eigenvalues.                *
!*                                                                    *
!*   Output parameter:                                                *
!*   ================                                                 *
!*      eivec    n x n matrix, whose columns are the eigenvectors     *
!*                                                                    *
!*   Return value :                                                   *
!*   =============                                                    *
!*      =  0     all ok                                               *
!*      =  1     h is the zero matrix.                                *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Subroutine in use:                                               *
!*   =================                                                *
!*                                                                    *
!*      comdiv(): complex division                                    *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Constant used:    XMACH_EPS                                      *
!*   =============                                                    *
!*                                                                    *
!*====================================================================*
  integer code, en,i,j,na,m,k,l
  real(kind=double)::  p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp

  r=ZERO; s=ZERO; z=ZERO; norm=ZERO

  do i = 0, n-1                               !find norm of h
    do j = i, n-1
      norm = norm + DABS(h(i,j))
    end do
  end do

  if (norm == ZERO) then
    rc = 1                                    !zero matrix
    return
  end if

  do en = n-1, 0, -1                          !transform back
    p = wr(en)
    q = wi(en)
    na = en - 1
    if (q == ZERO) then
      m = en
      h(en,en) = ONE
      do i = na, 0, -1
        w = h(i,i) - p
        r = h(i,en)
        do j = m, na 
                  r = r + h(i,j) * h(j,en)
        end do
        if (wi(i) < ZERO) then
          z = w
          s = r
        else
          m = i
          if (wi(i) == ZERO) then
            if (w.ne.ZERO) then 
                          temp = w 
                        else 
                          temp=XMACH_EPS * norm
            end if
            h(i,en) = -r/temp            
          else
            !Solve the linear system:
            !| w   x |  | h[i][en]   |   | -r |
            !|       |  |            | = |    |
            !| y   z |  | h[i+1][en] |   | -s |
            x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p)**2 + wi(i)**2
            h(i,en) = (x * s - z * r) / q
            t = h(i,en)
            if (DABS(x) > DABS(z)) then 
                          temp = (-r -w * t) / x 
                        else 
                          temp = (-s -y * t) / z
            end if
            h(i+1,en) = temp
          end if
        end if !wi[i] < 0
      end do !i loop
    else if (q < ZERO) then
      m = na
      if (DABS(h(en,na)) > DABS(h(na,en))) then
        h(na,na) = - (h(en,en) - p) / h(en,na)
        h(na,en) = - q / h(en,na)
      else
        call Comdiv(-h(na,en), 0.d0, h(na,na)-p, q, h(na,na), h(na,en),code)
      end if
      h(en,na) = ONE
      h(en,en) = ZERO
      do i = na - 1, 0, -1
        w = h(i,i) - p
        ra = h(i,en)
        sa = ZERO
        do j = m, na
          ra = ra + h(i,j) * h(j,na)
          sa = sa + h(i,j) * h(j,en)
        end do

        if (wi(i) < ZERO) then
          z = w
          r = ra
          s = sa
        else
          m = i
          if (wi(i) == ZERO) then
            call Comdiv(-ra, -sa, w, q, h(i,na), h(i,en),code)
          else
            !  solve complex linear system:
            !| w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa |
            !|             | |                        | = |          |
            !|   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   |
            x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p)**2 + wi(i)**2 - q*q
            vi = TWO * q * (wr(i) - p)
            if (vr == ZERO.AND.vi == ZERO) then
              vr = XMACH_EPS * norm * (DABS(w) + DABS(q) + DABS(x) + DABS(y) + DABS(z))
            end if

            call Comdiv (x * r - z * ra + q * sa, x * s - z * sa -q * ra, vr, vi, h(i,na), h(i,en),code)
            if (DABS(x) > DABS(z) + DABS(q)) then
              h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
              h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            else
              call Comdiv (-r - y * h(i,na), -s - y * h(i,en), z, q, h(i+1,na), h(i+1,en),code)
            end if
          end if !wi[i] = 0
        end if !wi[i] < 0
      end do !i loop
    end if !else if q < 0
  end do !en loop

  do i = 0, n-1                        !Eigenvectors for the evalues for
    if (i < low.or.i > high) then      !rows < low  and rows  > high
      do k = i + 1, n-1
        eivec(i,k) = h(i,k)
      end do
    end if
  end do                

  j = n-1
  do while (j>=low)
    if (j<=high) then
          m = j 
        else 
          j = high
    end if 
    if (wi(j) < ZERO) then
      l=j-1
      do i = low, high
        y=ZERO; z=ZERO
        do k = low, m
          y = y + eivec(i,k) * h(k,l)
          z = z + eivec(i,k) * h(k,j)
        end do
        eivec(i,l) = y
        eivec(i,j) = z
      end do
    else
      if (wi(j) == ZERO) then
        do i = low, high
          z = ZERO
          do k = low, m
            z = z + eivec(i,k) * h(k,j)
          end do
          eivec(i,j) = z
        end do
      end if
    end if
        j = j - 1
  end do !j loop

  rc = 0
  return
End Subroutine


!*******************************************************
!utility subroutine used only for debugging file
Subroutine WriteMat(fp,caption,A,n)
implicit none

integer fp,n,i,j
character*(*) caption
real(kind=double),dimension(0:n,0:n):: A
  write(fp,*)  caption
  do i=0, n-1
    write(fp,10) (A(i,j), j=0,n-1)
  end do
  return
10 format(10F10.6)
End Subroutine

!**********************************************************************
Subroutine hqr2        &  !compute eigenvalues .........
                (vec,  &  !switch for computing evectors
                 n,    &  !Dimension of matrix .........
                 low,  &  !first nonzero row ...........
                 high, &  !last nonzero row ............
                 h,    &  !Hessenberg matrix ...........
                 wr,   &  !Real parts of eigenvalues ...
                 wi,   &  !Imaginary parts of evalues ..
                 eivec,&  !Matrix of eigenvectors ......
                 cnt,  &  !Iteration counter ...........
                 rc    &  !return code .................
                )

implicit none
integer,intent(in):: n, low, high,vec
integer,intent(out)::rc
integer,intent(out),dimension(0:n):: cnt

real(kind=double),dimension(0:n,0:n) ::h 
real(kind=double),dimension(0:n),intent(out) ::wr,wi
real(kind=double),dimension(0:n,0:n),intent(out) ::eivec


!**********************************************************************
!*                                                                    *
!* hqr2 computes the eigenvalues and (if vec = True) the eigenvectors *
!* of an  n * n upper Hessenberg matrix.                              *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Control parameter:                                               *
!*   -----------------                                                *
!*      vec      integer;                                             *
!*       = 0     compute eigenvalues only                             *
!*       = 1     compute all eigenvalues and eigenvectors             *
!*                                                                    *
!*   Input parameters:                                                *
!*   ----------------                                                 *
!*      n        integer;  ( n > 0 )                                  *
!*               Dimension of  h and eivec,                           *
!*               length of the real parts vector  wr and of the       *
!*               imaginary parts vector  wi of the eigenvalues.       *
!*      low      integer;                                             *
!*      high     integer;  see balance                                *
!*      h        n x n matrix;                                        *
!*               upper Hessenberg matrix as output of Elmhes          *
!*               (destroyed in the process).                          *
!*                                                                    *
!*   Output parameters:                                               *
!*   -----------------                                                *
!*      eivec    n x n matrix;  (only if vec = 1)                     *
!*               Matrix, which for vec = 1 contains the               *
!*               eigenvectors as follows:                             *
!*               For real eigebvalues the corresponding column        *
!*               contains the corresponding eigenvactor, while for    *
!*               complex eigenvalues the corresponding column contains*
!*               the real part of the eigenvactor with its imaginary  *
!*               part is stored in the subsequent column of eivec.    *
!*               The eigenvactor for the complex conjugate eigenvactor*
!*               is given by the complex conjugate eigenvactor.       *
!*      wr       vector of size n;                                    *
!*               Real part of the n eigenvalues.                      *
!*      wi       vector of size n;                                    *
!*               Imaginary parts of the eigenvalues                   *
!*      cnt      Integer vector of size n;                            *
!*               vector of iterations used for each eigenvalue.       *
!*               For a complex conjugate eigenvalue pair the second   *
!*               entry is negative.                                   *
!*                                                                    *
!*   Return code:                                                     *
!*   ------------                                                     *
!*      =   0    all ok                                               *
!*      = 4xx    Iteration maximum exceeded when computing evalue xx  *
!*      =  99    zero  matrix                                         *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Subroutine in use:                                               *
!*   -----------------                                                *
!*                                                                    *
!*      hqrvec():  reverse transform for eigenvectors                 *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Constants used:     XMACH_EPS, MAXIT                             *
!*   --------------                                                   *
!*                                                                    *
!**********************************************************************
!Labels: 10, 12, 15, 20, 30
  integer debug, en, fp,i,iter,na,ll,l,j,m,k
  real(kind=double):: p, q, r, s, t, w, x, y, z
  
  debug=0; fp=3
  if (debug.ne.0) then
        open(unit=fp,file='hqr2.lst',status='unknown')
    call WriteMat(fp,'Matrix h in begin hqr2:',h,n)
    write(fp,*) ' low=',low,' high=',high
        write(fp,*) ' 1e15*MACH_EPS=', 1.d15*XMACH_EPS, ' ONE=', ONE
  end if

  p=ZERO; q=ZERO; r=ZERO 
  do i = 0, n-1
    if (i < low.or.i > high) then
      wr(i) = h(i,i)
      wi(i) = ZERO
      cnt(i) = 0
    end if
  end do

  en = high
  t = ZERO

  do while (en >= low)
    iter = 0
    na = en - 1

    do while(1<2)
      ll=999                          
      do l = en, low+1, -1                      !search for small
                                                !subdiagonal element
        if (debug.ne.0) then
          write(fp,*) 'l=',l,' a=',DABS(h(l,l-1)),' b=',XMACH_EPS*(DABS(h(l-1,l-1))+DABS(h(l,l)))
        end if
        if (DABS(h(l,l-1)) <= XMACH_EPS * (DABS(h(l-1,l-1)) + DABS(h(l,l)))) then
          ll=l;      !save current index
          goto 10    !exit l loop
        end if
      end do
10    if (ll.ne.999) then 
        l=ll 
          else 
            l=0          !restore l
      end if
      if (debug.ne.0) then
        write(fp,*) ' iter=',iter,' l=',l,' en=',en,' na=',na
        call WriteMat(fp,'Matrix h in hqr2 label 10:',h,n)
      end if

      x = h(en,en)
      if (l == en) then                         !found one evalue
        wr(en) = x + t
        h(en,en) = x + t
        wi(en) = ZERO
        cnt(en) = iter
        en = en - 1
        goto 15                                 !exit from loop while(True)
      end if

      y = h(na,na)
      w = h(en,na) * h(na,en)

      if (l == na) then                         !found two evalues
        p = (y - x) * 0.5d0
        q = p * p + w
        z = DSQRT(DABS(q))
        x = x + t
        h(en,en) = x + t
        h(na,na) = y + t
        cnt(en) = -iter
        cnt(na) = iter
        if (q >= ZERO) then                     !real eigenvalues
          if (p<ZERO) then 
                    z=p-z 
                  else 
                    z=p+z
          end if
          wr(na) = x + z
          wr(en) = x - w / z
          s = w - w / z
          wi(na) = ZERO
          wi(en) = ZERO
          x = h(en,na)
          r = DSQRT (x * x + z * z)

          if (vec.ne.0) then
            p = x / r
            q = z / r
            do j = na, n-1
              z = h(na,j)
              h(na,j) = q * z + p * h(en,j)
              h(en,j) = q * h(en,j) - p * z
            end do

            do i = 0, en
              z = h(i,na)
              h(i,na) = q * z + p * h(i,en)
              h(i,en) = q * h(i,en) - p * z
            end do

            do i = low, high
              z = eivec(i,na)
              eivec(i,na) = q * z + p * eivec(i,en)
              eivec(i,en) = q * eivec(i,en) - p * z
            end do
          end if !if vec <> ZERO
        else                                  !pair of complex
          wr(na) = x + p
          wr(en) = x + p
          wi(na) =   z
          wi(en) = - z
        end if !if q>=ZERO

        en = en - 2
        goto 15                               !exit while(1<2)
      end if !if l = na

      if (iter >= MAXIT) then
        cnt(en) = MAXIT + 1
        rc = en
        write(*,*) ' stop at iter >= MAXIT.'
        goto 20                               !MAXIT Iterations: close debugging file if open,
      end if                                  !and return

      if (iter.ne.0.and.MOD(iter,10) == 0) then
        t = t + x
        do i = low, en 
                  h(i,i) = h(i,i) - x
        end do
        s = DABS(h(en,na)) + DABS(h(na,en-2))
        x = 0.75d0 * s; y = x
        w = -0.4375d0 * s * s
      end if

      iter = iter + 1

      do m = en - 2, l, -1
        z = h(m,m)
        r = x - z
        s = y - z
        p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
        q = h(m + 1,m + 1) - z - r - s
        r = h(m + 2,m + 1)
        s = DABS(p) + DABS(q) + DABS (r)
        p = p / s
        q = q / s
        r = r / s
        if (m == l)  goto 12
        if (DABS(h(m,m-1)) * (DABS(q) + DABS(r)) <= XMACH_EPS * DABS(p)        &
                 * (DABS(h(m-1,m-1)) + DABS(z) + DABS(h(m+1,m+1)))) then
          goto 12                !exit m loop
        end if
      end do

12    do i = m + 2, en 
        h(i,i-2) = ZERO
      end do
      do i = m + 3, en 
            h(i,i-3) = ZERO
      end do

      do k = m, na
        if (k.ne.m) then         !double  QR step, for rows l to en
                                 !and columns m to en
          p = h(k,k-1)
          q = h(k+1,k-1)
          if (k.ne.na) then 
                    r = h(k+2,k-1) 
                  else 
                    r = ZERO
          end if
          x = DABS(p) + DABS(q) + DABS(r)
          if (x == ZERO) goto 30                  !next k
          p = p / x
          q = q / x
          r = r / x
        end if
        s = DSQRT(p * p + q * q + r * r)
        if (p < ZERO) s = -s

        if (k.ne.m) then
          h(k,k-1) = -s * x
        else if (l.ne.m) then
          h(k,k-1) = -h(k,k-1)
        end if
        p = p + s
        x = p / s
        y = q / s
        z = r / s
        q = q / p
        r = r / p

        do j = k, n-1                          !modify rows
          p = h(k,j) + q * h(k+1,j)
          if (k.ne.na) then
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * z
          end if
          h(k+1,j) = h(k+1,j) - p * y
          h(k,j)   = h(k,j) - p * x
        end do

        if (k+3 < en) then 
                  j=k+3 
                else 
                  j=en
        end if
        do i = 0, j                            !modify columns
          p = x * h(i,k) + y * h(i,k+1)
          if (k.ne.na) then
            p = p + z * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
          end if
          h(i,k+1) = h(i,k+1) - p * q
          h(i,k)   = h(i,k) - p
        end do

        if (vec.ne.0) then     !if eigenvectors are needed ..................
          do i = low, high
            p = x * eivec(i,k) + y * eivec(i,k+1)
            if (k.ne.na) then
              p = p + z * eivec(i,k+2)
              eivec(i,k+2) = eivec(i,k+2) - p * r
            end if
            eivec(i,k+1) = eivec(i,k+1) - p * q
            eivec(i,k)   = eivec(i,k) - p
          end do
        end if
30    end do !k loop

    end do !while(1<2)

15 end do !while en >= low                         All evalues found

  if (vec.ne.0) then                            !transform evectors back
    call hqrvec (n, low, high, h, wr, wi, eivec,rc)
  else
    rc=0
  end if

20 if (debug.ne.0)  close(fp)
  return
End Subroutine !hqr2


!**********************************************************************
Subroutine norm_1      &    !normalize eigenvectors to have one norm 1
                  (n,  &    !Dimension of matrix ...........
                   v,  &    !Matrix with eigenvectors ......
                   wi, &    !Imaginary parts of evalues ....
                   rc  &    !return code ...................
                  )
implicit none
integer,intent(in):: n
real(kind=double),dimension(0:n),intent(in)::  wi
real(kind=double),dimension(0:n,0:n),intent(inout)::  v

integer,intent(out):: rc
!**********************************************************************
!*                                                                    *
!*  norm_1 normalizes the one norm of the column vectors in v.        *
!*  (special attention to complex vectors in v  is given)             *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Input parameters:                                                *
!*   ----------------                                                 *
!*      n        integer; ( n > 0 )                                   *
!*               Dimension of matrix v                                *
!*      v        n x n matrix;                                        *
!*               Matrix of (not normalized) eigenvectors              *
!*      wi       vector of size n;                                    *
!*               Imaginary parts of the eigenvalues                   *
!*                                                                    *
!*   Output parameter:                                                *
!*   ----------------                                                 *
!*      v        n x n matrix;                                        *
!*               Matrix with normalized eigenvectors                  *
!*                                                                    *
!*   Return code:                                                     *
!*   -----------                                                      *
!*      = 0      all ok                                               *
!*      = 1      n < 1                                                *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   function or subroutine used:                                     *
!*   ---------------------------                                      *
!*      real*8   comabs():  complex absolute value                    *
!*               comdiv():  complex division                          *
!*                                                                    *
!**********************************************************************
  integer:: i,j
  real(kind=double):: maxi, tr, ti

  if (n < 1) then
    rc = 1
    return
  end if

  j = 0
  do while (j < n)
    if (wi(j) == ZERO) then
      maxi = v(0,j)
      do i = 1, n-1
        if (DABS(v(i,j)) > DABS(maxi))  maxi = v(i,j)
      end do

      if (maxi.ne.ZERO) then
        maxi = ONE / maxi
        do i = 0, n-1 
                  v(i,j) = v(i,j) * maxi
        end do
      end if
    else
      tr = v(0,j)
      ti = v(0,j+1)
      do i = 1, n-1
        if (Comabs(v(i,j), v(i,j+1)) > Comabs(tr, ti)) then
          tr = v(i,j)
          ti = v(i,j+1)
        end if
      end do

      if (tr.ne.ZERO.or.ti.ne.ZERO) then
        do i = 0, n-1
          call Comdiv(v(i,j), v(i,j+1), tr, ti, v(i,j), v(i,j+1),rc)
        end do
      end if
      j = j + 1                                        !raise j by two
    end if
        j = j + 1
  end do !j loop
  rc = 0
  return
End Subroutine

!**********************************************************************
!utility subroutine used only for debugging file
Subroutine WriteVec(fp,caption,V,n)
implicit none

integer,intent(in):: fp,n
integer,dimension(0:n),intent(in):: V
integer:: i


character*(*) caption


  write(fp,*)  caption
  write(fp,10) (V(i), i=0,n-1)
  write(fp,*) ' '
10 format(10I5)
End Subroutine


!**********************************************************************
Subroutine eigen (  &      !Compute all evalues/evectors of a matrix
           vec,     &      !switch for computing evectors ...
           ev_norm, &      !normalize Eigenvectors? .........
           n,       &      !size of matrix ..................
           mat,     &      !input matrix ....................
           eivec,   &      !Eigenvectors ....................
           valre,   &      !real parts of eigenvalues .......
           valim,   &      !imaginary parts of eigenvalues ..
           cnt,     &      !Iteration counter ...............
           rc       &      !return code .....................
          )

implicit none

integer,intent(in):: vec,ev_norm,n
real(kind=double),dimension(0:n,0:n):: mat

integer,intent(out):: rc
integer,dimension(0:n),intent(out):: cnt
real(kind=double),dimension(0:n),intent(out):: valre,valim
real(kind=double),dimension(0:n,0:n),intent(out):: eivec
!**********************************************************************
!*                                                                    *
!* The subroutine eigen  determines all eigenvalues and (if desired)  *
!* all eigenvectors of a real square  n * n  matrix via the QR method *
!* in the version of Martin, Parlett, Peters, Reinsch and Wilkinson.  *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Litterature:                                                     *
!*   -----------                                                      *
!*      1) Peters, Wilkinson: Eigenvectors of real and complex        *
!*         matrices by LR and QR triangularisations,                  *
!*         Num. Math. 16, p.184-204, (1970); [PETE70]; contribution   *
!*         II/15, p. 372 - 395 in [WILK71].                           *
!*      2) Martin, Wilkinson: Similarity reductions of a general      *
!*         matrix to Hessenberg form, Num. Math. 12, p. 349-368,(1968)*
!*         [MART 68]; contribution II,13, p. 339 - 358 in [WILK71].   *
!*      3) Parlett, Reinsch: Balancing a matrix for calculations of   *
!*         eigenvalues and eigenvectors, Num. Math. 13, p. 293-304,   *
!*         (1969); [PARL69]; contribution II/11, p.315 - 326 in       *
!*         [WILK71].                                                  *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Control parameters:                                              *
!*   ------------------                                               *
!*      vec      integer;                                             *
!*               call for eigen :                                     *
!*      = 0      compute eigenvalues only                             *
!*      = 1      compute all eigenvalues and eigenvectors             *
!*      ortho    boolean flag that shows if transformation of mat to  *
!*               Hessenberg form shall be done orthogonally by        *
!*               `orthes' (flag True) or elementarily by `elmhes'     *
!*               (flag False). The Householder matrices used in       *
!*               orthogonal transformation have the advantage of      *
!*               preserving the symmetry of input matrices.           *
!*               (only elmhes is implemented here).                   *
!*      ev_norm  integer flag that shows if Eigenvectors shall be     *
!*               normalized (flag=1) or not (flag=0).                 *
!*                                                                    *
!*   Input parameters:                                                *
!*   ----------------                                                 *
!*      n        integer; ( n > 0 )                                   *
!*               size of matrix, number of eigenvalues                *
!*      mat      n x n matrix;                                        *
!*               input matrix                                         *
!*                                                                    *
!*   Output parameters:                                               *
!*   -----------------                                                *
!*      eivec    n x n matrix;     (only if vec = 1)                  *
!*               matrix, if  vec = 1  that holds the eigenvectors     *
!*               thus :                                               *
!*               If the jth eigenvalue of the matrix is real then the *
!*               jth column is the corresponding real eigenvector;    *
!*               if the jth eigenvalue is complex then the jth column *
!*               of eivec contains the real part of the eigenvector   *
!*               while its imaginary part is in column j+1.           *
!*               (the j+1st eigenvector is the complex conjugate      *
!*               vector.)                                             *
!*      valre    vector of size n;                                    *
!*               Real parts of the eigenvalues.                       *
!*      valim    vector of size n;                                    *
!*               Imaginary parts of the eigenvalues                   *
!*      cnt      Integer vector of size n;                            *
!*               vector containing the number of iterations for each  *
!*               eigenvalue. (for a complex conjugate pair the second *
!*               entry is negative).                                  *
!*                                                                    *
!*   Return code:                                                     *
!*   ------------                                                     *
!*      =   0    all ok                                               *
!*      =   1    n < 1 or other invalid input parameter               *
!*      =   2    insufficient memory                                  *
!*      = 10x    error x from balance()                               *
!*      = 20x    error x from elmh()                                  *
!*      = 30x    error x from elmtrans()   (for vec = Trie only)      *
!*      = 4xx    error xx from hqr2()                                 *
!*      = 50x    error x from balback()    (for vec = True only)      *
!*      = 60x    error x from norm_1()     (for vec = True only)      *
!*                                                                    *
!* ------------------------------------------------------------------ *
!*                                                                    *
!*   Subroutines used:                                                *
!*   ----------------                                                 *
!*                                                                    *
!*   balance (): Balancing of an  n x n  matrix                       *
!*   elmhes ():  Transformation to upper Hessenberg form              *
!*   elmtrans(): intialize eigenvectors                               *
!*   hqr2 ():    compute eigenvalues/eigenvectors                     *
!*   balback (): Reverse balancing to obtain eigenvectors             *
!*   norm_1 ():  Normalize eigenvectors                               *
!*                                                                    *
!**********************************************************************
  integer high,i,low
  real(kind=double),dimension(0:n)::  scale 
  integer perm(0:n)
  integer fp  !debugging file number 
  integer debug

  debug = 0; fp=4
  
  if (n < 1) then                            !case n < 1 .............
    rc = 1
    return
  end if                                        

  do i = 0, n-1
    cnt(i) = 0
  end do

  if (n == 1) then                            !case n = 1 .............
    eivec(0,0) = ONE
    valre(0)   = mat(0,0)
    valim(0)   = ZERO
    rc = 0
    return
  end if
                                              !balance mat for nearly
  call Balance(n, mat, scale, low, high)      !equal row and column
                                              !one norms

  call Elmhes(n, low, high, mat, perm)        !reduce mat to upper
                                              !Hessenberg form

  if (vec.ne.0) then                          !initialize eivec
    call elmtrans (n, low, high, mat, perm, eivec)
  end if

  if (debug.ne.0) then
    open(unit=fp,file='eigen.lst',status='unknown')
    call WriteMat(fp,' Matrix A after Elmtrans:',mat,n)
    call WriteVec(fp,' Vector perm after Elmtrans:',perm,n)
    call WriteMat(fp,' Matrix Eivec after Elmtrans:',eivec,n)
  end if

  call hqr2 (vec, n, low, high, mat,        & !execute Francis QR algorithm
             valre, valim, eivec, cnt,rc)     !to obtain the eigenvalues and
                                              !the eigenvectors of input matrix
  if (rc.ne.0) then
    rc=2
    return
  end if

  if (vec.ne.0) then

    call balback(n, low, high, scale, eivec)  !reverse balancing if
                                              !eigenvaectors are to
                                              !be determined
    if (ev_norm.ne.0) then
      call norm_1 (n, eivec, valim, rc)       !normalize eigenvectors
    end if

    if (rc.ne.0) then
      rc=3
      return
    end if
  end if

  rc = 0
  if (debug.ne.0) close(fp)
End Subroutine !eigen

END MODULE
! ------------------------ END feigen0.f90 ----------------------------



