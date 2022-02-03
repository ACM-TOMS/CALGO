program DCONGEN
!
!  Generator for the BLAS/LAPACK constants module
!  E. Anderson
!  May 5, 2016
!
!  Change the following 4 lines and the format statement with the
!  label 88 for a different data type.  Also modify the setting of
!  wp in the two subroutines below.
!
   integer, parameter :: wp = 8
   character*(*), parameter :: modnam = "LA_CONSTANTS"
   character*1, parameter :: sprefix = 'D'
   character*1, parameter :: cprefix = 'Z'
!
   real(wp), parameter :: half = 0.5_wp, one = 1.0_wp
   character*2 :: day, month
   character*4 :: year
   character*8 :: dfield
   integer :: nmax
   real(wp) :: eps, ulp, safmin, safmax, smlnum, bignum, rtmin, rtmax 
   real(wp) :: b1,b2,s1m,s2m,overfl,rbig,relerr
!
!  Header
!
   call DATE_AND_TIME( dfield )
   year = dfield( 1:4 )
   month = dfield( 5:6 )
   day = dfield( 7:8 )
   write(*,92) 'module ', modnam
   write(*,91) '!'
   write(*,91) '!  -- BLAS/LAPACK module --'
   if( month == "01" ) then
      write(*,99) 'January ', day, year
   else if( month == "02" ) then
      write(*,99) 'February ', day, year
   else if( month == "03" ) then
      write(*,99) 'March ', day, year
   else if( month == "04" ) then
      write(*,99) 'April ', day, year
   else if( month == "05" ) then
      write(*,99) 'May ', day, year
   else if( month == "06" ) then
      write(*,99) 'June ', day, year
   else if( month == "07" ) then
      write(*,99) 'July ', day, year
   else if( month == "08" ) then
      write(*,99) 'August ', day, year
   else if( month == "09" ) then
      write(*,99) 'September ', day, year
   else if( month == "10" ) then
      write(*,99) 'October ', day, year
   else if( month == "11" ) then
      write(*,99) 'November ', day, year
   else if( month == "12" ) then
      write(*,99) 'December ', day, year
   else
      write(*,98) month, day, year
   end if
!
!  Standard constants
!
   write(*,91) '!'
   write(*,91) '!  Standard constants'
   write(*,91) '!'
   write(*,97) 'integer', 'wp = ', wp
   write(*,93) 'real(wp)', 'zero = 0.0_wp'
   write(*,93) 'real(wp)', 'half = 0.5_wp'
   write(*,93) 'real(wp)', 'one = 1.0_wp'
   write(*,93) 'real(wp)', 'two = 2.0_wp'
   write(*,93) 'real(wp)', 'three = 3.0_wp'
   write(*,93) 'real(wp)', 'four = 4.0_wp'
   write(*,93) 'real(wp)', 'eight = 8.0_wp'
   write(*,93) 'real(wp)', 'ten = 10.0_wp'
   write(*,93) 'complex(wp)', 'czero = ( 0.0_wp, 0.0_wp )'
   write(*,93) 'complex(wp)', 'chalf = ( 0.5_wp, 0.0_wp )'
   write(*,93) 'complex(wp)', 'cone = ( 1.0_wp, 0.0_wp )'
   write(*,94) 'character*1', 'sprefix = ', sprefix
   write(*,94) 'character*1', 'cprefix = ', cprefix
!
!  Model parameters
!
   write(*,91) '!'
   write(*,91) '!  Model parameters'
   write(*,91) '!'
   eps = half*epsilon( one )
   ulp = epsilon( one )
   safmin = tiny(one)
   safmax = one / safmin
   smlnum = safmin / ulp
   bignum = one / smlnum
   rtmin = sqrt( smlnum )
   rtmax = one / rtmin
   write(*,88) 'eps', eps
   write(*,88) 'ulp', ulp
   write(*,88) 'safmin', safmin
   write(*,88) 'safmax', safmax
   write(*,88) 'smlnum', smlnum
   write(*,88) 'bignum', bignum
   write(*,88) 'rtmin', rtmin
   write(*,88) 'rtmax', rtmax
!
!  Blue's scaling constants
!
   call x2init(nmax,b1,b2,s1m,s2m,overfl,rbig,relerr)
   write(*,91) '!'
   write(*,91) '!  Blue''s scaling constants'
   write(*,91) '!'
   write(*,88) 'tsml', b1
   write(*,88) 'tbig', b2
   write(*,88) 'ssml', s1m
   write(*,88) 'sbig', s2m
   write(*,92) 'end module ', modnam
99 format('!     ',a,a2,', ',a4)
98 format('!     ',a2,'/',a2,'/',a4)
97 format('   ',a,', parameter :: ',a,i4)
94 format('   ',a,', parameter :: ',a,'''',a,'''')
93 format('   ',a,', parameter :: ',a)
92 format(a,a)
91 format(a)
88 format('   real(wp), parameter :: ',a,' = ',e28.20e3,'_wp')
end program
subroutine x2init(nmax,b1,b2,s1m,s2m,overfl,rbig,relerr)
   integer, parameter :: wp = 8
   integer :: nmax,nbig,ibeta,it,iemin,iemax,iexp
   real(wp) :: bexp,abig,b1,b2,s1m,s2m,eps,relerr,overfl,rbig
   integer, parameter :: ione = 1
   real(wp), parameter :: one = 1.0_wp

!  x2init calculates the machine-dependent constants
!     b1, b2, s1m, s2m, relerr, overfl, nmax
!  from the integer and real data model parameters
!     nbig, ibeta, it, iemin, iemax, rbig.
!
!  Source: James L. Blue, "A Portable Fortran Program to Find
!  the Euclidean Norm of a Vector", ACM TOMS, vol. 4, no. 1,
!  March 1978.
!
!  2016-04-29:  Adjusted scaling factor for lower range to shift
!  denormals into normal range (E. Anderson)

   nbig   = huge(ione)         ! largest integer
   ibeta  = radix(one)         ! base for floating-point numbers
   it     = digits(one)        ! number of base-beta digits in mantissa
   iemin  = minexponent(one)   ! minimum exponent
   iemax  = maxexponent(one)   ! maximum exponent
   rbig   = huge(one)          ! largest floating-point number

   iexp   = -((1-iemin)/2)     ! or (iemin-1)/2
   b1     = bexp(ibeta,iexp)   ! lower boundary of midrange
   iexp   = (iemax+1-it)/2
   b2     = bexp(ibeta,iexp)   ! upper boundary of midrange

   iexp   = (1+it-iemin)/2
   s1m    = bexp(ibeta,iexp)   ! scaling factor for lower range
   iexp   = -((iemax+it)/2)
   s2m    = bexp(ibeta,iexp)   ! scaling factor for upper range

   overfl = rbig*s2m           ! overflow boundary for abig
   eps    = bexp(ibeta,1-it)
   relerr = sqrt(eps)          ! tolerance for neglecting asml
   abig   = one/eps - one
   if (float(nbig)>abig) then
      nmax = abig              ! largest safe n
   else
      nmax = nbig
   end if
   return
end subroutine x2init
function bexp(ibeta,iexp)
   integer, parameter :: wp = 8
   integer :: ibeta,iexp
   real(wp) :: bexp

!  bexp=ibeta**iexp by binary expansion of iexp,
!  exact if ibeta is the machine base
!
!  Source: James L. Blue, "A Portable Fortran Program to Find
!  the Euclidean Norm of a Vector", ACM TOMS, vol. 4, no. 1,
!  May 5, 2016

   integer :: n
   real(wp) :: tbeta

   tbeta = float(ibeta)
   bexp = 1.0
   n = iexp
   if (n<0) then
      n = -n
      tbeta = 1.0/tbeta
   end if
   do while( n >= 0 )
      if( mod(n,2) .NE. 0) bexp = bexp*tbeta
      n = n / 2
      if( n==0 ) exit
      tbeta = tbeta*tbeta
   end do
   return
end function bexp
