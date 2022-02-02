program in_gamma_test

use , intrinsic :: iso_fortran_env , only : iostat_end

use in_gamma_sequence , only : scaled_in_gamma

implicit none

integer , parameter :: dp = selected_real_kind( 15 )
integer , parameter :: li = selected_int_kind( 18 )
integer , parameter :: dg = digits( 1.0_dp )
integer , parameter :: u  = 10

!-----

logical      :: display_all

integer      :: run_opt , entry_opt
integer      :: ios
integer      :: n
integer      :: n0 , n1
integer      :: correct_digits

integer (li) :: total_evals
integer (li) :: xi , si

real    (dp) :: x

integer      , dimension (0:dg) :: res

real    (dp) , allocatable , dimension (:) :: s , s_maple

!-----

write (*,fmt='( "1 - compute S( n , x ) for a single value of x" )')
write (*,fmt='( "2 - use maple data file for result comparison" )')
write (*,fmt='( "Enter choice: " )',advance='no')

read (*,*,iostat=ios) run_opt

if ( ios == 0 ) then

  if ( run_opt < 1 .or. run_opt > 2 ) then
    write (*,*)
    write (*,fmt='( "Error: invalid input" )')
    stop
  end if
else
  write (*,*)
  write (*,fmt='( "Error: invalid input" )')
  stop
end if

!-----

select case( run_opt )

case( 1 )

  do
    write (*,fmt='( "Enter n0:" )',advance='no')
    read (*,*,iostat=ios) n0
    if ( ios == 0 ) exit
  end do

  do
    write (*,fmt='( "Enter n1 [>=" , i3 , "] :" )',advance='no') n0
    read (*,*,iostat=ios) n1
    if ( ios == 0 .and. n1 >= n0 ) exit
  end do

  do

    write (*,fmt='( "1 - enter decimal value for x" )')
    write (*,fmt='( "2 - enter transferred integer" )')
    write (*,fmt='( "Enter choice: " )',advance='no')

    read (*,*,iostat=ios) entry_opt
    if ( ios == 0 .and. entry_opt == 1 .or. entry_opt == 2 ) exit

  end do
  
  if ( entry_opt == 1 ) then

    do
      write (*,fmt='( "Enter x [>= 0]:" )',advance='no')
      read (*,*,iostat=ios) x
      if ( ios == 0 .and. x >= 0.0_dp ) exit
    end do

  else

    do
      write (*,fmt='( "Enter transferred integer representation for x:" )',advance='no')
      read (*,*,iostat=ios) xi
      if ( ios == 0 ) exit
    end do

    x = transfer( xi , 1.0_dp )
    write (*,fmt='( "x = " , f30.20 )') x

  end if

  !-----

  call scaled_in_gamma( n0 , n1 , x , s )

  do n = n0 , n1
    write (*,fmt='( i3 , " " , e30.20 )') n , s(n)
  end do

  !-----

case( 2 )

  do

    write (*,fmt='( "1 - display summary results (recommended)" )')
    write (*,fmt='( "2 - display all results" )')
    write (*,fmt='( "Enter choice: " )',advance='no')

    read (*,*,iostat=ios) n
    if ( ios == 0 .and. n == 1 .or. n == 2 ) exit

  end do

  display_all = ( n == 2 )

  res  = 0

  !-----

  open ( unit = u , file = "in_gamma_test.dat" , action = 'read' )

  !Skip three header lines
  read (u,*)
  read (u,*)
  read (u,*)

  read (u,*) n0 , n1

  !Skip next line
  read (u,*)

  allocate( s_maple(n0:n1) )

  !-----

  write (*,*)

  total_evals = 0

  x_vals : do

    read (u,*,iostat=ios) xi

    if ( ios == iostat_end ) exit x_vals

    x = transfer( xi , 1.0_dp )

    call scaled_in_gamma( n0 , n1 , x , s )

    if ( display_all ) write (*,fmt='( " x = " , f25.20 )') x

    do n = n0 , n1

      read (u,*) si
      s_maple(n) = transfer( si , 1.0_dp )

      if ( s(n) == s_maple(n) ) then
        correct_digits = dg 
      else
        correct_digits = exponent( s_maple(n) ) - exponent( s(n) - s_maple(n) )
      end if

      res(correct_digits) = res(correct_digits) + 1

      if ( display_all ) then
        write (*,fmt='( " S(" , i3 , ") = " , 2( " " , e30.20 ) , i3 )') n , s(n) , s_maple(n) , correct_digits
      else if ( correct_digits < dg - 3 ) then
        write (*,fmt='( i3 , " digits lost: n = " , i3 , "  x = " , f20.15 , "  trans. rep = " , i24 )') &
              & dg - correct_digits , n , x , xi
      end if

    end do

    total_evals = total_evals + ( n1 - n0 + 1 )

  end do x_vals

  !-----

  write (*,fmt='( "A total of " , i8 , " evaluations were performed." )') total_evals

  do n = 1 , dg
  
    if ( res(n) > 0 ) then
      write (*,fmt='( i3 , " correct digits : " , i8 , " occurrences (" , f7.4 , "%)" )') &
              & n , res(n) , real( 100 * res(n) ) / real( total_evals )
    end if

  end do

  write (*,fmt='( "(The floating point representation uses a " , i3 , " binary digit mantissa.)" )') digits( 1.0_dp )

  !-----

end select

!-----

end program
