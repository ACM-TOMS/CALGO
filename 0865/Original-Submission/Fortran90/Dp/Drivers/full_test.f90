! -- A Fully Portable High Performance Minimal --
! -- Storage Hybrid Format Cholesky Algorithm  --
!    F.G. Gustavson - IBM, USA;
!    J.K. Reid - RAL, UK;
!    J. Wasniewski - DTU, Denmark

program testprogram
! .. Use statements ..
  use block_hybrid_Cholesky, only: PpHpp, HppTf, &
                   HppPp, HppTs, HppTs1, default_nb
  implicit none
  ! .. Locals ..
  integer, parameter :: wp = kind(1.0d0) ! Precision parameter
  real(wp) :: anorm ! Infinity norm of matrix
  character*1 :: UpLo
  integer :: i  ! Temporary 
  integer :: ii ! Either 1 for nb present or 2 for nb absent
  integer :: info ! Error indicator variable
  integer :: iUpLo ! Either 1 for UpLo = 'l' or 2 for UpLo = 'u'  
  integer :: j ! Interior loop counter
  integer :: k ! Position counter of the array ap
  integer :: n ! Matrix order
  integer :: mb ! RHS block size
  integer :: nb ! Matrix block size
  integer :: nrhs ! Number of right-hand sides
  integer,parameter :: sing(10)=(/1,66,89,133,262,263,264,265,266,267/)
  real(wp), allocatable :: ap(:), app(:), b(:,:), bb(:,:), bbb(:,:)
  real(wp), allocatable :: v(:) ! Holds the column vector needed
           ! for multiplication a packed matrix
  real(wp) :: sf, ss, ssd, ss1 ! Hold the largest residuals
  ! .. Intrinsic ..
  intrinsic epsilon ! The model precision epsilon
  intrinsic huge ! The model precision largest value
  intrinsic tiny ! The model precision smallest value
  intrinsic precision ! The model precision maximum number of digits
  intrinsic max ! The maximum value of the sequence
  intrinsic maxval ! The maximum value of an array

  ! .. Executable Statements ..
    write(*,"(/,a)" ) '-- test example program results --'

  ! The base machine constants show the machine precision used 
    write( *,"(/,a,/,a)" ) 'The base machine constants:', &
               '    epsilon       huge       tiny    precision'
    write(*,"(3es12.4,i6)") epsilon(1.0_wp), huge(1.0_wp), &
               tiny(1.0_wp), precision(1.0_wp)

    write(*,"(/,a)") 'Checking default_nb'
    if(default_nb<0)write(*,"(a)")'default_nb is not positive'

    write(*,"(/,a)") 'Testing for invalid subroutine arguments:'
    nb = 40; n = 267; mb=40; nrhs = 10;

    allocate( ap(n*(n+1)/2), b(n,nrhs),stat=info )
    call infoer( 'allocate for ap and b', info )
    write(*,"(a)") '  _PpHpp:'
    call PpHpp( 'a', n, ap, info, nb )
    call infoer( '_PpHpp', info, .false. )
    call PpHpp( 'l', -2, ap, info, nb )
    call infoer( '_PpHpp', info, .false. )
    call PpHpp( 'l', n, ap, info, -2 )
    call infoer( '_PpHpp', info, .false. )
    
    write(*,"(/,a)") '  _HppPp:'
    call HppPp( 'a', n, ap, info, nb )
    call infoer( '_HppPp', info, .false. )
    call HppPp( 'l', -2, ap, info, nb )
    call infoer( '_HppPp', info, .false. )
    call HppPp( 'l', n, ap, info, -2 )
    call infoer( '_HppPp', info, .false. )
    
    write(*,"(/,a)") '  _HppTf:'
    call HppTf( 'a', n, ap, info, nb )
    call infoer( '_HppTf', info, .false. )
    call HppTf( 'l', -2, ap, info, nb )
    call infoer( '_HppTf', info, .false. )
    call HppTf( 'l', n, ap, info, -2 )
    call infoer( '_HppTf', info, .false. )

    write(*,"(/,a)") '  _HppTs:'
    call HppTs('a', n, nrhs, ap, b, n, info, nb, mb )
    call infoer( '_HppTs', info, .false. )
    call HppTs('l', -5, nrhs, ap, b, n, info, nb, mb )
    call infoer( '_HppTs', info, .false. )
    call HppTs('l', n, -5, ap, b, n, info, nb, mb )
    call infoer( '_HppTs', info, .false. )
    call HppTs('l', n, nrhs, ap, b, n-1, info, nb, mb )
    call infoer( '_HppTs', info, .false. )
    call HppTs('l', n, nrhs, ap, b, n, info, 0, mb )
    call infoer( '_HppTs', info, .false. )
    call HppTs('l', n, nrhs, ap, b, n, info, nb, 0 )
    call infoer( '_HppTs', info, .false. )

    write(*,"(/,a)") '  _HppTs1:'
    call HppTs1( 'a', n, ap, b(:,1), info, nb )
    call infoer( '_HppTs1', info, .false. )
    call HppTs1( 'l', -3, ap, b(:,1), info, nb )
    call infoer( '_HppTs1', info, .false. )
    call HppTs1( 'l', n, ap, b(:,1), info, 0 )
    call infoer( '_HppTs1', info, .false. )

    write(*,"(/,a)") 'Cases that are not positive definite:'
    do j = 1,2
      if( j == 1 )UpLo = 'l'
      if( j == 2 )UpLo = 'u'
      write(*,"(3a)") 'UpLo = ', UpLo, ':'
      do k = 1,size(sing)
        call getmatr( UpLo, n )
        if( UpLo == 'l' .or. UpLo == 'L' )then
          i = n + 1 - sing(k)
          i = (n*(n+1) - i*(i+1))/2 +1
        else
          i = sing(k)
          i = i*(i+1)/2
        end if
        ap(i) = 0.0_wp
        call PpHpp( UpLo, n, ap, info )
        call HppTf( UpLo, n, ap, info )
        call infoer( '_HppTf', info, .false. )
      end do
    end do
    deallocate(ap, b, stat=i)

    write( *,"(/,2a)" ) 'Valid cases for different', &
                     ' UpLo, n, nrhs, nb, and mb:'
    do
      deallocate( ap, app, b, bb, bbb, stat = info )
      read(*,*,iostat=info) n, nrhs, nb, mb
      if( info /= 0 ) exit
      write(*,"(/,a)")'    n nrhs   nb   mb'
      write(*,"(4i5,/)")n, nrhs, nb, mb

      allocate( ap(n*(n+1)/2), app(n*(n+1)/2), b(n,nrhs), &
                bb(n,nrhs), bbb(n,nrhs), stat=info )
      call infoer('allocate: ap, app, b, bb and bbb', info, .false. )

      if( info == 0 )then

        call random_number(bb)

        do iUpLo = 1,2
        if(info /= 0)exit
          if( iUpLo == 1 )UpLo = 'l'
          if( iUpLo == 2 )UpLo = 'u'
          write( *,"(2a)" ) 'UpLo = ', UpLo

          do ii = 1,2
          if(info /= 0)exit

!         A triangular of symmetric positive definite matrix
!         is generated, saved in Lapack packed format.
          call getmatr( UpLo, n ); app = ap

!           The matrix is transformed to the packed blocked hybrid format
            if( ii == 1 )call PpHpp( UpLo, n, ap, info, nb )
            if( ii == 2 )call PpHpp( UpLo, n, ap, info )
            call infoer( '_PpHpp', info, .false. )
            if(info /= 0)exit

!           The matrix is factorized in packed blocked hybrid format
            if( ii == 1 )call HppTf( UpLo, n, ap, info, nb )
            if( ii == 2 )call HppTf( UpLo, n, ap, info )
            call infoer( '_HppTf', info, .false. )
            if(info /= 0)exit

!           The system is solved in packed blocked hybrid format
!           mb has a specified value
            b = bb
            if( ii == 1 )call HppTs(UpLo, n, nrhs, ap, b, n, info, nb, mb )
            if( ii == 2 )call HppTs(UpLo, n, nrhs, ap, b, n, info, mb=mb )
            call infoer( '_HppTs', info, .false. )
            if(info /= 0)exit
!           The solution residual is calculated
            call get_res_rhs( ss )

!           The system is solved in packed blocked hybrid format
!           mb has a default value
            b = bb
            if( ii == 1 )call HppTs(UpLo, n, nrhs, ap, b, n, info, nb )
            if( ii == 2 )call HppTs(UpLo, n, nrhs, ap, b, n, info )
            call infoer( '_HppTs', info, .false. )
            if(info /= 0)exit
!           The solution residual is calculated
            call get_res_rhs( ssd )

!           The system is solved in packed blocked hybrid format with
!           a single rhs. The solution residual is calculated.
            if( nrhs == 1 )then
              b = bb
              if( ii == 1 )call HppTs1(UpLo, n, ap, b(:,1), info, nb )
              if( ii == 2 )call HppTs1(UpLo, n, ap, b(:,1), info )
              call infoer( '_HppTs1', info, .false. )
              if(info /= 0)exit
              call get_res_rhs( ss1 )
            end if

!           The matrix L is transformed to Lapack packed format
            if( ii == 1 )call HppPp( UpLo, n, ap, info, nb )
            if( ii == 2 )call HppPp( UpLo, n, ap, info )
            call infoer( '_HppTf', info, .false. )
            if(info /= 0)exit
!           The factorized matrix residual is calculated
            call get_res_fac( UpLo, n, sf )

!           Printing the residuals:
            if( ii == 1 )write( *,"(2(a,i5),a)") 'n = ', n, &
                              '  nb = ', nb, ' (nb specified)'
            if( ii == 2 )write(*,"(2(a,i5),a)") '  n = ', n, &
                              '  nb = ', default_nb, ' (nb is a default)'
            call print_res_fac( sf )

            write(*,"(2(a,i5),a)") 'nrhs = ', nrhs, '  mb = ', &
                    mb, ' (mb specified)'
            call print_res_rhs( ss )

            write(*,"(2(a,i5),a)") 'nrhs = ', nrhs, '  mb = ', &
                    nb, ' (mb is a defaulted to nb)'
            call print_res_rhs( ssd )

            if( nrhs == 1 )then
              write(*,"(a)") 'Single rhs'
              call print_res_rhs( ss1 )
            end if

          end do ! of ii
        end do ! of UpLo
      end if
    end do

contains
  function matr(n,i,j)
    integer, intent(in) :: n, i, j
    real(wp) :: matr
    matr = (n+1-i)/100.00_wp+j
  end function matr

! A symmetric positive definite matrix is generated and saved in 
! Lapack packed format.
  subroutine getmatr( UpLo, n )
    character, intent(in) :: UpLo
    integer, intent(in) :: n
    integer :: i, j, k
    if (allocated(v)) deallocate (v)
    allocate( v(n), stat=i )
    call infoer( 'allocate; v', i )
    v(:) = 0
    k = 1
    do j = 1,n
      if( UpLo == 'u' .or. UpLo == 'U' )k = j*(j+1)/2
      do i = j,n
        ap(k) = matr(n,i,j)
        v(i) = v(i) + ap(k)
        if (i>j) v(j) = v(j) + ap(k)
        if( UpLo == 'u' .or. UpLo == 'U' )then
          k = k + i
        else
          k = k + 1
        end if
      end do
    end do
    anorm = maxval(v(:))
  end subroutine getmatr

! The information error subroutine
  subroutine infoer( name, info, terminate )
    character(len=*), intent(in) :: name
    integer, intent(in) :: info
    logical, intent(in), optional :: terminate
    logical :: lterminate
    if( present(terminate) )then
      lterminate = terminate
    else
      lterminate = .true.
    end if  
    if( info /= 0 )then
      write(*,"(3a, i3)")'Failure in ', name, ' with info = ',info
      if( lterminate )then
        stop 'stop in infoer'
      end if
    end if
  end subroutine infoer  

! The residuum matrix of abs(A-L*L**T) is calculated
  subroutine get_res_fac( UpLo, n, s )
    character, intent(in) :: UpLo
    integer, intent(in) :: n
    real(wp), intent(out) :: s
    character*1 trans
    integer :: i
    if (allocated(v)) deallocate (v)
    allocate( v(n), stat=i )
    call infoer( 'allocate; v', i )
    
    if( UpLo == 'l' .or. UpLo == 'L' )then
      trans = 'n'
    else
      trans = 't'; k = 1
    end if
    s = 0.0_wp
    do i = 1,n
      v = 0.0_wp
      if( UpLo == 'l' .or. UpLo == 'L' )then
        k = i
      end if
      do j = 1,i
        v(j) = ap(k)
        if( UpLo == 'l' .or. UpLo == 'L' )then
          k = k+n-j
        else
          k = k +1
        end if
      end do
      call dtpmv ( UpLo , trans, 'n', n, ap, v, 1 )    
      do j = 1,i-1
        v(j) = v(j) - matr(n,i,j)
      end do
      do j = i,n
        v(j) = v(j) - matr(n,j,i)
      end do
      s = max(s,maxval(abs(v)))/anorm
    end do
  end subroutine get_res_fac

! The factorisation residual is printed
  subroutine print_res_fac( s )
    real(wp), intent(in) :: s
    
! Check the residual
     if (s< 25*epsilon(1.0_wp)) then
        write(*,"(a)") &
         'The scaled factorization residual is less then 25*epsilon'
     else
        write(*,"(a,f7.1,a)") &
         'The scaled factorization residual is', s/epsilon(1.0_wp), '*epsilon'
     end if
  end subroutine print_res_fac

! The solution residual is calculated
  subroutine get_res_rhs( ss )
    real(wp), intent(out) :: ss
    integer :: i
      bbb = bb
      ss = 0.0_wp
      do i = 1,nrhs
        call dspmv ( UpLo, n, 1.0_wp, app, b(:,i), 1, -1.0_wp, bbb(:,i), 1 )
        ss = max(ss,maxval(abs(bbb(:,i))))/anorm
      end do
  end subroutine get_res_rhs

! The solution residual is printed
  subroutine print_res_rhs( ss )
    real(wp), intent(in) :: ss
        if (ss< 25*epsilon(1.0_wp)) then  !jkr
          write(*,"(a)") &
             'The scaled solution residual is less then 25*anorm*epsilon'
        else
          write(*,"(a,f9.1,a)") 'The scaled solution residual is', &
               ss/epsilon(1.0_wp), '*epsilon'
        end if
  end subroutine print_res_rhs

end program testprogram
