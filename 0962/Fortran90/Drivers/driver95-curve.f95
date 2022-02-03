! Copyright (c) 2013, Paul Muir, Jack Pew
! Paul Muir, Mathematics and Computing Science, Saint Mary's University.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! This file contains the BACOLI source code (including the new code
! that implements the interpolation based spatial error estimates.)

!-----------------------------------------------------------------------
! This is an example driver for running BACOLI through the F95 wrapper.
! See the comments within bacoli95.f95 for details regarding use.

! This driver hardcodes values for npde, xa, xb, nderiv (the number of 
! PDEs, the endpoints of the spatial domain, the number of derivatives 
! of the solution to be output) and problem parameters, and uses 
! default values for many other parameters through the wrapper.
!
! Output is printed as columns to the standard output unit and
! optionally, interpolant data is written to the file 'Bsplines95'.
!
! This driver should be linked with with bacoli95.f95, bacoli.f,
! bacoli-aux.f, and a problem definition file, e.g., rcdsys.f
!
!---------------------------------------------------------------------
! Example Python plotting code using serialized B-spline data:
! (Useful if one wants to use Python to obtain plots of approximate 
! solutions computed by BACOLI.)
! 
!   import matplotlib as mpl
!   mpl.use('AGG')  # for systems not running a GUI
!   import matplotlib.pyplot as plt
!   import numpy as np
!   from scipy.interpolate import splev
!
!   with open('Bsplines') as f: lines = f.readlines()
!   npde = int(lines[0])
!   xbs = np.fromstring(lines[1], sep=' ')
!   y   = np.fromstring(lines[2], sep=' ').reshape((npde,-1))
!   p = int(lines[3])
!
!   x = np.linspace(0, 1, 1000)
!   for i in range(npde):
!     u = splev(x, (xbs,y[i],p))
!     plt.plot(x, u, label='$u_{%d}$'%(i+1))
!
!   plt.xlabel('$x$')
!   plt.ylabel('$u$')
!   plt.title('$u(t,x)$ at $t=$tout')
!   plt.legend()
!   plt.grid()
!   plt.savefig('curve.png')
!
! scipy.interpolate.splev is an interface to SPLEV of FITPACK.
!-----------------------------------------------------------------------

program curve_example_driver

    use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals
    use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown
    use bacoli95_mod, only: bacoli95_sol_to_splines
    use bacoli95_mod, only: bacoli95_splines, bacoli95_splines_teardown

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(bacoli95_sol)     :: sol
    type(bacoli95_splines) :: splines

    ! Choose npde to be consistent with the problem specific source 
    ! code; see burg1.f, burg2.f, CahnAllen.f, rcdsys.f, sinmads.f, 
    ! and steady.f.
    integer,  parameter    :: npde = 4, nderiv = 1
    real(dp), parameter    :: xa = 0, xb = 1

    integer                :: nout
    real(dp), allocatable  :: xout(:), uout(:,:,:)
    real(dp)               :: tout, atol(npde), rtol(npde)
    logical                :: serialize
    !character(1)           :: choice 
	!Uncomment above line if user gets to choose to save B-spline 
	!info. See below.

    integer                :: i, j, ier

    external f, bndxa, bndxb, uinit

    ! problem specific variables:
    ! eps is used in the burg1.f, burg2.f and Cahn_Allen.f examples.
    real(dp) :: eps
    common /burger/ eps
    ! The following parameters are used in the rcdsys.f example.
    real(dp) :: a1,a2,d1,d2,r,c,n,pe1,pe2
    common /rcdstuff/ a1,a2,d1,d2,r,c,n,pe1,pe2
    eps = 1e-3_dp
    a1 = 30._dp
    a2 = 30._dp
    d1 = 1.5_dp
    d2 = 1.2_dp
    r = 1000._dp
    c = .96_dp
    n = 1._dp
    pe1 = 1e4_dp
    pe2 = 1e4_dp

    !-------------------------------------------------------------------
    ! Write out value of npde to allow user to confirm that its value
    ! is appropriate for the problem to be solved.
    
    write(6,*) 'The number of PDEs is assumed to be ', npde
    write(6,*)
    
    ! Set output time
    tout = 1.0d0

    ! Set tolerance
    atol = 1d-6
    rtol = atol
    atol(4) = 1d-1

    ! The following allows the user to decide whether to save solution
    ! as a spline for use with above python code.
    !serialize = .false.
    !print*, "Would you like B-spline information saved to disk? [y/N]"
    !read(*,*,err=600) choice
    !if (choice .eq. 'y' .or. choice .eq. 'Y') serialize = .true.

    ! Here we simply set serialize
    serialize = .true.

    !-------------------------------------------------------------------
    ! Initialization: Allocate storage and set problem parameters.
    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol)

    ! Set ouput at 11 uniformly spaced points across the spatial domain
    nout = 11
    allocate(xout(nout), uout(npde,nout,0:nderiv), stat=ier)
    if (ier /= 0 .or. sol%idid == -1000) goto 700

    !-------------------------------------------------------------------
    ! Integrate solution from t=0 to t=tout.
    print '(/"THE INPUT IS")'
    print 900, sol%kcol, sol%nint, sol%npde, tout
    print 901, sol%atol(1), sol%rtol(1), "LOI"

    ! Compute the solution at tout
    call bacoli95(sol, tout, f, bndxa, bndxb, uinit)

    ! Output idid to check for a successful computation
    print '("idid=",i5)',sol%idid
    if (sol%idid > 0) print '("idid > 0 => Successful computation")'

    !-------------------------------------------------------------------
    ! Output results.
    if (sol%idid > 0) then

            ! Set output points and evaluate solution at these points
            xout = (/xa,(xa+i*(xb-xa)/(nout-1), i=1,nout-1)/)
            call bacoli95_vals(sol, xout, uout, nderiv=nderiv)

            print '(/"THE OUTPUT IS")'
            print '(a13,a27)', "XOUT", "UOUT"
            do i = 1, nout
                print*, xout(i), uout(:,i,0)
            end do

            do j = 1, nderiv
                print '(/"SPATIAL PARTIAL DERIVATIVE",i2)', j
                print '(a13,a27)', "XOUT", "UOUT"
                do i = 1, nout
                    print*, xout(i), uout(:,i,j)
                end do
            end do

        ! Write B-spline information to file.
        if (serialize) then

            call bacoli95_sol_to_splines(sol, splines, ier)
            if (ier /= 0) goto 700

            open(unit=20,file='Bsplines')
            ! line 1: number of PDEs in the system
            write(20,*) splines%npde
            ! line 2: breakpoint/knot sequence
            write(20,*) splines%knots
            ! line 3: coefficient values at time t0
            write(20,*) splines%y
            ! line 4: degree of the B-spline interpolant
            write(20,*) splines%p
            close(20)

            call bacoli95_splines_teardown(splines)
        end if
    end if

    !-------------------------------------------------------------------
    ! The end.
    call bacoli95_sol_teardown(sol) ; stop
600 print '("Error: Improperly formatted input")' ; stop
700 print '("Error: Could not allocate storage")' ; stop

    !-------------------------------------------------------------------
    ! Formats
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)
902 format("Number of subintervals in the current mesh:",i8)
end program
