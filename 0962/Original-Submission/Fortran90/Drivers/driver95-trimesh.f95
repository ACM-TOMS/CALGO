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

! This driver writes output for each PDE using BACOLI's spatial mesh
! points at ntout uniform points along the temporal domain. tstart is 0,
! tstop is the end time, and the spatial domain is [xa,xb]. A number of
! parameters are hard coded and should be edited as needed: npde (the
! number of PDEs), xa, xb. Several other parameters take default values
! through the F95 wrapper.
!
! The output format used here allows the user to take advantage of
! BACOLI's spatial adaptivity in their data visualisation, to some
! extent, but the trade-off is that data points are slightly more
! difficult to work with than those over a uniform grid.
!
! If returning after uniform time steps is not desirable,
! it is possible to modify this code to have bacoli return after
! a specific number of time steps (as performed by DASSL).
! See the documentation inside bacoli95.f95 for details.
!
! After a successful computation, this program will have written files
! named Points001, ..., PointsNPDE containing data points as ordered
! triples, one per line. The coordinates are ordered (X, T, U).
!
! This driver should be linked with with bacoli.f, bacoli-aux.f,
! bacoli95.f95 and a problem definition file, eg, burg1.f.
!
!--------------------------------------------------------------------
! Example Python plotting code:
!
!   import matplotlib as mpl
!   mpl.use('AGG')  # for systems not running a GUI
!   from mpl_toolkits.mplot3d import Axes3D
!   from matplotlib import cm
!   import matplotlib.pyplot as plt
!   import numpy as np
!
!   styling = {
!     'cmap': cm.coolwarm,
!     'linewidth': 0,
!     'antialiased': True
!   }
!
!   x, t, u = np.loadtxt('Points1', unpack=True)
!   fig = plt.figure()
!   ax = fig.add_subplot(111, projection='3d')
!   ax.plot_trisurf(x, t, u, **styling)
!
!   ax.set_xlabel('$x$')
!   ax.set_ylabel('$t$')
!   ax.set_zlabel('$u(t,x)$')
!
!   plt.savefig('trimesh.png')
!
! axes3d.plot_trisurf() was added to matplotlib in version 1.2.0.
! It drapes a surface over a triangularization of the data.
!-----------------------------------------------------------------------

program trimesh_example_driver

    use bacoli95_mod, only: bacoli95_init, bacoli95, bacoli95_vals
    use bacoli95_mod, only: bacoli95_sol, bacoli95_sol_teardown

    implicit none
    integer, parameter     :: dp = kind(0d0)
    type(bacoli95_sol)     :: sol

    ! Choose npde to be consistent with the problem specific source 
    ! code; see burg1.f, burg2.f, CahnAllen.f, rcdsys.f, sinmads.f, 
    ! and steady.f.

    integer,  parameter    :: npde = 1
    real(dp), parameter    :: xa = 0, xb = 1
    real(dp), allocatable  :: uout(:,:)
    real(dp)               :: tout, tstart, tstop, atol(npde), rtol(npde)

    integer                :: i, j, k, ier, ntout
    character(len=32)      :: fname, npde_str

    external f, bndxa, bndxb, uinit

!----------------------------------------------------------------------
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
    
    ! Get some user input
    !print*, "Enter tstop, the end of the temporal domain."
    !read(*,*,err=600) tstop
    tstop = 1.0d0

    !print*, "At how many equally-spaced points along the time domain" &
    !     // " is output desired?"
    !read(*,*,err=600) ntout
    ntout = 4

    !print*, "Please choose an error tolerance"
    !read(*,*,err=600) atol(1)
    atol(1) = 1d-6
    rtol(1) = atol(1)

    !-------------------------------------------------------------------
    ! Initialization: Allocate storage and set problem parameters.
!    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol, dirichlet=1)
    call bacoli95_init(sol, npde, (/xa,xb/), atol=atol, rtol=rtol)

    allocate(uout(npde,sol%nint_max+1), stat=ier)
    if (ier /= 0 .or. sol%idid == -1000) goto 700
    tstart = sol%t0

    !-------------------------------------------------------------------
    ! Open files for output.
    do k = 1, npde
        write(npde_str,*) k
        fname = 'Points' // adjustl(trim(npde_str))
        open(unit=10+k,file=fname)
    end do

    !-------------------------------------------------------------------
    ! Integrate solution from t=0 to t=tout.
    print '(/"THE INPUT IS")'
    print 900, sol%kcol, sol%nint, sol%npde, tstop
    print 901, sol%atol(1), sol%rtol(1), "LOI"

    do j = 2, ntout

        tout = tstart + (j-1)*(tstop-tstart)/(ntout-1)

        call bacoli95(sol, tout, f, bndxa, bndxb, uinit)
        if (sol%idid <= 0) goto 800
        !print 902, sol%nint

        if (j == 2) then
            do i = 1, sol%nint+1
                call uinit(sol%x(i), uout(1,i), npde)
            end do
            do k = 1, npde
                do i = 1, sol%nint+1
                    write(10+k,*) sol%x(i), tstart, uout(k,i)
                end do
            end do
        end if

        call bacoli95_vals(sol, sol%x(1:sol%nint+1), uout)

        do k = 1, npde
            do i = 1, sol%nint+1
                write(10+k,*) sol%x(i), sol%t0, uout(k,i)
            end do
        end do
    end do

    print '("IDID       = ",i10)', sol%idid
    print '("nsteps     = ",i10)', sol%num_accepted_time_steps

    !-------------------------------------------------------------------

    call bacoli95_sol_teardown(sol) ; stop
600 print '("Error: Improperly formatted input")' ; stop
700 print '("Error: Could not allocate storage")' ; stop
800 print '("Error: Was not able to integrate to tsop")' ; stop

    !-------------------------------------------------------------------
    ! Formats!
900 format("kcol = ",i2,", nint0 = ",i4,", npde = ",i3,", tout = ",es7.1)
901 format("atol = ",es7.1,", rtol = ",es7.1,",",17x,a3)
902 format("Number of subintervals in the current mesh:",i8)

end program
