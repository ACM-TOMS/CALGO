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


! This file (bacoli95.f95) contains the BACOLI95 (Fortran95) source code
! that provides a convenient user interface to the BACOLI (Fortran77)
! solver whose source code is contained in the separate file, bacoli.f.

! In order to run an example, this source code must be linked with
! four other compilation units:
! 1. The BACOLI solver routines in bacoli.f
! 2. The supporting routines collected in bacoli-aux.f.
! 3. The supporting routines in d1mach_i1mach.f95.
! 4. A driver program, e.g., driver95-curve.f95.
! 5. PDE system definition subroutines, e.g., rcdsys.f.
! Note that the number of PDES, npde, specified in 4. and 5. must agree.

module bacoli95_mod
! This module provides a convenient Fortran 95 wrapper for BACOLI-3,
! the current version of bacoli.f, written in Fortran 77. In this
! source code BACOLI-3 will be referred to simply as BACOLI.

    implicit none
    integer, parameter :: dp = kind(0d0)
    integer, public    :: eunit = 6  ! Error output unit for this file.

    private

    type, public :: bacoli95_sol
    ! This structured type encapsulates problem parameters and solution
    ! information for a PDE system to be solved by BACOLI.
    !
    ! constant problem parameters
        integer           :: npde
    !                        npde is the number of partial differential
    !                        equations in the system being solved.
    !                        Set by bacoli95_init.
        integer           :: nint_max
    !                        A maximum number of spatial mesh
    !                        subintervals is used to (indirectly) limit
    !                        the amount of work and storage employed by
    !                        BACOLI. This is that maximum.
    !                        Set by bacoli95_init.
        integer           :: kcol
    !                        BACOLI uses collocation methods. kcol is
    !                        the number of collocation points per mesh
    !                        subinterval and determines the order of
    !                        piecewise polynomials used.
    !                        Set by bacoli95_init.
        integer           :: estimator
    !                        BACOLI has two options for spatial error
    !                        estimation schemes, LOI/LE (estimator=0)
    !                        and SCI/ST (estimator=1).
    !                        Set by bacoli95_init.
        integer           :: maxord
    !                        The underlying time integrator, DASSL, uses
    !                        multistep BDF methods. maxord may be used
    !                        to limit the maximum order of BDF method
    !                        used, which can be useful since higher
    !                        order methods have smaller stability
    !                        regions.
    !                        Set by bacoli95_init.
        real(dp), pointer :: atol(:)
    !                        The absolute error tolerance, either
    !                        of length 1 or length npde. BACOLI
    !                        and DASSL work to (independently)
    !                        control errors associated with
    !                        spatial discretization and time
    !                        stepping using these tolerances.
    !                        Set by bacoli95_init.
        real(dp), pointer :: rtol(:)
    !                        The relative error tolerance, either
    !                        of length 1 or length npde. BACOLI
    !                        and DASSL work to (independently)
    !                        control errors associated with
    !                        spatial discretization and time
    !                        stepping using these tolerances.
    !                        Set by bacoli95_init.
        real(dp)          :: ini_ss
    !                        The user may optionally specify an initial
    !                        stepsize for the underlying time integrator
    !                        to attempt. This is that initial stepsize.
    !                        Set by bacoli95_init.
    !
    ! variable problem parameters (input-output)
        real(dp)          :: t0
    !                        BACOLI constructs B-spline interpolants
    !                        of solution approximations in space for
    !                        specific values of the time variable, t.
    !                        t0 is the starting value of t on input,
    !                        and on output (after a call to bacoli95)
    !                        it represents the value of t associated
    !                        with the current B-spline interpolant.
    !                        Set by bacoli95_init and bacoli95.
        real(dp), pointer :: x(:)
    !                        On input, x is the initial spatial mesh on
    !                        which to attempt time stepping. As BACOLI
    !                        adapts the mesh as part of its error
    !                        control, this is modified on output.
    !                        Set by bacoli95_init and bacoli95.
        integer           :: nint
    !                        nint is the number of subintervals in the
    !                        spatial mesh, x. As BACOLI adapts the
    !                        mesh as part of its error control, this
    !                        number is variable.
    !                        Set by bacoli95_init and bacoli95.
        integer           :: mflag(12)
    !                        A vector of control flags for BACOLI.
    !                        These flags are handled by this wrapper
    !                        module. Consult BACOLI (bacoli.f)
    !                        documentation for details.
    !                        Set by bacoli95_init and bacoli95.
    !
    ! variable problem parameters (output)
        integer           :: idid
    !                        idid is the status / error code which
    !                        indicates if solution information is
    !                        available or the nature of an error.
    !                        It is used independently by bacoli95_init
    !                        and bacoli95/BACOLI to report errors in
    !                        initializations and integration,
    !                        respectively.
    !                        idid = 0 is the initial value set by a
    !                        successful call to bacoli95_init.
    !                        idid > 0 indicates that time integration
    !                        to the current value of t0 was successful.
    !                        -1000 < idid < 0 indicates an error arose
    !                        within the call to BACOLI, and a message
    !                        was printed to stdout (unit=6).
    !                        idid <= -1000 indicates that an error
    !                        occurred during the bacoli95_init call,
    !                        and a message was printed to the unit
    !                        eunit. eunit is a module variable,
    !                        initially 6 (for stdout).
    !                        For more details regarding idid, refer
    !                        to BACOLI (bacoli.f) documentation.
    !                        Set by bacoli95_init and bacoli95.
        real(dp), pointer :: y(:)
    !                        BACOLI approximates solutions using
    !                        B-splines, and y is the vector of
    !                        coefficients which, when combined with
    !                        a B-spline basis, provides solution
    !                        information. To extract solution
    !                        information, use bacoli95_vals or an
    !                        alternative described along with the
    !                        bacoli95_splines type (see below).
    !                        Set by bacoli95.
        real(dp), pointer :: rpar(:)
    !                        rpar is a large work array for BACOLI.
    !                        Care should be taken when accessing or
    !                        modifying this memory.
    !                        Set by bacoli95_init and bacoli95.
        integer,  pointer :: ipar(:)
    !                        ipar is a work array for BACOLI.
    !                        Care should be taken when accessing or
    !                        modifying this memory.
    !                        Set by bacoli95_init and bacoli95.
    !
    ! miscellaneous counters and values (output by bacoli95)
        integer           :: num_remeshings
    !                        Number of spatial remeshings.
        integer           :: num_ini_remeshings
    !                        Number of spatial remeshings on the
    !                        initial step.
        integer           :: num_cold_restarts
    !                        Number of cold restarts of DASSL.
    !                        When BACOLI adapts the spatial mesh,
    !                        the DAE system passed to DASSL is
    !                        different from the previous one and
    !                        a restart of some sort is required.
    !                        BACOLI attempts to do warm restarts
    !                        several times before resorting to
    !                        a cold restart, so these are a sign
    !                        of difficulties in the spatial error
    !                        estimation.
        integer           :: num_accepted_time_steps
    !                        Number of time steps taken by DASSL
    !                        whose spatial error estimate was
    !                        found to be within tolerances.
        integer           :: min_len_ipar
    !                        The number of locations within ipar
    !                        currently being used.
        integer           :: min_len_rpar
    !                        The number of locations within rpar
    !                        currently being used.
        integer           :: prev_bdf_order
    !                        The order of the BDF method employed
    !                        by DASSL on the most recent time step.
        real(dp)          :: prev_time_step_size
    !                        The size of the most recent time step
    !                        taken by DASSL.
    end type

    public :: bacoli95_init
    ! Initializes a bacoli95_sol type for the problem specified by
    ! input parameters.
    !
    !   subroutine bacoli95_init(sol, npde, x, tstart, atol, rtol, &
    !           kcol, nint_max, estimator, dirichlet, maxord, ini_ss)
    !
    ! Required Input:
    ! sol:       bacoli95_sol type to be initialized.
    !              If sol is being reinitialized, be sure to call
    !              bacoli95_sol_teardown on it prior to this subroutine.
    ! npde:      Number of partial differential equations in the system.
    ! x:         Initial spatial mesh OR left and right endpoints for
    !              which a uniform inital mesh should be generated.
    !
    ! Optional Input:
    ! tstart:    Initial time. Default tstart=0
    ! atol:      Error tolerance (absolute). May be a scalar or a vector
    !              of length npde. Default scalar atol=1e-6
    ! rtol:      Error tolerance (relative). May be a scalar or a vector
    !              of length npde. Default scalar rtol=1e-6
    ! kcol:      Number of collocation points used in the spatial
    !              discretization scheme. (2 < kcol < 11)
    !              kcol=4 is generally a good choice (and the default),
    !              but higher values may be more reliable for difficult
    !              problems, but come with a high cost in efficiency.
    !              The spatial error schemes frequently encounter
    !              difficulties when kcol=3, so that value is not
    !              recommended (though it can in some cases be the most
    !              efficient choice).
    ! nint_max:  Maximum number of subintervals the spatial mesh may use.
    !              The default is 500.
    ! estimator: Integer to specify whether the LOI/LE or the SCI/ST
    !              spatial error estimation scheme should be used to
    !              control mesh adaptivity.
    !              0 indicates LOI/LE (the default),
    !              1 indicates SCI/ST.
    !              Both schemes experience difficulty with kcol=3.
    ! dirichlet: Integer to specify whether or not both left and right
    !              boundary conditions of the PDE system are Dirichlet.
    !              Nonzero indicates that they are. The default is 0.
    !              If the problem has Dirichlet boundary conditions, it
    !              is slightly more straightforward to find consistent
    !              initial input for the underlying time integrator.
    ! maxord:    Maximum order of the Backward Differentiation Formula
    !              employed by the underlying time integrator, DASSL.
    !              The default is 5, and should rarely be changed, but
    !              there are some problems for which BDF methods are not
    !              stable, such as Schrodinger type problems with purely
    !              imaginary eigenvalues. Contained by 1 <= maxord <= 5.
    ! ini_ss:     Initial step size for the underlying time integrator.
    !              The default is to allow DASSL to pick the initial
    !              step size.

    public :: bacoli95
    ! Integrates the solution stored in sol forward in time until either
    ! a stop condition or an error condition is reached.
    !
    !   subroutine bacoli95(sol, tout, f, bndxa, bndxb, uinit, &
    !           derivf, difbxa, difbxb, tstop, nsteps)
    !
    ! Required Input:
    ! sol:    The initialized solution type, possibly modified by a
    !           previous call to bacoli95 in the case of continuation
    !           calls.
    ! tout:   Next specific output time. bacoli95 will return with
    !           sol%t0 <= tout so that solution information may be
    !           obtained (unless an error condition arose).
    !           bacoli95 may also successfully return due to other
    !           criteria specified by optional parameters.
    ! f:      User subroutine describing PDE system.
    ! bndxa:  User subroutine describing boundary conditions at left
    !           spatial endpoint.
    ! bndxb:  User subroutine describing boundary conditions at right
    !           spatial endpoint.
    ! uinit:  User subroutine describing initial values of the PDE
    !           system, u(tstart).
    !
    ! Optional Input:
    ! derivf: User subroutine describing partial derivatives of the
    !           PDE system. If omitted, approximate values are obtained
    !           through finite difference approximation. Providing
    !           analytic values allows more efficient execution, but is
    !           frequently unnecessary.
    ! difbxa: User subroutine describing partial derivatives of the left
    !           boundary conditions. May be approximated if omitted.
    ! difbxb: User subroutine describing partial derivatives of the right
    !           boundary conditions. May be approximated if omitted.
    ! tstop:  For efficient integration with dense output, the underlying
    !           time integrator (DASSL) will usually step past tout and
    !           interpolate back to provide solution information.
    !           If the problem, for example, has a singularity at tout,
    !           it may be necessary to prevent this behaviour and step
    !           to precisely tstop.
    ! nsteps: To return with solution information after DASSL has taken
    !           at most nsteps time steps, set this parameter.
    !
    ! Output:
    ! sol:      The fields of sol are updated in place. Refer to
    !             bacoli95_sol documentation for detailed comments.
    ! sol%t0:   Value of time variable for which approximated solution
    !             information is available.
    ! sol%idid: Status / Error code. A message should be printed in the
    !             event of an error.
    !             If sol%idid=1, the code returned before reaching tout.
    !             If sol%idid=2 or 3, the code returned with sol%t0=tout.
    !             If sol%idid<0, an error occurred and a message should
    !             have been printed.
    !             See BACOLI documentation for details.

    public :: bacoli95_vals
    ! Extract solution values at specific spatial output points for
    ! t=sol%t0. If nderiv is provided, that many spatial partial
    ! derivatives are also output.
    !
    !   subroutine bacoli95_vals(sol, xout, uout, nderiv)
    !
    ! Required Input:
    ! sol:    Solution data returned by a successful call to bacoli95.
    ! xout:   Output points within the bounds of the spatial mesh.
    !
    ! Optional Input:
    ! nderiv: Number of spatial partial derivatives to also output
    !           approximated values for. Default values is 0.
    !
    ! Output:
    ! uout:   Output values. May be of rank 1, 2 or 3, depending upon
    !           whether the number of PDEs in the system is greater
    !           than one and upon whether spatial partial derivative
    !           output is requested.
    !           If npde=1 and nderiv=0, use uout(1:size(xout)).
    !           uout(1:npde, 1:size(xout), 1:nderiv+1) is also a useful
    !           shape.

    public :: bacoli95_sol_teardown
    ! Deallocate and nullify pointers inside sol.
    !
    !   subroutine bacoli95_sol_teardown(sol)

! BACOLI represents the approximate solution that it computes as a
! piecewise polynomial (a B-spline). While this wrapper and BACOLI both
! provide routines for evaluating this spline at points of the users
! choice along the spatial domain, x, it can also be convenient to
! have access to the B-spline itself to perform the evaluation in other
! programs.
!
! The bacoli95_splines type contains only the information from the
! bacoli95_sol type which pertains to the B-spline at a specific time.
! This time is the t0 field of the bacoli95_sol type, and the t field
! of the bacoli95_splines type. It is hoped that this provides the user
! with more flexibility in manipulating solutions returned by bacoli95.
!
! For an example making use of the bacoli95_sol_to_splines routine,
! see driver95-curve.f95. In that example, B-spline information is
! saved to a file to be later read in and evaluated within a trivial
! Python script by splev of FITPACK by Paul Dierckx (which scipy
! provides an interface to through scipy.interpolate.splev). The example
! Python script is included within the comments of driver95-curve.f95.

    type, public :: bacoli95_splines
    ! This structured type is provided for conveniently exporting
    ! B-spline information produced by BACOLI.
        integer           :: npde
    !                        npde is the number of B-splines contained
    !                        within this type. There is a separate
    !                        spline for each solution component.
        integer           :: p
    !                        p is the order of the piecewise polynomial
    !                        used. This is equal to sol%kcol+1
        integer           :: nc
    !                        NC is the number of B-spline coefficients
    !                        and basis functions used for each B-spline.
        real(dp)          :: t
    !                        The B-spline coeffients, y, correspond to
    !                        a specific value of the time variable
    !                        (sol%t0). This is that value.
        real(dp), pointer :: y(:,:)
    !                        These are the B-spline coefficients,
    !                        taken from sol%y and reordered into a
    !                        rank 2 form for convenience. The first
    !                        index ranges 1:nc, and the second 1:npde.
        real(dp), pointer :: knots(:)
    !                        B-spline basis functions are tied to knots.
    !                        BACOLI uses spatial mesh points as knots,
    !                        so these are obtained by repeating values
    !                        of sol%x.
    end type

    public :: bacoli95_sol_to_splines
    ! Extracts solution information at t=sol%t0 from the solution type
    ! into a possibly more convenient structure containing B-spline
    ! information for the npde PDE solution components.
    !
    !   subroutine bacoli95_sol_to_splines(sol, splines, iflag)
    !
    ! Input:
    ! sol:     Solution data returned by a successful call to bacoli95.
    !
    ! Output:
    ! splines: Contains copies of information defining the B-splines if
    !            iflag=0.
    ! iflag:   Error flag for this subroutine.
    !            If iflag=0, there was no error.
    !            If iflag=-1, then sol does not appear to be usable.
    !            If iflag=-2, then there was a memory allocation error.

    public :: bacoli95_splines_teardown
    ! Deallocate and nullify pointers inside splines.
    !
    !   subroutine bacoli95_splines_teardown(splines)

!=======================================================================

    contains
       subroutine bacoli95_init(sol, npde, x, tstart, atol, rtol, &
                kcol, nint_max, estimator, dirichlet, maxord, ini_ss)
            ! parameters
            type(bacoli95_sol)              :: sol
            integer,  intent(in)            :: npde
            real(dp), intent(in)            :: x(:)
            real(dp), intent(in), optional  :: tstart
            real(dp), intent(in), optional  :: atol(:)
            real(dp), intent(in), optional  :: rtol(:)
            integer,  intent(in), optional  :: kcol
            integer,  intent(in), optional  :: nint_max
            integer,  intent(in), optional  :: estimator
            integer,  intent(in), optional  :: dirichlet
            integer,  intent(in), optional  :: maxord
            real(dp), intent(in), optional  :: ini_ss
            ! locals
            integer,  parameter :: dflt_nint     = 10
            integer,  parameter :: dflt_nint_max = 500
            integer,  parameter :: dflt_kcol     = 4
            real(dp), parameter :: dflt_tol      = 1d-6
            integer :: i, ier, maxlrp, maxlip, maxy !,ntol

            ! Initialise sol.
            sol%idid = 0
            sol%mflag = 0
            sol%npde = npde
            sol%estimator = 0
            nullify(sol%x, sol%y, sol%atol, sol%rtol, sol%rpar, sol%ipar)

            ! Verify input and check for scalar or vector tolerance use.
            if (size(x) < 2) then
               sol%idid = -1001 ;                           goto 99
            end if
            if (present(nint_max)) then
                if (nint_max < 1) then
                   sol%idid = -1002 ;                       goto 99
                end if
                if (size(x) > nint_max+1) then
                   sol%idid = -1003 ;                       goto 99
                end if
            end if
            !ntol = 1
            !if (present(atol)) then
            !    ntol = size(atol)
            !end if
            !if (present(rtol)) then
            !    ntol = max(ntol, size(rtol))
            !end if
            if (present(estimator)) then
                if (estimator /= 0 .and. estimator /= 1) then
                   sol%idid = -1004 ;                       goto 99
                 end if
            end if
            if (present(maxord)) then
                if (maxord < 1 .or. 5 < maxord) then
                   sol%idid = -1005 ;                       goto 99
                end if
            end if

            ! Set fields of sol according to input parameters.
            if (present(kcol)) then
                sol%kcol = kcol
            else
                sol%kcol = dflt_kcol
            end if

            if (present(nint_max)) then
                sol%nint_max = nint_max
            else
                sol%nint_max = max(dflt_nint_max, size(x)-1)
            end if

            allocate(sol%x(sol%nint_max+1), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if
            if (size(x) == 2) then
                ! Create a uniform initial mesh of dflt_nint
                ! subintervals from x(1) to x(2).
                sol%x(1:dflt_nint+1) = &
                        (/ (i*(x(2)-x(1))/dflt_nint, i=0,dflt_nint) /)
                sol%nint = dflt_nint
            else
                sol%x(1:size(x)) = x
                sol%nint = size(x)-1
            end if

            if (present(tstart)) then
                sol%t0 = tstart
            else
                sol%t0 = 0
            end if

            !if (present(atol)) then
            !    allocate(sol%atol(ntol), stat=ier)
            !    if (ier /= 0) then
            !       sol%idid = -1000 ;                       goto 99
            !    end if
            !    if (size(atol) == 1) then
            !        sol%atol = atol(1)
            !    else
            !        sol%atol = atol
            !    end if
            !else
            !    allocate(sol%atol(ntol), stat=ier)
            !    if (ier /= 0) then
            !        sol%idid = -1000 ;                      goto 99
            !    end if
            !    sol%atol = dflt_tol
            !end if
            !
            !if (present(rtol)) then
            !    allocate(sol%rtol(ntol), stat=ier)
            !    if (ier /= 0) then
            !        sol%idid = -1000 ;                      goto 99
            !    end if
            !    if (size(rtol) == 1) then
            !        sol%rtol = rtol(1)
            !    else
            !        sol%rtol = rtol
            !    end if
            !else
            !    allocate(sol%rtol(ntol), stat=ier)
            !    if (ier /= 0) then
            !        sol%idid = -1000 ;                      goto 99
            !    end if
            !    sol%rtol = dflt_tol
            !end if

            allocate(sol%atol(npde), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if
            if (present(atol)) then
                sol%atol = atol(1)
                if (any(atol /= atol(1))) sol%mflag(2) = 1
            else
                sol%atol = dflt_tol
            end if

            allocate(sol%rtol(npde), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if
            if (present(rtol)) then
                sol%rtol = rtol(1)
                if (any(rtol /= rtol(1))) sol%mflag(2) = 1
            else
                sol%rtol = dflt_tol
            end if

            !if (ntol > 1) sol%mflag(2) = 1

            if (present(estimator)) then
                sol%estimator = estimator
                sol%mflag(8) = estimator
            end if

            if (present(dirichlet)) then
                if (dirichlet /= 0) sol%mflag(5) = 1
            end if

            if (present(maxord)) then
                sol%mflag(7) = 1
                sol%maxord = maxord
            end if

            if (present(ini_ss)) then
                if (ini_ss > 0) then
                    sol%mflag(6) = 1
                    sol%ini_ss = ini_ss
                end if
            end if

            call calc_array_sizes(sol%npde, sol%nint_max, sol%kcol, &
                    sol%estimator, maxlrp, maxlip, maxy)
            allocate(sol%y(maxy), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if
            allocate(sol%ipar(maxlip), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if
            allocate(sol%rpar(maxlrp), stat=ier)
            if (ier /= 0) then
                sol%idid = -1000 ;                          goto 99
            end if

            ! Initialise sol%rpar and sol%ipar to zero.
            sol%rpar = 0 ; sol%ipar = 0

            return
         99 continue
            ! Error section; free associated pointers and complain.
            call bacoli95_sol_teardown(sol)
            select case (sol%idid)
                case(-1000)
                    write(eunit,*) 'allocation error in bacoli95_init'
                case(-1001)
                    write(eunit,*) 'x must contain left and right ends'
                case(-1002)
                    write(eunit,*) 'nint_max < 1'
                case(-1003)
                    write(eunit,*) 'size(x) > nint_max+1'
                case(-1004)
                    write(eunit,*) 'estimator not 0 or 1'
                    write(eunit,*) ' use 0 to choose LOI/LE error control'
                    write(eunit,*) ' use 1 to choose SCI/ST error control'
                case(-1005)
                    write(eunit,*) 'maxord outside [1,5]'
                case default
                    write(eunit,*) 'error of unknown type'
            end select
        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95(sol, tout, f, bndxa, bndxb, uinit, &
                derivf, difbxa, difbxb, tstop, nsteps)
            ! parameters
            type(bacoli95_sol)                :: sol
            real(dp), intent(in)              :: tout
            real(dp), intent(in), optional    :: tstop
            integer,  intent(in), optional    :: nsteps
            external                             f
            external                             bndxa
            external                             bndxb
            external                             uinit
            external                             derivf
            external                             difbxa
            external                             difbxb
            optional  derivf, difbxa, difbxb
            external                             bacoli
            ! constants (names for indices and offsets within sol%ipar)
            integer, parameter :: iipstp = 6
            integer, parameter :: irpstp = 7
            integer, parameter :: irshin = 9
            integer, parameter :: isteps = 10
            integer, parameter :: irmesh = 11
            integer, parameter :: icolds = 14
            integer, parameter :: irwork = 39
            integer, parameter :: iiwork = 95
            integer, parameter :: ipbdf  = 8
            integer, parameter :: ihprv  = 7
            ! locals
            real(dp) :: outt
            integer  :: iflag
            logical  :: use_finite_differences_both
            logical  :: use_finite_differences_pdes
            logical  :: use_finite_differences_bconds

            if (sol%idid < 0) then
                write(eunit,*) 'BACOLI previously returned with an error'
                write(eunit,*) 'If it was rectified, set idid to 1'
                write(eunit,*) 'and try again.'
                return
            end if

            ! Check the lengths allocated to arrays in sol
            call check_array_sizes(sol, iflag)
            if (iflag == 1) then
                write(eunit,*) 'WARNING: array lengths in sol are short'
            else if (iflag /= 0) then
                write(eunit,*) 'ERROR: array lengths too short to call bacoli'
                return
            end if

            ! BACOLI allows for the finite difference approximation
            ! of Jacobians of the PDEs and the boundary conditions
            ! independently. sol%mflag(9) communicates what
            ! approximations are required like so:
            ! sol%mflag(9) = 0: approximate for both
            ! sol%mflag(9) = 1: approximate for b.conds.
            ! sol%mflag(9) = 2: approximate for pdes
            ! sol%mflag(9) = 3: use analytic Jacobians
            use_finite_differences_both   = .true.
            use_finite_differences_pdes   = .true.
            use_finite_differences_bconds = .true.
            if (present(derivf)) then
                use_finite_differences_pdes = .false.
                use_finite_differences_both = .false.
            end if
            if (present(difbxa)) then
                if (present(difbxb)) then
                    use_finite_differences_bconds = .false.
                    use_finite_differences_both   = .false.
                end if
            end if

            if (sol%mflag(1) == 1) sol%idid = 1
            if (present(tstop)) then
                ! BACOLI complains if tstop < tout
                outt = min(tout, tstop)
                sol%mflag(3) = 1
            else
                outt = tout
                sol%mflag(3) = 0
            end if
            if (present(nsteps)) then
                sol%mflag(4) = 1
                sol%ipar(8) = nsteps
            else
                sol%mflag(4) = 0
            end if

            if (present(tstop)) then
                if (sol%mflag(3) == 1) sol%rpar(1)  = tstop
            end if
            if (sol%mflag(6) == 1) sol%rpar(2)  = sol%ini_ss
            if (sol%mflag(7) == 1) sol%ipar(15) = sol%maxord

            if (use_finite_differences_both) then
                sol%mflag(9) = 0
                call bacoli(sol%t0, outt, sol%atol, sol%rtol, sol%npde,&
                            sol%kcol, sol%nint_max, sol%nint, sol%x,   &
                            sol%mflag, sol%rpar, size(sol%rpar),       &
                            sol%ipar, size(sol%ipar), sol%y, sol%idid, &
                            f, dummy_derivf, bndxa, dummy_difbx, bndxb,&
                            dummy_difbx, uinit)
            else if (use_finite_differences_bconds) then
                sol%mflag(9) = 1
                call bacoli(sol%t0, outt, sol%atol, sol%rtol, sol%npde,&
                            sol%kcol, sol%nint_max, sol%nint, sol%x,   &
                            sol%mflag, sol%rpar, size(sol%rpar),       &
                            sol%ipar, size(sol%ipar), sol%y, sol%idid, &
                            f, derivf, bndxa, dummy_difbx, bndxb,      &
                            dummy_difbx, uinit)
            else if (use_finite_differences_pdes) then
                sol%mflag(9) = 2
                call bacoli(sol%t0, outt, sol%atol, sol%rtol, sol%npde,&
                            sol%kcol, sol%nint_max, sol%nint, sol%x,   &
                            sol%mflag, sol%rpar, size(sol%rpar),       &
                            sol%ipar, size(sol%ipar), sol%y, sol%idid, &
                            f, dummy_derivf, bndxa, difbxa, bndxb,     &
                            difbxb, uinit)
            else
                sol%mflag(9) = 3
                call bacoli(sol%t0, outt, sol%atol, sol%rtol, sol%npde,&
                            sol%kcol, sol%nint_max, sol%nint, sol%x,   &
                            sol%mflag, sol%rpar, size(sol%rpar),       &
                            sol%ipar, size(sol%ipar), sol%y, sol%idid, &
                            f, derivf, bndxa, difbxa, bndxb, difbxb,   &
                            uinit)
            end if

            ! Update counters. This is merely for convenience; the truly
            ! determined user can, alternatively, read BACOLI source code
            ! and access appropriate indices of sol%rpar and sol%ipar.
            sol%num_remeshings          = sol%ipar(irmesh) + sol%ipar(irshin)
            sol%num_ini_remeshings      = sol%ipar(irshin)
            sol%num_cold_restarts       = sol%ipar(icolds)
            sol%num_accepted_time_steps = sol%ipar(isteps)
            sol%min_len_ipar            = sol%ipar(iipstp)
            sol%min_len_rpar            = sol%ipar(irpstp)
            sol%prev_bdf_order          = sol%ipar(iiwork-1+ipbdf)
            sol%prev_time_step_size     = sol%rpar(sol%ipar(irwork)-1+ihprv)

        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95_vals(sol, xout, uout, nderiv)
            type(bacoli95_sol)             :: sol
            real(dp), intent(in)           :: xout(:)
            real(dp), intent(out)          :: uout(*)
            integer,  intent(in), optional :: nderiv
            external values
            integer  :: nd, ier, lenwrk
            real(dp), allocatable :: work(:)

            if (present(nderiv)) then
                nd = nderiv
            else
                nd = 0
            end if

            lenwrk = (sol%kcol+2)*(nd+1) + sol%kcol*(sol%nint+1) + 4
            allocate(work(lenwrk), stat=ier)
            if (ier /= 0) then
                write(eunit,*) 'unable to allocate work storage for values'
                return
            end if

            call values(sol%kcol, xout, sol%nint, sol%x, sol%npde, &
                        size(xout), nd, uout, sol%y, work)
        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95_sol_to_splines(sol, splines, iflag)
            type(bacoli95_sol)     :: sol
            type(bacoli95_splines) :: splines
            integer, intent(out)   :: iflag
            integer :: npde, nint, kcol, ier, i, ii, j

            iflag = 0

            npde = sol%npde
            nint = sol%nint
            kcol = sol%kcol

            nullify(splines%y, splines%knots)

            ! Check sol%idid to determine if the values within sol
            ! are usable.
            if (sol%idid < 1) then
                iflag = -1 ;                                goto 99
            end if

            splines%npde = npde
            splines%p    = kcol+1
            splines%nc   = nint*kcol+2
            splines%t    = sol%t0

            allocate(splines%y(nint*kcol+2, npde), stat=ier)
            if (ier /= 0) then
                iflag = -2 ;                                goto 99
            end if

            ! Reorder the B-spline coefficients for convenient access
            ! across individual splines (corresponding to PDE components).
            splines%y = transpose(reshape(sol%y(1:npde*(nint*kcol+2)), &
                                          (/npde, nint*kcol+2/)))

            allocate(splines%knots((nint+1)*kcol+4), stat=ier)
            if (ier /= 0) then
                iflag = -2 ;                                goto 99
            end if

            ! The knot sequence (also known as a breakpoint sequence)
            ! is simply a repetition of the mesh points.
            ! There are extra copies of the leftmost mesh point.
            do i = 1, kcol+2
                splines%knots(i)             = sol%x(1)
                splines%knots(nint*kcol+2+i) = sol%x(nint+1)
            end do
            do i = 2, nint
                ii = (i-2)*kcol+kcol+2
                do j = 1, kcol
                    splines%knots(ii+j) = sol%x(i)
                end do
            end do

            return
         99 continue
            ! Error section; free associated pointers and complain.
            select case (iflag)
                case(-1)
                    write(eunit,*) 'sol does not appear to be valid'
                case(-2)
                    call bacoli95_splines_teardown(splines)
                    write(eunit,*) 'allocation error in bacoli95_sol_to_splines'
                case default
                    write(eunit,*) 'error of unknown type'
            end select
        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95_sol_teardown(sol)
            type(bacoli95_sol) :: sol
            if (associated(sol%x))    then
                deallocate(sol%x) ;    nullify(sol%x)
            end if
            if (associated(sol%y))    then
                deallocate(sol%y) ;    nullify(sol%y)
            end if
            if (associated(sol%atol)) then
                deallocate(sol%atol) ; nullify(sol%atol)
            end if
            if (associated(sol%rtol)) then
                deallocate(sol%rtol) ; nullify(sol%rtol)
            end if
            if (associated(sol%rpar)) then
                deallocate(sol%rpar) ; nullify(sol%rpar)
            end if
            if (associated(sol%ipar)) then
                deallocate(sol%ipar) ; nullify(sol%ipar)
            end if
        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95_splines_teardown(splines)
            type(bacoli95_splines) :: splines
            if (associated(splines%y)) then
                deallocate(splines%y) ;     nullify(splines%y)
            end if
            if (associated(splines%knots)) then
                deallocate(splines%knots) ; nullify(splines%knots)
            end if
        end subroutine

!-----------------------------------------------------------------------

        subroutine bacoli95_sol_trim(sol)
            ! Currently unused.
            type(bacoli95_sol) :: sol
            integer            :: minx, miny, minlrp, minlip, ier
            real(dp), pointer  :: tmp_x(:), tmp_y(:), tmp_rpar(:)
            integer,  pointer  :: tmp_ipar(:)

            if  ((.not.associated(sol%x))    &
            .or. (.not.associated(sol%y))    &
            .or. (.not.associated(sol%rpar)) &
            .or. (.not.associated(sol%ipar)))               return

            minx = sol%nint+1
            call calc_array_sizes(sol%npde, sol%nint, sol%kcol, &
                    sol%estimator, minlrp, minlip, miny)

            allocate(tmp_x(minx), stat=ier)
            if (ier /= 0)                                   return
            tmp_x = sol%x(1:minx)
            deallocate(sol%x)
            sol%x => tmp_x
            nullify(tmp_x)

            allocate(tmp_ipar(minlip), stat=ier)
            if (ier /= 0)                                   return
            tmp_ipar = sol%ipar(1:minlip)
            deallocate(sol%ipar)
            sol%ipar => tmp_ipar
            nullify(tmp_ipar)

            allocate(tmp_y(miny), stat=ier)
            if (ier /= 0)                                   return
            tmp_y = sol%y(1:miny)
            deallocate(sol%y)
            sol%y => tmp_y
            nullify(tmp_y)

            allocate(tmp_rpar(minlrp), stat=ier)
            if (ier /= 0)                                   return
            tmp_rpar = sol%rpar(1:minlrp)
            deallocate(sol%rpar)
            sol%rpar => tmp_rpar
            nullify(tmp_rpar)
        end subroutine

!-----------------------------------------------------------------------

        subroutine calc_array_sizes(npde, nint, kcol, est, &
                lrp, lip, ly)
            ! Helper routine.
            integer, intent(in)  :: npde, nint, kcol, est
            integer, intent(out) :: lrp, lip, ly

            ly = npde*(kcol*nint+2)

            lip = 115 + ly

            lrp = 113 + 59*npde + 27*nint + 13*npde*npde   &
                + 9*kcol + 24*kcol*nint + 6*nint*kcol*kcol &
                + 27*npde*nint*kcol + 7*nint*npde          &
                + 2*npde*npde*nint*kcol*kcol               &
                + 4*npde*npde*kcol*nint
            ! The LOI/LE scheme requires slightly less storage.
            if (est == 0) then
                lrp = lrp - 15*nint + 3*kcol - 8*kcol*nint &
                    + kcol*kcol - nint*kcol*kcol           &
                    - 3*nint*npde + 2
            end if
        end subroutine

!-----------------------------------------------------------------------

        subroutine check_array_sizes(sol, iflag)
            ! Helper routine.
            ! parameters
            type(bacoli95_sol)   :: sol
            integer, intent(out) :: iflag
            ! iflag = 0 : array lengths okay
            ! iflag = 1 : array lengths dangerously short for sol%nint_max
            ! iflag = 2 : array lengths below min (may segfault in bacoli)

            ! locals
            integer  :: minx,   midx,   maxx
            integer  :: miny,   midy,   maxy
            integer  :: minlrp, midlrp, maxlrp
            integer  :: minlip, midlip, maxlip

            ! Check the sizes of sol%{x,y,rpar,ipar}.
            ! Warn if they are less than min(1.5*minsize, maxsize),
            ! and error if they are smaller than minsize.
            minx = sol%nint+1 ; maxx = sol%nint_max+1
            call calc_array_sizes(sol%npde, sol%nint, sol%kcol, &
                    sol%estimator, minlrp, minlip, miny)
            call calc_array_sizes(sol%npde, sol%nint_max, sol%kcol, &
                    sol%estimator, maxlrp, maxlip, maxy)
            midx   = min(int(1.5*minx),   maxx)
            midy   = min(int(1.5*miny),   maxy)
            midlrp = min(int(1.5*minlrp), maxlrp)
            midlip = min(int(1.5*minlip), maxlip)

            iflag = 0
            if      ((size(sol%x)    < minx)    &
                .or. (size(sol%y)    < miny)    &
                .or. (size(sol%rpar) < minlrp)  &
                .or. (size(sol%ipar) < minlip)) then
                iflag = 1
            else if ((size(sol%x)    < midx)    &
                .or. (size(sol%y)    < midy)    &
                .or. (size(sol%rpar) < midlrp)  &
                .or. (size(sol%ipar) < midlip)) then
                iflag = 2
            end if
        end subroutine

!-----------------------------------------------------------------------

        subroutine dummy_derivf(t, x, u, ux, uxx, dfdu, dfdux, &
              dfduxx, npde)
           ! Dummy - used if user does not provide derivf.
           integer,  intent(in)  :: npde
           real(dp), intent(in)  :: t, x, u(npde)
           real(dp), intent(in)  :: ux(npde), uxx(npde)
           real(dp), intent(out) :: dfdu(npde,npde)
           real(dp), intent(out) :: dfdux(npde,npde)
           real(dp), intent(out) :: dfduxx(npde,npde)

           dfdu = 0 ; dfdux = 0 ; dfduxx = 0
        end subroutine

!-----------------------------------------------------------------------

        subroutine dummy_difbx(t, u, ux, dbdu, dbdux, dbdt, npde)
           ! Dummy - used if user does not provide difbxa, difbxb.
           integer,  intent(in)  :: npde
           real(dp), intent(in)  :: t, u(npde), ux(npde)
           real(dp), intent(out) :: dbdu(npde,npde)
           real(dp), intent(out) :: dbdux(npde,npde)
           real(dp), intent(out) :: dbdt(npde)

           dbdu = 0 ; dbdux = 0 ; dbdt = 0
        end subroutine

end module


