c BACOLI-3 (Modification of BACOL with interpolation-based spatial 
c error estimation)

c (BACOL is Copyright (c) 2002, Rong Wang, Pat Keast, Paul Muir
c Rong Wang, School of Mathematics and Statistics, Wuhan University,
c Pat Keast, Department of Mathematics and Statistics,
c                            Dalhousie University (retired),
c Paul Muir, Mathematics and Computing Science, Saint Mary's University.)

c BACOLI-3 is Copyright (c) 2013, Paul Muir, Jack Pew
c Paul Muir, Mathematics and Computing Science, Saint Mary's University.
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions
c are met:
c * Redistributions of source code must retain the above copyright
c   notice, this list of conditions and the following disclaimer.
c * Redistributions in binary form must reproduce the above copyright
c   notice, this list of conditions and the following disclaimer in the
c   documentation and/or other materials provided with the distribution.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
c HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
c LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
c THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
c OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c This file contains the BACOLI (Fortran77) source code (including the 
c new code that implements interpolation based spatial error estimates.)
c


      subroutine bacoli(t0, tout, atol, rtol, npde, kcol, nintmx, nint,
     &                  x, mflag, rpar, lrp, ipar, lip, y, idid,
     &                  f, derivf, bndxa, difbxa, bndxb, difbxb, uinit)

c-----------------------------------------------------------------------
c Purpose:
c       The purpose of BACOLI is to solve NPDE dimensional systems of
c       second order parabolic partial differential equations (PDEs)
c       in one space variable of the form:
c
c            dU            /                                \
c            -- (t,x) = f | t, x, U(t,x), U (t,x), U   (t,x) | ,
c            dt            \               x         xx     /
c
c       where x_a < x < x_b and t > t0, with initial conditions at
c       time t = t0 are given by:
c
c                          u(t0,x) = u_0(x),
c
c       for x_a <= x <= x_b, subject to separated boundary conditions
c       given by:
c
c                         /                      \
c                   b    | t, U(t,x_a), U (t,x_a) | = 0,
c                    x_a  \              x       /
c
c                         /                      \
c                   b    | t, U(t,x_b), U (t,x_b) | = 0,
c                    x_b  \              x       /
c
c       for t > t0 and x = x_a, x = x_b, respectively.
c
c       Guide to the above notation:
c          dU
c          -- (t,x) = denotes the first partial derivative of U(t,x)
c          dt         with respect to the time variable t.
c
c          U (t,x) = denotes the first partial derivative of U(t,x)
c           x        with respect to space variable x.
c
c          U  (t,x) = denotes the second partial derivative of U(t,x)
c           xx       with respect to space variable x.
c
c       Furthermore, the above functions are NPDE dimensional vector 
c       functions.
c
c       BACOLI is a method of lines algorithm which uses bspline
c       collocation to discretize the spatial domain [x_a,x_b].
c       The output is a vector of bspline coefficients which
c       can be used to calculate the approximate solution U(t,x) and
c       its spatial derivatives at (tout,x) where x_a <= x <= x_b
c       and t0 < tout.
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Setup of BACOLI:
c       BACOLI requires that the user specifies the system of PDEs and
c       the related initial and boundary conditions as well as setting
c       input parameters (which define the bspline space and the
c       requested error tolerances) and allocating work storage.
c
c       The calling sequence of BACOLI is:
c
c       call bacoli(t0, tout, atol, rtol, npde, kcol, nint, nintmx, x,
c    &             mflag, rpar, lrp, ipar, lip, y, idid,
c    &             f, derivf, bndxa, difbxa, bndxb, difbxb, uinit)
c
c       which will generate the vector y = Y(tout) upon successful
c       completion. Generally, the call to BACOLI will be followed by a
c       call to VALUES to calculate the solution at a set of points:
c
c       call values(kcol, xsol, nint, x, npde, npts, nderiv,
c     &             usol, y, work)
c
c       The details of the parameters to VALUES are documented within
c       the source code for that routine. The input parameters for
c       BACOLI are dealt with in detail below, but a quick summary is:
c
c       [t0, tout] is the time domain of the problem.
c       NPDE is the number of components in the PDE system.
c       atol is the absolute error tolerance.
c       rtol is the relative error tolerance.
c       kcol, nint, and x define the bspline space.
c       nintmx is the maximum number of subintervals allowed.
c       mflag(1:9) is used to control the operation of BACOLI.
c       rpar(lrp) is a floating point work array.
c       ipar(lip) is an integer work array.
c
c       The user must check idid to determine what further action needs
c       to be taken.
c
c       PDE system definition subroutines are to be set up as follows:
c
c       f(t, x, u, ux, uxx, fval, npde)
c       npde is an integer, the rest are double precision.
c       t, x are scalars, u, ux, uxx are vectors of length npde.
c       This subroutine defines the right hand side of the PDE system,
c       and ut = f(t, x, u, ux, uxx) should be returned in fval.
c
c       derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx, npde)
c       npde is an integer, the rest are double precision.
c       t, x are scalars, u, ux, uxx are vectors of length npde.
c       dfdu, dfdux, dfduxx are Jacobians of f evaluated at (t, x).
c
c       bndxa(t, u, ux, bval, npde)
c       npde is an integer, the rest are double precision.
c       t is scalar, u, ux, bval are vectors of length npde.
c       This subroutine defines the left boundary conditions,
c       where x = xa, by b(t, u, ux) = 0.
c       Return the residual b(t, u, ux) in bval.
c
c       difbxa(t, u, ux, dbdu, dbdux, dbdt, npde)
c       npde is an integer, the rest are double precision.
c       t is scalar, u, ux are vectors of length npde.
c       This subroutine defines the differentiated left boundary,
c       where x = xa.
c       Return the values of the Jacobians of bndxa evaluated at t in
c       dbdu, dbdux and dbdt.
c
c       bndxb and difbxb are the same as bndxa and difbxa, but for the
c       right boundary conditions where x = xb.
c
c       uinit(x, u, npde)
c       npde is an integer, x is a double precision scalar and
c       u is a double precision vector of length npde.
c       This subroutine defines the initial condition of the PDE system
c       by returning the value of u0(x) = u(t0, x) in u.
c       This initial condition should be C1-continuous.
c
c-----------------------------------------------------------------------
        implicit none
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               On input, t0 < tout is the initial
c                               time. On output, t0 is the current time,
c                               t0 <= tout.
c
        double precision        tout
c                               tout is the desired final output time.
c
c                               After a successful return from BACOLI,
c                               the time stepping may be resumed by
c                               changing tout so that t0 < tout and
c                               setting mflag(1) = 1 to indicate a
c                               continuation of the previous problem.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
c
c                               If PDE components vary in importance,
c                               then vector error tolerances may be
c                               used by setting mflag(2) = 1. In this
c                               case, the dimension of atol must be
c                               npde. The user will define atol(1),
c                               atol(2), ..., atol(npde) appropriately.
c                               Note that a change from scalar to vector
c                               tolerances (or vice versa) constitutes
c                               a new problem, and BACOLI will have to
c                               be reinitialized.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag(2) = 0.
c
c                               If PDE components vary in importance,
c                               then vector error tolerances may be
c                               used by setting mflag(2) = 1. In this
c                               case, the dimension of rtol must be
c                               npde. The user will define rtol(1),
c                               rtol(2), ..., rtol(npde) appropriately.
c                               Note that a change from scalar to vector
c                               tolerances (or vice versa) constitutes
c                               a new problem, and BACOLI will have to
c                               be reinitialized.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                                  2 < kcol <= mxkcol
c                               The degree of the piecewise polynomial
c                               is (kcol+1).
c
        integer                 nint
c                               on input, nint is the number of
c                               subintervals defined by the spatial
c                               mesh x at the initial time t0.
c                               on output, nint is the number of
c                               subintervals at tout.
c                               nint >= 1.
c
        integer                 nintmx
c                               the maximum number of subintervals that
c                               the user requires.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a,x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c                               at input, x(1:nint+1) stores the mesh
c                               points at the initial time t0.
c                               at output, x(1:nint+1) stores the mesh
c                               points at tout.
c
        integer                 mflag(12)
c                               This vector determines the interaction
c                               of BACOLI with DASSL and which error
c                               estimation scheme BACOLI will employ.
c
c                                         How to set mflag(1):
c                               On the initial call to BACOLI with a
c                               new problem, set mflag(1) = 0, which
c                               indicates that BACOLI and DASSL should
c                               perform the initialization steps that
c                               are required by each code, respectively.
c
c                               In order to continue time stepping in
c                               the current problem after a successful
c                               return from BACOLI, set mflag(1) = 1,
c                               idid = 1, and ensure that t0 < tout.
c
c                                         How to set mflag(2):
c                               If scalar absolute and relative error
c                               tolerances (atol and rtol) are desired,
c                               then set mflag(2) = 0.
c
c                               For vector absolute and relative error
c                               tolerances, set mflag(2) = 1, define
c                               atol(1), ..., atol(npde), and
c                               rtol(1), ..., rtol(npde), as described
c                               above, ensuring that the dimension of
c                               each of atol and rtol is at least npde.
c
c                                         How to set mflag(3):
c                               If there are no restrictions on t0,
c                               then set mflag(3) = 0.
c                               Since DASSL may actually "step over"
c                               tout and then interpolate, there is the
c                               option to enforce tstop >= tout.
c                               If this is desirable, set mflag(3) = 1,
c                               and define rpar(1) = tstop.
c
c                                         How to set mflag(4):
c                               If the user wishes, BACOLI will return
c                               the computed solution and derivative
c                               after a certain number of accepted time
c                               steps or at TOUT, whichever comes first.
c                               This is a good way to proceed if the
c                               user wants to see the behavior of the
c                               solution.
c                               If the user only wants the solution at
c                               TOUT, set mflag(4) = 0;
c                               else, set mflag(4) = 1, and assign a
c                               positive integer for ipar(8) that
c                               defines the number of time steps before
c                               BACOLI is stopped.
c
c                                         How to set mflag(5):
c                               If both boundary conditions are
c                               dirichlet, set mflag(5) = 1;
c                                    else, set mflag(5) = 0.
c
c                                         How to set mflag(6):
c                               If the user wants to specify an initial
c                               stepsize, set mflag(6) = 1,
c                                             and define rpar(2) = the
c                                             initial stepsize;
c                                   else, set mflag(6) = 0.
c
c                                         How to set mflag(7):
c                               If the user wants to use DASSL with the
c                               BDF methods of maximum order to default
c                               to 5, set mflag(7) = 0;
c                               else, set mflag(7) = 1, and define
c                                     ipar(15) = the maximum order.
c
c                                         How to set mflag(8):
c                               If the user wants to run BACOLI with the
c                               LOI error estimation scheme,
c                                  set mflag(8) = 0;
c                               SCI error estimation scheme,
c                                  set mflag(8) = 1.
c
c                                         How to set mflag(9):
c                               If the user wants to BACOLI use finite
c                               difference approximations in place of
c                               the user-provided analytic partial
c                               derivative subroutines derivf, difbxa
c                               and difbxb,
c                                  set mflag(9) = 0;
c                               If the user implements derivf, but not
c                               difbxa or difbxb,
c                                  set mflag(9) = 1;
c                               If the user implements difbxa and
c                               difbxb, but not derivf,
c                                  set mflag(9) = 2;
c                               If the user implements all of derivf,
c                               difbxa and difbxb,
c                                  set mflag(9) = 3.
c
c                               mflag(10:12): reserved for future use.
c
        integer                 lrp
c                               lrp is the size of the rpar storage
c                               array and must satisfy:
c                               lrp>=113+59*npde+27*nintmx
c    +                               +13*npde*npde+9*kcol+24*kcol*nintmx
c    +                               +6*nintmx*kcol*kcol+7*nintmx*npde
c    +                               +27*npde*nintmx*kcol
c    +                               +2*npde*npde*nintmx*kcol*kcol
c    +                               +4*npde*npde*kcol*nintmx
c    -                               -15*nintmx+3*kcol+kcol*kcol
c    -                               -8*kcol*nintmx-3*nintmx*npde
c    -                               -nintmx*kcol*kcol+2
c                               where the group of lines starting with
c                               a `-' should not be included when
c                               mflag(8) = 1. The SCI uses slightly
c                               more storage than the LOI.
c
        integer                 lip
c                               lip is the size of the ipar integer
c                               work array and must satisfy:
c                               lip>=115+npde*(nintmx*kcol+2)
c
        external                f
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
        external                uinit
c                               See "Setup of BACOLI" above.
c
c       Work Storage:
        double precision        rpar(lrp)
c                               rpar is a floating point work array
c                               of size lrp.
c
        integer                 ipar(lip)
c                               ipar is an integer work array
c                               of size lip.
c
c       Output:
        double precision        y(npde*(kcol*nintmx+2))
c                               On successful return from BACOLI,
c                               y(1:npde*(kcol*nint+2)) is
c                               the vector of bspline coefficients at
c                               the current time.
c
        integer                 idid
c                               idid is the BACOLI exit status flag
c                               which is based on the exit status from
c                               DASSL plus some additional status codes
c                               based on error checking performed by
c                               BACOLI on initialization. Positive
c                               values of idid indicate a successful
c                               return. Negative values of idid indicate
c                               an error which may or may not be fatal.
c                               The exact descriptions of idid return
c                               values will be discussed below.
c
c                               For calls other than the first call
c                               (mflag(1) = 1), the user must check
c                               idid, set idid = 1 (if necessary), and
c                               take other actions which are necessary
c                               such as defining a new value of tout.
c
c                               An excerpt from the DASSL source code
c                               documentation is included to define
c                               the idid return codes and clarify
c                               the operation of DASSL within BACOLI.
c
c-----------------------------------------------------------------------
c  The following is an excerpt from the DASSL source code documentation:
c
C  -------- OUTPUT -- AFTER ANY RETURN FROM DASSL ---------------------
C
C  The principal aim of the code is to return a computed solution at
C  TOUT, although it is also possible to obtain intermediate results
C  along the way. To find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the IDID parameter.
C
C  IDID -- Reports what the code did.
C
C                     *** Task completed ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- A step was successfully taken in the
C                   intermediate-output mode. The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- The integration to TSTOP was successfully
C                   completed (T=TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- The integration to TOUT was successfully
C                   completed (T=TOUT) by stepping past TOUT.
C                   Y(*) is obtained by interpolation.
C                   YPRIME(*) is obtained by interpolation.
C
C                    *** Task interrupted ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- A large amount of work has been expended.
C                   (About 500 steps)
C
C           IDID = -2 -- The error tolerances are too stringent.
C
C           IDID = -3 -- The local error test cannot be satisfied
C                   because you specified a zero component in ATOL
C                   and the corresponding computed solution
C                   component is zero. Thus, a pure relative error
C                   test is impossible for this component.
C
C           IDID = -6 -- DASSL had repeated error test
C                   failures on the last attempted step.
C
C           IDID = -7 -- The corrector could not converge.
C
C           IDID = -8 -- The matrix of partial derivatives
C                   is singular.
C
C           IDID = -9 -- The corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           IDID =-10 -- The corrector could not converge
C                   because IRES was equal to minus one.
C
C           IDID =-11 -- IRES equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           IDID =-12 -- DASSL failed to compute the initial
C                   YPRIME.
C
C                    *** Task terminated ***
C                Reported by the value of IDID=-33
C
C           IDID = -33 -- The code has encountered trouble from which
C                   it cannot recover. A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. For example, this occurs
C                   when invalid input is detected.
C
C  RTOL, ATOL -- These quantities remain unchanged except when
C               IDID = -2. In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration. However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
C                    (CALLS AFTER THE FIRST)
C
C  This code is organized so that subsequent calls to continue the
C  integration involve little (if any) additional effort on your
C  part. You must monitor the IDID parameter in order to determine
C  what to do next.
C
C  Recalling that the principal task of the code is to integrate
C  from T to TOUT (the interval mode), usually all you will need
C  to do is specify a new TOUT upon reaching the current TOUT.
C
C  Do not alter any quantity not specifically permitted below,
C  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
C  or the differential equation in subroutine RES. Any such
C  alteration constitutes a new problem and must be treated as such,
C  i.e., you must start afresh.
C
C  You cannot change from vector to scalar error control or vice
C  versa (mflag(2)), but you can change the size of the entries of
C  RTOL, ATOL. Increasing a tolerance makes the equation easier
C  to integrate. Decreasing a tolerance will make the equation
C  harder to integrate and should generally be avoided.
C
C  If it has been necessary to prevent the integration from going
C  past a point TSTOP (mflag(3), rpar(itstop)), keep in mind that
C  the code will not integrate to any TOUT beyond the currently
C  specified TSTOP. Once TSTOP has been reached you must change
C  the value of TSTOP or set mflag(3)=0. You may change mflag(3)
C  or TSTOP at any time but you must supply the value of TSTOP in
C  rpar(itstop) whenever you set mflag(3)=1.
C
C  -------- ERROR MESSAGES ---------------------------------------------
C
C      The SLATEC error print routine XERMSG is called in the event of
C   unsuccessful completion of a task.  Most of these are treated as
C   "recoverable errors", which means that (unless the user has directed
C   otherwise) control will be returned to the calling program for
C   possible action after the message has been printed.
C
C   In the event of a negative value of IDID other than -33, an appro-
C   priate message is printed and the "error number" printed by XERMSG
C   is the value of IDID.  There are quite a number of illegal input
C   errors that can lead to a returned value IDID=-33.  The conditions
C   and their printed "error numbers" are as follows:
C
C   Error number       Condition
C
C        1       Some element of INFO vector is not zero or one.
C        2       NEQ .le. 0
C        3       MAXORD not in range.
C        4       LRW is less than the required length for RWORK.
C        5       LIW is less than the required length for IWORK.
C        6       Some element of RTOL is .lt. 0
C        7       Some element of ATOL is .lt. 0
C        8       All elements of RTOL and ATOL are zero.
C        9       INFO(4)=1 and TSTOP is behind TOUT.
C       10       HMAX .lt. 0.0
C       11       TOUT is behind T.
C       12       INFO(8)=1 and H0=0.0
C       13       Some element of WT is .le. 0.0
C       14       TOUT is too close to T to start integration.
C       15       INFO(4)=1 and TSTOP is behind T.
C       16       --( Not used in this version )--
C       17       ML illegal.  Either .lt. 0 or .gt. NEQ
C       18       MU illegal.  Either .lt. 0 or .gt. NEQ
C       19       TOUT = T.
c----------------------------------------------------------------------
C
C   If DASSL is called again without any action taken to remove the
C   cause of an unsuccessful return, XERMSG will be called with a fatal
C   error flag, which will cause unconditional termination of the
C   program.  There are two such fatal errors:
C
C   Error number -998:  The last step was terminated with a negative
C       value of IDID other than -33, and no appropriate action was
C       taken.
C
C   Error number -999:  The previous call was terminated because of
C       illegal input (IDID=-33) and there is illegal input in the
C       present call, as well.  (Suspect infinite loop.)
C
c-----------------------------------------------------------------------
c Notes regarding modifications (Jack Pew):
c
c This version of BACOL, BACOLI, replaces the computation of a
c second collocation solution, previously required for error estimation,
c with one of two interpolants constructed from the first solution.
c
c The first scheme (LOI) replaces the low-order collocation solution
c with an interpolant of the same order and compares it against the
c high-order collocation solution. This interpolant is designed so that
c its error term is asymptotically the same as that of the collocation
c solution it replaces.
c
c The second scheme (SCI) replaces the high-order collocation solution
c computed by BACOL with an interpolant of the same order. The
c interpolant harvests superconvergent data from the low-order solution.
c
c Here is a list of changes made to BACOL:
c
c   BACOL -> BACOLI:
c
c An additional flag was added to the mflags vector to allow the user
c to select which interpolation scheme to use:
c mflag(8) = 0 for the LOI scheme, mflag(8) = 1 for the SCI scheme.
c
c The structure of and several pointers into both rpar and ipar were
c altered. There is no longer a space in these work arrays for 'Ubar',
c and rpar(ipar(iebas2)) was replaced with rpar(ipar(iecoef)), which
c stores basis function values and Hermite-Birkhoff coefficients used
c in the interpolation-based error estimates between remeshings.
c The principle of the old iebas2 and the new iecoef is the same.
c
c To remove the computation of the second collocation solution,
c DASSL was given a smaller number of equations.
c That is, neq is now what was previously neq1 - the number of equations
c in the primary computed solution.
c Additionally, JAC, RES, DDASLV and DDAJAC were modified to only call
c CALJAC, CALRES, CRSLVE and CRDCMP once each (this is where the two
c solutions were previously decoupled.)
c
c SUCSTP is now called to store the computed B-spline coefficients for
c the single collocation solution that is computed instead of only the
c higher order solution (it no longer makes sense to refer to a higher
c order and lower order solution.)
c
c REINIT is now called with kcol instead of kcol+1, as there is now
c only one computed solution for which we need to perform interpolation
c of history data after a warm start. Also, it is now only called once,
c as there is no longer a second solution to reinitialize.
c
c INIY and INIYP are no longer called a second time to initialize
c B-spline coefficients for a second solution.
c
c ERREST was heavily modified, replacing the second call to ERRVAL with
c a call to either LOWINT or SCINT, and computing the values of several
c new pointers into the section of rpar that was added to store B-spline
c basis function evaluations and Hermite-Birkhoff coefficients that are
c used in the interpolation subroutines.
c
c A minor bug in BACOL was fixed where, if mflag(3) were set,
c the value of tstop in rwork(itstop) became garbage on a cold start
c which could be caught as invalid input for DASSL.
c
c   BACOLI -> BACOLI-2
c
c BACOLI now takes the user subroutine names, which define the PDE system,
c as external parameters. These parameters had to be added to several
c BACOLI and DASSL subroutines to have them passed down to CALRES,
c CALJAC, INIY and INIYP.
c
c The D1MACH implementation from SLATEC was replaced with an alternative
c implementation from http://www.netlib.org/blas/d1mach.f which attempts
c to discern appropriate machine constants at runtime.
c
c Subroutines were shuffled around within this file in an attempt to
c better organise them.
c
c   BACOLI-2 -> BACOLI-3
c
c BACOLI can now use finite differences to approximate iteration
c matrices by default (mflag(9) = 0). The dummy subroutines derivf,
c difbxa and difbxb will only be called if mflag(9) is set to indicate
c that they should be.
c mflag(9) = 1 indicates that derivf is available,
c mflag(9) = 2 indicates that difbxa and difbxb are available, and
c mflag(9) = 3 indicates that all derivative subroutines are available.
c For the purposes of INIY and INIYP, the subroutine FDBNDX was
c developed to replace difbxa and difbxb via central differences on
c bndxa and bndxb. DDASSL was modified further to take a new MTYPE = 6
c which signals that approximated ABD linear algebra should be done
c (MTYPE = 3 is for dense ABD linear algebra), and DDAJAC was modified 
c to implement this through a minimal number of calls to RES.
c (This number is the number of columns in a central block of the ABD
c system: npde*(kcol+2).)
c
c A few calls to DCOPY within the BACOLI subroutine were replaced with
c loops to avoid potential aliasing issues.
c
c More space was added to mflag and ipar (reserving it for future use).
c This space in ipar exists in the two exclusive ranges of
c (imflg9=16, ih=31) and (iey=56, idasi=60).
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 MAXORD
        parameter              (MAXORD = 5)
c                               MAXORD is the maximum order of the
c                               backward differentiation formula (BDF)
c                               methods used by DASSL. MAXORD = 5 is
c                               the default used by DASSL.
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
        integer                 maxrsh
        parameter              (maxrsh = 20)
c                               maxrsh is the maximum number of
c                               remesh times at one time step,
c                               i.e., icount must less than or equal
c                               to maxrsh
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Local variables:
c
        integer                 kerr
c                               kerr is what kcol would be for the
c                               solution whose error is being estimated.
c                               This helps explain how the LOI error
c                               estimate is working.
c
        integer                 neq
c                               neq=npde*ncpts is the number of
c                               bspline coefficients (or DAEs) when
c                               using dassl_{kcol}.
c
        integer                 leniw
c                               leniw = 20 + neq is the length of the
c                               integer work array required by dassl.
c
        integer                 lenpd
c                               lenpd is the size of the Almost Block
c                               Diagonal (ABD) Jacobian.
c                               lenpd=npde*npde*(2*nconti
c                                      +kcol*(kcol+nconti)*nint)
c
        integer                 lenrw
c                               lenrw = 40+(MAXORD+4)*neq+lenpd
c                               is the total size of the floating point
c                               work array required by dassl.
c
        integer                 lenin
c                               lenin is the size of the floating
c                               point work array used by INIY and INIYP.
c                               lenin>=lenpd+2*neq+npde*4+2*npde*npde
c
        integer                 lenri
c                               lenri is the size of the floating
c                               point work array used by REINIT.
c
        integer                 lenrj
c                               lenrj is the size of the floating
c                               point work array used by RES and JAC.
c                               lenrj>=6*npde+5*npde*npde.
c
        integer                 lencof
c                               Length of the work storage array for
c                               Hermite-Birkhoff coefficients and
c                               values of B-spline basis functions used
c                               inside errest's interpolation
c                               subroutines. They only need to be
c                               recalculated when the mesh changes.
c                               The SCI estimate requires more storage
c                               than the LOI estimate here.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the total
c                               number of collocation points when using
c                               dassl_{kcol}.
c
        integer                 necpts
c                               necpts is the total number
c                               of quadrature points used for
c                               error estimate.
c
        integer                 icflag
c                               This is the status flag from the almost
c                               block diagonal factorization routine,
c                               CRDCMP.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        torign
c                               torign is the initial time, i.e. = t0
c                               at the beginning.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 isstep
c                               isstep is the number of accepted time
c                               steps since we restart BACOLI in the
c                               case that mflag(4) = 1.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before the current remeshing.
c
        integer                 ninpre
c                               ninpre is the number of subintervals
c                               when icount = 0 before remeshing.
c
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial step or any step
c                                           without needing remesh;
c                                      = 1, remesh with a hot start.
c                                      = 2, remesh with a cold start.
c
        integer                 neqpre
c                               neqpre is the number of bspline
c                               coefficients when icount = 0 before
c                               remeshing when using dassl_{kcol}.
c
        integer                 irold
c                               irold is the value of ipar(ixold) before
c                               remeshing.
c
        integer                 nstep
c                               nstep is the number of steps which need
c                               to be interpolated for a warm start.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 itstop
        parameter              (itstop =  1)
c                               rpar(itstop) = tstop as defined when
c                               mflag(3) = 1.
c
        integer                 iiniss
        parameter              (iiniss =  2)
c                               rpar(iiniss) = the initial stepsize when
c                               mflag(6) = 1.
c
        integer                 irpstr
        parameter              (irpstr = 11)
c                               rpar(1:irpstr-1) are reserved to store
c                               floating point scalar quantities.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpts
        parameter              (incpts =  4)
c                               ipar(incpts) = ncpts.
c
        integer                 ineq
        parameter              (ineq   =  5)
c                               ipar(ineq) = neq.
c
        integer                 iipstp
        parameter              (iipstp =  6)
c                               ipar(iipstp) = the minimum size of ipar.
c
        integer                 irpstp
        parameter              (irpstp =  7)
c                               ipar(irpstp) = the minimum size of rpar.
c
        integer                 iinstp
        parameter              (iinstp =  8)
c                               As input, ipar(iinstp) is the number of
c                               intermediate time steps before the user
c                               asks BACOLI to return the computed
c                               solution.
c                               As output, ipar(iinstp) will keep the
c                               same value unless TOUT is reached before
c                               ipar(iinstp) time steps. If TOUT is
c                               reached, ipar(iinstp) will be the number
c                               of time steps before TOUT is reached.
c                               See mflag(4) for more details.
c
        integer                 irshin
        parameter              (irshin =  9)
c                               ipar(irshin) is the number of remeshing
c                               times at the initial step.
c
        integer                 isteps
        parameter              (isteps = 10)
c                               ipar(isteps) is the number of time steps
c                               on the current problem.
c
        integer                 irmesh
        parameter              (irmesh = 11)
c                               ipar(irmesh) is the number of remeshing
c                               times after BACOLI starts the initial
c                               step.
c
        integer                 istalr
        parameter              (istalr = 12)
c                               ipar(istalr) is the number of accepted
c                               steps after the last successful
c                               remeshing.
c
        integer                 istblc
        parameter              (istblc = 13)
c                               ipar(istblc) is the number of steps
c                               BACOLI has taken before the latest cold
c                               start.
c
        integer                 icolds
        parameter              (icolds = 14)
c                               ipar(icolds) is the number of the times
c                               when BACOLI performs a cold start after
c                               remeshing.
c
        integer                 imxord
        parameter              (imxord = 15)
c                               ipar(imxord) is the maximum order for
c                               the BDF methods employed in DASSL.
c                               The default value is 5.
c
        integer                 imflg9
        parameter              (imflg9 = 16)
c                               ipar(imflg9) stores the value of
c                               mflag(9) in order for it to be passed
c                               down to caljac.
c
        integer                 idasi
        parameter              (idasi  = 60)
c                               ipar(idasi) stores, before remeshing,
c                               the first 20 elements of the integer
c                               point work array in dassl.
c
        integer                 iinfo
        parameter              (iinfo  = 80)
c                               ipar(iinfo) is an integer array
c                               required by dassl_{kcol}.
c                               ipar(iinfo)   = mflag(1),
c                               ipar(iinfo+1) = mflag(2),
c                               ipar(iinfo+3) = mflag(3),
c                               ipar(iinfo+7) = mflag(6).
c                               ipar(iinfo+8) = mflag(7).
c                               (See the documentation of DDASSL for
c                               details.)
c
        integer                 iiwork
        parameter              (iiwork = 95)
c                               ipar(iiwork) is the integer work array
c                               for dassl.
c
        integer                 ipivot
        parameter              (ipivot = 115)
c                               ipar(ipivot-1+i), i = 1, neq, contains
c                               the pivoting information from the
c                               factorization of the temporary matrix
c                               for dassl.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ih
        parameter              (ih     = 31)
c                               rpar(ipar(ih)) stores the mesh step
c                               size sequence.
c
        integer                 ixcol
        parameter              (ixcol  = 32)
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               dassl_{kcol}.
c
        integer                 ixbs
        parameter              (ixbs   = 33)
c                               rpar(ipar(ixbs)) stores the breakpoint
c                               sequence when using dassl_{kcol}.
c
        integer                 iy
        parameter              (iy     = 34)
c                               rpar(ipar(iy)) stores the vector of
c                               solution components to the DAE system
c                               when using dassl_{kcol}.
c
        integer                 iyp
        parameter              (iyp    = 35)
c                               rpar(ipar(iyp)) stores the vector of
c                               solution component derivatives of the
c                               DAE system when using dassl_{kcol}.
c
        integer                 iabtop
        parameter              (iabtop = 36)
c                               rpar(ipar(iabtop)) stores the top block
c                               of the ABD collocation matrices when
c                               using dassl_{kcol}.
c
        integer                 iabblk
        parameter              (iabblk = 37)
c                               rpar(ipar(iabblk)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_{kcol}.
c
        integer                 iabbot
        parameter              (iabbot = 38)
c                               rpar(ipar(iabbot)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using dassl_{kcol}.
c
        integer                 irwork
        parameter              (irwork = 39)
c                               rpar(ipar(irwork)) stores the floating
c                               point work array for DASSL, and is
c                               also used as work storage for the
c                               subroutines iniy, iniyp and reinit.
c
        integer                 iwkrj
        parameter              (iwkrj  = 40)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by RES and JAC.
c
        integer                 ibasi
        parameter              (ibasi  = 41)
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using dassl_{kcol}.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) stores
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 itatol
        parameter              (itatol = 42)
c                               rpar(ipar(itatol)) is the absolute error
c                               tolerance request on the integration
c                               error.
c
        integer                 itrtol
        parameter              (itrtol = 43)
c                               rpar(ipar(itrtol)) is the relative error
c                               tolerance request on the integration
c                               error.
c
        integer                 iexcol
        parameter              (iexcol = 44)
c                               rpar(ipar(iexcol)) stores the Gaussian
c                               quadrature points which are used for
c                               the L2-norm error estimate.
c
        integer                 iewts
        parameter              (iewts  = 45)
c                               rpar(ipar(iewts)) stores the Gaussian
c                               weights which are used for error
c                               estimate.
c
        integer                 iebasi
        parameter              (iebasi = 46)
c                               rpar(ipar(iebasi)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol}.
c
        integer                 iecoef
        parameter              (iecoef = 47)
c                               rpar(ipar(iecoef)) stores the values
c                               Hermite-Birkhoff coefficients and
c                               nonzero B-spline basis functions at
c                               points used in the error estimate.
c
        integer                 iercom
        parameter              (iercom = 48)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 49)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 50)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 iwkdnm
        parameter              (iwkdnm = 51)
c                               rpar(ipar(iwkdnm)) is the work storage
c                               for the modification version of the
c                               subroutine DDANRM.
c
        integer                 ixold
        parameter              (ixold  = 52)
c                               rpar(ipar(ixold)) stores the mesh point
c                               sequence when icount = 0 before
c                               remeshing.
c
        integer                 idasr
        parameter              (idasr  = 53)
c                               rpar(ipar(idasr)) stores, before
c                               remeshing, the first 40 elements of the
c                               floating point work array in dassl.
c
        integer                 iypre
        parameter              (iypre  = 54)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy2)) at the previous 6 steps.
c                               It is required for a hot restart after
c                               remeshing.
c
        integer                 iyprer
        parameter              (iyprer = 55)
c                               rpar(ipar(iyprer)) stores the
c                               information of at the previous steps
c                               after remeshing.
c
        integer                 iey
        parameter              (iey    = 56)
c                               rpar(ipar(iey)) stores the bspline
c                               coefficients at the farthest point that
c                               integration has reached.
c
c-----------------------------------------------------------------------
c Subroutines Passed to DASSL:
        external                jac
        external                res
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               colpnt
c                               ddassl
c                               ddatrp
c                               divdif
c                               errest
c                               iniy
c                               iniyp
c                               meshsq
c                               reinit
c                               remesh
c                               sucstp
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c
c Last modified by Jack Pew, August 9, 2013.
c
c-----------------------------------------------------------------------

c     Check validity of the mflag vector.
      do 1 i = 1, 8
         if ((mflag(i) .lt. 0) .or. (mflag(i) .gt. 1)) goto 710
   1  continue
      if ((mflag(9) .lt. 0) .or. (mflag(9) .gt. 3)) goto 711

      if (mflag(4) .eq. 1) then
         isstep = 0
         if (ipar(iinstp) .le. 0) goto 715
      endif

      irshfg = 0
      istart = 0
      icount = 0
      torign = t0
      ninpre = nint

c     The SCI scheme estimates the error for the computed collocation
c     solution, while the LOI scheme estimates the error for a
c     collocationsolution of one order lower. (In LOI mode, the code
c     computes an estimate of the error that would be obtained were the
c     code to compute a collocation solution that was one order lower.
c     This is similar to local extrapolation for initial value
c     solutions.)
c     kerr keeps track of this for calls to meshsq and remesh.
      if (mflag(8) .eq. 0) then
c        LOI scheme
         kerr = kcol - 1
      elseif (mflag(8) .eq. 1) then
c        SCI scheme
         kerr = kcol
      endif

c     mflag(9) is loaded into ipar so that it may be communicated to
c     caljac. This value can in theory change between calls to BACOLI.
      ipar(imflg9) = mflag(9)

c     Check for continuation of a previous problem.
      if (mflag(1) .eq. 1) then

         istart = 1

         neq   = ipar(ineq)
         lenpd = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)
         leniw = 20 + neq
         lenrw = 40 + (MAXORD + 4) * neq + lenpd
         if (mflag(8) .eq. 0) then
            necpts = (kcol + 2) * nint
            lenerr = (2 * necpts + nint) * npde
     &             + ((kcol - 3) * nint + (nint + 1) * 2) * npde
            lencof = (kcol + 1) * (kcol + 2)
     &             + (kcol + nconti) * (kcol - 3) * nint
     &             + (kcol + nconti) * 2 * (nint + 1)
         elseif (mflag(8) .eq. 1) then
            necpts = (kcol + 3) * nint
            lenerr = (2 * necpts + nint) * npde
     &             + ((kcol - 2) * nint + (nint + 1) * 2) * npde
            lencof = (kcol + 4) * (kcol + 3) * nint
     &             + (kcol + nconti) * (kcol - 2) * nint
     &             + (kcol + nconti) * 2 * (nint + 1)
         endif

         goto 200
      else

c        Check if the user specifies an initial stepsize
         if (mflag(6) .eq. 1) then
            if (((tout-t0)*rpar(iiniss)) .lt. zero) goto 720
            if (rpar(iiniss) .eq. zero) goto 725
         endif

         ipar(irmesh) = 0
         ipar(irshin) = 0
         ipar(istalr) = 0
         ipar(istblc) = 0
         ipar(icolds) = 0
      endif

c-----------------------------------------------------------------------
c     On the initial call or after remeshing, check for valid input and
c     initialize the workspace.
c-----------------------------------------------------------------------

  100 continue

c     Check validity of npde, kcol, and nint.
      if (npde .le. 0) goto 730
      if ((kcol .le. 2) .or. (kcol .gt. mxkcol)) goto 740
      if ((nint .le. 0) .or. (nint .gt. nintmx)) goto 640

c     Check for a monotone mesh.
      do 110 i = 1, nint
         if (x(i) .ge. x(i+1)) goto 760
  110 continue

c-----------------------------------------------------------------------
c     Calculate the extra storage requirements of res and jac.
      lenrj = (6 + 5 * npde) * npde

c     Calculate the number of collocation points.
      ncpts = nint * kcol + nconti

c     Calculate the number of DAEs given to dassl.
      neq = npde * ncpts

c     Size of the ABD iteration matrix in dassl.
      lenpd = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)

c     Calculate the extra storage requirements of iniy and iniyp.
      lenin = lenpd + 2 * neq + 2 * npde * (2 + npde)

c     Calculate the extra storage requirements of reinit.
      lenri = lenpd + kcol + nconti + kcol * (ninpre + 1)
     &         + 2 * nconti

c-----------------------------------------------------------------------
c     Total size of the DASSL floating point work array.
      lenrw = 40 + (MAXORD + 4) * neq + lenpd

c     Total size of the DASSL integer work array.
      leniw = 20 + neq

c-----------------------------------------------------------------------
c     Calculate the number of quadrature points used for error
c     estimate and the extra storage requirements of errest.
c     The SCI scheme needs slightly more storage here.

      if (mflag(8) .eq. 0) then
         necpts = (kcol + 2) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 3) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 1) * (kcol + 2)
     &          + (kcol + nconti) * (kcol - 3) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      elseif (mflag(8) .eq. 1) then
         necpts = (kcol + 3) * nint
         lenerr = (2 * necpts + nint) * npde
     &          + ((kcol - 2) * nint + (nint + 1) * 2) * npde
         lencof = (kcol + 4) * (kcol + 3) * nint
     &          + (kcol + nconti) * (kcol - 2) * nint
     &          + (kcol + nconti) * 2 * (nint + 1)
      endif

c-----------------------------------------------------------------------
c     Save the input parameters in the ipar integer communication
c     storage array.
      ipar(inpde)  = npde
      ipar(ikcol)  = kcol
      ipar(inint)  = nint
      ipar(incpts) = ncpts
      ipar(ineq)   = neq

c-----------------------------------------------------------------------
c     Calculate the offsets into the rpar floating point storage array.
c-----------------------------------------------------------------------
      ipar(itatol) = irpstr
      ipar(itrtol) = ipar(itatol) + npde
      ipar(ih)     = ipar(itrtol) + npde

      ipar(iy)     = ipar(ih)     + nint
      ipar(iyp)    = ipar(iy)     + neq

      ipar(ixcol)  = ipar(iyp)    + neq
      ipar(ixbs)   = ipar(ixcol)  + ncpts
      ipar(iabtop) = ipar(ixbs)   + ncpts + kcol + nconti
      ipar(iabblk) = ipar(iabtop) + npde * npde * nconti
      ipar(iabbot) = ipar(iabblk) + npde * npde * nint * kcol
     &                            * (kcol + nconti)
      ipar(ibasi)  = ipar(iabbot) + npde * npde * nconti

      ipar(irwork) = ipar(ibasi)  + (kcol + nconti) * 3 * ncpts
      ipar(iwkrj)  = ipar(irwork) + lenrw

      ipar(iexcol) = ipar(iwkrj)  + lenrj
      ipar(iewts)  = ipar(iexcol) + necpts
      ipar(ierint) = ipar(iewts)  + necpts
      ipar(iercom) = ipar(ierint) + nint
      ipar(iebasi) = ipar(iercom) + npde
      ipar(iecoef) = ipar(iebasi) + (kcol + nconti) * necpts
      ipar(iework) = ipar(iecoef) + lencof

      ipar(iwkdnm) = ipar(iework) + lenerr

      ipar(iyprer) = ipar(iwkdnm) + neq
      ipar(ixold)  = ipar(iyprer) + 6 * neq
      ipar(idasr)  = ipar(ixold)  + nintmx + 1
      ipar(iypre)  = ipar(idasr)  + 40

      ipar(iey)    = ipar(irwork) + 40 + 3 * neq

c     This offset is different between the initial call and remeshing.
      if ((irshfg .ne. 0) .and. (istart .eq. 1)) then
         ipar(irpstp) = ipar(iypre) + 6 * neqpre - 1
      else
         ipar(irpstp) = ipar(iypre) + 6 * neq - 1
      endif

c     Check for a sufficiently large rpar floating point work array.
      if (lrp .lt. ipar(irpstp)) goto 770

c     Calculate the offsets into the integer storage array.
      ipar(iipstp) = ipivot + neq - 1

c     Check for a sufficiently large ipar integer work array.
      if (lip .lt. ipar(iipstp)) goto 780

c     Check whether it is initial call or for remeshing.
c     A different set of initializations is done in each case.
      if ((irshfg .ne. 0) .and. (istart .ne. 0)) goto 300

c-----------------------------------------------------------------------
c     Perform initializations.
c-----------------------------------------------------------------------
c     Check that this is not just a remeshing at the initial step.
      if (icount .eq. 0) then

c        Set the info vector required by DASSL.
         ipar(iinfo)   = mflag(1)
         ipar(iinfo+1) = mflag(2)
         ipar(iinfo+2) = 1
         ipar(iinfo+3) = mflag(3)

         do 120 i = 6, 14
            ipar(iinfo+i-1) = 0
  120    continue

c        Indicate whether the maximum order of BDF method is restricted.
         ipar(iinfo+8) = mflag(7)

c        Indicate whether an user-supplied initial stepsize is used.
         ipar(iinfo+7) = mflag(6)

c        Indicate to DASSL whether or not an analytic Jacobian
c        (iteration) matrix is supplied. For the purposes of DASSL,
c        the Jacobian matrix is analytic if the user implements derivf.
         if (mflag(9).eq.1 .or. mflag(9).eq.3) then
            ipar(iinfo+4) = 1
         else
            ipar(iinfo+4) = 0
         endif

c        Indicate an ABD Jacobian (Iteration) matrix.
         ipar(iinfo+14) = 1

      else
c        This is a remeshing at the initial step. Increment the counter.
         ipar(irshin) = ipar(irshin) + 1
      endif

      call meshsq(kerr, nint, x, rpar(ipar(irwork)), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

      call colpnt(kcol, nint, ncpts, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol)),
     &            rpar(ipar(ixbs)))

      call iniy(t0, npde, kcol, nint, neq, ncpts, mflag(5),
     &          rpar(ipar(ixcol)), rpar(ipar(ixbs)),
     &          rpar(ipar(iabblk)), rpar(ipar(ibasi)), rpar(ipar(iy)),
     &          ipar(ipivot), rpar(ipar(irwork)), lenin, icflag,
     &          mflag(9), uinit, bndxa, difbxa, bndxb, difbxb)

      if (icflag .ne. 0) then
         goto 620
      endif

      call iniyp(t0, npde, kcol, nint, neq, ncpts,
     &           rpar(ipar(ixcol)), rpar(ipar(iabtop)),
     &           rpar(ipar(iabblk)), rpar(ipar(iabbot)),
     &           rpar(ipar(ibasi)), rpar(ipar(iy)), rpar(ipar(iyp)),
     &           ipar(ipivot), rpar(ipar(irwork)), lenin, icflag,
     &           mflag(9), f, bndxa, difbxa, bndxb, difbxb)

      if (icflag .ne. 0) then
         goto 620
      endif

      irshfg = 0

c     Set the parameters defining the ABD Jacobian (iteration) matrix.
c     This is inside the integer work array of dassl. Those parameter
c     include npde, kcol, nint.
      ipar(iiwork-1+17) = npde
      ipar(iiwork-1+18) = kcol
      ipar(iiwork-1+19) = nint

c     Set initial idid to be zero.
      idid = 0

c     Copy rtol and atol to be the relative and absolute error request
c     for the integration error.
      if (irshfg .eq. 0) then
         if (mflag(2) .eq. 0) then
            rpar(ipar(itatol)) = atol(1)
            rpar(ipar(itrtol)) = rtol(1)
         else
            do 130 i = 1, npde
               rpar(ipar(itatol)-1+i) = atol(i)
  130       continue
            do 140 i = 1, npde
               rpar(ipar(itrtol)-1+i) = rtol(i)
  140       continue
         endif
      endif

c     Second set of initializations removed.

c     Copy rpar(ipar(iy)) to rpar(ipar(iypre)).
      call dcopy(neq, rpar(ipar(iy)), 1, rpar(ipar(iypre)), 1)

      if (mflag(3) .eq. 1) then
         rpar(ipar(irwork)) = rpar(itstop)
      endif

c     Set the initial stepsize if applicable.
      if (mflag(6) .eq. 1) rpar(ipar(irwork)-1+3) = rpar(iiniss)

c     Set the maximum BDF order for DASSL if applicable.
      if (mflag(7) .eq. 1) ipar(iiwork-1+3) = ipar(imxord)

      goto 400

c-----------------------------------------------------------------------
c     When an adaptive mesh is used, this is not the first call for
c     the problem, and integration is to continue.
c-----------------------------------------------------------------------
  200 continue

c     Examine idid to determine if DASSL can be called again.
      if (idid .ne. 1) goto 790

c     If t0 must not go beyond tstop in the time stepping, then update
c     the DASSL floating point work array with the value of tstop.
      if (mflag(3) .eq. 1) then
         rpar(ipar(irwork)) = rpar(itstop)
         rpar(ipar(idasr))  = rpar(itstop)
      endif

c     Reset t0. If necessary, update rpar(ipar(iy)), rpar(ipar(iyp)).
      if (rpar(ipar(irwork)-1+4) .lt. tout) then
         t0 = rpar(ipar(irwork)-1+4)
         goto 400
      else
         t0 = tout
         call ddatrp(rpar(ipar(irwork)-1+4), tout, rpar(ipar(iy)),
     &               rpar(ipar(iyp)), neq, ipar(iiwork-1+8),
     &               rpar(ipar(irwork)+40+3*neq),
     &               rpar(ipar(irwork)-1+29))
         idid = 3
         goto 500
      endif

c-----------------------------------------------------------------------
c     Initialization after remeshing.
c-----------------------------------------------------------------------
  300 continue

      ipar(irmesh) = ipar(irmesh) + 1

      do 310 i = 1, 20
         ipar(iiwork-1+i) = ipar(idasi-1+i)
  310 continue

c     Tell DASSL to calculate the new iteration matrix and reset the
c     size of the iteration matrix.
      ipar(iiwork-1+5) = -1
      ipar(iiwork-1+16) = lenpd

c     Reset the ABD information in DASSL.
      ipar(iiwork-1+19) = nint

c     Move the last three sections of rpar (indexed by ipar(ixold),
c     ipar(idasr) and ipar(iypre)) to their new positions after a
c     remeshing in such a way as to avoid overwriting what you just
c     copied. The new and old ranges may overlap, so dcopy should
c     not be used -- aliasing is illegal in Fortran.
      if (nint .lt. ninold) then
c        rpar is shrinking; copy from the front
         do i = 1, nintmx+1+40+6*neqpre
            rpar(ipar(ixold)-1+i) = rpar(irold-1+i)
         end do
      else if (nint .gt. ninold) then
c        rpar is growing; copy in reverse
         do i = nintmx+1+40+6*neqpre, 1, -1
            rpar(ipar(ixold)-1+i) = rpar(irold-1+i)
         end do
      endif

c     Reset the first 40 elements of the floating point work array in
c     DASSL.
      call dcopy(40, rpar(ipar(idasr)), 1, rpar(ipar(irwork)), 1)

      call meshsq(kerr, nint, x, rpar(ipar(irwork)+40), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

      if (irshfg .ne. 2) then
         nstep = ipar(idasi-1+8) + 1
      else
         nstep = 1
      endif

c     reinit uses rpar(ipar(irwork)+40+6*neq) as its work array and
c     outputs y to rpar(ipar(irwork)+40).
      call reinit(npde, kcol, kcol, nint, ninpre, ncpts, neq,
     &            neqpre, icount, ipar(idasi-1+11), nstep, x,
     &            rpar(ipar(ixold)), rpar(ipar(iypre)), 0,
     &            rpar(ipar(irwork)+40+6*neq), lenri, ipar(ipivot),
     &            rpar(ipar(ih)), rpar(ipar(ixbs)), rpar(ipar(ixcol)),
     &            rpar(ipar(ibasi)), rpar(ipar(irwork)+40),
     &            rpar(ipar(iabblk)), icflag)

      if (icflag .ne. 0) then
         goto 630
      endif

c     cold start.
      if (irshfg .eq. 2) then
         call dcopy(neq, rpar(ipar(irwork)+40), 1, rpar(ipar(iy)), 1)
         call iniyp(t0, npde, kcol, nint, neq, ncpts,
     &              rpar(ipar(ixcol)), rpar(ipar(iabtop)),
     &              rpar(ipar(iabblk)), rpar(ipar(iabbot)),
     &              rpar(ipar(ibasi)), rpar(ipar(iy)),
     &              rpar(ipar(iyp)), ipar(ipivot),
     &              rpar(ipar(irwork)), lenin, icflag,
     &              mflag(9), f, bndxa, difbxa, bndxb, difbxb)
         if (icflag .ne. 0) then
            goto 630
         endif
c        If t0 must not go beyond tstop in the time stepping, update
c        the DASSL floating point work array with the value of tstop.
         if (mflag(3) .eq. 1) rpar(ipar(irwork)) = rpar(itstop)
         goto 320
      endif

c     Second set of initializations removed.

c     Back up the (possibly interpolated) B-spline coefficients
c     generated by reinit for a warm start.
      call dcopy(ipar(idasi-1+8)*neq, rpar(ipar(irwork)+40),
     &           1, rpar(ipar(iyprer)), 1)

c     Do divided differences to convert the interpolated B-spline
c     coefficients to the form used by DASSL - several derivatives
c     of the solution at the current time rather than the values at
c     the current and previous times.
c    (For space efficiency, this subroutine makes use of the work
c     array that is normally used by ERREST.)
      call divdif(neq, ipar(idasi-1+8)+1, rpar(ipar(idasr)-1+29),
     &            rpar(ipar(iework)), rpar(ipar(irwork)+40))

  320 continue

c     Move the B-spline coefficients from where reinit put them
c     to where DASSL expects to find them.
c     Copy in reverse since the ranges likely overlap, and the new
c     position is ahead of the old one.
      do i = neq*(ipar(iiwork-1+8)+1), 1, -1
         rpar(ipar(irwork)+40+3*neq-1+i) = rpar(ipar(irwork)+40-1+i)
      end do

c-----------------------------------------------------------------------
c     Time integration loop for DASSL.
c-----------------------------------------------------------------------

  400 continue

      call ddassl(res, neq, t0, rpar(ipar(iy)), rpar(ipar(iyp)),
     &            tout, ipar(iinfo), rpar(ipar(itrtol)),
     &            rpar(ipar(itatol)), idid, rpar(ipar(irwork)), lenrw,
     &            ipar(iiwork), leniw, rpar, ipar, jac,
     &            f, derivf, bndxa, difbxa, bndxb, difbxb)

c-----------------------------------------------------------------------
c     Check for a successful time step and decide whether to continue
c     integration or to perform a remeshing.
c-----------------------------------------------------------------------

      if (idid .le. 0) goto 600

      call errest(kcol, nint, npde, neq, necpts, icount,
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)),
     &            rpar(ipar(ixbs)), rpar(ipar(iey)), istart,
     &            mflag(2), atol, rtol, lenerr, rpar(ipar(iework)),
     &            rpar(ipar(iebasi)), rpar(ipar(iecoef)), lencof,
     &            errrat, rpar(ipar(ierint)), rpar(ipar(iercom)),
     &            ieflag, x, rpar(ipar(ih)), mflag(8))

      if (ieflag .eq. 0) then

c        The current step is accepted.
         if (icount .ne. 0) then

            ipar(istalr) = 1

            ipar(irpstp) = ipar(iypre) + 6 * neq - 1
c           Check for a sufficiently large rpar floating point
c           work array.
            if (lrp .lt. ipar(irpstp)) goto 770

         else
            ipar(istalr) = ipar(istalr) + 1
         endif

c        Update the backup information.
         call sucstp(ipar(iiwork-1+11), ipar(iiwork-1+8)+1, icount,
     &               neq, ipar(iiwork), rpar(ipar(irwork)),
     &               rpar(ipar(iey)), rpar(ipar(iyprer)), ipar(idasi),
     &               rpar(ipar(idasr)), rpar(ipar(iypre)))
         icount = 0
         istart = 1
         irshfg = 0

c        Check whether the integration is done or not.
         if (mflag(4) .eq. 0) then
            if (t0 .lt. tout) then
               goto 400
            else
               goto 500
            endif
         else
            isstep = isstep + 1
            if ((t0 .lt. tout) .and. (isstep .lt. ipar(iinstp))) then
               goto 400
            else
               goto 500
            endif
         endif

      else

c        The current step is rejected.
         if (icount .eq. maxrsh) goto 610

c        For the first remeshing at the current step, save nintpre and
c        neqpre at the last successful step.
         if (icount .eq. 0) then
            ninpre = nint
            neqpre = neq
         endif

         ninold = nint
         irold = ipar(ixold)

         call remesh(istart, icount, nintmx, ninpre, ninold,
     &               errrat, rpar(ipar(ierint)), irshfg,
     &               rpar(ipar(ixold)), nint, kerr, x,
     &               rpar(ipar(iework)))

         if (istart .eq. 1) then

c           This is not the initial step.
            t0 = rpar(ipar(idasr)-1+4)

c           In the first step after a remeshing, we do not allow DASSL
c           to increase the step size.
            if (rpar(ipar(idasr)-1+3) .gt. rpar(ipar(idasr)-1+7)) then
               rpar(ipar(idasr)-1+3) = rpar(ipar(idasr)-1+7)
            endif

c           In the first step after a remeshing, we do not allow DASSL
c           to increase the order of BDF method.
            if (ipar(idasi-1+7) .gt. ipar(idasi-1+8)) then
               ipar(idasi-1+7) = ipar(idasi-1+8)
            endif

            if (irshfg .eq. 2) then
c              This is a cold start
               ipar(istblc) = ipar(istblc) + ipar(idasi-1+11)
               ipar(icolds) = ipar(icolds) + 1
               ipar(iinfo) = 0
            endif

         else

c           This is the initial step.
            t0 = torign
            ipar(iinfo) = 0

         endif
c        print *, 't0 = ', t0, ' nint = ', nint
c        print *, 'istep = ', ipar(irmesh)
c        print *, (x(i),i=1,nint+1)
         goto 100

      endif

c-----------------------------------------------------------------------
c     Successful return section.
c-----------------------------------------------------------------------
  500 continue

c     Retrieve the value of mflag(1).
      mflag(1) = ipar(iinfo)

c     Retrieve the output vector y from the rpar communication array.
      do 510 i = 1, neq
         y(i) = rpar(ipar(iy)-1+i)
  510 continue

c     Retrieve information on the time stepping from the ipar array.
      ipar(isteps) = ipar(istblc) + ipar(iiwork-1+11)

c     Retrieve the value of ipar(iinstp) when mflag(4) = 1.
      if (mflag(4) .eq. 1) ipar(iinstp) = isstep

      return

c-----------------------------------------------------------------------
c     Unsuccessful return section.
c-----------------------------------------------------------------------
  600 continue
      write(6,9999) 'ERROR: BACOLI runtime error in time stepping.'
      write(6,9999) '       An error code and message should have'
      write(6,9999) '       been issued by DASSL.'
      return
  610 continue
      if (istart .eq. 1) then
         write(6,9998) 'ERROR: BACOLI has remeshed ', maxrsh
     &                 , ' times at', ' t0 =', rpar(ipar(idasr)-1+4)
      else
         write(6,9998) 'ERROR: BACOLI has remeshed ', maxrsh
     &                 , ' times at', ' t0 =', torign
      endif
      idid = -41
      return
  620 continue
      write(6,9999) 'ERROR: A singular matrix arises at the initial'
      write(6,9999) '       step. '
      idid = -42
      return
  630 continue
      write(6,9999) 'ERROR: A singular matrix arises during remeshing'
      write(6,9997) '       at t0 =', rpar(ipar(idasr)-1+4)
      idid = -43
      return
  640 continue
      if (istart .eq. 1) then
         write(6,9998) 'ERROR: nint >', nintmx, ' at'
     &                 , ' t0 =', rpar(ipar(idasr)-1+4)
      else
         write(6,9998) 'ERROR: nint >', nintmx, ' at'
     &                 , ' t0 =', torign
      endif
      idid = -44
      return

c-----------------------------------------------------------------------
c     The following section is the return point for invalid input.
c-----------------------------------------------------------------------

  710 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require:  0 <= mflag(i) <= 1, i = 1, 2, ..., 8.'
      idid = -51
      return
  711 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require:  0 <= mflag(9) <= 3'
      idid = -51
      return
  715 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require:  if mflag(4) = 1, ipar(8) must be set to'
      write(6,9999) 'be a positive integer.'
      idid = -52
      return
  720 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require:  if mflag(6) = 1, tout must be in front'
      write(6,9999) 'of t0.'
      idid = -53
      return
  725 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require:  if mflag(6) = 1, rpar(2) must be the'
      write(6,9999) 'initial stepsize, thus nonzero.'
      idid = -54
      return
  730 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require: npde > 0.'
      idid = -55
      return
  740 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require: 2 < kcol <=', mxkcol
      idid = -56
      return
  760 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require: x(1) < x(2) < ... < x(nint+1).'
      idid = -58
      return
  770 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require: lrp >= ', ipar(irpstp), '.'
      idid = -59
      return
  780 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'Require: lip >= ', ipar(iipstp), '.'
      idid = -60
      return
  790 continue
      write(6,9999) 'ERROR: BACOLI input violation.'
      write(6,9999) 'IDID .ne. 1, on a continuation call of BACOLI'
      write(6,9999) 'If IDID > 1, set idid = 1 and tout (t0 < tout)'
      write(6,9999) 'If IDID < -1, the code cannot be continued due to'
      write(6,9999) '              a previous error.'
      idid = -61
      return

c-----------------------------------------------------------------------
 9997 format(a,e12.5)
 9998 format(a,i4,a,a,e12.5)
 9999 format(a,i10,a,i4,a,i4,a,i4,a,i4)
c-----------------------------------------------------------------------
      end

      subroutine values(kcol, xsol, nint, x, npde, npts, nderiv,
     &                  usol, y, work)

c-----------------------------------------------------------------------
c Purpose:
c     This routine computes the solution u and the first nderv
c     derivatives of u at the npts points xsol. It then returns the
c     values in the array usol.
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------

c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c       kcol is the number of collocation points to be used in
c       each subinterval.
c
        integer                 npts
c       npts is the number of points in the x vector.
c
        double precision        xsol(npts)
c       xsol is an arbitrary set of spatial points at which the solution
c       and the first nderv derivative values are to be calculated.
c
        integer                 nint
c       nint >= 1 is the number of subintervals defined by the spatial
c       mesh x.
c
        double precision        x(nint+1)
c       x is the spatial mesh which divides the interval [x_a,x_b] into
c       x_a = x(1) < x(2) < x(3) < ... < x(nint+1) = x_b.
c
        integer                 npde
c       npde >= 1 is the number of components in the system of PDEs.
c
        integer                 nderiv
c       nderiv is the number of derivatives of the solution which are
c       to be calculated.
c
        double precision        y(npde*(nint*kcol+nconti))
c       y is the vector of B-spline coefficients at the output time.
c
c       output:
        double precision        usol(npde, npts, nderiv+1)
c       usol is the solution and the spatial partial derivatives (up to
c       the nderiv-th derivative) at xsol and at BACOLI's output time.
c
c       Work Storage:
        double precision        work((kcol+nconti)*(nderiv+1)
     *                                 +kcol*(nint+1)+2*nconti)
c       work is a floating point work storage array.
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 mflag
c                               mflag is required by subroutine
c                               interv.
c
        integer                 ilo
c                               ilo is required by subroutine
c                               interv.
c       Pointers into the floating point work array:
        integer                 ixbs
c                               work(ixbs) contains the breakpoint
c                               sequence.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 ii
        integer                 mj
        integer                 mm
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               interv
c
c-----------------------------------------------------------------------

c     set up the value for ileft, mflag and ncpts.
      ileft = 0
      mflag = -2
      ncpts = nint * kcol + nconti

c     set the pointer into the floating point work array
      ixbs  = (kcol+nconti)*(nderiv+1) + 1

c     Store the piecewise polynomial space breakpoint sequence in
c     work(ixbs).
c
      do 10 i = 1, kcol + nconti
         work(ixbs-1+i) = x(1)
         work(ixbs-1+i+ncpts) = x(nint+1)
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
   20 continue
   30 continue

      do 70 i = 1, npts
c
c     interv is called to compute ileft. bsplvd is called to compute
c     the values of the basis function at the required point.
         call interv(work(ixbs), ncpts, xsol(i), ileft, mflag, ilo)
         call bsplvd(work(ixbs),kcol+nconti,xsol(i),ileft,work,
     &               nderiv+1)
         ii = ileft - kcol - nconti
         do 60 j = 1, nderiv + 1
            do 50 k = 1, npde
               usol(k,i,j) = zero
               do 40 m = 1, kcol + nconti
                  mm = (m + ii - 1) * npde
                  mj = (j - 1) * (kcol + nconti) + m
                  usol(k,i,j) = usol(k,i,j) + y(mm+k) * work(mj)
   40          continue
   50       continue
   60    continue
   70 continue
      return
      end

      subroutine errest(kcol, nint, npde, neq, npts, icount,
     &                  xsol, wts, xbs, y, istart, mflag2, atol,
     &                  rtol, lenwk, work, errbas, ercoef, lencof,
     &                  errrat, errint, errcom, ieflag, x, h, est)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the error estimate for each subinterval
c       and for each component of the PDE system, and decides whether a
c       remeshing is necessary or not. If a remeshing is deemed
c       necessary, the distribution of mesh points is determined by this
c       error estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Jack Pew, May 17, 2013.
c
c Explanation of Modifications: This version of errest has replaced the
c evaluation of a second collocation solution with interpolated values
c from the (now only) computed collocation solution. For a choice of
c est, it will either use local extrapolation (the LOI scheme) or
c superconvergent interpolation (the SCI scheme) as the replacement.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 nintsm
        parameter              (nintsm = 15)
c                               when the current step is the first step
c                               after remeshing, we require
c                               if nint <= nintsm
c                                    errrat < saffa2
c                               else
c                                    saffa1 < errrat < saffa2.
c                               endif
c
        double precision        zero
        parameter              (zero = 0.0d0)
c
        double precision        one
        parameter              (one = 1.0d0)
c
        double precision        two
        parameter              (two = 2.0d0)
c
        double precision        saffa1
        parameter              (saffa1 = 0.1d0)
c
        double precision        saffa2
        parameter              (saffa2 = 0.4d0)
c                               These safety factors are heuristic
c                               values that bound the acceptable global
c                               scaled error estimate after a remeshing.
c                               They are used in Spatial Error Test II.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval in the
c                               computed solution.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(nint*kcol+nconti) is the
c                               number of bspline coefficients (or
c                               DAEs) when using dassl_kcol.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector, which is equal to
c                               nint*quad.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xsol(npts)
c                               xsol is the npts Gauss-Legendre
c                               points at which the solution are
c                               to be calculated for the L2-norm
c                               of the error.
c
        double precision        wts(npts)
c                               wts is the npts Gauss-Legendre
c                               weights at the corresponding xsol.
c
        double precision        x(nint+1)
c                               The current mesh.
c
        double precision        h(nint)
c                               The mesh step size sequence.
c
        double precision        xbs((nint+1)*kcol+2*nconti)
c                               xbs is the breakpoint sequence when
c                               using dassl_kcol.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients when using dassl_kcol.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 mflag2
c                               mflag2 = 0, scalar atol and rtol.;
c                               mflag2 = 1, vector atol and rtol.
c
        integer                 est
c                               est=0 => use LOI error estimate
c                               est=1 => use SCI error estimate

        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        integer                 lenwk
c                               lenwk is the size of the work storage
c                               array and must satisfy:
c                               lenwk >= 2*npde*nint*(kcol+2)
c                                        +npde*nint
c                                        +(kcol-3)*nint*npde
c                                        +2*(nint+1)*npde
c                               if using the LOI scheme, or:
c                               lenwk >= 2*npde*nint*(kcol+3)
c                                        +npde*nint
c                                        +(kcol-2)*nint*npde
c                                        +2*(nint+1)*npde
c                               if using the SCI scheme.
c
        integer                 lencof
c                               lencof is the size of the coefficient
c                               storage array and must satisfy:
c                               lencof >= (kcol+1)*(kcol+2)
c                                         +(kcol+nconti)*(kcol-3)*nint
c                                         +(kcol+nconti)*2*(nint+1)
c                               if using the LOI scheme, or:
c                               lencof >= (kcol+4)*(kcol+3)*nint
c                                         +(kcol+nconti)*(kcol-2)*nint
c                                         +(kcol+nconti)*2*(nint+1)
c                               if using the SCI scheme.
c
c       Work Storage:
        double precision        work(lenwk)
c                               work is a floating point work storage
c                               array of size lenwk.
c
c       output:
        double precision        errbas((kcol+nconti)*npts)
c                               errbas holds the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol.
c                               They are reusable between remeshings.
c
        double precision        ercoef(lencof)
c                               ercoef holds the values of the nonzero
c                               B-spline basis functions and the
c                               Hermite-Birkhoff coefficients used in
c                               the selected interpolation subroutine.
c                               They are reusable provide the spatial
c                               mesh does not change.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of errcom.
c
        double precision        errint(nint)
c                               errint is the error estimate at each
c                               subinterval.
c
        double precision        errcom(npde)
c                               errcom is the error estimate for
c                               each component of pdes at the whole
c                               range, i.e. from x_a to x_b.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
C-----------------------------------------------------------------------
c Local Variables:
        double precision        errsum
c                               errsum is the sum of errint.
c
        double precision        errmax
c                               errmax is the maximum value of
c                               errint(i), i = 1, nint.
c
        double precision        aerr
c                               aerr is the average value of errint(i),
c                               i = 1, nint.
c
        double precision        disind
c                               disind is equal to errmax/aerr, and it
c                               indicates the error distribution over
c                               the mesh.
c
        integer                 quad
c                               quad is the number of Gaussian
c                               quadrature points used in the L2-norm.
c
        double precision        power
c                               The exponent applied to each errint(i)
c                               component. This is one/dble(kcol+2) for
c                               the SCI and one/dble(kcol+1) for the LOI
c                               scheme.
c
        double precision        rL, rR
c                               Mesh subinterval size ratios
c                              (for scaling the SCI estimate.)
c
c       Pointers into the floating point work arrays:
        integer                 iusol1
c                               work(iusol1) stores the evaluations of
c                               the collocation solution at the npts
c                               Gaussian quadrature points.
c
        integer                 iusol2
c                               work(iusol2) stores the evaluations of
c                               interpolant at the npts Gaussian
c                               quadrature points.
c                               Under the LOI scheme, they are one
c                               order of accuracy lower than those in
c                               work(iusol1), and under the SCI scheme
c                               they are one order higher.
c
        integer                 ierrci
c                               work(ierrci) stores the error estimate
c                               for each subinterval for each component.
c
        integer                 iliu
c                               work(iliu) stores evaluations of the
c                               computed solution at the superconvergent
c                               points internal to each subinterval.
c
        integer                 ilium
c                               work(ilium) stores evaluations of the
c                               computed solution at the mesh points.
c
        integer                 ih
c                               ercoef(ih) stores the H coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ihd
c                               ercoef(ihd) stores the Hbar coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ig
c                               ercoef(ig) stores the G coefficients
c                               for the Hermite-Birkhoff interpolant.
c
        integer                 ibassc
c                               ercoef(ibassc) stores the evaluations
c                               of the B-spline basis functions at
c                               the superconvergent points internal
c                               to each subinterval.
c
        integer                 ibasm
c                               ercoef(ibasm) stores the evaluations of
c                               the B-spline basis functions at the
c                               mesh points. It and ercoef(ibassc) are
c                               analogs to the errbas vector, which
c                               stores basis function evaluations at
c                               quadrature points for the subroutine
c                               errval.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j, m, ij, im, mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               errval
c                               scint
c                               lowint
c
c-----------------------------------------------------------------------

c     quad is the number of quadrature points used per subinterval
      if (est .eq. 0) then
         quad = kcol + 2
      elseif (est .eq. 1) then
         quad = kcol + 3
      endif

c     Set pointers into two work arrays known to errest, work and ercoef
c     The SCI scheme requires more storage here, as its H-B coefficients
c     are different from one subinterval to the next due to the
c     dependence on adjacent subinterval sizes.
c     Also, the SCI uses a greater number of quadrature points.
      iusol1 = 1
      iusol2 = iusol1 + npde * nint * quad
      ierrci = iusol2 + npde * nint * quad
      iliu   = ierrci + npde * nint
      if (est .eq. 0) then
         ilium  = iliu   + npde * nint * (kcol - 3)
         ih     = 1
         ihd    = ih     + 2 * quad
         ig     = ihd    + 2 * quad
         ibassc = ig     + (kcol - 3) * quad
         ibasm  = ibassc + (kcol + nconti) * (kcol - 3) * nint
      elseif (est .eq. 1) then
         ilium  = iliu   + npde * nint * (kcol - 2)
         ih     = 1
         ihd    = ih     + 2 * quad * nint
         ig     = ihd    + 2 * quad * nint
         ibassc = ig     + kcol * quad * nint
         ibasm  = ibassc + (kcol + nconti) * (kcol - 2) * nint
      endif

c-----------------------------------------------------------------------
c     Generate the values of the collocation solution and the
c     interpolant at the Gaussian quadrature points, stored in xsol, and
c     save them in work(iusol1) and work(iusol2), respectively.
c     The LOI and the SCI schemes use quadrature rules of different
c     degrees.
c     errval evaluates the computed solution, while lowint and scint
c     evaluate the interpolant.

      call errval(kcol, nint, npde, neq, quad, istart, icount,
     &            xbs, xsol, y, errbas, work(iusol1))

      if (est .eq. 0) then
         call lowint(kcol, nint, npde, istart, icount, xbs, x, y,
     &         ercoef(ih), ercoef(ihd), ercoef(ig), ercoef(ibassc),
     &         ercoef(ibasm), work(iliu), work(ilium), work(iusol2))
      elseif (est .eq. 1) then
         call scint(kcol, nint, npde, istart, icount, xbs, x, y,
     &         ercoef(ih), ercoef(ihd), ercoef(ig), ercoef(ibassc),
     &         ercoef(ibasm), work(iliu), work(ilium), work(iusol2))
      endif

c-----------------------------------------------------------------------
c     Initialization task.
      do 10 i = 1, nint
         errint(i) = zero
   10 continue

      do 20 i = 1, npde
         errcom(i) = zero
   20 continue

      do 30 i = 1, npde * nint
         work(ierrci - 1 + i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Calculate the error estimate at each subinterval for each
c     component of PDEs.

      if (mflag2 .eq. 0) then
c        Use scalar error tolerance.
         do 60 m = 1, npde
            do 50 i = 1, nint
               do 40 j = 1, quad
                  ij = (i - 1) * quad + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im)
     &                     + ((work(iusol1-1+mm) - work(iusol2-1+mm))
     &                     / (atol(1) + rtol(1)*abs(work(iusol1-1+mm))))
     &                     **2 * wts(ij)
   40          continue
   50       continue
   60    continue
      else
c        Use vector error tolerance (to weight PDEs in the system.)
         do 90 m = 1, npde
            do 80 i = 1, nint
               do 70 j = 1, quad
                  ij = (i - 1) * quad + j
                  mm = npde * (ij - 1) + m
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im)
     &                     + ((work(iusol1-1+mm) - work(iusol2-1+mm))
     &                     / (atol(m) + rtol(m)*abs(work(iusol1-1+mm))))
     &                     **2 * wts(ij)
   70          continue
   80       continue
   90    continue
      endif

c-----------------------------------------------------------------------
c     Scale work(ierrci) to `correct' for overestimations by the
c     SCI due to mesh ratios. We divide by the product of the ratios
c     instead of the only the larger one because for internal
c     subintervals, the underestimation of errors in layer regions is an
c     issue. (This is also the case in the original BACOLI.)
      if (EST .eq. 1) then
         do j = 1, npde
c           i = 1
            rL = zero
            rR = h(2) / h(1)
            ij = ierrci + (j - 1) * nint
            work(ij) = work(ij) * rR * rR
            do i = 2, nint-1
               rL = h(i-1) / h(i)
               rR = h(i+1) / h(i)
               ij = ij + 1
               work(ij) = work(ij) * rL * rR
            end do
c           i = nint
            rL = h(nint-1) / h(nint)
            rR = zero
            ij = ij + 1
            work(ij) = work(ij) * rL * rL
         end do
      end if

c     Calculate errint and errcom.
      do 110 j = 1, npde
         do 100 i = 1, nint
            ij = ierrci - 1 + (j - 1) * nint + i
            errint(i) = errint(i) + work(ij)
            errcom(j) = errcom(j) + work(ij)
  100    continue
  110 continue

c     When using the LOI scheme, error is actually estimated for a
c     solution that is one order lower than that which is computed.
      if (est .eq. 0) power = one/dble(kcol+1)
      if (est .eq. 1) power = one/dble(kcol+2)

c     Take the square root and update errint and errcom.
      do 120 i = 1, nint
         errint(i) = sqrt(errint(i))
         errint(i) = errint(i) ** power
  120 continue

      do 130 i = 1, npde
         errcom(i) = sqrt(errcom(i))
  130 continue

c-----------------------------------------------------------------------
c     Decide whether remeshing is needed.
      ieflag = 0

c     update errrat (the max of errcom.)
      errrat = zero
      do 140 i = 1, npde
         if (errcom(i) .gt. errrat) then
            errrat = errcom(i)
         endif
  140 continue

c     Calculate errsum to be the sum of the errint. Find the maximum
c     errint(i) and save it in errmax.
      errsum = errint(1)
      errmax = errint(1)
      do 150 i = 2, nint
         if (errmax .lt. errint(i)) errmax = errint(i)
         errsum = errint(i) + errsum
  150 continue

c     Let aerr be the mean value of errint(i).
      aerr = errsum/dble(nint)

c     Calculate disind (which is a measure of error distribution.)
      disind = errmax/aerr

c     This is Spatial Error Test I.
      if (disind .gt. two) then
c        Mesh not well distributed. Remeshing needed.
         ieflag = 1
      else
c        Passed Test I, now do Test II.
         if ((istart .eq. 1) .and. (icount .eq. 0)) then
            if (errrat .ge. one) ieflag = 1
         else
            if (nint .gt. nintsm) then
               if ((errrat .ge. saffa2) .or. (errrat .le. saffa1))
     &            ieflag = 1
            else
               if (errrat .ge. saffa2) ieflag = 1
            endif
         endif
      endif

      return
      end

      subroutine errval(kcol, nint, npde, neq, nptse, istart, icount,
     &                  xbs, xsol, y, errbas, usol)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the values of the (kcol+nconti) nonzero
c       bspline basis function at each Gaussian point of xsol.
c       Then determine the solution usol, which is used for error
c       estimate, at xsol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 29, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(kcol*nint+2) is the number of
c                               bspline coefficients.
c
        integer                 nptse
c                               nptse is the number of Gaussian points
c                               in each subinterval for the error
c                               estimate.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xbs((nint+1)*kcol+2*nconti)
c                               The breakpoint sequence.
c
        double precision        xsol(nptse*nint)
c                               xsol is a set of spatial points at which
c                               the solution are to be calculated for
c                               error estimate.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients.
c
c       output:
        double precision        errbas(kcol+nconti, nptse*nint)
c                               errbas is the values of the nonzero
c                               basis functions at xsol.
c
        double precision        usol(npde, nptse*nint)
c                               uval is the solution at xsol.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c
c-----------------------------------------------------------------------

c     check whether errbas is necessary to be calculated.
      if ((istart .eq. 1) .and. (icount .eq. 0)) goto 30

c     calculate errbas.
      do 20 i = 1, nint
         ileft = kcol + nconti + (i - 1) * kcol
         do 10 j = 1, nptse
            jj = (i - 1) * nptse + j
            call bsplvd(xbs, kcol+nconti, xsol(jj), ileft, errbas(1,jj),
     &                  1)
   10    continue
   20 continue

   30 continue

c     compute the values of usol at xsol.
      do 70 i = 1, nint
         do 60 j = 1, nptse
            jj = (i - 1) * nptse + j
            do 50 k = 1, npde
               usol(k,jj) = zero
               do 40 m = 1, kcol + nconti
                  mm = npde * (m + (i - 1) * kcol - 1) + k
                  usol(k,jj) = usol(k,jj) + y(mm) * errbas(m,jj)
   40          continue
   50       continue
   60    continue
   70 continue

      return
      end

      subroutine sucstp(istep, nstep, icount, neq, icdas, cdasr, cypre,
     &                  cyprer, idas, dasr, ypre)

c-----------------------------------------------------------------------
c Purpose:
c       This routine stores the necessary information after each
c       accepted time step (i.e. no need for remeshing). This info
c       is needed if a remeshing is required next step.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 17, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Subroutine Parameters:
c       input:
        integer                 istep
c                               istep is the number of time steps that
c                               DASSL has taken when using dassl.
c
        integer                 nstep
c                               nstep is the number of previous steps
c                               necessary.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 neq
c                               neq is the number of bspline
c                               coefficients when using dassl_kcol.
c
        integer                 icdas(20)
c                               icdas stores the first 20 elements of
c                               the integer work array in dassl.
c
        double precision        cdasr(40)
c                               cdasr stores the first 40 elements of
c                               the floating point work array in dassl.
c
        double precision        cypre(neq)
c                               cypre is the vector of bspline
c                               coefficients at the current step.
c
        double precision        cyprer(6*neq)
c                               cyprer is the vector of bspline
c                               coefficients at the previous steps.
c
c       output:
        integer                 idas(20)
c                               idas is a copy of icdas.
c
        double precision        dasr(40)
c                               dasr is a copy of cdasr.
c
        double precision        ypre(6*neq)
c                               ypre1 stores the bspline coefficients
c                               at the past 6 steps of dassl_kcol.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 imod
        integer                 itemp
c
c-----------------------------------------------------------------------
c Loop Indices:
        integer                 i
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c Fortran Functions used:
c                               mod
c
c-----------------------------------------------------------------------

      do 10 i = 1, 20
         idas(i) = icdas(i)
   10 continue

      call dcopy(40, cdasr, 1, dasr, 1)

      imod = mod(istep, 6)
      if (imod .eq. 0) then
         imod = 6
      endif
      imod = 7 - imod

      call dcopy(neq, cypre, 1, ypre((imod-1)*neq+1), 1)

      if (icount .eq. 0) goto 99

c-----------------------------------------------------------------------
c     This modulus operations in this section are to decide which 'slot'
c     of the solution history section to put the current step's solution
c     information. Slots are cycled through, rather than putting the
c     current step's solution at the front after every step.
c
c     The following statement setting itemp = nstep is necessary to
c     avoid confusion for some compilers.
      itemp = nstep
      if (itemp .eq. 6) then
         itemp = 5
      endif

      if (imod .eq. 6) then
         call dcopy(itemp*neq, cyprer, 1, ypre, 1)
      else
         if ((imod+itemp) .le. 6) then
            call dcopy(itemp*neq, cyprer, 1, ypre(imod*neq+1), 1)
         else
            call dcopy((6-imod)*neq, cyprer, 1, ypre(imod*neq+1), 1)
            call dcopy((itemp+imod-6)*neq, cyprer((6-imod)*neq+1), 1,
     &                 ypre, 1)
         endif
      endif

c-----------------------------------------------------------------------
   99 continue

      return
      end

      subroutine remesh(istart, icount, nintmx, ninpre, ninold, errrat,
     &                  errint, irshfg, xold, nint, kcol, x, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates a new mesh by equidistributing the error
c       in each subinterval.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 22, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        double precision        point5
        parameter              (point5 = 0.5d0)
c
        double precision        one
        parameter              (one    = 1.0d0)
c
        double precision        two
        parameter              (two    = 2.0d0)
c
        double precision        saffac
        parameter              (saffac = 0.2d0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 nintmx
c                               the maximal number of subintervals that
c                               the user requires.
c
        integer                 ninpre
c                               ninpre is the number of subintervals
c                               when icount = 0 before remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        errint(ninold)
c                               errint is the error estimate at
c                               each subintervals.
c
c       Output:
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial call or continuation
c                                           calls;
c                                      = 1, remesh with a hot start.
c                                      = 2, remesh with a cold start.
c
        double precision        xold(ninpre+1)
c                               xold is the spatial mesh when icount = 0
c                               before remeshing.
c
c       In-output:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               ninmx >= nint >= 1.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh. As input, it is
c                               the value before remeshing; as output,
c                               it is the value after remeshing.
c
c       Work storage:
        double precision        work(2*ninold+1)
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        aerr
        double precision        berr
c
c       Pointers into the floating point work array:
        integer                 ierror
c                               work(ierror-1+i) is the L2-norm error
c                               estimate at the first i subintervals.
c
        integer                 ixold
c                               work(ixold) contains a copy of mesh
c                               points before remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
c
c-----------------------------------------------------------------------
c Functions used:
c                               dble
c                               int
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      ierror = 1
      ixold = ierror + ninold

c-----------------------------------------------------------------------
c     Update xold.
      if (icount .eq. 0) then
         do 10 i = 1, ninpre + 1
            xold(i) = x(i)
   10    continue
      endif

c-----------------------------------------------------------------------
c     Update icount, irshfg and nint.
      icount = icount + 1

c     If this is the first remesh at the current step which is not the
c     initial step.
      if ((icount .eq. 1) .and. (istart .eq. 1)) then
         irshfg = 1
         goto 20
      endif

c     If after four hot start the code still can not satisfy the error
c     requirement, a cold start will take place.
      If ((icount .eq. 5) .and. (istart .eq. 1)) then
         irshfg = 2
         nint = ninpre
         goto 20
      endif

c     Update errrat.
      errrat = (errrat/saffac) ** (one/dble(kcol+2))

c     Set the upper bound and lower bound of the ratio of nint over
c     ninold.
      if (errrat .gt. two) then
         errrat = two
      else
         if (errrat .lt. point5) then
            errrat = point5
         endif
      endif

      nint = int(ninold * errrat)

c     The code does not allow nint = ninold.
      if (nint .eq. ninold) then
         nint = nint + 1
      endif

c     Stop now if nint > nintmx and let bacoli pick up the error
c     when the code loops back up to the 100 label.
      if (nint .gt. nintmx) return

   20 continue

c-----------------------------------------------------------------------
c     Update work(ixold) to be the mesh before remeshing.
      do 30 i = 1, ninold + 1
         work(ixold-1+i) = x(i)
   30 continue

c-----------------------------------------------------------------------
c     Store work(i) to be the sum of the error at the first i
c     subintervals.
      work(ierror) = errint(1)
      do 40 i = ierror-1+2, ninold
         work(i) = errint(i) + work(i-1)
   40 continue

c     Let aerr to be the mean value of errint(i).
      aerr = work(ninold)/dble(nint)

c     Equidistribute the mesh points.
      berr = aerr
      j = 1

      do 60 i = 2, nint
   50    continue
         if (berr .gt. work(j)) then
            j = j + 1
            goto 50
         else
            if (j .eq. 1) then
               x(i) = work(ixold) + (work(ixold-1+2) - work(ixold))
     &                * berr/work(1)
            else
               x(i) = work(ixold-1+j) + (work(ixold-1+j+1) -
     &                work(ixold-1+j)) * (berr - work(j-1))/errint(j)
            endif
         endif
         berr = berr + aerr
   60 continue

      x(1) = work(ixold)
      x(nint+1) = work(ixold-1+ninold+1)

      return
      end

      subroutine meshsq(kcol, nint, x, work, h, excol, ewts)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mesh size sequence, then generates
c       the collocation points and Gaussian weights for error
c       estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 5, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
c       Work Storage:
        double precision        work((kcol+3)*(kcol+3))
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        double precision        excol(nint*(kcol+3))
c                               excol is the Gaussian quadrature point
c                               sequence which is used for the L2-norm
c                               error estimate.
c
        double precision        ewts(nint*(kcol+3))
c                               ewts is the Gaussian weight sequence
c                               which is used for error estimate.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        rho(mxkcol+3)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+3)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Calculate the mesh step size sequence.
      do 10 i = 1, nint
         h(i) = x(i+1)-x(i)
   10 continue

c     Compute the Gaussian points and Gaussian weights.
      call gauleg(kcol+3, (kcol+3)*(kcol+3), rho, wts,
     &            work, 4)

c     Define the Gaussian quadrature point sequence.
      do 30 i = 1, nint
         ii = (i - 1) * (kcol + 3)
         do 20 j = 1, kcol+3
            excol(ii + j) = x(i) + h(i) * rho(j)
            ewts(ii + j) = h(i) * wts(j)
   20    continue
   30 continue

      return
      end

      subroutine colpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the piecewise polynomial space breakpoint
c       sequence, and calculates the collocation point sequence.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 3, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
c       Work Storage:
        double precision        work(kcol*kcol)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [a,b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i),
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        rho(mxkcol+1)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+1)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Generate the piecewise polynomial space breakpoint sequence.
      do 10 i = 1, kcol + nconti
         xbs(i) = x(1)
         xbs(i + ncpts) = x(nint + 1)
   10 continue
      do 30 i = 2, nint
         ii = (i - 2) * kcol + kcol + nconti
         do 20 j = 1, kcol
            xbs(ii + j) = x(i)
   20    continue
   30 continue

c-----------------------------------------------------------------------
c     Compute the Gaussian points.
      call gauleg(kcol, kcol*kcol, rho, wts, work, 2)

c     Define the collocation point sequence.
      xcol(1) = x(1)
      do 50 i = 1, nint
         ii = (i - 1) * kcol + 1
         do 40 j = 1, kcol
            xcol(ii + j) = x(i) + h(i) * rho(j)
   40    continue
   50 continue
      xcol(ncpts) = x(nint + 1)

      return
      end

      subroutine iniy(t0, npde, kcol, nint, neq, ncpts, ifglin, xcol,
     &                xbs, abdblk, fbasis, y, ipivot, work, lw, icflag,
     &                ifgfdj, uinit, bndxa, difbxa, bndxb, difbxb)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs initialization tasks required on the
c       initial step, which are:
c
c               calculating the B-spline basis functions,
c               constructing abdblk of the collocation matrices, and
c               determining y(t0).
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, November 8, 2001.
c Last modified by Jack Pew, July 11, 2013.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        negone
        parameter              (negone = -1.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               t0 is the initial time.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bsplines
c                               coefficients (or DAEs).
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ifglin
c                               ifglin is a flag for the boundary
c                               conditions.
c                               ifglin = 1, if both boundary conditions
c                                           are Dirichlet;
c                                      = 0, otherwise.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i),
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
        integer                 lw
c                               lw is the size of the work storage
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +2*neq+2*npde+2*npde*npde
c
        integer                 ifgfdj
c                               Are finite differences being used to
c                               approximate boundary function Jacobi?
c                               Set to 0 or 1 if so,
c                               or to 2 or 3 to use user-provided
c                               analytic partial derivatives.
c
        external                uinit
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
c                               refer to the preamble of BACOLI.
c
c       Work Storage:
        integer                 ipivot(neq)
c                               pivoting information from the
c                               factorization of the temporary matrix.
c
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. fbasis(k,j,i) contains the
c                               values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        y(neq)
c                               y = y(t0) is the vector of B-spline
c                               coefficients at the initial time.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint index information.
c
        integer                 nels
c                               the number of elements in one
c                               collocation block of work.
c
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the
c                               bottom block which is required since
c                               crdcmp overwrites the input collocation
c                               matrix.
c
        integer                 idelta
c                               work(idelta) contains the residual which
c                               indicates how well y satisfies to the
c                               boundary condition and the initial
c                               condition at the internal collocation
c                               points.
c
        integer                 ivcol
c                               work(ivcol) contains the values of u
c                               at the internal collocation points.
c
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde). That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
        integer                 ifdwrk
c                               work(ifdwrk-1+i), i=1, 2*npde, is used
c                               as work storage for fdbndx.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 jj
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               bsplvd
c                               difbxa
c                               difbxb
c                               eval
c                               uinit
c                               crdcmp
c                               crslve
c                               fdbndx
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c                               dscal
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      idelta = iabdbt + npde*npde*nconti
      ivcol  = idelta + neq
      iu     = ivcol  + neq-2*npde
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde
      ifdwrk = idbdt  + npde

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 20 i = 1, npde*npde*nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue
      do 30 i = 1, nint*nels
         abdblk(i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j)
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Uinit is called to evaluate the first npde components at the
c     left boundary point, and save in y.
      call uinit(xcol(1), y(1), npde)

c     Making use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 40 i = 1, npde
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   40 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 80 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 70 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is the position in the y vector where the values for the
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
            jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bsplvd(xbs,kcol+nconti,xcol(ii),ileft,fbasis(1,1,ii),3)
            call uinit(xcol(ii), y(jj), npde)

            do 60 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 50 m = 1, npde
                  mm = ll + (m -1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   50          continue
   60       continue
   70    continue
   80 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      ii = neq - npde + 1
      call uinit(xcol(ncpts), y(ii), npde)
      do 90 i = 1, npde
         ii = ((i-1)+npde)*npde + i
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   90 continue

c-----------------------------------------------------------------------
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)

c     Check whether both boundary conditions are derichlet boundary
c     conditions. If no, copy the values at the internal collocation
c     points to work(ivcol), which will be used for newton iterations.
      if (ifglin .eq. 0) then
         call dcopy(neq-2*npde,y(npde+1),1,work(ivcol),1)
         call dscal(neq-2*npde,negone,work(ivcol),1)
      endif

c-----------------------------------------------------------------------
c     Generate the initial vector y(t0).
c-----------------------------------------------------------------------

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999

c     Solve the linear system. If Dirichlet boundary conditions are
c     given, this gives the basis function coefficients for the initial
c     conditions, i.e. y(t0). If not, this gives the predictor of y(t0).
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,y,0)

      if (icflag .ne. 0) goto 999

c     Check whether both boundary conditions are derichlet boundary
c     conditions.
      if (ifglin .eq. 1) goto 999

c-----------------------------------------------------------------------
c     Newton iteration loop.

c     Calculate (work(idelta-1+i), i = npde+1, neq-npde), which depends
c     on the nint blocks in the middle of the collocation matrix A.
      call dcopy(neq-2*npde,work(ivcol),1,work(idelta+npde),1)
      do 130 i = 1, nint
         do 120 j = 1, kcol + nconti
            do 110 l = 1, kcol
               ll = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(l-1)*npde
               do 100 m = 1, npde
                  ii = idelta-1+npde+(i-1)*npde*kcol+(l-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  work(ii) = work(ii) + abdblk(ll) * y(mm)
  100          continue
  110       continue
  120    continue
  130 continue

c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)

c     Update the values at the left boundary.
      call eval(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y)
      call bndxa(t0, work(iu), work(iux), work(idelta), npde)
      if (ifgfdj .lt. 2) then
         call fdbndx(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxa, work(ifdwrk))
      else
         call difbxa(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Set up the top block and save in work(iabdtp).
      do 150 j = 1, npde
         do 140 i = 1, npde
            ii = iabdtp - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            work(ii) = work(idbdu-1+mm) - work(jj)
  140    continue
  150 continue

c     Update the values at the right boundary.
      call eval(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y)
      call bndxb(t0, work(iu), work(iux), work(idelta+neq-npde), npde)
      if (ifgfdj .lt. 2) then
         call fdbndx(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxb, work(ifdwrk))
      else
         call difbxb(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Set up the bottom block and save in work(iabdbt).
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = iabdbt - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            work(jj) = work(idbdu-1+mm) - work(ii)
  160    continue
  170 continue

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999

c     Solve the corrector equation.
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            work(idelta),0)

      if (icflag .ne. 0) goto 999

c     Now apply the corrector of y(t0).
      do 180 i = 1, neq
         y(i) = y(i) - work(idelta-1+i)
  180 continue

c-----------------------------------------------------------------------

  999 return
      end

      subroutine iniyp(t0, npde, kcol, nint, neq, ncpts, xcol,
     &                 abdtop, abdblk, abdbot, fbasis, y, yprime,
     &                 ipivot, work, lw, icflag,
     &                 ifgfdj, f, bndxa, difbxa, bndxb, difbxb)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks required by
c       bacoli including:
c
c               determining yprime(t0).
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, November 8, 2001.
c Last modified by Jack Pew, July 11, 2013.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               t0 is the initial time.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bsplines
c                               coefficients (or DAEs).
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a,x_b].
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. fbasis(k,j,i) contains the
c                               values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        y(neq)
c                               y = y(t0) is the initial vector of
c                               bspline coefficients.
c
        integer                 ifgfdj
c                               Are finite differences being used to
c                               approximate boundary function Jacobi?
c                               Set to 0 or 1 if so,
c                               or to 2 or 3 to use user-provided
c                               analytic partial derivatives.
c
        external                f
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
c                               refer to the preamble of BACOLI.
c
c       Output:
        double precision        abdtop(npde*npde*nconti)
c                               The first block of the matrix A.
c
        double precision        abdbot(npde*npde*nconti)
c                               The last block of the matrix A.
c
        double precision        yprime(neq)
c                               yprime = yprime(t0) is the initial
c                               vector of bspline coefficients
c                               for the first temporal derivative.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        integer                 lw
c                               lw is the size of the work storage
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint+
c                                     npde*3+2*npde*npde+npde.
c
c       Work Storage:
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
        integer                 ipivot(neq)
c                               pivoting information from the
c                               factorization of the temporary matrix.
c
c-----------------------------------------------------------------------
c Local Variables:
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of abdtop
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of abdbot
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde), That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
        integer                 ifdwrk
c                               work(ifdwrk-1+1), i=1, 2*npde, is used
c                               as work storage for fdbndx.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               crdcmp
c                               crslve
c                               eval
c                               f
c                               difbxa
c                               difbxb
c                               fdbndx
c                               bndxa
c                               bndxb
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + npde*npde*kcol*(kcol+nconti)*nint
      iu     = iabdbt + npde*npde*nconti
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde
      ifdwrk = idbdt  + npde

c-----------------------------------------------------------------------
c     Initialize abdtop, abdbot and abdblk to zero.
      do 20 i = 1, npde * npde * nconti
         abdtop(i) = zero
         abdbot(i) = zero
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue

c-----------------------------------------------------------------------

c     Update the values at the left boundary.
      call eval(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y)
      if (ifgfdj .lt. 2) then
         call fdbndx(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxa, work(ifdwrk))
      else
         call difbxa(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Store -work(idbdt), which is the right side of the left boundary
c     conditions, into yprime.
      do 100 i = 1, npde
         yprime(i) = - work(idbdt-1+i)
  100 continue

c     Set up the top block and save in abdtop.
      do 120 j = 1, npde
         do 110 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            abdtop(ii) = work(idbdu-1+mm) - abdtop(jj)
  110    continue
  120 continue

c-----------------------------------------------------------------------
c     Generate the right side of ODEs at the collocation points
c     and save in yprime(i), i = npde + 1, neq - npde.
      do 140 i = 1, nint

c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol

         do 130 j = 1, kcol

c           jj is the index of the current collocation point.
            jj = (i - 1) * kcol + j + 1

c           mm is the pointer of yprime.
            mm = (jj - 1) * npde + 1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1,1,jj),y)

c           Evaluate the function f defining the PDE at the current
c           collocation point, storing the result in yprime.
            call f(t0,xcol(jj),work(iu),work(iux),work(iuxx),yprime(mm),
     &             npde)

  130    continue
  140 continue

c-----------------------------------------------------------------------
c     Update the values at the right boundary.
      call eval(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y)
      if (ifgfdj .lt. 2) then
         call fdbndx(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxb, work(ifdwrk))
      else
         call difbxb(t0, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Store -work(idbdt), which is the right side of the right boundary
c     conditions, into yprime.
      do 150 i = 1, npde
         ii = neq - npde + i
         yprime(ii) = - work(idbdt-1+i)
  150 continue

c     Set up the bottom block and save in abdbot.
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            abdbot(jj) = work(idbdu-1+mm) - abdbot(ii)
  160    continue
  170 continue

c-----------------------------------------------------------------------
c     Copy the collocation matrix into temporary storage.
      call dcopy(npde*npde*nconti, abdtop, 1, work(iabdtp), 1)
c
      call dcopy(npde*npde*kcol*(kcol+nconti)*nint, abdblk, 1,
     &           work(iabdbk), 1)
c
      call dcopy(npde*npde*nconti, abdbot, 1, work(iabdbt), 1)

c-----------------------------------------------------------------------
c     Generate the initial vector yp(t0).

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) go to 999
c
c     Solve the linear system. This gives yprime(t0)
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            yprime,0)

  999 continue
      return
      end

      subroutine reinit(npde, kcol, kold, nint, ninold, ncpts, neq,
     &                  neqold, icount, istep, nstep, x, xold, yold,
     &                  iflag, work, lw, ipivot, h, xbs, xcol,
     &                  fbasis, y, abdblk, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks after remeshing:
c
c               calculating the mesh step size sequence,
c               generating the piecewise polynomial space breakpoint
c               sequence,
c               calculating the collocation point sequence,
c               calculating the B-spline basis functions,
c               constructing abdblk of the collocation matrices and
c               calculating the bspline coefficients at the last nstep
c               steps which is needed for a warm start.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c Last modified by Paul Muir, Mar. 6, 2014: Commented out jj.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval after
c                               remeshing.
c
        integer                 kold
c                               kold is the number of collocation points
c                               to be used in each subinterval before
c                               remeshing.
c
        integer                 nint
c                               nint is the number of subintervals after
c                               remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bspline
c                               coefficients after remeshing.
c
        integer                 neqold
c                               neqold=npde*(kold*ninold+nconti) is
c                               the number of bspline
c                               coefficients before remeshing.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 istep
c                               istep is the number of time steps that
c                               DASSL has taken when kcol collocation
c                               points are used in each subinterval.
c
        integer                 nstep
c                               nstep is the number of time steps
c                               on which the remeshing is needed.
c
        double precision        x(nint+1)
c                               x is the spatial mesh after remeshing.
c
        double precision        xold(ninold+1)
c                               xold is the spatial mesh before
c                               remeshing.
c
        double precision        yold(6*neqold)
c                               yold is the vector of bspline
c                               coefficients at the last nstep time
c                               steps before remeshing.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        integer                 iflag
c                               iflag is a flag.
c                               iflag = 0, initialization is done for
c                                          dassl_kcol.
c                                       1, initialization is done for
c                                          dassl_kcol+1.
c
        integer                 lw
c                               lw is the size of the work storage
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +(kold+nconti)+kold*(ninold+1)
c                                     +2*nconti
c                               Since nint >= ninold/2 and kcol >=
c                               kold+1, it implies that lw >= 3*neqold.
c
c       Work Storage:
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
        integer                 ipivot(neq)
c                               pivoting information from the
c                               factorization of the temporary matrix.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points.
c
        double precision        y(nstep*neq)
c                               y is the vector of bspline coefficients
c                               at the last nstep time steps after
c                               remeshing.
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one
c                               collocation block of work.
c
        integer                 imod
c
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the
c                               bottom block which is required since
c                               crdcmp overwrites the input collocation
c                               matrix.
c
        integer                 ivwork
c                               work(ivwork) is the work storage
c                               required by revalu.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
c       integer                 jj  ! No longer used.
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               colpnt
c                               crdcmp
c                               crslve
c                               revalu
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------
c Fortran Functions Used:
c                               mod
c
c-----------------------------------------------------------------------

c     Generate the piecewise polynomial space breakpoint sequence,
c     and calculates the collocation point sequences.
      call colpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c     Update yold.
c     The modulus operations and dcopys in this section are to order the
c     steps in order. They undo the shuffling that was done by the
c     subroutine sucstp in the name of efficiency.
      if ((iflag .ne. 0) .or. (icount .ne. 1)) goto 5

      imod = mod(istep, 6)
      if (imod .eq. 0) then
         imod = 6
      endif
      imod = 7 - imod
      if (imod .eq. 1) goto 5

      if ((imod+nstep-1) .gt. 6) then
         if (imod .le. 3) then
            call dcopy((imod+nstep-1-6)*neqold, yold, 1, work, 1)
            call dcopy((6-imod+1)*neqold, yold((imod-1)*neqold+1), 1,
     &                 yold, 1)
            call dcopy((imod+nstep-1-6)*neqold, work, 1,
     &                 yold((6-imod+1)*neqold+1), 1)
         else
            call dcopy((6-imod+1)*neqold, yold((imod-1)*neqold+1), 1,
     &                 work, 1)
            call dcopy((imod+nstep-1-6)*neqold, yold, -1,
     &                 yold((6-imod+1)*neqold+1), -1)
            call dcopy((6-imod+1)*neqold, work, 1, yold, 1)
         endif
      else
         call dcopy(nstep*neqold, yold((imod-1)*neqold+1), 1, yold, 1)
      endif

    5 continue

c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      ivwork = iabdbt + npde*npde*nconti

c-----------------------------------------------------------------------
c     Call revalu to calculate the values at xcol and at the last nstep
c     time step. Then save in y.
      call revalu(kold, xcol, ninold, xold, npde, ncpts, nstep,
     &            y, yold, work(ivwork))

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 10 i = 1, npde * npde * nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   10 continue
      do 20 i = 1, nint*nels
         abdblk(i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j)
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Making use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 30 i = 1, npde
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   30 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 70 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is no longer used so we comment out this assignment:
c     jj is the position in the y vector where the values for the
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
c           jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bsplvd(xbs,kcol+nconti,xcol(ii),ileft,fbasis(1,1,ii),3)

            do 50 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 40 m = 1, npde
                  mm = ll + (m-1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      do 80 i = 1, npde
         ii = ((i-1)+npde)*npde + i
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   80 continue

c-----------------------------------------------------------------------
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nels*nint,abdblk,1,work(iabdbk),1)

c-----------------------------------------------------------------------
c     Generate the vector y.

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) go to 999

c     Solve the linear system. This gives the basis function
c     coefficients for the initial conditions, i.e. y(t0).
      do 90 i = 1, nstep
         ii = (i - 1) * neq + 1
         call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &               (kcol+nconti)*npde,nint,work(iabdbt),npde,
     &               ipivot,y(ii),0)
         if (icflag .ne. 0) go to 999
   90 continue

  999 return
      end

      subroutine revalu(kcol, xsol, nint, x, npde, npts, nstep, usol,
     &                  y, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the solution u at the npts points xsol
c       and at the current and previous nstep-1 time step. Then return
c       them in the array usol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Mar 31, 2006.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector.
c
        double precision        xsol(npts)
c                               xsol is an arbitrary set of spatial
c                               points at which the solution and the
c                               first nderv derivative values are
c                               to be calculated.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a,x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 nstep
c                               nstep-1 is the number of previous steps.
c                               When user wants to calculate solution
c                               at tout, let nstep = 1.
c
        double precision        y(npde*(nint*kcol+nconti), nstep)
c                               y is the vector of bspline
c                               coefficients at the current time step
c                               and previous nstep-1 steps.
c
c       output:
        double precision        usol(npde, npts, nstep)
c                               usol is the solution at the given
c                               points and at the current time step
c                               and previous nstep-1 steps.
c
c       Work Storage:
        double precision        work((kcol+nconti)+kcol*(nint+1)
     *                               +2*nconti)
c                               work is a floating point work storage
c                               array.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 mflag
c                               mflag is required by subroutine
c                               interv.
c
        integer                 ilo
c                               ilo is required by subroutine
c                               interv.
c       Pointers into the floating point work array:
        integer                 ixbs
c                               work(ixbs) contains the breakpoint
c                               sequence.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 k
        integer                 m
        integer                 n
        integer                 ii
        integer                 mm
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               interv
c
c-----------------------------------------------------------------------

c     set up the value for ileft and ncpts.
      ileft = 0
      ncpts = nint * kcol + nconti

c     set the pointer into the floating point work array
      ixbs  = (kcol+nconti) + 1

c     Store the piecewise polynomial space breakpoint sequence in
c     work(ixbs).
c
      do 10 i = 1, kcol + nconti
         work(ixbs-1+i) = x(1)
         work(ixbs-1+i+ncpts) = x(nint+1)
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
   20 continue
   30 continue

      do 70 n = 1, nstep
c
c     set up the value for ileft and ncpts.
         mflag = -2
         do 60 i = 1, npts
c
c     interv is called to compute ileft. bsplvd is called to compute
c     the values of the basis function at the required point.
            call interv(work(ixbs), ncpts, xsol(i), ileft, mflag, ilo)
            call bsplvd(work(ixbs),kcol+nconti,xsol(i),ileft,work,1)
            ii = ileft - kcol - nconti
            do 50 k = 1, npde
               usol(k,i,n) = zero
               do 40 m = 1, kcol + nconti
                  mm = (m + ii - 1) * npde
                  usol(k,i,n) = usol(k,i,n) + y(mm+k,n) * work(m)
   40          continue
   50       continue
   60    continue
   70 continue
      return
      end

      subroutine eval(npde,kcol,ileft,icpt,ncpts,uval,uxval,uxxval,
     &                fbasis,y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine evaluates u(k), ux(k), and uxx(k), k=1 to npde,
c       at the icpt-th collocation point.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Feb. 11, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 icpt
c                               the index of the collocation point.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        fbasis((kcol+nconti)*3)
c                               Basis function values at the icpt-th
c                               collocation point.
c                               fbasis(k+(j-1)*(kcol+nconti)) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti).
c
        double precision        y(ncpts*npde)
c                               y is the vector of bspline coefficients.
c
c       Output:
        double precision        uval(npde)
c                               uval gives the approximation to
c                               u(t,x).
c
        double precision        uxval(npde)
c                               uxval gives the approximation to
c                               the first spatial derivative of u(t,x).
c
        double precision        uxxval(npde)
c                               uxxval gives the approximation to
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 j
        integer                 m
        integer                 mj
        integer                 mj0
c-----------------------------------------------------------------------
      do 10 j = 1, npde
         uval(j)   = zero
         uxval(j)  = zero
         uxxval(j) = zero
   10 continue
      if (icpt .ne. 1 .and. icpt .ne. ncpts) then
         mj0 = (ileft-kcol-3) * npde
         do 30 j = 1, npde
            mj0 = mj0 + 1
            mj = mj0
            do 20 m = 1, kcol + nconti
               mj = mj + npde
               uval(j)   = uval(j)   + fbasis(m)                 * y(mj)
               uxval(j)  = uxval(j)  + fbasis(m+kcol+nconti)     * y(mj)
               uxxval(j) = uxxval(j) + fbasis(m+2*(kcol+nconti)) * y(mj)
   20          continue
   30       continue
      else
         if (icpt .eq. 1) then
            do 40 j = 1, npde
               uval(j)   = uval(j)   + fbasis(1) * y(j)
               uxval(j)  = uxval(j)  + fbasis(1+kcol+nconti) * y(j)
     &                     + fbasis(2+kcol+nconti) * y(npde + j)
               uxxval(j) = uxxval(j) + fbasis(1+2*(kcol+nconti)) * y(j)
     &                     + fbasis(2+2*(kcol+nconti)) * y(npde + j)
     &                     + fbasis(3+2*(kcol+nconti)) * y(2*npde + j)
   40       continue
         else
            do 50 j = 1, npde
               uval(j)   = uval(j)   + fbasis(kcol+nconti)
     &                     * y((ncpts - 1) * npde + j)
               uxval(j)  = uxval(j)  + fbasis((kcol+nconti)*2)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*2-1)
     &                     * y((ncpts - 2) * npde + j)
               uxxval(j) = uxxval(j) + fbasis((kcol+nconti)*3)
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*3-1)
     &                     * y((ncpts - 2) * npde + j)
     &                     + fbasis((kcol+nconti)*3-2)
     &                     * y((ncpts - 3) * npde + j)
   50       continue
         endif
      endif
      return
      end

      subroutine divdif(neq, nstep, psi, work, y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the divided difference, which is required
c       by DASSL for a hot start, after calculating the bspline
c       coefficients at the last nstep steps.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 4, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of bspline
c                               coefficients after remeshing.
c
        integer                 nstep
c                               nstep is the number of time steps
c                               on which the remeshing is needed.
c
        double precision        psi(6)
c                               psi is the stepsize vector of the
c                               previous 6 time steps.
c
c       Work Storage:
        double precision        work(6)
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        y(nstep*neq)
c                               y is the vector of bspline coefficients
c                               at the last nstep time steps after
c                               remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 m
c
c-----------------------------------------------------------------------

c     Update y to be the divide difference.

      do 40 i = 1, nstep - 1
         work(1) = psi(i)
         if (nstep .gt. 2) then
            do 10 j = 2, nstep-i
               work(j) = psi(j+i-1) - psi(j-1)
   10       continue
         endif
         do 30 j = nstep, i+1, -1
            do 20 m = 1, neq
               y((j-1)*neq+m) = (y((j-2)*neq+m) - y((j-1)*neq+m))
     *                          / work(j-i)
   20       continue
   30    continue
   40 continue

      work(1) = psi(1)
      do 50 i = 2, nstep - 1
         work(i) = work(i-1) * psi(i)
   50 continue
      do 70 i = 2, nstep
         do 60 m = 1, neq
            y((i-1)*neq+m) = y((i-1)*neq+m) * work(i-1)
   60    continue
   70 continue

c-----------------------------------------------------------------------
      return
      end

      subroutine jac(t, y, yprime, pd, cj, rpar, ipar,
     &               derivf, bndxa, difbxa, bndxb, difbxb)
c-----------------------------------------------------------------------
c Purpose:
c       This is the subroutine which defines the Jacobian of the
c       differential/algebraic system to be solved by DASSL. It returns:
c                       PD := dG/dY + cj * dG/dY'
c       To be precise, the (i,j)th element of PD involves the partial
c       derivative of equation i with respect to the variable j (or its
c       time derivative). Or in pseudo-code we have:
c
c               PD(i,j) = dG(i)/dY(j) + cj * dG(i)/dYprime(j).
c
c       The DAE G(t, Y, Y') = 0 arises from applying the method-of-lines
c       and bspline collocation to the system of NPDE PDES of
c       the form:
c                       u_t = f(t, x, u, u_x, u_xx)
c       In the discretized form this yields:
c                      G(t, Y, Y') = A*Y' - F~
c       The abd matrix A contains the collocation equations and some
c       boundary condition information, the vector F~ contains the rhs
c       of the collocation equations and the corresponding boundary
c       conditions.
c
c       In view of this, we have:
c                      PD = cj * A - dF~/dY.
c       Now by the product rule we can express df/dY as:
c                   df/dY = df/du * du/dY +
c                           df/du_x * du_x/dY +
c                           df/du_xx * du_xx/dY.
c       So, in this fashion the elements of dF~/dY can be calculated.
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, August 13, 2001.
c Last modified by Jack Pew, May 17, 2013.
c Last modified by Paul Muir, Mar. 6, 2014: Commented out nconti.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
c       integer                 nconti
c       parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t
c                               t is the current time.
c
        double precision        y(*)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(*)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        double precision        cj
c                               cj is a scalar chosen by DASSL to
c                               accelerate convergence of the modified
c                               Newton iteration used to solve the
c                               implicit equations resulting from the
c                               BDF methods.
c
        double precision        rpar(*)
c                               rpar is the BACOLI floating point work
c                               array.
c
        integer                 ipar(*)
c                               rpar is the BACOLI integer work array.
c
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
c                               refer to the preamble of BACOLI.
c
c       Output:
        double precision        pd(*)
c                               pd is the ABD Jacobian (iteration)
c                               matrix of the residual of the DAE
c                               system defined by RES.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
c                               ipar(inpde) = npde
c
        integer                 ikcol
c                               ipar(ikcol) = kcol.
c
        integer                 inint
c                               ipar(inint) = nint.
c
        integer                 incpts
c                               ipar(incpts) = ncpts.
c
        integer                 ineq
c                               ipar(ineq) = neq.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ixcol
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               dassl_kcol.
c
        integer                 iabtop
c                               rpar(ipar(iabtop)) stores the top block
c                               of the ABD collocation matrices when
c                               using dassl_kcol.
c
        integer                 iabblk
c                               rpar(ipar(iabblk)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_kcol.
c
        integer                 iabbot
c                               rpar(ipar(iabbot)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using dassl_kcol.
c
        integer                 iwkrj
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using dassl_kcol.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 imflg9
c                               ipar(imflg9) stores the value of
c                               mflag(9) in order for it to be passed
c                               down to caljac. If this is 0 or 1, then
c                               the user has not implemented difbxa and
c                               difbxb and caljac should use fdbndx to
c                               replace them with finite differences.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts
        integer                 neq
c
c-----------------------------------------------------------------------
c       Direct IPAR indices:
        parameter              (inpde  =  1)
        parameter              (ikcol  =  2)
        parameter              (inint  =  3)
        parameter              (incpts =  4)
        parameter              (ineq   =  5)
c
c-----------------------------------------------------------------------
c       IPAR indices which serve as an indirect pointer into RPAR:
        parameter              (imflg9 = 16)
        parameter              (ixcol  = 32)
        parameter              (iabtop = 36)
        parameter              (iabblk = 37)
        parameter              (iabbot = 38)
        parameter              (iwkrj  = 40)
        parameter              (ibasi  = 41)
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                              caljac
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts  = ipar(incpts)
      neq    = ipar(ineq)

c     Calculate jacobian for dassl_kcol.
      call caljac(npde, kcol, nint, ncpts, neq, rpar(ipar(ixcol)),
     &            rpar(ipar(iabtop)), rpar(ipar(iabblk)),
     &            rpar(ipar(iabbot)), rpar(ipar(ibasi)), t, y,
     &            yprime, cj, rpar(ipar(iwkrj)), pd,
     &            ipar(imflg9), derivf, bndxa, difbxa, bndxb, difbxb)

c     Second call to caljac removed.

      return
      end

      subroutine caljac(npde, kcol, nint, ncpts, neq, xcol, abdtop,
     &                  abdblk, abdbot, fbasis, t, y, yprime, cj,
     &                  work, pd,
     &                  ifgfdj, derivf, bndxa, difbxa, bndxb, difbxb)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by jac. It provides a lower-level
c       interface to generate the iteration matrix.
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, August 13, 2001.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bspline
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdtop(npde*npde*nconti)
c                               abdtop stores the top block of the ABD
c                               matrices.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
        double precision        abdbot(npde*npde*nconti)
c                               abdbot stores the bottom block of the
c                               ABD matrices.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(neq)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        double precision        cj
c                               cj is a scalar chosen by DASSL to
c                               accelerate convergence of the modified
c                               Newton iteration used to solve the
c                               implicit equations resulting from the
c                               BDF methods.
c
        integer                 ifgfdj
c                               Are finite differences being used to
c                               approximate boundary function Jacobi?
c                               Set to 0 or 1 if so,
c                               or to 2 or 3 to use user-provided
c                               analytic partial derivatives.
c
        external                derivf
        external                bndxa
        external                difbxa
        external                bndxb
        external                difbxb
c                               refer to the preamble of BACOLI.
c
c       Work storage:
        double precision        work(4*npde+5*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+5*npde*npde.
c
c       Output:
        double precision        pd(npde*npde*(2*nconti
     *                             +nint*kcol*(kcol+nconti)))
c                               pd is the ABD Jacobian (iteration)
c                               matrix of the residual of the DAE
c                               system defined by RES.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ipdtop
c                               ipdtop is the pointer into pd where the
c                               top block of the ABD Jacobian is stored.
c
        integer                 ipdblk
c                               ipdblk is the pointer into pd where the
c                               nint blocks in the middle of the ABD
c                               Jacobian are stored.
c
        integer                 ipdbot
c                               ipdbot is the pointer into pd where the
c                               bottom block of the ABD Jacobian is
c                               stored.
c
        integer                 nsiztb
c                               nsiztb is the size of the top block
c                               as same as the bottom block of the ABD
c                               Jacobian.
c
        integer                 nsizbk
c                               nsizbk is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 n
c
        integer                 ii
        integer                 ij
        integer                 jj
        integer                 kk
        integer                 nn
        integer                 mm
        integer                 jk
        integer                 jk2
        integer                 jk3
        integer                 mn
        integer                 mn2
        integer                 mn3
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idfdu
c                               work(idfdu) stores the Jacobian of f
c                               with respect to u.
c
        integer                 idfdux
c                               work(idfdux) stores the Jacobian of f
c                               with respect to u_x.
c
        integer                 idfuxx
c                               work(idfuxx) stores the Jacobian of f
c                               with respect to u_xx.
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde), That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
        integer                 ifdwrk
c                               work(ifdwrk-1+1), i=1, 2*npde, is used
c                               as work storage for fdbndx.
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               eval
c                               derivf
c                               difbxa
c                               difbxb
c                               fdbndx
c                               bndxa
c                               bndxb
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               daxpy
c
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde
      idfdu  = iuxx   + npde
      idfdux = idfdu  + npde * npde
      idfuxx = idfdux + npde * npde
      idbdu  = idfuxx + npde * npde
      idbdux = idbdu  + npde * npde
      idbdt  = idbdux + npde * npde
      ifdwrk = idbdt  + npde

c     Set the indices into pd which define the ABD Jacobian.
      ipdtop = 1
      ipdblk = ipdtop + nconti * npde * npde
      ipdbot = ipdblk + nint * npde * npde * kcol * (kcol + nconti)

c-----------------------------------------------------------------------
c     Calculate the size of top (or bottom) block and the size of a
c     subblock in the middle.
      nsiztb = npde * npde * nconti
      nsizbk = npde * npde * kcol * (kcol + nconti)

c     Initialize pdtop, pdblk and pdbot to zero.
      do 10 i = 1, nsiztb
         pd(ipdtop-1+i) = zero
         pd(ipdbot-1+i) = zero
   10 continue
      do 20 i = 1, nint * nsizbk
         pd(ipdblk-1+i) = zero
   20 continue
c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations and
c     calculate the portion of dG/dY which depends on them.

      do 70 i = 1, nint

c        ii+1 is the pointer to the first element at the i-th subblock
c        of the jacobian matrix, i = 1, nint.
         ii = ipdblk - 1 + (i - 1) * nsizbk

c        ij is the value of ileft for the current collocation point.
         ij = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c           jj+1 is the pointer to the first element corresponding to
c           the j-th collocation point in the i-th interval.
            jj = ii + (j - 1) * npde

c           mm is the index of the current collocation point.
            mm = (i - 1) * kcol + j + 1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ij,mm,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1+(mm-1)*(kcol+nconti)*3),y)

c           Generate dfdu, dfdux, and dfdux at the current
c           collocation point (the j-th point of the i-th
c           subinterval).
            call derivf(t, xcol(1+(i-1)*kcol+j), work(iu),
     &                  work(iux), work(iuxx), work(idfdu),
     &                  work(idfdux), work(idfuxx), npde)

            do 50 k = 1, kcol + nconti

c              kk+1 is the pointer to the first element of a npde by
c              npde submatrix, which is corresponding to the j-th
c              collocation point in the i-th interval, and the k-th
c              nonzero basis function.
               kk = jj + (k-1) * npde * npde * kcol

c              jk is the pointer to the k-th nonzero function at the
c              mm-th collocation point in the basis function,
c              fbasis(1).
               jk = (mm - 1) * (kcol + nconti) * 3 + k

c              jk2 is the pointer to the first derivative for the
c              above basis function.
               jk2 = jk + kcol + nconti

c              jk3 is the pointer to the second derivative for the
c              above basis function.
               jk3 = jk2 + kcol + nconti

               do 40 m = 1, npde
                  do 30 n = 1, npde

c                    nn is the pointer to the (n, m) element of the
c                    npde by npde submatrix.
                     nn = kk + (m-1)*npde*kcol + n

c                    mn is the pointer to the (n, m) element of dfdu.
                     mn = idfdu - 1 + (m - 1) * npde + n

c                    mn2 is the pointer to the (n, m) element of dfdux.
                     mn2 = mn + npde * npde

c                    mn3 is the pointer to the (n, m) element of dfduxx.
                     mn3 = mn2 + npde * npde

c                    now set up the value in pd at the place nn.
                     pd(nn) = - work(mn) * fbasis(jk)
     &                        - work(mn2) * fbasis(jk2)
     &                        - work(mn3) * fbasis(jk3)

   30             continue
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Update the values at the left boundary.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      if (ifgfdj .lt. 2) then
         call fdbndx(t, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxa, work(ifdwrk))
      else
         call difbxa(t, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Update the top block of the collocation matrix dG/dY'.
      do 90 j = 1, npde
         do 80 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) =
     &            fbasis(2+kcol+nconti) * work(idbdux-1+mm)
            abdtop(ii) =
     &            work(idbdu-1+mm) - abdtop(jj)
   80    continue
   90 continue

c-----------------------------------------------------------------------
c     Update the values at the right boundary.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      if (ifgfdj .lt. 2) then
         call fdbndx(t, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde,
     &               bndxb, work(ifdwrk))
      else
         call difbxb(t, work(iu), work(iux), work(idbdu),
     &               work(idbdux), work(idbdt), npde)
      end if

c     Update the bottom block of the collocation matrix.
      do 110 j = 1, npde
         do 100 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) =
     &            fbasis(1+kcol+kcol+nconti+(ncpts-1)*(kcol+nconti)*3)
     &            * work(idbdux-1+mm)
            abdbot(jj) =
     &            work(idbdu-1+mm) - abdbot(ii)
  100    continue
  110 continue

c-----------------------------------------------------------------------
c     Add cj * A to pdmat. This is cj * dG/dY'.

      call daxpy(nsiztb, cj, abdtop, 1, pd(ipdtop), 1)

      call daxpy(nint*nsizbk, cj, abdblk, 1, pd(ipdblk), 1)

      call daxpy(nsiztb, cj, abdbot, 1, pd(ipdbot), 1)

c-----------------------------------------------------------------------
      return
      end

      subroutine res(t, y, yprime, delta, ires, rpar, ipar,
     &               f, bndxa, bndxb)

c-----------------------------------------------------------------------
c Purpose:
c       This is the subroutine which defines the differential/algebraic
c       system to be solved by DASSL. It returns the residual (delta) of
c       the DAE system:
c                          delta := G(t, Y, Y')
c       The DAE G(t, Y, Y') = 0 arises from applying the method-of-lines
c       and bspline collocation to the system of NPDE PDES of
c       the form:
c                       u_t = f(t, x, u, u_x, u_xx)
c       In the discretized form this yields:
c                      G(t, Y, Y') = A*Y' - F~
c       The abd matrix A contains the collocation equations and some
c       boundary condition information, the vector F~ contains the rhs
c       of the collocation equations and the corresponding boundary
c       conditions.
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, November 8, 2001.
c Last modified by Jack Pew, May 17, 2013.
c
c-----------------------------------------------------------------------
        implicit none
c Subroutine Parameters:
c       Input:
        double precision        t
c                               T is the current time.
c
        double precision        y(*)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(*)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        double precision        rpar(*)
c                               rpar is the BACOLI floating point work
c                               array.
c
        integer                 ipar(*)
c                               rpar is the BACOLI integer work array.
c
        external                f
        external                bndxa
        external                bndxb
c                               refer to the preamble of BACOLI.
c
c       Output:
        integer                 ires
c                               ires is a user flag set to alert DASSL
c                               of an illegal y vector. ires is not
c                               used in this subroutine since there is
c                               no restriction on the vector of
c                               bspline coefficients.
c
        double precision        delta(*)
c                               delta is the residual of the DAE system.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
c                               ipar(inpde) = npde
c
        integer                 ikcol
c                               ipar(ikcol) = kcol.
c
        integer                 inint
c                               ipar(inint) = nint.
c
        integer                 incpts
c                               ipar(incpts) = ncpts.
c
        integer                 ineq
c                               ipar(ineq) = neq.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ixcol
c                               rpar(ipar(ixcol)) stores the
c                               collocation points when using
c                               dassl_kcol.
c
        integer                 iabblk
c                               rpar(ipar(iabblk)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               dassl_kcol.
c
        integer                 iwkrj
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi
c                               rpar(ipar(ibasi)) stores the basis
c                               function values at the collocation
c                               points when using dassl_kcol.
c                               rpar(ipar(ibasi)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts
        integer                 neq
c
c-----------------------------------------------------------------------
c       Direct IPAR indices:
        parameter              (inpde  =  1)
        parameter              (ikcol  =  2)
        parameter              (inint  =  3)
        parameter              (incpts =  4)
        parameter              (ineq   =  5)
c
c-----------------------------------------------------------------------
c       IPAR indices which serve as an indirect pointer into RPAR:
        parameter              (ixcol  = 32)
        parameter              (iabblk = 37)
        parameter              (iwkrj  = 40)
        parameter              (ibasi  = 41)
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                              calres
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts  = ipar(incpts)
      neq    = ipar(ineq)

c     Calculate residual for dassl.
      call calres(npde, kcol, nint, ncpts, neq, rpar(ipar(ixcol)),
     &            rpar(ipar(iabblk)), rpar(ipar(ibasi)), t, y,
     &            yprime, rpar(ipar(iwkrj)), delta,
     &            f, bndxa, bndxb)

c     Second call to calres removed.

      return
      end

      subroutine calres(npde, kcol, nint, ncpts, neq, xcol, abdblk,
     &                  fbasis, t, y, yprime, work, delta,
     &                  f, bndxa, bndxb)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by res. It provides a lower-level
c       interface to generate the residue at the current time t.
c
c-----------------------------------------------------------------------
c
c Modified by Rong Wang, November 8, 2001.
c Last modified by Jack Pew, May 17, 2013.
c
c-----------------------------------------------------------------------
        implicit none
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        negone
        parameter              (negone = -1.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bsplines
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               T is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        double precision        yprime(neq)
c                               yprime is the derivative of y with
c                               respect to time at the current time.
c
        external                f
        external                bndxa
        external                bndxb
c                               refer to the preamble of BACOLI.
c
c       Work storage:
        double precision        work(4*npde+2*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+2*npde*npde.
c
c       Output:
        double precision        delta(neq)
c                               delta is the residual of the DAE system.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 m
        integer                 k
c
c       Indices:
        integer                 ii
        integer                 jj
        integer                 mm
        integer                 kk
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               f
c                               eval
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c BLAS Subroutines:
c       double precision:
c                               dscal
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde

c-----------------------------------------------------------------------
c     Initialize the residual to the zero vector.
      do 10 i = 1, neq
         delta(i) = zero
   10 continue

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations.

      do 30 i = 1, nint

c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol

         do 20 j = 1, kcol

c           jj is the pointer of collocation point.
            jj = (i - 1) * kcol + j + 1

c           mm is the pointer of delta.
            mm = (jj - 1) * npde + 1

c           kk is the pointer of the basis function values at
c           the current collocation point.
            kk =(jj-1)*(kcol+nconti)*3+1

c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(kk),y)

c           Evaluate the function f defining the PDE at the current
c           collocation point, storing the result in delta.
            call f(t, xcol(jj), work(iu), work(iux),
     &              work(iuxx), delta(mm), npde)

   20    continue
   30 continue

c     Scale (delta(i), i=npde+1,npde*(ncpts-1)) with negative one.
      call dscal(npde*kcol*nint, negone, delta(npde+1), 1)

c-----------------------------------------------------------------------
c     Calculate the portion of the residual vector which depends on the
c     collocation matrix. delta := delta + A*Yprime.

c     Calculate (delta(i), i=1, npde), which depend on the left
c     boundary point.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call bndxa(t, work(iu), work(iux), delta(1), npde)

c-----------------------------------------------------------------------
c     Calculate (delta(i), i=1, npde), which depend on the right
c     boundary point.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call bndxb(t, work(iu), work(iux), delta(neq-npde+1), npde)

c-----------------------------------------------------------------------
c     Calculate (delta(i), i = npde+1, (ncpts-1)*npde), which depend
c     on the nint blocks in the middle of the collocation matrix A.
      do 70 i = 1, nint
         do 60 j = 1, kcol + nconti
            do 50 k = 1, kcol
               kk = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(k-1)*npde
               do 40 m = 1, npde
                  ii = npde+(i-1)*npde*kcol+(k-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  delta(ii) = delta(ii) + abdblk(kk) * yprime(mm)
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
      return
      end

c--- Interpolation based spatial error estimation code.

      subroutine lowint(kcol, nint, npde, istart, icount, xbs,
     &                  x, y, h, hd, g, bassc, basm, u, um, usol)
c-----------------------------------------------------------------------
c     This subroutine evaluates a Hermite-Birkhoff (H-B) interpolant of
c     order kcol (degree kcol-1) at kcol+2 Gaussian quadrature points on
c     the interval (0,1). (These values will be used elsewhere to obtain
c     an L2-norm error estimate.)
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
      integer                 nconti
      parameter              (nconti = 2)
c     nconti continuity conditions are imposed at the internal mesh
c     points of the collocation solution.
c
      double precision        zero
      parameter              (zero = 0.0d0)
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     Number of collocation points per subinterval for the collocation
c     solution. The collocation solution is a piecewise polynomial of
c     order kcol+1 (degree kcol.) The values obtained from this
c     interpolant are of order kcol. The interpolant is a piecewise
c     polynomial of degree kcol-1.
c
      integer                 nint
c     Number of subintervals on the current mesh.
c
      integer                 npde
c     Number of PDEs in the system being solved.
c
      integer                 istart
c     istart indicates whether or not the integration has just begun.
c         istart = 0, it is the initial step;
c                = 1, it is not the initial step.
c
      integer                 icount
c     icount is the number of remeshings that have been performed at
c     the current time step.
c     istart and icount are used to determine whether certain
c     coefficients need to be calculated, depending upon whether the
c     mesh has changed.
c
      double precision        xbs((nint+1)*kcol+nconti*2)
c     The breakpoint sequence. This is required by BSPLVD to compute the
c     values of the B-spline basis functions.
c
      double precision        x(nint+1)
c     The current mesh point sequence.
c
      double precision        y(npde*(nint*kcol+nconti))
c     B-spline coefficients computed by DASSL.
c
c-----------------------------------------------------------------------
c Inoutput
      double precision        h(2, kcol+2)
      double precision        hd(2, kcol+2)
      double precision        g(max(1,kcol-3), kcol+2)
c     Coefficients of the Hermite-Birkhoff interpolant.
c     h and hd coefficients multiply the values and derivatives of the
c     collocation solution at mesh points, g coefficients multiply the
c     values of the collocation solution at the internal points.
c     For the LOI scheme, no points external to the subinterval are
c     used, and these coefficient do not depend upon the mesh.
c     This means they really only need to be computed once,
c     but then the values would need to be copied from place to place
c     whenever nint changed and rpar grew or shrank. This copying would
c     potentially by of some confusion, as it would have to be done in
c     the main BACOLI loop. So while there is unnecessary work done here
c     (to recompute these coefficients,) we do it anyway since it takes
c     an insignificant fraction of the total CPU time.
c
      double precision        bassc(kcol+nconti, max(1,(kcol-3)*nint))
c     Values of all nonzero B-spline basis functions at the
c     interpolated points inside the subintervals.
c     These are stored between remeshings.
c
      double precision        basm(kcol+nconti, 2, nint+1)
c     Values of all nonzero B-spline basis functions and their first
c     derivatives at the mesh points. kcol of every kcol+nconti of these
c     are actually zero, but the trouble to reduce the storage used by
c     this vector seemed too great.
c     These are stored between remeshings.
c
c-----------------------------------------------------------------------
c Output
      double precision        usol(npde, nint*(kcol+2))
c     Array of evaluations of the Hermite-Birkhoff interpolant at the
c     quadrature points used for the computation of an L2-norm error
c     estimate elsewhere.
c
c-----------------------------------------------------------------------
c Work storage
      double precision        u(npde, max(1,(kcol-3)*nint))
c     Work array for storage of values of the collocation solution at
c     the internal points used in the Hermite-Birkhoff interpolant.
c
      double precision        um(2, npde, nint+1)
c     Work array for storage of values of the collocation solution, and
c     its first derivative, at the mesh points used in the Hermite-
c     Birkhoff interpolant.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        scpt(mxkcol-2)
c     The kcol-2 internal interpolation points, over (0,1), used in the
c     construction of the Hermite-Birkhoff interpolant.
c     These points are not actually superconvergent for the
c     computed solution, as they would be in the SCI case. They are what
c     would be superconvergent points, if kcol were one lower.
c     This choice of points is important in order to have the leading
c     terms of the interpolation and collocation errors match.
c
      double precision        theta(mxkcol+2)
c     Quadrature points for a subinterval, which are the Gauss points
c     on (0,1).
c
      double precision        gwts(mxkcol+2)
c     Gaussian weights. Unused here.
c
      double precision        gwork((mxkcol+2)*(mxkcol+2))
c     Work storage for gauleg.
c
      integer                 ileft
c     Breakpoint information. ileft is used in conjunction with the
c     breakpoint sequence for the evaluation of B-spline basis functions
c     by BSPLVD.
c
      double precision        hx
c     Width of the ith subinterval, x(i+1) - x(i).
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j, k, m, ii, jj, mm
c
c-----------------------------------------------------------------------
c Subroutines Called
c                             bsplvd
c                             gauleg
c                             lowcoeff
c                             setpts
c
c-----------------------------------------------------------------------

c     Check whether the H-B coefficients and B-spline basis function
c     evaluations need to be computed. Previous values may be reused if
c     the mesh has not changed since the last step.
      if ((icount .ne. 0) .or. (istart .eq. 0)) then

c        Set thetas (weights will not be computed since last parameter
c        in the call to gauleg is set to 2.)
         call gauleg(kcol+2, (kcol+2)*(kcol+2), theta, gwts, gwork, 2)

c        H-B coefficients are computed for the Gaussian quadrature
c        points. They are the same on every subinterval, so only one set
c        is needed.
c        Note that LOWCOEFF is called with kcol-1 instead of kcol. (This
c        is since for a given kcol and a collocation solution
c        corresponding to this given kcol, we actually want to construct
c        a LOI scheme that is one order lower.)
         do j = 1, kcol+2
            call lowcoeff(theta(j), kcol-1, h(1,j), hd(1,j), g(1,j))
         end do

c        Call BSPLVD to obtain values of B-spline basis functions at the
c        internal interpolated points.
c        The subroutine INTERV does not need to be called to find ileft,
c        since we know which subinterval we are in already.
c
c        Note that setpts is called with kcol-1, so the points returned
c        are not in fact superconvergent for the computed collocation
c        solution. It is important that these points be used however so
c        that the leading terms of the interpolation and collocation
c        errors match.
c        There are no such points for the LOI scheme when kcol = 3.
         if (kcol .ne. 3) then
            ileft = nconti
            do i = 1, nint
               call setpts(scpt, kcol-1, i, x, nint)
               ileft = ileft + kcol
               do j = 1, kcol-3
                  jj = (i-1) * (kcol-3) + j
                  call bsplvd(xbs, kcol+nconti, scpt(j), ileft,
     &                        bassc(1, jj), 1)
               end do
            end do
         endif

c        Call BSPLVD for the values and first derivatives of the
c        B-spline basis functions at the nint+1 mesh points.
c       (ileft does not follow the pattern at the rightmost mesh point,
c        so we need a separate call for the last mesh point.)
         ileft = nconti
         do i = 1, nint
            ileft = ileft + kcol
            call bsplvd(xbs, kcol+nconti, x(i), ileft,
     &                  basm(1, 1, i), 2)
         end do
c        ileft = kcol+nconti+kcol*(nint-1)
         call bsplvd(xbs, kcol+nconti, x(nint+1), ileft,
     &               basm(1, 1, nint+1), 2)

      endif

c-----------------------------------------------------------------------
c     Compute the values of the collocation solution at the internal
c     interpolated points and mesh points by multiplying the values of
c     the basis functions at those points by the B-spline coefficients
c     computed by DASSL.
c
c     The end result of these two steps (evaluating the basis functions
c     and later multiplying by the coefficients) is the same as calling
c     the subroutine VALUES, but this approach avoids a large number
c     of unnecessary calls to BSPLVD by saving some of the required
c     values.
      ii = 0
      do i = 1, nint
         do j = 1, kcol-3
            jj = (i-1) * (kcol-3) + j
            do k = 1, npde
               u(k, jj) = zero
               mm = ii
               do m = 1, kcol+nconti
c                 mm = (m + (i-1) * kcol - 1) * npde
                  u(k, jj) = u(k, jj) + y(mm+k) * bassc(m, jj)
                  mm = mm + npde
               end do
            end do
         end do
         ii = ii + kcol * npde
      end do
c     All but two of the values of the basis functions at the mesh
c     points are zero, so um gets special treatment.
      mm = 0
      do i = 1, nint
         do k = 1, npde
            um(1, k, i) = y(mm+k) * basm(1, 1, i)
     &             + y(npde+mm+k) * basm(2, 1, i)
            um(2, k, i) = y(mm+k) * basm(1, 2, i)
     &             + y(npde+mm+k) * basm(2, 2, i)
         end do
         mm = mm + kcol * npde
      end do
c     It is the final two, rather than the first two, values in the
c     first dimension of basm that are nonzero at the right endpoint.
c     Thus, this case needs to be handled separately.
      do k = 1, npde
c        mm = nint*kcol*npde from the previous loop
         um(1, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 1, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   1, nint+1)
         um(2, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 2, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   2, nint+1)
      end do

c-----------------------------------------------------------------------
c     With the solution values at our interpolation points extracted
c     from the B-spline interpolant (i.e. the collocation solution,) and
c     the Hermite-Birkhoff coefficients precomputed, we can now evaluate
c     Hermite-Birkhoff interpolant at the quadrature points for an error
c     estimate. (Recall that the Hermite-Birkhoff interpolant uses
c     collocation solution values as its interpolation data.)
c
c     The values are returned in the usol array.
      do i = 1, nint
         hx = (x(i+1) - x(i))
         ii = (i-1) * (kcol-3)
         do j = 1, kcol+2
            jj = (i-1) * (kcol+2) + j
            do m = 1, npde

               usol(m, jj) =  h(1, j) * um(1, m, i)
     &                     +  h(2, j) * um(1, m, i+1)
     &                     + hd(1, j) * um(2, m, i)   * hx
     &                     + hd(2, j) * um(2, m, i+1) * hx

               do k = 1, kcol-3
                  usol(m, jj) = usol(m, jj) + g(k, j) * u(m, ii+k)
               end do

            end do
         end do
      end do

      return
      end

      subroutine lowcoeff(x, kcol, h, hd, g)
c-----------------------------------------------------------------------
c     This subroutine computes the Hermite-Birkhoff coefficients for a
c     given point at which the interpolant is to be evaluated.
c     They do not depend on the current subinterval.
c     This code uses a barycentric technique for phi, and exploits the
c     symmetry of the interpolation points.
c     Horner's rule is not used for phi, since that formulation lead to
c     some loss of precision.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c Input
      double precision        x
c     A point at which to evaluate the interpolant.
c     Typically a Gaussian quadrature point.
c
      integer                 kcol
c     The order of the polynomial is kcol+1.
c     Note that a value of kcol-1 is passed to this subroutine, however,
c     and therefore the interpolation points used here are not in fact
c     superconvergent (aside from the mesh points.) They would be
c     superconvergent for kcol-1, but are in effect arbitrary points
c     for kcol.
c
c-----------------------------------------------------------------------
c Output
      double precision        h(2)
      double precision        hd(2)
      double precision        g(kcol-2)
c     Hermite-Birkhoff coefficients for the given x.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        phi
c     phi(x) = product(x-w_i), i=1,kcol-2
c     where w_i are the interpolation points.
c     A barycentric technique is applied to replace
c     phi(x)_j with (phi(x))/(x-w_j).
c
      double precision        eta1, eta2, eta
c     These factors are common between the coefficients.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i
c
c-----------------------------------------------------------------------

c     precompute etas, which correspond to eta_1^2, eta_2^2 and eta^2
         eta1 = (x-1.d0)*(x-1.d0)
         eta2 = x*x
         eta = eta1*eta2

c     select which polynomial to evaluate for the H-B coefficients
      if     (kcol .eq. 2) then
         h(1) = (1d0+2d0*x)*eta1
         h(2) = (1d0-2d0*(x-1))*eta2
         hd(1) = x*eta1
         hd(2) = (x-1d0)*eta2

      elseif (kcol .eq. 3) then
         phi = x-0.5d0
         h(1) = -2.d0*(1.d0+4.d0*x)*eta1*phi
         h(2) = 2.d0*(5.d0-4.d0*x)*eta2*phi
         hd(1) = -2.d0*x*eta1*phi
         hd(2) = 2.d0*(x-1.d0)*eta2*phi
         g(1) = 16.d0*eta

      elseif (kcol .eq. 4) then
         phi = (x-0.3110177634953864d0)*(x-0.6889822365046136d0)
         h(1) = (0.4666666666666667d1+0.3111111111111111d2*x)*eta1*phi
         h(2) = (0.3577777777777778d2-0.3111111111111111d2*x)*eta2*phi
         hd(1) = 0.4666666666666667d1*x*eta1*phi
         hd(2) = (0.4666666666666667d1*x-0.4666666666666667d1)*eta2*phi
         g(1) = -0.5761858410762886d2*eta
         g(2) = -g(1)*(x-0.3110177634953864d0)
         g(1) = g(1)*(x-0.6889822365046136d0)

      elseif (kcol .eq. 5) then
         phi = (x-0.2113248654051871d0)*(x-0.5d0)
     #         *(x-0.7886751345948129d0)
         h(1) = (-0.12d2-0.12d3*x)*eta1*phi
         h(2) = (0.132d3-0.12d3*x)*eta2*phi
         hd(1) = -0.12d2*x*eta1*phi
         hd(2) = (0.12d2*x-0.12d2)*eta2*phi
         g(1) = 0.216d3*eta*phi
         g(2) = -192.d0*eta*phi/(x-0.5d0)
         g(3) = g(1)/(x-0.7886751345948129d0)
         g(1) = g(1)/(x-0.2113248654051871d0)

      elseif (kcol .eq. 6) then
         phi = (x-0.1526267046965671d0)*(x-0.3747185964571342d0)
     #         *(x-0.6252814035428658d0)*(x-0.8473732953034329d0)
         h(1) = (0.33d2+0.462d3*x)*eta1*phi
         h(2) = (0.495d3-0.462d3*x)*eta2*phi
         hd(1) = 0.33d2*x*eta1*phi
         hd(2) = (-0.33d2+0.33d2*x)*eta2*phi
         g(1) = -0.8197591840300222d3*eta*phi
         g(2) = 0.6925405260332655d3*eta*phi
         g(3) = -g(2)/(x-0.6252814035428658d0)
         g(4) = -g(1)/(x-0.8473732953034329d0)
         g(1) = g(1)/(x-0.1526267046965671d0)
         g(2) = g(2)/(x-0.3747185964571342d0)

      elseif (kcol .eq. 7) then
         phi = (x-0.1152723378341063d0)*(x-0.2895425974880943d0)
     #         *(x-0.5d0)*(x-0.7104574025119057d0)
     #         *(x-0.8847276621658937d0)
         h(1) = (-0.9533333333333333d2-0.1779555555555556d4*x)*eta1*phi
         h(2) = (0.1874888888888889d4-0.1779555555555556d4*x)*eta2*phi
         hd(1) = -0.9533333333333333d2*x*eta1*phi
         hd(2) = (-0.9533333333333333d2+0.9533333333333333d2*x)*eta2*phi
         g(1) = 0.3131255164938139d4*eta*phi
         g(2) = -0.257196627604925d4*eta*phi
         g(3) = 0.2440533333333333d4*eta*phi/(x-0.5d0)
         g(4) = g(2)/(x-0.7104574025119057d0)
         g(5) = g(1)/(x-0.8847276621658937d0)
         g(1) = g(1)/(x-0.1152723378341063d0)
         g(2) = g(2)/(x-0.2895425974880943d0)

      elseif (kcol .eq. 8) then
         phi = (x-0.9007700226825652d-1)*(x-0.2296976813063206d0)
     #         *(x-0.405661288754607d0)*(x-0.594338711245393d0)
     #         *(x-0.7703023186936794d0)*(x-0.9099229977317435d0)
         h(1) = (0.286d3+0.6864d4*x)*eta1*phi
         h(2) = (0.715d4-0.6864d4*x)*eta2*phi
         hd(1) = 0.286d3*x*eta1*phi
         hd(2) = (-0.286d3+0.286d3*x)*eta2*phi
         g(1) = -0.1201314415336355d5*eta*phi
         g(2) = 0.9696023778080566d4*eta*phi
         g(3) = -0.8929458911121293d4*eta*phi
         g(4) = -g(3)/(x-0.594338711245393d0)
         g(5) = -g(2)/(x-0.7703023186936794d0)
         g(6) = -g(1)/(x-0.9099229977317435d0)
         g(1) = g(1)/(x-0.9007700226825652d-1)
         g(2) = g(2)/(x-0.2296976813063206d0)
         g(3) = g(3)/(x-0.405661288754607d0)

      elseif (kcol .eq. 9) then
         phi = (x-0.7229898685756272d-1)*(x-0.1863109301186906d0)
     #         *(x-0.3341852231986051d0)*(x-0.5d0)*(x-0.66581477680139
     #         49d0)*(x-0.8136890698813094d0)*(x-0.9277010131424373d0)
         h(1) = (-0.884d3-0.2652d5*x)*eta1*phi
         h(2) = (0.27404d5-0.2652d5*x)*eta2*phi
         hd(1) = -0.884d3*x*eta1*phi
         hd(2) = (-0.884d3+0.884d3*x)*eta2*phi
         g(1) = 0.4624522420172525d5*eta*phi
         g(2) = -0.3688889968687672d5*eta*phi
         g(3) = 0.3332824691372289d5*eta*phi
         g(4) = -0.3232914285714286d5*eta*phi/(x-0.5d0)
         g(5) = g(3)/(x-0.6658147768013949d0)
         g(6) = g(2)/(x-0.8136890698813094d0)
         g(7) = g(1)/(x-0.9277010131424373d0)
         g(1) = g(1)/(x-0.7229898685756272d-1)
         g(2) = g(2)/(x-0.1863109301186906d0)
         g(3) = g(3)/(x-0.3341852231986051d0)

      else
         h(1) = 0.d0
         h(2) = 0.d0
         hd(1) = 0.d0
         hd(2) = 0.d0
         do i =1,7
           g(i) = 0.d0
         end do
      endif

      return
      end

      subroutine setpts(scpt, kcol, i, x, nint)
c-----------------------------------------------------------------------
c     This subroutine fills the scpt array with the superconvergent
c     points for the given mesh subinterval and kcol.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     The value of kcol for which superconvergent points are located.
c
      integer                 nint
c     Number of subintervals in the current mesh.
c
      double precision        x(nint+1)
c     The current mesh.
c
      integer                 i
c     Subinterval to which the superconvergent points for the given kcol
c     over [0,1] are mapped to.
c
c-----------------------------------------------------------------------
c Output
      double precision        scpt(mxkcol-2)
c     Superconvergent points for the given kcol mapped to the given
c     mesh subinterval.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        h
c     The width of the given mesh subinterval.
c
c-----------------------------------------------------------------------

      h = (x(i+1) - x(i))

      if     (kcol .eq. 3) then
         scpt(1) = x(i) + 0.5000000000000000d0 * h

      elseif (kcol .eq. 4) then
         scpt(1) = x(i) + 0.3110177634953864d0 * h
         scpt(2) = x(i) + 0.6889822365046136d0 * h

      elseif (kcol .eq. 5) then
         scpt(1) = x(i) + 0.2113248654051871d0 * h
         scpt(2) = x(i) + 0.5000000000000000d0 * h
         scpt(3) = x(i) + 0.7886751345948129d0 * h

      elseif (kcol .eq. 6) then
         scpt(1) = x(i) + 0.1526267046965671d0 * h
         scpt(2) = x(i) + 0.3747185964571342d0 * h
         scpt(3) = x(i) + 0.6252814035428658d0 * h
         scpt(4) = x(i) + 0.8473732953034329d0 * h

      elseif (kcol .eq. 7) then
         scpt(1) = x(i) + 0.1152723378341063d0 * h
         scpt(2) = x(i) + 0.2895425974880943d0 * h
         scpt(3) = x(i) + 0.5000000000000000d0 * h
         scpt(4) = x(i) + 0.7104574025119057d0 * h
         scpt(5) = x(i) + 0.8847276621658937d0 * h

      elseif (kcol .eq. 8) then
         scpt(1) = x(i) + 0.9007700226825652d-1* h
         scpt(2) = x(i) + 0.2296976813063206d0 * h
         scpt(3) = x(i) + 0.4056612887546070d0 * h
         scpt(4) = x(i) + 0.5943387112453930d0 * h
         scpt(5) = x(i) + 0.7703023186936794d0 * h
         scpt(6) = x(i) + 0.9099229977317435d0 * h

      elseif (kcol .eq. 9) then
         scpt(1) = x(i) + 0.7229898685756272d-1* h
         scpt(2) = x(i) + 0.1863109301186906d0 * h
         scpt(3) = x(i) + 0.3341852231986051d0 * h
         scpt(4) = x(i) + 0.5000000000000000d0 * h
         scpt(5) = x(i) + 0.6658147768013949d0 * h
         scpt(6) = x(i) + 0.8136890698813094d0 * h
         scpt(7) = x(i) + 0.9277010131424373d0 * h

      elseif (kcol .eq. 10) then
         scpt(1) = x(i) + 0.5929571219129399d-1* h
         scpt(2) = x(i) + 0.1539696908715823d0 * h
         scpt(3) = x(i) + 0.2792835119457421d0 * h
         scpt(4) = x(i) + 0.4241841678533667d0 * h
         scpt(5) = x(i) + 0.5758158321466333d0 * h
         scpt(6) = x(i) + 0.7207164880542579d0 * h
         scpt(7) = x(i) + 0.8460303091284177d0 * h
         scpt(8) = x(i) + 0.9407042878087060d0 * h

      endif

      return
      end

      subroutine scint(kcol, nint, npde, istart, icount, xbs,
     &                 x, y, h, hd, g, bassc, basm, u, um, usol)
c-----------------------------------------------------------------------
c     This subroutine evaluates a Hermite-Birkhoff (H-B) interpolant of
c     order kcol+2 (degree kcol+1) at kcol+3 Gaussian quadrature points
c     on the interval (0,1). (These values will be used elsewhere to
c     obtain an L2-norm error estimate.)
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     bacoli can have up to 10 collocation points per subinterval.
c
      integer                 nconti
      parameter              (nconti = 2)
c     nconti continuity conditions are imposed at the internal mesh
c     points of the collocation solution.
c
      double precision        zero
      parameter              (zero = 0.0d0)
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     Number of collocation points per subinterval for the collocation
c     solution. The collocation solution is a piecewise polynomial of
c     order kcol+1 (degree kcol.) The values obtained from this
c     interpolant are of order kcol+2. The interpolant is a piecewise
c     polynomial of degree kcol+1.
c
      integer                 nint
c     Number of subintervals in the current mesh.
c
      integer                 npde
c     Number of PDEs in the system being solved.
c
      integer                 istart
c     istart indicates whether or not the integration has just begun.
c         istart = 0, it is the initial step;
c                = 1, it is not the initial step.
c
      integer                 icount
c     icount is the number of remeshings that have been performed at
c     the current time step.
c     istart and icount are used to determine whether certain
c     coefficients need to be calculated, depending upon whether the
c     mesh has changed.
c
      double precision        xbs((nint+1)*kcol+nconti*2)
c     The breakpoint sequence. This is required by BSPLVD to compute the
c     values of the B-spline basis functions.
c
      double precision        x(nint+1)
c     The current mesh point sequence.
c
      double precision        y(npde*(nint*kcol+nconti))
c     B-spline coefficients computed by DASSL.
c
c-----------------------------------------------------------------------
c Inoutput
      double precision        h(2, kcol+3, nint)
      double precision        hd(2, kcol+3, nint)
      double precision        g(kcol, kcol+3, nint)
c     Coefficients of the Hermite-Birkhoff interpolant.
c     h and hd coefficients multiply the values and derivatives of the
c     collocation solution at mesh points, g coefficients multiply the
c     values of the collocation solution at the internal points.
c     These are stored between remeshings and depend upon mesh ratios
c     in the SCI case, since points external to the subinterval are
c     used. Thus, nint sets of these coefficients must be kept.
c
      double precision        bassc(kcol+nconti, (kcol-2)*nint)
c     Values of all nonzero B-spline basis functions at the
c     superconvergent points inside the subintervals.
c     These are stored between remeshings.
c
      double precision        basm(kcol+nconti, 2, nint+1)
c     Values of all nonzero B-spline basis functions and their first
c     derivatives at the mesh points. kcol of every kcol+nconti of these
c     are actually zero, but the trouble to reduce the storage used by
c     this vector seemed too great. These are stored between remeshings.
c
c-----------------------------------------------------------------------
c Output
      double precision        usol(npde, nint*(kcol+3))
c     Array of evaluations of the Hermite-Birkhoff interpolant at the
c     quadrature points used for the computation of an L2-norm error
c     estimate in elsewhere.
c
c-----------------------------------------------------------------------
c Work storage
      double precision        u(npde, (kcol-2)*nint)
c     Work array for storage of values of the collocation solution at
c     the internal points used in the Hermite-Birkhoff interpolant.
c
      double precision        um(2, npde, nint+1)
c     Work array for storage of values of the collocation solution, and
c     its first derivative, at the mesh points used in the Hermite-
c     Birkhoff interpolant.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        scpt(mxkcol-2)
c     The kcol interpolated superconvergent points internal to each
c     subinterval. Some are used in more than one subinterval of the
c     Hermite-Birkhoff interpolant.
c
      double precision        theta(mxkcol+3)
c     Quadrature points for a subinterval, which are the Gauss points
c     on (0,1).
c
      double precision        gwts(mxkcol+3)
c     Gaussian weights. Unused here.
c
      double precision        gwork((mxkcol+3)*(mxkcol+3))
c     Work storage for gauleg.
c
      integer                 ileft
c     Breakpoint information. ileft is used in conjunction with the
c     breakpoint sequence for the evaluation of B-spline basis functions
c     by BSPLVD.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j, k, m, jj, mm, ii
c
c-----------------------------------------------------------------------
c Subroutines Called
c                             bsplvd
c                             gauleg
c                             scicoeff
c                             setpts
c
c-----------------------------------------------------------------------

c     Check whether the H-B coefficients and B-spline basis function
c     evaluations need to be computed.
c     Previous values may be reused if the mesh has not changed
c     since the last step.
      if ((icount .ne. 0) .or. (istart .eq. 0)) then

c        Set thetas (weights will not be computed since last parameter
c        in the call to gauleg is set to 2.)
         call gauleg(kcol+3, (kcol+3)*(kcol+3), theta, gwts, gwork, 2)

c        H-B coefficients are computed for the Gaussian quadrature
c        points at which the H-B interpolant is evaluated.
c        Since two points external to the subinterval are used in the
c        interpolant piece for each subinterval, all the coefficients
c        end up depending upon mesh subinteval ratios.
         call scicoeff(theta, kcol, nint, x, h, hd, g)

c        Call BSPLVD to obtain values of the B-spline basis functions at
c        the internal superconvergent points used in the interpolant.
c        The subroutine INTERV does not need to be called to find ileft,
c        since we know which subinterval we are in already.
         ileft = nconti
         do i = 1, nint
            call setpts(scpt, kcol, i, x, nint)
            ileft = ileft + kcol
            do j = 1, kcol-2
               jj = (i-1) * (kcol-2) + j
               call bsplvd(xbs, kcol+nconti, scpt(j), ileft,
     &                     bassc(1, jj), 1)
            end do
         end do

c        Call BSPLVD for the values and first derivatives of the
c        B-spline basis functions at the nint+1 mesh points.
c       (ileft does not follow the pattern at the rightmost mesh point,
c        so we need a separate call for the last mesh point.)
         ileft = nconti
         do i = 1, nint
            ileft = ileft + kcol
            call bsplvd(xbs, kcol+nconti, x(i), ileft,
     &                  basm(1, 1, i), 2)
         end do
c        ileft = kcol+nconti+kcol*(nint-1)
         call bsplvd(xbs, kcol+nconti, x(nint+1), ileft,
     &               basm(1, 1, nint+1), 2)

      endif

c-----------------------------------------------------------------------
c     Compute the values of the collocation solution at the
c     superconvergent and mesh points by multiplying the values of the
c     basis functions at those points by the B-spline coefficients
c     computed by DASSL.
c
c     The end result of these two steps (evaluating the basis functions
c     and later multiplying by the coefficients) is the same as calling
c     the subroutine VALUES, but this approach avoids a large number
c     of unnecessary calls to BSPLVD by saving some of the required
c     values between remeshings.
      ii = 0
      do i = 1, nint
         do j = 1, kcol-2
            jj = (i-1) * (kcol-2) + j
            do k = 1, npde
               u(k, jj) = zero
               mm = ii
               do m = 1, kcol+nconti
c                 mm = (m + (i-1) * kcol - 1) * npde
                  u(k, jj) = u(k, jj) + y(mm+k) * bassc(m, jj)
                  mm = mm + npde
               end do
            end do
         end do
         ii = ii + kcol * npde
      end do
c     All but two of the values of the basis functions at the mesh
c     points are zero, so um gets special treatment.
      mm = 0
      do i = 1, nint
         do k = 1, npde
            um(1, k, i) = y(mm+k) * basm(1, 1, i)
     &             + y(npde+mm+k) * basm(2, 1, i)
            um(2, k, i) = y(mm+k) * basm(1, 2, i)
     &             + y(npde+mm+k) * basm(2, 2, i)
         end do
         mm = mm + kcol * npde
      end do
c     It is the final two, rather than the first two, values in the
c     first dimension of basm that are nonzero at the right endpoint.
c     This case needs to be handled separately.
      do k = 1, npde
c        mm = nint*kcol*npde from the previous loop
         um(1, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 1, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   1, nint+1)
         um(2, k, nint+1) = y(mm+k) * basm(kcol+nconti-1, 2, nint+1)
     &               + y(npde+mm+k) * basm(kcol+nconti,   2, nint+1)
      end do

c-----------------------------------------------------------------------
c     With the solution values at our interpolation points extracted
c     from the B-spline interpolant (i.e. the collocation solution,) and
c     the Hermite-Birkhoff coefficients precomputed, we can now evaluate
c     Hermite-Birkhoff interpolant at the quadrature points for an error
c     estimate. (Recall that the Hermite-Birkhoff interpolant uses
c     collocation solution values as its interpolation data.)
c
c     The values are returned in the usol array.
      if (kcol .ne. 3) then
         do i = 1, nint
            if (i .ne. 1 .and. i .ne. nint) then
               ii = (i-1) * (kcol-2) - 1
            elseif (i .eq. 1) then
               ii = (i-1) * (kcol-2)
            else
               ii = (i-1) * (kcol-2) - 2
            endif
            do j = 1, kcol+3
               jj = (i-1) * (kcol+3) + j
               do m = 1, npde
                  usol(m, jj) =  h(1, j, i) * um(1, m, i)
     &                        +  h(2, j, i) * um(1, m, i+1)
     &                        + hd(1, j, i) * um(2, m, i)
     &                        + hd(2, j, i) * um(2, m, i+1)
                  do k = 1, kcol
                     usol(m, jj) = usol(m, jj)
     &                           + g(k, j, i) * u(m, ii + k)
                  end do
               end do
            end do
         end do
c     This has a different form with kcol=3, which uses a third mesh
c     point when interpolating at boundary subintervals.
      else
         do i = 1, nint
            do j = 1, kcol+3
               jj = (i-1) * (kcol+3) + j
               do m = 1, npde
                  if (i .ne. 1 .and. i .ne. nint) then
                     usol(m, jj) = g(1, j, i) * u(m, i-1)
     &                           + g(2, j, i) * u(m, i)
     &                           + g(3, j, i) * u(m, i+1)
                  elseif (i .eq. 1) then
                     usol(m, jj) = g(1, j, 1) * u(m, 1)
     &                           + g(2, j, 1) * u(m, 2)
     &                           + g(3, j, 1) * um(1, m, 3)
                  else
                     usol(m, jj) = g(1, j, nint) * um(1, m, nint-1)
     &                           + g(2, j, nint) * u(m, nint-1)
     &                           + g(3, j, nint) * u(m, nint)
                  endif
                  usol(m, jj) = usol(m, jj)
     &                        +  h(1, j, i) * um(1, m, i)
     &                        +  h(2, j, i) * um(1, m, i+1)
     &                        + hd(1, j, i) * um(2, m, i)
     &                        + hd(2, j, i) * um(2, m, i+1)
               end do
            end do
         end do
      endif

      return
      end

      subroutine scicoeff(x, kcol, nint, mesh, h, hd, g)
c-----------------------------------------------------------------------
c     This subroutine computes the Hermite-Birkhoff coefficients for the
c     kcol+3 quadrature points per subinterval at which the interpolant
c     is to be evaluated. They differ between subintervals since points
c     from adjacent subintervals are used to construct the interpolant.
c
c     Note that the hd coefficients are also scaled by the subinterval
c     width here.
c
c     Last modified by Jack Pew, September 1, 2011.
c
c-----------------------------------------------------------------------
      implicit none
c Constants
      integer                 mxkcol
      parameter              (mxkcol = 10)
c     Maximum value kcol may take.
c
c-----------------------------------------------------------------------
c Input
      integer                 kcol
c     The value of kcol for which superconvergent points are located.
c
      double precision        x(kcol+3)
c     Points relative to each subinterval at which to evaluate the
c     interpolant.
c
      integer                 nint
c     Number of mesh subintervals.
c
      double precision        mesh(nint+1)
c     The mesh array.
c
c-----------------------------------------------------------------------
c Output
      double precision        h(2, kcol+3, nint)
      double precision        hd(2, kcol+3, nint)
      double precision        g(kcol, kcol+3, nint)
c     H-B coefficients evaluated for the current mesh at the points x
c     on each subinterval.
c
c-----------------------------------------------------------------------
c Local variables
      double precision        phi
c     phi(x) = product(x-w_i), i=1,kcol
c     where w_i are the superconvergent points.
c     A barycentric technique is applied to replace
c     phi(x)_j with (phi(x))/(x-w_j).
c
      double precision        a, aa, b, bb, hx
c     Superconvergent points from adjacent intervals,
c     and the subinterval size.
c     a and aa are the nearest two such points from the left,
c     b and bb are the nearest two such points from the right,
c     and hx is the current subinterval size.
c     While hx is not strictly a part of any coefficient, it is
c     multiplied here to save a little work.
c
      double precision        eta(3, mxkcol+3), gamma1, gamma2
c     Other precomputed factors.
c
c-----------------------------------------------------------------------
c Loop indices
      integer                 i, j
c
c-----------------------------------------------------------------------

c     Precompute etas, which correspond to eta_1^2, eta_2^2 and eta^2
      do j = 1, kcol+3
         eta(1, j) = (x(j)-1.d0) * (x(j)-1.d0)
         eta(2, j) = x(j) * x(j)
         eta(3, j) = eta(1, j) * eta(2, j)
      end do

c     Evaluate the coefficients for a given kcol and mesh at the
c     kcol+3 points in x.
c     This is greatly complicated over the LOI case by the use of points
c     external to a subinterval, bringing mesh ratios into the mix, and
c     necessitating a distinction between internal and boundary
c     subintervals.
      if (kcol .eq. 3) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.5d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-4.d0-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+4.d0+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -2.d0/a/b
               h(2,kcol+3,i) = 2.d0/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.5d0)*(a-b)*a*a
     #            *(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 16.d0/(0.5d0-a)/(0.5d0-b)
               g(3,kcol+3,i) = 1.d0/((b-a)*(b-0.5d0)*b*b*(b-1.d0)
     #            *(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.5d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -4.d0-1.d0/b-1.d0/bb
               gamma2 = 4.d0+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -2.d0/b/bb
               h(2,kcol+3,i) = 2.d0/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 16.d0/(0.5d0-bb)/(0.5d0-b)
               g(2,kcol+3,i) = 1.d0/((b-bb)*(b-0.5d0)*b*b*(b-1.d0)
     #            *(b-1.d0))
               g(3,kcol+3,i) = 1.d0/((bb-b)*(bb-0.5d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.5d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -(mesh(i)-mesh(i-1))/hx
               a = -0.5d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-4.d0
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+4.d0
               h(1,kcol+3,i) = -2.d0/aa/a
               h(2,kcol+3,i) = 2.d0/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.5d0)*(aa-a)*aa*aa
     #            *(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.5d0)*(a-aa)*a*a
     #            *(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 16.d0/(0.5d0-a)/(0.5d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.5d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 4) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.3110177634953864d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.3110177634953864d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.6666666666666667d1-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.6666666666666667d1+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = 0.4666666666666667d1/a/b
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/(a-0.3110177634953864d0)/(a-0.6889
     #            822365046136d0)/(a-b)/a/a/(a-1.d0)/(a-1.d0)
               g(2,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-a)/(0.3110177634953864d0-b)
               g(3,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-a)/(0.6889822365046136d0-b)
               g(4,kcol+3,i) = 1.d0/(b-a)/(b-0.3110177634953864d0)
     #            /(b-0.6889822365046136d0)/b/b/(b-1.d0)/(b-1.d0)

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.3110177634953864d0)*(x(j)
     #                  -0.6889822365046136d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.3110177634953864d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.6889822365046136d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.6666666666666667d1-1.d0/b-1.d0/bb
               gamma2 = 0.6666666666666667d1+1.d0/(1.d0-b)+1.d0
     #            /(1.d0-bb)
               h(1,kcol+3,i) = 0.4666666666666667d1/b/bb
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-b)/(0.3110177634953864d0-bb)
               g(2,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-b)/(0.6889822365046136d0-bb)
               g(3,kcol+3,i) = 1.d0/(b-bb)/(b-0.3110177634953864d0)
     #            /(b-0.6889822365046136d0)/b/b/(b-1.d0)/(b-1.d0)
               g(4,kcol+3,i) = 1.d0/(bb-b)/(bb-0.3110177634953864d0)
     #            /(bb-0.6889822365046136d0)/bb/bb/(bb-1.d0)/(bb-1.d0)

               do j = 1, kcol+3
                  phi = (x(j)-0.3110177634953864d0)*(x(j)
     #                  -0.6889822365046136d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.6889822365046136d0*(mesh(i)-mesh(i-1))/hx
               a = -0.3110177634953864d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.6666666666666667d1
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.6666666666666667d1
               h(1,kcol+3,i) = 0.4666666666666667d1/aa/a
               h(2,kcol+3,i) = 0.4666666666666667d1/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/(aa-a)/(aa-0.3110177634953864d0)
     #            /(aa-0.6889822365046136d0)/aa/aa/(aa-1.d0)/(aa-1.d0)
               g(2,kcol+3,i) = 1.d0/(a-aa)/(a-0.3110177634953864d0)
     #            /(a-0.6889822365046136d0)/a/a/(a-1.d0)/(a-1.d0)
               g(3,kcol+3,i) = -0.5761858410762886d2/(0.3110177634953
     #            864d0-a)/(0.3110177634953864d0-aa)
               g(4,kcol+3,i) = 0.5761858410762886d2/(0.68898223650461
     #            36d0-a)/(0.6889822365046136d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.3110177634953864d0)
     #                  *(x(j)-0.6889822365046136d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3110177634953864d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6889822365046136d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 5) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.2113248654051871d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.2113248654051871d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.1d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.1d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -0.12d2/a/b
               h(2,kcol+3,i) = 0.12d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.2113248654051871d0)*(a-0.5d0)
     #            *(a-0.7886751345948129d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.216d3/(0.2113248654051871d0-a)
     #            /(0.2113248654051871d0-b)
               g(3,kcol+3,i) = -192.d0/(0.5d0-a)/(0.5d0-b)
               g(4,kcol+3,i) = 0.216d3/(0.7886751345948129d0-a)
     #            /(0.7886751345948129d0-b)
               g(5,kcol+3,i) = 1.d0/((b-a)*(b-0.2113248654051871d0)*(b-
     #            0.5d0)*(b-0.7886751345948129d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.2113248654051871d0)*(x(j)
     #                  -0.5d0)*(x(j)-0.7886751345948129d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.2113248654051871d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.5d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.1d2-1.d0/b-1.d0/bb
               gamma2 = 0.1d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.12d2/b/bb
               h(2,kcol+3,i) = 0.12d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.216d3/(0.2113248654051871d0-b)
     #            /(0.2113248654051871d0-bb)
               g(2,kcol+3,i) = -192.d0/(0.5d0-b)/(0.5d0-bb)
               g(3,kcol+3,i) = 0.216d3/(0.7886751345948129d0-b)
     #            /(0.7886751345948129d0-bb)
               g(4,kcol+3,i) = 1.d0/((b-0.2113248654051871d0)*(b-0.5d0)*
     #            (b-0.7886751345948129d0)*(b-bb)*b*b*(b-1.d0)*(b-1.d0))
               g(5,kcol+3,i) = 1.d0/((bb-0.2113248654051871d0)*(bb
     #            -0.5d0)*(bb-0.7886751345948129d0)*(bb-b)*bb*bb
     #            *(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.2113248654051871d0)*(x(j)-0.5d0)*(x(j)
     #                  -0.7886751345948129d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.5d0*(mesh(i)-mesh(i-1))/hx
               a = -0.2113248654051871d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.1d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.1d2
               h(1,kcol+3,i) = -0.12d2/aa/a
               h(2,kcol+3,i) = 0.12d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.2113248654051871d0)*(aa
     #            -0.5d0)*(aa-0.7886751345948129d0)*(aa-a)*aa*aa
     #            *(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.2113248654051871d0)*(a-0.5d0)*
     #            (a-0.7886751345948129d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.216d3/(0.2113248654051871d0-a)
     #            /(0.2113248654051871d0-aa)
               g(4,kcol+3,i) = -192.d0/(0.5d0-a)/(0.5d0-aa)
               g(5,kcol+3,i) = 0.216d3/(0.7886751345948129d0-a)
     #            /(0.7886751345948129d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.2113248654051871d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.7886751345948129d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2113248654051871d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7886751345948129d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 6) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.1526267046965671d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.1526267046965671d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.14d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.14d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = 0.33d2/a/b
               h(2,kcol+3,i) = 0.33d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.1526267046965671d0)*(a-0.37471
     #            85964571342d0)*(a-0.6252814035428658d0)*(a-0.84737329
     #            53034329d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-a)/(0.1526267046965671d0-b)
               g(3,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-a)/(0.3747185964571342d0-b)
               g(4,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-a)/(0.6252814035428658d0-b)
               g(5,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-a)/(0.8473732953034329d0-b)
               g(6,kcol+3,i) = 1.d0/((b-a)*(b-0.1526267046965671d0)*(b
     #            -0.3747185964571342d0)*(b-0.6252814035428658d0)*(b
     #            -0.8473732953034329d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.1526267046965671d0)*(x(j)-0.3
     #                  747185964571342d0)*(x(j)-0.6252814035428658d0)
     #                  *(x(j)-0.8473732953034329d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.1526267046965671d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.3747185964571342d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.14d2-1.d0/b-1.d0/bb
               gamma2 = 0.14d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.33d2/b/bb
               h(2,kcol+3,i) = 0.33d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-b)/(0.1526267046965671d0-bb)
               g(2,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-b)/(0.3747185964571342d0-bb)
               g(3,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-b)/(0.6252814035428658d0-bb)
               g(4,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-b)/(0.8473732953034329d0-bb)
               g(5,kcol+3,i) = 1.d0/((b-bb)*(b-0.1526267046965671d0)*(b
     #            -0.3747185964571342d0)*(b-0.6252814035428658d0)*(b
     #            -0.8473732953034329d0)*b*b*(b-1.d0)*(b-1.d0))
               g(6,kcol+3,i) = 1.d0/((bb-b)*(bb-0.1526267046965671d0)*
     #            (bb-0.3747185964571342d0)*(bb-0.6252814035428658d0)*
     #            (bb-0.8473732953034329d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.1526267046965671d0)*(x(j)-0.374718
     #                  5964571342d0)*(x(j)-0.6252814035428658d0)
     #                  *(x(j)-0.8473732953034329d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.3747185964571342d0*(mesh(i)-mesh(i-1))/hx
               a = -0.1526267046965671d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.14d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.14d2
               h(1,kcol+3,i) = 0.33d2/aa/a
               h(2,kcol+3,i) = 0.33d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.1526267046965671d0)*(aa-0.374
     #            7185964571342d0)*(aa-0.6252814035428658d0)*(aa-0.8473
     #            732953034329d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.1526267046965671d0)*(a-0.37471
     #            85964571342d0)*(a-0.6252814035428658d0)*(a-0.84737329
     #            53034329d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.8197591840300222d3/(0.152626704696567
     #            1d0-a)/(0.1526267046965671d0-aa)
               g(4,kcol+3,i) = 0.6925405260332655d3/(0.37471859645713
     #            42d0-a)/(0.3747185964571342d0-aa)
               g(5,kcol+3,i) = -0.6925405260332655d3/(0.625281403542865
     #            8d0-a)/(0.6252814035428658d0-aa)
               g(6,kcol+3,i) = 0.8197591840300222d3/(0.847373295303432
     #            9d0-a)/(0.8473732953034329d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.1526267046965671d0)
     #                  *(x(j)-0.3747185964571342d0)*(x(j)-0.625281403
     #                  5428658d0)*(x(j)-0.8473732953034329d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1526267046965671d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3747185964571342d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6252814035428658d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8473732953034329d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 7) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.1152723378341063d0*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.1152723378341063d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.1866666666666667d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.1866666666666667d2+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = -0.9533333333333333d2/a/b
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.1152723378341063d0)*(a-0.28954
     #            25974880943d0)*(a-0.5d0)*(a-0.7104574025119057d0)
     #            *(a-0.8847276621658937d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-a)/(0.1152723378341063d0-b)
               g(3,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-a)/(0.2895425974880943d0-b)
               g(4,kcol+3,i) = 0.2440533333333333d4/(0.5d0-a)/(0.5d0-b)
               g(5,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-a)/(0.7104574025119057d0-b)
               g(6,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-a)/(0.8847276621658937d0-b)
               g(7,kcol+3,i) = 1.d0/((b-a)*(b-0.1152723378341063d0)*(b-
     #            0.2895425974880943d0)*(b-0.5d0)*(b-0.710457402511905
     #            7d0)*(b-0.8847276621658937d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.1152723378341063d0)*(x(j)-0.28
     #                  95425974880943d0)*(x(j)-0.5d0)*(x(j)-0.71045740
     #                  25119057d0)*(x(j)-0.8847276621658937d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.1152723378341063d0*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.2895425974880943d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.1866666666666667d2-1.d0/b-1.d0/bb
               gamma2 = 0.1866666666666667d2+1.d0/(1.d0-b)
     #            +1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.9533333333333333d2/b/bb
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-bb)/(0.1152723378341063d0-b)
               g(2,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-bb)/(0.2895425974880943d0-b)
               g(3,kcol+3,i) = 0.2440533333333333d4/(0.5d0-bb)/(0.5d0-b)
               g(4,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-bb)/(0.7104574025119057d0-b)
               g(5,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-bb)/(0.8847276621658937d0-b)
               g(6,kcol+3,i) = 1.d0/((b-bb)*(b-0.1152723378341063d0)*(b-
     #            0.2895425974880943d0)*(b-0.5d0)*(b-0.710457402511905
     #            7d0)*(b-0.8847276621658937d0)*b*b*(b-1.d0)*(b-1.d0))
               g(7,kcol+3,i) = 1.d0/((bb-b)*(bb-0.1152723378341063d0)
     #            *(bb-0.2895425974880943d0)*(bb-0.5d0)*(bb-0.71045740
     #            25119057d0)*(bb-0.8847276621658937d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.1152723378341063d0)*(x(j)-0.289542597488
     #                  0943d0)*(x(j)-0.5d0)*(x(j)-0.7104574025119057d0)
     #                  *(x(j)-0.8847276621658937d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.2895425974880943d0*(mesh(i)-mesh(i-1))/hx
               a = -0.1152723378341063d0*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.1866666666666667d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.1866666666666667d2
               h(1,kcol+3,i) = -0.9533333333333333d2/aa/a
               h(2,kcol+3,i) = 0.9533333333333333d2/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.1152723378341063d0)*(aa-0.28
     #            95425974880943d0)*(aa-0.5d0)*(aa-0.7104574025119057d0)
     #            *(aa-0.8847276621658937d0)*(aa-a)*aa*aa*(aa-1.d0)
     #            *(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.1152723378341063d0)*(a-0.28954
     #            25974880943d0)*(a-0.5d0)*(a-0.7104574025119057d0)*(a
     #            -0.8847276621658937d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.3131255164938139d4/(0.115272337834106
     #            3d0-a)/(0.1152723378341063d0-aa)
               g(4,kcol+3,i) = -0.257196627604925d4/(0.289542597488094
     #            3d0-a)/(0.2895425974880943d0-aa)
               g(5,kcol+3,i) = 0.2440533333333333d4/(0.5d0-a)/(0.5d0-aa)
               g(6,kcol+3,i) = -0.257196627604925d4/(0.710457402511905
     #            7d0-a)/(0.7104574025119057d0-aa)
               g(7,kcol+3,i) = 0.3131255164938139d4/(0.884727662165893
     #            7d0-a)/(0.8847276621658937d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.1152723378341063d0)*(
     #                  x(j)-0.2895425974880943d0)*(x(j)-0.5d0)*(x(j)-0.
     #                  7104574025119057d0)*(x(j)-0.8847276621658937d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1152723378341063d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2895425974880943d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7104574025119057d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8847276621658937d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 8) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.9007700226825652d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.9007700226825652d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.24d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.24d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = 0.286d3/a/b
               h(2,kcol+3,i) = 0.286d3/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.9007700226825652d-1)*(a-0.229
     #            6976813063206d0)*(a-0.405661288754607d0)*(a-0.5943387
     #            11245393d0)*(a-0.7703023186936794d0)*(a-0.90992299773
     #            17435d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-a)/(0.9007700226825652d-1-b)
               g(3,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-a)/(0.2296976813063206d0-b)
               g(4,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-a)/(0.405661288754607d0-b)
               g(5,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -a)/(0.594338711245393d0-b)
               g(6,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-a)/(0.7703023186936794d0-b)
               g(7,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-a)/(0.9099229977317435d0-b)
               g(8,kcol+3,i) = 1.d0/((b-a)*(b-0.9007700226825652d-1)*(b
     #            -0.2296976813063206d0)*(b-0.405661288754607d0)*(b-0.5
     #            94338711245393d0)*(b-0.7703023186936794d0)*(b-0.90992
     #            29977317435d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.9007700226825652d-1)*(x(j)-0.2
     #                  296976813063206d0)*(x(j)-0.405661288754607d0)
     #                  *(x(j)-0.594338711245393d0)*(x(j)-0.7703023186
     #                  936794d0)*(x(j)-0.9099229977317435d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.9007700226825652d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.2296976813063206d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.24d2-1.d0/b-1.d0/bb
               gamma2 = 0.24d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.286d3/b/bb
               h(2,kcol+3,i) = 0.286d3/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-bb)/(0.9007700226825652d-1-b)
               g(2,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-bb)/(0.2296976813063206d0-b)
               g(3,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-bb)/(0.405661288754607d0-b)
               g(4,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -bb)/(0.594338711245393d0-b)
               g(5,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-bb)/(0.7703023186936794d0-b)
               g(6,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-bb)/(0.9099229977317435d0-b)
               g(7,kcol+3,i) = 1.d0/((b-bb)*(b-0.9007700226825652d-1)*(b
     #            -0.2296976813063206d0)*(b-0.405661288754607d0)*(b-0.5
     #            94338711245393d0)*(b-0.7703023186936794d0)*(b-0.90992
     #            29977317435d0)*b*b*(b-1.d0)*(b-1.d0))
               g(8,kcol+3,i) = 1.d0/((bb-b)*(bb-0.9007700226825652d-1)*
     #            (bb-0.2296976813063206d0)*(bb-0.405661288754607d0)*(bb
     #            -0.594338711245393d0)*(bb-0.7703023186936794d0)*(bb-0.
     #            9099229977317435d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.9007700226825652d-1)*(x(j)-0.2296976813
     #                  063206d0)*(x(j)-0.405661288754607d0)*(x(j)-0.59
     #                  4338711245393d0)*(x(j)-0.7703023186936794d0)
     #                  *(x(j)-0.9099229977317435d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.2296976813063206d0*(mesh(i)-mesh(i-1))/hx
               a = -0.9007700226825652d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.24d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.24d2
               h(1,kcol+3,i) = 0.286d3/aa/a
               h(2,kcol+3,i) = 0.286d3/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.9007700226825652d-1)*(aa-0.2
     #            296976813063206d0)*(aa-0.405661288754607d0)*(aa-0.594
     #            338711245393d0)*(aa-0.7703023186936794d0)*(aa-0.90992
     #            29977317435d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.9007700226825652d-1)*(a-0.229
     #            6976813063206d0)*(a-0.405661288754607d0)*(a-0.5943387
     #            11245393d0)*(a-0.7703023186936794d0)*(a-0.90992299773
     #            17435d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.1201314415336355d5/(0.900770022682565
     #            2d-1-a)/(0.9007700226825652d-1-aa)
               g(4,kcol+3,i) = 0.9696023778080566d4/(0.229697681306320
     #            6d0-a)/(0.2296976813063206d0-aa)
               g(5,kcol+3,i) = -0.8929458911121293d4/(0.40566128875460
     #            7d0-a)/(0.405661288754607d0-aa)
               g(6,kcol+3,i) = 0.8929458911121293d4/(0.594338711245393d0
     #            -a)/(0.594338711245393d0-aa)
               g(7,kcol+3,i) = -0.9696023778080565d4/(0.770302318693679
     #            4d0-a)/(0.7703023186936794d0-aa)
               g(8,kcol+3,i) = 0.1201314415336355d5/(0.909922997731743
     #            5d0-a)/(0.9099229977317435d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.9007700226825652d-1)
     #                  *(x(j)-0.2296976813063206d0)*(x(j)-0.4056612887
     #                  54607d0)*(x(j)-0.594338711245393d0)*(x(j)-0.770
     #                  3023186936794d0)*(x(j)-0.9099229977317435d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9007700226825652d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2296976813063206d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.405661288754607d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.594338711245393d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7703023186936794d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9099229977317435d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 9) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.7229898685756272d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.7229898685756272d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.3d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.3d2+1.d0/(1.d0-b)
               h(1,kcol+3,i) = -0.884d3/a/b
               h(2,kcol+3,i) = 0.884d3/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.7229898685756272d-1)*(a-0.186
     #            3109301186906d0)*(a-0.3341852231986051d0)*(a-0.5d0)*
     #            (a-0.6658147768013949d0)*(a-0.8136890698813094d0)*(a
     #            -0.9277010131424373d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-a)/(0.7229898685756272d-1-b)
               g(3,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-a)/(0.1863109301186906d0-b)
               g(4,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-a)/(0.3341852231986051d0-b)
               g(5,kcol+3,i) = -0.3232914285714286d5/(0.5d0-a)/(0.5d0-b)
               g(6,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-a)/(0.6658147768013949d0-b)
               g(7,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-a)/(0.8136890698813094d0-b)
               g(8,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-a)/(0.9277010131424373d0-b)
               g(9,kcol+3,i) = 1.d0/((b-a)*(b-0.7229898685756272d-1)*(b
     #            -0.1863109301186906d0)*(b-0.3341852231986051d0)*(b-0.
     #            5d0)*(b-0.6658147768013949d0)*(b-0.8136890698813094d0)
     #            *(b-0.9277010131424373d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.7229898685756272d-1)*(x(j)-0.1
     #                  863109301186906d0)*(x(j)-0.3341852231986051d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.6658147768013949d0)*(x(j)
     #                  -0.8136890698813094d0)*(x(j)-0.927701013142437
     #                  3d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.7229898685756272d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.1863109301186906d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.3d2-1.d0/b-1.d0/bb
               gamma2 = 0.3d2+1.d0/(1.d0-b)+1.d0/(1.d0-bb)
               h(1,kcol+3,i) = -0.884d3/b/bb
               h(2,kcol+3,i) = 0.884d3/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-bb)/(0.7229898685756272d-1-b)
               g(2,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-bb)/(0.1863109301186906d0-b)
               g(3,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-bb)/(0.3341852231986051d0-b)
               g(4,kcol+3,i) = -0.3232914285714286d5/(0.5d0-bb)
     #            /(0.5d0-b)
               g(5,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-bb)/(0.6658147768013949d0-b)
               g(6,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-bb)/(0.8136890698813094d0-b)
               g(7,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-bb)/(0.9277010131424373d0-b)
               g(8,kcol+3,i) = 1.d0/((b-bb)*(b-0.7229898685756272d-1)*(b
     #            -0.1863109301186906d0)*(b-0.3341852231986051d0)*(b-0.
     #            5d0)*(b-0.6658147768013949d0)*(b-0.8136890698813094d0)
     #            *(b-0.9277010131424373d0)*b*b*(b-1.d0)*(b-1.d0))
               g(9,kcol+3,i) = 1.d0/((bb-b)*(bb-0.7229898685756272d-1)*
     #            (bb-0.1863109301186906d0)*(bb-0.3341852231986051d0)*
     #            (bb-0.5d0)*(bb-0.6658147768013949d0)*(bb-0.8136890698
     #            813094d0)*(bb-0.9277010131424373d0)*bb*bb*(bb-1.d0)
     #            *(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.7229898685756272d-1)*(x(j)-0.1
     #                  863109301186906d0)*(x(j)-0.3341852231986051d0)
     #                  *(x(j)-0.5d0)*(x(j)-0.6658147768013949d0)*(x(j)
     #                  -0.8136890698813094d0)*(x(j)-0.927701013142437
     #                  3d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.1863109301186906d0*(mesh(i)-mesh(i-1))/hx
               a = -0.7229898685756272d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.3d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)+0.3d2
               h(1,kcol+3,i) = -0.884d3/aa/a
               h(2,kcol+3,i) = 0.884d3/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.7229898685756272d-1)*(aa-0.1
     #            863109301186906d0)*(aa-0.3341852231986051d0)*(aa-0.
     #            5d0)*(aa-0.6658147768013949d0)*(aa-0.813689069881309
     #            4d0)*(aa-0.9277010131424373d0)*(aa-a)*aa*aa*(aa-1.d0)
     #            *(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.7229898685756272d-1)*(a-0.186
     #            3109301186906d0)*(a-0.3341852231986051d0)*(a-0.5d0)*
     #            (a-0.6658147768013949d0)*(a-0.8136890698813094d0)*(a
     #            -0.9277010131424373d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = 0.4624522420172525d5/(0.722989868575627
     #            2d-1-a)/(0.7229898685756272d-1-aa)
               g(4,kcol+3,i) = -0.3688889968687672d5/(0.186310930118690
     #            6d0-a)/(0.1863109301186906d0-aa)
               g(5,kcol+3,i) = 0.3332824691372289d5/(0.3341852231986051
     #            d0-a)/(0.3341852231986051d0-aa)
               g(6,kcol+3,i) = -0.3232914285714286d5/(0.5d0-a)
     #            /(0.5d0-aa)
               g(7,kcol+3,i) = 0.3332824691372289d5/(0.665814776801394
     #            9d0-a)/(0.6658147768013949d0-aa)
               g(8,kcol+3,i) = -0.3688889968687672d5/(0.813689069881309
     #            4d0-a)/(0.8136890698813094d0-aa)
               g(9,kcol+3,i) = 0.4624522420172525d5/(0.927701013142437
     #            3d0-a)/(0.9277010131424373d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.7229898685756272d-1)
     #                  *(x(j)-0.1863109301186906d0)*(x(j)-0.3341852231
     #                  986051d0)*(x(j)-0.5d0)*(x(j)-0.665814776801394
     #                  9d0)*(x(j)-0.8136890698813094d0)
     #                  *(x(j)-0.9277010131424373d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7229898685756272d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1863109301186906d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.3341852231986051d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.6658147768013949d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8136890698813094d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.9277010131424373d0)
               end do
            endif
         end do

c-----------------------------------------------------------------------
      elseif (kcol .eq. 10) then
         do i = 1, nint
            hx = (mesh(i+1)-mesh(i))
            if (i .ne. 1 .and. i .ne. nint) then
               a = -0.5929571219129399d-1*(mesh(i)-mesh(i-1))/hx
               b = 1.d0+0.5929571219129399d-1*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -1.d0/a-0.3666666666666667d2-1.d0/b
               gamma2 = 1.d0/(1.d0-a)+0.3666666666666667d2+1.d0
     #            /(1.d0-b)
               h(1,kcol+3,i) = 0.2799333333333333d4/a/b
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-a)/(1.d0-b)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((a-0.5929571219129399d-1)*(a-0.1539
     #            696908715823d0)*(a-0.2792835119457421d0)*(a-0.4241841
     #            678533667d0)*(a-0.5758158321466333d0)*(a-0.7207164880
     #            542579d0)*(a-0.8460303091284177d0)*(a-0.9407042878087
     #            06d0)*(a-b)*a*a*(a-1.d0)*(a-1.d0))
               g(2,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-a)/(0.5929571219129399d-1-b)
               g(3,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-a)/(0.1539696908715823d0-b)
               g(4,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-a)/(0.2792835119457421d0-b)
               g(5,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-a)/(0.4241841678533667d0-b)
               g(6,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-a)/(0.5758158321466333d0-b)
               g(7,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-a)/(0.7207164880542579d0-b)
               g(8,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-a)/(0.8460303091284177d0-b)
               g(9,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-a)/(0.940704287808706d0-b)
               g(10,kcol+3,i) = 1.d0/((b-a)*(b-0.5929571219129399d-1)*(b
     #            -0.1539696908715823d0)*(b-0.2792835119457421d0)*(b-0.
     #            4241841678533667d0)*(b-0.5758158321466333d0)*(b-0.720
     #            7164880542579d0)*(b-0.8460303091284177d0)*(b-0.940704
     #            287808706d0)*b*b*(b-1.d0)*(b-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-a)*(x(j)-0.5929571219129399d-1)*(x(j)-0.1
     #               539696908715823d0)*(x(j)-0.2792835119457421d0)
     #               *(x(j)-0.4241841678533667d0)*(x(j)-0.5758158321466
     #               333d0)*(x(j)-0.7207164880542579d0)*(x(j)-0.8460303
     #               091284177d0)*(x(j)-0.940704287808706d0)*(x(j)-b)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
               end do

            elseif (i .eq. 1) then
               b = 1.d0+0.5929571219129399d-1*(mesh(i+2)-mesh(i+1))/hx
               bb = 1.d0+0.1539696908715823d0*(mesh(i+2)-mesh(i+1))/hx
               gamma1 = -0.3666666666666667d2-1.d0/b-1.d0/bb
               gamma2 = 0.3666666666666667d2+1.d0/(1.d0-b)
     #            +1.d0/(1.d0-bb)
               h(1,kcol+3,i) = 0.2799333333333333d4/b/bb
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-b)/(1.d0-bb)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-bb)/(0.5929571219129399d-1-b)
               g(2,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-bb)/(0.1539696908715823d0-b)
               g(3,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-bb)/(0.2792835119457421d0-b)
               g(4,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-bb)/(0.4241841678533667d0-b)
               g(5,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-bb)/(0.5758158321466333d0-b)
               g(6,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-bb)/(0.7207164880542579d0-b)
               g(7,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-bb)/(0.8460303091284177d0-b)
               g(8,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-bb)/(0.940704287808706d0-b)
               g(9,kcol+3,i) = 1.d0/((b-bb)*(b-0.5929571219129399d-1)*(b
     #            -0.1539696908715823d0)*(b-0.2792835119457421d0)*(b-0.
     #            4241841678533667d0)*(b-0.5758158321466333d0)*(b-0.720
     #            7164880542579d0)*(b-0.8460303091284177d0)*(b-0.940704
     #            287808706d0)*b*b*(b-1.d0)*(b-1.d0))
               g(10,kcol+3,i) = 1.d0/((bb-b)*(bb-0.5929571219129399d-1)
     #            *(bb-0.1539696908715823d0)*(bb-0.2792835119457421d0)
     #            *(bb-0.4241841678533667d0)*(bb-0.5758158321466333d0)
     #            *(bb-0.7207164880542579d0)*(bb-0.8460303091284177d0)
     #            *(bb-0.940704287808706d0)*bb*bb*(bb-1.d0)*(bb-1.d0))

               do j = 1, kcol+3
                  phi = (x(j)-0.5929571219129399d-1)*(x(j)-0.1539696908
     #               715823d0)*(x(j)-0.2792835119457421d0)*(x(j)-0.4241
     #               841678533667d0)*(x(j)-0.5758158321466333d0)*(x(j)
     #               -0.7207164880542579d0)*(x(j)-0.8460303091284177d0)
     #               *(x(j)-0.940704287808706d0)*(x(j)-b)*(x(j)-bb)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-b)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-bb)
               end do

            else
               aa = -0.1539696908715823d0*(mesh(i)-mesh(i-1))/hx
               a = -0.5929571219129399d-1*(mesh(i)-mesh(i-1))/hx
               gamma1 = -1.d0/aa-1.d0/a-0.3666666666666667d2
               gamma2 = 1.d0/(1.d0-aa)+1.d0/(1.d0-a)
     #            +0.3666666666666667d2
               h(1,kcol+3,i) = 0.2799333333333333d4/aa/a
               h(2,kcol+3,i) = 0.2799333333333333d4/(1.d0-aa)/(1.d0-a)
               hd(1,kcol+3,i) = h(1,kcol+3,i)*hx
               hd(2,kcol+3,i) = h(2,kcol+3,i)*hx
               g(1,kcol+3,i) = 1.d0/((aa-0.5929571219129399d-1)*(aa-0.1
     #            539696908715823d0)*(aa-0.2792835119457421d0)*(aa-0.42
     #            41841678533667d0)*(aa-0.5758158321466333d0)*(aa-0.720
     #            7164880542579d0)*(aa-0.8460303091284177d0)*(aa-0.9407
     #            04287808706d0)*(aa-a)*aa*aa*(aa-1.d0)*(aa-1.d0))
               g(2,kcol+3,i) = 1.d0/((a-0.5929571219129399d-1)*(a-0.1539
     #            696908715823d0)*(a-0.2792835119457421d0)*(a-0.4241841
     #            678533667d0)*(a-0.5758158321466333d0)*(a-0.7207164880
     #            542579d0)*(a-0.8460303091284177d0)*(a-0.9407042878087
     #            06d0)*(a-aa)*a*a*(a-1.d0)*(a-1.d0))
               g(3,kcol+3,i) = -0.1785201442945021d6/(0.592957121912939
     #            9d-1-a)/(0.5929571219129399d-1-aa)
               g(4,kcol+3,i) = 0.1412225097382635d6/(0.153969690871582
     #            3d0-a)/(0.1539696908715823d0-aa)
               g(5,kcol+3,i) = -0.1259240998337037d6/(0.279283511945742
     #            1d0-a)/(0.2792835119457421d0-aa)
               g(6,kcol+3,i) = 0.1197516586185479d6/(0.424184167853366
     #            7d0-a)/(0.4241841678533667d0-aa)
               g(7,kcol+3,i) = -0.1197516586185479d6/(0.575815832146633
     #            3d0-a)/(0.5758158321466333d0-aa)
               g(8,kcol+3,i) = 0.1259240998337037d6/(0.720716488054257
     #            9d0-a)/(0.7207164880542579d0-aa)
               g(9,kcol+3,i) = -0.1412225097382635d6/(0.846030309128417
     #            7d0-a)/(0.8460303091284177d0-aa)
               g(10,kcol+3,i) = 0.1785201442945021d6/(0.94070428780870
     #            6d0-a)/(0.940704287808706d0-aa)

               do j = 1, kcol+3
                  phi = (x(j)-aa)*(x(j)-a)*(x(j)-0.5929571219129399d-1)
     #               *(x(j)-0.1539696908715823d0)*(x(j)-0.2792835119457
     #               421d0)*(x(j)-0.4241841678533667d0)*(x(j)-0.5758158
     #               321466333d0)*(x(j)-0.7207164880542579d0)*(x(j)
     #               -0.8460303091284177d0)*(x(j)-0.940704287808706d0)
                  h(1,j,i) = h(1,kcol+3,i)*(1.d0-x(j)*gamma1)
     #                  *eta(1,j)*phi
                  h(2,j,i) = h(2,kcol+3,i)*(1.d0+gamma2-x(j)
     #                  *gamma2)*eta(2,j)*phi
                  hd(1,j,i) = hd(1,kcol+3,i)*x(j)*eta(1,j)*phi
                  hd(2,j,i) = hd(2,kcol+3,i)*(x(j)-1.d0)*eta(2,j)*phi
                  g(1,j,i) = g(1,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-aa)
                  g(2,j,i) = g(2,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-a)
                  g(3,j,i) = g(3,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5929571219129399d-1)
                  g(4,j,i) = g(4,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.1539696908715823d0)
                  g(5,j,i) = g(5,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.2792835119457421d0)
                  g(6,j,i) = g(6,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.4241841678533667d0)
                  g(7,j,i) = g(7,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.5758158321466333d0)
                  g(8,j,i) = g(8,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.7207164880542579d0)
                  g(9,j,i) = g(9,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.8460303091284177d0)
                  g(10,j,i) = g(10,kcol+3,i)*eta(3,j)*phi
     #                  /(x(j)-0.940704287808706d0)
               end do
            endif
         end do
      endif

      return
      end

      subroutine fdbndx(t, u, ux, dbdu, dbdux, dbdt, npde, bndx, work)
c-----------------------------------------------------------------------
c        This subroutine produces central difference approximations for
c        first partial derivatives of the left and right boundary
c        conditions used in BACOLI models, b(t,u,ux) = 0
c        This subroutine is called by INIY and INIYP in place of
c        user subroutines DIFBXA and DIFBXB in order to simplify
c        the user interface when finite difference approximated Jacobi
c        are used (ie, when mflag(9) = 0).
c-----------------------------------------------------------------------
         implicit none
c        input
         integer npde
         double precision t, u(npde), ux(npde)
         external bndx
c        work storage
         double precision work(npde,2)
c        output
         double precision dbdu(npde,npde), dbdux(npde,npde), dbdt(npde)
c        locals
         external d1mach
         double precision del, d1, d2, delinv, oldval
         double precision squr, uround, d1mach
         integer i, j
c-----------------------------------------------------------------------

c        we will apply increments on the order of size sqrt uround
         uround = d1mach(4)
         squr = sqrt(uround)

c        first, dbdu
         do j = 1, npde
            oldval = u(j)
            del = squr * max(abs(u(j)), uround)
c           d1 and d2 are what is *actually* added or subtracted
            d1 = (u(j) + del) - u(j)
            d2 = (u(j) - del) - u(j)
            u(j) = u(j) + d1
            call bndx(t, u, ux, work(1,1), npde)
            u(j) = oldval + d2
            call bndx(t, u, ux, work(1,2), npde)
            u(j) = oldval
            delinv = 1d0/(d1-d2)
            do i = 1, npde
               dbdu(i,j) = (work(i,1) - work(i,2)) * delinv
            end do
         end do

c        second, dbdux
         do j = 1, npde
            oldval = ux(j)
            del = squr * max(abs(ux(j)), uround)
c           d1 and d2 are what is *actually* added or subtracted
            d1 = (ux(j) + del) - ux(j)
            d2 = (ux(j) - del) - ux(j)
            ux(j) = ux(j) + d1
            call bndx(t, u, ux, work(1,1), npde)
            ux(j) = oldval + d2
            call bndx(t, u, ux, work(1,2), npde)
            ux(j) = oldval
            delinv = 1d0/(d1-d2)
            do i = 1, npde
               dbdux(i,j) = (work(i,1) - work(i,2)) * delinv
            end do
         end do

c        last, dbdt
         del = squr * max(abs(t), uround)
c        d1 and d2 are what is *actually* added or subtracted
         d1 = (t + del) - t
         d2 = (t - del) - t
         call bndx(t+d1, u, ux, work(1,1), npde)
         call bndx(t+d2, u, ux, work(1,2), npde)
         delinv = 1d0/(d1-d2)
         do i = 1, npde
            dbdt(i) = (work(i,1) - work(i,2)) * delinv
         end do

      end
