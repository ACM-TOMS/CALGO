***********************************************************************
*                                                                     *
*                     Routines for SLEDGE                             *
*     (Sturm-Liouville Estimates Determined by Global Errors)         *
*                                                                     *
***********************************************************************
C
C     Release 2.2     12/02/94
C
C     Steven Pruess,  Colorado School of Mines
C                     spruess@mines.colorado.edu
C     Charles Fulton, Florida Institute of Technology
C                     fulton@zach.fit.edu
C
***********************************************************************
*      These routines estimate eigenvalues, eigenfunctions and/or     *
*      spectral density functions for Sturm-Liouville problems.       *
*      The differential equation has the form:                        *
*                                                                     *
*           -(p(x)u')' + q(x)u  =  EV*r(x)u       for x in [A,B]      *
*                                                                     *
*      with boundary conditions (at regular points)                   *
*                                                                     *
*           A1*u - A2*(pu')  =  EV*(A1'*u - A2'*(pu'))    at A        *
*           B1*u + B2*(pu')  =  0                         at B .      *
*                                                                     *
*      The functions p(x) and r(x) are assumed to be positive in      *
*      the open interval (A,B).                                       *
***********************************************************************
C
C      Possible outputs are:
C         a set of eigenvalues;
C         a set of eigenvalues and tables of values for their eigen-
C            functions;
C         a table of the spectral density function (for cases with
C            continuous spectrum).
C         a classification of the problem (regular or singular; if
C            singular then limit point or limit circle, oscillatory
C            or nonoscillatory).
C
C      The code can find eigenvalues and eigenfunctions for problems
C      in spectral category 1 (both endpoints NONOSC), spectral
C      category 2 (one endpoint NONOSC and the other O-NO), and
C      those discrete eigenvalues below the essential spectrum in
C      spectral category 10 (both endpoints O-NO).  Here OSC at an
C      endpoint means the Sturm-Liouville equation is oscillatory for
C      all real values of EV at that endpoint, NONOSC at an endpoint
C      means the equation is nonoscillatory for all real values of EV
C      at that endpoint, and O-NO means there is a `cutoff' value EV'
C      such that the equation is nonoscillatory for real values of
C      EV < EV' and oscillatory for real values of EV > EV'.  For
C      problems in other spectral categories an error return will
C      be generated.  The manner in which SLEDGE classifies singular
C      endpoints of Sturm-Liouville problems as LP/LC (Limit Point/
C      Limit Circle), OSC/NONOSC/O-NO, and uses this information to
C      determine the spectral category is explained in detail in
C      reference [2].
C
C      There is one subroutine called SLEDGE of direct interest to the
C      user; additionally, a secondary routine INTERV is available which
C      determines the indices of eigenvalues located in a specified
C      subinterval of the real line.
C
C      The names of other routines in this package are AITKEN, ASYMEV,
C      ASYMR, BRCKET, CLASS, CLSEND, DENSEF, DSCRIP, EXTRAP, GETEF,
C      GETRN, MESH, POWER, PQRINT, REGULR, SHOOT, START, STEP, and
C      ZZERO.
C
C      There are 4 blocks of labeled COMMON with the names SLREAL,
C      SLINT, SLLOG, and SLCLSS.
C
C      This is the double precision version of the code; all floating
C      point variables should be declared DOUBLE PRECISION in the
C      calling program.  In these subprograms all such local
C      variables and constants have been explicitly declared; also,
C      FORTRAN77 generic intrinsic functions have been used, so
C      conversion to single precision should be straightforward, if
C      desired.
C
C      ACKNOWLEDGMENT:  This work was partially supported by the
C      National Science Foundation under grants DMS-8813113 and DMS-
C      8905202 to Florida Institute of Technology and DMS-8800839 and
C      DMS-8905232 to the Colorado School of Mines.
C
C      References
C
C      The following papers are available from the authors on request:
C
C      [1]. Pruess & Fulton, Mathematical software for Sturm-Liouville
C           problems, ACM Trans. on Math. Software, 19 (1993), 360-376.
C
C      [2]. Fulton, Pruess & Xie, The automatic classification of Sturm-
C           Liouville problems, submitted, 1992.
C
C      [3]. Pruess, Fulton & Xie, An asymptotic numerical method for a
C           class of singular Sturm-Liouville problems, to appear in
C           SIAM J. Numer. Anal.
C
C      [4]  Fulton and Pruess, Eigenvalue and eigenfunction asymptotics
C           for regular Sturm-Liouville problems, Jour. Math. Anal. and
C           Appls., 188 (1994), 297-340.
C
C      [5]  Fulton and Pruess, Numerical Approximation of singular
C           spectral functions arising from the Fourier-Jacobi problem
C           on a half line with continuous spectra, Sixth International
C           Workshop in Analysis and its Applications, June, 1992.

C      [6]. Pruess, Fulton & Xie, Performance of the Sturm-Liouville
C           software package SLEDGE, Colo. School of Mines, Dept. of
C           Math. and Comp. Sci., MCS-91-19, 1991. Revision 12/92.
C-----------------------------------------------------------------------
C      Brief overview of algorithms:
C
C         The code constructs (or takes from input) an initial mesh,
C      called the level 0 mesh.  Subsequent meshes (for level 1,2,...)
C      are unions of the previous level's mesh with its midpoints.  A
C      sequence of estimates for desired eigenvalues and eigenfunctions
C      is constructed, one set for each level.  These estimates (the
C      eigenvalue is called EvHat) are exact solutions (up to the
C      requested tolerance) of a Sturm-Liouville problem which is an
C      approximation to the original one; this approximation results
C      from replacing the given coefficient functions with step function
C      approximations relative to the current level's mesh.  The eigen-
C      functions of the resulting ODE's are piecewise trigonometric
C      (circular or hyperbolic) functions.
C         If estimates for the spectral density function are reqested,
C      these are computed as limits of a sequence of spectral density
C      functions of approximating regular problems.  For these regular
C      problems the spectral density function is a step function, and
C      is computed directly from the definition making use of computed
C      eigenvalues and the norm reciprocals of the corresponding eigen-
C      functions.  If verbose output is rquested by the user, there
C      will be displayed iterations (corresponding to the sequence of
C      approximating regular intervals which the code automatically
C      selects) and within each iteration there will be levels
C      (corresponding to increasingly finer meshes as described above).
C      A step spectral density function will be printed at each level
C      of each iteration.  The spectral density function displayed at
C      the end of each iteration is the result of an h-squared extra-
C      polation over the regular step functions generated at each level
C      of this iteration.  The condition for stopping at a given iter-
C      ation is a straightforward comparison of the spectral function
C      data for the current iteration with the previous iteration. There
C      is no extrapolation over the sequence of regular approximating
C      intervals as no extrapolation theory for the approximation of
C      the singular spectral function by regular step spectral functions
C      is known.  (To achieve closer approximation of the regular step
C      spectral functions to the singular spectral function, it is
C      actually the piecewise linear function obtained by joining the
C      midpoints of successive steps by a straight line which is used as
C      the `regular' spectral function for the purpose of generating the
C      actual data used for the h-squared extrapolation.)
C         The classification is determined by applying standard theory
C      to an approximating problem, each of whose coefficient functions,
C      in a small neighborhood of each endpoint, consists of the leading
C      term in a power-like asymptotic development.  For this reason
C      there are many problems, particularly those with oscillatory
C      coefficient functions, for which the code's output for the
C      classification information is labelled `uncertain'.  For further
C      information on the theory used by the code to generate endpoint
C      classifications and spectral category information see [2] above.
C-----------------------------------------------------------------------
C  Usage (simple explanation) -
C      The subroutine SLEDGE is called in the following manner:
C
C      SUBROUTINE SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,
C                        PDEF,T,RHO,IFLAG,STORE)
C
C      If k eigenvalues (no eigenvectors) are sought then set
C      (a) the logical 5-vector JOB to (True, False, False, False, True);
C      (b) the real 8-vector CONS to the values of A1,A1',A2,A2',B1,B2,
C          A,B for the boundary condition information.  It does not
C          matter what values are used for infinite endpoints, nor for
C          the boundary constants at a singular endpoint; the code
C          automatically selects the Friedrichs' boundary condition at
C          NONOSC singular endpoints, overriding user input for the
C          boundary condition constants; for infinite endpoints the code
C          also automatically selects these constants.
C      (c) the logical 2-vector ENDFIN to (True, True) if both endpoints
C          are finite, (True, False) if A is finite but B infinite, etc.;
C      (d) the integer vector INVEC should have
C              INVEC(1) = 0 (no internal printing)
C              INVEC(2) = 0
C              INVEC(3) = k, the number of eigenvalues sought
C              INVEC(3+i) = index of ith eigenvalue sought,i = 1,...,k;
C      (e) the real 6-vector TOL should have
C              TOL(1) = absolute error tolerance desired,
C              TOL(2) = relative error tolerance desired,
C          the remaining 4 entries of TOL are ignored;
C      (f) the output estimate for the ith eigenvalue is returned
C          in EV(i), i = 1,...,k;
C      (g) the output integer k-vector IFLAG(*) should have all entries
C          zero; nonzero values indicate warnings or error returns and
C          are explained in the detailed usage section below;
C      (h) the auxiliary vector STORE(*) should be dimensioned at least
C          155 in the calling program;
C      (i) the logical 4 by 2 vector TYPE, the real vectors XEF(1),
C          EF(1), PDEF(1), T(1), and RHO(1) can be ignored except that
C          they need to be declared in the calling program.  The integer
C          scalar NUMX can also be ignored.
C
C      If k eigenfunctions are also desired, then follow the above
C      pattern except make JOB(1) False and JOB(2) True.  The values of
C      TOL(3) and TOL(4) control the absolute and relative errors in
C      each u(x); TOL(5) and TOL(6) control the absolute and relative
C      errors in each (pu')(x).  It is usually appropriate to set TOL(5)
C      = TOL(3) = TOL(1) and TOL(6) = TOL(4) = TOL(2), but the user has
C      the option of entering all six tolerance parameters as desired.
C      The output eigenfunction information is returned in the three
C      real vectors X(*) for the independent variable x , EF(*) for u(x),
C      PDEF(*) for (pu')(x).  The code automatically chooses the x
C      values; the number of values is returned in NUMX.  If you prefer
C      another choice of output points, see the detailed explanation
C      below on usage of the code.  The values for the first requested
C      u(x) are returned in the first NUMX locations of EF(*), those for
C      the second are in the next NUMX locations, etc.  PDEF(*) is part-
C      itioned similarly; X(*) must be dimensioned at least 31 in the
C      calling program while EF(*) and PDEF(*) must be dimensioned at
C      least 31*k.  The auxiliary vector STORE(*) should be dimensioned
C      at least 420.
C
C      For other possibilities, see the detailed description which
C      follows.
C-----------------------------------------------------------------------
C  Usage (detailed explanation) -
C
C  SUBROUTINE SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,PDEF,
C                    T,RHO,IFLAG,STORE)
C
C  Input parameters;
C     JOB(*)    = logical 5-vector,
C                 JOB(1) = .True. iff a set of eigenvalues are to be
C                                 computed but not their eigenfunctions.
C                 JOB(2) = .True. iff a set of eigenvalue and eigenfunc-
C                                 tion pairs are to be calculated.
C                 JOB(3) = .True. iff the spectral function is to be
C                                 computed over some subinterval of the
C                                 essential spectrum.
C                 JOB(4) = .True. iff the normal call to the routines for
C                                 classification (regular/singular, etc.)
C                                 is OVERRIDDEN.  If JOB(4) is True then
C                                 TYPE(*,*) discussed below must be
C                                 INPUT correctly!  Most users will not
C                                 want to override the classification
C                                 routines, but it would, of course, be
C                                 appropriate for users experimenting with
C                                 problems for which the coefficient
C                                 functions do not have power-like
C                                 behavior near the singular endpoints.
C                                 Note: the code may perform poorly if
C                                 the classification information is
C                                 incorrect; since the cost is usually
C                                 negligible, it is strongly recommended
C                                 that JOB(4) be False.  The classifica-
C                                 tion is deemed sufficiently important
C                                 for spectral density function calcul-
C                                 ations that JOB(4) is ignored when
C                                 the input JOB(3) is True.
C                 JOB(5) = .True. iff mesh distribution is to be chosen
C                                 by SLEDGE.  If JOB(5) is True and NUMX
C                                 is zero, then the number of mesh
C                                 points is also chosen by SLEDGE; if
C                                 NUMX > 0 then NUMX mesh points will be
C                                 used.  If JOB(5) is False, then the
C                                 number (NUMX) and distribution
C                                 (XEF(*)) must be input by the user.
C                                 If JOB(3) is True and JOB(5) False
C                                 then the user must set BOTH the number
C                                 NUMX and distribution.  In this case,
C                                 NO global error estimates are made.
C     CONS(*)    = real vector of length 8, values are the boundary
C                  condition constants A1, A1', A2, A2', B1, B2, A, B.
C                  In the case of a NONOSC singular endpoint, the class-
C                  ification routine uses the Friedrichs' boundary
C                  condition constants.  The code cannot automatically
C                  choose a non-Friedrichs' boundary condition; however,
C                  interval truncation in the user's calling program can
C                  be used, together with many calls to SLEDGE, to
C                  compute singular eigenvalues associated with a non-
C                  Friedrichs' boundary condition at a NONOSC endpoint
C                  (see remark 12 below).
C     ENDFIN(*) = logical 2-vector, values are
C        ENDFIN(1) = .True. iff endpoint A is finite.
C        ENDFIN(2) = .True. iff endpoint B is finite.
C     INVEC(*)  = integer vector of length 3+(number of eigenvalues
C                 desired).  This vector contains a variety of input
C                 information.
C        INVEC(1) controls the amount of internal printing: values are
C                 from 0 (no printing) to 5 (much printing).
C                 For INVEC(1) > 0 much of the output will be to a file
C                 attached to unit #21 which should be named in the
C                 user's calling program via an OPEN statement.
C                 Output for the various cases is, when INVEC(1) =
C                    0   no printing.
C                 When JOB(1) or JOB(2) is True
C                    1   initial mesh (the first 51 or fewer points),
C                        eigenvalue estimate at each level,
C                    4   the above,
C                        at each level
C                          matching point for eigenfunction shooting,
C                          X(*), EF(*), PDE(*) values,
C                    5   all the above,
C                        at each level
C                          brackets for the eigenvalue search,
C                          intermediate shooting info for the eigen-
C                          function and eigenfunction norm.
C                 When JOB(3) is True
C                    1   the actual (a,b) used at each iteration,
C                        the total number of eigenvalues computed,
C                    2   the above,
C                        switchover points to the asymptotic formulas,
C                        some intermediate Rho(t) approximations,
C                    3   all the above,
C                        initial meshes for each iteration,
C                        index of the largest EV which may be computed,
C                        various Ev and RsubN values,
C                    4   all of the above,
C                        RhoHat values at each level,
C                    5   all of the above,
C                        all Ev and RsubN values below switchover point.
C                 When JOB(4) is False
C                    2   output a description of the spectrum,
C                    3   the above plus the constants for the
C                        Friedrichs' boundary condition(s),
C                    5   all the above plus intermediate details of
C                        classification calculation.
C                 Some of the output may go to the default output device
C                 (screen or printer), but all information requested is
C                 also directed to the file attached to unit #21.
C        INVEC(2) gives the number (positive) of output values desired
C                 for the array RHO(*) (not referenced if JOB(3) is
C                 False).
C        INVEC(3) is the total number of eigenvalues to be output in
C                 EV(*).
C        INVEC(J) for J = 4, 5, ..., 3+INVEC(3) contains the indices for
C                 the eigenvalues sought. If JOB(1) and JOB(2) are
C                 False, this part of INVEC(*) is not referenced.
C     TOL(*)    = real vector of from 2 to 6 tolerances.
C                 If JOB(1) or JOB(2) is True then
C                   TOL(1) is the absolute error tolerance for e-values,
C                   TOL(2) is the relative error tolerance for e-values,
C                   TOL(3) is the abs. error tolerance for e-functions,
C                   TOL(4) is the rel. error tolerance for e-functions,
C                   TOL(5) is the abs. error tolerance for eigenfunction
C                          derivatives,
C                   TOL(6) is the rel. error tolerance for eigenfunction
C                          derivatives.
C                   Eigenfunction tolerances need not be set if JOB(2)
C                   is False.
C                 If JOB(3) is True then
C                   TOL(1) is the absolute error tolerance,
C                   TOL(2) is the relative error tolerance;
C                   the output RHO values are NOT required to satisfy
C                   these tolerances when JOB(5) is False.
C                 All absolute error tolerances must be positive; all
C                 relative error tolerances must be at least 100 times
C                 the unit roundoff.
C     NUMX      = integer whose value is
C                    the number of output points where each eigen-
C                    function is to be evaluated (the number of entries
C                    in XEF(*)) when JOB(2) is True,
C                 or
C                    the number of points in the initial mesh used when
C                    JOB(5) is False and NUMX>0.
C                 If JOB(5) is False, the points in XEF(*) should be
C                 chosen to have a reasonable distribution.  Since the
C                 endpoints A and B must be part of any mesh, NUMX
C                 cannot be 1 in this case.  If JOB(5) is FALSE and
C                 JOB(3) is True, then NUMX must be positive.
C     XEF(*)    = real vector of points where
C                    eigenfunction estimates are desired (JOB(2) True)
C                 or
C                    where user's initial mesh is entered (JOB(5) False
C                    and NUMX>0).
C                 The values must satisfy
C                      A = XEF(1) < XEF(2) < ... < XEF(NUMX) = B .
C                 When JOB(2) is True the initial mesh corresponds to
C                 the set of points where eigenfunction output is
C                 desired. If JOB(2) is False and NUMX = 0, then this
C                 vector is not referenced.  When A and/or B are
C                 infinite (as indicated through ENDFIN(*)), the
C                 entries XEF(1) and/or XEF(NUMX) are ignored; however,
C                 it is required that XEF(2) be negative when ENDFIN(1)
C                 is False, and XEF(NUMX-1) be positive when ENDFIN(2)
C                 is False (otherwise, IFLAG = -39 will result).
C     T(*)      = real vector of INVEC(2) values where the spectral
C                 function RHO(*) is desired (the existence and location
C                 of continuous spectrum can be found by first calling
C                 SLEDGE with JOB(J) False, J=1,...,4 and INVEC(1) = 1).
C                 Vector T(*) is not referenced if JOB(3) is False.  Its
C                 entries must be in increasing order.
C
C   Output parameters:
C     TYPE(*,*) = 4 by 2 logical array; column 1 carries information
C                 about endpoint A while column 2 refers to B.
C                 TYPE(1,*) = True  iff the endpoint is regular,
C                 TYPE(2,*) = True  iff it is limit circle,
C                 TYPE(3,*) = True  iff it is nonoscillatory for all EV,
C                 TYPE(4,*) = True  iff it is oscillatory for all EV,
C                 Important note: all of these must be correctly INPUT
C                 if JOB(4) is True!
C     EV(*)     = real vector containing the computed approximations to
C                 the eigenvalues whose indices are specified in
C                 INVEC(*); if JOB(1) and JOB(2) are False, then the
C                 output has no meaning.
C     NUMX      = the number of output points for eigenfunctions when
C                 input NUMX = 0, and JOB(2) or JOB(5) is True.
C     XEF(*)    = input values (if any) are changed only if JOB(2) and
C                 JOB(5) are True; in this case, the output values
C                 are chosen by the code.  If JOB(2) is False then this
C                 vector is not referenced; if JOB(2) is True and NUMX>0
C                 on input then XEF(*) should be dimensioned at least
C                 NUMX+16 in the calling program.  If JOB(2) is True and
C                 NUMX=0 on input (so that the code chooses NUMX), then
C                 dimension XEF(*) at least 31 in the calling program.
C     EF(*)     = real vector of eigenfunction values: EF((k-1)*NUMX+i)
C                 is the estimate of u(XEF(i)) corresponding to the
C                 eigenvalue in EV(k).  If JOB(2) is False then this
C                 vector is not referenced.  Otherwise, if JOB(2) is
C                 True and NUMX>0 on input then EF(*) should be
C                 dimensioned at least NUMX*INVEC(3) in the calling
C                 program.  If JOB(2) is True and NUMX=0 on input (so
C                 that the code chooses NUMX), then dimension XEF(*)
C                 at least 31*INVEC(3) in the calling program.
C     PDEF(*)   = real vector of eigenfunction derivative values:
C                 PDEF((k-1)*NUMX+i) is the estimate of (pu')(XEF(i))
C                 corresponding to the eigenvalue in EV(k).  If JOB(2)
C                 is False then this vector is not referenced; otherwise,
C                 it must be dimensioned as is EF(*).
C     RHO(*)    = real vector of values for the spectral density
C                 function rho(t), RHO(I) = rho(T(I)).  RHO(*) must be
C                 dimensioned at least INVEC(2); this vector is not
C                 referenced if JOB(3) is False.
C     IFLAG(*)   = integer vector carrying information about the output.
C       Declared length must be at least max(1,INVEC(3)).  For the Kth
C       requested eigenvalue (when JOB(1) or JOB(2) is true; otherwise,
C       only IFLAG(1) is used):
C       IFLAG(K) =  0, normal return, output should be reliable.
C                <  0, fatal error, calculations ceased: if
C                = -1, too many levels needed for the eigenvalue
C                      calculation; problem seems too difficult for
C                      this algorithm at this tolerance. Are the
C                      coefficient functions nonsmooth?
C                = -2, too many levels needed for the eigenfunction
C                      calculation; problem seems too difficult for
C                      this algorithm at this tolerance.  Are the
C                      eigenfunctions ill-conditioned?
C                = -3, too many levels needed for the spectral density
C                      calculation; problem seems too difficult for
C                      this algorithm at this tolerance.
C                = -4, the user has requested the spectral density
C                      function for a problem which has no continuous
C                      spectrum.
C                = -5, the user has requested the spectral density
C                      function for a problem with both endpoints
C                      generating essential spectrum, i.e., both
C                      endpoints being either OSC or O-NO.  The spectral
C                      density function calculation has not been
C                      implemented for such cases.  For spectral
C                      category 10 (both endpoints O-NO) the spectral
C                      multiplicity is generally two, proper normal-
C                      izations for the solutions against which the
C                      spectral functions will be normalized will
C                      depend on how the user wants to express the
C                      eigenfunction expansion.  Users having problems
C                      in spectral category 10 are encouraged to supply
C                      them to the authors, and if possible, recommend
C                      normalizations of the two solutions to be used in
C                      writing the associated eigenfunction expansion.
C                = -6, the user has requested the spectral density
C                      function for a problem in spectral category 2 for
C                      which a proper normalization of solution at the
C                      NONOSC endpoint is not known; for example,
C                      problems with an irregular singular point or
C                      infinite endpoint at one end and continuous
C                      spectrum generated at the other.  Users with
C                      problems of this type are encouraged to supply
C                      them to the authors, and if possible, recommend a
C                      normalization of solution at the NONOSC endpoint
C                      which they would like to see implemented. As a
C                      rule it is best to pick a normalization which
C                      ensures that the solution is uniquely fixed and
C                      entire in the eigenvalue parameter EV for all x
C                      in the Sturm-Liouville interval; for further
C                      mathematical information on NONOSC endpoints we
C                      refer to paper [2] above.
C                = -7, problems encountered in obtaining a bracket.
C                = -8, too small a step used in the integration;
C                      TOL(*) values may be too small for this problem.
C                = -9, too small a step used in a spectral density
C                      function calculation for which the continuous
C                      spectrum is generated by a finite endpoint.  Try
C                      transforming to Liouville (or some other) form.
C                = -10, an argument to the circular trig functions is
C                       too large.  Try rerunning with a finer initial
C                       mesh, or, on singular problems, use interval
C                       truncation (see remark (12)).
C                = -15, p(x) and r(x) not positive in (A,B).
C                = -20, eigenvalues/functions were requested for a
C                       problem with an OSC singular endpoint.
C                       Interval truncation (see remark (12)) must be
C                       used on such problems.
C                = -3?, illegal input, viz.
C                  -30,  NUMX = 1 when JOB(5) is True,
C                        or NUMX = 0 when JOB(3) is True and JOB(5) is
C                        False,
C                  -31,  B1 = B2 = 0 (at a regular endpoint),
C                  -32,  A1'*A2-A1*A2' .le. 0 when A1' or A2' nonzero,
C                  -33,  A1 = A2 = A1'= A2'= 0 (at a regular endpoint),
C                  -34,  A .ge. B (when both are finite),
C                  -35,  TOL(odd) .le. 0 ,
C                  -36,  TOL(even)  <  100*unit roundoff,
C                  -37,  INVEC(k) < 0  for some k>3 when INVEC(3)>0,
C                  -38,  INVEC(2) .le. 0 when JOB(3) is True ,
C                  -39,  XEF(*) entries out of order or not in [A,B].
C                        or XEF(2), XEF(NUMX-1) have the wrong sign in
C                           infinite interval cases,
C                        or T(*) entries are out of order.
C                >  0,  indicates some kind of warning, in this case the
C                       value may contain ANY of the following digits:
C                =  1,  failure in routine BRCKET probably due to a
C                       cluster of eigenvalues which the code cannot
C                       separate.  Calculations have continued as best
C                       as possible, but any eigenfunction results are
C                       suspect.  Try rerunning with tighter input
C                       tolerances to separate the cluster.
C                =  2,  there is uncertainty in the classification for
C                       this problem.  Because of the limitations of the
C                       floating point arithmetic on the computer used,
C                       and the nature of the finite sampling, the
C                       routine is cannot be decisive about the
C                       classification information at the requested
C                       tolerance.
C                =  3,  there may be some eigenvalues imbedded in the
C                       essential spectrum; using IPRINT greater than
C                       zero will result in additional output giving
C                       the location of the approximating eigenvalues
C                       for the step function problem.  These could be
C                       extrapolated to estimate the actual eigenvalue
C                       embedded in the essential spectrum.
C                =  5,  a change of variables was made to avoid poten-
C                       tial slow convergence; however, the global
C                       error estimates may not be as reliable.  Some
C                       experimentation using different tolerances is
C                       recommended.
C                =  6,  there were problems with eigenfunction conver-
C                       gence in a spectral density calculation; the
C                       output Rho(t) may not be accurate.
C
C   Auxiliary storage:
C     STORE(*) = real vector of auxiliary storage, must be dimensioned
C                at least
C             max(155,NUMX+16)     in general;
C               26*(NUMX+16)       for any eigenfunction calculation;
C             2400+13*INVEC(2)     for any spectral density calculation.
C-----------------------------------------------------------------------
C  SUBROUTINE INTERV(FIRST,ALPHA,BETA,CONS,ENDFIN,NFIRST,NTOTAL,
C                    IFLAG,STORE)
C
C    Input parameters:
C     FIRST      = logical; value is True if various internal variables
C                  have not yet been set.  If a prior call has been made
C                  to INTERV with FIRST True, then a little time can
C                  be saved by letting FIRST be False.
C                  IMPORTANT NOTE: setting FIRST = True will clobber any
C                  initial mesh the user has input (when NUMX > 0 or
C                  JOB(5) is False);  also, INTERV will classify the
C                  problem irregardless of what JOB(4) is set to
C                  for SLEDGE.
C      ALPHA     = real value of left end point of search interval.
C      BETA      = real value of right end point of search interval.
C      CONS(* )  = real vector of 8 input constants: A1, A1', A2, A2',
C                  B1, B2, A, B.
C      ENDFIN(*) = logical 2-vector, same meaning as in SLEDGE.
C      STORE(*)  = real vector holding initial mesh.
C
C    Output parameters:
C      NFIRST = index of first eigenvalue > ALPHA.
C      NTOTAL = total number of eigenvalues in the interval.
C      IFLAG  = integer status indicator.
C               IFLAG =   0 , normal return, output should be reliable,
C                     =  11 , there are no eigenvalues in [alpha, beta],
C                     =  12 , low confidence in NFIRST or NTOTAL or both,
C                     =  13 , BETA and/or ALPHA exceed the cutoff for
C                             the continuous spectrum.  If only BETA
C                             is too big then NFIRST may be OK, but
C                             NTOTAL is meaningless.
C                     = -11 , ALPHA .ge. BETA,
C                     = -25 , oscillatory endpoint, output meaningless,
C                     = -3? , illegal CONS(*) values (see above comments
C                             on SLEDGE for an explanation).
C--------------------------------------------------------------------------
C         In addition, a subroutine subprogram must be provided for the
C      coefficient functions p(x), q(x), and r(x); the form of this
C      routine is
C
C          SUBROUTINE COEFF(X,PX,QX,RX)
C          DOUBLE PRECISION X,PX,QX,RX
C               ...
C          PX = ...
C          QX = ...
C          RX = ...
C          RETURN
C          END
C
C      The subroutine name MUST be COEFF, though of course the names of
C      arguments only need follow the usual FORTRAN77 rules. X is the
C      independent variable; PX, QX, and RX are the output values of the
C      respective coefficient functions p(x), q(x), and r(x) at X.
C-----------------------------------------------------------------------
C     This is a simple sample driver for SLEDGE.
CC
CC     Declare all variables:
CC
C      INTEGER IFLAG(1),INVEC(4),NUMX, I,J,K
C      LOGICAL JOB(5),TYPE(4,2),ENDFIN(2)
C      DOUBLE PRECISION CONS(8),TOL(6),EV(1),T(3),RHO(3),STORE(2450),
C     &                 XEF(5),EF(5),PDEF(5)
CC
CC     Load the boundary condition information into CONS(*).
CC     This example has a Neumann condition at A = 1, and a
CC     singular point at B = +infinity.
CC
C      DATA CONS/0.0, 0.0, 1.0, 0.0,   0.0, 0.0,   1.0, 0.0/
C      DATA ENDFIN/.TRUE., .FALSE./
CC
CC     The eigenfunctions will be estimated at 5 points.
CC
C      DATA NUMX,XEF/5, 1.0, 1.5D0, 2.0, 4.0, 100.0/
CC
CC     Initialize the vector INVEC(*):
CC        little printing,
CC        3 output points for the density function Rho(t),
CC        estimates for the first (index 0) eigenvalue/function.
CC
C      DATA INVEC/1, 3, 1, 0/
CC
CC     Set the JOB(*) vector:
CC        estimate both eigenvalues and eigenvectors,
CC        estimate the spectral density function,
CC        classify,
CC        force the initial mesh to be the output points.
CC
C      DATA JOB/.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE./
CC
CC     Set the tolerances:
CC
C      DATA TOL/1.D-5,1.D-4,  1.D-5,1.D-4,  1.D-5,1.D-4/
CC
CC     Initialize the 3 output points for the density function.
CC
C      DATA T/0.0, 0.5, 2.0/
CC
CC     Open file for output.
CC
C      OPEN(21,FILE = 'sample.out')
C      CALL SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,PDEF,
C     &            T,RHO,IFLAG,STORE)
CC
CC     Print results:
CC
C      DO 30 I = 1,INVEC(3)
C         WRITE (*,10) INVEC(3+I),EV(I),IFLAG(I)
C         WRITE (21,10) INVEC(3+I),EV(I),IFLAG(I)
C   10    FORMAT(' Nev =',I6,';   Ev =',D25.15,';     Flag = ',I3)
C         IF (IFLAG(I) .GT. -10) THEN
C            WRITE (*,15)
C            WRITE (21,15)
C   15       FORMAT(13X,'x',23X,'u(x)',18X,'(pu`)(x)')
C            K = NUMX*(I-1)
C            DO 25 J = 1,NUMX
C               WRITE (21,20) XEF(J),EF(J+K),PDEF(J+K)
C   20          FORMAT(3D25.15)
C   25       CONTINUE
C         ENDIF
C   30    CONTINUE
C      WRITE (*,35)
C      WRITE (21,35)
C   35 FORMAT(/,8X,'t',21X,'Rho(t)')
C      DO 45 I = 1,INVEC(2)
C         WRITE (*,40) T(I),RHO(I)
C         WRITE (21,40) T(I),RHO(I)
C   40    FORMAT(F11.3,D32.15)
C   45 CONTINUE
C      CLOSE(21)
C      STOP
C      END
CC
C      SUBROUTINE COEFF(X,PX,QX,RX)
CC
CC     Define the coefficient functions; here a Yukawa potential.
CC
C      DOUBLE PRECISION X,PX,QX,RX, T
CC
CC     Be careful with potential over/underflows; here we assume the
CC     IEEE double precision exponent range.
CC
C      IF (X .LT. 650.0) THEN
C         T = EXP(-X)
C      ELSE
C         T = 0.0
C      ENDIF
C      PX = 1.0
C      QX = -T/X
C      RX = PX
C      RETURN
C      END
CC
CC     End of sample driver for SLEDGE.
C-----------------------------------------------------------------------
C      General remarks:
C      (1) Two machine dependent constants must be set in a DATA
C          statement in routine START (in part 4 of the package):
C          URN   -  an estimate of the unit roundoff; infinite output
C                   values are assigned the value 1/URN.
C          UFLOW -  a number somewhat smaller than -ln(underflow level).
C                   Values of certain variables z for which
C                             ln(abs(z)) < -under
C                   will be set to zero.
C      (2) A value of IFLAG = -1, -2, or -3 may be the result of a
C          lack of smoothness in the coefficient functions.  In such
C          cases a user input mesh may perform better (see (4) below).
C      (3) The heuristics for generating the initial mesh distribution
C          work reasonably well over a wide range of examples, but
C          occasionally they are far from optimal.  The code's choice
C          can be over-ridden by setting JOB(5) False, setting NUMX
C          appropriately and supplying a mesh in XEF(*).
C      (4) If any of the coefficient functions p,q, or r (or their first
C          few derivatives) have finite jump discontinuities at points
C          in the interior of (A,B), then it is advantageous to have
C          these points in SLEDGE's mesh.  Currently, this can only be
C          accomplished by setting JOB(5) False and supplying an
C          appropriate mesh using NUMX and XEF(*).
C      (5) In general, eigenvalue convergence is observed to be more
C          rapid than eigenfunction convergence; hence, it is
C          recommended that JOB(2) be False unless eigenfunction
C          information really is necessary.
C      (6) When eigenfunction output is sought, unless some knowledge
C          of the eigenfunction is known in advance, it is recommended
C          that JOB(5) be True so that the code will attempt to choose
C          a reasonable distribution for the initial mesh points.
C      (7) Computing the spectral density function for problems having
C          continuous spectrum can be very expensive; it is recommended
C          that initially, relatively crude tolerances (0.001 or so) be
C          used to get some idea of the effort required.
C      (8) It is recommended that every problem be classified (JOB(4)
C          False) by the code before any calculation of spectral
C          quantities occurs.  Only if the user is certain as to what
C          the classification is (and describes it correctly through
C          INVEC and TYPE) should the classification option be bypassed.
C      (9) If the code does the classification of singular problems, it
C          will automatically choose the Friedrichs' boundary condition
C          at NONOSC endpoints.  If another boundary condition is
C          desired, the user must use interval truncation in the
C          calling program (see remark (12)).
C     (10) While all parts of the code should function on machines
C          with a fairly narrow exponent range (such as IEEE single
C          precision), it is better to have a relatively wide exponent
C          range (IEEE double precision).  The classification algorithm,
C          in particular, is far more reliable if done on a machine with
C          a fairly wide exponent range.
C     (11) Care must be taken in writing the subroutine COEFF for the
C          evaluation of p(x), q(x), and r(x) to avoid arithmetic
C          exceptions such as overflow and underflow (or trig function
C          arguments too large).  This can be especially delicate on
C          machines with a small exponent range.
C     (12) In some cases `interval truncation' is recommended.  By this
C          is meant the user should call SLEDGE several times using a
C          sequence of regular endpoints (with appropriate boundary
C          conditions) converging to the singular endpoint.  The eigen-
C          values of the regular problems selected by the user should be
C          arranged so as to converge to those of the desired singular
C          problem. For example, if the user wishes to compute eigen-
C          values associated with a non-Friedrichs' boundary condition
C          for problems in spectral category 1, the user can experiment
C          with choosing a sequence of regular approximating intervals,
C          and vary the boundary conditions appropriately by means of a
C          `boundary condition function' or known solution of the
C          equation for a real value of EV on the sequence of regular
C          intervals until convergence of the regular eigenvalues to
C          the desired singular one is observed.  Similarly, for
C          problems in spectral category 3 or 5 which involve one or two
C          endpoints which are LC and OSC, the (necessarily discrete)
C          spectrum is known to be unbounded below and above. To
C          implement a given LC boundary condition at a singular LC
C          endpoint one may choose a `boundary condition function' or
C          known solution of the equation for a real value of EV and
C          make use of it on a sequence of regular approximating
C          intervals to vary the boundary condition on successive calls
C          to SLEDGE for the sequence of regular intervals until
C          convergence to the desired singular eigenvalue is observed.
C          At present these methods are highly experimental and problem-
C          dependent as good heuristics for the choice of the rate of
C          convergence of the regular intervals to the singular one
C          which work well over a wide class of problems are not known.
C          (The only case in which SLEDGE automatically selects regular
C          approximating subintervals is for spectral density function
C          calculations for problems in spectral category 2; but
C          here the singular endpoint is of LP type, so no singular
C          boundary condition is required to be implemented.)
C     (13) Problems of slow convergence can sometimes be avoided by a
C          judicious change of either dependent or independent variable
C          (or both).
C     (14) If the Liouville normal form potential Q(t) has a minimum
C          far from zero, then the heuristics for generating the initial
C          mesh may well miss it.  In this case, it is advisable to
C          shift the independent variable.
C     (15) The determination of the total number of eigenvalues is the
C          most difficult part of the classification process.  When the
C          theory provides this number, of course, there is no problem;
C          otherwise, it should be viewed with some skepticism.  A more
C          reliable count of the eigenvalues below the cutoff point of
C          the essential spectrum can be gained (at some expense) by
C          trying to compute many eigenvalues near that point.
C-----------------------------------------------------------------------
C     Changes since version 2.1:
C     (1) bug fix in MESH when a or b infinite            09/12/91
C     (2) new mesh heuristics for infinite intervals      09/13/91
C     (3) bug fix in CLSEND: undefined CP(*)              09/15/91
C     (4) bug fix in AA, BB definition of SLEDGE          09/19/91
C     (5) added NUMEV to DENSEF                           09/22/91
C     (6) fixed NADD initialization in MESH               09/24/91
C     (7) fixed uninitialized ASYMEV value                09/27/91
C     (8) change eigenvalue count in CLASS                10/01/91
C     (9) improved NZERO calculation in STEP              10/07/91
C    (10) bug fix on testing input bc                     10/10/91
C    (11) relaxed error tests on Rho(t) in EXTRAP         10/13/91
C    (12) altered printing options in DENSEF,MESH,SLEDGE  10/16/91
C    (13) altered KCLASS values                           10/16/91
C    (14) fixed RLOW, RHIGH bugs in EXTRAP                10/17/91
C    (15) add KCLASS print to CLSEND                      10/19/91
C    (16) finished heuristics for density calculation     10/21/91
C    (17) minor change to SYMM part of GETEF              11/06/91
C    (18) changed printing format in GETEF                11/06/91
C    (19) updated INTERV to conform to earlier changes    11/06/91
C    (20) moved initialization of U, UNDER back to START  11/15/91
C    (21) wholesale bug fixes in INTERV                   11/15/91
C    (22) added print of spectral category to DSCRIP      11/21/91
C    (23) delete FLAG = 0 in DENSEF                       12/02/91
C    (24) tightened up oscillatory test in POWER          04/19/92
C    (25) altered printouts for oscillatory coeff. case   04/19/92
C    (26) increased max NxInit for DENSEF                 04/25/92
C    (27) pass iteration index to DENSEF                  04/25/92
C    (28) minor change to POWER                           04/30/92
C    (29) allow user to input mesh to DENSEF              05/09/92
C    (30) renumbered labels in SLEDGE                     05/09/92
C    (31) minor changes to CLASS, DSCRIP                  05/10/92
C    (32) tightened switchover test in DENSEF             05/17/92
C    (33) altered tests for IRREG in CLSEND               06/24/92
C    (34) avoid integer overflow in POWER                 06/24/92
C    (35) added print option to EXTRAP                    06/25/92
C    (36) NUMX reset by SLEDGE; change comments           07/21/92
C    (37) further tinkering with OSC in POWER, CLSEND     07/22/92
C    (38) bug fix in KCLASS = 10 case in CLSEND           08/27/92
C    (39) added print options to POWER                    11/03/92
C    (40) Aitken used instead of Wynn in POWER            11/11/92
C    (41) redo eigenfunction pointers in REGULR           12/23/92
C    (42) new mesh heuristics for KCLASS = 3              12/31/92
C    (43) redo user interface                             01/28/93
C    (44) don't extrapolate infinity from GETEF           02/13/93
C    (45) more changes to user interface                  03/14/93
C    (46) update format 65 in SLEDGE                      04/12/93
C    (47) bug fix in END calculation in CLASS             01/24/94
C    (48) eliminate undefined FLAG test in DENSEF         12/02/94
C///////////////////////////////////////////////////////////////////////
      module SLEDGEMD
      private
      public:: SLEDGE,INTERV

      contains

      subroutine SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,
     +                  PDEF,T,RHO,IFLAG,STORE)
C
C     This is the interface routine between the user and other routines
C     which carry out most of the actual calculations.
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,TWO,TOLMAX
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,TWO=2.0,TOLMAX=1.D-4)
C     ..
C     .. Scalar Arguments ..
      integer NUMX
C     ..
C     .. Array Arguments ..
      double precision CONS(*),EF(*),EV(*),PDEF(*),RHO(*),STORE(*),T(*),
     +                 TOL(*),XEF(*)
      integer IFLAG(*),INVEC(*)
      logical ENDFIN(*),JOB(*),TYPE(4,*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision AA,ALPHA,BB,DENS,DENSHI,DENSLO,DENSOP,ENDFAC,
     +                 ERROR,FZ,HMIN,RHOTOL,SGN,TOL1,XTOL,ZETA
      integer I,IBASE,IEV,IPRINT,J,JTOL,K,KCL1,KCL2,LASTEV,MAXITS,MAXT,
     +        MU1,MU2,NEV,NEXTRP,NUMEV,NUMT
      logical AAFIN,BBFIN,DOMESH,DONE,EDONE,LBASE,LMESH,OSCILL
C     ..
C     .. Local Arrays ..
      double precision CEV(2),ENDI(5),ZETAI(5)
      logical CSPEC(2),JOBST(3),LPLC(2)
C     ..
C     .. External Subroutines ..
cc      external CLASS,DENSEF,DSCRIP,MESH,REGULR,SHOOT,START
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,LOG10,MAX,MIN,MOD,SIGN
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Data statements ..
      data DENSLO,DENSOP,DENSHI/4.0,6.0,12.0/
      data ENDI/12.0,20.0,85.0,240.0,500.0/
      data ZETAI/2.2,2.0,1.5,1.4,1.3/
C     ..
C
C     Initialize.
C
      AFIN = ENDFIN(1)
      BFIN = ENDFIN(2)
      IPRINT = INVEC(1)
      NUMT = INVEC(2)
      NEV = INVEC(3)
      LNF = .FALSE.
      DOMESH = .TRUE.
      LMESH = .FALSE.
      FLAG = 0
      IFLAG(1) = 0
      IBASE = 1
      LBASE = .FALSE.
      do 5 I = 1,6
         LFLAG(I) = .FALSE.
    5 continue
      TOL1 = MIN(TOL(1)+TOL(2),TOLMAX)
      JOBST(1) = JOB(2)
      if ((NUMX.gt.0) .and. (.not.JOB(5))) then
         JOBST(2) = .TRUE.
      else
         JOBST(2) = .FALSE.
      end if

      JOBST(3) = JOB(3)
      if ((.not.JOB(1)) .and. (.not.JOB(2))) NEV = 0
      call START(JOBST,CONS,TOL,NEV,INVEC(4),NUMX,XEF,NUMT,T,NEXTRP,
     +           STORE)
      if (JOB(4)) then
         ALPHA = A2*A1P - A1*A2P
         if ((A1P.ne.ZERO) .or. (A2P.ne.ZERO)) then
            if (ALPHA.le.ZERO) FLAG = -32
         else
            if ((A1.eq.ZERO) .and. (A2.eq.ZERO)) FLAG = -33
         end if

         if ((B1.eq.ZERO) .and. (B2.eq.ZERO)) FLAG = -31
      end if

      if (FLAG.lt.0) go to 120
      if (JOB(1) .or. JOB(2)) then
         do 10 K = 1,NEV
            EV(K) = ZERO
   10    continue
      end if

      if (JOB(3)) then
         do 15 K = 1,NUMT
            RHO(K) = ZERO
   15    continue
      end if

      if ((.not.JOB(4)) .or. JOB(3)) then
         call CLASS(IPRINT,TOL1,JOBST(2),CSPEC,CEV,LASTEV,LPLC,STORE,
     +              JOB(5),HMIN,DOMESH)
         ALPHA = A2*A1P - A1*A2P
         if ((A1P.ne.ZERO) .or. (A2P.ne.ZERO)) then
            if (ALPHA.le.ZERO) FLAG = -32
         else
            if ((A1.eq.ZERO) .and. (A2.eq.ZERO)) FLAG = -33
         end if

         if ((B1.eq.ZERO) .and. (B2.eq.ZERO)) FLAG = -31
         if (FLAG.lt.0) go to 120
         do 20 K = 1,2
            TYPE(1,K) = REG(K)
            TYPE(2,K) = LC(K)
            TYPE(3,K) = .not. OSC(K)
            TYPE(4,K) = OSC(K)
            if (CSPEC(K)) then
               TYPE(3,K) = .FALSE.
               TYPE(4,K) = .FALSE.
            end if

   20    continue
         if (IPRINT.gt.2) call DSCRIP(LC,LPLC,TYPE,REG,CSPEC,CEV,CUTOFF,
     +                                LASTEV,A1,A1P,A2,A2P,B1,B2)
      else
         LNF = .FALSE.
         KCLASS(1) = 0
         KCLASS(2) = 0
         do 25 K = 1,2
            REG(K) = TYPE(1,K)
            LC(K) = TYPE(2,K)
            OSC(K) = TYPE(4,K)
            CSPEC(K) = .not. (TYPE(3,K) .or. TYPE(4,K))
   25    continue
         if (.not.AFIN) STORE(1) = -99999.0
         if (.not.BFIN) STORE(NXINIT) = 99999.0
      end if
C
C     Use NSGNF to hold the sign of F when EV is large negative.
C
      SGN = A2P*B2
      if (SGN.ne.ZERO) then
         NSGNF = SIGN(ONE,SGN)
      else
         SGN = A1P*B2 + A2P*B1
         if (SGN.ne.ZERO) then
            NSGNF = SIGN(ONE,SGN)
         else
            SGN = A1P*B1 + A2*B2
            if (SGN.ne.ZERO) then
               NSGNF = SIGN(ONE,SGN)
            else
               SGN = A1*B2 + A2*B1
               if (SGN.ne.ZERO) then
                  NSGNF = SIGN(ONE,SGN)
               else
                  NSGNF = SIGN(ONE,A1*B1)
               end if

            end if

         end if

      end if

      OSCILL = ((.not.TYPE(1,1)) .and. TYPE(4,1)) .or.
     +         ((.not.TYPE(1,2)) .and. TYPE(4,2))
      TOL1 = TOL(1) + TOL(2)
      if (JOB(1) .or. JOB(2)) then
C
C        Set up approximating regular problems for eigenvalues.
C
         if (OSCILL) then
            if (IPRINT.ge.1) then
               write (*,FMT=30)
               write (21,FMT=30)

   30          format (' This problem is oscillatory, you must use ',
     +                'interval truncation.')

            end if

            FLAG = -20
            go to 120

         end if

         if (DOMESH) then
C
C           Calculate the initial mesh.
C
            K = NXINIT + 16
            call MESH(JOB(5),-1,STORE,STORE(K),STORE(2*K+1),
     +                STORE(3*K+1),STORE(4*K+1),TOL1,HMIN)
            if (FLAG.lt.0) go to 120
         end if

         if (((KCLASS(1).eq.3).or. (KCLASS(2).eq.3)) .and.
     +       JOB(5)) LMESH = .TRUE.
         if ((.not.LMESH) .and. (IPRINT.ge.1)) then
            write (*,FMT=35) (STORE(I),I=1,NXINIT)
            write (21,FMT=35) (STORE(I),I=1,NXINIT)

   35       format (' Level 0 mesh:',/ (5g15.6))

         end if

         if (JOB(5)) NUMX = NXINIT
C
C        Set MAXLVL, the maximum number of levels (mesh bisections).
C
C        IMPORTANT NOTE: the size of various fixed arrays in this
C        package depends on the value of MAXLVL in this FORTRAN77
C        implementation.  If MAXLVL is increased, then more storage
C        may have to be allocated to these arrays.  In particular,
C        check RATIO(*), R(*,*), and W(*,*) in EXTRAP; EVEXT(*)
C        in REGULR.
C
         MAXLVL = 10
C

         do 45 K = 1,NEV
            EV(K) = ZERO
            IFLAG(K) = 0
            FLAG = 0
            call REGULR(JOB(2),LMESH,TOL,INVEC(3+K),EV(K),IPRINT,NEXTRP,
     +                  XEF,EF(1+NUMX* (K-1)),PDEF(1+NUMX* (K-1)),HMIN,
     +                  STORE)
            if ((CSPEC(1).or.CSPEC(2)) .and. (IPRINT.ge.1) .and.
     +          (.not.JOB(4)) .and. (FLAG.gt.-5)) then
               if ((EV(K).ge.CUTOFF) .or. ((LASTEV.ne.-5).and.
     +             (INVEC(3+K).ge.LASTEV))) then
                  write (*,FMT=40) INVEC(3+K)
                  write (21,FMT=40) INVEC(3+K)

   40             format (' WARNING: Requested eigenvalue ',i6,
     +                   ' may not be below the continuous spectrum.')

               end if

            end if

            if (LFLAG(1)) then
               IFLAG(K) = IFLAG(K) + IBASE
               LFLAG(1) = .FALSE.
               LBASE = .TRUE.
            end if

            if (FLAG.lt.0) IFLAG(K) = FLAG
   45    continue
         if (LBASE) then
            IBASE = 10*IBASE
            LBASE = .FALSE.
         end if

      end if

      if (JOB(3)) then
         if (CSPEC(1) .and. CSPEC(2)) then
            IFLAG(1) = -5
            if (IPRINT.gt.0) write (*,FMT=50)

   50       format (' This problem has continuous spectrum generated by'
     +             ,' both endpoints.  The',/
     +             ' calculation of the spectral density',
     +             ' function has not yet been implemented',/
     +             ' for such cases.',/)

            go to 120

         end if

         if (.not. (CSPEC(1).or.CSPEC(2))) then
            IFLAG(1) = -4
            if (IPRINT.gt.0) write (*,FMT=55)

   55       format (' This problem has no continuous spectrum.')

            go to 120

         end if

         if ((CSPEC(1).and. ((KCLASS(2).eq.5).or.
     +       (KCLASS(2).eq.9))) .or. (CSPEC(2).and.
     +       ((KCLASS(1).eq.5).or. (KCLASS(1).eq.9)))) then
            IFLAG(1) = -6
            if (IPRINT.gt.0) write (*,FMT=60)

   60       format (
     +             ' The normalization of the spectral density function'
     +             ,' is unknown for this problem.')

            go to 120

         end if

         if ((CSPEC(1).and. (.not.BFIN)) .or.
     +       (CSPEC(2).and. (.not.AFIN))) then
            IFLAG(1) = -6
            if (IPRINT.gt.0) write (*,FMT=60)
            go to 120

         end if

         if (OSCILL) then
            FLAG = -25
            if (IPRINT.gt.0) then
               write (*,FMT=30)
               write (21,FMT=30)
            end if

            go to 120

         end if

         XTOL = -LOG10(MAX(TOL(1),TOL(2)))
         JTOL = XTOL - HALF
         JTOL = MIN(MAX(JTOL,1),5)
         DENSOP = 3*JTOL
         MAXITS = (15-JTOL)/3
C
C        Set Maxlvl for the density function calculation; see above
C        "IMPORTANT NOTE" if this is to be increased.
C
         MAXLVL = (7+JTOL)/2
         AAFIN = AFIN
         AA = A
         KCL1 = KCLASS(1)
         BBFIN = BFIN
         BB = B
         KCL2 = KCLASS(2)
         if (JOB(5)) then
C
C           Use interval truncation in this oscillatory regime.
C
            OSCILL = .FALSE.
            if ((.not.JOB(4)) .and. ((KCLASS(1).eq.1).or.
     +          (KCLASS(2).eq.1))) OSCILL = .TRUE.
            if (.not.OSCILL) then
               NXINIT = 4*JTOL + 5
               ENDFAC = ENDI(JTOL)
            else
               NXINIT = 24*JTOL + 36
               ENDFAC = 48.0
            end if

            if (CSPEC(1)) then
               if (AFIN) then
                  KCLASS(1) = 7
                  ENDFAC = 4.0*ENDFAC
                  if (BFIN) then
                     A = AA + (BB-AA)/ENDFAC
                  else
                     A = AA + ABS(AA)/ENDFAC
                  end if

               else
                  AFIN = .TRUE.
                  KCLASS(1) = 0
                  if (BFIN) then
                     A = -ENDFAC - MIN(-B,ZERO)
                  else
                     A = -ENDFAC
                  end if

               end if

            else
               if (BFIN) then
                  KCLASS(2) = 7
                  ENDFAC = 4.0*ENDFAC
                  if (AFIN) then
                     B = BB - (BB-AA)/ENDFAC
                  else
                     B = BB - ABS(BB)/ENDFAC
                  end if

               else
                  BFIN = .TRUE.
                  KCLASS(2) = 0
                  if (AFIN) then
                     B = ENDFAC + MAX(A,ZERO)
                  else
                     B = ENDFAC
                  end if

               end if

            end if

         else
            if (CSPEC(1)) AFIN = .TRUE.
            if (CSPEC(2)) BFIN = .TRUE.
            if (NUMX.eq.0) then
               IFLAG(1) = -30
               return

            end if

         end if

         MAXT = NUMT
C
C        Loop over the choices of intervals.
C
         NUMEV = 0
         do 105 K = 1,MAXITS
            STORE(1) = A
            STORE(NXINIT) = B
            LFLAG(3) = .FALSE.
            FLAG = 0
            if (IPRINT.ge.1) then
               write (*,FMT=65) K
               write (21,FMT=65) K

   65          format (60 ('-'),/' Iteration ',i2)

               write (21,FMT=70) A,B,NXINIT

   70          format (/' For a, b =',2f15.8,/' Nxinit = ',i4,/)

            end if

            if (JOB(5)) then
               I = NXINIT + 16
               call MESH(.TRUE.,-1,STORE,STORE(I+1),STORE(2*I+1),
     +                   STORE(3*I+1),STORE(4*I+1),TOL1,HMIN)
            end if

            if (IPRINT.ge.3) then
               write (*,FMT=75) (STORE(I),I=1,NXINIT)
               write (21,FMT=75) (STORE(I),I=1,NXINIT)

   75          format (' Level 0 mesh:',/ (5g15.6))

            end if

            call DENSEF(TOL,CSPEC,IPRINT,K,NEXTRP,MAXT,T,RHO,IEV,HMIN,
     +                  NUMEV,STORE)
            if (FLAG.eq.-3) then
               LFLAG(6) = .FALSE.
               FLAG = 0
            end if

            if (FLAG.lt.0) go to 120
            if (.not.JOB(5)) go to 110
            if (K.gt.1) then
               DONE = .TRUE.
               J = MAXT
               do 80 I = 1,J
                  RHOTOL = TWO*ZETA*MAX(TOL(1),TOL(2)*RHO(I))
                  ERROR = RHO(I) - STORE(2320+ (MAXLVL+2)*NUMT+I)
                  if (ABS(ERROR).le.RHOTOL) then
                     EDONE = .TRUE.
                  else
                     EDONE = .FALSE.
                     MAXT = I
                  end if

                  DONE = DONE .and. EDONE
   80          continue
               if (DONE) go to 110
            end if

            if (IPRINT.ge.2) then
               write (*,FMT=85)
               write (21,FMT=85)

   85          format (9x,'t',15x,'Truncated Rho(t)')

               do 95 I = 1,NUMT
                  write (*,FMT=90) T(I),RHO(I)
                  write (21,FMT=90) T(I),RHO(I)

   90             format (f12.4,d31.15)

   95          continue
            end if

            do 100 I = 1,MAXT
               STORE(2320+ (MAXLVL+2)*NUMT+I) = RHO(I)
  100       continue
            COUNTZ = .TRUE.
            call SHOOT(CUTOFF,STORE,MU1,FZ)
            call SHOOT(T(NUMT),STORE,MU2,FZ)
            COUNTZ = .FALSE.
            if (T(NUMT).gt.CUTOFF) then
               DENS = (MU2-MU1)/ (K* (T(NUMT)-CUTOFF))
            else
               DENS = DENSOP
            end if

            if (.not.OSCILL) then
               NXINIT = NXINIT + 10
               ZETA = ZETAI(JTOL)
            else
               ZETA = 2.0
               NXINIT = ZETA*NXINIT
            end if

            if (CSPEC(1)) then
               if (AAFIN) then
                  ENDFAC = 5.0*ZETA
                  if (DENS.lt.DENSLO) ENDFAC = 75.0
                  if ((DENS.gt.DENSHI) .and. (.not.OSCILL)) ENDFAC = 8.0
                  A = AA + (A-AA)/ENDFAC
                  if ((AA-A)**2.lt.U) then
                     FLAG = -9
                     go to 110

                  end if

               else
                  if (DENS.lt.DENSLO) ZETA = 2.0
                  if ((DENS.gt.DENSHI) .and. (.not.OSCILL)) ZETA = 1.4
                  A = ZETA*A
               end if

            end if

            if (CSPEC(2)) then
               if (BBFIN) then
                  ENDFAC = 5.0*ZETA
                  if (DENS.lt.DENSLO) ENDFAC = 75.0
                  if ((DENS.gt.DENSHI) .and. (.not.OSCILL)) ENDFAC = 8.0
                  B = BB - (BB-B)/ENDFAC
                  if ((B-BB)**2.lt.U) then
                     FLAG = -9
                     go to 110

                  end if

               else
                  if (DENS.lt.DENSLO) ZETA = 2.0
                  if ((DENS.gt.DENSHI) .and. (.not.OSCILL)) ZETA = 1.4
                  B = ZETA*B
               end if

            end if

            if (MOD(NXINIT,2).eq.0) NXINIT = NXINIT + 1
            NXINIT = MIN(464,NXINIT)
  105    continue
         FLAG = -3
  110    if (IPRINT.ge.1) write (21,FMT=115) NUMEV

  115    format (' The total number of eigenvalues computed was ',i10)

         if (CSPEC(1)) then
            if (AAFIN) then
               A = AA
               KCLASS(1) = KCL1
            else
               AFIN = .FALSE.
            end if

         end if

         if (CSPEC(2)) then
            if (BBFIN) then
               B = BB
               KCLASS(2) = KCL2
            else
               BFIN = .FALSE.
            end if

         end if

      end if
C
C     Set fatal output flags.
C
  120 if (FLAG.lt.-9) then
         do 125 K = 1,MAX(NEV,1)
            IFLAG(K) = FLAG
  125    continue
         return

      else
         if ((FLAG.lt.0) .and. (.not. (JOB(1).or.
     +       JOB(2)))) IFLAG(1) = FLAG
      end if
C
C     Set warning flags.
C
      do 135 I = 2,5
         do 130 K = 1,MAX(NEV,1)
            if (LFLAG(I) .and. (IFLAG(K).ge.0)) IFLAG(K) = IFLAG(K) +
     +          I*IBASE
  130    continue
         if (LFLAG(I)) IBASE = 10*IBASE
  135 continue
      return

      end subroutine SLEDGE
C=======================================================================
      subroutine INTERV(FIRST,ALPHA,BETA,CONS,ENDFIN,NFIRST,NTOTAL,
     +                  IFLAG,X)
***********************************************************************
*                                                                     *
*     INTERV calculates the indices of eigenvalues found in a         *
*     specified interval.                                             *
*                                                                     *
***********************************************************************
C     Local variables:
C
C
C     .. Scalar Arguments ..
      double precision ALPHA,BETA
      integer IFLAG,NFIRST,NTOTAL
      logical FIRST
C     ..
C     .. Array Arguments ..
      double precision CONS(*),X(*)
      logical ENDFIN(*)
C     ..
C     .. Scalars in Common ..
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision CUTOFF,HMIN,V
      integer I1,I2,I3,J1,J2,J3,K,LASTEV,LEVEL0,MU,NEXTRP,NLAST,NUMX
      logical DOMESH
C     ..
C     .. Local Arrays ..
      double precision CEV(2),TOL(6)
      integer IDUMMY(1)
      logical CSPEC(2),JOBST(3),LPLC(2)
C     ..
C     .. External Subroutines ..
cc      external CLASS,MESH,SHOOT,START
C     ..
C     .. Intrinsic Functions ..
      intrinsic MIN
C     ..
C     .. Common blocks ..
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
C      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Save statement ..
      save CUTOFF
C     ..
      NFIRST = -5
      NLAST = -5
      IFLAG = 0
      if (ALPHA.ge.BETA) then
         IFLAG = -11
         return

      end if

      TOL(1) = 0.000001
      TOL(2) = 0.000001
      if (FIRST) then
         JOBST(1) = .FALSE.
         JOBST(2) = .FALSE.
         JOBST(3) = .FALSE.
         AFIN = ENDFIN(1)
         BFIN = ENDFIN(2)
         MAXLVL = 10
         NUMX = 0
         call START(JOBST,CONS,TOL,0,IDUMMY,NUMX,X,0,X,NEXTRP,X)
         call CLASS(0,TOL(1),JOBST(2),CSPEC,CEV,LASTEV,LPLC,X,.TRUE.,
     +              HMIN,DOMESH)
         CUTOFF = MIN(CEV(1),CEV(2))
         if (FLAG.lt.0) then
            IFLAG = FLAG
            return

         end if

         if (OSC(1) .or. OSC(2)) then
            IFLAG = -25
            return

         end if

         if (DOMESH) then
            K = NUMX + 16
            call MESH(.TRUE.,-1,X,X(K),X(2*K+1),X(3*K+1),X(4*K+1),
     +                TOL(1),HMIN)
         end if

      end if

      if (FLAG.lt.0) then
         IFLAG = FLAG
         return

      end if

      LEVEL0 = LEVEL
      COUNTZ = .TRUE.
      LEVEL = 3
      call SHOOT(ALPHA,X,MU,V)
      I1 = MU
      call SHOOT(BETA,X,MU,V)
      J1 = MU
      LEVEL = LEVEL + 1
      call SHOOT(ALPHA,X,MU,V)
      I2 = MU
      call SHOOT(BETA,X,MU,V)
      J2 = MU
   10 LEVEL = LEVEL + 1
      if (NFIRST.eq.-5) then
         call SHOOT(ALPHA,X,MU,V)
         I3 = MU
         if ((I1.eq.I2) .and. (I2.eq.I3)) then
            NFIRST = I1
            go to 15

         end if

         I1 = I2
         I2 = I3
      end if

   15 if (NLAST.eq.-5) then
         call SHOOT(BETA,X,MU,V)
         J3 = MU
         if ((J1.eq.J2) .and. (J2.eq.J3)) then
            NLAST = J1 - 1
            if (NFIRST.ne.-5) go to 20
         end if

         J1 = J2
         J2 = J3
      end if

      if (LEVEL.lt.MAXLVL) go to 10
   20 if (NFIRST.eq.-5) then
         NFIRST = I3
         IFLAG = 12
      end if

      if (NLAST.eq.-5) then
         NLAST = J3
         IFLAG = 12
      end if

      NTOTAL = NLAST + 1 - NFIRST
      if (NTOTAL.eq.0) IFLAG = 11
      if (BETA.gt.CUTOFF) IFLAG = 13
      COUNTZ = .FALSE.
      LEVEL = LEVEL0
      return

      end subroutine INTERV
C========================= End of Part 1 ===============================
*************************   Start of Part 2 ****************************
C///////////////////////////////////////////////////////////////////////
      subroutine AITKEN(XLIM,TOL,N,X,ERROR)
C
C     Use Aitken's algorithm to accelerate convergence of the sequence
C     in X(*).
C
C
C     .. Parameters ..
      double precision ZERO,ONE,TWO
      parameter (ZERO=0.0,ONE=1.0,TWO=2.0)
C     ..
C     .. Scalar Arguments ..
      double precision ERROR,TOL,XLIM
      integer N
C     ..
C     .. Array Arguments ..
      double precision X(*)
C     ..
C     .. Local Scalars ..
      double precision DENOM,XOLD
      integer I
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX
C     ..
      if (N.le.2) then
         XLIM = X(N)
         ERROR = ZERO
         return

      end if

      XOLD = 1.D30
      do 10 I = 1,N - 2
         DENOM = X(I+2) - TWO*X(I+1) + X(I)
         if (DENOM.ne.ZERO) then
            XLIM = X(I) - (X(I+1)-X(I))**2/DENOM
            ERROR = XLIM - XOLD
         else
            ERROR = X(I+2) - X(I+1)
            XLIM = X(I+2)
         end if

         if (ABS(ERROR).lt.MAX(ONE,ABS(XLIM))*TOL) return
         XOLD = XLIM
   10 continue
      return

      end subroutine AITKEN
C----------------------------------------------------------------------
      double precision function ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,
     +                 BETA1,BETA2)
C
C
C     Evaluate the asymptotic formula for eigenvalue NEV.
C        Note: not all cases have been implemented yet.
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,TWO,PI
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,TWO=2.0,
     +          PI=3.14159265358979324D0)
C     ..
C     .. Scalar Arguments ..
      double precision ALPHA1,ALPHA2,BETA1,BETA2,QINT,RPINT
      integer NEV
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision FNEV
C     ..
C     .. Common blocks ..
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      ASYMEV = -999999.0
      FNEV = NEV
      if (REG(1)) then
         if ((A1P.ne.ZERO) .or. (A2P.ne.ZERO)) then
            if (A2P.ne.ZERO) then
               ASYMEV = (TWO*A1P/A2P+QINT)/RPINT
               if (B2.ne.ZERO) then
C                 Case 1
                  ASYMEV = ASYMEV + TWO*B1/ (B2*RPINT) + ((FNEV-ONE)*PI/
     +                     RPINT)**2
               else
C                 Case 2
                  ASYMEV = ASYMEV + ((FNEV-HALF)*PI/RPINT)**2
               end if

            else
               ASYMEV = (TWO*A2/A1P+QINT)/RPINT
               if (B2.ne.ZERO) then
C                 Case 3
                  ASYMEV = ASYMEV + TWO*B1/ (B2*RPINT) +
     +                     ((FNEV-HALF)*PI/RPINT)**2
               else
C                 Case 4
                  ASYMEV = ASYMEV + (FNEV*PI/RPINT)**2
               end if

            end if

         else
            if (A2.ne.ZERO) then
               if (B2.ne.ZERO) then
C                 Case 1
                  ASYMEV = (FNEV*PI/RPINT)**2 +
     +                     (TWO* (BETA1/BETA2+ALPHA1/ALPHA2)+QINT)/RPINT
               else
C                 Case 2   (Dirichlet at B)
                  ASYMEV = ((FNEV+HALF)*PI/RPINT)**2 + (TWO*ALPHA1/
     +                     ALPHA2+QINT)/RPINT
               end if

            else
               if (B2.ne.ZERO) then
C                 Case 3   (Dirichlet at A)
                  ASYMEV = ((FNEV+HALF)*PI/RPINT)**2 + (TWO*BETA1/
     +                     BETA2+QINT)/RPINT
               else
C                 Case 4   (Dirichlet at A and at B)
                  ASYMEV = ((FNEV+ONE)*PI/RPINT)**2 + QINT/RPINT
               end if

            end if

         end if

         return

      end if

      if (REG(2)) then
         if (B2.ne.ZERO) then
            if (A2.ne.ZERO) then
C              Case 1
               ASYMEV = (FNEV*PI/RPINT)**2 +
     +                  (TWO* (ALPHA1/ALPHA2+BETA1/BETA2)+QINT)/RPINT
            else
C              Case 2   (Dirichlet at A)
               ASYMEV = ((FNEV+HALF)*PI/RPINT)**2 + (TWO*BETA1/
     +                  BETA2+QINT)/RPINT
            end if

         else
            if (B2.ne.ZERO) then
C              Case 3   (Dirichlet at B)
               ASYMEV = ((FNEV+HALF)*PI/RPINT)**2 + (TWO*ALPHA1/
     +                  ALPHA2+QINT)/RPINT
            else
C              Case 4   (Dirichlet at A and at B)
               ASYMEV = ((FNEV+ONE)*PI/RPINT)**2 + QINT/RPINT
            end if

         end if

      end if

      return

      end function ASYMEV
C----------------------------------------------------------------------
      double precision function ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
C
C
C     Evaluate the asymptotic formula for RsubNEV.
C        Note: not all cases have been implemented yet.  See the note
C        above in ASYMEV.
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,TWO,PI
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,TWO=2.0,
     +          PI=3.14159265358979324D0)
C     ..
C     .. Scalar Arguments ..
      double precision RPATA,RPATB,RPINT,SCALE
      integer NEV
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision FNEV
C     ..
C     .. Intrinsic Functions ..
      intrinsic MAX
C     ..
C     .. Common blocks ..
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      ASYMR = ZERO
      FNEV = MAX(NEV,2)
      if (REG(1)) then
         if ((A1P.ne.ZERO) .or. (A2P.ne.ZERO)) then
            if (A2P.ne.ZERO) then
               if (B2.ne.ZERO) then
                  ASYMR = TWO*RPINT**3/ (RPATA*A2P**2*
     +                    ((FNEV-ONE)*PI)**4)
               else
                  ASYMR = TWO*RPINT**3/ (RPATA*A2P**2*
     +                    ((FNEV-HALF)*PI)**4)
               end if

            else
               if (B2.ne.ZERO) then
                  ASYMR = RPATA*TWO*RPINT/ (A1P* (FNEV-HALF)*PI)**2
               else
                  ASYMR = RPATA*TWO*RPINT/ (A1P*FNEV*PI)**2
               end if

            end if

         else
            if (A2.ne.ZERO) then
               ASYMR = TWO/ (RPATA*A2*A2*RPINT)
            else
               if (B2.ne.ZERO) then
                  ASYMR = TWO*RPATA* ((FNEV+HALF)*PI/A1)**2/RPINT**3
               else
                  ASYMR = TWO*RPATA* ((FNEV+ONE)*PI/A1)**2/RPINT**3
               end if

            end if

         end if

         return

      end if

      if (REG(2)) then
         if (A2.ne.ZERO) then
            ASYMR = TWO/ (RPATB*B2*B2*RPINT)
         else
            if (B2.ne.ZERO) then
               ASYMR = TWO*RPATB* ((FNEV+HALF)*PI/B1)**2/RPINT**3
            else
               ASYMR = TWO*RPATB* ((FNEV+ONE)*PI/B1)**2/RPINT**3
            end if

         end if

      end if

      ASYMR = ASYMR*SCALE**2
      return

      end function ASYMR
C-----------------------------------------------------------------------
      subroutine BRCKET(N,EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,X)
C
C     Find values for EVLOW and EVHIGH which bracket the Nth eigenvalue;
C     in particular,
C           EV(N-1) < EVLOW < EV(N) < EVHIGH < EV(N+1)  .
C     It is assumed that if U(X,LAMBDA) has NZ zeros in (A,B) then
C           EV(MU-1) < LAMBDA < EV(MU)
C     where MU is a function of NZ, LAMBDA, and the constants in the
C     boundary conditions.  The value of MU for a given LAMBDA is
C     returned by a call to subprogram SHOOT.
C
C
C     Set COUNTZ so that zeros are counted in SHOOT.
C
C     .. Parameters ..
      double precision ZERO,TWO
      parameter (ZERO=0.0,TWO=2.0)
C     ..
C     .. Scalar Arguments ..
      double precision ABSERR,EVHIGH,EVLOW,FHIGH,FLOW,RELERR
      integer N
C     ..
C     .. Array Arguments ..
      double precision X(*)
C     ..
C     .. Scalars in Common ..
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision DIFF,EV,EVSIGN,FEV
      integer K,MU
      logical HIGH,LOW
C     ..
C     .. External Subroutines ..
cc      external SHOOT
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MOD
C     ..
C     .. Common blocks ..
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
C      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      COUNTZ = .TRUE.
      EVSIGN = NSGNF
C
C     SHOOT with Ev = Evlow should return FEV having sign EVSIGN.
C
      if (N.ne.2* (N/2)) EVSIGN = -EVSIGN
      LOW = .FALSE.
      HIGH = .FALSE.
C
C     Make EVLOW a lower bound for EV(N).
C
      EV = EVLOW
      DIFF = ABS(EVHIGH-EVLOW)
      if (DIFF.eq.ZERO) DIFF = ABSERR + RELERR
   10 call SHOOT(EV,X,MU,FEV)
      if (FLAG.lt.0) then
         COUNTZ = .FALSE.
         return

      end if

      if (MU.gt.N) then
         EVHIGH = EV
         FHIGH = FEV
         EV = EV - DIFF
         DIFF = TWO*DIFF
         if ((MU.eq.N+1) .and. (EVSIGN*FEV.le.ZERO)) HIGH = .TRUE.
         go to 10

      else
         EVLOW = EV
         FLOW = FEV
         if ((MU.eq.N) .and. (EVSIGN*FEV.ge.ZERO)) LOW = .TRUE.
      end if
C
C     Make EVHIGH an upper bound for EV(N).
C
      if (.not.HIGH) then
         EV = EVHIGH
         DIFF = ABS(EVHIGH-EVLOW)
   20    call SHOOT(EV,X,MU,FEV)
         if (FLAG.lt.0) return
         K = NSGNF* (-1)**MOD(MU,2)
         if (K*FEV.lt.ZERO) then
            EV = EV + DIFF
            DIFF = TWO*DIFF
            go to 20

         else
            if (MU.le.N) then
               EVLOW = EV
               FLOW = FEV
               EV = EV + DIFF
               DIFF = TWO*DIFF
               go to 20

            else
               EVHIGH = EV
               FHIGH = FEV
               if (MU.eq.N+1) HIGH = .TRUE.
            end if

         end if

      end if
C
C     Refine the interval [EVLOW,EVHIGH] to include only the Nth
C     eigenvalue.
C
   30 if ((.not.LOW) .or. (.not.HIGH)) then
         DIFF = EVHIGH - EVLOW
         EV = EVLOW + DIFF/TWO
C
C        Check for a cluster of eigenvalues within user's tolerance.
C
         if (TWO*DIFF.lt.MAX(ABSERR,RELERR* (MAX(ABS(EVLOW),
     +       ABS(EVHIGH))))) then
            LFLAG(1) = .TRUE.
            COUNTZ = .FALSE.
            return

         end if

         call SHOOT(EV,X,MU,FEV)
         if (FLAG.lt.0) then
            COUNTZ = .FALSE.
            return

         end if
C
C        Update EVLOW and EVHIGH.
C
         if (MU.eq.N) then
            EVLOW = EV
            LOW = .TRUE.
            FLOW = FEV
         else
            if (MU.eq.N+1) then
               EVHIGH = EV
               HIGH = .TRUE.
               FHIGH = FEV
            else
               if (MU.lt.N) then
                  EVLOW = EV
                  FLOW = FEV
               else
                  EVHIGH = EV
                  FHIGH = FEV
               end if

            end if

         end if

         go to 30

      end if

      COUNTZ = .FALSE.
      return

      end subroutine BRCKET
C-----------------------------------------------------------------------
      subroutine CLASS(IPRINT,TOL,JOB,CSPEC,CEV,LASTEV,LPLC,X,JMESH,
     +                 HMIN,DOMESH)
C
C     This routine classifies the Sturm-Liouville problem.  Note:
C     (1) any computational algorithm must be based on a finite
C     amount of information; hence, there will always be cases that
C     any algorithm misclassifies.  In addition, some problems are
C     inherently ill-conditioned, in that a small change in the
C     coefficients can produce a totally different classification.
C     (2)  The maximum number of points sampled for singular problems
C     is given by the variable KMAX.  By increasing this number, the
C     reliability of the classification may increase; however, the
C     computing time may also increase.  The values we have chosen
C     seem to be a reasonable balance for most problems.
C     (3) The algorithms apply standard theorems involving limits of
C     the Liouville normal form potential.  When this is not available,
C     each coefficient function is approximated by a power function
C     (c*x^r) and classified according to the properties of the
C     resulting Liouville approximation.
C
C
C     Sample the coefficient functions; determine if the problem as
C     given is in Liouville normal form.
C
C     .. Parameters ..
      double precision ZERO,TENTH,ONE,TWO,EIGHT
      parameter (ZERO=0.0,TENTH=0.1D0,ONE=1.0,TWO=2.0,EIGHT=8.0)
C     ..
C     .. Scalar Arguments ..
      double precision HMIN,TOL
      integer IPRINT,LASTEV
      logical DOMESH,JMESH,JOB
C     ..
C     .. Array Arguments ..
      double precision CEV(*),X(*)
      logical CSPEC(*),LPLC(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision BASE,END,EV,FEV,OVER,S,SGN
      integer IFLAG,J,K,KMAX,M
C     ..
C     .. Local Arrays ..
      double precision BC(2),PZ(40,2),QZ(40,2),RZ(40,2),Y(40),Z(40,2)
      integer KUSED(2),LAST(3)
      logical ENDFIN(2)
C     ..
C     .. External Subroutines ..
cc      external CLSEND,COEFF,MESH,SHOOT
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,INT,LOG,LOG10,MAX,MIN
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      LNF = .TRUE.
      DOMESH = .TRUE.
      do 40 J = 1,2
         if (J.eq.1) then
            END = A
            ENDFIN(1) = AFIN
            SGN = ONE
         else
            END = B
            ENDFIN(2) = BFIN
            SGN = -ONE
         end if

         K = TENTH - LOG10(TOL)
         KMAX = MIN(MAX(4*K,10),40)
         M = 0
         BASE = ONE/MIN(MAX(K,4),8)
         if (ENDFIN(J)) then
            if (END.ne.ZERO) KMAX = MIN(-INT(LOG10(U))-1,KMAX)
         else
            KMAX = MIN(KMAX,20)
            BASE = EIGHT
         end if

         do 10 K = 1,KMAX
            Z(K,J) = BASE**K
   10    continue
         OVER = UNDER/TWO
         do 30 K = 1,KMAX
            if (ENDFIN(J)) then
               S = END + SGN*Z(K,J)
            else
               S = -SGN*Z(K,J)
            end if

            call COEFF(S,PZ(K,J),QZ(K,J),RZ(K,J))
            NCOEFF = NCOEFF + 1
            if ((PZ(K,J).le.ZERO) .or. (RZ(K,J).le.ZERO)) then
               FLAG = -15
               return

            end if

            if (LNF .and. (K.gt.1)) then
               if (PZ(K,J).ne.PZ(K-1,J)) LNF = .FALSE.
               if (RZ(K,J).ne.RZ(K-1,J)) LNF = .FALSE.
            end if

            S = LOG(PZ(K,J))
            if (ABS(S).gt.OVER) M = 1
            S = LOG(ONE+ABS(QZ(K,J)))
            if (ABS(S).gt.OVER) M = 1
            S = LOG(RZ(K,J))
            if (ABS(S).gt.OVER) M = 1
            if (M.ne.0) go to 35
   30    continue
         K = KMAX
   35    KUSED(J) = K - 1
   40 continue
      do 50 J = 1,2
         call CLSEND(Z(1,J),PZ(1,J),QZ(1,J),RZ(1,J),KUSED(J),IPRINT,END,
     +               ENDFIN(J),J,TOL,CEV(J),CSPEC(J),BC,Y,LPLC,IFLAG)
         if (.not.JOB) then
            if ((.not.AFIN) .and. (J.eq.1)) X(1) = -END
            if ((.not.BFIN) .and. (J.eq.2)) X(NXINIT) = END
         end if

         if ((CSPEC(J).or. (.not.OSC(J))) .and. (.not.REG(J))) then
            if (J.eq.1) then
               A1 = BC(1)
               A1P = ZERO
               A2 = BC(2)
               A2P = ZERO
            else
               B1 = BC(1)
               B2 = BC(2)
            end if

         end if

         if (IFLAG.eq.1) LFLAG(2) = .TRUE.
   50 continue
      CUTOFF = MIN(CEV(1),CEV(2))
C
C     Find the number of eigenvalues below the start of the
C     continuous spectrum.
C
      LASTEV = -5
      if ((.not.OSC(1)) .and. (.not.LC(2)) .and. OSC(2)) LASTEV = 0
      if ((.not.OSC(2)) .and. (.not.LC(1)) .and. OSC(1)) LASTEV = 0
      if (OSC(1) .and. (.not.LC(1)) .and. OSC(2) .and. LC(2)) LASTEV = 0
      if (OSC(2) .and. (.not.LC(2)) .and. OSC(1) .and. LC(1)) LASTEV = 0
      if ((CSPEC(1).and.OSC(2)) .or. (CSPEC(2).and.OSC(1))) LASTEV = 0
      if ((CSPEC(1).and. (.not.OSC(2))) .or.
     +    (CSPEC(2).and. (.not.OSC(1))) .or.
     +    (CSPEC(1).and.CSPEC(2))) then
         K = NXINIT + 16
         call MESH(JMESH,-1,X,X(K),X(2*K+1),X(3*K+1),X(4*K+1),TOL,HMIN)
         DOMESH = .FALSE.
         COUNTZ = .TRUE.
         do 75 J = 1,2
            LEVEL = 3*J
            EV = CUTOFF
            call SHOOT(EV,X,LAST(J),FEV)
            if (FLAG.lt.0) return
            if (IPRINT.ge.5) write (21,FMT=70) LEVEL,LAST(J)

   70       format (' When level = ',i2,', Ev index at cutoff is ',i12)

   75    continue
         COUNTZ = .FALSE.
         if (LAST(1).ge.LAST(2)) then
            LASTEV = LAST(2)
            if (LAST(1).ne.LAST(2)) then
               LFLAG(2) = .TRUE.
               if (IPRINT.ge.3) write (21,FMT=80)

   80          format (' The eigenvalue count is uncertain.')

               if (LAST(1).gt.2*LAST(2)) LASTEV = 0
            end if

         else
            LASTEV = -5
         end if

      end if

      return

      end subroutine CLASS
C-----------------------------------------------------------------------
      subroutine CLSEND(Z,PZ,QZ,RZ,KMAX,IPRINT,END,ENDFIN,IEND,TOL,CEV,
     +                  CSPEC,BC,Y,LPLC,IFLAG)
C
C     Iflag = 0    if reasonably certain of the classification;
C           = 1    if not sure.
C
C     Information about the nature of the problem at singular point
C     IEND is passed through the variable KCLASS(IEND):
C         KCLASS(*) = 0    normal;
C                   = 1    oscillatory coefficient function;
C                   = 2    regular, but 1/p, q, or r unbounded;
C                   = 3    infinite endpoint, Eqlnf = -1 ;
C                   = 4    finite singular endpoint, Tau unbounded,
C                          (not 8-10);
C                   = 5    not "hard", irregular;
C                   = 6    "hard" irregular with Eta(1) < 0;
C                   = 7    finite end which generates Cspectrum;
C                   = 8    Q is unbounded (< 1/t^2) near a nonoscill-
C                          atory finite end;
C                   = 9    Q is unbounded (like 1/t^2) near a nonosc-
C                          illatory finite end;
C                   = 10   "hard", irregular, Eta(1) > 0.
C                    Note: "hard" means Tau goes to +infinity at a
C                          finite nonoscillatory endpoint.
C         REG(*) = .True.  iff endpoint is regular.
C         LC(*)  = .True.  iff endpoint is limit circle.
C         OSC(*) = .True.  iff endpoint is oscillatory for all Ev.
C         CSPEC  = .True.  iff endpoint generates continuous spectrum.
C         LPLC(*)= .True.  iff theory yields Lp/Lc classification.
C
C
C     .. Parameters ..
      double precision ZERO,QUART,HALF,QUART3,ONE,TWO,FOUR
      parameter (ZERO=0.0,QUART=0.25D0,HALF=0.5D0,QUART3=0.75D0,ONE=1.0,
     +          TWO=2.0,FOUR=4.0)
C     ..
C     .. Scalar Arguments ..
      double precision CEV,END,TOL
      integer IEND,IFLAG,IPRINT,KMAX
      logical CSPEC,ENDFIN
C     ..
C     .. Array Arguments ..
      double precision BC(*),PZ(*),QZ(*),RZ(*),Y(*),Z(*)
      logical LPLC(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CPT(2),CRT(2),D(4,2),EMU(2),EPT(2),EQLNF(2),
     +                 ERT(2),ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision C1,C2,C3,CP,CQ,CR,DELTA,EP,EQ,ER,GAMMA,SGN,TOL4,
     +                 ZZ
      integer I,IQLNF,K
      logical EX,EXACT,IRREG,POSC,QOSC,ROSC
C     ..
C     .. External Subroutines ..
cc      external COEFF,POWER
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN,SIGN,SQRT
C     ..
C     .. Common blocks ..
      common /SLCLSS/CPT,CRT,CUTOFF,D,EMU,EPT,EQLNF,ERT,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      IFLAG = 0
      if (IPRINT.ge.3) then
         if (IEND.eq.1) write (21,FMT=5)

    5    format (/' For endpoint A')

         if (IEND.eq.2) write (21,FMT=6)

    6    format (/' For endpoint B:')

         write (21,FMT=7) KMAX

    7    format ('  Kmax = ',i3)

      end if

      CSPEC = .FALSE.
      IRREG = .FALSE.
      LPLC(IEND) = .TRUE.
      KCLASS(IEND) = 0
      PNU(IEND) = ZERO
      CEV = ONE/U
      EX = .TRUE.
      C1 = ZERO
      C2 = ZERO
      TOL4 = FOUR*TOL
C
C     Seek monomial approximations to each coefficient function.
C
      call POWER(Z,QZ,KMAX,TOL,IPRINT,EQ,CQ,QOSC,EXACT,Y,IFLAG)
      if (ABS(CQ).le.TOL) CQ = ZERO
      if (CQ.eq.ZERO) EQ = ZERO
      if (LNF) then
         EP = ZERO
         CP = PZ(1)
         ER = ZERO
         CR = RZ(1)
         POSC = .FALSE.
         ROSC = .FALSE.
      else
         call POWER(Z,PZ,KMAX,TOL,IPRINT,EP,CP,POSC,EXACT,Y,IFLAG)
         if (ABS(CP).le.TOL) EP = ZERO
         EX = EX .and. EXACT
         call POWER(Z,RZ,KMAX,TOL,IPRINT,ER,CR,ROSC,EXACT,Y,IFLAG)
         if (ABS(CR).le.TOL) ER = ZERO
      end if

      if (POSC .or. ROSC) then
         if (ENDFIN) then
            REG(IEND) = .TRUE.
         else
            IFLAG = 1
            if (IPRINT.ge.3) write (21,FMT=10)

   10       format (
     +             '  WARNING: p(x) or r(x) is not well-approximated by'
     +             ,' a power potential.',/
     +             '  Classification is uncertain.')

            REG(IEND) = .FALSE.
            KCLASS(IEND) = 1
         end if

         LC(IEND) = .TRUE.
         OSC(IEND) = .FALSE.
      end if

      if (QOSC) then
         IFLAG = 1
         if (IPRINT.ge.3) write (21,FMT=20)

   20    format ('  WARNING: q(x) is not well-approximated by a power ',
     +          'potential.',/'  Classification is uncertain.')

         if (ENDFIN) then
            REG(IEND) = .TRUE.
            LC(IEND) = .TRUE.
            OSC(IEND) = .FALSE.
         else
            KCLASS(IEND) = 1
            REG(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CSPEC = .TRUE.
            OSC(IEND) = .FALSE.
            BC(1) = ONE
            BC(2) = ZERO
            CEV = QZ(KMAX-1)
            K = 40
            DELTA = (Z(KMAX)-Z(KMAX-1))/ (K+1)
            do 30 I = 0,K
               ZZ = Z(KMAX) - I*DELTA
               call COEFF(ZZ,CP,CQ,CR)
               NCOEFF = NCOEFF + 1
               CEV = MIN(CEV,CQ)
   30       continue
            if (ABS(CEV).lt.TOL4) CEV = ZERO
            if (U*ABS(CEV).ge.ONE) CSPEC = .FALSE.
         end if

         EQLNF(IEND) = ZERO
      end if

      if (IPRINT.ge.3) then
         write (21,FMT=40) CP,EP,CQ,EQ,CR,ER

   40    format ('  Cp, Ep; Cq, Eq; Cr, Er =',/3 (2d25.12,/))

      end if

      CPT(IEND) = CP
      EPT(IEND) = EP
C
C     Analyze this endpoint.
C
      if ((EP.lt.ONE) .and. (EQ.gt.-ONE) .and. (ER.gt.-ONE) .and.
     +    ENDFIN) then
         REG(IEND) = .TRUE.
         LC(IEND) = .TRUE.
         OSC(IEND) = .FALSE.
         CSPEC = .FALSE.
         if ((EP.gt.ZERO) .or. (EQ.lt.ZERO) .or.
     +       (ER.lt.ZERO)) KCLASS(IEND) = 2
         return

      end if

      REG(IEND) = .FALSE.
      ETA(1,IEND) = HALF* (ER-EP+TWO)
      if (ABS(ETA(1,IEND)).le.TOL) ETA(1,IEND) = ZERO
      if (ETA(1,IEND).ne.ZERO) then
         EQLNF(IEND) = (EQ-ER)/ETA(1,IEND)
         IQLNF = EQLNF(IEND) + SIGN(HALF,EQLNF(IEND))
         if (ABS(IQLNF-EQLNF(IEND)).lt.TOL4) EQLNF(IEND) = IQLNF
         C1 = (CQ/CR)* (ABS(ETA(1,IEND))*SQRT(CP/CR))**EQLNF(IEND)
         if (C1.eq.ZERO) EQLNF(IEND) = ZERO
         C2 = (EP+ER)/ (FOUR*ETA(1,IEND))
         C2 = C2* (C2-ONE)
         if (IPRINT.ge.5) then
            write (21,FMT=50) C1,C2,EQLNF(IEND),ETA(1,IEND)

   50       format ('  C1, C2      =',2d20.10,/
     +             '  Eqlnf, Eta1 =',2d20.10)

         end if

      else
         C3 = (CP/CR)* (QUART* (EP+ER))**2
         if (IPRINT.ge.5) then
            write (21,FMT=60) C3

   60       format ('  C3 = ',d19.10)

         end if

      end if

      if (.not.ENDFIN) then
C
C        Make an initial estimate for "infinity" (used in MESH).
C
         if ((ER.gt.EQ) .or. (CQ.eq.ZERO)) then
            if (ER.ne.EP) then
               GAMMA = ER - EP
               DELTA = CR/CP
            else
               GAMMA = EQ - EP
               DELTA = ABS(CQ)/CP
               if (DELTA.eq.ZERO) DELTA = ONE
            end if

         else
            if (EQ.gt.ER) then
               GAMMA = EQ - EP
               DELTA = ABS(CQ)/CP
            else
               if (ER.gt.EP) then
                  GAMMA = ER - EP
                  DELTA = CR/CP
               else
                  GAMMA = ZERO
                  if (ER.eq.EP) then
                     DELTA = ABS(CR-CQ)/CP
                  else
                     DELTA = ABS(CQ)/CP
                  end if

               end if

            end if

         end if

         if (GAMMA.gt.HALF) then
            if (GAMMA.lt.TWO) then
               END = 80.0
            else
               END = MIN(MAX(64.0/ ((TWO*GAMMA-3.0)*DELTA** (ONE/
     +               (GAMMA+TWO))),ONE),80.D0)
               if (GAMMA.gt.24.0) END = 12.0
            end if

         else
            if (GAMMA.lt.-HALF) then
               END = MAX(MIN(600.0*DELTA** (ONE/
     +               GAMMA)*5.0**GAMMA,120.0D0),TWO)
            else
               END = 12.0
            end if

            if ((GAMMA.eq.ZERO) .and. (CQ.ne.ZERO)) END = 40.0
         end if

      end if
C
C     Test for finite irregular singular points.
C
      if (ENDFIN) then
         SGN = ONE
         I = ER - EP + SIGN(HALF,ER-EP)
         K = EQ - EP + SIGN(HALF,EQ-EP)
         IRREG = .TRUE.
         if (CQ.eq.ZERO) then
            if ((I.ge.-2) .and. (ABS(ER-EP-I).le.TOL4)) IRREG = .FALSE.
         else
            if ((ER.le.EQ) .and. (I.ge.-2) .and.
     +          (ABS(ER-EP-I).le.TOL4)) IRREG = .FALSE.
            if ((ER.gt.EQ) .and. (K.ge.-2) .and.
     +          (ABS(EQ-EP-K).le.TOL4)) IRREG = .FALSE.
         end if

         EMU(IEND) = HALF* (ONE-EP)
         if (IRREG) then
            if (IPRINT.ge.3) write (21,FMT=70)

   70       format ('  This is an irregular singular point.')

            KCLASS(IEND) = 5
         else
C
C           Compute the principal Frobenius root.
C
            if (ETA(1,IEND).ne.ZERO) then
               if ((CQ.ne.ZERO) .and. (ER.gt.EQ) .and. (K.eq.-2)) then
                  PNU(IEND) = EMU(IEND)**2 + CQ/CP
                  if (ABS(PNU(IEND)).le.TOL4) PNU(IEND) = ZERO
                  if (PNU(IEND).ge.ZERO) then
                     PNU(IEND) = EMU(IEND) + SQRT(PNU(IEND))
                  else
                     PNU(IEND) = -EP
                  end if

               else
                  PNU(IEND) = MAX(ONE-EP,ZERO)
               end if

               if (PNU(IEND).gt.-EP) then
                  if (IPRINT.ge.5) write (21,FMT=75) PNU(IEND)

   75             format ('  The principal Frobenius root is ',e20.8)

               end if

            end if

         end if

      else
         SGN = -ONE
      end if

      if (SGN*ETA(1,IEND).gt.ZERO) then
C
C        Carry out the Case 1 tests.
C
         K = 0
         if (EQLNF(IEND).lt.-TWO) then
            if (CQ.lt.ZERO) K = 1
            if (CQ.gt.ZERO) K = -1
         end if

         if (EQLNF(IEND).eq.-TWO) then
            if (ABS(C1+C2+QUART).le.TOL4) then
               if (IPRINT.ge.3) write (21,FMT=80)

   80          format ('  WARNING: borderline nonoscillatory/oscillato',
     +                'ry classification.')

               K = -1
               IFLAG = 1
            else
               if (C1+C2.lt.-QUART-TOL4) K = 1
               if (C1+C2.gt.-QUART) K = -1
            end if

         end if

         if (EQLNF(IEND).gt.-TWO) then
            if (ABS(C2+QUART).le.TOL4) then
               C2 = -QUART
               if (IPRINT.ge.3) write (21,FMT=80)
               IFLAG = 1
            end if

            if (C2.ge.-QUART) K = -1
         end if

         if (K.eq.1) then
            OSC(IEND) = .TRUE.
         else
            if (K.eq.-1) then
               OSC(IEND) = .FALSE.
            else
               if (IPRINT.ge.3) write (21,FMT=85)

   85          format ('  NO INFORMATION on osc/nonosc class.')

            end if

         end if

         K = 0
         if (EQLNF(IEND).lt.-TWO) then
            if (CQ.gt.ZERO) K = -1
            if (CQ.lt.ZERO) K = 1
         end if

         if (EQLNF(IEND).eq.-TWO) then
            if (ABS(C1+C2-QUART3).le.TOL4) then
               K = -1
               if (IPRINT.ge.3) write (21,FMT=90)

   90          format ('  WARNING: borderline Lc/Lp classification.')

               IFLAG = 1
            end if

            if (C1+C2.ge.QUART3) K = -1
            if (ABS(C1+C2).lt.QUART3-TOL4) K = 1
            if (C1+C2.lt.-TOL4) K = 1
         end if

         if (EQLNF(IEND).gt.-TWO) then
            if (ABS(C2-QUART3).le.TOL4) then
               K = -1
               if (IPRINT.ge.3) write (21,FMT=90)
               IFLAG = 1
            end if

            if (C2.ge.QUART3) K = -1
            if (ABS(C2).lt.QUART3-TOL4) K = 1
            if (C2.lt.-TOL4) K = 1
         end if

         if (K.eq.1) then
            LC(IEND) = .TRUE.
         else
            if (K.eq.-1) then
               LC(IEND) = .FALSE.
            else
               write (21,FMT=95)

   95          format ('  NO INFORMATION on Lp/Lc class.')

            end if

         end if

      end if

      if (SGN*ETA(1,IEND).lt.ZERO) then
C
C        Carry out the Case 2 tests.
C
         K = 0
         if ((EQLNF(IEND).gt.ZERO) .and. (CQ.lt.ZERO)) K = 1
         if ((EQLNF(IEND).gt.ZERO) .and. (CQ.gt.ZERO)) K = -1
         if (EQLNF(IEND).eq.ZERO) then
            K = -1
            CEV = CQ/CR
            CSPEC = .TRUE.
            if (U*ABS(CEV).ge.ONE) CSPEC = .FALSE.
         end if

         if (EQLNF(IEND).lt.ZERO) then
            K = -1
            CEV = ZERO
            CSPEC = .TRUE.
         end if

         if (K.eq.1) then
            OSC(IEND) = .TRUE.
         else
            if (K.eq.-1) then
               OSC(IEND) = .FALSE.
            else
               write (21,FMT=100)

  100          format ('  NO INFORMATION on Osc/Nonosc class.')

            end if

         end if

         K = 0
         if ((EQLNF(IEND).gt.TWO) .and. (CQ.gt.ZERO)) K = -1
         if (EQLNF(IEND).le.TWO) K = -1
         if ((EQLNF(IEND).gt.TWO) .and. (CQ.lt.ZERO)) K = 1
         if (K.eq.1) then
            LC(IEND) = .TRUE.
         else
            if (K.eq.-1) then
               LC(IEND) = .FALSE.
            else
               write (21,FMT=105)

  105          format ('  NO INFORMATION on Lp/Lc class.')

            end if

         end if

      end if

      if (ETA(1,IEND).eq.ZERO) then
C
C        Carry out the Case 3 and 4 tests.
C
         if ((SGN* (EQ-ER).lt.ZERO) .and. (CQ.lt.ZERO)) then
            OSC(IEND) = .TRUE.
            LC(IEND) = .TRUE.
         end if

         if ((SGN* (EQ-ER).lt.ZERO) .and. (CQ.gt.ZERO)) then
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
         end if

         if (EQ.eq.ER) then
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CEV = CQ/CR + C3
            CSPEC = .TRUE.
            if (U*ABS(CEV).ge.ONE) CSPEC = .FALSE.
         end if

         if ((SGN* (EQ-ER).gt.ZERO) .or. (CQ.eq.ZERO)) then
            OSC(IEND) = .FALSE.
            LC(IEND) = .FALSE.
            CEV = C3
            CSPEC = .TRUE.
            if (U*ABS(CEV).ge.ONE) CSPEC = .FALSE.
         end if

      end if

      if (ABS(CEV).le.TOL4) CEV = ZERO
C
C     Calculate the Friedrichs boundary condition (if appropriate).
C
      if (CSPEC) then
         BC(1) = ONE
         BC(2) = ZERO
      end if

      if ((.not.CSPEC) .and. (.not.OSC(IEND))) then
         if ((SGN* (ER+EP).gt.ZERO) .and. (SGN* (EQ+EP).gt.ZERO)) then
            BC(1) = ZERO
            BC(2) = ONE
         else
            if ((SGN* (ER+EP).gt.ZERO) .and. (EQ+EP.eq.ZERO)) then
               BC(1) = SQRT(CP*ABS(CQ))
               BC(2) = ONE
               if (BC(1).gt.ONE) then
                  BC(2) = ONE/BC(1)
                  BC(1) = ONE
               end if

            else
               if ((SGN* (ER+EP).lt.ZERO) .or.
     +             (SGN* (EQ+EP).lt.ZERO)) then
                  BC(1) = ONE
                  BC(2) = ZERO
               end if

            end if

         end if

      end if

      if (.not.OSC(IEND)) then
         if ((.not.ENDFIN) .and. (EQLNF(IEND).eq.-ONE)) KCLASS(IEND) = 3
         I = ER - EP
         if (CQ.ne.ZERO) then
            K = EQ - EP
            if (CQ.gt.ZERO) then
               if (K.lt.I) I = 0
            else
               I = MIN(I,K)
            end if

         end if

         if (ENDFIN .and. (I.lt.0)) then
            if (IRREG) KCLASS(IEND) = 6
            if (ETA(1,IEND).gt.ZERO) then
C
C              Transform some nonoscillatory problems for which Tau
C              is unbounded near a finite singular endpoint.
C
               if (IRREG) KCLASS(IEND) = 10
               CPT(IEND) = CP
               CRT(IEND) = CR
               EPT(IEND) = EP
               ERT(IEND) = ER
               EMU(IEND) = ZERO
               D(1,IEND) = (ETA(1,IEND)*SQRT(CP/CR))** (ONE/ETA(1,IEND))
               D(2,IEND) = ONE/SQRT(SQRT(CP*CR*D(1,IEND)** (EP+ER)))
               if ((EQLNF(IEND).eq.-TWO) .or. (C1.eq.ZERO)) then
                  if (.not.IRREG) KCLASS(IEND) = 9
                  EMU(IEND) = ABS(QUART+C1+C2)
                  if (EMU(IEND).lt.TOL4) then
                     EMU(IEND) = HALF
                  else
                     EMU(IEND) = HALF + SQRT(EMU(IEND))
                  end if

               else
                  if (.not.IRREG) KCLASS(IEND) = 8
               end if

               ETA(2,IEND) = EMU(IEND) - QUART* (EP+ER)/ETA(1,IEND)
               if ((KCLASS(IEND).eq.10) .and. (EMU(IEND).eq.ZERO)) then
                  ETA(2,IEND) = HALF* (ONE-EP)/ETA(1,IEND)
                  EMU(IEND) = ETA(2,IEND) + QUART* (EP+ER)/ETA(1,IEND)
               end if

               D(3,IEND) = ETA(2,IEND)* (ETA(2,IEND)+ (EP-ONE)/
     +                     ETA(1,IEND))
               D(4,IEND) = D(3,IEND)
               if (EQLNF(IEND).eq.-TWO) D(4,IEND) = D(4,IEND) - C1
               if (ABS(D(4,IEND)).le.TOL4) D(4,IEND) = ZERO
               D(4,IEND) = SQRT(ABS(D(4,IEND)))
               if (IPRINT.ge.5) then
                  write (21,FMT=110) EMU(IEND),ETA(2,IEND)
                  write (21,FMT=115) D(3,IEND),D(4,IEND)

  110             format ('  Mu =',d20.6,';  Eta2 =',d20.6)
  115             format ('  D3 =',d20.6,';    D4 =',d20.6)

               end if

            end if

            if (KCLASS(IEND).ge.9) then
               if (.not.EX) LFLAG(5) = .TRUE.
               if (IPRINT.ge.3) then
                  write (21,FMT=120)

  120             format ('  This problem has unbounded ',
     +                   '[Ev*r(x)-q(x)]/p(x).')

                  if ((EMU(IEND).gt.ZERO) .or.
     +                (EP+ER.ne.ZERO)) write (21,FMT=125)

  125             format ('  A change of variables will be used near',
     +                   ' this endpoint.')

               end if

            end if

         end if

         if (ENDFIN .and. (.not.REG(IEND)) .and.
     +       (KCLASS(IEND).eq.0)) KCLASS(IEND) = 4
      end if

      if ((POSC.or.QOSC.or.ROSC) .and. (.not.ENDFIN)) END = 99.0
      if (IPRINT.ge.5) then
         write (21,FMT=130) KCLASS(IEND)

  130    format ('  Classification type (KCLASS) is: ',i2)

      end if

      return

      end subroutine CLSEND
C-----------------------------------------------------------------------
      subroutine DENSEF(TOL,CSPEC,IPRINT,ITER,NEXTRP,NUMT,T,RHO,NEV,
     +                  HMIN,NUMEV,STORE)
***********************************************************************
*                                                                     *
*     This routine computes the spectral density function rho(t).     *
*                                                                     *
***********************************************************************
C
C    Input parameters:
C      TOL(*) as in SLEDGE.
C      CSPEC(*)  = logical 2-vector; CSPEC(i) = .true. iff endpoint i
C                  (1 = A, 2 = B) generates continuous spectrum.
C      IPRINT    = integer controlling printing.
C      ITER      = iteration from SLEDGE.
C      NEXTRP    = integer giving maximum no. of extrapolations.
C      NUMT      = integer equalling number of T(*) points.
C      T(*)      = real vector of abcissae for spectral function rho(t).
C      HMIN      = minimum stepsize in Level 0 mesh.
C
C    Output parameters:
C      RHO(*)    = real vector of values for spectral density function
C                  rho(t),  RHO(I) = rho(T(I)).
C      NEV       = integer pointer to eigenvalue.  On a normal return
C                  (FLAG = 0) this is set to the index of the last
C                  eigenvalue computed; if FLAG is not zero, then NEV
C                  gives the index of the eigenvalue where the problem
C                  occurred.
C      NUMEV     = cumulative number of eigenvalues computed.
C
C    Auxiliary storage:
C      STORE(*) = real vector of auxiliary storage, must be dimensioned
C                 at least 5*Nxinit+(Maxlvl+2)*NUMT.  The value of
C                 Nxinit is either the input NUMX or Maxint.  Currently,
C                 Maxlvl = 8 and Maxint = 235.
C            1    ->     Nxinit         vector of mesh points X(*),
C       Nxinit+1  ->    5*Nxinit        intermediate RsubN calculations,
C     5*Nxinit+1  -> 5*Nxinit+10*NUMT   intermediate RHO values.
C-----------------------------------------------------------------------
C     The definition of a spectral density function assumes a certain
C     normalization on the eigenfunctions.  For the case when x = b
C     generates the continuous spectrum, the normalization used here
C     (and in routine GETRN below) is:
C     (1) when x = a is regular  u(a) = (A2-A2P*Ev)/SCALE
C                            (pu')(a) = (A1-A1P*Ev)/SCALE,
C         with SCALE = sqrt(A1**2+A2**2) when A1' = A2' = 0, and
C              SCALE = sqrt(ALPHA) otherwise.
C     (2) When x = a is a regular singular point then u(x) is taken to
C         be asymptotic to the principal Frobenius solution, i.e.,
C         near x = a   u(x) ~ (x-a)**Nu    with Nu the larger
C         root of the indicial equation.
C     Analogous normalizations hold at x = b when the endpoint x = a
C     generates the continuous spectrum.
C----------------------------------------------------------------------
C     Local variables:
C
C
C     Initialization:
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,TWO,FOUR,SIX,TEN,TOLMIN
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,TWO=2.0,FOUR=4.0,SIX=6.0,
     +          TEN=10.0,TOLMIN=5.D-3)
C     ..
C     .. Scalar Arguments ..
      double precision HMIN
      integer IPRINT,ITER,NEV,NEXTRP,NUMEV,NUMT
C     ..
C     .. Array Arguments ..
      double precision RHO(*),STORE(*),T(*),TOL(*)
      logical CSPEC(2)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETAT(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision ABSERR,AEV,ALPHA,ALPHA1,ALPHA2,ARN,AVGRHO,BETA1,
     +                 BETA2,BIG,DELTA,DENOM,DX,EFNORM,ERROR,ETA,ETAOLD,
     +                 EV,EVHAT,EVHIGH,EVLOW,EVSAVE,FHIGH,FLOW,H,HALFH,
     +                 PDU,PSRHO,PX,QINT,QLNF,QX,RELERR,RHOSUM,ROLD,
     +                 RPATA,RPATB,RPINT,RSUBN,RX,SCALE,SQRTRP,TENU,
     +                 TOL1,TOL2,TOLMAX,UX,XLEFT,XTOL,Z,ZABS,ZREL
      integer I,IASYMP,J,JS,KLVL,KRHO,LPRINT,MAXNEV,NRHO,NSAVE
      logical DONE,JUMP,RDONE
C     ..
C     .. Local Arrays ..
      double precision EVOLD(200),RSAVE(5)
      logical EFIN(2,2)
C     ..
C     .. External Functions ..
cc      double precision ASYMEV,ASYMR
cc      external ASYMEV,ASYMR
C     ..
C     .. External Subroutines ..
cc      external BRCKET,COEFF,EXTRAP,GETEF,GETRN,PQRINT,ZZERO
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN,SQRT
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETAT,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      NSAVE = 200
      JS = 5*NXINIT
      TENU = TEN*U
      BIG = ONE/U
      AVGRHO = BIG
      LPRINT = MIN(IPRINT,3)
      MAXNEV = 0
      if (REG(1)) then
         DX = SQRT(U)*MAX(ONE,ABS(A))
         call COEFF(A+DX,PX,QX,RX)
         if (FLAG.lt.0) return
         RPATA = RX*PX
         ALPHA1 = ONE/RX
         call COEFF(A+TWO*DX,PX,QX,RX)
         ALPHA1 = ALPHA1* (RPATA-PX*RX)/ (DX*FOUR)
         RPATA = SQRT(RPATA)
         ALPHA2 = SQRT(RPATA)
         ALPHA1 = (A1+A2*ALPHA1)/ALPHA2
         ALPHA2 = ALPHA2*A2
         NCOEFF = NCOEFF + 2
      else
         ALPHA1 = ZERO
         ALPHA2 = ONE
         RPATA = ONE
      end if

      if (REG(2)) then
         DX = SQRT(U)*MAX(ONE,ABS(B))
         call COEFF(B-DX,PX,QX,RX)
         if (FLAG.lt.0) return
         RPATB = RX*PX
         BETA1 = ONE/RX
         call COEFF(B-TWO*DX,PX,QX,RX)
         BETA1 = BETA1* (RPATB-PX*RX)/ (DX*FOUR)
         RPATB = SQRT(RPATB)
         BETA2 = SQRT(RPATB)
         BETA1 = (B1-B2*BETA1)/BETA2
         BETA2 = BETA2*B2
         NCOEFF = NCOEFF + 2
      else
         BETA1 = ZERO
         BETA2 = ONE
      end if

      ALPHA = A1P*A2 - A1*A2P
      if (CSPEC(2)) then
         if (ALPHA.eq.ZERO) then
            SCALE = SQRT(A1**2+A2**2)
         else
            SCALE = SQRT(ALPHA)
         end if

      end if

      if (CSPEC(1)) SCALE = SQRT(B1**2+B2**2)
      TOLMAX = MAX(TOL(1),TOL(2))
      TOL1 = MIN(TOL(1),TOLMIN)
      TOL2 = MIN(TOL(2),TOLMIN)
      ABSERR = TOL1
      RELERR = TOL2
      KLVL = 1
      DELTA = HALF
C
C     Begin the Main loop over LEVEL.
C
      do 120 LEVEL = 0,MAXLVL
         if (IPRINT.ge.2) then
            write (*,FMT=10) LEVEL,ITER
            write (21,FMT=10) LEVEL,ITER

   10       format (' Level',i3,' of iteration',i3)

         end if

         NRHO = 1
         KRHO = 0
         PSRHO = ZERO
         RHOSUM = ZERO
         ROLD = ZERO
         IASYMP = 0
         EVSAVE = -BIG
         NEV = 0
         JUMP = .FALSE.
C
C        Compute integrals needed in asymptotic formulas.
C
         QINT = ZERO
         RPINT = ZERO
         do 20 I = 2,NXINIT
            XLEFT = STORE(I-1)
            H = (STORE(I)-XLEFT)/KLVL
            HALFH = HALF*H
            do 15 J = 1,KLVL
               Z = XLEFT + HALFH
               XLEFT = XLEFT + H
               call PQRINT(Z,SQRTRP,QLNF)
               if (FLAG.lt.0) return
               QINT = QINT + H*QLNF
               RPINT = RPINT + H*SQRTRP
   15       continue
   20    continue
         if (QINT.gt.ONE/U) QINT = ZERO
         ZABS = MAX(MIN(ABSERR/100.0,RELERR)/10.0,TENU)
         ZREL = RELERR/10.0
         DELTA = MAX(DELTA/SIX,TOL1+TOL2)
C
C        Begin the secondary loop over NEV.
C
   25    if (IASYMP.ge.2) then
            AEV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
            EVHAT = AEV
            if (IASYMP.le.3) go to 35
            RSUBN = ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
            go to 45

         end if

         if (HMIN/KLVL.le.TENU) FLAG = -8
         if (FLAG.lt.0) return
         if ((LEVEL.gt.0) .and. (NEV.lt.MIN(MAXNEV,NSAVE))) then
            EV = MAX(HALF*TOL1,DELTA,HALF*TOL2*ABS(EVOLD(NEV+1)))
            EVLOW = EVOLD(NEV+1) - EV
            EVHIGH = EVOLD(NEV+1) + EV
         else
            if (LEVEL.eq.0) then
               EV = ZERO
            else
               EV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
            end if

            ETA = MAX(HALF*TOL1,DELTA,HALF*TOL2*ABS(EV))
            EVLOW = EV - ETA
            EVHIGH = EV + ETA
         end if

         call BRCKET(NEV,EVLOW,EVHIGH,FLOW,FHIGH,ZABS,ZREL,STORE)
         if ((LEVEL.eq.0) .and. (NEV.eq.0)) DELTA = EVHIGH - EVLOW
         if (FLAG.lt.0) return
         if (ABS(EVHIGH-EVLOW).gt.MAX(ZABS,
     +       ZREL*ABS(EVHIGH))) call ZZERO(EVLOW,EVHIGH,FLOW,FHIGH,ZABS,
     +       ZREL,J,STORE)
         EVHAT = MIN(EVLOW,EVHIGH)
         call GETEF(EVHAT,EFNORM,LPRINT,STORE,EFIN)
         if (FLAG.lt.0) return
         if (CSPEC(2)) then
            if (REG(1) .or. (PNU(1).eq.ZERO) .or.
     +          (PNU(1).eq.ONE-EP(1))) then
               UX = ABS(STORE(NXINIT+1))
               PDU = ABS(STORE(2*NXINIT+1))
               if (A2-A2P*EVHAT.ne.ZERO) then
                  DENOM = SCALE*UX/ABS(A2-A2P*EVHAT)
               else
                  DENOM = SCALE*PDU/ABS(A1-A1P*EVHAT)
               end if

            else
               H = STORE(2) - STORE(1)
               UX = ABS(STORE(NXINIT+2))
               PDU = ABS(STORE(2*NXINIT+2))
               if (UX.ge.PDU) then
                  DENOM = UX/H**PNU(1)
               else
                  DENOM = PDU/ (CP(1)*ABS(PNU(1))*
     +                    H** (EP(1)+PNU(1)-ONE))
               end if

            end if

         else
            if (REG(2) .or. (PNU(2).eq.ZERO) .or.
     +          (PNU(2).eq.ONE-EP(2))) then
               UX = ABS(STORE(2*NXINIT))
               PDU = ABS(STORE(3*NXINIT))
               if (B2.ne.ZERO) then
                  DENOM = SCALE*UX/ABS(B2)
               else
                  DENOM = SCALE*PDU/ABS(B1)
               end if

            else
               H = STORE(NXINIT) - STORE(NXINIT-1)
               UX = ABS(STORE(2*NXINIT-1))
               PDU = ABS(STORE(3*NXINIT-1))
               if (UX.ge.PDU) then
                  DENOM = UX/H**PNU(2)
               else
                  DENOM = PDU/ (CP(2)*ABS(PNU(2))*
     +                    H** (EP(2)+PNU(2)-ONE))
               end if

            end if

         end if

         if (BIG*DENOM.ge.ONE) then
            EFNORM = ONE/DENOM**2
         else
            EFNORM = ONE/U**2
         end if
C
C           Test for asymptotic EV.
C
         AEV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
         if (TEN*ABS(AEV-EVHAT).gt.MAX(ABSERR,RELERR*ABS(EVHAT))) then
            IASYMP = 0
         else
            if (IASYMP.lt.1) then
               IASYMP = 1
            else
               if (IPRINT.ge.2) then
                  write (*,FMT=30) NEV
                  write (21,FMT=30) NEV

   30             format (' Switchover to asymptotic eigenvalues at',
     +                   ' Nev =',i8)

               end if

               IASYMP = 2
            end if

         end if

         if (EFNORM.lt.BIG) then
            RSUBN = ONE/ (ALPHA+EFNORM)
         else
            RSUBN = ZERO
         end if
C
C           Test for asymptotic RsubN.
C
   35    ARN = ASYMR(NEV,RPINT,RPATA,RPATB,SCALE)
         if (IASYMP.ge.2) then
C
C           Eigenvalues from asymptotic formulas; produce current RsubN.
C
            call GETRN(AEV,ALPHA,CSPEC,SCALE,RSUBN,STORE)
            if (FLAG.lt.0) return
            if (ABS(ARN-RSUBN).gt.MAX(ABSERR/100.0,RELERR*ARN)) then
               IASYMP = 2
            else
               if (IASYMP.lt.3) then
                  IASYMP = 3
               else
                  if (IPRINT.ge.2) then
                     write (*,FMT=40) NEV
                     write (21,FMT=40) NEV

   40                format (' Switchover to asymptotic RsubN at ',
     +                      'Nev = ',i8)

                  end if

                  IASYMP = 4
               end if

            end if

         end if

   45    if (NEV.lt.NSAVE) EVOLD(NEV+1) = EVHAT
         if (NEV.gt.0) ETA = HALF* (EVHAT+EVSAVE)
         if (EVHAT.lt.CUTOFF) PSRHO = PSRHO + RSUBN
   50    if (T(NRHO).le.CUTOFF) then
C
C              Use step functions for Rho(t).
C
            if (T(NRHO).le.EVHAT) then
               RHO(NRHO) = RHOSUM
               NRHO = NRHO + 1
               if (NRHO.le.NUMT) then
                  go to 50

               else
                  go to 85

               end if

            end if

         else
            if (NEV.eq.0) go to 78
C
C              Use linear interpolation for Rho(t).
C
   55       if (T(NRHO).le.ETA) then
               if (JUMP) then
   60             if (T(NRHO).le.EVSAVE) then
                     RHO(NRHO) = RHOSUM - ROLD
                     NRHO = NRHO + 1
                     if (NRHO.le.NUMT) then
                        go to 60

                     else
                        JUMP = .FALSE.
                        go to 85

                     end if

                  else
                     JUMP = .FALSE.
   65                if (T(NRHO).le.ETA) then
                        RHO(NRHO) = RHOSUM
                        NRHO = NRHO + 1
                        if (NRHO.gt.NUMT) go to 85
                        go to 65

                     else
                        go to 70

                     end if

                  end if

               end if

               if ((NEV.le.1) .or. (CUTOFF.gt.EVSAVE)) then
                  UX = CUTOFF
               else
                  UX = ETAOLD
               end if

               RHO(NRHO) = RHOSUM - ROLD* (ETA-T(NRHO))/ (ETA-UX)
               NRHO = NRHO + 1
               if (NRHO.le.NUMT) then
                  go to 55

               else
                  go to 85

               end if

            end if

   70       if ((RSUBN.gt.MAX(SIX*ROLD,TEN*TOLMAX,FOUR*AVGRHO)) .and.
     +          (EVHAT.gt.CUTOFF) .and. (KRHO.gt.4)) then
C
C                 Possible eigenvalue in the continuous spectrum.
C
               LFLAG(3) = .TRUE.
               if (IPRINT.ge.1) then
                  write (*,FMT=75) EVHAT
                  write (21,FMT=75) EVHAT

   75             format (' Large jump in the step spectral density',
     +                   ' function at',d17.10)

                  write (*,FMT=76) ITER,LEVEL,RSUBN
                  write (21,FMT=76) ITER,LEVEL,RSUBN

   76             format (18x,'Iteration =',i2,', level = ',i2,
     +                   ', jump = ',d17.10)

               end if

               JUMP = .TRUE.
            end if

         end if

         if (NEV.gt.0) ETAOLD = ETA
   78    RHOSUM = RHOSUM + RSUBN
         ROLD = RSUBN
         EVSAVE = EVHAT
C
C           Output requested information.
C
         if (IPRINT.ge.3) then
            if ((NEV.le.25) .or. ((IASYMP.le.1).and.
     +          (IPRINT.ge.5))) then
               write (21,FMT=80) NEV,EVHAT,RSUBN
               write (*,FMT=80) NEV,EVHAT,RSUBN

   80          format (' Nev =',i7,', EvHat =',d15.6,', RHat =',d15.6)

            else
               if ((NEV.lt.100) .and. (10* (NEV/10).eq.NEV)) then
                  write (21,FMT=80) NEV,EVHAT,RSUBN
                  write (*,FMT=80) NEV,EVHAT,RSUBN
               else
                  if ((NEV.lt.1000) .and. (100* (NEV/100).eq.NEV)) then
                     write (21,FMT=80) NEV,EVHAT,RSUBN
                     write (*,FMT=80) NEV,EVHAT,RSUBN
                  else
                     if (1000* (NEV/1000).eq.NEV) then
                        write (21,FMT=80) NEV,EVHAT,RSUBN
                        write (*,FMT=80) NEV,EVHAT,RSUBN
                     end if

                  end if

               end if

            end if

         end if

         NEV = NEV + 1
         if (EVHAT.ge.CUTOFF) then
            if (KRHO.eq.5) then
               RSAVE(1) = RSAVE(2)
               RSAVE(2) = RSAVE(3)
               RSAVE(3) = RSAVE(4)
               RSAVE(4) = RSAVE(5)
            else
               KRHO = KRHO + 1
            end if

            RSAVE(KRHO) = RSUBN
            AVGRHO = ZERO
            do 84 I = 1,KRHO
               AVGRHO = AVGRHO + RSAVE(I)
   84       continue
            if (KRHO.gt.0) AVGRHO = AVGRHO/KRHO
         end if

         go to 25
C
C        End of Nev loop ------------------
C
   85    MAXNEV = NEV
         NUMEV = NUMEV + MAXNEV
         if (IPRINT.ge.3) then
            write (*,FMT=90) MAXNEV
            write (21,FMT=90) MAXNEV

   90       format (' MaxNev = ',i8)

         end if

         if (IPRINT.ge.4) then
            write (*,FMT=95)
            write (21,FMT=95)

   95       format (9x,'t',18x,'RhoHat(t)')

            do 105 J = 1,NUMT
               write (*,FMT=100) T(J),RHO(J)
               write (21,FMT=100) T(J),RHO(J)

  100          format (f12.4,e31.15)

  105       continue
         end if
C
C        Extrapolate interpolated approximations.
C
         DONE = .TRUE.
         do 110 J = 1,NUMT
            XTOL = MAX(TOL1,ABS(RHO(J))*TOL2)
            call EXTRAP(RHO(J),XTOL,LEVEL+1,NEXTRP,.TRUE.,.FALSE.,1,
     +                  STORE((MAXLVL+1)*J+JS),IPRINT,ERROR,RDONE)
            if (RHO(J).lt.ZERO) RHO(J) = ZERO
            if (J.gt.1) RHO(J) = MAX(RHO(J),RHO(J-1))
            if (ERROR.le.HALF*XTOL) RDONE = .TRUE.
            DONE = RDONE .and. DONE
  110    continue
         if (DONE) return
         ABSERR = MAX(HALF*ABSERR,TENU)
         RELERR = MAX(HALF*RELERR,TENU)
  120 continue
      FLAG = -3
      return

      end subroutine DENSEF
****************************** End of Part 2 ***************************
******************************* Start of Part 3 ************************
C-----------------------------------------------------------------------
      subroutine DSCRIP(LC,LPLC,TYPE,REG,CSPEC,CEV,CUTOFF,LASTEV,A1,A1P,
     +                  A2,A2P,B1,B2)
C
C     Output (if requested) a description of the spectrum.
C
C     .. Scalar Arguments ..
      double precision A1,A1P,A2,A2P,B1,B2,CUTOFF
      integer LASTEV
C     ..
C     .. Array Arguments ..
      double precision CEV(*)
      logical CSPEC(*),LC(*),LPLC(*),REG(*),TYPE(4,*)
C     ..
      write (21,FMT=*)
      if (TYPE(3,1) .and. TYPE(3,2)) then
C
C     Category 1
C
         write (21,FMT=123) 1
         write (21,FMT=100)
         write (21,FMT=121)
         write (21,FMT=102)
         write (21,FMT=110)
         if (REG(1)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            write (21,FMT=114)
            if (LC(1)) then
               if (LPLC(1)) write (21,FMT=117)
            else
               if (LPLC(1)) write (21,FMT=118)
            end if

            write (21,FMT=119) A1,A1P,A2,A2P
         end if

         write (21,FMT=111)
         if (REG(2)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            write (21,FMT=114)
            if (LC(2)) then
               if (LPLC(2)) write (21,FMT=117)
            else
               if (LPLC(2)) write (21,FMT=118)
            end if

            write (21,FMT=119) B1,B2
         end if

      end if
C
      if ((TYPE(3,1).and.CSPEC(2)) .or. (TYPE(3,2).and.CSPEC(1))) then
C
C        Category 2
C
         write (21,FMT=123) 2
         write (21,FMT=100)
         write (21,FMT=103) CUTOFF
         if (LASTEV.eq.-5) then
            write (21,FMT=105)
            write (21,FMT=120)
         else
            if (LASTEV.eq.0) then
               write (21,FMT=107)
            else
               if (LASTEV.eq.1) then
                  write (21,FMT=108)
               else
                  write (21,FMT=109) LASTEV
               end if

            end if

         end if

         write (21,FMT=110)
         if (REG(1)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (CSPEC(1)) then
               write (21,FMT=116) CEV(1)
            else
               write (21,FMT=114)
               write (21,FMT=119) A1,A1P,A2,A2P
            end if

            if (LC(1)) then
               if (LPLC(1)) write (21,FMT=117)
            else
               if (LPLC(1)) write (21,FMT=118)
            end if

         end if

         write (21,FMT=111)
         if (REG(2)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (CSPEC(2)) then
               write (21,FMT=116) CEV(2)
            else
               write (21,FMT=114)
               write (21,FMT=119) B1,B2
            end if

            if (LC(2)) then
               if (LPLC(2)) write (21,FMT=117)
            else
               if (LPLC(2)) write (21,FMT=118)
            end if

         end if

      end if
C
      if ((TYPE(3,1).and. (TYPE(4,2).and.LC(2).and.
     +    (.not.CSPEC(2)))) .or. (TYPE(3,2).and. (TYPE(4,
     +    1).and.LC(1).and. (.not.CSPEC(1))))) then
C
C        Category 3
C
         write (21,FMT=123) 3
         write (21,FMT=100)
         write (21,FMT=106)
         write (21,FMT=102)
         write (21,FMT=110)
         if (REG(1)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (TYPE(4,1)) then
               write (21,FMT=115)
            else
               write (21,FMT=114)
               write (21,FMT=119) A1,A1P,A2,A2P
            end if

            if (LC(1)) then
               if (LPLC(1)) write (21,FMT=117)
            else
               if (LPLC(1)) write (21,FMT=118)
            end if

         end if

         write (21,FMT=111)
         if (REG(2)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (TYPE(4,2)) then
               write (21,FMT=115)
            else
               write (21,FMT=114)
               write (21,FMT=119) B1,B2
            end if

            if (LC(2)) then
               if (LPLC(2)) write (21,FMT=117)
            else
               if (LPLC(2)) write (21,FMT=118)
            end if

         end if

      end if
C
      if ((TYPE(3,1).and. ((.not.LC(2)).and.TYPE(4,
     +    2).and. (.not.CSPEC(2)))) .or.
     +    (TYPE(3,2).and. ((.not.LC(1)).and.TYPE(4,
     +    1).and. (.not.CSPEC(1))))) then
C
C        Category 4
C
         write (21,FMT=123) 4
         write (21,FMT=100)
         write (21,FMT=104)
         write (21,FMT=110)
         if (REG(1)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (TYPE(4,1)) then
               write (21,FMT=115)
            else
               write (21,FMT=114)
               write (21,FMT=119) A1,A1P,A2,A2P
            end if

            if (LC(1)) then
               if (LPLC(1)) write (21,FMT=117)
            else
               if (LPLC(1)) write (21,FMT=118)
            end if

         end if

         write (21,FMT=111)
         if (REG(2)) then
            write (21,FMT=112)
         else
            write (21,FMT=113)
            if (TYPE(4,2)) then
               write (21,FMT=115)
            else
               write (21,FMT=114)
               write (21,FMT=119) B1,B2
            end if

            if (LC(2)) then
               if (LPLC(2)) write (21,FMT=117)
            else
               if (LPLC(2)) write (21,FMT=118)
            end if

         end if

      end if
C
      if ((LC(1).and.TYPE(4,1)) .and. (LC(2).and.TYPE(4,2))) then
C
C        Category 5
C
         write (21,FMT=123) 5
         write (21,FMT=100)
         write (21,FMT=106)
         write (21,FMT=102)
         write (21,FMT=110)
         write (21,FMT=113)
         write (21,FMT=115)
         write (21,FMT=117)
         write (21,FMT=111)
         write (21,FMT=113)
         write (21,FMT=115)
         write (21,FMT=117)
      end if
C
      if ((TYPE(4,1).and. (.not.LC(1)).and. (.not.CSPEC(1))) .and.
     +    (TYPE(4,2).and. (.not.LC(2)).and. (.not.CSPEC(2)))) then
C
C        Category 6
C
         write (21,FMT=123) 6
         write (21,FMT=122)
         write (21,FMT=110)
         write (21,FMT=113)
         write (21,FMT=115)
         write (21,FMT=118)
         write (21,FMT=111)
         write (21,FMT=113)
         write (21,FMT=115)
         write (21,FMT=118)
      end if
C
      if ((TYPE(4,1).and.TYPE(4,2).and..not. (CSPEC(1).or.
     +    CSPEC(2))) .and. ((LC(1).and. (.not.LC(2))).or.
     +    ((.not.LC(1)).and.LC(2)))) then
C
C        Category 7
C
         write (21,FMT=123) 7
         write (21,FMT=104)
         write (21,FMT=110)
         write (21,FMT=113)
         write (21,FMT=115)
         if (LC(1)) then
            if (LPLC(1)) write (21,FMT=117)
         else
            if (LPLC(1)) write (21,FMT=118)
         end if

         write (21,FMT=111)
         write (21,FMT=113)
         write (21,FMT=115)
         if (LC(2)) then
            if (LPLC(2)) write (21,FMT=117)
         else
            if (LPLC(2)) write (21,FMT=118)
         end if

      end if
C
      if ((LC(1).and.TYPE(4,1).and.CSPEC(2)) .or.
     +    (LC(2).and.TYPE(4,2).and.CSPEC(1))) then
C
C        Category 8
C
         write (21,FMT=123) 8
         write (21,FMT=100)
         write (21,FMT=103) CUTOFF
         write (21,FMT=110)
         write (21,FMT=113)
         if (CSPEC(1)) then
            write (21,FMT=116) CEV(1)
         else
            write (21,FMT=115)
         end if

         if (LC(1)) then
            if (LPLC(1)) write (21,FMT=117)
         else
            if (LPLC(1)) write (21,FMT=118)
         end if

         write (21,FMT=111)
         write (21,FMT=113)
         if (CSPEC(2)) then
            write (21,FMT=116) CEV(2)
         else
            write (21,FMT=115)
         end if

         if (LC(2)) then
            if (LPLC(2)) write (21,FMT=117)
         else
            if (LPLC(2)) write (21,FMT=118)
         end if

      end if
C
      if ((((.not.LC(1)).and.TYPE(4,1).and. (.not.CSPEC(1))).and.
     +    CSPEC(2)) .or. (((.not.LC(2)).and.TYPE(4,
     +    2).and. (.not.CSPEC(2))).and.CSPEC(1))) then
C
C        Category 9
C
         write (21,FMT=123) 9
         write (21,FMT=101)
         write (21,FMT=110)
         write (21,FMT=113)
         if (CSPEC(1)) then
            write (21,FMT=116) CEV(1)
         else
            write (21,FMT=115)
         end if

         if (LC(1)) then
            if (LPLC(1)) write (21,FMT=117)
         else
            if (LPLC(1)) write (21,FMT=118)
         end if

         write (21,FMT=111)
         write (21,FMT=113)
         if (CSPEC(2)) then
            write (21,FMT=116) CEV(2)
         else
            write (21,FMT=115)
         end if

         if (LC(2)) then
            if (LPLC(2)) write (21,FMT=117)
         else
            if (LPLC(2)) write (21,FMT=118)
         end if

      end if
C
      if (CSPEC(1) .and. CSPEC(2)) then
C
C        Category 10
C
         write (21,FMT=123) 10
         write (21,FMT=101)
         write (21,FMT=103) CUTOFF
         if (LASTEV.eq.-5) then
            write (21,FMT=120)
         else
            if (LASTEV.eq.0) then
               write (21,FMT=107)
            else
               if (LASTEV.eq.1) then
                  write (21,FMT=108)
               else
                  write (21,FMT=109) LASTEV
               end if

            end if

         end if

         write (21,FMT=110)
         write (21,FMT=113)
         write (21,FMT=116) CEV(1)
         if (LC(1)) then
            if (LPLC(1)) write (21,FMT=117)
         else
            if (LPLC(1)) write (21,FMT=118)
         end if

         write (21,FMT=111)
         write (21,FMT=113)
         write (21,FMT=116) CEV(2)
         if (LC(2)) then
            if (LPLC(2)) write (21,FMT=117)
         else
            if (LPLC(2)) write (21,FMT=118)
         end if

      end if

      write (21,FMT=*)
      return

  100 format (' This problem has simple spectrum.')
  101 format (' This problem may have non-simple spectrum.')
  102 format (' There is no continuous spectrum.')
  103 format (' There is continuous spectrum in [Ev, infinity) where',
     +       ' Ev =',g15.6)
  104 format (' There is continuous spectrum consisting of the entire ',
     +       'real line.')
  105 format (' The set of eigenvalues is bounded below.')
  106 format (' There are infinitely many negative and infinitely many',
     +       ' positive',/'    eigenvalues (unbounded in either',
     +       ' direction).')
  107 format (' There appear to be no eigenvalues below the start of',
     +       ' the continuous spectrum.')
  108 format (' There appears to be 1 eigenvalue below the start of the'
     +       ,' continuous spectrum.')
  109 format (' There appear to be ',i12,' eigenvalues below the start',
     +       /'    of the continuous spectrum.')
  110 format (' At endpoint A')
  111 format (' At endpoint B')
  112 format ('    the problem is regular;')
  113 format ('    the problem is singular;')
  114 format ('    it is nonoscillatory for all Ev.')
  115 format ('    it is oscillatory for all Ev.')
  116 format ('    it is nonoscillatory for Ev <',g15.6,
     +       ' and oscillatory otherwise.')
  117 format ('    It is limit circle.')
  118 format ('    It is limit point.')
  119 format ('    The constants for the Friedrichs boundary conditions'
     +       ,' are',/4e18.8)
  120 format (' There appear to be infinitely many eigenvalues below',
     +       ' the start',/'    of the continuous spectrum.')
  121 format (' There are infinitely many eigenvalues, bounded below.')
  122 format (' The nature of the spectrum is unknown; there is likely',
     +       ' to be ',/' continuous spectrum.')
  123 format (' The spectral category is',i3,'.')

      end subroutine DSCRIP
C-----------------------------------------------------------------------
      subroutine EXTRAP(V,TOL,IROW,MAXCOL,FULL,TIGHT,MODE,VSAVE,IPRINT,
     +                  ERROR,DONE)
C
C     Use Richardson's h**2 extrapolation (based on doubling) when
C     suitable, otherwise use Wynn's acceleration scheme.
C
C     Input:
C        V       = real value at current level.
C        TOL     = real tolerance.
C        IROW    = integer giving current row index (1 .le. IROW).
C        MAXCOL  = integer giving maximum number of columns in table.
C        FULL    = logical, True iff entire table is to be computed.
C        TIGHT   = logical, True for conservative convergence tests.
C        MODE    = integer, value is
C                  0   both Richardson and Wynn algorithms can be used;
C                  1   only Richardson is used;
C                  2   only Wynn is used.
C        IPRINT  = integer controlling amount of printing.
C     Output:
C        V       = real, best output estimate.
C        VSAVE   = real vector, holds previous level values.
C        DONE    = logical, True iff Error is sufficiently small.
C
C     If FULL is True, then the entire acceleration array is produced
C     (through row IROW); if False, then only the next row is appended.
C     Hence, the choice of FULL = True requires more work, but it may
C     save some global storage.  The vector VSAVE contains the V values
C     for levels 0 through max(IROW-1,MAXCOL).
C
C
C     The local arrays RATIO(*), R(*,*), and W(*,*) must be declared to
C     have at least as many rows as the value of MAXLVL initialized in
C     routine START.
C
C
C     .. Parameters ..
      double precision ZERO,TENTH,ONE,TWO
      parameter (ZERO=0.0,TENTH=0.1D0,ONE=1.0,TWO=2.0)
C     ..
C     .. Scalar Arguments ..
      double precision ERROR,TOL,V
      integer IPRINT,IROW,MAXCOL,MODE
      logical DONE,FULL,TIGHT
C     ..
C     .. Array Arguments ..
      double precision VSAVE(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
C     ..
C     .. Local Scalars ..
      double precision DIFF,EPS,ETEMP,RHIGH,RLOW,RTOL,T,TOL1,TOL2,VTEMP
      integer I,IMIN,J,MAXJ,NCOL
C     ..
C     .. Local Arrays ..
      double precision R(12,8),RATIO(12),W(40,11)
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,LOG,MAX,MIN,MOD
C     ..
C     .. Common blocks ..
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Save statement ..
      save R,W
C     ..
      EPS = 0.2
      DONE = .FALSE.
      ETEMP = ONE/U
      VTEMP = V
      VSAVE(IROW) = V
      TOL1 = TOL
      TOL2 = TOL/3.0
      if (MODE.eq.2) then
         MAXJ = MAXCOL
         go to 40

      else
         MAXJ = 11
      end if

      RLOW = MAX(3.0,4.2-0.2* (IROW-1))
      RHIGH = MIN(5.0,3.8+0.2* (IROW-1))
      RTOL = U
C
C     Analyze the rate of convergence to determine NCOL and tolerances.
C
      NCOL = 2
      do 10 I = 1,IROW
         R(I,1) = VSAVE(I)
         RATIO(I) = ONE/U
         if (I.ge.3) then
            T = R(I,1) - R(I-1,1)
            if (T.ne.ZERO) then
               RATIO(I) = (R(I-1,1)-R(I-2,1))/T
            else
               V = R(I,1)
               DONE = .TRUE.
               return

            end if

            if (((RATIO(I).ge.RLOW).and. (RATIO(I).le.RHIGH)) .or.
     +          (.not.TIGHT)) then
               RTOL = TOL1
               NCOL = MIN(MAXCOL,NCOL+1)
            else
               RTOL = TOL2
               if (RATIO(I).lt.ZERO) RTOL = TENTH*TOL2
               if (RATIO(I).lt.TWO) NCOL = 2
            end if

         end if

   10 continue
      if (FULL) then
         IMIN = 2
      else
         IMIN = IROW
      end if
C
C     Use Richardson's h^2 extrapolation.  The number of columns used
C     is a function of the amount of data (IROW), the requested order
C     (MAXCOL), the observed rate of convergence (NCOL), and the amount
C     of storage allocated to R(*,*).
C
      do 30 I = IMIN,IROW
         do 20 J = 2,MIN(I,NCOL,8)
            DIFF = (R(I,J-1)-R(I-1,J-1))/ (4** (J-1)-1)
            R(I,J) = R(I,J-1) + DIFF
            if ((.not.FULL) .or. (I.eq.IROW)) then
               T = ABS(DIFF)
               if (T.le.ETEMP) then
                  ETEMP = T
                  VTEMP = R(I,J)
                  if (T.le.RTOL) then
                     DONE = .TRUE.
                     V = VTEMP
                     ERROR = ETEMP
                     return

                  end if

               end if

            end if

   20    continue
   30 continue
      V = VTEMP
      ERROR = ETEMP
      if (IROW.lt.4) return
C
C     Test for rate of convergence other than second order.
C
      if (ABS(RATIO(IROW)-RATIO(IROW-1))+
     +    ABS(RATIO(IROW-1)-RATIO(IROW-2)).gt.EPS) return
      if ((3.5.lt.RATIO(IROW)) .and. (RATIO(IROW).lt.4.5)) return
      if (RATIO(IROW).lt.ONE) return
      if (MODE.ne.0) return
C
C     Use Wynn's algorithm.
C
   40 if (IPRINT.ge.4) then
         DIFF = LOG(ABS(RATIO(IROW)))/LOG(TWO)
         write (21,FMT=50) DIFF

   50    format (' In EXTRAP: using Wynn`s acceleration; rate = ',f8.5)

      end if

      W(1,1) = VSAVE(1)
      ETEMP = ONE/U
      do 60 I = 2,IROW
         W(I,1) = VSAVE(I)
         DIFF = W(I,1) - W(I-1,1)
         if ((I.eq.IROW) .or. (MODE.eq.2)) then
            T = ABS(DIFF)
            if (T.le.ETEMP) then
               ETEMP = T
               VTEMP = W(I,1)
               if (T.le.TOL2) then
                  V = VTEMP
                  ERROR = ETEMP
                  DONE = .TRUE.
                  return

               end if

            end if

         end if

         if (DIFF.ne.ZERO) then
            W(I,2) = ONE/DIFF
         else
            V = W(I,1)
            ERROR = ZERO
            DONE = .TRUE.
            return

         end if

   60 continue
      do 80 J = 3,MIN(IROW,MAXJ)
         do 70 I = J,IROW
            DIFF = W(I,J-1) - W(I-1,J-1)
            if (DIFF.ne.ZERO) then
               DIFF = ONE/DIFF
            else
               if (MOD(J,2).eq.0) then
                  V = W(I,J-1)
                  ERROR = ZERO
                  DONE = .TRUE.
                  return

               else
                  DIFF = ONE/U**2
               end if

            end if

            W(I,J) = W(I-1,J-2) + DIFF
            if ((MOD(J,2).eq.1) .and. ((I.eq.IROW).or.
     +          (MODE.eq.2))) then
               T = ABS(DIFF)
               if (T.le.ETEMP) then
                  ETEMP = T
                  VTEMP = W(I,J)
                  if (T.le.TOL2) then
                     V = VTEMP
                     ERROR = ETEMP
                     DONE = .TRUE.
                     return

                  end if

               end if

            end if

   70    continue
   80 continue
      V = VTEMP
      ERROR = ETEMP
      return

      end subroutine EXTRAP
C-----------------------------------------------------------------------
      subroutine GETEF(EV,EFNORM,IPRINT,X,EFIN)
C
C     Compute an eigenfunction for one fixed mesh.
C
C
C     .. Parameters ..
      double precision ZERO,C10M4,HALF,ONE,TWO,THREE,FIVE,C15,C21,TINY
      parameter (ZERO=0.0,C10M4=1.D-4,HALF=0.5D0,ONE=1.0,TWO=2.0,
     +          THREE=3.0,FIVE=5.0,C15=15.0,C21=21.0,TINY=1.D-38)
C     ..
C     .. Scalar Arguments ..
      double precision EFNORM,EV
      integer IPRINT
C     ..
C     .. Array Arguments ..
      double precision X(*)
      logical EFIN(2,2)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision CHI,DPSI,DV,DW,FNORM,FSCALE,FSUM,H,HALFH,HOM,OM,
     +                 OSCALE,PDV,PDW,PN,PROD,PSI,RATIO,RN,RNORM,RSCALE,
     +                 RSUM,SCALE,T,TAU,TAUHH,TAUMAX,V,VNEW,W,WNEW,
     +                 XLEFT,XRIGHT,Z,
     +                 FTERM,RTERM
      integer I,J,JDU,JLAST,JS,JU,JX,KLVL,MODE
      logical ALLOK,SYMM
C     ..
C     .. Local Arrays ..
      integer MIDDLE(2)
C     ..
C     .. External Subroutines ..
cc      external STEP
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,EXP,LOG,MAX,SIGN
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      KLVL = 2**LEVEL
C
C     For this EV and mesh, calculate FSCALE, RSCALE, and MIDDLE(*).
C         FSCALE = sum of logs of scale factors 1 through m
C         RSCALE = sum of logs of scale factors m+1 through N
C         MIDDLE(1), MIDDLE(2) describe the coordinates of the matching
C            point M for the shooting; in particular,
C            M = X(MIDDLE(1)-1) + MIDDLE(2)*H(MIDDLE(1)-1)/2**LEVEL
C     The matching point is chosen to be roughly (a+b)/2 if either
C     a = -b and Tau(x) > 0 near 0, or if Tau(x) > 0 for all x;
C     otherwise, it is chosen to roughly maximize Tau(x).
C
      FSCALE = ZERO
      RSCALE = ZERO
      TAUMAX = -ONE/U
      ALLOK = .TRUE.
      SYMM = .FALSE.
      if (A.eq.-B) SYMM = .TRUE.
      EFIN(1,1) = .TRUE.
      EFIN(2,1) = .TRUE.
      EFIN(1,2) = .TRUE.
      EFIN(2,2) = .TRUE.
      JX = (NXINIT+1)/2
      do 20 I = 2,NXINIT
         MODE = 0
         XLEFT = X(I-1)
         H = X(I) - XLEFT
         if (I.eq.2) then
            if (KCLASS(1).ge.9) MODE = 1
            if (.not.AFIN) then
               MODE = 3
               XLEFT = ZERO
               H = -ONE/X(2)
            end if

         end if

         if (I.eq.NXINIT) then
            if (KCLASS(2).ge.9) MODE = 2
            if (.not.BFIN) then
               MODE = 4
               H = ONE/X(I-1)
               XLEFT = -H
            end if

         end if

         H = H/KLVL
         HALFH = HALF*H
         do 10 J = 1,KLVL
            Z = XLEFT + HALFH
            XLEFT = XLEFT + H
            call STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            if (TAU.lt.ZERO) ALLOK = .FALSE.
            if (TAU.gt.TAUMAX) then
               TAUMAX = TAU
               MIDDLE(1) = I
               MIDDLE(2) = J
               FSCALE = FSCALE + RSCALE + SCALE
               OSCALE = SCALE
               RSCALE = ZERO
            else
               RSCALE = RSCALE + SCALE
            end if

   10    continue
         if ((I.eq.JX) .and. (TAU.lt.ZERO)) SYMM = .FALSE.
   20 continue
      if (ALLOK .or. SYMM) then
         MIDDLE(1) = JX
         MIDDLE(2) = MAX(KLVL-1,1)
      end if

      if ((.not.AFIN) .or. (.not.BFIN)) then
C
C        Don't split near infinity!
C
         if ((.not.AFIN) .and. (MIDDLE(1).eq.2)) then
            if (REG(2)) then
               MIDDLE(1) = NXINIT
            else
               MIDDLE(1) = JX
            end if

            MIDDLE(2) = MAX(KLVL-1,1)
         end if

         if ((.not.BFIN) .and. (MIDDLE(1).eq.NXINIT)) then
            if (REG(1)) then
               MIDDLE(1) = 2
            else
               MIDDLE(1) = JX
            end if

            MIDDLE(2) = 1
         end if

      end if

      if ((LEVEL.gt.1) .and. (MIDDLE(2).eq.KLVL)) then
         MIDDLE(2) = KLVL - 1
         FSCALE = FSCALE - OSCALE
         RSCALE = RSCALE + OSCALE
      end if

      if (IPRINT.ge.4) write (21,FMT=21) MIDDLE(1),MIDDLE(2)

   21 format (' Coordinates of matching point =',2i6)

      JU = NXINIT
      JDU = 2*NXINIT
      JS = 3*NXINIT
C
C     Shoot from x=A to the middle.
C
      V = A2 - A2P*EV
      PDV = A1 - A1P*EV
      FNORM = ZERO
      if (KCLASS(1).lt.9) then
         SCALE = MAX(ABS(V),ABS(PDV))
         V = V/SCALE
         PDV = PDV/SCALE
         MODE = 0
      else
         V = ONE
         PDV = D(4,1)
         MODE = 1
      end if

      if (.not.AFIN) MODE = 3
      FSUM = -FSCALE
      if ((IPRINT.ge.5) .and. (MODE.eq.0)) write (21,FMT=35) X(1),V,PDV,
     +    FSUM
      do 40 I = 2,MIDDLE(1)
         X(JU+I-1) = V
         X(JDU+I-1) = PDV
         X(JS+I-1) = FSUM
         if (MODE.eq.0) then
            XLEFT = X(I-1)
            H = X(I) - XLEFT
         else
            XLEFT = ZERO
            if (MODE.eq.1) then
               H = ((X(2)-X(1))/D(1,1))**ETA(1,1)
            else
               H = -ONE/X(2)
            end if

         end if

         H = H/KLVL
         HALFH = HALF*H
         JLAST = KLVL
         if (I.eq.MIDDLE(1)) JLAST = MIDDLE(2)
         do 30 J = 1,JLAST
            Z = XLEFT + HALFH
            XLEFT = XLEFT + H
            call STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            DV = PDV/PN
            if (ABS(TAU)*H*H.ge.C10M4) then
               if (TWO* (FSUM+SCALE).gt.-UNDER) then
                  FNORM = FNORM + RN* (PSI*DPSI* (V*V-DV*DV/TAU)/
     +                    TWO+V*PSI*DV*PSI)*EXP(TWO* (FSUM+SCALE))
                  if (TWO*FSUM.gt.-UNDER) FNORM = FNORM +
     +                RN*EXP(TWO*FSUM)*H* (V*V+DV*DV/TAU)/TWO
               end if

            else
               if (TWO*FSUM.gt.-UNDER) then
                  TAUHH = TAU*H*H
                  FNORM = FNORM + RN*H*EXP(TWO*FSUM)*
     +                    (V*V* (ONE+TAUHH* (TAUHH/FIVE-ONE)/
     +                    THREE)+H*V*DV* (ONE+TAUHH* (TWO*TAUHH/
     +                    C15-ONE)/THREE)+ (H*DV)**2*
     +                    (ONE+TAUHH* (TWO*TAUHH/C21-ONE)/FIVE)/THREE)
               end if

            end if

            FSUM = FSUM + SCALE
            VNEW = DPSI*V + PSI*DV
            PDV = -PN*TAU*PSI*V + DPSI*PDV
            V = VNEW
   30    continue
         if (MODE.eq.1) then
C
C           Convert from V(t) to u(x).
C
            if (ETA(2,1).lt.ZERO) then
               X(JU+1) = ONE/U
               EFIN(1,1) = .FALSE.
            else
               if (ETA(2,1).eq.ZERO) then
                  X(JU+1) = D(2,1)
               else
                  X(JU+1) = ZERO
               end if

            end if

            Z = ETA(1,1)*ETA(2,1) + EP(1) - ONE
            if (Z.lt.ZERO) then
               X(JDU+1) = SIGN(ONE/U,ETA(2,1))
               EFIN(2,1) = .FALSE.
            else
               if (Z.gt.ZERO) then
                  X(JDU+1) = ZERO
               else
                  X(JDU+1) = D(2,1)*CP(1)*ETA(1,1)*ETA(2,1)/
     +                       D(1,1)** (ETA(1,1)*ETA(2,1))
               end if

            end if

            H = X(2) - A
            PN = CP(1)*H** (EP(1)-ONE)
            T = (H/D(1,1))**ETA(1,1)
            CHI = D(2,1)*T**ETA(2,1)
            V = V*CHI
            PDV = PN*ETA(1,1)* (ETA(2,1)*V+CHI*PDV*T** (ONE-TWO*EMU(1)))
            if (IPRINT.ge.5) write (21,FMT=35) X(1),X(JU+1),X(JDU+1)
         end if

         MODE = 0
         if (IPRINT.ge.5) write (21,FMT=35) XLEFT,V,PDV,FSUM

   35    format (g16.6,3d15.6)

   40 continue
C
C     Shoot from x=B to the middle.
C
      RNORM = ZERO
      MODE = 0
      if (KCLASS(2).lt.9) then
         SCALE = MAX(ABS(B1),ABS(B2))
         W = -B2/SCALE
         PDW = B1/SCALE
         MODE = 0
      else
         W = -ONE
         PDW = -D(4,2)
         MODE = 2
      end if

      if (.not.BFIN) MODE = 4
      RSUM = -RSCALE
      if ((IPRINT.ge.5) .and. (MODE.eq.0)) write (21,FMT=35) X(NXINIT),
     +    W,PDW,RSUM
      do 60 I = NXINIT,MIDDLE(1),-1
         X(JU+I) = W
         X(JDU+I) = PDW
         X(JS+I) = RSUM
         if (MODE.eq.0) then
            XRIGHT = X(I)
            H = XRIGHT - X(I-1)
         else
            XRIGHT = ZERO
            if (MODE.eq.2) then
               H = ((X(NXINIT)-X(NXINIT-1))/D(1,2))**ETA(1,2)
            else
               H = ONE/X(I-1)
            end if

         end if

         H = H/KLVL
         HALFH = HALF*H
         JLAST = KLVL
         if (I.eq.MIDDLE(1)) JLAST = JLAST - MIDDLE(2)
         do 50 J = 1,JLAST
            Z = XRIGHT - HALFH
            XRIGHT = XRIGHT - H
            call STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,MODE)
            DW = PDW/PN
            if (ABS(TAU)*H*H.ge.C10M4) then
               if (TWO* (RSUM+SCALE).gt.-UNDER) then
                  RNORM = RNORM + RN* (PSI*DPSI* (W*W-DW*DW/TAU)/
     +                    TWO-W*PSI*DW*PSI)*EXP(TWO* (RSUM+SCALE))
                  if (TWO*RSUM.gt.-UNDER) RNORM = RNORM +
     +                RN*EXP(TWO*RSUM)*H* (W*W+DW*DW/TAU)/TWO
               end if

            else
               if (TWO*RSUM.gt.-UNDER) then
                  TAUHH = TAU*H*H
                  RNORM = RNORM + RN*H*EXP(TWO*RSUM)*
     +                    (W*W* (ONE+TAUHH* (TAUHH/FIVE-ONE)/
     +                    THREE)-H*W*DW* (ONE+TAUHH* (TWO*TAUHH/
     +                    C15-ONE)/THREE)+ (H*DW)**2*
     +                    (ONE+TAUHH* (TWO*TAUHH/C21-ONE)/FIVE)/THREE)
               end if

            end if

            RSUM = RSUM + SCALE
            WNEW = DPSI*W - PSI*DW
            PDW = PN*TAU*PSI*W + DPSI*PDW
            W = WNEW
   50    continue
         if (MODE.eq.2) then
C
C           Convert from V(t) to u(x).
C
            if (ETA(2,2).lt.ZERO) then
               X(JU+NXINIT) = -ONE/U
               EFIN(1,2) = .FALSE.
            else
               if (ETA(2,2).eq.ZERO) then
                  X(JU+NXINIT) = -D(2,2)
               else
                  X(JU+NXINIT) = ZERO
               end if

            end if

            Z = ETA(1,2)*ETA(2,2) + EP(2) - ONE
            if (Z.lt.ZERO) then
               X(JDU+NXINIT) = SIGN(ONE/U,ETA(2,2))
               EFIN(2,2) = .FALSE.
            else
               if (Z.gt.ZERO) then
                  X(JDU+NXINIT) = ZERO
               else
                  X(JDU+NXINIT) = D(2,2)*CP(2)*ETA(1,2)*ETA(2,2)/
     +                            D(1,2)** (ETA(1,2)*ETA(2,2))
               end if

            end if

            if (IPRINT.ge.5) write (21,FMT=35) X(NXINIT),X(JU+NXINIT),
     +          X(JDU+NXINIT)
            H = X(NXINIT) - X(NXINIT-1)
            PN = CP(2)*H** (EP(2)-ONE)
            T = (H/D(1,2))**ETA(1,2)
            CHI = D(2,2)*T**ETA(2,2)
            W = CHI*W
            PDW = PN*ETA(1,2)* (CHI*PDW*T** (ONE-TWO*EMU(2))-ETA(2,2)*W)
         end if

         if (MODE.eq.4) then
            X(JU+NXINIT) = ZERO
            X(JDU+NXINIT) = ONE
            if (IPRINT.ge.5) write (21,FMT=35) X(NXINIT),X(JU+NXINIT),
     +          X(JDU+NXINIT)
         end if

         MODE = 0
         if ((JLAST.ne.0) .and. (IPRINT.ge.5)) write (21,FMT=35) XRIGHT,
     +       W,PDW,RSUM
   60 continue
      if (ABS(W).ge.ABS(PDW)) then
         RATIO = V/W
         if (IPRINT.ge.5) write (21,FMT=61) RATIO*PDW - PDV,RATIO

   61    format ('  DuHat jump, ratio =',2d24.15)

      else
         RATIO = PDV/PDW
         if (V*W*RATIO.lt.ZERO) RATIO = -RATIO
         if (IPRINT.ge.5) write (21,FMT=62) RATIO*W - V,RATIO

   62    format ('  UHat jump, ratio =',2d24.15)

      end if
C
C     Calculate weighted 2-norm and scale approximate eigenfunction.
C
CJDP The expression for EFNORM caused overflow unnecessarily so
C    next 4 lines altered. Note EFNORM, SCALE both used later
C    but FSCALE, RSCALE are not.
C      FSCALE = EXP(-FSUM)
C      RSCALE = EXP(-RSUM)
C      EFNORM = SQRT(FNORM*FSCALE**2+RNORM* (RATIO*RSCALE)**2)
C      SCALE = LOG(EFNORM)
C    Compute the logs of the 2 terms under the square root above,
C    with precaution against log of zero.
      FTERM = HALF*LOG(MAX(TINY,FNORM)) - FSUM
      RTERM = HALF*LOG(MAX(TINY,RNORM*RATIO**2)) - RSUM
C    SCALE can never cause overflow, and EFNORM will only cause
C    overflow if it is actually too large to store:
      if (FTERM.ge.RTERM) then
        SCALE = FTERM + HALF*LOG(ONE+EXP(TWO*(RTERM-FTERM)))
      else
        SCALE = RTERM + HALF*LOG(ONE+EXP(TWO*(FTERM-RTERM)))
      end if
      EFNORM = EXP(SCALE)
CJDP end of amendment

      if (IPRINT.ge.5) write (21,FMT=65) EFNORM

   65 format ('  EFnorm =',d24.15)

      do 70 I = 1,NXINIT
         TAU = X(JS+I) - SCALE
         if (TAU.le.-UNDER) then
            X(JU+I) = ZERO
            X(JDU+I) = ZERO
         else
            PROD = EXP(TAU)
            if (I.ge.MIDDLE(1)) PROD = PROD*RATIO
            X(JU+I) = X(JU+I)*PROD
            X(JDU+I) = X(JDU+I)*PROD
         end if

   70 continue
      if (IPRINT.ge.4) then
         write (21,FMT=75)

   75    format (10x,'x',15x,'Uhat(x)',13x,'PUhat`(x)')

         do 85 I = 1,NXINIT
            write (21,FMT=80) X(I),X(JU+I),X(JDU+I)

   80       format (g16.6,2d20.8)

   85    continue
      end if

      return

      end subroutine GETEF
C----------------------------------------------------------------------
      subroutine GETRN(EV,ALPHA,CSPEC,DENOM,RSUBN,X)
C
C     Compute the RsubN value from the weighted eigenfunction 2-norm
C     when standard shooting is stable and an accurate eigenvalue is
C     available (from the asymptotic formulas).
C
C
C     .. Parameters ..
      double precision ZERO,C10M4,HALF,ONE,TWO,THREE,FIVE,C15,C21
      parameter (ZERO=0.0,C10M4=1.D-4,HALF=0.5D0,ONE=1.D0,TWO=2.0,
     +          THREE=3.0,FIVE=5.0,C15=15.0,C21=21.0)
C     ..
C     .. Scalar Arguments ..
      double precision ALPHA,DENOM,EV,RSUBN
C     ..
C     .. Array Arguments ..
      double precision X(*)
      logical CSPEC(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,UNDER,URN
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision DPHI,DPSI,DU,DUSAVE,FNORM,H,HALFH,HOM,HSAVE,OM,
     +                 PDU,PHI,PN,PSI,RN,SCALE,TAU,TAUHH,U,UNEW,USAVE,
     +                 XLEFT,Z
      integer I,J,KLVL
C     ..
C     .. External Subroutines ..
cc      external STEP
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,EXP
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,URN,UNDER
C     ..
      U = A2 - A2P*EV
      PDU = A1 - A1P*EV
      FNORM = ZERO
      KLVL = 2**LEVEL
C
C     Shoot from x=A to x=B.
C
      do 20 I = 2,NXINIT
         XLEFT = X(I-1)
         H = (X(I)-XLEFT)/KLVL
         HALFH = HALF*H
         do 10 J = 1,KLVL
            Z = XLEFT + HALFH
            XLEFT = XLEFT + H
            call STEP(Z,H,EV,PN,RN,TAU,OM,HOM,PSI,DPSI,SCALE,0)
            SCALE = EXP(SCALE)
            PHI = PSI*SCALE
            DPHI = DPSI*SCALE
            DU = PDU/PN
            if (ABS(TAU)*H*H.ge.C10M4) then
               FNORM = FNORM + RN* (PHI*DPHI* (U*U-DU*DU/TAU)/
     +                 TWO+U*PHI*DU*PHI+H* (U*U+DU*DU/TAU)/TWO)
            else
               TAUHH = TAU*H*H
               FNORM = FNORM + RN*H* (U*U* (ONE+TAUHH* (TAUHH/FIVE-ONE)/
     +                 THREE)+H*U*DU* (ONE+TAUHH* (TWO*TAUHH/C15-ONE)/
     +                 THREE)+ (H*DU)**2* (ONE+TAUHH* (TWO*TAUHH/
     +                 C21-ONE)/FIVE)/THREE)
            end if

            if ((I.eq.NXINIT) .and. (J.eq.KLVL) .and. CSPEC(1)) then
               HSAVE = H
               USAVE = ABS(U)
               DUSAVE = ABS(PDU)
            end if

            UNEW = DPHI*U + PHI*DU
            PDU = -PN*TAU*PHI*U + DPHI*PDU
            U = UNEW
            if ((I.eq.2) .and. (J.eq.1) .and. CSPEC(2)) then
               HSAVE = H
               USAVE = ABS(U)
               DUSAVE = ABS(PDU)
            end if

   10    continue
   20 continue
      if (CSPEC(2)) then
         if (REG(1) .or. (PNU(1).eq.ZERO) .or.
     +       (PNU(1).eq.ONE-EP(1))) then
            PHI = DENOM
         else
            if (USAVE.ge.DUSAVE) then
               PHI = USAVE/HSAVE**PNU(1)
            else
               PHI = DUSAVE/ (CP(1)*ABS(PNU(1))*
     +               HSAVE** (EP(1)+PNU(1)-ONE))
            end if

         end if

      else
         if (REG(2) .or. (PNU(2).eq.ZERO) .or.
     +       (PNU(2).eq.ONE-EP(2))) then
            PHI = DENOM
         else
            if (USAVE.ge.DUSAVE) then
               PHI = USAVE/HSAVE**PNU(2)
            else
               PHI = DUSAVE/ (CP(2)*ABS(PNU(2))*
     +               HSAVE** (EP(2)+PNU(2)-ONE))
            end if

         end if

      end if

      RSUBN = ONE/ (ALPHA+FNORM/PHI**2)
      return

      end subroutine GETRN
******************************** End of Part 3 *************************
******************************* Start of Part 4 ************************
C-----------------------------------------------------------------------
      subroutine MESH(JOB,NEV,X,G,H,QLNF,Z,TOL,HMIN)
C
C     If JOB = True then calculate the initial mesh; redistribute so
C     that H(*) is approximately equidistributed.  If JOB = False
C     then use the mesh input by the user.
C
C     .. Parameters ..
      double precision ZERO,TENTH,HALF,ONE,TWO,FOUR
      parameter (ZERO=0.0,TENTH=0.1D0,HALF=0.5D0,ONE=1.0,TWO=2.0,
     +          FOUR=4.0)
C     ..
C     .. Scalar Arguments ..
      double precision HMIN,TOL
      integer NEV
      logical JOB
C     ..
C     .. Array Arguments ..
      double precision G(*),H(*),QLNF(*),X(*),Z(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision DX,ENDA,ENDB,EPS,EQMAX,EQMIN,EV,GAMMA,P1,P2,P3,
     +                 Q1,Q2,Q3,QMAX,QMIN,R1,R2,R3,WEIGHT,Y,Y1,Y2,Y3
      integer I,ITS,J,JTOL,K,MAXITS,N,NADD
      logical DONE
C     ..
C     .. External Subroutines ..
cc      external COEFF
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,EXP,LOG10,MAX,MIN,SQRT
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Save statement ..
      save ENDA,ENDB
C     ..
C     .. Data statements ..
      data EPS/0.0001/
C     ..
C
      if (.not.JOB) then
         HMIN = B - A
         do 5 I = 2,NXINIT
            HMIN = MIN(HMIN,X(I)-X(I-1))
    5    continue
         return

      end if

      N = NXINIT - 1
      EV = NEV
C
C     Find an appropriate initial mesh.
C
      if (AFIN) then
         X(1) = A
         if (BFIN) then
            X(NXINIT) = B
            DX = (B-A)/N
            do 10 I = 2,N
               X(I) = X(1) + (I-1)*DX
   10       continue
         else
            if (NEV.lt.0) then
               ENDB = X(NXINIT)
            else
               X(NXINIT) = (1+4*NEV)*ENDB
            end if

            Y1 = X(1)/ (ONE+ABS(X(1)))
            DX = (X(NXINIT)/ (ONE+X(NXINIT))-Y1)/N
            do 11 I = 2,N
               Y = Y1 + (I-1)*DX
               X(I) = Y/ (ONE-ABS(Y))
   11       continue
         end if

      else
         if (NEV.lt.0) then
            ENDA = X(1)
         else
            X(1) = (1+4*NEV)*ENDA
         end if

         if (BFIN) then
            X(NXINIT) = B
            Y1 = X(NXINIT)/ (ONE+ABS(X(NXINIT)))
            DX = (Y1-X(1)/ (ONE-X(1)))/N
            do 12 I = 2,N
               Y = Y1 - (I-1)*DX
               X(NXINIT+1-I) = Y/ (ONE-ABS(Y))
   12       continue
         else
            Y1 = X(1)/ (ONE-X(1))
            if (NEV.lt.0) then
               ENDB = X(NXINIT)
            else
               X(NXINIT) = (1+4*NEV)*ENDB
            end if

            Y2 = X(NXINIT)/ (ONE+X(NXINIT))
            DX = (Y2-Y1)/N
            do 13 I = 2,N
               Y = Y1 + (I-1)*DX
               X(I) = Y/ (ONE-ABS(Y))
   13       continue
            if (ABS(X((NXINIT+1)/2)).lt.TOL) X((NXINIT+1)/2) = ZERO
         end if

      end if

      JTOL = -LOG10(TOL) + HALF
      if (REG(1) .and. REG(2)) then
         MAXITS = 6
      else
         MAXITS = 3
      end if
C
C     Calculate H(*) and G(*).
C
      ITS = 1
      if (.not.AFIN) then
         if (X(2).ge.ZERO) X(2) = MAX(HALF*X(1),-ONE)
      end if

      if (.not.BFIN) then
         if (X(N).le.ZERO) X(N) = MIN(HALF*X(NXINIT),ONE)
      end if

      QMIN = 1.E31
      QMAX = -QMIN
   20 GAMMA = ZERO
C
C    Equidistribute { [Qmax - Q]^2 * max[abs(p') , abs(q') , abs(r')] }.
C
      EQMAX = ZERO
      EQMIN = 1.E31
      do 25 J = 1,N
         DX = X(J+1) - X(J)
         Y2 = X(J) + HALF*DX
         call COEFF(Y2,P2,Q2,R2)
         if ((P2.eq.ZERO) .or. (R2.eq.ZERO) .or. (P2*R2.lt.ZERO)) then
            FLAG = -15
            return

         end if

         Y1 = MAX(Y2-EPS,X(1)+TWO*U*ABS(X(1)))
         call COEFF(Y1,P1,Q1,R1)
         Y3 = MIN(Y2+EPS,X(NXINIT)-TWO*U*ABS(X(NXINIT)))
         call COEFF(Y3,P3,Q3,R3)
         if (LNF) then
            H(J) = ABS(Q3-Q1)/ (Y3-Y1)
            QLNF(J) = Q2/R2
         else
            H(J) = MAX(ABS(P3-P1),ABS(Q3-Q1),ABS(R3-R1))/ (Y3-Y1)
            Y1 = SQRT(SQRT(R1*P1))
            Y2 = SQRT(SQRT(R2*P2))
            Y3 = SQRT(SQRT(R3*P3))
            Y = SQRT(P2/R2)
            QLNF(J) = Q2/R2 + Y* ((Y3-Y1)* (SQRT(P3/R3)-SQRT(P1/R1))/
     +                FOUR+ (Y3-TWO*Y2+Y1)*Y)/ (Y2*EPS**2)
            if (ABS(QLNF(J)).le.EPS) QLNF(J) = ZERO
         end if

         QMAX = MAX(QMAX,QLNF(J))
         QMIN = MIN(QMIN,QLNF(J))
   25 continue
      Y = MAX(QMAX-QMIN,ONE)
      EV = 100.0 + MAX(ZERO,QMIN)
      do 30 J = 1,N
         DX = X(J+1) - X(J)
         if (QLNF(J).le.EV) then
            WEIGHT = 3.0* ((QMAX-QLNF(J))/Y)**2 + ONE
            H(J) = MAX(H(J)*WEIGHT,U)
         else
            Y2 = TWO*DX*SQRT(QLNF(J)-EV)
            if (Y2.le.UNDER) then
               WEIGHT = EXP(-Y2)
               H(J) = MAX(WEIGHT*H(J),U)
            else
               H(J) = U
            end if

         end if

         if ((.not.AFIN) .and. (X(J+1).lt.ZERO)) then
            H(J) = MAX(H(J)*EXP(X(J+1)),U)
            if ((J.eq.1) .and. (H(1).eq.U)) X(1) = HALF* (X(1)+X(2))
         end if

         if ((.not.BFIN) .and. (X(J).gt.ZERO)) then
            H(J) = MAX(H(J)*EXP(-X(J)),U)
            if ((J.eq.N) .and. (H(N).eq.U)) X(NXINIT) = HALF*
     +          (X(N)+X(NXINIT))
         end if

         EQMIN = MIN(EQMIN,H(J)*DX)
         EQMAX = MAX(EQMAX,H(J)*DX)
   30 continue
      NCOEFF = NCOEFF + 3*N
      if (EQMAX-EQMIN.le.MAX(TENTH*EQMAX,U/TENTH)) go to 75
C
C     Use a roughly locally quasi-uniform mesh.
C
      GAMMA = ZERO
      do 35 I = 1,N
         GAMMA = GAMMA + H(I)
   35 continue
      GAMMA = ONE/GAMMA
      do 45 I = 1,N
         Y = ZERO
         do 40 J = 1,N
            Y = MAX(Y,H(J)/ (ONE+GAMMA*ABS(X(I)+X(I+1)-X(J)-X(J+
     +          1))*H(J)))
   40    continue
         Z(I) = Y
   45 continue
      do 50 I = 1,N
         H(I) = Z(I)
   50 continue
      G(1) = ZERO
      do 55 J = 1,N
         G(J+1) = G(J) + H(J)* (X(J+1)-X(J))
   55 continue
      GAMMA = G(N+1)/N
C
C     Redistribution algorithm:
C
      Y = GAMMA
      I = 1
      do 65 J = 1,N
   60    if (Y.le.G(J+1)) then
            I = I + 1
            Z(I) = X(J) + (Y-G(J))/H(J)
            Y = Y + GAMMA
            go to 60

         end if

   65 continue
      Z(1) = X(1)
      Z(NXINIT) = X(NXINIT)
      DONE = .TRUE.
      do 70 J = 2,N
         if (ABS(Z(J)-X(J)).gt.TENTH* (Z(J+1)-Z(J-1))) DONE = .FALSE.
         X(J) = Z(J)
   70 continue
      if (.not.DONE) then
         if (ITS.lt.MAXITS) then
            ITS = ITS + 1
            go to 20

         end if

      end if

   75 if (.not.AFIN) then
         if (X(2).ge.ZERO) X(2) = -ONE
      end if

      if (.not.BFIN) then
         if (X(N).le.ZERO) X(N) = ONE
      end if

      do 95 K = 1,2
         NADD = 0
         if (KCLASS(K).gt.0) then
C
C           Add Nadd extra points near endpoint K.
C
            if (KCLASS(K).eq.1) then
               NADD = MIN(MAX((JTOL+2)/3,1),4)
               Y = TENTH
            end if

            if (KCLASS(K).eq.2) then
               NADD = MIN(MAX(JTOL,3),4)
               Y = MIN(MAX(TENTH** (JTOL/3),1.D-3),1.D-2)
            end if

            if (KCLASS(K).eq.3) NADD = MIN(MAX(JTOL/2,2),4)
            if (KCLASS(K).eq.4) then
               NADD = MIN(MAX((5+JTOL)/3,2),5)
               Y = TENTH** (NADD-1)
            end if

            if (KCLASS(K).eq.5) then
               NADD = (2**JTOL)**0.4
               NADD = MIN(MAX(NADD,3),6)
               Y = MIN(MAX((0.1**JTOL)** (0.333),0.001),0.1)
            end if

            if (KCLASS(K).eq.6) then
               NADD = MIN(MAX(JTOL,2),8)
               Y = 0.005
            end if

            if ((KCLASS(K).eq.7) .or. (KCLASS(K).eq.10)) then
               NADD = MIN(MAX((2*JTOL+6)/3,3),8)
               Y = MIN(MAX(SQRT(TENTH*TOL),1.D-6),1.D-2)
            end if

            if (KCLASS(K).eq.8) then
               NADD = (2**JTOL)**0.4
               NADD = MIN(MAX(NADD,2),6)
               Y = MIN(MAX((0.1**JTOL)** (0.333),0.001),0.05)
            end if

            if (KCLASS(K).eq.9) then
               if (LFLAG(5)) then
                  NADD = MIN(MAX((JTOL+4)/3,2),5)
                  Y = TENTH** (NADD-1)
               else
                  NADD = MIN(MAX(2+JTOL* (JTOL-3)/40,2),4)
                  Y = 0.25
               end if

            end if

            if (K.eq.1) then
               do 80 I = NXINIT,2,-1
                  X(I+NADD) = X(I)
   80          continue
               NXINIT = NXINIT + NADD
               if (AFIN) DX = X(2) - A
               do 85 I = 1,NADD
                  if (AFIN) then
                     X(I+1) = A + DX*Y** ((NADD-I+ONE)/NADD)
                  else
                     if (KCLASS(1).ne.1) then
                        X(NADD+2-I) = X(NADD+3-I) -
     +                                (X(NADD+4-I)-X(NADD+3-I))*2.4
                     else
                        X(NADD+2-I) = X(NADD+3-I) -
     +                                (X(NADD+4-I)-X(NADD+3-I))
                     end if

                  end if

   85          continue
            else
               if (BFIN) DX = B - X(NXINIT-1)
               N = NXINIT - 1
               NXINIT = NXINIT + NADD
               X(NXINIT) = B
               do 90 I = 1,NADD
                  if (BFIN) then
                     X(NXINIT-I) = B - DX*Y** ((NADD-I+ONE)/NADD)
                  else
                     if (KCLASS(2).ne.1) then
                        X(N+I) = X(N+I-1) + (X(N+I-1)-X(N+I-2))*2.4
                     else
                        X(N+I) = X(N+I-1) + (X(N+I-1)-X(N+I-2))
                     end if

                  end if

   90          continue
            end if

         end if

   95 continue
      if (.not.AFIN) X(1) = -ONE/U
      if (.not.BFIN) X(NXINIT) = ONE/U
      HMIN = X(NXINIT) - X(1)
      do 100 I = 2,NXINIT
         HMIN = MIN(HMIN,X(I)-X(I-1))
  100 continue
      return

      end subroutine MESH
C-----------------------------------------------------------------------------
      subroutine POWER(X,F,N,TOL,IPRINT,EF,CF,OSC,EXACT,Y,IFLAG)
C
C     Find the power function which "dominates" the tabled
C     coefficient function.  The output is Cf and Ef such that
C           f(x)  is asymptotic to  Cf*x^Ef .
C     The vectors X(*) and F(*) hold the N input points:
C           F(I) = f(X(I)) I = 1,...,N.
C     Set IFLAG = 0 for normal return; 1 for uncertainty in Ef;
C     2 if uncertain about Cf (oscillatory).
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0)
C     ..
C     .. Scalar Arguments ..
      double precision CF,EF,TOL
      integer IFLAG,IPRINT,N
      logical EXACT,OSC
C     ..
C     .. Array Arguments ..
      double precision F(*),X(*),Y(*)
C     ..
C     .. Local Scalars ..
      double precision ERROR,TOLAIT,TOLMIN
      integer K,NY
C     ..
C     .. External Subroutines ..
cc      external AITKEN
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,LOG,MAX,MIN,SIGN,SQRT
C     ..
C     .. Data statements ..
      data TOLMIN/1.D-6/
C     ..
C
C     Estimate the exponent.
C
      OSC = .FALSE.
      NY = N - 1
      ERROR = 1.E30
      TOLAIT = MIN(TOLMIN,TOL)
      do 10 K = 1,NY
         if ((F(K).ne.ZERO) .and. (F(K+1).ne.ZERO)) then
            Y(K) = LOG(ABS(F(K+1)/F(K)))/LOG(ABS(X(K+1)/X(K)))
         else
            Y(K) = ZERO
         end if

   10 continue
      EF = Y(NY)
      if (IPRINT.gt.4) then
         write (21,FMT=*) ' From POWER; E_k and c_k sequences:'
         write (21,FMT=15) (Y(K),K=1,NY)

   15    format (4d19.10)

         write (21,FMT=*)
      end if

      call AITKEN(EF,TOLAIT,NY,Y,ERROR)
      K = EF + SIGN(HALF,EF)
      if (ABS(K-EF).le.SQRT(TOL)) EF = K
      if (ABS(ERROR).gt.TOL*MAX(ONE,ABS(EF))) then
C
C        There is uncertainty in the exponent.
C
         IFLAG = 1
      end if

      if (ABS(EF).le.TOL) EF = ZERO
C
C     Estimate the coefficient.
C
      do 20 K = 1,N - 1
         Y(K) = F(K)/ABS(X(K))**EF
   20 continue
      CF = Y(N-1)
      if (IPRINT.ge.5) write (21,FMT=15) (Y(K),K=1,N-1)
      call AITKEN(CF,TOLAIT,N-1,Y,ERROR)
      if ((EF.gt.20.) .and. (ABS(CF).le.TOL)) then
C
C        Coefficient probably has exponential behavior.
C
         CF = SIGN(ONE,Y(N-1))
      else
         if ((ABS(ERROR).gt.TOL*MAX(ONE,ABS(CF))) .or.
     +       ((ABS(F(N)-CF*X(N)**EF).gt.20.0*TOL*ABS(F(N))).and.
     +       (EF.ne.ZERO))) then
C
C           There is uncertainty in the coefficient; call such
C           cases oscillatory.
C
            IFLAG = 2
            OSC = .TRUE.
         end if

      end if

      if (ABS(CF).gt.1.D7) then
         EXACT = .FALSE.
         return

      end if

      K = CF + HALF
      if ((ABS(K-CF).le.SQRT(TOL)) .and. (K.ne.0)) CF = K
      EXACT = .TRUE.
      do 30 K = 1,N
         if (ABS(F(K)-CF*X(K)**EF).gt.TOL*ABS(F(K))) EXACT = .FALSE.
   30 continue
      return

      end subroutine POWER
C----------------------------------------------------------------------
      subroutine PQRINT(X,SQRTRP,QLNF)
C
C     Evaluate the integrands needed for the asymptotic formulas.
C       (1) The Liouville normal form potential Qlnf:
C             Qlnf(t) =  q/r + f"(t)/f   with   f = (pr)**.25  .
C       (2) The term in the change of independent variable:
C                        sqrt(r/p) .
C
C     .. Parameters ..
      double precision ZERO,TWO,FOUR
      parameter (ZERO=0.0,TWO=2.0,FOUR=4.0)
C     ..
C     .. Scalar Arguments ..
      double precision QLNF,SQRTRP,X
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision EPS,FL,FM,FR,PX,QX,RX,XDOTL,XDOTR,Z
C     ..
C     .. External Subroutines ..
cc      external COEFF
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MIN,SQRT
C     ..
C     .. Common blocks ..
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      if (LNF) then
         call COEFF(X,PX,QX,RX)
         if ((PX.eq.ZERO) .or. (RX.eq.ZERO) .or. (PX*RX.lt.ZERO)) then
            FLAG = -15
            return

         end if

         NCOEFF = NCOEFF + 1
         QLNF = QX/RX
         SQRTRP = SQRT(RX/PX)
      else
         EPS = MIN(1.D-4,MIN(ABS(B-X),ABS(X-A))/TWO)
         Z = X - EPS
         call COEFF(Z,PX,QX,RX)
         if ((PX.eq.ZERO) .or. (RX.eq.ZERO) .or. (PX*RX.lt.ZERO)) then
            FLAG = -15
            return

         end if

         XDOTL = SQRT(PX/RX)
         FL = SQRT(SQRT(PX*RX))
         Z = X + EPS
         call COEFF(Z,PX,QX,RX)
         XDOTR = SQRT(PX/RX)
         FR = SQRT(SQRT(PX*RX))
         call COEFF(X,PX,QX,RX)
         SQRTRP = SQRT(RX/PX)
         FM = SQRT(SQRT(PX*RX))
         NCOEFF = NCOEFF + 3
         QLNF = QX/RX + ((FR-FM)* (XDOTR-XDOTL)/FOUR+ (FR-TWO*FM+FL)/
     +          SQRTRP)/ (EPS*EPS*FM*SQRTRP)
         if (ABS(QLNF).le.EPS) QLNF = ZERO
      end if

      return

      end subroutine PQRINT
C-----------------------------------------------------------------------
      subroutine REGULR(JOB,JOBMSH,TOL,NEV,EV,IPRINT,NEXTRP,XEF,EF,PDEF,
     +                  HMIN,STORE)
***********************************************************************
*                                                                     *
*     REGULR calculates Sturm-Liouville eigenvalue and (optionally)   *
*     eigenfunction estimates for the problem described initially.    *
*                                                                     *
***********************************************************************
C
C    Input parameters:
C      JOB       = logical variable describing tasks to be carried out.
C                  JOB = .True. iff an eigenfunction is to be calculated.
C      JOBMSH    = logical variable, JOBMSH = .True. iff initial mesh
C                  is a function of the eigenvalue index.
C      TOL(*)    = real vector of 6 tolerances.
C                  TOL(1) is the absolute error tolerance for e-values,
C                  TOL(2) is the relative error tolerance for e-values,
C                  TOL(3) is the abs. error tolerance for e-functions,
C                  TOL(4) is the rel. error tolerance for e-functions,
C                  TOL(5) is the abs. error tolerance for e-function
C                         derivatives,
C                  TOL(6) is the rel. error tolerance for e-function
C                         derivatives.
C                  Eigenfunction tolerances need not be set if JOB is
C                  False.  All absolute error tolerances must be
C                  positive; all relative must be at least 100 times
C                  the unit roundoff.
C      NEV       = integer index for the eigenvalue sought; NEV .GE. 0 .
C      EV        = real initial guess for eigenvalue NEV; accuracy is
C                  not at all critical, but if a good estimate is
C                  available some time may be saved.
C      IPRINT    = integer controlling amount of internal printing done.
C
C    Output parameters:
C      EV        = real computed approximation to NEVth eigenvalue.
C      XEF(*)    = real vector of points for eigenfunction output.
C      EF(*)     = real vector of eigenfunction values: EF(i) is the
C                  estimate of u(XEF(i)).  If JOB is False then this
C                  vector is not referenced.
C      PDEF(*)   = real vector of eigenfunction derivative values:
C                  PDEF(i) is the estimate of (pu')(XEF(i)).  If JOB is
C                  False then this vector is not referenced.
C
C    Auxiliary storage:
C      STORE(*) = real vector of auxiliary storage, must be dimensioned
C                 at least max[100,26N]. (N the number of mesh points)
C
C     Storage allocation in auxiliary vector (currently Maxlvl = 10):
C         STORE(*)
C       1   ->      N           vector of mesh points X(*),
C     N+1   ->     2N           best current eigenfunction values,
C    2N+1   ->     3N           best current derivative values,
C    3N+1   ->     4N           scale factors in GETEF,
C    4N+1   ->  (6+2*Maxlvl)N   intermediate eigenfunction values.
C-----------------------------------------------------------------------
C     Local variables:
C
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,THREE,FIVE,TEN,TOLMIN
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,THREE=3.0,FIVE=5.0,
     +          TEN=10.0,TOLMIN=1.D-3)
C     ..
C     .. Scalar Arguments ..
      double precision EV,HMIN
      integer IPRINT,NEV,NEXTRP
      logical JOB,JOBMSH
C     ..
C     .. Array Arguments ..
      double precision EF(*),PDEF(*),STORE(*),TOL(*),XEF(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision ABSERR,ALPHA1,ALPHA2,BETA1,BETA2,DELTA,EFNORM,
     +                 ERROR,EVHAT,EVHIGH,EVLOW,FHIGH,FLOW,H,PDUMAX,
     +                 QINT,QLNF,RELERR,RPINT,SQRTRP,TOL1,TOL2,TOL3,
     +                 TOL4,TOL5,TOL6,TOLEXT,TOLPDU,TOLSUM,UMAX,Z
      integer I,II,J,JDU,JU,KDU,KK,KU
      logical DONE,EFDONE,EVDONE,EXFULL
C     ..
C     .. Local Arrays ..
      double precision EVEXT(20)
      logical EFIN(2,2)
C     ..
C     .. External Functions ..
cc      double precision ASYMEV
cc      external ASYMEV
C     ..
C     .. External Subroutines ..
cc      external BRCKET,EXTRAP,GETEF,MESH,PQRINT,ZZERO
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,MIN
C     ..
C     .. Common blocks ..
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
      ALPHA1 = ZERO
      ALPHA2 = ONE
      BETA1 = ZERO
      BETA2 = ONE
      if (NEV.lt.0) then
         FLAG = -37
         return

      end if

      EVDONE = .FALSE.
      TOL1 = MIN(TOL(1),TOLMIN)
      TOL2 = MIN(TOL(2),TOLMIN)
      if (.not. (REG(1).and.REG(2))) then
         TOL1 = TOL1/THREE
         TOL2 = TOL2/THREE
      end if

      if (JOB) then
         TOL3 = MIN(TOL(3),TOLMIN)
         TOL4 = MIN(TOL(4),TOLMIN)
         TOL5 = MIN(TOL(5),TOLMIN)
         TOL6 = MIN(TOL(6),TOLMIN)
         ABSERR = TOL1/FIVE
         RELERR = TOL2/FIVE
         EXFULL = .TRUE.
      else
         EFDONE = .TRUE.
         ABSERR = TOL1/TEN
         RELERR = TOL2/TEN
         EXFULL = .FALSE.
      end if

      TOLSUM = TOL1 + TOL2
      if (JOBMSH) then
         if (NEV.ge.0) then
            II = NXINIT + 16
            call MESH(.TRUE.,NEV,STORE,STORE(II),STORE(2*II+1),
     +                STORE(3*II+1),STORE(4*II+1),TOLSUM,HMIN)
         end if

         if (IPRINT.ge.1) then
            write (*,FMT=15) (STORE(I),I=1,NXINIT)
            write (21,FMT=15) (STORE(I),I=1,NXINIT)

   15       format (' Level 0 mesh:',/ (5g15.6))

         end if

      end if
C
C     Compute estimates for integrals in asymptotic formulas (accuracy
C     is not all that critical).
C
      QINT = ZERO
      RPINT = ZERO
      do 20 I = 2,NXINIT
         H = STORE(I) - STORE(I-1)
         Z = STORE(I-1) + HALF*H
         if ((.not.AFIN) .and. (I.eq.2)) then
            H = STORE(3) - STORE(2)
            Z = STORE(2)
         end if

         if ((.not.BFIN) .and. (I.eq.NXINIT)) then
            H = STORE(NXINIT-1) - STORE(NXINIT-2)
            Z = STORE(NXINIT-1)
         end if

         call PQRINT(Z,SQRTRP,QLNF)
         if (FLAG.lt.0) return
         QINT = QINT + H*QLNF
         RPINT = RPINT + H*SQRTRP
   20 continue
      if (QINT.gt.ONE/U) QINT = ZERO
      if (RPINT.gt.ONE/U) RPINT = ZERO
C
C     Loop over the levels.
C
      do 60 LEVEL = 0,MAXLVL
         if (HMIN/2**LEVEL.le.TEN*U) then
            FLAG = -8
            go to 70

         end if
C
C        Find a bracket for the Nevth eigenvalue.
C
         if (LEVEL.eq.0) then
            EV = ASYMEV(NEV,QINT,RPINT,ALPHA1,ALPHA2,BETA1,BETA2)
            EV = MAX(EV,ZERO)
            DELTA = HALF
         else
            DELTA = MAX(TOLSUM*ABS(EVHAT),HALF*HALF*DELTA)
            if (LEVEL.gt.1) then
               ERROR = (EVEXT(LEVEL)-EVEXT(LEVEL-1))/THREE
               if (ABS(ERROR).le.100.) then
                  EV = EVHAT + ERROR
               else
                  DELTA = ONE
                  EV = EVHAT
               end if

            else
               EV = EVHAT
            end if

         end if

         EVLOW = EV - DELTA
         EVHIGH = EV + DELTA
         if (IPRINT.ge.4) write (21,FMT=25) EVLOW,EVHIGH

   25    format ('      In  bracket:',2d24.15)

         call BRCKET(NEV,EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,STORE)
         if (IPRINT.ge.4) write (21,FMT=30) EVLOW,EVHIGH

   30    format ('      Out bracket:',2d24.15)

         DELTA = HALF* (EVHIGH-EVLOW)
         if (FLAG.lt.0) return
         if (ABS(EVHIGH-EVLOW).gt.MAX(ABSERR,RELERR*ABS(EVHIGH))) then
            call ZZERO(EVLOW,EVHIGH,FLOW,FHIGH,ABSERR,RELERR,J,STORE)
            if (J.ne.0) then
               FLAG = -7
               return

            end if

         end if

         EVHAT = MIN(EVLOW,EVHIGH)
         if (IPRINT.ge.1) then
            write (*,FMT=40) LEVEL,EVHAT
            write (21,FMT=40) LEVEL,EVHAT

   40       format (' Level ',i3,' ;    EvHat = ',d24.15)

         end if

         EV = EVHAT
         TOLEXT = MAX(TOL1,ABS(EV)*TOL2)
         call EXTRAP(EV,TOLEXT,LEVEL+1,NEXTRP,EXFULL,.TRUE.,0,EVEXT,
     +               IPRINT,ERROR,EVDONE)
         if (JOB) then
            call GETEF(EVHAT,EFNORM,IPRINT,STORE,EFIN)
            if (LEVEL.eq.0) then
               UMAX = ONE
               PDUMAX = ONE
C
C              Set pointers to STORE(*).
C
               JU = NXINIT
               JDU = 2*NXINIT
               KU = 4*NXINIT
               KDU = (MAXLVL+5)*NXINIT
               KK = MAXLVL + 1
            end if
C
C           Extrapolate eigenfunction values.
C
            TOLEXT = MAX(TOL3,UMAX*TOL4)
            TOLPDU = MAX(TOL5,PDUMAX*TOL6)
            EFDONE = .TRUE.
            if (AFIN) then
               UMAX = ABS(STORE(JU+1))
               PDUMAX = ABS(STORE(JDU+1))
            else
               UMAX = ZERO
               PDUMAX = ZERO
            end if

            if (EFIN(1,1)) then
               call EXTRAP(STORE(JU+1),TOLEXT,LEVEL+1,NEXTRP,EXFULL,
     +                     .TRUE.,0,STORE(KU+1),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
            end if

            if (EFIN(2,1)) then
               call EXTRAP(STORE(JDU+1),TOLPDU,LEVEL+1,NEXTRP,EXFULL,
     +                     .TRUE.,0,STORE(KDU+1),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
            end if

            do 50 I = 2,NXINIT - 1
               II = KK* (I-1) + 1
               call EXTRAP(STORE(JU+I),TOLEXT,LEVEL+1,NEXTRP,EXFULL,
     +                     .TRUE.,0,STORE(KU+II),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
               call EXTRAP(STORE(JDU+I),TOLPDU,LEVEL+1,NEXTRP,EXFULL,
     +                     .TRUE.,0,STORE(KDU+II),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
               UMAX = MAX(UMAX,ABS(STORE(JU+I)))
               PDUMAX = MAX(PDUMAX,ABS(STORE(JDU+I)))
   50       continue
            II = KK* (NXINIT-1) + 1
            if (EFIN(1,2)) then
               call EXTRAP(STORE(JU+NXINIT),TOLEXT,LEVEL+1,NEXTRP,
     +                     EXFULL,.TRUE.,0,STORE(KU+II),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
            end if

            if (EFIN(2,2)) then
               call EXTRAP(STORE(JDU+NXINIT),TOLPDU,LEVEL+1,NEXTRP,
     +                     EXFULL,.TRUE.,0,STORE(KDU+II),0,ERROR,DONE)
               EFDONE = EFDONE .and. DONE
            end if

            if (BFIN) then
               UMAX = MAX(UMAX,ABS(STORE(JU+NXINIT)))
               PDUMAX = MAX(PDUMAX,ABS(STORE(JDU+NXINIT)))
            end if

            ABSERR = MAX(HALF*ABSERR,TEN*U)
            RELERR = MAX(HALF*RELERR,TEN*U)
         end if

         if (EVDONE .and. (LEVEL.ge.2) .and. EFDONE) go to 70
   60 continue
      if (.not.EVDONE) then
         FLAG = -1
         return

      end if
C
C     Unload eigenfunction values.
C
   70 if (JOB) then
         do 80 I = 1,NXINIT
            XEF(I) = STORE(I)
            EF(I) = STORE(JU+I)
            PDEF(I) = STORE(JDU+I)
   80    continue
         if ((FLAG.ge.0) .and. (.not.EFDONE)) FLAG = -2
      end if

      return

      end subroutine REGULR
C---------------------------------------------------------------------
      subroutine SHOOT(EV,X,MU,FEV)
C
C     .. Parameters ..
      double precision ZERO,HALF,ONE,TWO,PI
      parameter (ZERO=0.0,HALF=0.5D0,ONE=1.0,TWO=2.0,
     +          PI=3.141592653589793D0)
C     ..
C     .. Scalar Arguments ..
      double precision EV,FEV
      integer MU
C     ..
C     .. Array Arguments ..
      double precision X(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision CHI,DPSI,DV,H,HALFH,HOMEGA,OMEGA,PDV,PHASE,PN,
     +                 PSI,RN,SA1,SA2,SB1,SB2,SCALE,SGN,T,TAU,V,VNEW,X2,
     +                 XLEFT,Z
      integer I,IMAX,J,K1,K2,KLVL,MODE,NSAVE,NZERO
C     ..
C     .. External Subroutines ..
cc      external STEP
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,ATAN,INT,MAX,MIN,MOD,SIGN
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Data statements ..
      data IMAX/1000000/
C     ..
C
C     Make one shot across (A,B) for the current mesh using the scaled
C     variable v(x). Count zeros if COUNTZ is True.
C
C     Shoot from x=A to x=B.
C
      V = A2 - A2P*EV
      PDV = A1 - A1P*EV
      NZERO = 0
      NSAVE = NSGNF
      SA1 = A1
      SA2 = A2
      SB1 = B1
      SB2 = B2
      KLVL = 2**LEVEL
      MODE = 0
C
C     Modify base sign count if necessary.
C
      if (((KCLASS(1).gt.8).or. (KCLASS(2).gt.8)) .and.
     +    (EV.lt.CUTOFF)) then
         if (KCLASS(1).gt.8) then
            SA1 = D(4,1)
            SA2 = ONE
            V = SA2
            PDV = SA1
            MODE = 1
         end if

         if (KCLASS(2).gt.8) then
            SB1 = D(4,2)
            SB2 = ONE
         end if

         SGN = A2*B2
         if (SGN.eq.ZERO) then
            SGN = SA1*SB2 + SA2*SB1
            if (SGN.eq.ZERO) SGN = A1*B1
         end if

         NSGNF = SIGN(ONE,SGN)
      end if

      if (.not.AFIN) MODE = 3
      SCALE = MAX(ABS(V),ABS(PDV))
      V = V/SCALE
      PDV = PDV/SCALE
      do 20 I = 2,NXINIT
         XLEFT = X(I-1)
         H = X(I) - XLEFT
         if (MODE.eq.1) then
            H = (H/D(1,1))**ETA(1,1)
            XLEFT = ZERO
         end if

         if (MODE.eq.3) then
            XLEFT = ZERO
            H = -ONE/X(2)
         end if

         if ((KCLASS(2).gt.8) .and. (I.eq.NXINIT) .and.
     +       (EV.lt.CUTOFF)) then
C
C           Convert from u(x) to V(t) near x=b.
C
            MODE = 2
            T = (H/D(1,2))**ETA(1,2)
            CHI = D(2,2)*T**ETA(2,2)
            V = V/CHI
            PN = CP(2)*H** (EP(2)-ONE)
            PDV = (PDV/ (PN*CHI*ETA(1,2))+ETA(2,2)*V)*
     +            T** (TWO*EMU(2)-ONE)
            H = T
            XLEFT = -H
         end if

         if ((.not.BFIN) .and. (I.eq.NXINIT)) then
            MODE = 4
            H = ONE/X(I-1)
            XLEFT = -H
         end if

         H = H/KLVL
         HALFH = HALF*H
         do 10 J = 1,KLVL
            Z = XLEFT + HALFH
            call STEP(Z,H,EV,PN,RN,TAU,OMEGA,HOMEGA,PSI,DPSI,SCALE,MODE)
            if (FLAG.lt.0) then
               FEV = ZERO
               return

            end if

            DV = PDV/PN
            VNEW = DPSI*V + PSI*DV
            XLEFT = XLEFT + H
            if (COUNTZ) then
C
C              Count zeros of v(x).
C
               if (TAU.le.ZERO) then
                  if (VNEW*V.lt.ZERO) NZERO = NZERO + 1
               else
                  if (DV.eq.ZERO) then
                     NZERO = NZERO + INT(HALF+HOMEGA/PI)
                  else
                     PHASE = ATAN(V*OMEGA/DV)
                     K1 = PHASE/PI
                     X2 = (PHASE+HOMEGA)/PI
                     if (X2.lt.IMAX) then
                        K2 = X2
                        NZERO = NZERO + K2 - K1
                     else
                        NZERO = IMAX
                     end if

                     if (PHASE* (PHASE+HOMEGA).lt.ZERO) NZERO = NZERO +
     +                   1
                  end if

                  NZERO = MIN(IMAX,NZERO)
               end if

            end if

            PDV = -PN*TAU*PSI*V + DPSI*PDV
            V = VNEW
   10    continue
         if (MODE.eq.1) then
C
C           Convert from V(t) back to u(x) near x=a.
C
            PN = CP(1)* (X(2)-A)** (EP(1)-ONE)
            T = ((X(2)-A)/D(1,1))**ETA(1,1)
            CHI = D(2,1)*T**ETA(2,1)
            PDV = PN*ETA(1,1)*CHI* (ETA(2,1)*V+PDV*T** (ONE-TWO*EMU(1)))
            V = CHI*V
         end if

         MODE = 0
   20 continue
      FEV = SB1*V + SB2*PDV
      if (COUNTZ) then
C
C        Adjust zero count.
C
         MU = NZERO
         if (A2P.ne.ZERO) then
            if (EV.ge.SA2/A2P) MU = NZERO + 1
         end if

         if (SB2.ne.ZERO) then
            SGN = SIGN(ONE,NSGNF*FEV)
            if (MOD(MU,2).eq.1) SGN = -SGN
            if (SGN.lt.ZERO) MU = MU + 1
         end if

      end if

      NSGNF = NSAVE
      return

      end subroutine SHOOT
C-----------------------------------------------------------------------
      subroutine START(JOB,CONS,TOL,NEV,INDXEV,N,XEF,NUMT,T,NEXTRP,X)
C
C     This routine tests the input data, initializes the labeled
C     common blocks, and generates the first mesh.  Check
C        eigenfunction tolerances iff JOB(1) is True ,
C        XEF(*) iff JOB(2) is True,
C        NUMT, T(*) iff JOB(3) is True ,
C
C     .. Parameters ..
      double precision ZERO,HUNDRD
      parameter (ZERO=0.0,HUNDRD=100.0)
C     ..
C     .. Scalar Arguments ..
      integer N,NEV,NEXTRP,NUMT
C     ..
C     .. Array Arguments ..
      double precision CONS(*),T(*),TOL(*),X(*),XEF(*)
      integer INDXEV(*)
      logical JOB(*)
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision UFLOW,URN
      integer I,K
C     ..
C     .. Intrinsic Functions ..
      intrinsic LOG10,MAX,MIN
C     ..
C     .. Common blocks ..
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Data statements ..
************************************************************************
*   In the following DATA statement, initialize URN to an estimate for *
*   the unit roundoff and UFLOW to a value somewhat less than          *
*   -ln(underflow).  E.g.,                                             *
*      for IEEE double precision use URN = 2.D-16, UFLOW = 650;        *
*      for VAXen double precision use URN = 1.D-17, UFLOW = 85;        *
*      for Crays (single precision) use URN = 7.D-15, UFLOW = 5000.    *
*   Exact values are not at all critical.                              *
*   Here, we assume IEEE double precision:                             *
*                                                                      *
      data URN/2.D-16/,UFLOW/650.0/
C     ..
C***********************************************************************
      U = URN
      UNDER = UFLOW
C
C     Initialize.
C
      A1 = CONS(1)
      A1P = CONS(2)
      A2 = CONS(3)
      A2P = CONS(4)
      A = CONS(7)
      B1 = CONS(5)
      B2 = CONS(6)
      B = CONS(8)
      NCOEFF = 0
      FLAG = 0
C
C     Test input.
C
      if (AFIN .and. BFIN .and. (A.ge.B)) FLAG = -34
      if (TOL(1).le.ZERO) FLAG = -35
      if (TOL(2).lt.HUNDRD*U) FLAG = -36
      if (JOB(1)) then
         if (TOL(3).le.ZERO) FLAG = -35
         if (TOL(4).lt.HUNDRD*U) FLAG = -36
         if (TOL(5).le.ZERO) FLAG = -35
         if (TOL(6).lt.HUNDRD*U) FLAG = -36
      end if

      if (JOB(2)) then
         if (N.eq.1) then
            FLAG = -30
         else
            if (AFIN .and. (XEF(2).le.A)) FLAG = -39
            do 10 I = 3,N - 1
               if (XEF(I-1).ge.XEF(I)) FLAG = -39
   10       continue
            if (BFIN .and. (XEF(N-1).gt.B)) FLAG = -39
         end if

         if ((.not.AFIN) .and. (XEF(2).ge.ZERO)) FLAG = -39
         if ((.not.BFIN) .and. (XEF(N-1).le.ZERO)) FLAG = -39
      end if

      if ((JOB(2).or.JOB(3)) .and. (NEV.gt.0)) then
         do 20 I = 1,NEV
            if (INDXEV(I).lt.0) FLAG = -37
   20    continue
      end if

      if (JOB(3)) then
         if (NUMT.le.0) FLAG = -38
         do 30 I = 2,NUMT
            if (T(I).le.T(I-1)) FLAG = -39
   30    continue
      end if

      if (FLAG.lt.0) return
C
C     Set MAXEXT, the maximum number of extrapolations allowed, and
C         MAXINT, the maximum number of intervals in X(*) allowed when
C         mesh is chosen by START.
C
C     IMPORTANT NOTE: the size of various fixed arrays in this package
C     depends on the value of MAXEXT in this FORTRAN77 implementation.
C     If MAXEXT is increased, then more storage may have to be allocated
C     to the columns of R(*,*) in EXTRAP.
C
      MAXEXT = 6
      MAXINT = 31
C
C     Calculate maximum number of columns in extrapolation table
C     and the maximum number of levels allowed.
C
      K = -LOG10(TOL(2))
      I = MAX(K+3,0)
      NEXTRP = MIN(MAX(3,I/2),MAXEXT)
C
C     Calculate the initial mesh.
C
      if (N.gt.0) then
         NXINIT = N
      else
         NXINIT = MIN(2*NEXTRP+3,MAXINT)
         if (JOB(1)) N = NXINIT
      end if

      if (JOB(2)) then
         A = XEF(1)
         B = XEF(NXINIT)
         do 40 I = 1,NXINIT
            X(I) = XEF(I)
   40    continue
      end if

      return

      end subroutine START
C-----------------------------------------------------------------------
      subroutine STEP(X,H,EV,PX,RX,TAU,OMEGA,HOMEGA,PSI,DPSI,SCLOG,MODE)
C
C    Evaluate the coefficient functions, the scaled basis function PSI,
C    its derivative DPSI, and the log of the scale factor SCLOG.
C
C
C     .. Parameters ..
      double precision ZERO,HNDRTH,HALF,ONE,TWO,SIX,TWELVE,TWENTY
      parameter (ZERO=0.0,HNDRTH=.01,HALF=0.5D0,ONE=1.0,TWO=2.0,SIX=6.0,
     +          TWELVE=12.0,TWENTY=20.0)
C     ..
C     .. Scalar Arguments ..
      double precision DPSI,EV,H,HOMEGA,OMEGA,PSI,PX,RX,SCLOG,TAU,X
      integer MODE
C     ..
C     .. Scalars in Common ..
      double precision A,A1,A1P,A2,A2P,B,B1,B2,CUTOFF,U,UNDER
      integer FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
C      logical AFIN,BFIN,COUNTZ,LNF
C     ..
C     .. Arrays in Common ..
      double precision CP(2),CR(2),D(4,2),EMU(2),EP(2),EQLNF(2),ER(2),
     +                 ETA(2,2),PNU(2)
      integer KCLASS(2)
C      logical LC(2),LFLAG(6),OSC(2),REG(2)
C     ..
C     .. Local Scalars ..
      double precision DX,FP,FR,OVER,QX,T,TMU,Z
C     ..
C     .. External Subroutines ..
cc      external COEFF
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,COS,LOG,MAX,MIN,SIN,SQRT,TANH
C     ..
C     .. Common blocks ..
      common /SLCLSS/CP,CR,CUTOFF,D,EMU,EP,EQLNF,ER,ETA,PNU,KCLASS
      common /SLINT/FLAG,LEVEL,MAXEXT,MAXINT,MAXLVL,NCOEFF,NSGNF,NXINIT
C      common /SLLOG/AFIN,BFIN,COUNTZ,LFLAG,LNF,LC,OSC,REG
      common /SLREAL/A1,A1P,A2,A2P,B1,B2,A,B,U,UNDER
C     ..
C     .. Data statements ..
      data OVER/1.D8/
C     ..
C
C    Evaluate the coefficient functions at X and calculate TAU.  The
C    error flag FLAG is zero for a successful calculation; if p(x)
C    or r(x) are zero, then FLAG is set to -15.  If the argument for
C    a trig function exceeds OVER then FLAG is set to -10.
C         Proceed normally when Mode = 0; when Mode = 1 or 2 use the
C    change of variable for "hard" problems; when Mode = 3 or 4 use
C    the change of variable t = -1/x near infinity.
C
      if (MODE.eq.0) then
         call COEFF(X,PX,QX,RX)
      else
         if (MODE.gt.2) then
            T = X
            X = -ONE/T
            call COEFF(X,PX,QX,RX)
            T = T*T
            PX = T*PX
            QX = QX/T
            RX = RX/T
         else
            T = X
            if (ETA(1,MODE).eq.ONE) then
               DX = D(1,MODE)*ABS(T)
            else
               DX = D(1,MODE)* (ABS(T))** (ONE/ETA(1,MODE))
            end if

            if (MODE.eq.1) then
               Z = A + DX
            else
               Z = B - DX
            end if

            call COEFF(Z,PX,QX,RX)
         end if

      end if

      NCOEFF = NCOEFF + 1
      if ((PX.eq.ZERO) .or. (RX.eq.ZERO)) then
         FLAG = -15
         return

      end if

      if ((MODE.eq.1) .or. (MODE.eq.2)) then
         if (EMU(MODE).eq.HALF) then
            TMU = ABS(T)
         else
            TMU = (T*T)**EMU(MODE)
         end if

         FP = CP(MODE)
         if (EP(MODE).ne.ZERO) FP = FP*DX**EP(MODE)
         FR = CR(MODE)
         if (ER(MODE).ne.ZERO) FR = FR*DX**ER(MODE)
         PX = PX*TMU/FP
         QX = (QX/FR-D(3,MODE)/T**2)*TMU
         RX = RX*TMU/FR
      end if

      TAU = (EV*RX-QX)/PX
      OMEGA = SQRT(ABS(TAU))
      HOMEGA = H*OMEGA
      SCLOG = ZERO
C
C     Evaluate the scaled basis functions.
C
      if (HOMEGA.gt.HNDRTH) then
         if (TAU.gt.ZERO) then
            if (HOMEGA.gt.OVER) then
               FLAG = -10
               return

            end if

            DPSI = COS(HOMEGA)
            PSI = SIN(HOMEGA)/OMEGA
         else
            SCLOG = HOMEGA
            if (HOMEGA.lt.UNDER) then
               T = TANH(HOMEGA)
               DPSI = ONE/ (ONE+T)
               PSI = T*DPSI/OMEGA
            else
               SCLOG = MIN(SCLOG,TWO*UNDER)
               DPSI = HALF
               PSI = DPSI/OMEGA
            end if

         end if

      else
         T = TAU*H*H
         DPSI = ONE + T* (T/TWELVE-ONE)/TWO
         PSI = H* (ONE+T* (T/TWENTY-ONE)/SIX)
         if (T.lt.ZERO) then
            T = MAX(ABS(PSI),ABS(DPSI))
            SCLOG = LOG(T)
            PSI = PSI/T
            DPSI = DPSI/T
         end if

      end if

      return

      end subroutine STEP
C-----------------------------------------------------------------------
      subroutine ZZERO(B,C,FB,FC,ABSERR,RELERR,IFLAG,X)
C
C
C  ZZERO computes a root of F.  The method used is a combination of
C  bisection and the secant rule.  This code is adapted from one in
C  the text "Foundations of Numerical Computing" written by Allen,
C  Pruess, and Shampine.
C
C  Input parameters:
C     B,C   = values of X such that F(B)*F(C) .LE. 0.
C     FB,FC = values of F at input B and C, resp.
C     ABSERR,RELERR = absolute and relative error tolerances.  The
C             stopping criterion is:
C               ABS(B-C) .LE. 2.0*MAX(ABSERR,ABS(B)*RELERR).
C  Output parameters:
C     B,C   = see IFLAG returns.
C     FB    = value of final residual F(B).
C     IFLAG = 0 for normal return; F(B)*F(C) .LT. 0 and the
C               stopping criterion is met (or F(B)=0).  B always
C               satisfies ABS(F(B)) .LE. ABS(F(C)).
C           = 1 if too many function evaluations were made; in this version
C               200 are allowed.
C           =-2 if F(B)*F(C) is positive on input.
C
C  Local variables:
C
C  Internal constants
C
C
C  Initialization.
C
C     .. Parameters ..
      double precision ZERO,ONE,TWO,EIGHT
      integer MAXF
      parameter (ZERO=0.0,ONE=1.0,TWO=2.0,EIGHT=8.0,MAXF=200)
C     ..
C     .. Scalar Arguments ..
      double precision ABSERR,B,C,FB,FC,RELERR
      integer IFLAG
C     ..
C     .. Array Arguments ..
      double precision X(*)
C     ..
C     .. Local Scalars ..
      double precision A,ACMB,CMB,FA,P,Q,TOL,WIDTH
      integer KOUNT,MU,NF
C     ..
C     .. External Subroutines ..
cc      external SHOOT
C     ..
C     .. Intrinsic Functions ..
      intrinsic ABS,MAX,SIGN
C     ..
      KOUNT = 0
      WIDTH = ABS(B-C)
      A = C
      FA = FC
      if (SIGN(ONE,FA).eq.SIGN(ONE,FB)) then
         IFLAG = -2
         return

      end if

      FC = FA
      NF = 2
   20 if (ABS(FC).lt.ABS(FB)) then
C
C  Interchange B and C so that ABS(F(B)) .LE. ABS(F(C)).
C
         A = B
         FA = FB
         B = C
         FB = FC
         C = A
         FC = FA
      end if

      CMB = (C-B)/TWO
      ACMB = ABS(CMB)
      TOL = MAX(ABSERR,ABS(B)*RELERR)
C
C  Test stopping criterion and function count.
C
      if (ACMB.le.TOL) then
         IFLAG = 0
         return

      end if

      if (NF.ge.MAXF) then
         IFLAG = 1
         return

      end if
C
C  Calculate new iterate implicitly as B+P/Q where we arrange
C     P .GE. 0.  The implicit form is used to prevent overflow.
C
      P = (B-A)*FB
      Q = FA - FB
      if (P.lt.ZERO) then
         P = -P
         Q = -Q
      end if
C
C  Update A; check if reduction in the size of bracketing interval is
C     satisfactory.  If not, bisect until it is.
C
      A = B
      FA = FB
      KOUNT = KOUNT + 1
      if (KOUNT.ge.4) then
         if (EIGHT*ACMB.ge.WIDTH) then
            B = B + CMB
            go to 30

         end if

         KOUNT = 0
         WIDTH = ACMB
      end if
C
C  Test for too small a change.
C
      if (P.le.ABS(Q)*TOL) then
C
C  Increment by tolerance.
C
         B = B + SIGN(TOL,CMB)
      else
C
C  Root ought to be between B and (C+B)/2.
C
         if (P.lt.CMB*Q) then
C
C  Use secant rule.
C
            B = B + P/Q
         else
C
C  Use bisection.
C
            B = B + CMB
         end if

      end if
C
C  Have completed computation for new iterate B.
C
   30 call SHOOT(B,X,MU,FB)
      NF = NF + 1
      if (ABS(FB).eq.ZERO) then
         IFLAG = 0
         C = B
         FC = FB
         return

      end if

      if (SIGN(ONE,FB).eq.SIGN(ONE,FC)) then
         C = A
         FC = FA
      end if

      go to 20

      end subroutine ZZERO

      end module SLEDGEMD
