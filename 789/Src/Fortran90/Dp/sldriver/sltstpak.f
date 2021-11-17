C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
C                                                                      *
C                    SLTSTPAK, SLDRIVER by JDP 1993-1998               *
C                          DOCUMENTATION SECTION                       *
C                         version 4.1, June 1998                       *
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C BUG NOTES etc:
C
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C REVISION HISTORY
C Spring 96. As part of implementation of SLDRIVER v3.0
C - Problem set (TSTSET routine) put in a separate file
C - GETBCS had ISING argument removed. This was 'a concession to SL02F'
C   and my instinct that it was cack-handed proved correct when it
C   was part of the reason for a subtle bug that only affected problems
C   with LC, non-Friedrichs, BCs applied at a truncated endpoint
C Spring 97. As part of implementation of SLDRIVER v4.0
C - Converted to a F90 module.
C June 98. As part of tidying up for TOMS, documentation corrected to
C   describe the variables that were formerly in COMMON and are now
C   module variables of each Problem Set module.
C   These put in a separate module SLTSTVAR, and changes made so that
C   XINFTY is a global constant rather than a variable
C
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                                                                      *
C                          GENERAL DESCRIPTION                         *
C                                                                      *
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C
C This package of 13 routines implements an interface between a Test Set
C of Sturm-Liouville Problems, (SLPs), and a Driver program which
C applies a particular Sturm-Liouville solver to the Test Set.
C
C The Test Set is entirely contained within a routine TSTSET which is in
C a separate module TESTMOD. If a Problem-writer wishes to alter the
C Test Set or replace it completely, only TSTSET needs to be altered.
C
C My naming convention is that TSTSETs live in files whose name reflects
C the contents, e.g the set of 60 problems from my book is in
C STANDARD.FOR and a small sample collection is in SAMPLE.FOR.
C
C Communication between SLTSTPAK & TESTMOD is by global (module)
C variables of TESTMOD. I chose to put them in TESTMOD rather than in a
C separate module in the hope that this will be more efficient: TSTSET
C may be called millions of times during the solution of a difficult
C SLP.
C The alternative of putting them in SLTSTPAK is not available as it
C would lead to recursive USE statements. (I think such a structure is
C OK in an OO language like Java where each test set would extend a base
C class.)
C
C The documentation uses a few LaTeX conventions, e.g.: \ starts a
C special symbol such as \lambda; underscore _ and caret ^ stand
C for subscript and superscript; { } enclose a group of symbols.
C
C DESIGN CONSIDERATIONS:
C Functional considerations.
C 1.  The class of problem catered for is the classical regular or
C     singular SLP
C         -(p(x)u')' + q(x)u = \lambda.w(x)u,      a<x<b,           (1)
C     with separated boundary conditions (BCs) as follows:
C       at a regular endpoint (say x=a)
C         BC   A1.u(a) = A2.pu'(a)          with A1,A2 constant;   (2a)
C       at a limit-point (LP) endpoint
C         though no BC is needed, a regular BC is available in case
C         one wishes to solve a truncated problem.
C       at a limit-circle (LC) singular (possibly oscillatory) endpoint
C         BC   0 = [u,h](a) = \lim_{x->a+0} [u,h](x)               (2b)
C         where [u,h](x) = u(x)ph'(x) - pu'(x)h(x);
C         and   h = A1.f + A2.g             with A1,A2 constant;
C         and   f, g are appropriate 'admissible' functions defining
C               independent BCs at the endpoint, such that, if the
C               endpoint is non-oscillatory, f specifies the
C               Friedrichs BC.
C
C 2.  The SLP may involve parameters whose values need to be set by the
C     user.
C
C 3.  The definition of a 'Problem' in the TSTSET code includes
C     the DE (1), and 'default' BCs (2a) or (2b) as appropriate, but
C     NOT
C        the values of any parameters
C        eigenvalue index k
C        tolerance TOL
C     These are to be supplied by the calling main program.
C
C 4.  Brief titling information should be provided
C     (a) to identify the selected problem
C     (b) to number and name parameters (for use in prompt for input)
C     It should be made convenient for the driver program to obtain the
C     titles of as many Problems as it wishes, e.g. to display a menu of
C     Problems.
C
C 5.  The endpoints may depend on parameters. So may their type (limit-
C     circle, regular etc) and hence the kind of BC required there.
C
C 6.  About Boundary Conditions
C 6a. The BCs may involve parameters as well as x and \lambda.
C     These therefore cannot be determined till the parameters have
C     been set.
C
C 6b. A default BC shall be provided at each endpoint that needs one
C     but this may be overridden.
C     (We support this by providing default values of the BC
C     coefficients (A1, A2 at x=a, B1, B2 at x=b) and a routine
C     to provide new values for them.)
C
C 6c. The typical automatic solution at a LC endpoint (say x=a) entails
C     replacing the given BC by a sequence of regular BCs of the form
C     [u,h](a*)=0 at points a* converging to a. This is to be supported.
C
C 7.  It is very desirable that for a problem with parameters, it be
C     possible to take a parameter other than \lambda to be the
C     eigenparameter (provided it doesn't occur in the problem in too
C     complicated a way, what this means is solver-dependent).
C
C     This has been done for the sort of parameter-dependence which 
C     the NAG code D02KEF can handle, and an interface routine COEFF2
C     for D02KEF, though not strictly part of the package, is provided 
C     with it.
C
C 8.  The typical pattern of usage of this test set is
C     (a) Select the problem, set the values of any parameters and
C         evaluate 'constants' that may depend on these parameters.
C     (b) Possibly set up BC information for later use.
C     (c) Set up information needed by the SL solver.
C     (d) Invoke the solver, which will usually perform a number of
C         integrations with trial values of \lambda, each of which
C         typically (depending on the solver used) evaluates the
C         coefficient functions and invokes the BC information.
C     (e) Report results.
C
C Non-functional considerations.
C 9.   During the solution process the coefficient functions p,q,w
C      will be evaluated many times (500 to 50000 is typical).
C      Evaluation of p,q,w should therefore be reasonably rapid.
C
C10.  Any global storage used by the package is to be invisible to the
C     calling program.
C
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                  MODULE VARIABLES & CONSTANTS                        *
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Usage:
C NPROB, XINFTY, PI
C   hold global constants (no. of problems in set, 'infinity', and
C   pi=3.14... respectively).
C   NPROB depends on the particular Problem Set being used, accordingly
C   it is set inside the particular TSTSET routine.
C   XINFTY, PI are PARAMETERs.
C
C IPROB, NPARM, PARM, A, B, ATYPE, BTYPE, SYMC
C                                , A1, A2, B1, B2
C   hold the State Information for the current Problem namely:
C   current Problem no.;
C   no. and values of Parameters for current Problem;
C   endpoints and their type-info:'regular', 'limit-circle' etc.;
C   whether Problem is symmetric (=unchanged when reflected in (a+b)/2);
C   the boundary condition coefficients (set to defaults in SETUP1 code
C   of TSTSET but alterable by call of SETBCS).
C
C NEPRM, IPARM
C   hold further state information only used by a routine like D02KEF
C   that can choose other parameters as eigenparameter.
C   NEPRM indicates params 1 to NEPRM are eligible as eigenparameter
C     in addition to standard lambda which is parameter 0.
C   IPARM is index of currently selected eigenparameter.
C     (So 0<=IPARM<=NEPRM<=NPARMA)
C
C TITLE, PARNM
C   hold Problem title and names of its Parameters
C PARM(0) (equivalenced to EIGA inside TSTSET), X, P, Q, W, PDU,
C   U are communication variables between TSTSET and other package
C   routines.
C   X transfers input data of COEFFN to corresponding sections of
C   TSTSET code.
C   P, Q, W transfer output data from TSTSET back to COEFFN.
C   PARM(0) aka EIGA, X transfer input data of GETBCS to corresponding
C   sections of TSTSET code.
C   PDU, U transfer output data from TSTSET back to GETBCS.
C NFEVAL
C   is a statistics variable. It is reset to 0 by calling NEVAL with
C   IFLAG=0 and increased by 1 at each call to COEFFN. Its value is
C   returned by the function NEVAL.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                                                                      *
C                        END OF DOCUMENTATION                          *
C                           START OF CODE                              *
C                                                                      *
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      module SLTSTPAK
      use SLTSTVAR
      use TESTMOD
      implicit none
      private
      public:: TSTINI,SETUP0,GTDAT0,SETUP1,COEFFN,GETBCS,SETBCS,NEVAL,
     +         SETIPM,SETEIG,CPU,COEFF2
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      subroutine TSTINI(TSETNMA,NPROBA)
C Global initialization routine for Test Set, to be called once per run.
C There are no input arguments.
C Output arguments:
C  Title of the test set
C  The number of Problems in the package,
C  The value to be taken as 'infinity'.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      integer NPROBA
      character*8 TSETNMA
C     .. Intrinsic Functions ..
cc      intrinsic ATAN2

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Call TSTSET to set Test-set Title and no. of Problems into global
C storage:
      call TSTSET(-1)
C Copy information into output arguments for use by Driver program:
      TSETNMA = TITLE(1:8)
      NPROBA = NPROB
      end subroutine TSTINI

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETUP0(IPROBA,TITLEA,NPARMA,NEPRMA,PARNMA)
C Problem initialization routine stage 0, to be called once per Problem.
C Input is IPROBA. The remaining arguments are output and
C their meanings are as in the instructions on coding up a Problem.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      integer IPROBA,NEPRMA,NPARMA
      character PARNMA*72,TITLEA*72

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Copy problem number into global storage:
      IPROB = IPROBA

C Call the problem-specific set-up code to set title and info about
C parameters:
      call TSTSET(0)
C Copy resulting info from global storage into output arguments:
      TITLEA = TITLE
      NPARMA = NPARM
      NEPRMA = NEPRM
      PARNMA = PARNM

      end subroutine SETUP0

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine GTDAT0(IPROBA,TITLEA,NPARMA,NEPRMA,PARNMA)
C This returns the same information as does SETUP0, but *without*
C altering the internal state of SLTSTPAK.
C Input is IPROBA. The remaining arguments are output and
C their meanings are as in the instructions on coding up a Problem.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      integer IPROBA,NEPRMA,NPARMA
      character PARNMA*72,TITLEA*72
C     .. Local Scalars ..
      integer IPRSAV

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Save the value of global variable IPROB
      IPRSAV = IPROB
C Copy problem number into global storage
      IPROB = IPROBA

C Call the problem-specific set-up code to set title and info about
C parameters:
      call TSTSET(0)
C Copy resulting info from global variables into output arguments:
      TITLEA = TITLE
      NPARMA = NPARM
      NEPRMA = NEPRM
      PARNMA = PARNM
C Restore IPROB, then restore TITLE,NPARM,NEPRM,PARNM to previous values
      IPROB = IPRSAV
      call TSTSET(0)

      end subroutine GTDAT0

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETUP1(PARMA,AA,BA,ATYPEA,BTYPEA,A1A,A2A,B1A,B2A,SYMA)
C Problem initialization routine stage 1, to be called once per Problem
C after setting PARMs if any.
C Must not be skipped even if NPARMA is 0.
C Input is PARMA. The remaining arguments are output and
C their meanings are as in the instructions on coding up a Problem.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      double precision AA,A1A,A2A,BA,B1A,B2A
      logical SYMA
      character ATYPEA*4,BTYPEA*4
C     .. Array Arguments ..
      double precision PARMA(0:10)
C     .. Local Scalars ..
      integer I

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Copy PARMA values into global storage
      do 100 I = 0,NPARM
         PARM(I) = PARMA(I)
  100 continue

C Call the 2nd section of set-up code which may depend on PARMA values
C It sets endpoints, their type etc., and may set certain constants to
C be saved inside  TSTSET, used by COEFFN and GETBCS for this problem.
      call TSTSET(1)
C Copy resulting info from global storage into output arguments:
      AA = A
      BA = B
      ATYPEA = ATYPE
      BTYPEA = BTYPE
      SYMA = SYM
      A1A = A1
      A2A = A2
      B1A = B1
      B2A = B2

      end subroutine SETUP1

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine COEFFN(XA,PA,QA,WA)
C Routine to compute coefficient functions p,q,w
C Input is XA, output is PA, QA, WA.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      double precision PA,QA,WA,XA

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Copy input argument into global storage
      X = XA
      call TSTSET(2)
C Copy output back from global storage
      PA = P
      QA = Q
      WA = W
C Increment function evaluation count
      NFEVAL = NFEVAL + 1
C Put a dot to screen every 1000 evals
      if (mod(NFEVAL,1000).eq.0) then
        write(*,fmt='(".")',advance='NO')
        if (mod(NFEVAL,50000).eq.0) write(*,*)
      end if
      end subroutine COEFFN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine GETBCS(IEND,XA,EIGA,PDUA,UA)
C Routine to evaluate boundary conditions
C Input is IEND, XA, EIGA. Output is PDUA, UA.

C GETBCS takes as input
C   IEND     0 or 1 for left or right endpoint
C   XA        current x value
C   EIGA      current lambda value, in case BC is lambda-dependent
C and produces as output
C   PDUA,UA    defining BC at x for this lambda
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      double precision PDUA,UA,XA,EIGA
      integer IEND
C     .. External Subroutines ..
cc      external SETEIG,TSTSET

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Copy input arguments into global storage before TSTSET call.
C Note SETEIG copies into current eigenparameter, i.e. PARM(IPARM):
      X = XA
      call SETEIG(EIGA)
C Call TSTSET with switch for left or right B to extract
C information in global variables; copy to output.
      if (IEND.eq.0) then
         call TSTSET(3)
      else
         call TSTSET(4)
      end if

      PDUA = PDU
      UA = U
      end subroutine GETBCS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETBCS(IEND,C1,C2)
C SETBCS overrides the given stored BC coefficients with new ones
C As in GETBCS, IEND=0 for left end, IEND=1 for right end
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      double precision C1,C2
      integer IEND

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      if (IEND.eq.0) then
         A1 = C1
         A2 = C2
      else
         B1 = C1
         B2 = C2
      end if

      end subroutine SETBCS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      integer function NEVAL(IFLAG)
C NEVAL returns the number of COEFFN calls since it was last reset.
C IFLAG is input. If IFLAG=0 then the internal counter is reset to
C zero (after returning the current value). If IFLAG is not 0 the
C counter is not reset.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      integer IFLAG

      NEVAL = NFEVAL
      if (IFLAG.eq.0) NFEVAL = 0
      end function NEVAL

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function TYPEOK(TYPE)
C TYPEOK verifies that TYPE is one of the allowed endpoint type codes
C 'R   ', 'WR  ', 'LCN ', 'LCO ', 'LPN ' or 'LPNO'.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
C     .. Scalar Arguments ..
      character TYPE*4
C     .. Intrinsic Functions ..
      intrinsic INDEX

      TYPEOK = INDEX('$R   $WR  $LCN $LCO $LPN $LPNO','$'//TYPE) .ne. 0
      end function TYPEOK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETIPM(I)
C SETIPM sets the global variable IPARM, which selects the parameter
C PARMA(IPARM) to be used as the eigenparameter, for use by the D02Kxx
C routines which can handle problems nonlinear in a parameter.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      implicit none
      integer I
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      if (I.ge.0 .and. I.le.NPARM) then
         IPARM = I
      else
         print*,'Serious SETIPM error, IPARM=',I,' out of range'
         stop
      end if
      end subroutine SETIPM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETEIG(EIGA)
C SETEIG copies its argument to the parameter in the SLTSTPAK global storage
C that is currently selected as the eigenparameter.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      implicit none
      double precision EIGA
      if (IPARM.ge.0 .and. IPARM.le.NPARM) then
         PARM(IPARM) = EIGA
      else
C        This should never happen because (I hope) IPARM is protected
C        from ever having a value outside the range
         print*,'Serious SETEIG error, IPARM=',IPARM,' out of range'
         stop
      end if
      end subroutine SETEIG

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      SUBROUTINE CPU(OPTION,CURREN,ELAPSE)
C Subroutine to initialise CPU time or return the current and
C elapsed CPU times in seconds
C A general F90 version is supplied, with a VAX/VMS F77 version in
C comments
C Input
C  OPTION - INTEGER - 0 initialises CPU time and sets
C                     elapsed to zero. Any other value returns the
C                     elapsed from the input current time and the
C                     current time.
C Input/Output
C  CURREN - DOUBLE PRECISION    - Current CPU time in seconds. On entry the time
C                     from some previous call. On exit the time at the
C                     moment with respect to zero.
C Output
C  ELAPSE - DOUBLE PRECISION    - Elapsed CPU time in seconds. The difference
C                     between the current time with respect to zero
C                     and the input current time.
C     .. Scalar Arguments ..
      INTEGER OPTION
      DOUBLE PRECISION CURREN,ELAPSE
      INTEGER ICURR, IRATE
C     .. Local Scalars ..
C(VAX)C  TIMECS - CPU time with respect to zero in centi-seconds
C(VAX)      INTEGER TIMECS
C(VAX)      INTEGER CPUTIM
C(VAX)
C(VAX)      IF (OPTION .EQ. 0) THEN
C(VAX)         CURREN = CPUTIM()/100.0
C(VAX)         ELAPSE = 0.0
C(VAX)      ELSE
C(VAX)         TIMECS = CPUTIM()
C(VAX)         ELAPSE = TIMECS / 100.0 - CURREN
C(VAX)         CURREN = TIMECS / 100.0
C(VAX)      END IF
      ELAPSE = CURREN
      call SYSTEM_CLOCK(ICURR,IRATE)
      CURREN = dble(ICURR)/dble(IRATE)
      if (OPTION.eq.0) then
        ELAPSE = 0D0
      else
        ELAPSE = CURREN - ELAPSE
      end if
      END SUBROUTINE CPU

C(VAX)      INTEGER FUNCTION CPUTIM()
C(VAX)C     .. Parameters ..
C(VAX)      INTEGER*2 NUMLST,MAXLST
C(VAX)      PARAMETER (NUMLST=1,MAXLST=NUMLST*6 + 2)
C(VAX)C     .. Local Scalars ..
C(VAX)      INTEGER*4 STATUS,FRED
C(VAX)C     .. Local Arrays ..
C(VAX)      INTEGER*2 ITMLST(MAXLST)
C(VAX)C     .. Data Statements ..
C(VAX)      DATA ITMLST / 4,'407'X,0,0,0,0,0,0 /
C(VAX)C     .. External Functions ..
C(VAX)      INTEGER*4 SYS$TRNLNM
C(VAX)      EXTERNAL SYS$TRNLNM
C(VAX)      INTEGER SYS$GETJPI
C(VAX)      INTEGER HARRY
C(VAX)      EQUIVALENCE (HARRY,ITMLST(3))
C(VAX)      SAVE
C(VAX)C     .. Executable Statements ..
C(VAX)      HARRY = %LOC(FRED)
C(VAX)      STATUS = SYS$GETJPI(,,,ITMLST,,,)
C(VAX)      IF (.NOT. STATUS) CALL LIB$STOP(%VAL(STATUS))
C(VAX)      CPUTIM = FRED
C(VAX)      END

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine COEFF2(PA, QQ, WA, XA, MU, JINT)
C COEFF2 is of the form needed by the NAG SL solver D02KEF.  It
C interfaces to SLTSTPAK in such a way that it can handle problems that
C are nonlinear in the eigenparameter, this being a facility of the
C D02K routines. The interface is a KLUGE but at least it can be done.
C
C Explanation of logic:
C Either the standard parameter lambda held in EIGA (if IPARM .eq. 0),
C or the IPARMth of the NPARMA parameters in PARMA(1:NPARMA) (if IPARM
C .gt. 0), has been nominated as eigenparameter. Temporarily NPARMA is
C called m, IPARM is called i, and the parameters are called mu(1) ..
C mu(m). It is assumed that in the standard definition 
C           -(p u')' + q u = lambda w u
C of a SLP,
C           p depends only on x;
C           q, w may depend on any of mu(1) .. mu(m) and x.
C COEFF2 must return QQ = lambda*w - q, and the partial derivative 
C WA = d(QQ)/d(mu(i)), through the arguments.
C If IPARM, that is i, = 0:
C   this is no problem since WA is just w. Also the input argument MU
C   can be ignored as far as COEFFN is concerned as the functions do not
C   depend on it. 
C Otherwise if i>0,
C   the value of mu is copied to the mu(i) in global storage currently
C   nominated as eigenparameter, by calling SETEIG; the writer of the
C   Problem Code in TSTSET is responsible for ensuring that the COEFFN
C   code sets QQ = lambda *w - q and WA = d(QQ)/d(mu(i)).

C Argument JINT is to do with D02KEF's 'breakpoint' facility. It is 
C ignored in this routine tho' could be used for Problem 53 of the
C standard test set.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      implicit none
      double precision XA,MU,PA,QQ,WA
      integer JINT
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      PARM(IPARM) = MU
      call COEFFN(XA,PA,QQ,WA)
      if (IPARM .eq. 0) QQ = MU*WA-QQ
      end subroutine COEFF2

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      end module SLTSTPAK
