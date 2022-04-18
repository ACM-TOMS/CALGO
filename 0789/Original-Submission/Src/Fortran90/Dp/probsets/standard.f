C***+****|****+****|* COPYRIGHT J D PRYCE 1997 **|****+****|****+****|**
      module TESTMOD
      use SLTSTVAR
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C STANDARD Problem Set v1.2 Apr 96 by John Pryce
C Revision History
C v1.0 with original SLTSTPAK 1994
C v1.1 Apr 96:
C - Minor changes to coefficient functions to reduce risk of overflow
C - Implement/correct almost all the BC functions for LC problems.
C   Still some that I haven't managed to solve, on LCN/LCO border!
C v1.2 Feb 97: converted to a F90 module.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine TSTSET(ISWTCH)

C   TSTSET contains the full definition of each problem.
C   The effect of a call is controlled by the module variable IPROB
C  (Problem number) and argument ISWTCH(section of code within the
C   Problem) via a number of computed GOTO statements at the top.
C   Case ISWTCH=-1: initialization call done once per run
C   Case ISWTCH=0: SETUP0 code, first setting-up stage
C   Case ISWTCH=1: SETUP1 code, second setting-up stage after reading
C                  parameter values
C   Case ISWTCH=2: COEFFN code, evaluate coefficient functions
C   Case ISWTCH=3: GETBCS(Left) code, evaluate left BC
C   Case ISWTCH=4: GETBCS(Right) code, evaluate right BC
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C   Uses module variables described below
C * NOTE. For readability the module variables are renamed inside
C   this routine(PARM becomes PARM etc); the description below uses
C   the 'standard' names.
C   The different types of call are:
C+------------------------------------------------------------------+
C|  Case ISWTCH=-1(called from TSTINI)                              |
C|  INITIALIZATION CALL done once per run:                          |
C+------------------------------------------------------------------+
C     Copies value of
C       NPROBP(the no. of Problems in the collection)
C       TSETNM(the name of the collection)
C     set in a PARAMETER statement, to module vars NPROB & TITLE, so
C     they can be transmitted to the Driver by TSTINI.
C  *  NOTE. TSTSET is the only package routine that 'knows' this
C     information. Encapsulating this data inside TSTSET makes it
C     replaceable by another TSTSET without changing *any* other
C     part of the package.
C
C+------------------------------------------------------------------+
C|  PROBLEM-SPECIFIC CALLS:                                         |
C|  Various computed GOTOs controlled by IPROB cause a jump to      |
C|  label 10*IPROB+ISWTCH, for instance calling TSTSET(3) for       |
C|  Problem 17 will jump to label 173.                              |
C+------------------------------------------------------------------+
C   Case ISWTCH=0(called from SETUP0).
C     Input is the value of IPROB, set by SETUP0.
C     The title of the problem, number of parameters and the names
C     of parameters are put in the module variables TITLE, NPARM
C     and PARNM respectively.
C   Case ISWTCH=1(called from SETUP1).
C     Input, as well as the data set up by SETUP0 call, is the
C     parameter array PARM, set by SETUP1.
C     The endpoints are put in module variables A and B.
C     Their types are put in module variables ATYPE and BTYPE.
C     Their default BC coefficients are put in module variables
C     A1, A2, B1, B2.
C     Symmetry information is put in SYM.
C   Case ISWTCH=2(called from COEFFN).
C     Input is the value of X, set by COEFFN.
C     Evaluate p,q,w at x=X and put the result in P,Q,W.
C   Case ISWTCH=3(called from GETBCS(Left), i.e. GETBCS with IEND=0).
C     Input is the value of X and PARM(0), aka EIG, set by GETBCS.
C     Evaluate the LEFT-hand BC information at x=X, lambda=PARM(0),
C     aka EIG, putting the results in module variables PDU, U.
C   Case ISWTCH=4(called from GETBCS(Right), i.e. GETBCS with IEND=1).
C     Input is the value of X and PARM(0), aka EIG, set by GETBCS.
C     Evaluate the RIGHT-hand BC information at x=X, lambda=PARM(0),
C     aka EIG, putting the results in module variables PDU, U.
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

C     .. Parameters ..
C+------------------------------------------------------------------+
C| TSETNM is name given to Test Set by its implementer.             |
C| NPROBP is total no. of Problems.                                 |
C| There must be a Problem number I for every I from 1 to NPROBP.   |
C| *NOTE* When altering the Problem Set, re-set TSETNM & NPROBP !!! |
C+------------------------------------------------------------------+
      implicit none
      integer NPROBP
      character*8 TSETNM
      parameter(TSETNM='standard', NPROBP=60)
C     .. Local Scalars ..
C     These are used in various problems. Some are given a value by the
C     TSTSET(1) call, which is held for later use, so they need to be
C     SAVEd.  Others are used as temporaries so don't need to be SAVEd:
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      double precision ALPHA,BETA,C,D,DE,DESCAL,EXPM2X,EXPMX,L,
     +       LLP1,N,NU,NUSQ,OMEGA,PDUF,PDUNF,R,RE,REX6,T,TEMP,THETA,UF,
     +       UNF,VEL
C     .. Scalar Arguments ..
      integer ISWTCH
C     .. Save statement ..
      save ALPHA,BETA,C,NU,NUSQ,N,
     +     OMEGA,T,VEL,DE,RE,L,DESCAL,LLP1,REX6
C     .. Data statements ..
      data DE/62D0/,RE/3.56D0/,L/7D0/

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Put ISWTCH.EQ.2 case first as this is the most frequently called:
      if (ISWTCH.eq.2) then
C     'COEFFN' call:
         go to(12,22,32,42,52,62,72,82,92,102,112,122,132,142,152,162,
     +          172,182,192,202,212,222,232,242,252,262,272,282,292,302,
     +          312,322,332,342,352,362,372,382,392,402,412,422,432,442,
     +          452,462,472,482,492,502,512,522,532,542,552,562,572,582,
     +          592,602) IPROB

      else if (ISWTCH.eq.-1) then
C     'TSTINI' call:
C        Copy title & total no. of Problems into module variables:
         TITLE(1:8) = TSETNM
         NPROB = NPROBP
         return

      else if (ISWTCH.eq.0) then
C     'SETUP0' call:
C     Start of new problem so set default values of IPARM, NEPRM:
         IPARM = 0
         NEPRM = 0

         go to(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     +          170,180,190,200,210,220,230,240,250,260,270,280,290,300,
     +          310,320,330,340,350,360,370,380,390,400,410,420,430,440,
     +          450,460,470,480,490,500,510,520,530,540,550,560,570,580,
     +          590,600) IPROB

      else if (ISWTCH.eq.1) then
C     'SETUP1' call:
         go to(11,21,31,41,51,61,71,81,91,101,111,121,131,141,151,161,
     +          171,181,191,201,211,221,231,241,251,261,271,281,291,301,
     +          311,321,331,341,351,361,371,381,391,401,411,421,431,441,
     +          451,461,471,481,491,501,511,521,531,541,551,561,571,581,
     +          591,601) IPROB

      else if (ISWTCH.eq.3) then
         go to(13,23,33,43,53,63,73,83,93,103,113,123,133,143,153,163,
     +          173,183,193,203,213,223,233,243,253,263,273,283,293,303,
     +          313,323,333,343,353,363,373,383,393,403,413,423,433,443,
     +          453,463,473,483,493,503,513,523,533,543,553,563,573,583,
     +          593,603) IPROB

      else if (ISWTCH.eq.4) then
         go to(14,24,34,44,54,64,74,84,94,104,114,124,134,144,154,164,
     +          174,184,194,204,214,224,234,244,254,264,274,284,294,304,
     +          314,324,334,344,354,364,374,384,394,404,414,424,434,444,
     +          454,464,474,484,494,504,514,524,534,544,554,564,574,584,
     +          594,604) IPROB

      else
C       This signals a serious coding error in the Package so:
         print *,'Invalid ISWTCH',ISWTCH,' and/or IPROB',IPROB,
     +     ' input to TSTSET call'
         stop

      end if

C***+****|****+****|     START OF PROBLEM SET    |****+****|****+****|**
C Problem #1
C This is parametrized to test the nonlinear-in-eigenparameter facility
C of NAG D02Kxx routines
   10 TITLE =  '-u" + u/(x+alpha)^2 = lambda u'
      NPARM = 1
      NEPRM = 1
      PARNM = 'alpha(>0, 0.1 gives standard problem in book)'
      go to 1000

   11 A = 0D0
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   12 ALPHA = PARM(1)
      P = 1D0
      Q = 1D0/(X+ALPHA)**2
      W = 1D0
      if (IPARM .eq. 1) then
         Q = EIG*W - Q
         W = 2d0/(X+ALPHA)**3
      end if
      go to 1002

   13 U = A2
      PDU = A1
      go to 1003

   14 U = B2
      PDU = B1
      go to 1004
C For $\lam_k$ see Appendix A

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #2
   20 TITLE = 'Mathieu equation -u" + 2r cos(2x) u = lambda u'
      NPARM = 1
      NEPRM = 1
      PARNM = 'r'
      go to 1000

   21 A = 0D0
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   22 P = 1D0
      Q = 2D0*PARM(1)*COS(2D0*X)
      W = 1D0
      if (IPARM .ne. 0) then
C       Only IPARM=1 allowed, taking PARM(1)=r as eigenparameter
         Q = PARM(0)*W - Q
         W = 4D0*SIN(2D0*X)
      end if
      go to 1002

   23 U = A2
      PDU = A1
      go to 1003

   24 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Regular problems with strongly varying coefficient behaviour
C Problem #3
   30 TITLE = 'Ref: Klotter, Technisch Schwingungslehre, I, p.12'
      NPARM = 0
      PARNM = ' '
      go to 1000

   31 A = 8D0/7D0
      B = 8D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   32 P = 1D0
      Q = 3D0/(4D0*X**2)
      W = 64*PI*PI/(9D0*X**6)
      go to 1002

   33 U = A2
      PDU = A1
      go to 1003

   34 U = B2
      PDU = B1
      go to 1004
C $a = 1$  \qquad  Regular  \qquad  $u(a) = 0$ \\
C $b = 2$  \qquad  Regular  \qquad  $u(b) = 0$
C $\lam_k =(k+1)^{2}$, $k = 0,1,\ldots$\\
C Transformation of $-d^2v/dt^2=\lam v$, $v(\pi/3)=0=v(4\pi/3)$ by
C $t={4\pi / 3x^2}$, $u=x^{3/2}v$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #4
   40 TITLE = 'Truncated Hydrogen equation'
      NPARM = 1
      PARNM = 'Right-hand endpoint b'
      go to 1000

   41 A = 0D0
      B = PARM(1)
      ATYPE = 'LPN'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   42 P = 1D0
      Q = -1D0/X + 2D0/X**2
      W = 1D0
      go to 1002

   43 U = A2
      PDU = A1
      go to 1003

   44 U = B2
      PDU = B1
      go to 1004
C As $b$ is increased some codes run into meshing difficulties, putting
C too few meshpoints near 0.
C The lower evs well approximate those of the infinite problem.
C E.g. with b=1000:
C \lambda_0=-6.250000000000E-2, \lambda_9=-2.066115702478E-3,
C \lambda_17=-2.5757359232E-4, \lambda_18=2.873901310E-5.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Regular problems with oscillatory coefficients
C Problem #5
   50 TITLE = 'Version of Mathieu equation, -u" + c cos(x) u = lambda u'
      NPARM = 1
      PARNM = 'c'
      go to 1000

   51 A = 0D0
      B = 40D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   52 P = 1D0
      Q = PARM(1)*COS(X)
      W = 1D0
      go to 1002

   53 U = A2
      PDU = A1
      go to 1003

   54 U = B2
      PDU = B1
      go to 1004
C The lower eigenvalues form clusters of 6; more and tighter clusters as
C $c$ increases.
C For $c=5$: $\lam_0=-3.484229$, $\lam_5=-3.484187$, $\lam_6=-0.599543$,
C $\lam_11=-0.595603$, $\lam_12=1.932915$, $\lam_17=1.995459$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #6
   60 TITLE = 'Truncated Gelfand-Levitan.'
      NPARM = 0
      PARNM = ' '
      go to 1000

   61 A = 0D0
      B = 100D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 1D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   62 P = 1D0
      T = 1D0 + X/2D0 + SIN(2D0*X)/4D0
      Q = 2D0*(T*SIN(2D0*X)+COS(X)**4)/T**2
      W = 1D0
      go to 1002

   63 U = A2
      PDU = A1
      go to 1003

   64 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad Regular \qquad $u(a)-(pu')(a) = 0$ \\
C $b = 100$ \qquad Regular \qquad $u(b) = 0$
C $\lam_0 = 0.00024681157  \qquad\lam_{99} = 9.77082852816$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Close-eigenvalue problems
C Problem #7
   70 TITLE =
     + 'Coffey-Evans eqn -u" + beta[beta sin(2x)^2 - 2cos(2x)]u = lam u'
      NPARM = 1
      PARNM = 'beta'
      go to 1000

   71 A = -PI/2D0
      B = -A
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      BETA = PARM(1)
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   72 P = 1D0
      Q = -2D0*BETA*COS(2.D0*X) + (BETA*SIN(2D0*X))**2
      W = 1D0
      go to 1002

   73 U = A2
      PDU = A1
      go to 1003

   74 U = B2
      PDU = B1
      go to 1004
C $a = -\pi/2 $ \qquad  Regular \qquad $u(a) = 0$ \\
C $b =  \pi/2 $ \qquad  Regular \qquad $u(b) = 0$
C As $\beta$ increases there are very close eigenvalue triplets
C $\{\lam_2,\lam_3,\lam_4\}$, $\{\lam_6,\lam_7,\lam_8\}$, \ldots, with
C the other evs well separated. $\lam_0$ is very close to zero. \\
C $(\beta=20)\lam_0 = 0.0000000000000 \qquad \lam_3 = 151.46322365766$\\
C $(\beta=30)\lam_0 = 0.0000000000000 \qquad \lam_3 = 231.66492931296$\\
C $(\beta=50)\lam_0 = 0.0000000000000 \qquad \lam_3 = 391.808191489$\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #8
   80 TITLE = 'Truncated Lennard-Jones LJ(12,6)'
      NPARM = 1
      PARNM = 'endpoint b'
      go to 1000

   81 A = 0D0
      B = PARM(1)
      ATYPE = 'LPN'
      BTYPE = 'R'
      SYM = .FALSE.
      DESCAL = 1.92D0/16.858056*DE
      LLP1 = L*(L+1D0)
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   82 P = 1D0
      REX6 =(RE/X)**6
      Q = DESCAL*REX6*(REX6-2D0) + LLP1/X**2
      W = 1D0
      go to 1002

   83 U = A2
      PDU = A1
      go to 1003

   84 U = B2
      PDU = B1
      go to 1004
C Shows close evs can happen even with highly asymmetric potentials.
C The effect is due to the resonance barrier.
C With constants as given, b=38.85 gives approximately the minimum
C splitting:
C \lam_0=0.0899594272, \lam_1=0.0899769187.
C q(x) gets large at 0: left BC u(0.8)=0 is acceptable substitute

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Regular problems that look singular
C Problem #9
   90 TITLE =
     +   'Regular with nasty w(x). Ref Fox, BVPs in DEs, Wisconsin 1960'
      NPARM = 0
      PARNM = ' '
      go to 1000

   91 A = -1D0
      B = 1D0
      ATYPE = 'WR'
      BTYPE = 'WR'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   92 P = 1D0/SQRT(1-X**2)
      Q = 0D0
      W = P
      go to 1002

   93 U = A2
      PDU = A1
      go to 1003

   94 U = B2
      PDU = B1
      go to 1004
C $a = -1$ \qquad  Regular \qquad $u(a) = 0$ \\
C $b = 1$ \qquad   Regular \qquad $u(b) = 0$
C $\lam_0 = 3.559279966  \qquad  \lam_9 = 258.8005854$
C $\lam_{24} = 1572.635284

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #10
  100 TITLE =
     +   'Regular with nasty 1/p(x). Ref: J D Pryce, Num Sol SLPs, 1993'
      NPARM = 0
      PARNM = ' '
      go to 1000

  101 A = -1D0
      B = 1D0
      ATYPE = 'WR'
      BTYPE = 'WR'
      SYM = .FALSE.
      A2 = 1D0
      A1 = 0D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  102 P = SQRT(1-X**2)
      Q = 0D0
      W = 1D0
      go to 1002

  103 U = A2
      PDU = A1
      go to 1003

  104 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty $  \qquad  \Spc: none
C $\lam_0 = 0.3856819?? \qquad \lam_9 = 1031.628??$(values uncertain)

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #11
  110 TITLE = 'Regular with nasty q(x). Ref: Pruess/Fulton 133'
      NPARM = 0
      PARNM = ' '
      go to 1000

  111 A = 0D0
      B = 4D0
      ATYPE = 'WR'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  112 P = 1D0
      Q = LOG(X)
      W = 1D0
      go to 1002

  113 U = A2
      PDU = A1
      go to 1003

  114 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad  Regular \qquad  $u(a) = 0$
C $b = 4$ \qquad  Regular \qquad  $u(b) = 0$
C $ \lambda _0 = 1.1248168097 \qquad  \lambda _{24} = 385.92821596$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Parametrized singular problems
C Problem #12
  120 TITLE = 'Bessel`s equation -(xu'')'' + (nu**2/x) u = lambda x u'
      NPARM = 1
      PARNM = 'nu**2(can be <0)'
      go to 1000

  121 A = 0D0
      B = 1D0
      NUSQ = PARM(1)
C Note usual NU is actually i*(this NU) when NUSQ<0:
      NU = SQRT(ABS(NUSQ))

      if (NUSQ.ge.1D0) then
         ATYPE = 'LPN'
      else if (NUSQ.ge.0D0) then
         ATYPE = 'LCN'
      else
         ATYPE = 'LCO'
      end if

      BTYPE = 'R'
      SYM = .FALSE.
C in all cases:
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

C Remaining bits assume NUSQ & NU have been set:
  122 P = X
      Q = NUSQ/X
      W = X
      go to 1002

  123 if (NUSQ.ge.1D0) then
         U = A2
         PDU = A1

      else if (NUSQ.gt.0D0) then
         U = A1*(X**NU) + A2/(X**NU)
         PDU = NU*A1*(X**NU) - NU*A2/(X**NU)

      else if (NUSQ.eq.0D0) then
         U = A1 + A2*LOG(X)
         PDU = A2

      else
         U = A1*COS(NU*LOG(X)) + A2*SIN(NU*LOG(X))
         PDU = NU*(-A1*SIN(NU*LOG(X))+A2*COS(NU*LOG(X)))
      end if

      go to 1003

  124 U = B2
      PDU = B1
      go to 1004
C $ a = 0 $ \qquad LcN if $0\leq\nu<1$, LpN if $\nu\geq1$
C \qquad Fried. bc: $\nu u(a)-(pu')(a) = 0$ \\
C $ b = 1 $ \qquad Regular \qquad  $u(b) = 0$
C Number of eigenvalues: $\infty $  \qquad  \Spc: none \qquad
C $\nu={1\over3}$: $\lam_0=8.4250069295$ \qquad $\lam_9=970.7184494$\\
C $\nu=2$: $\lam_0 = 26.374616427  \qquad  \lam_9 = 1136.8036878 $\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #13
  130 TITLE =
     +    'Bessel eqn in normal form. -u" + ((nu**2-1/4)/x**2)u = lam u'
      NPARM = 1
      PARNM = 'nu**2(can be <0)'
      go to 1000

  131 NUSQ = PARM(1)
      A = 0D0
      B = 1D0
C     settng NU for NUSQ < 0(as well as >= 0) is useful for LCO case:
      NU = SQRT(ABS(NUSQ))
      if (NUSQ.ge.1D0) then
         ATYPE = 'LPN'
C special case which doesn't occur in the non-normal form
      else if (NUSQ.eq.0.25D0) then
         ATYPE = 'R'

      else if (NUSQ.ge.0D0) then
         ATYPE = 'LCN'

      else
         ATYPE = 'LCO'
      end if
C In all cases take:
      A2 = 0D0
      A1 = 1D0

      BTYPE = 'R'
      SYM = .FALSE.
      B2 = 0D0
      B1 = 1D0
      go to 1001

C Remaining bits assume NUSQ & NU have been set:
  132 P = 1D0
      Q = 0D0
      if (NUSQ.ne.0.25D0) Q =(NUSQ-0.25D0)/(X*X)
      W = 1D0
      go to 1002
C BCs  handle all cases:
  133 if (NUSQ.ge.1D0) then
C        LPN case
         U = A2
         PDU = A1

      else if (NUSQ.eq.0.25D0) then
C         Regular case
         U = A1*X + A2
         PDU = A1

      else if (NUSQ.ge.0D0) then
C        LCN case, Fried. BC function=x**(1/2+nu), other=x**(1/2-nu)
         UF = X**(0.5D0+NU)
         UNF = X**(0.5D0-NU)

         PDUF =(0.5D0+NU)*UF/X
         PDUNF =(0.5D0-NU)*UNF/X

         U = A1*UF + A2*UNF
         PDU = A1*PDUF + A2*PDUNF
c      print*,'TSTSET#13:X,A1,A2,UF,UNF,PDUF,PDUNF,U,PDU',
c     + X,A1,A2,UF,UNF,PDUF,PDUNF,U,PDU
      else
C        LCO case, BC function RealPart(x**(1/2+i*nu)(a1-i*a2))
         TEMP = NU*LOG(X)
         U =(A1*COS(TEMP)+A2*SIN(TEMP))*SQRT(X)
         PDU =((0.5D0*A1+NU*A2)*COS(TEMP)-(NU*A1-0.5D0*A2)*SIN(TEMP))/
     +         SQRT(X)
      end if

      go to 1003

  134 U = B2
      PDU = B1
      go to 1004
C Same cases as Problem 12 except for extra Regular case.
C BC functions are sqrt(x) times those of Problem 12.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #14
  140 TITLE =
     +       'Ref: Bender & Orzsag. -u" - l(l+1)sech**2(x) u = lambda u'
      NPARM = 1
      PARNM = 'l'
      go to 1000

  141 L = PARM(1)
      A = -XINFTY
      B = XINFTY
      ATYPE = 'LPNO'
      BTYPE = 'LPNO'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  142 P = 1D0
C Compute sech^2 x without overflow risk:
      EXPM2X = EXP(-2D0*ABS(X))
      Q = -L*(L+1D0) * 4D0*EXPM2X/(1D0+EXPM2X)**2
      W = 1D0
      go to 1002

  143 U = A2
      PDU = A1
      go to 1003

  144 U = B2
      PDU = B1
      go to 1004
C $a = -\infty $  \qquad LpN/O \\
C $b = +\infty $ \qquad  LpN/O \\
C Number of eigenvalues: $l$ \qquad \Spc:(0,$\infty $)  \qquad
C $ \lam_k = -(l+1-k)^{2}$, $k = 1,2,\ldots,l$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #15
  150 TITLE =
     +   'Assoc. Legendre eqn. -((1-x^2)u'')'' + c u/(1-x^2) = lambda u'
      NPARM = 1
      PARNM = 'c(c=0 for usual Legendre, c=1/4 for Chebyshev)'
      go to 1000

  151 A = -1D0
      B = 1D0
      C = PARM(1)
C See note re NU in Prob 12:
      NU = SQRT(ABS(C))

      if (C.ge.1D0) then
         ATYPE = 'LPN'
         BTYPE = 'LPN'

      else if (C.ge.0D0) then
         ATYPE = 'LCN'
         BTYPE = 'LCN'

      else
         ATYPE = 'LCO'
         BTYPE = 'LCO'
      end if
C In all cases:
      A1 = 1D0
      A2 = 0D0
      B1 = 1D0
      B2 = 0D0
      SYM = .TRUE.
      go to 1001

C Remaining bits assume C & NU have been set:
  152 P = 1D0 - X**2
      Q = C/(1D0-X**2)
      W = 1D0
      go to 1002

  153 if (C.ge.1D0) then
C        LPN case:
         U = A2
         PDU = A1

      else if (C.gt.0D0) then
C        LCN case:
         TEMP =(1D0-X**2)**(NU*0.5D0)
         U = A1*TEMP + A2/TEMP
         PDU = NU*X*(-A1*TEMP+A2/TEMP)

      else if (C.eq.0D0) then
C        Also LCN:
         U = A1 + A2*LOG(1-X**2)
         PDU = -A2*2D0*X

      else
C        LCO case:
         TEMP = NU*0.5D0*LOG(1-X**2)
         U = A1*COS(TEMP) + A2*SIN(TEMP)
         PDU = NU*X*(A1*SIN(TEMP)-A2*COS(TEMP))
      end if

      go to 1003
C Right end is same as left end except A1 A2 become B1 B2:
  154 if (C.gt.1D0) then
C         LPN case:
         U = B2
         PDU = B1

      else if (C.gt.0D0) then
C         LCN case:
         TEMP =(1D0-X**2)**(NU*0.5D0)
         U = B1*TEMP + B2/TEMP
         PDU = NU*X*(-B1*TEMP+B2/TEMP)

      else if (C.eq.0D0) then
C         Also LCN:
         U = B1 + B2*LOG(1-X**2)
         PDU = -B2*2D0*X

      else
C         LCO case:
         TEMP = NU*0.5D0*LOG(1-X**2)
         U = B1*COS(TEMP) + B2*SIN(TEMP)
         PDU = NU*X*(B1*SIN(TEMP)-B2*COS(TEMP))
      end if

      go to 1004
C Number of eigenvalues: $\infty $ \qquad     \Spc: none \\
C $ \lambda _k =(k+\nu+1)(k+\nu),  k = 0,1,...$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #16
  160 TITLE =
     +'Assoc Legendre tfmed. -u" + [(c-1/4)sec(x)^2 - 1/4] u = lambda u'
      NPARM = 1
      PARNM =
     +       'c(c=0 for usual Legendre, c=1/4 reduces to regular case)'
      go to 1000

  161 A = -PI/2D0
      B = PI/2D0
      C = PARM(1)
      NU = SQRT(ABS(C))
      if (C.ge.1D0) then
         ATYPE = 'LPN'
         BTYPE = 'LPN'

      else if (C.eq.0.25D0) then
         ATYPE = 'R'
         BTYPE = 'R'

      else if (C.ge.0D0) then
         ATYPE = 'LCN'
         BTYPE = 'LCN'

      else
         ATYPE = 'LCO'
         BTYPE = 'LCO'
      end if
C In all cases:
      A1 = 1D0
      A2 = 0D0
      B1 = 1D0
      B2 = 0D0
      SYM = .TRUE.
      go to 1001

C Remaining bits assume C & NU have been set:
  162 P = 1D0
      Q =(C-0.25D0)/COS(X)**2 - 0.25D0
      W = 1D0
      go to 1002

  163 if (C.gt.1D0) then
C         LPN case:
         U = A2
         PDU = A1

      else if (C.gt.0D0) then
C         LCN case:
         TEMP = COS(X)**NU
         U =(A1*TEMP+A2/TEMP)*SQRT(COS(X))
         PDU = -((0.5D0+NU)*A1*TEMP+ (0.5D0-NU)*A2/TEMP)*SIN(X)/
     +         SQRT(COS(X))

      else if (C.eq.0D0) then
      if(cos(x).lt.1d-7) print*,'x,x-a,cos x=',x,x-a,cos(x)
C         Also LCN:
         U =(A1+A2*LOG(COS(X)))*SQRT(COS(X))
         PDU = -(A1+A2*(1D0+0.5D0*LOG(COS(X))))*SIN(X)/SQRT(COS(X))

      else
C         LCO case:
         TEMP = NU*LOG(COS(X))
         U =(A1*COS(TEMP)+A2*SIN(TEMP))*SQRT(COS(X))
         PDU =((-0.5D0*COS(TEMP)+NU*SIN(TEMP))*A1+
     +        (-0.5D0*SIN(TEMP)-NU*COS(TEMP))*A2)*SIN(X)/SQRT(COS(X))
      end if

      go to 1003
C Right end is same as left end except A1 A2 become B1 B2:
  164 if (C.gt.1D0) then
C         LPN case:
         U = B2
         PDU = B1

      else if (C.gt.0D0) then
C         LCN case:
         TEMP = COS(X)**NU
         U =(B1*TEMP+B2/TEMP)*SQRT(COS(X))
         PDU = -((0.5D0+NU)*B1*TEMP+ (0.5D0-NU)*B2/TEMP)*SIN(X)/
     +         SQRT(COS(X))

      else if (C.eq.0D0) then
      if(cos(x).lt.1d-7) print*,'x,b-x,cos x=',x,b-x,cos(x)
C         Also LCN:
         U =(B1+B2*LOG(COS(X)))*SQRT(COS(X))
         PDU = -(B1+B2*(1D0+0.5D0*LOG(COS(X))))*SIN(X)/SQRT(COS(X))

      else
C         LCO case:
         TEMP = NU*LOG(COS(X))
         U =(B1*COS(TEMP)+B2*SIN(TEMP))*SQRT(COS(X))
         PDU =((-0.5D0*COS(TEMP)+NU*SIN(TEMP))*B1+
     +        (-0.5D0*SIN(TEMP)-NU*COS(TEMP))*B2)*SIN(X)/SQRT(COS(X))
      end if

      go to 1004
C Spectrum as for previous problem

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #17
  170 TITLE =
     + 'Anharmonic oscillator(Marletta PhD) -u" + x^alpha u = lambda u'
      NPARM = 1
      PARNM = 'alpha'
      go to 1000

  171 ALPHA = PARM(1)
      A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPN'
      A1 = 1D0
      A2 = 0D0
      B1 = 1D0
      B2 = 0D0
      SYM = .FALSE.
      go to 1001

  172 P = 1D0
      Q = X**ALPHA
      W = 1D0
      go to 1002

  173 U = A2
      PDU = A1
      go to 1003

  174 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  Regular  \qquad $u(a) = 0$ \\
C $b = +\infty $ \qquad  LpN  \\
C Number of eigenvalues: $\infty $ \qquad    \Spc: none
C $\alpha=2$: $\lam_k = 4k+3$, $k = 0,1,2,\ldots$
C $\alpha=3$: $\lam_0 = 3.4505626899 \qquad  \lam_{24} = 228.52088139$
C $\alpha=4$: $\lam_0 = 3.7996730298 \qquad  \lam_{24} = 397.14132678$
C $\alpha=5$: $\lam_0 = 4.0891593149  \qquad  \lam_{24} = 588.17824969$
C As $\alpha$ increases, and for large $k$, choice of mesh at the RH
C truncation point seems to become more difficult for Pruess methods
C and stiffness makes initial-value methods expensive.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Finite interval singular: LPN, and LCN Friedrichs BC
C Problem #18
  180 TITLE = 'Bessel, order 0. -(xu'')'' = lambda x u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  181 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'R'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .FALSE.
      go to 1001

  182 P = X
      Q = 0D0
      W = X
      go to 1002

  183 U = A1 + A2*LOG(X)
      PDU = A2
      go to 1003

  184 U = B2
      PDU = B1
      go to 1004
C $\nu=0$: $\lam_0 = 5.78318596295 \qquad  \lam_{19} = 3850.01252885 $\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #19
  190 TITLE = 'Bessel, order 1/2. -(xu'')'' + 1/(4x) u = lambda x u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  191 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'R'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .FALSE.
      go to 1001

  192 P = X
      Q = 0.25D0/X
      W = X
      go to 1002

  193 U = A1*sqrt(X) + A2/sqrt(X)
      PDU = 0.5D0*(A1*sqrt(X) - A2/sqrt(X))
      go to 1003

  194 U = B2
      PDU = B1
      go to 1004
C $\lam_k=((k+1)\pi)^2$, this is $-v''=\lam v$ transformed
C by $v=x^{\half}u$.\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #20
  200 TITLE = 'Legendre eqn. -((1-x^2)u'')'' = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  201 A = -1D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'LCN'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .TRUE.
C Set C, NU for subsequent invocation of Problem 15
      C = 0D0
      NU = 0D0
      go to 1001

  202 P = 1D0 - X**2
      Q = 0D0
      W = 1D0
      go to 1002

  203 U = A1 + A2*LOG(1-X**2)
      PDU = -A2*2D0*X

      go to 1003
C Right end is same as left end except A1 A2 become B1 B2:
  204 U = B1 + B2*LOG(1-X**2)
      PDU = -B2*2D0*X

      go to 1004

C Number of eigenvalues: $\infty $ \qquad     \Spc: none \\
C $ \lambda _k =(k+1)k,  k = 0,1,...$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #21
  210 TITLE = 'Slightly tricky q(x). -u" -(10/x^1.5) u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  211 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  212 P = 1D0
      Q = -10D0/(X*SQRT(X))
      W = 1D0
      go to 1002
C The BC functions can be defined in terms of Bessel fns.
C However the following is a more elementary way to define them.
C Non-Friedrichs BC was provided by Prof Patrick Parks formerly of RMCS
C then at Oxford Centre for Industrial & Applied Mathematics(OCIAM),
C Oxford, UK. Sadly died Jan 1994
  213 UF = X -(40D0/3D0)*X*SQRT(X)
      PDUF = 1D0 - 20D0*SQRT(X)
      UNF = LOG(X)*UF - 0.0025D0 - 0.1D0*SQRT(X) + 320.D0/9.D0*X*SQRT(X)
      PDUNF = LOG(X)*PDUF + 1.D0 + 40.D0*SQRT(X) - 1.D0/(20.D0*SQRT(X))

      U = A1*UF + A2*UNF
      PDU = A1*PDUF + A2*PDUNF
      go to 1003

  214 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty $ \qquad     \Spc: none \\
C With f, g being UF, UNF above:
C Though f ~ const.x, g ~ const, neither 1 nor x is an admissible
C function! Check-calculation: with Lu := -u''-(10/x^1.5) u we have
C Lf = 40/3, Lg =(40/3)ln(x) - 3200/9
C so f, g, Lf, Lg *are* square-integrable, i.e. f, g are admissible.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #22
  220 TITLE =
     +  'Legendre eqn, Liouville form. -u" -(1/4)sec(x)^2 u = lambda u'
C but note Prob 16 can't be invoked for coeff functions
C as eigenvalue is shifted by 0.25(could invoke it for BCs though)
      NPARM = 0
      PARNM = ' '
      go to 1000

  221 A = -PI/2D0
      B = PI/2D0
      ATYPE = 'LCN'
      BTYPE = 'LCN'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .TRUE.
      go to 1001

  222 P = 1D0
      Q = -0.25D0/(COS(X)**2)
      W = 1D0
      go to 1002

  223 U =(A1+A2*LOG(COS(X)))*SQRT(COS(X))
      PDU = -(A1+A2*(1D0+0.5D0*LOG(COS(X))))*SIN(X)/SQRT(COS(X))
      go to 1003

  224 U =(B1+B2*LOG(COS(X)))*SQRT(COS(X))
      PDU = -(B1+B2*(1D0+0.5D0*LOG(COS(X))))*SIN(X)/SQRT(COS(X))
      go to 1004
C Number of eigenvalues: $\infty $  \qquad \Spc: none \qquad
C $\lambda _k = k(k+1)+0.25$, $k = 0,1,\ldots $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #23
  230 TITLE =
     +'Generalized hypergeom., -u"+ (-242ch x + 241)/(4sh(x)^2)u=lam u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  231 A = 0D0
      B = XINFTY
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .FALSE.
      go to 1001

  232 P = 1D0
C      Q =(-242D0*COSH(X)+241D0)/(4D0*SINH(X)**2)
      if (X.le.1D0) then
        Q =(-121D0*SINH(0.5D0*X)**2 - 0.25D0) /(SINH(X)**2)
      else
        EXPMX = EXP(-X)
        Q =(-121D0*(1D0-EXPMX)**2-EXPMX)*EXPMX /(1D0-EXPMX**2)**2
      end if
      W = 1D0
      go to 1002

  233 U = SQRT(X)*(A1+A2*LOG(X))
      PDU =(A1*0.5D0+A2*(0.5D0*LOG(X)+1D0))/SQRT(X)
      go to 1003

  234 U = B2
      PDU = B1
      go to 1004
C $ \lambda_k = -(5-k)^2, k=0,1,2,3,4$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #24
  240 TITLE = 'Latzko equation. -((1-x^7)u'')'' = lambda x^7 u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  241 A = 0D0
      B = 1D0
      ATYPE = 'R'
      BTYPE = 'LCN'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  242 P = 1D0 - X**7
      Q = 0D0
      W = X**7
      go to 1002

  243 U = A2
      PDU = A1
      go to 1003

  244 U = B1 - B2*LOG(1D0-X)
      PDU = B2*( (((((X+1D0)*X+1D0)*X+1D0)*X+1D0)*X+1D0)*X+1D0 )
      go to 1004
C $a = 0$  \qquad Regular  \qquad $u(a) = 0$ \\
C $b = 1$  \qquad LcN  \qquad BC fns f=1, g=ln(1-x)\\
C Number of eigenvalues: $\infty $ \qquad      \Spc: none  \qquad
C $\lam_0 = 8.7274703526  \qquad \lam_2 = 435.06333218$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #25
  250 TITLE = 'Transformed regular -(x^4 u'')'' - 2x^2 u = lambda x^4 u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  251 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  252 P = X**4
      Q = -2D0*X**2
      W = P
      go to 1002

  253 U = A1/X + A2/(X**2)
      PDU = -A1*X**2 - 2D0*A2*X
      go to 1003

  254 U = B2
      PDU = B1
      go to 1004
C This is -v" = lambda v transformed by v = x^2 u.
C The p=x^4 may cause difficulty in choosing mesh at x=0
C $\lam_k =((k+1)\pi)^2, k=0,1,...$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #26
  260 TITLE =
     +   'Mysterious exact lam_0=7. -(x^3 u'')'' + x^3 u = lambda x^2 u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  261 A = 0D0
      B = 1D0
      ATYPE = 'LPN'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 1D0
      A1 = 0D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  262 P = X**3
      Q = P
      W = X**2
      go to 1002

  263 U = A2
      PDU = A1
      go to 1003

  264 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty $  \qquad  \Spc: none \qquad
C $\lam_0 = 7.0000000000  \qquad  \lam_9 = 284.53608972 $\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Infinite interval singular: LPN, and LCN Friedrichs BC
C Problem #27
  270 TITLE = 'Airy equation. -u" + x u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  271 A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPN'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  272 P = 1D0
      Q = X
      W = 1D0
      go to 1002

  273 U = A2
      PDU = A1
      go to 1003

  274 U = B2
      PDU = B1
      go to 1004
C $ a = 0  $ \qquad Regular \qquad  $u(a) = 0$ \\
C $ b = +\infty $  \qquad  LpN
C Number of eigenvalues: $\infty $ \qquad  \Spc: none
C Eigenvalues are the zeros of Airy function $\mbox{Ai}(\lam)=
C(J_{1/3}+J_{-1/3})({2 / 3}\lam^{1/3})$.\\
C $\lam_0 = 2.3381074104 \qquad \lam_9 = 12.828776753$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #28
  280 TITLE = 'Harmonic oscillator. -u" + x^2 u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  281 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  282 P = 1D0
      Q = X**2
      W = 1D0
      go to 1002

  283 U = A2
      PDU = A1
      go to 1003

  284 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty$  \qquad \Spc: none \qquad
C $\lam_k = 2k+1$, $k = 0,1,\ldots  $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #29
  290 TITLE = 'Hydrogen atom. -u" + (2/x^2 - 1/x) u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  291 A = 0D0
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  292 P = 1D0
      Q = -1/X + 2/X**2
      W = 1D0
      go to 1002

  293 U = A2
      PDU = A1
      go to 1003

  294 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty $  \qquad \Spc:(0,$\infty $) \qquad
C $\lam_k = -1/(2k+4)^2$, $k = 0,1,\dots$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #30
  300 TITLE = 'Coulomb potential. -u" - u/x = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  301 A = 0D0
      B = XINFTY
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  302 P = 1D0
      Q = -1D0/X
      W = 1D0
      go to 1002

  303 U = A1*X + A2*(1D0-X*LOG(X))
      PDU = A1 + A2*(-1D0-LOG(X))
      go to 1003

  304 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  LcN  \\
C $b = +\infty \qquad  $ LpN/O
C Number of eigenvalues: $\infty $.  \qquad  \Spc:(0,$\infty $)\\
C $ \lam_k = -1/[4(k+1)^2]$, $k = 0,1,\ldots $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #31
  310 TITLE =
     +  'Transformed H-atom. -(x^2 u'')'' -(l(l+1)-x) u = lambda x^2 u'
      NPARM = 1
      PARNM = 'l(integer>=0)'
      go to 1000

  311 A = 0D0
      B = XINFTY
      L = PARM(1)
      if (L.lt.1D0) then
         ATYPE = 'LCN'
      else
         ATYPE = 'LPN'
      end if

      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 1D0
      A1 = 0D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  312 P = X*X
      Q = L*(L+1D0) - X
      W = P
      go to 1002
C I haven't checked what happens for l that are not integer >=0
C BEWARE!
  313 if (L.eq.0D0) then
         U = A1 + A2*(1D0/X-LOG(X))
         PDU = A2*(-1D0-X)

      else
         U = A2
         PDU = A1
      end if

      go to 1003

  314 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $\infty $.  \qquad  \Spc:(0,$\infty $)\\
C $ \lam_k = -1/[4(k+l+1)^2]$, $k = 0,1,\ldots $
C This is Hydrogen atom eqn -v" + (l(l+1)/x^2 - 1/x) v = lambda v
C transformed by u = x v.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #32
  320 TITLE = 'Laguerre eqn. -u" + (x^2 + (alpha-1/4)/x^2) u = lambda u'
      NPARM = 1
      PARNM = 'alpha(=1 for standard problem)'
      go to 1000

  321 A = 0D0
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  322 P = 1D0
      Q = X**2 + (PARM(1)-0.25)/(X**2)
      W = 1D0
      go to 1002

  323 U = A2
      PDU = A1
      go to 1003

  324 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad  LpN\qquad Fried. bc: $u(a) = 0$ \\
C $b = +\infty $ \qquad LpN \qquad Fried. bc: $u(b) = 0$
C Number of eigenvalues: $\infty
C $  \qquad      \Spc: none \\
C When $\alpha=1$, $ \lambda _k = 4(k+1),  k = 0,1,...$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #33
  330 TITLE =
     + 'Raman Scattering. -u" + (-7x^2 + 0.5x^3 + x^4) u = 0.5 lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  331 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .FALSE.
      go to 1001

  332 P = 1D0
      Q = X*X*(-7D0+X*(0.5D0+X))
      W = 0.5D0
      go to 1002

  333 U = 0D0
      PDU = 1D0
      go to 1003

  334 U = 0D0
      PDU = 1D0
      go to 1004
C $a = -\infty $ \qquad LpN  \\
C $b =  \infty $ \qquad LpN
C Number of eigenvalues: $\infty $ \qquad    \Spc: none
C $\lam_0 = -24.5175977072  \qquad \lam_5 = 8.10470769427 $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #34
  340 TITLE =
     +'Pruess LCN/LCO border. -u"/2 -(1/(8x^2)+sech(x)^2) u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  341 A = 0D0
      B = XINFTY
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  342 P = 0.5D0
C  Compute sech^2 x without overflow risk:
      EXPM2X = EXP(-2D0*ABS(X))
      Q = -0.125D0/X**2 - 4D0*EXPM2X/(1D0+EXPM2X)**2
      W = 1D0
      go to 1002

  343 U =(A1+A2*LOG(X))*SQRT(X)
      PDU = 0.5D0*(A1*0.5D0+A2*(0.5D0*LOG(X)+1D0))/SQRT(X)
      go to 1003

  344 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: 1  \qquad \Spc:(0,$\infty $) \\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Typical problems from chemical physics literature
C Problem #35
  350 TITLE =
     +   'Morse(1929) potential. -u" + (9e^(-x)-18e^(-2x))u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  351 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  352 P = 1D0
C     Put a 'ceiling' on T hopefully to stop incautious codes
C     from overflowing with x<<0
C      T = EXP(-X)
      T = EXP(min(22D0,-X))
      Q = 9D0*T*(T-2D0)
      W = 1D0
      go to 1002

  353 U = A2
      PDU = A1
      go to 1003

  354 U = B2
      PDU = B1
      go to 1004
C $ a = -\infty $ \qquad LpN  \\
C $ b = +\infty $ \qquad LpN/O \\
C Number of eigenvalues: 3  \qquad  \Spc:(0,$\infty $) \qquad
C $\lam_k = -0.25-(3-k)(2-k)$, $k = 0,1,2 $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #36
  360 TITLE = 'Morse potential, Secrest et al.(1962)'
      NPARM = 0
      PARNM = ' '
      go to 1000

  361 A = 0D0
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  362 P = 1D0
      T = EXP(-1.7D0*(X-1.3D0))
      Q = 2D0/(X**2) + 2000D0*T*(T-2D0)
      W = 1D0
      go to 1002

  363 U = A2
      PDU = A1
      go to 1003

  364 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  LpN  \\
C $b = +\infty $ \qquad LpN/O \\
C Number of eigenvalues: 26 \qquad \Spc:(0,$\infty$)
C $\lam_0 = -1923.529653 \qquad \lam_1 = -177.2908125 \qquad
C \lambda _{13}  = -473.29709743  \qquad\lam_{25} = -1.7670126867 $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #37
  370 TITLE =
     +   'Quartic anharmonic oscillator. -u" + (x^4 + x^2) u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  371 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  372 P = 1D0
      Q = X**2 + X**4
      W = 1D0
      go to 1002

  373 U = A2
      PDU = A1
      go to 1003

  374 U = B2
      PDU = B1
      go to 1004
C $a = -\infty $ \qquad LpN  \\
C $b =  \infty $ \qquad LpN  \\
C Number of eigenvalues: $\infty $  \qquad  \Spc: none
C $\lam_0 = 1.3923516415 \qquad \lam_9 = 46.965069501 $

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #38
  380 TITLE = 'Close-evs problem -u" + (x^4-25x^2) u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  381 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  382 P = 1D0
      Q =(X**2-25D0)*X**2
      W = 1D0
      go to 1002

  383 U = A2
      PDU = A1
      go to 1003

  384 U = B2
      PDU = B1
      go to 1004
C $a = -\infty $ \qquad LpN \qquad Fried. bc: $u(a)= 0$ \\
C $b = +\infty $ \qquad LpN \qquad Fried. bc: $u(b) = 0$
C Number of eigenvalues: $\infty $  \qquad      \Spc: none\\
C $ \lambda _0 = -149.2194561  \qquad \lambda _1 = -149.2194561 \qquad
C \lambda {_40} = 75.69072485$\\
C Double-well version of quartic anharmonic oscillator.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #39
  390 TITLE =
     + 'Morse(Marletta). -u" + 8000[e^(-3x) - 2e^(-3x/2)] u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  391 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  392 P = 1D0
      T = EXP(-1.5D0*X)
      Q = 8D3*T*(T-2D0)
      W = 1D0
      go to 1002

  393 U = A2
      PDU = A1
      go to 1003

  394 U = B2
      PDU = B1
      go to 1004
C $ a = -\infty $ \qquad LpN  \\
C $ b = +\infty $ \qquad LpN/O \\
C Number of eigenvalues: 60 \qquad  \Spc:(0,$\infty $)
C With this deep well a large truncated interval seems to be
C needed to give good approximations to higher \ev s.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #40
  400 TITLE = 'Wicke and Harris(1976), spike at bottom of well'
      NPARM = 0
      PARNM = ' '
      go to 1000

  401 A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  402 P = 1D0
      Q = 1250D0*EXP(-83.363D0*(X-2.47826D0)**2) +
     +    3906.25D0*(1D0-EXP(2.3237D0-X))**2
      W = 1D0
      go to 1002

  403 U = A2
      PDU = A1
      go to 1003

  404 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad regular \qquad $u(a) = 0$ \\
C $b = +\infty $\qquad  LpN/O \\
C Number of eigenvalues: 62  \qquad    \Spc:(3906.25,$\infty$)\\
C $ \lambda _0 = 163.2238021    \qquad    \lambda _9 = 1277.536963 $\\
C This has a spike at the bottom of the well: whether this causes any
C trouble to automatic codes is doubtful.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #41
  410 TITLE = 'Woods-Saxon potential(Vanden Berghe et al. 1989).'
      NPARM = 1
      PARNM = 'l(>=0)'
      go to 1000

  411 A = 0D0
      B = XINFTY
      L = PARM(1)
      LLP1 = L*(L+1D0)
      if (L.eq.0.D0) then
         ATYPE = 'R'
      else
         ATYPE = 'LPN'
      end if

      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  412 P = 1D0
C T is the reciprocal of t as given in SL Book
      T = EXP(-(X-7D0)/0.6D0)
      Q = LLP1/(X*X) - 50D0*T/(T+1D0)*(1D0-(5D0/3D0)/(T+1D0))
      W = 1D0
      go to 1002

  413 U = A2
      PDU = A1
      go to 1003

  414 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad Regular($l=0$), LpN($l=1$) \qquad  $u(a) = 0$ \\
C $b = +\infty $ \qquad LpN/O \\
C($l=0$):\\
C Number of eigenvalues: 14  \qquad \Spc:(0,$\infty $) \\
C $ \lambda _0 = -49.457788728 \qquad \lambda _{10} = -18.094688282$ \\
C($l=1$):\\
C Number of eigenvalues: 13  \qquad \Spc:(0,$\infty $) \\
C $ \lambda _0 = -48.349481052  \qquad  \lambda _{10} = -13.522303353$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #42
  420 TITLE = 'Another model potential(Vanden Berghe et al. 1989).'
      NPARM = 1
      PARNM = 'l'
      go to 1000

  421 A = 0D0
      B = XINFTY
      L = PARM(1)
      if (L.eq.0D0) then
         ATYPE = 'LCN'
      else
         ATYPE = 'LPN'
      end if

      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  422 P = 1D0
      Q = L*(L+1D0)/X**2 + (-1D0+5D0*EXP(-2D0*X))/X
      W = 1D0
      go to 1002

  423 if (L.eq.0D0) then
         U = A1*X + A2*(1D0+4D0*X*LOG(X))
         PDU = A1 + A2*(4D0*LOG(X)+4D0)
      else
         U = A2
         PDU = A1
      end if

      go to 1003

  424 U = B2
      PDU = B1
      go to 1004
C $a=0$ \qquad LcN($l=0$), LpN($l=1$) \qquad Fried. bc: $u(a)=0$\\
C $b = +\infty $  \qquad  LpN/O \\
C Number of eigenvalues: $\infty $  \qquad  \Spc:(0,$\infty $) \\
C($l=0$):\\
C $\lambda_0=-0.156358880971 \qquad \lambda _2 = -0.023484895664$\\
C($l=0$):\\
C $ \lambda _0 = -0.061681846633 \qquad \lambda _2 = -0.015501561691$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Borderline Problems
C Problem #43
  430 TITLE = 'Bessel, LCN/LCO border. -u" + ((alpha-1/4)/x^2)u = lam u'
      NPARM = 1
      PARNM = 'alpha, say in range -0.1 to 0.1'
      go to 1000
C This is just Problem 13 so code is an exact copy (apart from labels)
  431 NUSQ = PARM(1)
      A = 0D0
      B = 1D0
C     settng NU for NUSQ < 0(as well as >= 0) is useful for LCO case:
      NU = SQRT(ABS(NUSQ))
      if (NUSQ.ge.1D0) then
         ATYPE = 'LPN'
C special case which doesn't occur in the non-normal form
      else if (NUSQ.eq.0.25D0) then
         ATYPE = 'R'

      else if (NUSQ.ge.0D0) then
         ATYPE = 'LCN'

      else
         ATYPE = 'LCO'
      end if
C In all cases take:
      A2 = 0D0
      A1 = 1D0

      BTYPE = 'R'
      SYM = .FALSE.
      B2 = 0D0
      B1 = 1D0
      go to 1001

C Remaining bits assume NUSQ & NU have been set:
  432 P = 1D0
      Q =(NUSQ-0.25D0)/(X*X)
      W = 1D0
      go to 1002
C BCs  handle all cases:
  433 if (NUSQ.ge.1D0) then
C        LPN case
         U = A2
         PDU = A1

      else if (NUSQ.eq.0.25D0) then
C         Regular case
         U = A1*X + A2
         PDU = A1

      else if (NUSQ.ge.0D0) then
C        LCN case, Fried. BC function=x**(1/2+nu), other=x**(1/2-nu)
         UF = X**(0.5D0+NU)
         UNF = X**(0.5D0-NU)

         PDUF =(0.5D0+NU)*UF/X
         PDUNF =(0.5D0-NU)*UNF/X

         U = A1*UF + A2*UNF
         PDU = A1*PDUF + A2*PDUNF

      else
C        LCO case, BC function RealPart(x**(1/2+i*nu)(a1-i*a2))
         TEMP = NU*LOG(X)
         U =(A1*COS(TEMP)+A2*SIN(TEMP))*SQRT(X)
         PDU =((0.5D0*A1+NU*A2)*COS(TEMP)-(NU*A1-0.5D0*A2)*SIN(TEMP))/
     +         SQRT(X)
      end if

      go to 1003

  434 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad LcN $(0\leq\alpha<0.25)$, LcO $(\alpha<0)$ \\
C $b = 1$ \qquad Regular \qquad  $u(b) = 0$
C Number of eigenvalues: $\infty\quad(\alpha\geq0)$,
C   $\pm\infty\quad(\alpha<0)$, \qquad  \Spc: none  \qquad
C $\alpha=0.01$: $\lambda_0=6.540555712 \qquad\lam_{24}=6070.441468$\\
C $\alpha=0$: $\lambda_0=5.78318596295 \qquad\lam_{24}=6045.999522$\\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #44
  440 TITLE = 'Border of LPN and LCN. -u" + x^(alpha-2) u = lambda u'
      NPARM = 1
      PARNM = 'alpha(near 0)'
      go to 1000

  441 A = 0.D0
      B = 1.D0
      ALPHA = PARM(1)
      if (ALPHA.gt.0D0) then
         ATYPE = 'LCN'
      else
         ATYPE = 'LPN'
      end if

      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  442 P = 1D0
      Q = X**(ALPHA-2D0)
      W = 1D0
      go to 1002

C The BC function in the book is in terms of Bessel fns.
C This one is wrong: only the 1st 2 terms of a series that needs
C more & more as alpha decreases to 0
  443 if (ALPHA.gt.0D0) then
         UF = X*(1D0 + X**ALPHA/(ALPHA*(1D0+ALPHA)))
         PDUF = 1D0+X**ALPHA/ALPHA
         UNF =(1D0 + X**ALPHA/(ALPHA*(ALPHA-1D0)))
         PDUNF = X**(ALPHA-1D0)/ALPHA
         U = A1*UF + A2*UNF
         PDU = A1*PDUF + A2*PDUNF

      else
         U = 0D0
         PDU = 1D0
      end if

      go to 1003

  444 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad Changes from LcN to LpN as $\alpha$ decreases thru $-2$
C  \qquad Fried. bc: $u(a) = 0$ \\
C $b = 1$  \qquad Regular  \qquad  $u(b) = 0$
C Number of eigenvalues: $\infty $ \qquad  \Spc: none  \qquad
C $\alpha=-1.99$: $\lam_0=15.87305674 \qquad \lam_{24} = 6316.899940$\\
C $\alpha=-2.01$: $\lam_0=15.96808975 \qquad \lam_{24} = 6325.038047$

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #45
  450 TITLE = 'Border of LCN and LCO. -u" - x^(alpha-2) u = lambda u'
      NPARM = 1
      PARNM = 'alpha(near 0)'
      go to 1000

  451 A = 0D0
      B = 1D0
      ALPHA = PARM(1)
      if (ALPHA.gt.0D0) then
         ATYPE = 'LCN'
      else
         ATYPE = 'LCO'
      end if

      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  452 P = 1D0
      Q = -X**(ALPHA-2D0)
      W = 1D0
      go to 1002

C I'm not too sure of this one(JDP):
C It should be analysed in same way as #44
  453 if (ALPHA.gt.0D0) then
         U = A1*(X-X**(1D0+ALPHA)/(ALPHA*(1D0+ALPHA))) + A2*0D0
         PDU = A1*(1D0-X**ALPHA/ALPHA) + A2*0D0

      else if (ALPHA.eq.0D0) then
         U = A1*(SIN(SQRT(0.75)*LOG(X))*SQRT(X)) + A2*0D0
         PDU = A1*(SIN(SQRT(0.75)*LOG(X)+PI/6D0)/SQRT(X)) + A2*0D0
      else
         write(*,FMT=*) 'Problem',IPROB,': LCO BCs not implemented yet'
         stop

      end if

      go to 1003

  454 U = B2
      PDU = B1
      go to 1004
C $a=0$ \qquad Changes from LcN to LcO as $\alpha$ decreases thru $0$
C \\
C Number of eigenvalues: $\pm \infty $  \qquad  \Spc: none\\
C $\alpha=-1.9$: $\lam_0 = -6152??.  \qquad\lam_2 = 11.38??$\\
C $\alpha=-2.01$: requires LcO BC at $x=0$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Singular problems: LCN, Non-Friedrichs BC
C Problem #46
  460 TITLE =
     +'Bessel normal form, nonFried. BC, ref. Bailey Everitt Zettl 1991'
      NPARM = 0
      PARNM = ' '
      go to 1000

C This is just Problem 13 with \nu=3/4
  461 NU = 0.75D0
      A = 0D0
      B = 1D0
      ATYPE = 'LCN'
C But with non-default left-hand BC defined by g=x^(-\nu+1/2):
      A2 = 1D0
      A1 = 0D0

      BTYPE = 'R'
      SYM = .FALSE.
      B2 = 0D0
      B1 = 1D0
      go to 1001

  462 P = 1D0
      Q =(5D0/16D0)/(X*X)
      W = 1D0
      go to 1002

  463 continue
C     LCN case, Fried. BC function=x**(1/2+nu), other=x**(1/2-nu)
      UF = X**(1.25D0)
      UNF = X**(-0.25D0)

      PDUF =(1.25D0)*UF/X
      PDUNF =(-0.25D0)*UNF/X

      U = A1*UF + A2*UNF
      PDU = A1*PDUF + A2*PDUNF
      go to 1003

  464 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #47
  470 TITLE =
     +'NonFriedrichs, Pryce/Marletta. -(xu'')''+(b/x)^2 u = lam x^(2b)u'
      NPARM = 1
      PARNM = 'b(>0)'
      go to 1000

  471 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'R'
      SYM = .FALSE.
      BETA = PARM(1)
      A2 = 1.D0
      A1 = 0.D0
      B2 = 1.D0
      B1 = -BETA
      go to 1001

  472 P = X
      Q = BETA**2/X
      W = X**(2D0*BETA)
      go to 1002

  473 U = A1*X**(BETA) + A2/X**(BETA)
      PDU = BETA*(A1*X**(BETA)-A2/X**(BETA))
      go to 1003

  474 U = B2
      PDU = B1
      go to 1004
C $a=0$ \qquad LcN \qquad non-Fried. bc: $[x^{-\beta},u](a)=0$\\
C $b=1$ \qquad regular \qquad $u'(b)=-\beta u(b)$\\
C $\lam_0=0$ for all $\beta$. For $\beta>14$(approx.) this appears
C to cause severe problems to all codes, see \scref{sngexpts}.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #48
  480 TITLE = 'Spectral density fn. example. Ref: Pruess/Fulton 75'
      NPARM = 0
      PARNM = ' '
      go to 1000

  481 A = 0D0
      B = 1D0
      ATYPE = 'LPNO'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  482 P = X**4
      Q = 0D0
      W = X**2
      go to 1002

  483 U = A2
      PDU = A1
      go to 1003

  484 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: 0 \qquad \Spc:(2.25,$\infty $) \\
C Spectral density $\rho(t) = {2\over 3\pi}(t-2.25)^{1.5}$
C The Frobenius index at 0 depends on lambda, and the p=x^4 can
C cause meshing problems.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Singular problems: LCO
C Problem #49
  490 TITLE =
     + 'Modified Bessel(Bessel of order nu=i/2). Ref: Pruess/Fulton 33'
      NPARM = 0
      PARNM = ' '
      go to 1000

  491 A = 0D0
      B = 1D0
      ATYPE = 'LCO'
      BTYPE = 'R'
      SYM = .FALSE.
      A1 = 1D0
      A2 = 0D0
      B1 = 1D0
      B2 = 0D0
      go to 1001

  492 P = X
      Q = -1/(4D0*X)
      W = X
      go to 1002

  493 TEMP = 0.5D0*LOG(X)
      U = A1*SIN(TEMP) + A2*COS(TEMP)
      PDU = 0.5D0*(A1*COS(TEMP)-A2*SIN(TEMP))
      go to 1003

  494 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  LcO \\
C $b = 1$  \qquad  Regular  \qquad  $u(b) = 0$
C Number of eigenvalues: $\pm \infty $  \qquad \Spc: none \\
C General BC at 0 is $[x^{1/2}\sin(\half\ln x+\gamma),u](0)=0$ where
C$\gamma$ is constant.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #50
  500 TITLE =
     +'Weierstrass-Mandelbrot. -u" + (c-(omega^2+.25)/x^2) u = lambda u'
      NPARM = 2
      PARNM =
     +'omega >0 & offset c; c>0 vital for SLEIGN2 to find lambda_0'
      go to 1000

  501 OMEGA = PARM(1)
      A = 0D0
      B = XINFTY
      ATYPE = 'LCO'
      BTYPE = 'LPNO'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .FALSE.
      go to 1001

  502 P = 1D0
      Q = PARM(2)-(OMEGA**2+0.25D0)/X**2
      W = 1D0
      go to 1002

  503 TEMP = OMEGA*LOG(X)
      U =(A2*COS(TEMP)+A1*SIN(TEMP))*SQRT(X)
      PDU =((0.5D0*A2+OMEGA*A1)*COS(TEMP) -
     +      (OMEGA*A2-0.5D0*A1)*SIN(TEMP))/SQRT(X)
      go to 1003

  504 U = B2
      PDU = B1
      go to 1004
C Number of eigenvalues: $+-\infty$ \qquad \Spc:(0,$\infty $)  \qquad
C BC functions f=x^(1/2) sin(omega ln(x) + gamma), gamma arbitrary
C For c=0 & any gamma the evs are a 2-sided infinite geometric
C progression with common ratio exp(2 pi/omega).
C For ratio=10 take omega = 2.72875270768368
C For ratio=2  take omega = 9.06472028365439

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #51
  510 TITLE = 'Titchmarsh problem. -u" - e^x u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  511 A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LCO'
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .FALSE.
      go to 1001

  512 P = 1D0
      Q = -EXP(X)
      W = 1D0
      go to 1002

  513 U = A2
      PDU = A1
      go to 1003

  514 TEMP = 2D0*EXP(0.5D0*X)
      U = EXP(-0.25D0*X)*(B1*SIN(TEMP)+B2*COS(TEMP))
      PDU = -0.5D0*U + 2D0*EXP(0.25D0*X)*(B1*COS(TEMP)-B2*SIN(TEMP))
      go to 1004
C $a = 0$  \qquad  Regular \qquad $u(a) = 0$  \\
C $b = +\infty $  \qquad LcO
C Number of eigenvalues: $\pm \infty $  \qquad \Spc: none \\
C General BC at $\infty$ is of form $[f,u](\infty)=0$ where
C $f=e^{-x/4}\sin(2e^{x/2} + \gamma)$, which is an exact solution of
C the DE with $\lam=-1/16$ for any value of the constant $\gamma$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #52
  520 TITLE =
     +    'Limit-circle osc with cont spectrum. -u" -1/x^3 u = lambda u'
      NPARM = 0
      PARNM = ' '
      go to 1000

  521 A = 0D0
      B = XINFTY
      ATYPE = 'LCO'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  522 P = 1D0
      Q = -1D0/X**3
      W = 1D0
      go to 1002

  523 continue
C     This needed JWKB expansion up to R_3 term (JDP book p41)
C     to get approximations u,v s.t. Lu,Lv are square integrable.
      R = X**0.75D0*exp((3D0/64D0)*X)
      THETA = 2D0/sqrt(X) + (3D0/16D0)*sqrt(X)
      C = 0.75D0/X + 3D0/64D0
      D = (-1D0/X + (3D0/32D0))/sqrt(X)
      UF = R*cos(THETA)
      UNF = R*sin(THETA)
      PDUF = C*UF-D*UNF
      PDUNF = D*UF+C*UNF
      U = A1*UF + A2*UNF
      PDU = A1*PDUF + A2*PDUNF
      go to 1003

  524 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  LcO, BC functions x^{0.5) Z(2x^{-1/2})
C where Z is a combination of Bessel fns J1 and Y1.\\
C Elementary BC functions from real, imaginary parts of $w$ where
C \[w'/w = -ix^{-3/2} + {3\over4}x^{-1} + i{3\over32}x^{-1/2}
C    + {3\over64}x^0 \]
C obtained by terms up to $R_3$ in JWKB expansion (book p41).
C $b = +\infty $  \qquad LpN/O \\
C \Spc:(0,$\infty $)   \\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problems with discontinuous coefficients
C Problem #53
  530 TITLE = 'Approximate harmonic oscillator, Pryce 1993'
      NPARM = 1
      PARNM =
     +   'b (reg BCs u(-b)=u(b)=0 are used, b>=10^35 taken as Infinity)'
      go to 1000

  531 B = MIN(PARM(1),XINFTY)
      A = -B
      if (B.lt.XINFTY) then
         ATYPE = 'R'
         BTYPE = 'R'

      else
         ATYPE = 'LPN'
         BTYPE = 'LPN'
      end if

      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      SYM = .TRUE.
      go to 1001

  532 P = 1D0
      T = ABS(X)
      N = FLOAT(INT(T))
      Q = N*N + (2*N+1)*(T-N)
      W = 1D0
      go to 1002

  533 U = A2
      PDU = A1
      go to 1003

  534 U = B2
      PDU = B1
      go to 1004
C q is piecewise linear function joining the points(x,x^2) for x an
C integer. Codes seem to find automatic endpoint analysis hard,
C hence the PARM which allows one to set finite regular endpoints
C at which u=0 is imposed.
C Number of eigenvalues: $\infty $ \qquad     \Spc: none \\
C $\lambda_0=.3456792711$, $\lambda_9=18.079846195$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Extensions to the classical SLP
C Problem #54
  540 TITLE = 'Indefinite weight function, Marletta PhD 1991'
      NPARM = 0
      PARNM = ' '
      go to 1000

  541 A = 0D0
      B = 1D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  542 P = 1D0
      Q = 1D0
      if (X.le.0.5D0) then
         W = 0D0
      else
         W = 1D0
      end if

      go to 1002

  543 U = A2
      PDU = A1
      go to 1003

  544 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  Regular \qquad $u(a)=0$\\
C $b = 1$  \qquad  Regular \qquad $u(b)=0$\\
C Number of eigenvalues: $\infty $.
C $\lam_k=\mu_k^2+1$ where $\mu_k$ are roots of
C $\tan{\mu\over2}+\mu\tanh\half=0$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #55
  550 TITLE = 'Indefinite weight fn, Bailey et al. ACM TOMS 1978'
      NPARM = 0
      PARNM = ' '
      go to 1000

  551 A = -1D0
      B = 1D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  552 P = EXP(-X**2)
      Q = P
      W = -2D0*X*P
      go to 1002

  553 U = A2
      PDU = A1
      go to 1003

  554 U = B2
      PDU = B1
      go to 1004
C $a = -1$  \qquad  Regular \qquad $u(a)=0$\\
C $b = 1$  \qquad  Regular \qquad $u(b)=0$\\
C Number of eigenvalues: $\pm\infty $.
C The conversion to self-adjoint form of the problem
C u'' - x u' -(1+2\lambda x) u = 0, u(-1)=u(1)=0
C $\lambda_5=250.974149734$.
C By the antisymmetry, if lambda is an eigenvalue so is -lambda.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #56
  560 TITLE = 'lambda-dependent BC at x=a, Fulton/Pruess 1992'
      NPARM = 4
      PARNM = 'A1 A2 A1'' A2'' in SLEDGE-style BCs'//
     *', 1 0 0 1 for prob. in book'
      go to 1000

  561 A = 1D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPNO'
      SYM = .FALSE.
C The BC at x=a is not adjustable so A1, A2 not set:
      B2 = 0D0
      B1 = 1D0
      go to 1001

  562 P = 1D0
      Q = -0.25D0/X**2
      W = 1D0
      go to 1002

  563 U = PARM(1)-EIG*PARM(3)
      PDU = PARM(2)-EIG*PARM(4)
      go to 1003

  564 U = B2
      PDU = B1
      go to 1004
C $a = 1$  \qquad  Regular \qquad $\lambda u(a) -(pu')(a) = 0$\\
C $b = \infty$  \qquad  LPN/O\\
C Number of eigenvalues: $1$ \qquad     \Spc: $(0,\infty)$ \\
C $\lambda_0=-0.02229899486$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #57
  570 TITLE =
     +    'lambda-dependent BC at x=a, Fulton/Pruess private comm. 1992'
      NPARM = 4
      PARNM = 'A1 A2 A1'' A2'' in SLEDGE-style BCs'//
     *', 1 0 0 1 for prob. in book'
      go to 1000

  571 A = 1D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPNO'
      SYM = .FALSE.
C The BC at x=a is not adjustable so A1, A2 not set:
      B2 = 0D0
      B1 = 1D0
      go to 1001

  572 P = 1D0
      Q = -1D0/X
      W = 1D0
      go to 1002

  573 U = PARM(1)-EIG*PARM(3)
      PDU = PARM(2)-EIG*PARM(4)
      go to 1003

  574 U = B2
      PDU = B1
      go to 1004
C $a = 1$  \qquad  Regular \qquad $\lambda u(a) -(pu')(a) = 0$\\
C $b = \infty$  \qquad  LPN/O\\
C Number of eigenvalues: $\infty$ \qquad     \Spc: $(0,\infty)$ \\
C $\lambda_0=-0.30017359365$ \qquad $\lambda_9=-0.0025669706$.

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #58
  580 TITLE = 'Rossby wave equation, Drazin et al. J Fluid Mech 1982'
      NPARM = 2
      NEPRM = 2
      PARNM = 'alpha, c'
      go to 1000

  581 A = -PI
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

C This problem is quite far from the classical SLP and one cannot
C classify it for SLEIGN (q finite at a or b, etc) in a simple way.
C In the following section:
C BETA is beta, the standard 'linear' eigenparameter.
C VEL is called U(x) in the cited paper.
C To provide an interface for problems nonlinear in parameters,
C for a NAG D02KEF-style solver, if IPARM>0,
C writing QQ=beta*w-q, we compute WW=partial d(QQ)/d(alpha),
C WW = partial d(QQ)/d(c) according as IPARM is 1 or 2
C and return QQ, WW in variables Q and W.
  582 BETA = PARM(0)
      ALPHA = PARM(1)
      C = PARM(2)
      VEL =(X/B)**2
      P = 1D0
      Q = ALPHA**2 + VEL/(VEL-C)
      W = 1D0/(VEL-C)
      if (IPARM .ne. 0) then
         Q = BETA*W - Q
         if (IPARM .eq. 1) then
            W = -2D0*ALPHA
         else if (IPARM .eq. 2) then
            W =(BETA-VEL)/(VEL-C)**2
         end if
      end if
      go to 1002

  583 U = A2
      PDU = A1
      go to 1003

  584 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problems contrived to show special features
C Problem #59
  590 TITLE = 'Ref Gelfand & Levitan, AMS Translations 1955.'
      NPARM = 0
      PARNM = ' '
      go to 1000

  591 A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 1D0
      A1 = -1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  592 P = 1D0
      T = 1D0 + X/2D0 + SIN(2D0*X)/4D0
      Q = 2D0*(T*SIN(2D0*X)+COS(X)**4)/T**2
      W = 1D0
      go to 1002

  593 U = A2
      PDU = A1
      go to 1003

  594 U = B2
      PDU = B1
      go to 1004
C $a = 0$ \qquad Regular \qquad $u(a)+ (pu')(a) = 0$ \\
C $b = \infty$ \qquad LPN/O \\
C Number of(isolated) eigenvalues: none \qquad \Spc: $(0,\infty)$ \\
C There is an eigenvalue at $\lambda=1$ in the \Spc. \\

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #60
  600 TITLE =
     + 'Problem with pseudo-eigenvalue, Pryce 1989, ref Marletta PhD 91'
      NPARM = 0
      PARNM = ' '
      go to 1000

  601 A = 0D0
      B = XINFTY
      ATYPE = 'R'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = -8D0
      A1 = 5D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  602 P = 1D0
      Q =(3D0*(X-31D0))/(4D0*(1D0+X)*(4D0+X)**2)
      W = 1D0
      go to 1002

  603 U = A2
      PDU = A1
      go to 1003

  604 U = B2
      PDU = B1
      go to 1004
C $a = 0$  \qquad  Regular \qquad $5u(a)+8(pu')(a) = 0$ \\
C $b = +\infty $  \qquad LpN/O \\
C Number of eigenvalues: 1  \qquad \Spc:(0,$\infty $) \qquad
C $\lam_0 = -1.185214105$\\
C There is a solution for $\lam=0$ satisfying the BC at 0 and converging
C to 0 at $\infty$. However it is not square-integrable, so $0$ is not
C an eigenvalue.

C***+****|****+****|      END OF PROBLEM SET     |****+****|****+****|**

 1000 continue
 1001 continue
 1002 continue
 1003 continue
 1004 continue
      end subroutine TSTSET

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine ABINFO(P0ATA,QFATA,P0ATB,QFATB)
C This routine is used only by SLEIGN & SLEIGN2.
C Its functionality should be incorporated into TSTSET eventually.
C******* THIS ROUTINE APPLIES TO THE STANDARD TEST SET ONLY ************
      implicit none
      integer NPROB
      parameter(NPROB=60)
C     .. Scalar Arguments ..
      double precision P0ATA,QFATA,P0ATB,QFATB
C     .. Local Arrays ..
      integer P0A(1:NPROB),QFA(1:NPROB),P0B(1:NPROB),QFB(1:NPROB)

      data P0A/
C      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     + 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1,
C     21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
     + 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
C     41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
     + 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/

      data QFA/
C      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     + 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
C     21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
     + 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
C     41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
     + 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1/

      data P0B/
C      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     + 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 
C     21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
     + 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
C     41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
     + 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/

      data QFB/
C      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     + 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 
C     21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
     + 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 
C     41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
     + 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1/

C***+****|****+****|**** Special Cases Section **|****+****|****+****|**
C in several problems this info depends on values of parameters so:
      goto(999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     +    999, 120, 130, 999, 150, 160, 999, 999, 999, 999,
     +     999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     +     999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     +     410, 999, 430, 440, 450, 999, 999, 999, 999, 999,
     +     999, 999, 999, 999, 999, 999, 999, 999, 999, 999)
     +     IPROB

C Problem #12, Bessel: q finite at 0 if nusq=0
  120 if (PARM(1).eq.0D0) QFATA=1
      go to 999

C Problem #13, Bessel LNF: q finite at 0 if nusq=1/4
  130 if (PARM(1).eq.0.25D0) QFATA=1
      go to 999

C Problem #15, Assoc Legendre: q finite at -1, +1 if c=0
  150 if (PARM(1).eq.0D0) then
        QFATA=1
        QFATB=1
      end if
      go to 999

C Problem #16, Assoc Legendre LNF: q finite at -1, +1 if c=1/4
  160 if (PARM(1).eq.0.25D0) then
        QFATA=1
        QFATB=1
      end if
      go to 999

C Problem #41, Woods-Saxon: q finite if l = -1 or 0
  410 if (PARM(1).eq.-1D0 .or. PARM(1).eq.0D0) QFATA=1
      go to 999

C Problem #43, `borderline 1', is a copy of Prob 13
  430 if (PARM(1).eq.0.25D0) QFATA=1
      go to 999

C Problem #44, `borderline 2': q finite at 0 iff alpha>=2
  440 QFATA = 0
      if (PARM(1).ge.2D0) QFATA=1
      go to 999

C Problem #45, `borderline 3': q behaves like #44
  450 QFATA = 0
      if (PARM(1).ge.2D0) QFATA=1
      go to 999

C***+****|****+****|**** End of Special Cases  **|****+****|****+****|**
 999  continue
      P0ATA = 2*P0A(IPROB) - 1
      QFATA = 2*QFA(IPROB) - 1
      P0ATB = 2*P0B(IPROB) - 1
      QFATB = 2*QFB(IPROB) - 1

      end subroutine ABINFO

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine GTHUHV(IEND,XA,UA,VA,HU,HV)
C This routine is used only by SLEIGN2.
C Its functionality should be incorporated into TSTSET eventually.
C******* THIS ROUTINE APPLIES TO THE STANDARD TEST SET ONLY ************
      implicit none
      integer IEND
      double precision XA,UA,VA,HU,HV
      double precision AA,BB,ALPHA,NU,OM,OMSQ,TEMP

C Returns the HU,HV info needed by SLEIGN2 as part of its UV routine.
C (May be included in main TSTSET eventually)
C     Kluge, to set valid fl.pt. values: SLEIGN2 does funny things..
      HU = 0D0
      HV = 0D0

      if (ATYPE(1:2).eq.'LC' .or. BTYPE(1:2).eq.'LC') then
        go to(10,10,10,10,10,10,10,10,10, 10, 10,120,130, 10,150,160,
     +          10,180,190,200,210,220,230,240,250, 10, 10, 10, 10,300,
     +         310, 10, 10,340, 10, 10, 10, 10, 10, 10, 10,420, 10,440,
     +         450,460,470, 10,490,500,510,520, 10, 10, 10, 10, 10, 10,
     +          10, 10) IPROB
      else
C       If neither endpoint is LC, UV is not needed so do nothing
        return
      end if

   10 write(*,*) 'Error in UV, Problem',IPROB,' has LC endpoint but',
     +  ' no HU,HV info implemented for it'
      stop

  120 if (IEND.eq.0) then
C     Whatever the value of PARM(1)=NU**2:
      HU = 0D0
      HV = 0D0
      end if
      go to 1000

  130 if (IEND.eq.0) then
C     Whatever the value of PARM(1)=NU**2:
      HU = 0D0
      HV = 0D0
      end if
      go to 1000

  150 continue
C     Whether IEND is 0 or 1:
      if (PARM(1).ge.0D0) then
        NU = sqrt(PARM(1))
        HU =(NU*NU+NU)*UA
        HV =(NU*NU-NU)*VA
      else
C       Here UA,VA are Re,Im of(1-x*x)**{i*omega}, nu=i*omega
        OMSQ = PARM(1)
        OM = sqrt(-OMSQ)
        HU = OMSQ*UA-OM*VA
        HV = OMSQ*VA+OM*UA
      end if
      go to 1000

  160 continue
C     Whatever the value of PARM(1)=NU**2 and at both ends:
      HU = 0D0
      HV = 0D0
      go to 1000

  180 if (IEND.eq.0) then
      HU = 0D0
      HV = 0D0
      end if
      go to 1000

  190 if (IEND.eq.0) then
      HU = 0D0
      HV = 0D0
      end if
      go to 1000

  200 continue
C     At both ends:
      HU = 0D0
      HV = 0D0
      go to 1000

  210 continue
C for truncated series BC functions(see TSTSET) not the Bessel ones
      HU = 400D0/3D0
      HV = HU*log(XA) - 3200D0/9D0
      go to 1000

  220 continue
C     At both ends:
      HU = 0D0
      HV = 0D0
      go to 1000

  230 continue
      print*,'UV: Prob 23 not implemented yet, program terminated'
      stop
      go to 1000

  240 if (IEND.eq.1) then
        HU = 0D0
        HV = - (1D0+XA*(2D0+XA*(3D0+XA*(4D0+XA*(5D0+XA*6D0)))) )
      end if
      go to 1000

  250 if (IEND.eq.0) then
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  300 if (IEND.eq.0) then
        HU = -1D0
        HV = log(XA)
      end if
      go to 1000

  310 if (IEND.eq.0) then
        HU = -XA
        HV = log(XA)
        print*,'UV: Prob 31 Warning, only valid for l=0'
      end if
      go to 1000

  340 if (IEND.eq.0) then
C       Probably no need to worry about overflow in cosh x here:
        HU = -sqrt(XA)/(cosh(XA)**2)
        HV = HU*log(XA)
      end if
      go to 1000

  420 if (IEND.eq.0) then
        HU = 5D0*exp(-2D0*XA) - 1D0
        if (abs(XA).le.1D-4) then
          TEMP = -2D0+2D0*XA-(4D0/3D0)*XA*XA+ (2D0/3D0)*XA*XA*XA
        else
          TEMP = 5D0*(exp(-2D0*XA)-1D0)/XA
        end if
        HV = TEMP + (5D0*exp(-2D0*XA)-1D0)*4D0*log(XA)
        print*,'UV: Prob 42 Warning, only valid for l=0'
      end if
      go to 1000

  440 if (IEND.eq.0) then
        print*,'UV: Prob 44, BC fns are known NOT to be correct!'
        ALPHA = PARM(1)
        HU = XA**(2D0*ALPHA-1D0)/((ALPHA+1D0)*ALPHA)
        HV = XA**(2D0*ALPHA-2D0)/(ALPHA*(ALPHA-1D0))
      end if
      go to 1000

  450 if (IEND.eq.0) then
        print*,'UV: Prob 45, BC fns are known NOT to be correct!'
        ALPHA = PARM(1)
C this is false!
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  460 if (IEND.eq.0) then
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  470 if (IEND.eq.0) then
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  490 if (IEND.eq.0) then
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  500 if (IEND.eq.0) then
        HU = 0D0
        HV = 0D0
      end if
      go to 1000

  510 if (IEND.eq.1) then
        HU = UA/(-16D0)
        HV = VA/(-16D0)
      end if
      go to 1000

  520 if (IEND.eq.0) then
        AA = (63D0/1024D0)/XA + 9D0/4096D0
        BB = (9D0/1024D0)/sqrt(XA)
        HU = AA*UA - BB*VA
        HV = AA*VA + BB*UA
      end if
      go to 1000

 1000 end subroutine GTHUHV

      end module TESTMOD
