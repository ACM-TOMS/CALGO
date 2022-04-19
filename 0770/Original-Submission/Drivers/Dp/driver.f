
        PROGRAM DBVSDR

C  ABSTRACT
C
C  DBVSDR is a driver program which shows the use of the subroutines
C  DBVSIS, DBVSSC and DBVSSE with some examples.
C
C
C  DESCRIPTION
C
C  DBVSDR consists of five examples which handle the execution of
C  DBVSSC (Boundary-Valued Shape-preserving Spline Computation) and
C  DBVSSE (Boundary-Valued Shape-preserving Spline Evaluation) or, when
C  convenient, the support routine DBVSIS, which simply collects the
C  calls to the previous ones, and cover the main features of the
C  package. The outputs are written to the external file RES.
C  For an easier reading of the program, all the output statements
C  are collected in the subroutine DWRTO.
C  The first and third examples make use of the data file data1 the
C  second of data2; obviously, they must be stored in the same 
C  working directory as the executable version of DBVSDR.
C  The second and fourth examples make also use of external functions
C  which characterize the boundary conditions; although when the
C  non-separable boundary conditions are not imposed this function are
C  meaningless, vacuous function subprograms must be supplied for
C  linking purposes (for details see the description of the method and
C  of the input parameters BETA, BETAI, RHO and RHOI in DBVSIS).
C
C
C  DBVSDR makes use of:
C
C  a) variables and arrays used as actual arguments for DBVSIS.
C     They have the same names of the corresponding dummy arguments
C     and their meaning is explained in the above subroutine;
C
C  b) symbolic constant names, used for dimensioning the arrays. They
C     have the mnemonic form P-name, where P stands for PARAMETER and
C     'name' is the name of a constant or variable which is used in 
C     DBVSIS. So, their meaning can be easily obtained from the related
C     constant or variable 'name' described in the above subroutine;
C
C  c) symbolic constants and variable names, used for input operations.
C     They will be explained at their first occurrence in executable 
C     statements;
C
C  d) function names used as actual arguments for DBVSIS. They have one
C     of the following forms DBn , DBnI , DRn , DRnI or DDUMFN. The
C     number n refers to the example in which the actual function is
C     used, B and R refer to functions 'beta' or 'rho' described in
C     DBVSIS and I , when present, to their inverse. DDUMFN is finally
C     used for linking purposes when 'beta' and/or 'rho' functions are
C     not referenced.
C
C  The following are in class a)
C
C  X,Y,NP,N,K,OPT,D0.DNP,D20,D2NP,CONSTR,EPS,KMAX,MAXSTP,XTAB,NTAB,
C  SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,D,D2,Y0TAB,Y1TAB,Y2TAB,WORK,NWORK;
C
C  in class b)
C
C  PNP,PNTAB,PN,PCOMM,PPART,PNWORK;
C
C  in class c)
C
C  UN1,UN2,FNM1,FNM2;
C
C  and in class d)
C
C  DB2,DB2I,DR2,DR2I,DB4,DB4I,DDUMFN.
C
C  DBVSDR calls the subroutines DBVSSC, DBVSSE and DWRTO
C  (at the entry points DWRT1A, DWRT1B, DWRT2, DWRT3A, DWRT3B, DWRT4, 
C   DWRT5A, DWRT5B ).

        INTEGER PNP,PNTAB,PPART,PCOMM,PN,PNWORK
        PARAMETER (PNP=10,PNTAB=15,PPART=2,PCOMM=5,PN=6,
     *             PNWORK=PCOMM+(PPART+7)*PNP+PN*(PN+11)/2+9)

C  The following statements define the symbolic names of constants used
C  for the input operations.

        INTEGER UN1,UN2
        CHARACTER FNM1*9,FNM2*9
        PARAMETER (UN1=1,UN2=2,FNM1='data1',FNM2='data2')

        INTEGER NP,N,K,OPT,CONSTR(0:PNP-1),KMAX,MAXSTP,NTAB,SBOPT,
     *          Y0OPT,Y1OPT,Y2OPT,ERRE,ERRC,DIAGN(0:PNP-1),NWORK,I

        DOUBLE PRECISION X(0:PNP),Y(0:PNP),D0,DNP,D20,D2NP,EPS,
     *                   XTAB(0:PNTAB),D(0:PNP),D2(0:PNP),
     *                   Y0TAB(0:PNTAB),Y1TAB(0:PNTAB),
     *                   Y2TAB(0:PNTAB),WORK(1:PNWORK)

        DOUBLE PRECISION DB2,DB2I,DR2,DR2I,DB4,DB4I,DDUMFN

        EXTERNAL DB2,DB2I,DR2,DR2I,DB4,DB4I,DDUMFN

C  The following variables will be used in the output statements of
C  subroutine DWRTO.

        COMMON /DBVCOM/ NP,NTAB,ERRC,ERRE,X,Y,XTAB,Y0TAB,Y1TAB,Y2TAB,
     *                  D,D2,DIAGN,CONSTR

C  ---------------------------------------------------------------------

C  Call subroutine DWRTO for printing the heading in the output file.

        CALL DWRTO

C  ---------------------------------------------------------------------

C  Example n.1 .
C
C  We want to construct a spline interpolating the set of data points
C  (x(i),y(i)) , i=0,1,...,10 , with  x(i) = 2*i/20  and
C  y(i)= sin x(i) + ( cos 5x(i) )/5 , stored in the external file FNM1.
C  No constraints or boundary conditions are imposed on it and its
C  slopes are computed as the best approximation to the Bessel
C  estimates. The spline has continuity 1 and degree 3. Note that with
C  these selections we get the classical cubic Bessel interpolant.
C  We also want to estimate the spline at the set of tabulation points
C  xtab(i)=x(i) , i=0,...,10 . Then we want to evaluate the first and
C  the second derivative at the point  1.39 .
C
C
C  Read the data from the input file. It contains:
C  - the number NP of the data points,
C  - the coordinates (x(i),y(i)) of the i-th point, i=0,...,NP.

        OPEN  (UNIT=UN1,FILE=FNM1,STATUS='OLD',ACCESS='SEQUENTIAL',
     *        FORM='FORMATTED')
        READ  (UNIT=UN1,FMT=*) NP
        READ  (UNIT=UN1,FMT=*) (X(I),Y(I) , I=0,NP)
        CLOSE (UNIT=UN1,STATUS='KEEP')

C  Set the degree and the class of continuity of the spline.

        N=3
        K=1

C  Set the size of the work area

        NWORK=PCOMM+(PPART+7)*NP+(N*(N+11))/2+9

C  Set the options on the spline: slopes as the best approximation to
C  the Bessel estimates, no boundary conditions, no constraints.

        OPT=210

C  Set EPS, KMAX and MAXSTP to 0 for the automatic choice of their
C  values.

        EPS=0.D+00
        KMAX=0
        MAXSTP=0

C  Call DBVSSC to compute the spline's parameters.

        CALL DBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,DDUMFN,
     *               DDUMFN,DDUMFN,DDUMFN,KMAX,MAXSTP,ERRC,D,D2,DIAGN,
     *               WORK,NWORK)


C  Set the options for evaluating the spline at the tabulation points.

        Y0OPT=1
        Y1OPT=0
        Y2OPT=0

C  Define the number and the set of tabulation points.

        NTAB=10
        DO 10 I=0,NTAB
            XTAB(I)=X(I)
10      CONTINUE

C  Set SBOPT to 1 for selecting the sequential search algorithm.

        SBOPT=1

C  Call DBVSSE for evaluating the spline.

        CALL DBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *               ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

C  Call DWRTO at the entry point DWRT1A for printing the results in the
C  output file.

        CALL DWRT1A

C  Set the options for evaluating the spline at the tabulation point.

        Y0OPT=0
        Y1OPT=1
        Y2OPT=1

C  Define the new tabulation point.

        NTAB=0
        XTAB(0)=0.139D+01

C  Set SBOPT to 2 for selecting the binary search algorithm.

        SBOPT=2

C  Call DBVSSE for evaluating the spline.

        CALL DBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *               ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

C  Call DWRTO at the entry point DWRT1B for printing the results in the
C  output file.

        CALL DWRT1B

C  ---------------------------------------------------------------------

C  Example n.2 .
C
C  Given the set of points (x(i),y(i)) , i=0,...,8 where
C       {x(0),...,x(8)}={0.,1.,2.,3.,4.,5.,6.,7.,8.}   and
C    {y(0),...,y(8)}={0.,0.75,1.,0.75,0.,-0.75,-1.,-0.75,0.}
C  stored in the external file FNM2, we want to construct a
C  co-monotone spline subject to the additional boundary condition 
C  s'(8)=s'(0) ,  s''(8)=s''(0) .
C  We also ask the slopes to be the best approximation to the set of
C  third order monotonicity-preserving estimates described in DBVSSC
C  (see the input parameter OPT and the subroutine TODC).
C  The spline has degree 6 and continuity 2.
C  We evaluate the spline and its second derivative at the set of
C  tabulation points xtab(i), i=0,...,10 uniformly distributed in the
C  interval  (x(0),x(8)).
C  We use in this case the support routine DBVSISP.
C
C  Read the data from the input file. It contains:
C  - the number NP of the data points,
C  - the coordinates (x(i),y(i)) of the i-th point, i=0,...,8 .

        OPEN  (UNIT=UN2,FILE=FNM2,STATUS='OLD',
     *        ACCESS='SEQUENTIAL',FORM='FORMATTED')
        READ  (UNIT=UN2,FMT=*) NP
        READ  (UNIT=UN2,FMT=*) (X(I),Y(I),I=0,NP)
        CLOSE (UNIT=UN2,STATUS='KEEP')

C  Set the degree and the class of continuity of the spline.

        N=6
        K=2

C  Set the size of the work area

        NWORK=PCOMM+(PPART+7)*NP+(N*(N+11))/2+9

C  Set the options on the spline: slopes as the best approximation to
C  the set of third order accurate estimates, non-separable boundary
C  conditions, monotonicity constraints.

        OPT=321

C  Set EPS, KMAX and MAXSTP to 0 for the automatic choice of their
C  values.

        EPS=0.0D+00
        KMAX=0
        MAXSTP=0

C  Set the options for evaluating the spline and its second derivative
C  at the tabulation points.

        Y0OPT=1
        Y1OPT=0
        Y2OPT=1

C  Define the number and the set of tabulation points.

        NTAB=10
        DO 20 I=0,NTAB
            XTAB(I)=X(0)+I*(X(NP)-X(0))/NTAB
20      CONTINUE

C  Set SBOPT to 1 for the sequential search algorithm.

        SBOPT=1

C  Call DBVSIS for computing the spline.

        CALL DBVSIS (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,
     *               DB2,DB2I,DR2,DR2I,KMAX,MAXSTP,XTAB,NTAB,
     *               SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,D,D2,DIAGN,
     *               Y0TAB,Y1TAB,Y2TAB,WORK,NWORK)

C  Call DWRTO at the entry point DWRT2 for printing the results into
C  the output file.

        CALL DWRT2

C  ---------------------------------------------------------------------

C  Example n.3 .
C
C  Given the set of data points (x(i),y(i)) i=0,1,...,10 as in the
C  example 1, we want to construct a spline of degree 6 and continuity 2
C  which is co-convex and subject to separable boundary conditions:
C  s'(x(0))=y'(x(0)) , s'(x(10))=y'(x(10)) , s''(x(0))=y''(x(0)) ,
C  s''(x(10))=y''(x(10)) . The parameters are determined only by the 
C  convexity constraints. We evaluate the spline and its first 
C  derivative at the tabulation points xtab(i) , i=0,...,15 uniformly 
C  distributed in the interval  (x(0),x(10))  and then the second
C  derivative at the knots.
C
C
C  Read the data from the input file.

        OPEN  (UNIT=UN1,FILE=FNM1,STATUS='OLD',
     *        ACCESS='SEQUENTIAL',FORM='FORMATTED')
        READ  (UNIT=UN1,FMT=*) NP
        READ  (UNIT=UN1,FMT=*) (X(I),Y(I),I=0,NP)
        CLOSE (UNIT=UN1,STATUS='KEEP')

C  Set the degree and the class of continuity of the spline.

        N=6
        K=2

C  Set the size of the work area

        NWORK=PCOMM+(PPART+7)*NP+(N*(N+11))/2+9

C  Set the options on the spline: the slopes are computed only from
C  the constraints, non-separable boundary conditions,
C  convexity constraints.

        OPT=132

C  Set the boundary conditions.

        D0=COS(X(0))+SIN(5.0*X(0))
        DNP=COS(X(NP))+SIN(5.0*X(NP))
        D20=-SIN(X(0))-5.0*COS(5.0*X(0))
        D2NP=-SIN(X(NP))-5.0*COS(5.0*X(NP))

C  Set EPS, KMAX and MAXSTP for the automatic choice of their values.

        EPS=0.0D+00
        KMAX=0
        MAXSTP=0

C  Call DBVSSC for computing the spline.

        CALL DBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,DDUMFN,
     *               DDUMFN,DDUMFN,DDUMFN,KMAX,MAXSTP,ERRC,D,D2,DIAGN,
     *               WORK,NWORK)

C  Set the options for evaluating the spline and its first derivative
C  at the tabulation points.

        Y0OPT=1
        Y1OPT=1
        Y2OPT=0

C  Define the number and the set of tabulation points.

        NTAB=15
        DO 30 I=0,NTAB
            XTAB(I)=X(0)+I*(X(NP)-X(0))/NTAB
30      CONTINUE

C  Set SBOPT to 1 for selecting the sequential search algorithm.

        SBOPT=1

C  Call DBVSSE for evaluating the spline.

        CALL DBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *               ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

C  Call DWRTO at the entry point DWRT3A for printing the results in the
C  output file.

        CALL DWRT3A

C  Set the options for evaluating the second derivative at the knots.

        Y0OPT=0
        Y1OPT=0
        Y2OPT=1

C  Define the number and the set of tabulation points.

        NTAB=NP
        DO 40 I=0,NTAB
            XTAB(I)=X(I)
40      CONTINUE

C  Set SBOPT to 1 for selecting the sequential search algorithm.

        SBOPT=1

C  Call DBVSSE for evaluating the spline.

        CALL DBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *               ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

C  Call DWRTO at the entry point DWRT3B for printing the results in the
C  output file.

        CALL DWRT3B

C  ---------------------------------------------------------------------

C  Example n.4 .
C
C  Given the 'Runge function': f(x) = 1/(1+x**2) , we want to construct
C  a co-monotone and co-convex spline which interpolates the function at
C  x(i) = -1+2*i/4  , i=0,1,...,4 . We also impose the additional
C  boundary conditions  s'(1)=-s'(-1) and we compute the slopes as the
C  best approximation to the estimates  s'(x(i)) = f'(x(i)) , i=0,...,4.
C  The spline has degree 3 and continuity 1. We evaluate the spline and
C  its first derivative at the tabulation points xtab(i), i=0,...,15
C  uniformly distributed in the interval(x(0),x(4)).
C  It is used the support routine DBVSIS
C
C  Compute the input data

        NP=4
        DO 50 I=0,NP
            X(I)=-1.0+2.0*I/4.0
            Y(I)=1.0/(1.0+X(I)**2)
            D(I)=-2.0*X(I)/((1.0+X(I)**2)**2)
50      CONTINUE

C  Set the degree and the class of continuity of the spline.

        N=3
        K=1

C  Set the size of the work area

        NWORK=PCOMM+(PPART+7)*NP+(N*(N+11))/2+9

C  Set the options on the spline: best approximation to user supplied
C  derivatives, boundary conditions, monotonicity and convexity
C  constraints.

        OPT=423

C  Set EPS, KMAX and MAXSTP

        EPS=1.0D-05
        KMAX=5
        MAXSTP=15

C  Set the options for evaluating the spline and its first derivative
C  at the tabulation points.

        Y0OPT=1
        Y1OPT=1
        Y2OPT=0

C  Define the number and the set of tabulation points.

        NTAB=15
        DO 60 I=0,NTAB
            XTAB(I)=X(0)+I*(X(NP)-X(0))/NTAB
60      CONTINUE

C  Set SBOPT to 1 for selecting the sequential search algorithm.

        SBOPT=1

C  Call DBVSIS for computing the spline.

        CALL DBVSIS (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,DB4,DB4I,
     *               DDUMFN,DDUMFN,KMAX,MAXSTP,XTAB,NTAB,
     *               SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,D,D2,DIAGN,
     *               Y0TAB,Y1TAB,Y2TAB,WORK,NWORK)

C  Call DWRTO at the entry point DWRT4 for printing the results into
C  the output file.

        CALL DWRT4

C  ---------------------------------------------------------------------

C  Example n.5 .
C
C  Given the points { (0,0),(1,1),(2,2),(3,2),(4,1),(5,0) } we want to
C  construct a spline of degree 4 and continuity 1 which is co-convex
C  in the interval (2,3) and co-monotone elsewhere. We also require the
C  separable boundary conditions  s'(0)=1 , s'(5)=-1. The derivatives
C  are computed as the best approximation to Bessel estimates.
C  The spline is evaluated at four random points.
C  The routine DBVSIS is used.
C
C  Compute the input data

        NP=5
        X(0)=0.0D+00
        X(1)=1.0D+00
        X(2)=2.0D+00
        X(3)=3.0D+00
        X(4)=4.0D+00
        X(5)=5.0D+00
        Y(0)=0.0D+00
        Y(1)=1.0D+00
        Y(2)=2.0D+00
        Y(3)=2.0D+00
        Y(4)=1.0D+00
        Y(5)=0.0D+00

C  Set the degree and the class of continuity of the spline.

        N=4
        K=1

C  Set the size of the work area

        NWORK=PCOMM+(PPART+7)*NP+(N*(N+11))/2+9

C  Set the options on the spline: Bessel like derivatives,
C  separable boundary conditions, user specified shape constraints.

        OPT=234

C  Set the boundary conditions.

        D0=1.0D+00
        DNP=-1.0D+00

C  Select the shape constraints.

        CONSTR(0)=1
        CONSTR(1)=1
        CONSTR(2)=2
        CONSTR(3)=1
        CONSTR(4)=1

C  Call DWRTO at the entry point DWRT5A for printing the input
C  constraints values.

        CALL DWRT5A

C  Set EPS, KMAX and MAXSTP for the automatic choice of their values.

        EPS=0.
        KMAX=0
        MAXSTP=0

C  Set the options for evaluating the spline and its first derivative
C  at the tabulation points.

        Y0OPT=1
        Y1OPT=0
        Y2OPT=0

C  Define the number and the set of tabulation points.

        NTAB=3
        XTAB(0)=0.127D+01
        XTAB(1)=-0.5D+00
        XTAB(2)=0.57D+01
        XTAB(3)=0.35D+01

C  Set SBOPT to 2 for selecting the binary search algorithm.

        SBOPT=2

C  Call DBVSIS for computing the spline.

        CALL DBVSIS (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,DDUMFN,
     *               DDUMFN,DDUMFN,DDUMFN,KMAX,MAXSTP,XTAB,NTAB,
     *               SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,D,D2,DIAGN,
     *               Y0TAB,Y1TAB,Y2TAB,WORK,NWORK)

C  Call DWRTO at the entry point DWRT5B for printing the results into
C  the output file.

        CALL DWRT5B

        STOP
        END

C  ---------------------------------------------------------------------

        SUBROUTINE DWRTO

C  DWRTO is an auxiliary routine which prints the results of the driver
C  program DBVSDR in the external file RES. All the parameters
C  are shared with the main program in the common area and their meaning
C  is explained in DBVSDR.
C
C  Items of possible interest are:
C
C  OUN   : symbolic name of an integer constant, used for identifying
C          the logical output number;
C  FNAME : symbolic name of a character constant, used for identifying
C          the name of the output file.
C
C  The subroutine is called at the entry points DWRTO, DWRT1A, DWRT1B,
C  DWRT2, DWRT3A, DWRT3B, DWRT4, DWRT5A, DWRT5B.

        INTEGER PNP,PNTAB,PN
        PARAMETER (PNP=10,PNTAB=15,PN=6)

        INTEGER OUN
        INTEGER I
        CHARACTER FNAME*12
        PARAMETER (OUN=10,FNAME='RES')

        INTEGER NP,NTAB,CONSTR(0:PNP-1),ERRC,DIAGN(0:PNP-1),ERRE

        DOUBLE PRECISION X(0:PNP),Y(0:PNP),XTAB(0:PNTAB),Y0TAB(0:PNTAB),
     *                   Y1TAB(0:PNTAB),Y2TAB(0:PNTAB),D(0:PNP),
     *                   D2(0:PNP)

        COMMON /DBVCOM/ NP,NTAB,ERRC,ERRE,X,Y,XTAB,Y0TAB,Y1TAB,Y2TAB,
     *                  D,D2,DIAGN,CONSTR

C  Connect the output file FNAME to the logical unit number OUN.

        OPEN (UNIT=OUN,FILE=FNAME,STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL',FORM='FORMATTED')

C  Write the heading of the output file

        WRITE(UNIT=OUN,FMT=10)

        RETURN

C  Write the results of the example 1.

        ENTRY DWRT1A
        WRITE (UNIT=OUN,FMT=20) 'EXAMPLE  N. 1'
        WRITE (UNIT=OUN,FMT=30) ERRC,ERRE
        WRITE (UNIT=OUN,FMT=33)
        WRITE (UNIT=OUN,FMT=36) (I,DIAGN(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=40)
        WRITE (UNIT=OUN,FMT=50) (I,D(I),I=0,NP)
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=100) (I,XTAB(I),Y0TAB(I),I=0,NTAB)
        RETURN

        ENTRY DWRT1B
        WRITE (UNIT=OUN,FMT=25) 'FURTHER EVALUATIONS'
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=011) (I,XTAB(I),Y1TAB(I),Y2TAB(I),I=0,NTAB)
        RETURN

C  Write the results of the example 2.

        ENTRY DWRT2
        WRITE (UNIT=OUN,FMT=20) 'EXAMPLE  N. 2'
        WRITE (UNIT=OUN,FMT=30) ERRC,ERRE
        WRITE (UNIT=OUN,FMT=33)
        WRITE (UNIT=OUN,FMT=36) (I,DIAGN(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=45)
        WRITE (UNIT=OUN,FMT=55) (I,D(I),D2(I),I=0,NP)
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=101) (I,XTAB(I),Y0TAB(I),Y2TAB(I),I=0,NTAB)
        RETURN

C  Write the results of the example 3.

        ENTRY DWRT3A
        WRITE (UNIT=OUN,FMT=20) 'EXAMPLE  N. 3'
        WRITE (UNIT=OUN,FMT=30) ERRC,ERRE
        WRITE (UNIT=OUN,FMT=33)
        WRITE (UNIT=OUN,FMT=36) (I,DIAGN(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=45)
        WRITE (UNIT=OUN,FMT=55) (I,D(I),D2(I),I=0,NP)
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=110) (I,XTAB(I),Y0TAB(I),Y1TAB(I),I=0,NTAB)
        RETURN

        ENTRY DWRT3B
        WRITE (UNIT=OUN,FMT=25) 'FURTHER EVALUATIONS'
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=001) (I,XTAB(I),Y2TAB(I),I=0,NTAB)
        RETURN

C  Write the results of the example 4.

        ENTRY DWRT4
        WRITE (UNIT=OUN,FMT=20) 'EXAMPLE  N. 4'
        WRITE (UNIT=OUN,FMT=30) ERRC,ERRE
        WRITE (UNIT=OUN,FMT=33)
        WRITE (UNIT=OUN,FMT=36) (I,DIAGN(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=40)
        WRITE (UNIT=OUN,FMT=50) (I,D(I),I=0,NP)
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=110) (I,XTAB(I),Y0TAB(I),Y1TAB(I),I=0,NTAB)
        RETURN

C  Write the results of the example 5.

        ENTRY DWRT5A
        WRITE (UNIT=OUN,FMT=20) 'EXAMPLE  N. 5'
        WRITE (UNIT=OUN,FMT=70)
        WRITE (UNIT=OUN,FMT=90) (I,CONSTR(I),I=0,NP-1)
        RETURN

        ENTRY DWRT5B
        WRITE (UNIT=OUN,FMT=80)
        WRITE (UNIT=OUN,FMT=90) (I,CONSTR(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=30) ERRC,ERRE
        WRITE (UNIT=OUN,FMT=33)
        WRITE (UNIT=OUN,FMT=36) (I,DIAGN(I),I=0,NP-1)
        WRITE (UNIT=OUN,FMT=40)
        WRITE (UNIT=OUN,FMT=50) (I,D(I),I=0,NP)
        WRITE (UNIT=OUN,FMT=60)
        WRITE (UNIT=OUN,FMT=100) (I,XTAB(I),Y0TAB(I),I=0,NTAB)

C  Disconnect the output file.

        CLOSE (UNIT=OUN,STATUS='KEEP')

        RETURN

10      FORMAT(//22X,'DRESULTS.OUT'//
     *  '  THIS FILE CONTAINS THE COMPUTED RESULTS OF THE DRIVER'/
     *  '  PROGRAM DBVSDR, WHICH DEMONSTRATES THE USE OF THE'/
     *  '  SUBROUTINES DBVSIS, DBVSSC AND DBVSSE WITH SOME EXAMPLES')
20      FORMAT(//22X,A13/)
25      FORMAT(/22X,A19/)
30      FORMAT(/'  VALUES OF THE ERROR FLAGS :   ERRC=',
     *          I2,' ; ERRE=',I2)
33      FORMAT(/'  DIAGNOSTIC INFORMATION'//,
     *         10X,'I',7X,'DIAGN(I)'/)
36      FORMAT(9X,I2,11X,I1)
40      FORMAT(/'  VALUES OF THE FIRST DERIVATIVES AT THE KNOTS'//,
     *         10X,'I',11X,'D(I)'/)
45      FORMAT(/'  VALUES OF THE FIRST AND SECOND DERIVATIVES AT ',
     *         'THE KNOTS'//,10X,'I',11X,'D(I)',14X,'D2(I)'/)
50      FORMAT(9X,I2,5X,F13.9)
55      FORMAT(9X,I2,5X,F13.9,5X,F13.9)
60      FORMAT(/
     *  '  VALUES OF THE SPLINE AND ITS DERIVATIVES AT THE TABULATION',
     *  '  POINTS'//,
     *  7X,'I',8X,'XTAB(I)',8X,'Y0TAB(I)',8X,'Y1TAB(I)',8X,'Y2TAB(I)'/)
70      FORMAT(/'  VALUES OF THE INPUT CONSTRAINTS'//,
     *         10X,'I',7X,'CONSTR(I)'/)
80      FORMAT(/'  VALUES OF THE OUTPUT CONSTRAINTS'//,
     *         10X,'I',7X,'CONSTR(I)'/)
90      FORMAT(9X,I2,11X,I1)
001     FORMAT(5X,I3,3X,F13.9,35X,F13.9)
011     FORMAT(5X,I3,3X,F13.9,19X,F13.9,3X,F13.9)
100     FORMAT(5X,I3,3X,F13.9,3X,F13.9)
101     FORMAT(5X,I3,3X,F13.9,3X,F13.9,19X,F13.9)
110     FORMAT(5X,I3,3X,F13.9,3X,F13.9,3X,F13.9)

        END

C  ---------------------------------------------------------------------

C  The following functions define the non-separable boundary conditions.
C  For linking purposes, these must be given even when the conditions
C  are not required.
C  In this case the 'vacuous function' DDUMFN is used.

        DOUBLE PRECISION FUNCTION DB2(X)
        DOUBLE PRECISION X
        DB2=X
        RETURN
        END

        DOUBLE PRECISION FUNCTION DB2I(Y)
        DOUBLE PRECISION Y
        DB2I=Y
        RETURN
        END

        DOUBLE PRECISION FUNCTION DR2(X)
        DOUBLE PRECISION X
        DR2=X
        RETURN
        END

        DOUBLE PRECISION FUNCTION DR2I(Y)
        DOUBLE PRECISION Y
        DR2I=Y
        RETURN
        END

        DOUBLE PRECISION FUNCTION DB4(X)
        DOUBLE PRECISION X
        DB4=-X
        RETURN
        END

        DOUBLE PRECISION FUNCTION DB4I(Y)
        DOUBLE PRECISION Y
        DB4I=-Y
        RETURN
        END

        DOUBLE PRECISION FUNCTION DDUMFN(T)
        DOUBLE PRECISION T
        DDUMFN=0.0D+00*T
        RETURN
        END
