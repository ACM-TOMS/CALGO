        SUBROUTINE SBVSIS (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,
     *                     BETA,BETAI,RHO,RHOI,KMAX,MAXSTP,XTAB,NTAB,
     *                     SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,D,D2,
     *                     DIAGN,Y0TAB,Y1TAB,Y2TAB,WORK,NWORK)

C  SBVSIS is merely a support routine which calls SBVSSC and SBVSSE
C  for the computation of the needed parameters and for the evaluation
C  of a shape-preserving, C(k), k=1,2 , interpolating spline,
C  optionally subject to boundary conditions.
C  The use of SBVSIS is not recommended when more evaluations of the
C  same spline are required; in this case it is better to separately
C  call SBVSSC and then SBVSSE repeatedly.
C  For an explanation of input and output parameters, the user is
C  referred to the comments of SBVSSC and SBVSSE.


        EXTERNAL BETA,BETAI,RHO,RHOI

        INTEGER NP,N,K,OPT,CONSTR(0:NP-1),NTAB,KMAX,MAXSTP,SBOPT,
     *          Y0OPT,Y1OPT,Y2OPT,DIAGN(0:NP-1),ERRC,ERRE,NWORK

        REAL X(0:NP),Y(0:NP),D0,DNP,D20,D2NP,XTAB(0:NTAB),
     *       EPS,BETA,BETAI,RHO,RHOI,D(0:NP),D2(0:NP),
     *       Y0TAB(0:NTAB),Y1TAB(0:NTAB),Y2TAB(0:NTAB),
     *       WORK(1:NWORK)


        CALL SBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,BETA,
     *               BETAI,RHO,RHOI,KMAX,MAXSTP,ERRC,D,D2,DIAGN,
     *               WORK,NWORK)

        CALL SBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *               ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,
     *                     BETA,BETAI,RHO,RHOI,KMAX,MAXSTP,ERRC,D,D2,
     *                     DIAGN,WORK,NWORK)

C  -------------------------------------------------
C            Lines 48-548 are comment lines.
C            The actual code begins at line 554.
C  -------------------------------------------------

C  ABSTRACT:
C
C  SBVSSC is designed to compute the coefficients (first and, if
C  appropriate, second derivatives) of a shape-preserving spline, of
C  continuity class C(k), k=1,2 , which interpolates a set of data 
C  points and, if required, satisfies additional boundary conditions.
C  SBVSSC furnishes the input parameters for SBVSSE, which, in turn,
C  provides to evaluate the spline and its derivatives at a set of
C  tabulation points.
C
C  The user is allowed to use the following options:
C
C  - to compute a spline subject to:
C        - no constraint,
C        - monotonicity constraints,
C        - convexity constraints,
C        - monotonicity and convexity constraints,
C        - one of the above constraints in each subinterval;
C
C  - to impose separable or non-separable boundary conditions on the
C    spline,
C
C  - to assign the first derivatives d(i), i=0,1,...,np , in input or to
C    compute them from the constraints only or as the best approximation
C    to a set of optimal values. Although the final sequence of
C    derivatives does in any case satisfy the imposed restrictions on 
C    the shape, the resulting graphs may exhibit different behaviors.
C
C
C  REMARK:
C
C  In these comments variable and array names will be denoted with
C  capital letters, and their contents with small letters. Moreover:
C  .IN.   means belonging to;
C  .INT.  stands for intersection.
C
C
C  The code has the following structure:
C
C         SBVSSC
C              SBVC
C                   SSTINF
C                        SMSK1
C                        SMSK2
C                             SPRJ0
C                             SPRJ1
C                        STDC
C                             SMNMOD
C                             SMDIAN
C                   SALG3
C                        SPRJ0
C                        SALG1
C                             SPRJ1
C                        SINTRS
C                        STST
C                        SFPSVF
C                             SSL
C                             SALG1D
C                                  SPRJ1
C                        SAL2
C                             SPRJ2
C                        SAL2DP
C                             SMNIND
C                             SPRJ2
C                             SSL
C                        SSCDRC
C
C
C  CALLING SEQUENCE:
C
C       CALL SBVSSC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,EPS,BETA,
C    *               BETAI,RHO,RHOI,KMAX,MAXSTP,ERRC,D,D2,DIAGN,
C    *               WORK,NWORK)
C
C
C  INPUT PARAMETERS:
C
C  X       : floating array, of bounds 0:NP, containing the data
C            abscissas  x(i), i=0,1,...,np.
C            Restriction: x(i).LT.x(i+1), i=0,1,...,np.
C  Y       : floating array, of bounds 0:NP, containing the data
C            ordinates  y(i), i=0,1,...,np.
C  NP      : integer variable, defining the number of interpolation
C            points. Restriction: np.GE.2 .
C  N       : integer variable, containing the degree of s.
C            Restriction: n.GE.3 .
C  K       : integer variable, containing the class of continuity of s.
C            Restriction:  k.EQ.1  or  k.EQ.2  and  n.GE.3*k .
C  OPT     : integer variable, containing a control parameter. It is
C            a three-digit decimal of the form  pqr  (that is of
C            numerical value  p*100+q*10+r ) where:
C            r  controls the constraints on the shape. 
C            q  controls the boundary conditions and 
C            p  controls the computation of the derivatives, 
C            More specifically:
C            r=0 (opt=pq0) : no constraint on the shape is required;
C            r=1 (opt=pq1) : monotonicity constraints are required;
C            r=2 (opt=pq2) : convexity constraints are required;
C            r=3 (opt=pq3) : monotonicity and convexity constraints are
C                            required;
C            r=4 (opt=pq4) : local constraints for any subinterval are
C                            supplied by the user (see the description
C                            of the array CONSTR);
C            q=1 (opt=p1r) : no boundary condition is imposed;
C            q=2 (opt=p2r) : non-separable boundary conditions are 
C                            imposed (see the description of BETA,
C                            BETAI, RHO, RHOI);
C            q=3 (opt=p3r) : separable boundary conditions are imposed
C                            (see the description of D0, DNP, D20,
C                             D2NP);
C            p=1 (opt=1qr) : the sequence of first derivatives 
C                            d(0),....,d(np)  is computed using the 
C                            constraints only using subroutine SAL2;
C            p=2 (opt=2qr) : the sequence is computed as the constrained
C                            best approximation to Bessel derivatives
C                            using subroutine SAL2DP;
C            p=3 (opt=3qr) : the sequence is computed as the constrained
C                            best approximation to a set of third order
C                            accurate derivative estimates produced in
C                            subroutine STDC using subroutine SAL2DP 
C                            (since this estimates are inherently mono-
C                            tonicity preserving, it is not recommended
C                            to associate this option with the convexity
C                            constraints only);
C            p=4 (opt=4qr) : the sequence is computed as the constrained
C                            best approximation to a set of values given
C                            in input by the user using SAL2DP; note
C                            that opt.EQ.410 will provide the classical
C                            C(k) function interpolating the data and 
C                            the derivatives.
C         Restriction: ( p.GE.1 .AND. p.LE.4 ) .AND.
C                      ( q.GE.1.AND. q.LE.3 ) .AND.
C                      ( r.GE.0 .AND. r.LE.4 ) .AND.
C                      .NOT. ( r.EQ.0 .AND. p.EQ.1 ) .
C  D0      : floating variable containing the left separable boundary
C            condition for the first derivative (d(0)=d0).
C            D0 is referenced only when  q=3 .
C  DNP     : floating variable containing the right separable boundary
C            condition for the first derivative (d(np)=dnp).
C            DNP is referenced only when  q=3 .
C  D20     : floating variable containing the left separable boundary
C            condition for the second derivative (d2(0)=d20).
C            D20 is referenced only when  q=3  and  k=2 .
C  D2NP    : floating variable containing the right separable boundary
C            condition for the second derivative (d2(np)=d2np).
C            D2NP is referenced only when  q=3  and  k=2 .
C  EPS     : floating variable, containing a value for computing the
C            relative tolerance of the method. It should be set greater
C            or equal to the machine precision. However, if eps.LE.0, 
C            SBVSSC resets it to 0.0001 which has turned out to be a
C            good choice for graphical applications.
C  BETA    : user supplied function, which represents non-separable
C            boundary conditions for the first derivatives.
C            BETA is referenced only when  q=2 .
C  BETAI   : user supplied function, which is the inverse of BETA.
C            BETAI is referenced only when  q=2 .
C  RHO     : user supplied function, which represents non-separable
C            boundary conditions for the second derivatives.
C            RHO is referenced only when  q=2  and  k=2 .
C  RHOI    : user supplied function, which is the inverse of RHO.
C            RHOI is referenced only when  q=2  and  k=2 .
C  KMAX    : integer variable, containing the number of iterations
C            allowed for selecting the minimal set ASTAR described
C            below. If kmax.LE.0, SBVSSC resets it to 10 .
C            KMAX is referenced only when  q=2 .
C  MAXSTP  : integer variable, containing the number of steps allowed
C            to find dstar in the set of admissible values.
C            If maxstp.LE.0, SBVSSC resets it to 10 .
C            MAXSTP is referenced only when  q=2 .
C
C
C  INPUT / OUTPUT PARAMETERS:
C
C  CONSTR  : integer array, of bounds  0:NP , containing, in input the
C            desired constraints on the shape for each subinterval.
C            constr(i)=kind , kind=0,1,2,3 , means that none, monotoni-
C            city, convexity, monotonicity and convexity constraint is
C            imposed on the subinterval [x(i),x(i+1)]. If constr(i) is
C            not compatible with the data it is relaxed according to
C            their shape (see subroutine SMSK1 for details). So, on out-
C            put, CONSTR contains the constraints actually imposed.
C            For example, if the data are convex and on input we have
C            constr(i)=3 , the result in output will be  constr(i)=2.
C            Restriction: constr(i).GE.0 .AND. constr(i).LE.3 .
C            CONSTR is referenced only when  r=4 .
C  D       : floating array, of bounds 0:NP, containing the first
C            derivatives at  x(i), i=0,1,...,np . If  p=4 , d(i) is the
C            input value to be approximated by the computed derivative,
C            which is then stored in the same location.
C            On output, D is computed by the routine SAL2 if  p=1  and
C            is computed by the routine SAL2DP if  p=2  or  p=3 .
C
C
C  OUTPUT PARAMETERS
C
C  ERRC    : integer variable, containing an error flag which displays
C            the status of the code. The status is divided into: severe
C            error (error on the input data, no computation has been
C            done), error (some computation has been done and some
C            information or suggestions are available), warning (some
C            requirement is not fulfilled, but the spline's parameters
C            have been computed), success.
C            errc=0 : success, normal return of the code;
C            errc=1 : severe error, incorrect assignment for some of
C                     the values nwork, opt, np;
C            errc=2 : severe error, for some i the restriction
C                     0.LE.constr(i) .AND. constr(i).LE.3  is not 
C                     fulfilled;
C            errc=3 : severe error, incorrect assignment for some of 
C                     the values n,k;
C            errc=4 : severe error, the restriction x(i).LT.x(i+1) is 
C                     not fulfilled for some i;
C            errc=5 : error, the problem does not have any solution
C                     because the set
C                     betai ( phi(a(0,k)) .INT. beta(a(0,k)) )
C                     is empty for some k. In other words the boundary
C                     conditions cannot be satisfied and the output
C                     parameters are meaningless.
C                     The user is suggested to increase the value of n.
C            errc=6 : warning; for some i, the constraints on the
C                     interval  [x(i),x(i+1)]  are too strong and they 
C                     have not been considered. There is no guarantee
C                     that the spline is shape-preserving within all
C                     the intervals. More accurate diagnostic details
C                     can be found in the array DIAGN.
C                     The user is suggested to increase the value of n.
C            errc=7 : error, dstar such that beta(dstar).IN.phi(dstar)
C                     has not been found. The integer parameter maxstp 
C                     should be increased.
C                     The output parameters are meaningless.
C            errc=8 : error, both situations described in errc=6 and
C                     errc=7  have occurred.
C            errc=9 : warning, one of the separable boundary conditions
C                     d(0)=d0  and/or  d(np)=dnp  are not compatible
C                     with the constraints in  [x(0),x(1)]  and/or
C                     [x(np-1),x(np)]  which have consequently been
C                     relaxed. The user is suggested to increase the 
C                     value of n. More accurate diagnostic details can
C                     be found in the array DIAGN.
C            errc=10: warning, both situations described for errc=6 and
C                     errc=9 have occurred.
C            errc=11: warning, one of the separable boundary conditions
C                     d2(0)=d20  and/or  d2(np)=d2np  are not compatible
C                     with the constraints in  [x(0),x(1)]  and/or
C                     [x(np-1),x(np)] . The boundary conditions have
C                     consequently been approximated. The user is
C                     suggested to increase the value of n.
C            errc=12: warning, both situations described for errc=6 and
C                     errc=11 have occurred.
C            errc=13: warning, both situations described for errc=9 and
C                     errc=11 have occurred.
C            errc=14: warning, both situations described for errc=10 and
C                     errc=11 have occurred.
C            errc=15: warning, the non-separable boundary conditions
C                     d2(np)=rho(d2(0))  are not compatible with the 
C                     constraints. The boundary conditions have
C                     consequently been approximated. The user is
C                     suggested to increase the value of n.
C            errc=16: warning, both situations described for errc=6 and
C                     errc=15 have occurred.
C            errc=17: warning, both situations described for errc=9 and
C                     errc=15 have occurred.
C            errc=18: warning, both situations described for errc=10 and
C                     errc=15 have occurred.
C  D2      : floating array of bounds 0:NP containing the second
C            derivatives at knots. D2 is computed in subroutine DCDERC .
C            D2 is referenced only when  k=2 .
C  DIAGN   : integer array of bounds 0:NP-1 containing further
C            diagnostic information:
C            diagn(i)=0 : the constraints in the interval [x(i),x(i+1)]
C                         have been satisfied;
C            diagn(i)=1 : the constraints in the interval [x(i),x(i+1)]
C                         have not been satisfied;
C            
C
C
C  OTHER PARAMETERS:
C
C  WORK    : floating array, of bounds 1:NKORK, which is used as
C            a work area to store intermediate results.
C            The same array can be used to provide workspace for both
C            the main subroutines  SBVSSC and SBVSSE .
C  NWORK   : integer variable containing the size of the work area.
C            Restriction: nwork .GE. comm+(part+7)*np+(n*(n+11))/2+9
C                           that is
C                         nwork .GE. 5+(2+7)*np+(n*(n+11))/2+9
C
C
C  ------------------------------------------------
C
C  METHOD:
C
C  Let the integers n and k, such that k=1,2 ; n >= 3k , and the
C  sequence of points  (x(i),y(i)), i=0,1,...,np , with 
C  x(0) < x(1) < ... <x(np) , be given; let us denote with  BS(n;k)
C  the set of the splines s of degree n and continuity k whose second
C  derivative, in the case k=2 , vanishes at the knots. We are 
C  interested in the existence and construction, if possible, of a 
C  shape-preserving interpolating spline s of BS(n;k) such that
C
C            s(x(i)) = y(i) , i=0,1,...,np                          (1)
C
C  and optionally subject to constraints on the shape in each interval 
C  [x(i),x(i+1)] .
C
C  In the case k=2 the zero derivatives of the spline  s.IN.BS(n;k) are
C  then modified to assume non-zero values which are not in contrast
C  with the shape constraints and, if possible, satisfy the boundary
C  conditions eventually imposed by the user. For details we refer to 
C  the comments in subroutine SSCDRC.
C
C  Since any s.IN.BS(n;k) is determined by its values and slopes at
C  x(i) , i=0,1,...,np , we can reformulate the problem as follows: 
C  compute the values  d(i), i=0,1,...,np , such that the spline s, 
C  satisfying (1) and
C
C            Ds(x(i)) := d(i) , i=0,1,...,np                        (2)
C
C  is shape-preserving.
C  Setting  delta(i) := (y(i+1)-y(i))/(x(i+1)-x(i)) , we have that s is
C  increasing (I) ( decreasing (D) ) in [x(i),x(i+1)] if and only if
C  (d(i),d(i+1))  belongs to
C
C    D(i) := { (u,v).IN.RxR : u >= 0, v >= 0, v =< -u+ n/k delta(i) }
C                                                                    (3)
C  ( D(i) := { (u,v).IN.RxR : u =< 0, v =< 0, v >= -u+ n/k delta(i) } )
C
C  s is convex (CVX) ( concave (CNC) ) if and only if (d(i),d(i+1))
C  belongs to
C
C    D(i) := { (u,v).IN.RxR : v >= - (k/(n-k)) u + (n/(n-k)) delta(i) ,
C                             v =< - ((n-k)/k) u + (n/k) delta(i) }
C                                                                    (4)
C  ( D(i) := { (u,v).IN.RxR : v =< - (k/(n-k)) u + (n/(n-k)) delta(i) ,
C                             v >= - ((n-k)/k) u + (n/k) delta(i) }  )
C
C  and that s is I (D) and CVX (CNC) if and only if (d(i),d(i+1)) 
C  belongs to
C
C             D(i) := { (u,v) satisfying (3) and (4) } .
C
C  So, if we choose the family of sets D(i) , i=0,1,...,np-1 , according
C  to the shape of the data, we have to solve:
C
C  PROBLEM P1. Does a sequence ( d(0), d(1), ..., d(np) ) such that
C              (d(i),d(i+1)) .IN. D(i) , i=0,1,...,np-1 , exist ?
C
C  PROBLEM P2. If P1 is feasible, how can a (the best) solution be 
C              computed ?
C
C  Let SPRJ1: RxR -> R and SPRJ2: RxR -> R be, respectively, the 
C  projection maps from uv-plane onto the u-axis and v-axis and let us 
C  denote with  B(i) := SPRJ1(D(i)) :
C
C      ALGORITHM A1[B0].
C        1. Set A(0):=B(0); J:=np.
C        2. For i=1,2,...,np
C           2.1. Set A(i):= SPRJ2( D(i-1).INT.{ A(i-1) x B(i) } ) .
C           2.2. If A(i) is empty, set J:=i and stop.
C        3. Stop.
C
C  We have the following result:
C
C  THEOREM 1. P1 has a solution if, and only if, J=np, that is A(i) is
C             not empty , i=0,1,...,np . If ( d(0), d(1), ...,d(np) )
C             is a solution then  d(i).IN.A(i) , i=0,1,...,np .
C
C  A solution can be computed with the following algorithm:
C
C      ALGORITHM A2[A(np),B0].
C        1. Choose d(np).IN.A(np).
C        2. For i=np-1, np-2, ..., 0
C           2.1. Choose d(i).IN.SPRJ1( D(i).INT.{ A(i) x { d(i+1) }}).
C        3. Stop.
C
C  For more theoretical details about A1 and A2 see \1\ , and for
C  practical details see subprograms SALG1, SAL2, SAL2DP. In the latter
C  a dynamic programming scheme is used to find the best solution in 
C  the feasible set. More specifically, it is possible to compute the
C  values  d(i),i=0,..,np which satisfy the constraints and are as close 
C  as possible to another sequence which does not satisfy the 
C  constraints but is, in some sense, optimal.
C  
C  From a theoretical point of view, algs A1 and A2 give a complete 
C  answer to problems P1 and P2. However, it could be pointed out that,
C  for practical applications, we would like to have the best possible 
C  plot, whether or not P1 has solutions. Let us suppose that the 
C  problem is solvable from 0 to j and from j to np, but that alg A1 
C  applied to the whole family of sets  D(i), i=0,1,...,np-1  gives
C  J.eq.j.ne.np ; if we reset  D(j-1) := A(j-1) x B(j) , alg A1 applied
C  to this new family of sets will produce J=np . However, it must be
C  recalled that, in this way, we do not consider the constraints in the
C  interval [x(j-i),x(j)] and so there is no guarantee that the spline
C  is shape-preserving in this interval. Whenever this fact cannot be
C  accepted it is convenient to rerun the code with a larger value for
C  the degree n , because the domains of constraints enlarge as n 
C  increases (see (3) and (4)).
C
C  It is immediate to see that separable boundary conditions of the form
C
C            d(0) := d0 ; d(np) := dnp
C
C  can be easily inserted with a reduction of the corresponding
C  admissible sets which does not modify the above theory:
C
C       D(0) := D(0).INT.{d(0)=d0} ; D(np) := D(np).INT.{d(np)=dnp}
C
C  In the case k=2 the corresponding conditions  d2(0) = d20 ,
C  d2(np) = d2np  are imposed only if not in contrast with the shape of
C  the data; otherwise the admissible values for  d2(0) and d2(np) 
C  respectively closest to d20 and d2np are chosen.
C
C  Now, let beta be a continuous function from R to R, with continuous
C  inverse betai, we want to solve the following non-separable boundary
C  valued problem:
C
C  PROBLEM P3. Do sequences ( d(0), d(1), ..., d(np) ) , such that
C              (d(i),d(i+1)).IN.D(i), i=0,1,...,np-1    and
C              d(np) = beta ( d(0) ) , exist ?
C
C  It is obvious that a solution of this new problem, if it exists, can
C  be found amongst the solutions of P1. Let A(0), A(1),...,A(np) be the
C  sequence of sets given by alg A1 (we assume that A(i) is not empty,
C  i=0,1,...,np , that is P1 is solvable or, if this is not the case,
C  the constraints have been relaxed ), we can assume that 
C  A(np) = phi(A(0)) , where  phi: R -> R is a set valued function 
C  (see \1\ for details). It can be demonstrated that:
C
C  THEOREM 2. P1 is solvable if, and only if, there is  dstar.IN.A(0)
C             such that   beta(dstar).IN.phi({dstar}) .
C
C  It should be noted that if ( d(0), d(1), ..., d(np) ) satisfies P1,
C       d(0) .IN. betai(phi(A(0)).INT.beta(A(0))) =: gamma(A(0))
C  and, consequently, the set of admissible values is reduced. If we 
C  repeat this procedure, we get a gradually diminishing admissible set
C  for d(0). We define
C     ASTAR := lim A(0,m)  where
C     A(0,0) := A(0)   and   A(0,m) := gamma(A(0,m-1)) ;
C  ASTAR is the minimal admissible set for dstar. We can now combine the
C  various theorems and algorithms and give the general algorithm to
C  solve P3:
C
C      ALGORITHM A3.
C        1. Set A(0,0) := B0 ; m:=1.
C        2. Use A1[A(0,0)] for computing phi (A(0,0)).
C        3. Set A(0,1) := gamma(A(0,0))
C                       = betai(phi(A(0,0)).INT.beta(A(0,0))).
C        4. If A(0,1) is empty, stop (P1 is unsolvable).
C        5. While ( convergence test not satisfied ) do
C           5.1. Use A1[A(0,m)] for computing A(np,m) = phi (A(0,m)).
C           5.2. Set A(0,m+1) := gamma(A(0,m)).
C           5.3. Set m:=m+1.
C        6. Set ASTAR := A(0,m).
C        7. Use A1[{d(0)}] to find dstar.IN.ASTAR such that
C           beta(dstar).IN.phi(dstar).
C        8. Use A2[beta(dstar),dstar] for computing a sequence
C           ( d(0), d(1), ..., d(np) )  which solves P1.
C
C  In the case k=2 the corresponding condition  d2(np) = beta2(d2(0))
C  is imposed only if not in contrast with the shape of
C  the data; otherwise the admissible values for  d2(0) and d2(np) 
C  closest to the boundary condition are chosen.
C
C  References
C
C  \1\ P.Costantini: Boundary Valued Shape-Preserving Interpolating 
C      Splines, ACM Trans. on Math. Softw., companion paper.
C  \2\ R.Bellman, S.Dreyfus: Applied Dynamic Programming, Princeton
C      University Press, New York, 1962.
C  \3\ H.T.Huynh: Accurate Monotone Cubic Interpolation, SIAM J. Num.
C      Anal., 30 (1993), 57-100.
C
C  The ideas involved in Algorithm A3 have been implemented in the code
C  in a general form. Since Algorithm A3 resembles closely the abstract
C  formulation it could, therefore, be used for several practical
C  problems. The particular case actually treated is reflected in the
C  contents of the information array INFO (see its description in 
C  subroutine SSTINF) which contains all the data needed for the
C  construction of the operators SPRJ0, SPRJ1 and SPRJ2.
C
C  As a consequence, the user has the following options:
C
C  - to compute a Spline subject to:
C        - no constraint;
C        - monotonicity constraints,
C        - convexity constraints,
C        - monotonicity and convexity constraints,
C        - one of the above constraints in each subinterval, as 
C          specified in the corresponding array CONSTR;
C
C  - to impose separable or non-separable boundary conditions on the
C    spline. In the latter case, the external functions BETA, BETAI,
C    RHO and RHOI must be supplied,
C
C  - to assign the first derivatives d(i), i=0,1,...,np , in input or to
C    compute them from the only constraints or as the best approximation
C    to a set of optimal values. Although the final sequence of
C    derivatives does in any case satisfy the imposed restrictions on 
C    the shape, the resulting graphs may exhibit different behaviors.
C
C  See the description of the input parameter OPT for more details.

C  ------------------------------------------------
C            End of comments.
C  ------------------------------------------------

        INTEGER COMM,PART

C  COMM contains the number of global data referring to the initial
C  points  (x(i),y(i)) stored in the array INFO, described in
C  subroutine SSTINF.
C  PART contains the number of particular data referring to each
C  interval  (x(i),x(i+1)) , i=0,1,...,np , stored in the array INFO.

        PARAMETER (COMM=5, PART=2)

        EXTERNAL BETA,BETAI,RHO,RHOI

        INTEGER NP,N,K,OPT,CONSTR(0:NP-1),KMAX,MAXSTP,ERRC,
     *          DIAGN(0:NP-1),NWORK,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10

        REAL X(0:NP),Y(0:NP),D0,DNP,D20,D2NP,EPS,BETA,BETAI,
     *       RHO,RHOI,D(0:NP),D2(0:NP),WORK(1:NWORK)


C  Assign the success value to the error flag.

        ERRC=0

C  Check the size of the work area.

        IF ( NWORK.LT.COMM+(PART+7)*NP+(N*(N+11))/2+9 ) THEN
            ERRC=1
            RETURN
        END IF

C  Compute indices for the splitting of the work array WORK.

        I1=1
        I2=I1+(COMM+PART*NP+NP+1)
        I3=I2+2*(NP+1)
        I4=I3+2*(NP+1)
        I5=I4+(N*(N+1)/2+N)
        I6=I5+(N+1)
        I7=I6+(N+1)
        I8=I7+(N+1)
        I9=I8+(N+1)
        I10=I9+NP-1

C  SBVSSC is essentially an interfacing routine which relieves the
C  user of a longer calling sequence. The structure of the method can 
C  be seen in SBVC and in the subroutines called.

        CALL SBVC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,WORK(I1),
     *             COMM,PART,EPS,KMAX,MAXSTP,WORK(I2),WORK(I3),
     *             WORK(I9),WORK(I10),BETA,BETAI,RHO,RHOI,
     *             D,D2,ERRC,DIAGN)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *                     ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)

C  -------------------------------------------------
C            Lines 619-752 are comment lines.
C            The actual code begins at line 758.
C  -------------------------------------------------

C  ABSTRACT:
C
C  SBVSSE is designed to evaluate the interpolating, shape-preserving
C  spline computed in subroutine SBVSSC.
C
C
C  REMARK:
C
C  In these comments variable and array names will be denoted with
C  capital letters, and with small letters their contents.
C
C
C  METHOD:
C
C  Let a spline  s:=s(x)  of degree n and continuity k (k=1,2) ,
C  interpolating at the knots the point (x(i),y(i)) , i=0,1,...,np ,
C  be previously computed in subroutine SBVSSC. Then, given a set of
C  tabulation points  xtab(i) , i=0,1,...,ntab , SBVSSE computes the
C  values  y0tab(itab):=s(xtab(itab))  and/or 
C  y1tab(itab):=Ds(xtab(itab))  and/or  y2tab(itab):=DDs(xtab(itab)) ,
C  using, under user selection, a sequential or binary search scheme.
C
C  The code has the following structure:
C
C         SBVSSE
C             SBVE
C                 STRMB
C                 SSQ
C                     SLSPIS
C                     SBL
C                     SBL1
C                     SBL2
C                 SBNTAB
C                     SBSEAR
C                     SLSPIS
C                     SBL
C                     SBL1
C                     SBL2
C
C
C  CALLING SEQUENCE:
C
C        CALL SBVSSE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
C    *                ERRC,D,D2,Y0TAB,Y1TAB,Y2TAB,ERRE,WORK,NWORK)
C
C
C  INPUT PARAMETERS:
C
C  X       : floating array, of bounds 0:NP, containing the data
C            abscissas  x(i), i=0,1,...,np. Restriction:
C            x(i).LT.x(i+1), i=0,1,...,np , checked in SBVSSC.
C  Y       : floating array, of bounds 0:NP, containing the data
C            ordinates  y(i), i=0,1,...,np.
C  NP      : integer variable, defining the number of interpolation
C            points. Restriction: np.GE.2 , checked in SBVSSC.
C  N       : integer variable, containing the degree of s.
C            Restriction: n.GE.3 , checked in SBVSSC
C  K       : integer variable, containing the class of continuity of s.
C            Restriction:  k.EQ.1  or  k.EQ.2  and  n.GE.3*k , checked
C            in SBVSSC.
C  XTAB    : floating array, of bounds 0:NTAB, containing the abscissas
C            of tabulation points.
C            Restriction: xtab(i).LE.xtab(i+1), i=0,1,...,ntab-1 .
C  NTAB    : integer variable, defining the number of tabulation points.
C            Restriction: ntab.GE.0 .
C  SBOPT   : integer variable, containing a control parameter.
C            If sbopt=1 then the sequential search is used for selecting
C            the interval of interpolation points in which xtab(i)
C            falls. If sbopt=2, binary search is used.
C            Restriction: sbopt.EQ.1 .OR. sbopt.EQ.2 .
C  Y0OPT   : integer variable, containing a control parameter. 
C            If y0opt=1, the spline is evaluated at the points
C            xtab(i), i=0,1,...,ntab and the results are stored at the
C            array  Y0TAB.
C            Restriction: y0opt.EQ.0 .OR. y0opt.EQ.1 .
C  Y1OPT   : integer variable, containing a control parameter.
C            If y1opt=1 the first derivatives of the spline at points
C            xtab(i) i=0,1,...,ntab , are computed and the results are
C            stored in the array Y1TAB .
C            Restriction: y1opt.EQ.0 .OR. y1opt.EQ.1 .
C  Y2OPT   : integer variable, containing a control parameter. 
C            If y2opt=1 the second derivatives of the spline at points 
C            xtab(i), i=0,1,...,ntab  are computed and the results are
C            stored in the array Y2TAB.
C            Restriction: y2opt.EQ.0 .OR. y2opt.EQ.1 .
C  ERRC    : integer variable, containing the status of the last
C            execution of subroutine SBVSSC.
C  D       : floating array, of bounds 0:NP, containing the first
C            derivatives at the knots.
C  D2      : floating array of bounds 0:NP containing the second
C            derivatives at the knots.
C
C
C  OUTPUT PARAMETERS:
C
C
C  Y0TAB   : floating array, of bounds 0:NTAB, containing the values of
C            the spline at the tabulation points xtab(i) ,
C            i=0,1,...,ntab when the option  y0opt=1  is activated.
C  Y1TAB   : floating array, of bounds 0:NTAB, containing the values of 
C            the first derivative of the spline at the tabulation points
C            xtab(i) , i=0,1,...ntab , when the option y1opt=1 is 
C            activated.
C  Y2TAB   : floating array, of bounds 0:NTAB, containing the values of 
C            the second derivative of the spline at the tabulation 
C            points xtab(i) , i=0,1,...,ntab , when the option y2opt=1
C            is activated.
C  ERRE    : integer variable, containing an error flag which displays
C            the status of the code. SBVSSE has only two levels of error
C            (see SBVSSC for comparison): success and severe error,
C            which means that some incorrect assignment for input data
C            have been set.
C            ERRE=0:  success, normal return of the code;
C            ERRE=1:  severe error, the value errc gives a status of 
C                     error, which means that the output of SBVSSC is 
C                     meaningless. Check the input parameters of SBVSSC.
C            ERRE=2:  severe error, incorrect assignment for some of
C                     the values ntab, sbopt, y0opt, y1opt, y2opt ,
C                     nwork;
C            ERRE=3:  severe error, the restriction xtab(i).LT.xtab(i+1)
C                     is not fulfilled for some i when sequential search
C                     is required;
C 
C
C  OTHER PARAMETERS:
C
C  WORK    : floating array, of bounds 1:NKORK, which is used as
C            a work area to store intermediate results.
C            The same array can be used to provide workspace for both
C            the main subroutines  SBVSSC and SBVSSE .
C  NWORK   : integer variable containing the size of the work area.
C            Restriction: nwork .GE. comm+(part+7)*np+(n*(n+11))/2+9
C                           that is
C                         nwork .GE. 3+(2+7)*np+(n*(n+11))/2+9

C  -------------------------------------------------
C            End of comments.
C  -------------------------------------------------

        INTEGER COMM,PART
        PARAMETER (COMM=5, PART=2)

        INTEGER NP,N,K,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,NWORK,
     *          I1,I2,I3,I4,I5,I6,I7,I8,I9,I10

        REAL X(0:NP),Y(0:NP),XTAB(0:NTAB),Y0TAB(0:NTAB),
     *       Y1TAB(0:NTAB),Y2TAB(0:NTAB),D(0:NP),D2(0:NP),
     *       WORK(1:NWORK)


C  Assign the success value to the error flag.

        ERRE=0

C  Check the size of the work area.

        IF ( NWORK.LT.COMM+(PART+7)*NP+(N*(N+11))/2+9 ) THEN
            ERRE=2
            RETURN
        END IF

C  Compute indices for the splitting of the work array WORK.

        I1=1
        I2=I1+(COMM+PART*NP+NP+1)
        I3=I2+2*(NP+1)
        I4=I3+2*(NP+1)
        I5=I4+(N*(N+1)/2+N)
        I6=I5+(N+1)
        I7=I6+(N+1)
        I8=I7+(N+1)
        I9=I8+(N+1)
        I10=I9+NP-1

C  SBVSSE is essentially an interfacing routine which relieves the
C  user of a longer calling sequence. The structure of the method can
C  be seen in SBVE and in the subroutines called.

        CALL SBVE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
     *             D,D2,ERRC,WORK(I4),WORK(I5),WORK(I6),
     *             WORK(I7),WORK(I8),Y0TAB,Y1TAB,Y2TAB,ERRE)

        RETURN
        END


        SUBROUTINE SALG1 (A1,NP,INFO,COMM,PART,EPS,A2,ERRC,DIAGN)

C  SALG1 implements the algorithm A1[B(0)] described in subr. DBVSSC.
C
C  The input parameters NP,COMM,PART,EPS and the output parameters
C  ERRC, DIAGN are described in DBVSSC. The input parameter INFO is 
C  described in SSTINF.
C
C  Items of possible interest are:
C
C  A1: floating array, of bounds 1:2, 0:NP, containing the sequence of
C      the sets  B(i), i=0,1,...,np (see the comments in DBVSSC).
C      More precisely,  B(i) = [a1(1,i),a1(2,i)] .
C
C  A2: floating array, of bounds 1:2, 0:NP, containing the sequence of
C      the sets  A(i), i=0,1,...,np (see the comments in DBVSSC).
C      More precisely, A(i) = [a2(1,i),a2(2,i)] .


        INTEGER NP,COMM,PART,ERRC,DIAGN(0:NP-1),ERRC1,I

        REAL A1(1:2,0:NP),INFO(1:COMM+PART*NP+NP+1),EPS,
     *       A2(1:2,0:NP)

        SAVE FL0
        REAL FL0
        DATA FL0/0.0E00/

        ERRC1=0

C  Step 1.

        A2(1,0)=A1(1,0)
        A2(2,0)=A1(2,0)

C  Step 2.

        DO 20 I=1,NP
10          CONTINUE

C  Step 2.1.

                CALL SPRJ1 (A2(1,I-1),A2(2,I-1),A1(1,I),A1(2,I),I,INFO,
     *                      COMM,PART,NP,A2(1,I),A2(2,I))

C  Ignore the constraints in  [x(i),x(i+1)]  when A(i) is empty.

            IF (A2(1,I).GT.A2(2,I)+EPS) THEN
                INFO(COMM+1+(I-1))=FL0
                ERRC1=1
                DIAGN(I-1)=1
                GO TO 10
            END IF
20      CONTINUE

        IF (ERRC1.EQ.1 .AND. ERRC.EQ.9) THEN
            ERRC=10
        ELSE IF (ERRC1.EQ.1) THEN
            ERRC=6
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,A2,ERRC1)

C  SALG1D computes the sequence of sets A(i), i=0,1,...,np, implementing
C  the algorithm A1[{dstar}], that is with A(0)={dstar} (see the com-
C  ments in subroutine DBVSSC for details).
C
C  The input parameters NP,COMM,PART,EPS are described in DBVSSC; the
C  input parameter INFO is described in SSTINF; the input parameters A1
C  and A2 are described in subprogram SALG1.
C
C  Item of possible interest is:
C
C  ERRC1  : Integer parameter, containing a control variable which is
C           then used in subr. SFPSVF
C           errc1 = 0 - success, normal return of the subprogram;
C           errc1 = 1 - A(i) is empty for some i.


        INTEGER NP,COMM,PART,ERRC1,I

        REAL DSTAR,A1(1:2,0:NP),INFO(1:COMM+PART*NP+NP+1),
     *       EPS,A2(1:2,0:NP)

        ERRC1=0

C  Step 1.

        A2(1,0)=DSTAR
        A2(2,0)=DSTAR

C  Step 2.

        DO 10 I=1,NP
            CALL SPRJ1 (A2(1,I-1),A2(2,I-1),A1(1,I),A1(2,I),I,INFO,
     *                  COMM,PART,NP,A2(1,I),A2(2,I))
            IF (A2(1,I).GT.A2(2,I)+EPS) THEN
                ERRC1=1
                RETURN
            END IF
10      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SALG3 (INFO,NP,COMM,PART,OPT,D0,DNP,EPS,KMAX,MAXSTP,
     *                    BETA,BETAI,A1,A2,D,ERRC,DIAGN)

C  SALG3 computes a sequence of slopes ( d(0), d(1), ..., d(np) ) which
C  can be used to compute a shape-preserving interpolating spline with
C  or without boundary conditions, as requested by the user. It is an
C  implementation of the algorithm A3 described in subroutine DBVSSC.
C
C  The input parameters NP,COMM,PART,OPT,EPS,KMAX,MAXSTP,BETA,BETAI,D
C  and the output parameter ERRC are described in subprogram DBVSSC.
C  The input parameter INFO is described in subprogram SSTINF.


        EXTERNAL BETA,BETAI

        INTEGER NP,COMM,PART,OPT,KMAX,MAXSTP,ERRC,DIAGN(0:NP-1),
     *          I,K,P,Q

        REAL INFO(1:COMM+PART*NP+NP+1),D0,DNP,EPS,BETA,
     *       BETAI,A1(1:2,0:NP),A2(1:2,0:NP),D(0:NP),DSTAR,
     *       P1,P2

        LOGICAL STST

        INTEGER INK,INSTP
        PARAMETER (INK=10,INSTP=10)

        SAVE FL2
        REAL FL2
        DATA FL2/2.0E00/

        P=OPT/100
        Q=MOD(OPT,100)/10

C  If kmax.LE.0 it is reset to ink.

        IF(KMAX.LE.0) THEN
            KMAX=INK
        END IF

C  If maxstp.LE.0 it is reset to instp.

        IF(MAXSTP.LE.0) THEN
            MAXSTP=INSTP
        END IF

C  Start step 1: store the sets  B(i), i=0,1,...,np , into the array A1.

        DO 10 I=1,NP
            CALL SPRJ0 (I,INFO,COMM,PART,NP,A1(1,I-1),A1(2,I-1))
10      CONTINUE
        A1(1,NP)=INFO(4)
        A1(2,NP)=INFO(5)

C  Reset the first and the last interval if separable boundary condtions
C  are required

        IF (Q.EQ.3) THEN
            A1(1,0)=D0
            A1(2,0)=D0
            A1(1,NP)=DNP
            A1(2,NP)=DNP
        END IF

C  Start step 2. Call SALG1 to compute the array A2 containing the
C  sets A(i) , i=0,1,...,np.

        CALL SALG1 (A1,NP,INFO,COMM,PART,EPS,A2,ERRC,DIAGN)

C  Start step 3 (steps 3-7 are activated only if boundary conditions are
C  required).

        IF (Q.EQ.2) THEN

C  Compute  betai(phi(A(0).INT.beta(A(0))) .

            CALL SINTRS (MIN(BETA(A2(1,0)),BETA(A2(2,0))),
     *                   MAX(BETA(A2(1,0)),BETA(A2(2,0))),
     *                   A2(1,NP),A2(2,NP),P1,P2)
            A1(1,0)=MIN(BETAI(P1),BETAI(P2))
            A1(2,0)=MAX(BETAI(P1),BETAI(P2))

C  Start step 4.

            IF (P1.GT.P2+EPS) THEN
                ERRC=5
                RETURN
            END IF

C  Start step 5 : initialization

            K=1
            DO 20 I=1,NP
                A1(1,I)=A2(1,I)
                A1(2,I)=A2(2,I)
20          CONTINUE

C  Iteration. The loop is stopped if a convergence test is satisfied
C  or kmax iterations have already been done.

30          CONTINUE
            IF (STST(A1(1,0),A1(2,0),A2(1,0),A2(2,0),EPS).OR.
     *                                             K.EQ.KMAX)  GO TO 50

C  Step 5.1 .

                CALL SALG1 (A1,NP,INFO,COMM,PART,EPS,A2,ERRC,DIAGN)

C  Step 5.2 .

                CALL SINTRS (MIN(BETA(A2(1,0)),BETA(A2(2,0))),
     *                       MAX(BETA(A2(1,0)),BETA(A2(2,0))),
     *                       A2(1,NP),A2(2,NP),P1,P2)
                A1(1,0)=MIN(BETAI(P1),BETAI(P2))
                A1(2,0)=MAX(BETAI(P1),BETAI(P2))

C  If  gamma(A(0))  is empty for some k the problem does not have any
C  solution.

                IF (P1.GT.P2+EPS) THEN
                    ERRC=5
                    RETURN
                END IF

                DO 40 I=1,NP
                    A1(1,I)=A2(1,I)
                    A1(2,I)=A2(2,I)
40              CONTINUE

                K=K+1
                GO TO 30
50          CONTINUE

C  Start step 7.
C  Assign to dstar a suitable value

            IF(P.EQ.1) THEN
                DSTAR=(A1(1,0)+A1(2,0))/FL2
            ELSE
                DSTAR=INFO(COMM+PART*NP+1)
            END IF

C  Check if dstar solves the problem, that is,  beta(dstar)  belongs to
C  phi(dstar); if it is not the case, another value for dstar
C  is looked for.

            CALL SFPSVF (A1,A2,NP,INFO,COMM,PART,EPS,MAXSTP,
     *                   BETA,ERRC,DSTAR)
            IF(ERRC.NE.0.AND.ERRC.NE.6.AND.ERRC.NE.9
     *               .AND.ERRC.NE.10) RETURN
            INFO(COMM+PART*NP+NP+1)=BETA(DSTAR)
            A2(1,NP)=BETA(DSTAR)
            A2(2,NP)=BETA(DSTAR)

        END IF

C  Start step 8.

        IF(P.EQ.1) THEN
            CALL SAL2 (A2,NP,INFO,COMM,PART,D)
        ELSE
            CALL SAL2DP (A2,NP,INFO,COMM,PART,D)
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SAL2(A2,NP,INFO,COMM,PART,D)

C  SAL2 computes a sequence of slopes (d(0),d(1),...,d(np)) implementing
C  alg. A2  described in subr. DBVSSC. Each d(i),i=0,1,...,np , is
C  chosen as the midpoint of the interval of all feasible values .
C
C  The input parameters NP,COMM,PART and the output parameter D are
C  described in DBVSSC; the input parameter INFO is described in SSTINF.
C
C  Item of possible interest is:
C
C  A2   : floating array, of bounds 1:2, 0:NP; [a2(1,i),a2(2,i)]
C         is the feasible interval for d(i) .


        INTEGER NP,COMM,PART,I

        REAL A2(1:2,0:NP),INFO(1:COMM+PART*NP+NP+1),
     *       D(0:NP),P1,P2

        SAVE FL1D2
        REAL FL1D2
        DATA FL1D2/0.5E00/

        D(NP)=(A2(1,NP)+A2(2,NP))*FL1D2

        DO 10 I=NP,1,-1
            CALL SPRJ2(A2(1,I-1),A2(2,I-1),D(I),D(I),I,INFO,
     *                       COMM,PART,NP,P1,P2)
            D(I-1)=(P1+P2)*FL1D2
10      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SAL2DP (A2,NP,INFO,COMM,PART,D)

C  SAL2DP links algorithm A2 and a dynamic programming scheme
C  to select, among the set of all feasible solutions, the sequence
C  ( d(0),d(1), ..., d(np) ) which is the best 2-norm approximation to
C  a set of optimal values. More precisely, if (ds(0),ds(1), ...,ds(np))
C  is the sequence of optimal values, SAL2DP use the following dynamic
C  programming relations
C
C    psi(0;d(0)) := (d(0)-ds(0))**2
C    psi(j;d(j)) := (d(j)-ds(j))**2 + min(psi(j-1;d(j-1)))
C
C  for describing the objective function
C
C      SUM  ((d(j) - ds(j)) ** 2
C    j=0,np
C
C  For a complete comprehension of the algorithm see the book \2\
C  quoted in the references of subr. DBVSSC
C
C  The input parameters NP,COMM,PART and the output parameter D are
C  described in subprogram DBVSSC; the input parameter INFO is described
C  in subprogram SSTINF and the input parameter A2 is described in SAL2.
C  The constant NSUBD defined below is related to the discretization of
C  the admissible domain.


        INTEGER NSUBD
        PARAMETER (NSUBD=20)

        INTEGER NP,COMM,PART,IND,I,J,JD0,SMNIND

        REAL A2(1:2,0:NP),INFO(1:COMM+PART*NP+NP+1),D(0:NP),
     *       PSI0(0:NSUBD+1),PSI1(0:NSUBD+1),
     *       PART0(0:NSUBD+1),PART1(0:NSUBD+1),H0,H1,PSI1MN,
     *       P1,P2,D0,D1,SSL

        SAVE FL0,FLE30
        REAL FL0,FLE30
        DATA FL0,FLE30/0.0E00,1.0E30/

        IND=COMM+PART*NP+1

        H0=MAX(FL0,(A2(2,0)-A2(1,0))/NSUBD)
        DO 5 J=0,NSUBD
            PART0(J)=A2(1,0)+J*H0
5       CONTINUE
        PART0(NSUBD+1)=SSL(A2(1,0),A2(2,0),INFO(IND))
        D(0)=SSL(A2(1,0),A2(2,0),INFO(IND))
        DO 10 J=0,NSUBD+1
            PSI0(J)=(PART0(J)-D(0))**2
10      CONTINUE

        DO 40 I=1,NP
            H1=MAX(FL0,(A2(2,I)-A2(1,I))/NSUBD)
            DO 15 J=0,NSUBD
                PART1(J)=A2(1,I)+J*H1
15          CONTINUE
            PART1(NSUBD+1)=SSL(A2(1,I),A2(2,I),INFO(IND+I))
            PSI1MN=FLE30
            DO 20 J=0,NSUBD+1
                D1=PART1(J)
                CALL SPRJ2 (A2(1,I-1),A2(2,I-1),D1,D1,I,INFO,
     *                      COMM,PART,NP,P1,P2)
                D0=SSL(P1,P2,D(I-1))
                IF (H0.GT.FL0) THEN
                    JD0=SMNIND(D0,PART0)
                ELSE
                    JD0=0
                END IF
                PSI1(J)=(D1-INFO(IND+I))**2+PSI0(JD0)
                IF (PSI1(J).LT.PSI1MN) THEN
                    PSI1MN=PSI1(J)
                    D(I)=D1
                END IF
20            CONTINUE
            H0=H1
            DO 30 J=0,NSUBD+1
                PSI0(J)=PSI1(J)
                PART0(J)=PART1(J)
30          CONTINUE
40      CONTINUE

        DO 50 I=NP,1,-1
            CALL SPRJ2 (A2(1,I-1),A2(2,I-1),D(I),D(I),I,INFO,
     *                  COMM,PART,NP,P1,P2)
            D(I-1)=SSL(P1,P2,D(I-1))
50      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SBL(X,N,L,X0,XN,TB,FLAG,LAUX0)

C  SBL computes the value assumed by the n-degree Bernstein polynomial
C  of a function  l  in the interval  (x0,xn)  at the point  x .
C  The evaluation is made using a Horner scheme, and the instructions
C  which do not depend upon  x  are executed under the control of
C  the variable  FLAG , for avoiding useless computations in
C  subsequent calls.
C  The degree  n  is supposed greater or equal to  3 .
C
C
C  INPUT PARAMETERS
C
C  X     : floating variable, containing the evaluation point.
C  N     : integer variable, containing the degree of Bernstein
C          polynomial.
C  L     : floating array, of bounds  0:N , containing the values
C          of the function  l .
C  X0    : floating variable, containing the left extreme of the
C          interval.
C  XN    : floating variable, containing the right extreme of the
C          interval.
C  TB    : floating array, of bounds  0:N , containing the binomial
C          terms used for computing the Bernstein polynomial.
C  FLAG  : integer variable, containing a control parameter.
C          In the case  flag=0  SBL  assumes to perform the first
C          evaluation of the polynomial, and computes the values
C          tb(i)*l(i) , i=0,1,...,n . In the case  flag=1  SBL
C          assumes to perform subsequent evaluations, and uses the
C          values previously computed.
C
C
C  OTHER PARAMETERS
C
C  LAUX0 : floating array, of bounds 0:N used as a work area to store
C          intermediate results.


        INTEGER N,FLAG,I

        REAL X,L(0:N),X0,XN,TB(0:N),LAUX0(0:N),XNMX,XMX0,
     *       AUX,FL1

        SAVE FL1
        DATA FL1/1.E00/

        IF(FLAG.EQ.0) THEN

            DO 10 I=0,N
                LAUX0(I)=TB(I)*L(I)
10          CONTINUE

        END IF

        XNMX=XN-X
        XMX0=X-X0
        AUX=FL1
        SBL=LAUX0(N)

        DO 20 I=N-1,0,-1
            AUX=XNMX*AUX
            SBL=LAUX0(I)*AUX+XMX0*SBL
20      CONTINUE

        SBL=SBL/(XN-X0)**N

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SBL1(X,N,L,X0,XN,TB,FLAG,LAUX1)

C  SBL1 computes the value assumed by the first derivative of an
C  n-degree Bernstein polynomial of a function  l  in the interval
C  (x0,xn)  at the point  x .
C  The evaluation is made using a Horner scheme, and the instructions
C  which do not depend upon  x  are executed under the control of
C  the variable  FLAG , for avoiding useless computations in
C  subsequent calls.
C  The degree  n  is supposed greater or equal to  3 .
C
C  INPUT PARAMETERS
C
C  X     : floating variable, containing the evaluation point.
C  N     : integer variable, containing the degree of Bernstein
C          polynomial.
C  L     : floating array, of bounds  0:N , containing the values
C          of the function  l .
C  X0    : floating variable, containing the left extreme of the
C          interval.
C  XN    : floating variable, containing the right extreme of the
C          interval.
C  TB    : floating array, of bounds  0:N-1 , containing the binomial
C          terms used for computing the Bernstein polynomial.
C  FLAG  : integer variable, containing a control parameter.
C          In the case  flag=0  SBL1  assumes to perform the first
C          evaluation of the polynomial, and computes the values
C          tb(i)*(l(i+1)-l(i)) , i=0,1,...,n-1 . In the case  flag=1
C          SBL1 assumes to perform subsequent evaluations, and uses 
C          the values previously computed.
C
C
C  OTHER PARAMETERS
C
C  LAUX1 : floating array, of bounds 0:N-1 used as a work area to store
C          intermediate results.


        INTEGER N,FLAG,I

        REAL X,L(0:N),X0,XN,TB(0:N-1),LAUX1(0:N-1),XNMX,
     *       XMX0,AUX,FL1

        SAVE FL1
        DATA FL1/1.E00/

        IF(FLAG.EQ.0) THEN

            DO 10 I=0,N-1
                LAUX1(I)=TB(I)*(L(I+1)-L(I))
10          CONTINUE

        END IF

        XNMX=XN-X
        XMX0=X-X0
        AUX=FL1
        SBL1=LAUX1(N-1)

        DO 20 I=N-2,0,-1
            AUX=XNMX*AUX
            SBL1=LAUX1(I)*AUX+XMX0*SBL1
20      CONTINUE

        SBL1=N*SBL1/(XN-X0)**N

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SBL2(X,N,L,X0,XN,TB,FLAG,LAUX2)

C  SBL2 computes the value assumed by the second derivative of an
C  n-degree Bernstein polynomial of a function  l  in the interval
C  (x0,xn)  at the point  x .
C  The evaluation is made using a Horner scheme, and the instructions
C  which do not depend upon  x  are executed under the control of
C  the variable  FLAG , for avoiding useless computations in
C  subsequent calls.
C  The degree  n  is supposed greater or equal to  3 .
C
C  INPUT PARAMETERS
C
C  X     : floating variable, containing the evaluation point.
C  N     : integer variable, containing the degree of Bernstein
C          polynomial.
C  L     : floating array, of bounds  0:N , containing the values
C          of the function  l .
C  X0    : floating variable, containing the left extreme of the
C          interval.
C  XN    : floating variable, containing the right extreme of the
C          interval.
C  TB    : floating array, of bounds  0:N-2 , containing the binomial
C          terms used for computing the Bernstein polynomial.
C  FLAG  : integer variable, containing a control parameter.
C          In the case  flag=0  SBL2  assumes to perform the first
C          evaluation of the polynomial, and computes the values
C          tb(i)*(l(i+2)-2*l(i+1)+l(i)) , i=0,1,...,n-2 .
C          In the case  flag=1  SBL2 assumes to perform subsequent 
C          evaluations, and uses the values previously computed.
C
C
C  OTHER PARAMETERS
C
C  LAUX2 : floating array, of bounds 0:N-2 used as a work area to store
C          intermediate results.


        INTEGER N,FLAG,I

        REAL X,L(0:N),X0,XN,TB(0:N-2),LAUX2(0:N-2),XNMX,
     *       XMX0,AUX,FL1

        SAVE FL1
        DATA FL1/1.E00/

        IF(FLAG.EQ.0) THEN

            DO 10 I=0,N-2
                LAUX2(I)=TB(I)*(L(I+2)-L(I+1)-L(I+1)+L(I))
10          CONTINUE

        END IF

        XNMX=XN-X
        XMX0=X-X0
        AUX=FL1
        SBL2=LAUX2(N-2)

        DO 20 I=N-3,0,-1
            AUX=XNMX*AUX
            SBL2=LAUX2(I)*AUX+XMX0*SBL2
20      CONTINUE

        SBL2=N*(N-1)*SBL2/(XN-X0)**N

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBNTAB (X,Y,NP,XTAB,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,D,D2,
     *                     TB,L,LAUX0,LAUX1,LAUX2,Y0TAB,Y1TAB,Y2TAB)

C  SBNTAB evaluates the spline and/or its first derivative and/or its
C  second derivative at the points  xtab(j) , j=0,1,...,ntab  using
C  a binary search for finding the interval  [x(i),x(i+1)] in which
C  the tabulation point falls. The input (X,Y,NP,XTAB,NTAB,Y0OPT,
C  Y1OPT,Y2OPT,N,K,D,D2,TB) and the output (Y0TAB,Y1TAB,Y2TAB)
C  parameters have been explained in subroutine DBVSSE. For the others
C  see subroutines STRMB, SLSPIS.


        INTEGER NP,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,IND,J

        REAL X(0:NP),Y(0:NP),XTAB(0:NTAB),D(0:NP),D2(0:NP),
     *       TB(1:N*(N+1)/2+N),L(0:N),LAUX0(0:N),LAUX1(0:N),
     *       LAUX2(0:N),Y0TAB(0:NTAB),Y1TAB(0:NTAB),
     *       Y2TAB(0:NTAB),SBL,SBL1,SBL2

        DO 40 J=0,NTAB

C  Call subprogram  SBSEAR  to compute the index  ind  such that
C       x(ind).LE.xtab(j).LT.x(ind+1) .

            CALL SBSEAR(X,NP,XTAB(J),IND)

C  Call subprogram  SLSPIS  to compute the linear shape-preserving
C  interpolating spline  l  at
C      x(ind)+p*(x(ind+1)-x(ind))/n , p=0,1,...,n .

            CALL SLSPIS(X,Y,D,D2,NP,N,K,IND,L)

            IF(Y0OPT.EQ.1) THEN

C  Evaluate the spline at  xtab(j) .

                Y0TAB(J)=SBL(XTAB(J),N,L,X(IND),X(IND+1),
     *                       TB(N*(N+1)/2),0,LAUX0)
            END IF

            IF(Y1OPT.EQ.1) THEN

C  Evaluate the first derivative of the spline at  xtab(j) .

                Y1TAB(J)=SBL1(XTAB(J),N,L,X(IND),X(IND+1),
     *                        TB((N-1)*N/2),0,LAUX1)
            END IF

            IF(Y2OPT.EQ.1) THEN

C  Evaluate the second derivative of the spline at  xtab(j) .

                Y2TAB(J)=SBL2(XTAB(J),N,L,X(IND),X(IND+1),
     *                        TB((N-2)*(N-1)/2),0,LAUX2)
            END IF

40      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBSEAR(X,NP,XTAB,IND)

C  Given an ordered set of points  (x(i), i=0,1,...,np)  and the
C  point  xtab , SBSEAR finds the index  ind  such that
C
C              x(ind) .LE. xtab .LT. x(ind+1)
C
C  using a standard binary search. SBSEAR  sets  ind=0  or  ind=np-1
C  whenever  xtab.LT.x(0)  or  x(np).LE.xtab .
C
C
C  INPUT PARAMETERS
C
C  X     : floating array, of bounds  0:NP , containing the set of
C          ordered points.
C  XTAB  : floating variable, containing the point to be processed.
C  NP    : integer  variable, defining the number of points of the
C          ordered set.
C
C
C  OUTPUT PARAMETERS
C
C  IND   : integer variable, whose value selects the interval in
C          which the point  xtab  falls.


        INTEGER NP,IND,I1,I2,MED

        REAL X(0:NP),XTAB

        IF(XTAB.LE.X(0)) THEN
            IND=0
            RETURN
        END IF

        IF(XTAB.GE.X(NP)) THEN
            IND=NP-1
            RETURN
        END IF

        I1=0
        I2=NP

10      CONTINUE
        IF (.NOT.(I1.NE.I2-1)) GO TO 20

            MED=(I1+I2)/2
            IF (XTAB.LT.X(MED)) THEN
                I2=MED
            ELSE IF (XTAB.GT.X(MED)) THEN
                I1=MED
            ELSE
                IND=MED
                RETURN
            END IF

        GO TO 10

20      CONTINUE
        IND=I1

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBVC (X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,INFO,
     *                   COMM,PART,EPS,KMAX,MAXSTP,A1,A2,DAUX2,
     *                   DAUX3,BETA,BETAI,RHO,RHOI,D,D2,ERRC,DIAGN)

C  SBVC checks input parameters and computes the required spline.
C
C  The input parameters X,Y,NP,N,K,OPT,D0,DNP,D20,D2NP,CONSTR,COMM,PART,
C  EPS,KMAX,MAXSTP,BETA,BETAI,RHO,RHOI and the output parameters
C  D,D2,ERRC,DIAGN are described in subroutine DBVSSC.
C  The other parameters are described in the called subprograms.


        EXTERNAL BETA,BETAI,RHO,RHOI

        INTEGER NP,N,K,OPT,CONSTR(0:NP-1),COMM,PART,KMAX,MAXSTP,
     *          ERRC,DIAGN(0:NP-1),P,Q,R,I

        REAL X(0:NP),Y(0:NP),D0,DNP,D20,D2NP,
     *       INFO(1:COMM+PART*NP+NP+1),EPS,A1(1:2,0:NP),
     *       A2(1:2,0:NP),D(0:NP),D2(0:NP),BETA,BETAI,
     *       RHO,RHOI,DAUX2(1:NP-1),DAUX3(0:NP-1)

C  Check the input parameters NP and OPT.

        P=OPT/100
        Q=MOD(OPT,100)/10
        R=MOD(OPT,10)

        IF ( (NP.LT.2) .OR. (P.LT.1.OR.P.GT.4) .OR.
     *       (Q.LT.1.OR.Q.GT.3) .OR. (R.LT.0.OR.R.GT.4) .OR.
     *       (P.EQ.1.AND.R.EQ.0) )                              THEN
             ERRC=1
             RETURN
        END IF

C  Check the array CONSTR.

        IF (MOD(OPT,10).EQ.4) THEN
            DO 10 I=0,NP-1
                IF (CONSTR(I).LT.0 .OR. CONSTR(I).GT.3) THEN
                    ERRC=2
                    RETURN
                END IF
10          CONTINUE
        END IF

C  Check the input parameters N and K.

        IF ( (K.LT.1 .OR. K.GT.2) .OR. (N.LT.3*K) ) THEN
            ERRC=3
            RETURN
        END IF

C  Check the abscissas of the interpolation points.

        DO 20 I=0,NP-1
            IF(X(I+1).LE.X(I)) THEN
                ERRC=4
                RETURN
            END IF
20      CONTINUE

C  Call subprogram SSTINF to set the information array INFO.

C  Initialize the array DIAGN.

        DO 30 I=0,NP-1
            DIAGN(I)=0
30      CONTINUE

        CALL SSTINF (OPT,D0,DNP,CONSTR,N,K,X,Y,D,NP,COMM,PART,EPS,
     *                     BETA,BETAI,DAUX2,DAUX3,INFO,ERRC,DIAGN)

C  Call subprogram SALG3 to compute the array D containing the first
C  derivative at initial points.

        CALL SALG3 (INFO,NP,COMM,PART,OPT,D0,DNP,EPS,KMAX,MAXSTP,
     *              BETA,BETAI,A1,A2,D,ERRC,DIAGN)

        IF (K.EQ.2) THEN

C  A  C(2) spline is required. Compute the sequence of second derivati-
C  ves d2(i), i=0,...,np , according to the shape constraints and, if
C  possible, to boundary conditions.

            CALL SSCDRC (N,X,Y,D,OPT,NP,EPS,D20,D2NP,
     *                   RHO,RHOI,A1,A2,DAUX3,D2,ERRC)
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SBVE (X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,
     *                   Y2OPT,D,D2,ERRC,TB,L,LAUX0,LAUX1,LAUX2,
     *                   Y0TAB,Y1TAB,Y2TAB,ERRE)

C  SBVE checks input parameters and evaluates the required spline.
C
C  The input parameters X,Y,NP,N,K,XTAB,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,
C  D,D2,ERRC and the output parameters Y0TAB,Y1TAB,Y2TAB,ERRE are
C  described in subroutine DBVSSE. The others are used as work areas
C  and will be eventually described in the subsequent routines.

 
        INTEGER NP,N,K,NTAB,SBOPT,Y0OPT,Y1OPT,Y2OPT,ERRC,ERRE,I

        REAL X(0:NP),Y(0:NP),XTAB(0:NTAB),D(0:NP),D2(0:NP),
     *       L(0:N),LAUX0(0:N),LAUX1(0:N),LAUX2(0:N),
     *       TB(1:N*(N+1)/2+N),Y0TAB(0:NTAB),Y1TAB(0:NTAB),
     *       Y2TAB(0:NTAB)

C  Check the correctness of input data, that is if subroutine DBVSSC
C  has correctly run.

        IF (ERRC.EQ.1.OR.ERRC.EQ.2.OR.ERRC.EQ.3.OR.ERRC.EQ.4
     *      .OR.ERRC.EQ.5.OR.ERRC.EQ.7.OR.ERRC.EQ.8.OR.ERRC.EQ.10) THEN
            ERRE = 1
            RETURN
        END IF

C  Check the input parameters NTAB, SBOPT, Y0OPT, Y1OPT, Y2OPT.

        IF (NTAB.LT.0.OR.(SBOPT.NE.1.AND.SBOPT.NE.2).OR.
     *      (Y0OPT.NE.0.AND.Y0OPT.NE.1).OR.(Y1OPT.NE.0.AND.Y1OPT.NE.1)
     *      .OR.(Y2OPT.NE.0.AND.Y2OPT.NE.1)) THEN
            ERRE=2
            RETURN
        END IF

        IF(SBOPT.EQ.1.AND.NTAB.GT.0) THEN

C  Check the abscissas of the tabulation points when the sequential
c  search is required.

            DO 10 I=0,NTAB-1
                IF(XTAB(I+1).LE.XTAB(I)) THEN
                    ERRE=3
                    RETURN
                END IF
10          CONTINUE
        END IF

C  Call subprogram STRMB to compute the binomial terms needed
C  in the expression of Bernstein polynomials.

        CALL STRMB(N,TB)

        IF(SBOPT.EQ.1) THEN

C  sbopt=1:  sequential search is required.

            CALL SSQTAB (X,Y,NP,XTAB,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,D,D2,
     *                   TB,L,LAUX0,LAUX1,LAUX2,Y0TAB,Y1TAB,Y2TAB)
        ELSE

C  sbopt=2: binary search is required.

            CALL SBNTAB (X,Y,NP,XTAB,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,D,D2,
     *                   TB,L,LAUX0,LAUX1,LAUX2,Y0TAB,Y1TAB,Y2TAB)
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SFPSVF (A1,A2,NP,INFO,COMM,PART,EPS,MAXSTP,
     *                     BETA,ERRC,DSTAR)

C  SFPSVF finds, if possible, dstar.IN.[a1(1,0),a1(2,0)] such that
C               beta(dstar) .INT. phi(dstar)                     (1)
C
C  The input parameters NP,COMM,PART,EPS,MAXSTP,BETA, and the output
C  parameter ERRC are described in DBVSSC. The input parameters A1 and
C  A2 are described in  SALG1.


        EXTERNAL BETA

        INTEGER NP,COMM,PART,MAXSTP,ERRC,STEP,SUBSTP,ERRC1,I

        REAL A1(1:2,0:NP),A2(1:2,0:NP),
     *       INFO(1:COMM+PART*NP+NP+1),EPS,BETA,DSTAR,MID,
     *       H,SSL

        SAVE FL1D2
        REAL FL1D2
        DATA FL1D2/0.5E00/

C  If the optimum input value of dstar does not belong to 
C  [a1(1,0),a1(2,0)] , the nearest extreme of this interval
C  replaces the old value of dstar.

        DSTAR=SSL(A1(1,0),A1(2,0),DSTAR)

C  Compute phi(dstar).

        CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,A2,ERRC1)

C  If phi(dstar) is not empty and dstar satisfies (1), it is the desired
C  value.

        IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *      BETA(DSTAR).LE.A2(2,NP)+EPS))                      RETURN

C  If it is not the case, look for another value. First, check
C  if the midpoint of the interval of all possible values satisfies
C  condition (1).

        MID=(A1(1,0)+A1(2,0))*FL1D2
        DSTAR=MID
        CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,A2,ERRC1)
        IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *      BETA(DSTAR).LE.A2(2,NP)+EPS))                       RETURN

C  Second, check if any point of a tabulation of interval
C  [a1(1,0),a1(2,0)]  satisfies the condition (1). The tabulation 
C  points are given by jumps of decreasing lenghts with alternate 
C  direction with respect to the middle of the interval.

        H=(A1(2,0)-A1(1,0))*FL1D2
        DO 20 STEP=1,MAXSTP
            H=H*FL1D2
            SUBSTP=2**STEP-1
            DO 10 I=1,SUBSTP,2
                DSTAR=MID-I*H
                CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,
     *                       A2,ERRC1)
                IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *              BETA(DSTAR).LE.A2(2,NP)+EPS))               RETURN
                DSTAR=MID+I*H
                CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,
     *                       A2,ERRC1)
                IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *              BETA(DSTAR).LE.A2(2,NP)+EPS))               RETURN
10          CONTINUE
20      CONTINUE

C  Finally, check if condition (1) is satisfied by one of the
C  [a1(1,0),a1(2,0)] extremes.

        DSTAR=A1(1,0)
        CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,A2,ERRC1)
        IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *      BETA(DSTAR).LE.A2(2,NP)+EPS))                       RETURN
        DSTAR=A1(2,0)
        CALL SALG1D (DSTAR,A1,NP,INFO,COMM,PART,EPS,A2,ERRC1)
        IF (ERRC1.EQ.0.AND.(A2(1,NP)-EPS.LE.BETA(DSTAR).AND.
     *      BETA(DSTAR).LE.A2(2,NP)+EPS))                       RETURN

C  If dstar satisfying (1) has not been found, send a message resetting
C  the error flag errc.

        IF (ERRC.EQ.6) THEN
            ERRC=8
        ELSE
            ERRC=7
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SINTRS (A,B,C,D,P1,P2)

C  SINTRS computes the intersection of the two intervals  [a,b]
C  and [c,d]. [p1,p2] is the result. The output  p1.GT.p2 means that
C  the intervals are disjoint. SINTRS assumes  a.LE.b  and  c.LE.d .


        REAL A,B,C,D,P1,P2

        P1=MAX(A,C)
        P2=MIN(B,D)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SLSPIS(X,Y,D,D2,NP,N,K,IND,L)

C  SLSPIS   evaluates the control points of the Bernstein-Bezier net
C  l:=l(x) of the interpolating spline  s:=s(x) , s.IN.BS(n;k) in the
C  interval  [x(ind),x(ind+1)] . For a description of the function  l
C  see the comments in subroutines  DBVSSC and SSCDRC. Here we only
C  recall that the structure of the net is different for k=1 or k=2 .
C
C  The input parameters  X,Y,D,D2,NP,N,K  are explained in subroutine
C  SPISE.
C
C  OTHER PARAMETERS
C
C  IND   : integer variable, used to select the knots interval.
C  L     : floating array, of bounds  0:N , containing the ordinates
C          of the control points.


        INTEGER NP,N,K,IND,I

        REAL X(0:NP),Y(0:NP),D(0:NP),D2(0:NP),L(0:N),H,
     *       ALPHA,Q1,Q2,FL1,FL2,FL4

        SAVE FL1,FL2,FL4
        DATA FL1,FL2,FL4/1.0E00,2.0E00,4.0E00/

        H=(X(IND+1)-X(IND))

C  Compute the net in the case  k=1 .

        IF (K.EQ.1) THEN

            L(0)=Y(IND)
            L(1)=Y(IND)+H*D(IND)/N
            L(N)=Y(IND+1)
            L(N-1)=Y(IND+1)-H*D(IND+1)/N

            DO 10 I=2,N-2
                ALPHA=(I-FL1)/(N-FL2)
                L(I)=(FL1-ALPHA)*L(1)+ALPHA*L(N-1)
10          CONTINUE

        ELSE IF (K.EQ.2) THEN

C  Compute the net in the case  k=2 .

            L(0)=Y(IND)
            L(1)=L(0)+H*D(IND)/N
            L(2)=H**2*D2(IND)/(N*(N-1))+FL2*L(1)-L(0)
            L(N)=Y(IND+1)
            L(N-1)=L(N)-H*D(IND+1)/N
            L(N-2)=H**2*D2(IND+1)/(N*(N-1))+FL2*L(N-1)-L(N)

            ALPHA=((N-4)/2)/(N-FL4)
            Q1=ALPHA*(Y(IND+1)-FL2*H*D(IND+1)/N)+
     *         (FL1-ALPHA)*(Y(IND)+FL2*H*D(IND)/N)
            Q2=(FL1-ALPHA)*(Y(IND+1)-FL2*H*D(IND+1)/N)+
     *          ALPHA*(Y(IND)+FL2*H*D(IND)/N)

            DO 30 I=3,N/2
                ALPHA=(I-FL2)/((N/2)-FL2)
                L(I)=(FL1-ALPHA)*L(2)+ALPHA*Q1
30          CONTINUE

            DO 40 I=N/2+1,N-3
                ALPHA=(N-FL2-I)/((N/2)-FL2)
                L(I)=(FL1-ALPHA)*L(N-2)+ALPHA*Q2
40          CONTINUE

        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SMDIAN(A,B,C)

C  Given three numbers a,b,c , median  is the one which lies between
C  the other two.


        REAL A,B,C

        SMDIAN=MIN(MAX(A,B),MAX(B,C),MAX(C,A))

        RETURN
        END

C  ---------------------------------------------------------------------

        INTEGER FUNCTION SMNIND(D,PART)

C  SMNIND finds the index of the component of the array PART closest
C  to d .


        INTEGER NSUBD
        PARAMETER (NSUBD=20)

        INTEGER J

        REAL D,PART(0:NSUBD+1),AUX(0:NSUBD+1),MINDIS

        SAVE FL30
        REAL FL30
        DATA FL30/1.0E30/

        DO 10 J=0,NSUBD+1
            AUX(J)=ABS(D-PART(J))
10      CONTINUE

        MINDIS=FL30

        DO 20 J=0,NSUBD+1
            IF (AUX(J).LT.MINDIS) THEN
                SMNIND=J
                MINDIS=AUX(J)
            END IF
20      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SMNMOD(A,B)

C  Given two real numbers a and b, SMNMOD returns the number between
C  a and b which is closest to zero.


        REAL A,B,FL1,FL2

        SAVE FL1,FL2
        DATA FL1,FL2/1.0E00,2.0E00/

        SMNMOD=((SIGN(FL1,A)+SIGN(FL1,B))/FL2)*MIN(ABS(A),ABS(B))

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SMSK1(INFO,CONSTR,COMM,PART,IND1,NP)

C  SMSK1 compares the constraints required in input by the user and
C  stored in the array CONSTR with the shape of the data, stored by
C  SSTINF in the array INFO. If the required and real shapes do not
C  agree, SMSK1 resets both INFO and CONSTR with the 'intersection'
C  of the shapes. For example, if info(ind1+i)=11 , that is the data
C  are increasing and convex, and constr(i)=2 , that is only convexity
C  is required, then the output values will be  info(ind1+i)=10
C  (convexity) and constr(i)=2 (unchanged). If  info(ind1+i)=20
C  (concavity) and  constr(i)=1 (monotonicity) the output will be
C  info(ind1+i)=constr(i)=0 (no constraints). So, the computations made
C  in SALG3 will be based on these new values for selecting the domains
C  of admissible derivatives, and CONSTR will contain information on
C  the constraints effectively imposed.
C  Further details on the parameters INFO and IND1 can be found in sub-
C  routine SSTINF; CONSTR, COMM, PART, NP are explaained in subroutine
C  DBVSSC.


        INTEGER COMM,PART,IND1,NP,I
        INTEGER CONSTR(0:NP-1)

        REAL INFO(1:COMM+PART*NP+NP+1)

        DO 10 I=0,NP-1

            IF ( INFO(IND1+I).EQ.0 .OR. CONSTR(I).EQ.0 ) THEN
                INFO(IND1+I)=0
                CONSTR(I)=0
            ELSE IF ( INFO(IND1+I).EQ.1 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=1
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=0
                    CONSTR(I)=0
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=1
                    CONSTR(I)=1
                END IF
            ELSE IF ( INFO(IND1+I).EQ.2 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=2
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=0
                    CONSTR(I)=0
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=2
                    CONSTR(I)=1
                END IF
            ELSE IF ( INFO(IND1+I).EQ.10 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=0
                    CONSTR(I)=0
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=10
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=10
                    CONSTR(I)=2
                END IF
            ELSE IF ( INFO(IND1+I).EQ.20 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=0
                    CONSTR(I)=0
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=20
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=20
                    CONSTR(I)=2
                END IF
            ELSE IF ( INFO(IND1+I).EQ.11 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=1
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=10
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=11
                    CONSTR(I)=3
                END IF
            ELSE IF ( INFO(IND1+I).EQ.21 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=1
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=20
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=21
                    CONSTR(I)=3
                END IF
            ELSE IF ( INFO(IND1+I).EQ.12 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=2
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=10
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=12
                    CONSTR(I)=3
                END IF
            ELSE IF ( INFO(IND1+I).EQ.22 ) THEN
                IF ( CONSTR(I).EQ.1 ) THEN
                    INFO(IND1+I)=2
                    CONSTR(I)=1
                ELSE IF ( CONSTR(I).EQ.2 ) THEN
                    INFO(IND1+I)=20
                    CONSTR(I)=2
                ELSE IF ( CONSTR(I).EQ.3 ) THEN
                    INFO(IND1+I)=22
                    CONSTR(I)=3
                END IF
            END IF

10      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SMSK2(INFO,COMM,PART,IND1,NP,D0,DNP,EPS,
     *                   ERRC,DIAGN)

C  This routine controls if the separable boundary conditions d(0)=d0
C  and d(np)=dnp are compatible with the first and the last domain of
C  constraints. The error flag is reset correspondingly.
C  Details on the parameters INFO, IND1 and COMM, PART, NP, D0, DNP,
C  EPS, ERRC, DIAGN can be found in subroutines 
C  SSTINF and DBVSSC respectively.


        INTEGER COMM,PART,IND1,NP,ERRC,DIAGN(0:NP-1)

        REAL INFO(1:COMM+PART*NP+NP+1),D0,DNP,EPS,P1,P2,A,B

        CALL SPRJ0 (1,INFO,COMM,PART,NP,P1,P2)

        IF ( .NOT.(P1-EPS.LE.D0 .AND. D0.LE.P2+EPS) ) THEN
            INFO(IND1)=0
            ERRC=9
            DIAGN(0)=1
        END IF

        CALL SPRJ0(NP,INFO,COMM,PART,NP,A,B)
        CALL SPRJ1(A,B,INFO(4),INFO(5),NP,INFO,COMM,PART,NP,P1,P2)

        IF ( .NOT.(P1-EPS.LE.DNP .AND. DNP.LE.P2+EPS) ) THEN
            INFO(IND1+NP-1)=0
            ERRC=9
            DIAGN(NP-1)=1
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SPRJ0 (I,INFO,COMM,PART,NP,P1,P2)

C  Given the integer i , SPRJ0 computes the set B(i) performing the
C  projection of D(i) (subset of the (i-1)i-plane) onto the (i-1)-axis.
C
C  The input parameters COMM,PART,NP are described in DBVSSC; the input
C  parameter INFO is described in subroutine SSTINF.
C
C  OUTPUT PARAMETERS:
C
C  P1  : floating variable, containing the left extreme of the
C        resulting interval.
C
C  P2  : floating variable, containing the right extreme of the
C        resulting interval.


        INTEGER I,COMM,PART,NP,KIND

        REAL INFO(1:COMM+PART*NP+NP+1),P1,P2,N,K,DEL

        SAVE FL0
        REAL FL0
        DATA FL0/0.0E00/

        N=INFO(1)
        K=INFO(2)
        KIND=INFO(COMM+1+(I-1))
        DEL=INFO(COMM+NP+1+(I-1))

C  No constraint

        IF (KIND.EQ.0) THEN
            P1=INFO(4)
            P2=INFO(5)

C  Increase constraints

        ELSE IF (KIND.EQ.1) THEN
            P1=FL0
            P2=N*DEL/K

C  Decrease constraints

        ELSE IF (KIND.EQ.2) THEN
            P1=N*DEL/K
            P2=FL0

C  Convexity constraints

        ELSE IF (KIND.EQ.10) THEN
            P1=INFO(4)
            P2=DEL

C  Concavity constraints

        ELSE IF (KIND.EQ.20) THEN
            P1=DEL
            P2=INFO(5)

C  Increase and convexity

        ELSE IF (KIND.EQ.11) THEN
            P1=FL0
            P2=DEL

C  Increase and concavity

        ELSE IF (KIND.EQ.21) THEN
            P1=DEL
            P2=N*DEL/K

C  Decrease and convexity

        ELSE IF (KIND.EQ.12) THEN
            P1=N*DEL/K
            P2=DEL

C  Decrease and concavity

        ELSE IF (KIND.EQ.22) THEN
            P1=DEL
            P2=FL0

        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SPRJ1 (A,B,C,D,I,INFO,COMM,PART,NP,P1,P2)

C  Given the set S=[a,b]x[c,d] and the integer i , SPRJ1 performs the
C  intersection of S with the domain D(i) and the projection of the
C  resulting set (a subset of (i-1)i-plane) onto the i-axis .
C
C  The input parameters COMM,PART,NP are described in DBVSSC; the input
C  parameter INFO is described in SSTINF.
C
C  OUTPUT PARAMETERS:
C
C  P1  : floating variable, containing the left extreme of the
C        resulting interval.
C
C  P2  : floating variable, containing the right extreme of the
C        resulting interval.


        INTEGER I,COMM,PART,NP,KIND

        REAL A,B,C,D,INFO(1:COMM+PART*NP+NP+1),P1,P2,N,K,
     *       DEL,F1,F2,F3,X

        SAVE FL0
        REAL FL0
        DATA FL0/0.0E00/

        F1(X,N,K,DEL) = -X+N*DEL/K
        F2(X,N,K,DEL) = -K*X/(N-K)+N*DEL/(N-K)
        F3(X,N,K,DEL) = -(N-K)*X/K+N*DEL/K

        N=INFO(1)
        K=INFO(2)
        KIND=INFO(COMM+1+(I-1))
        DEL=INFO(COMM+NP+1+(I-1))

C  No constraint

        IF (KIND.EQ.0) THEN
            P1=C
            P2=D

C  Increase constraints

        ELSE IF (KIND.EQ.1) THEN
            P1=FL0
            P2=F1(A,N,K,DEL)

C  Decrease constraints

        ELSE IF (KIND.EQ.2) THEN
            P1=F1(B,N,K,DEL)
            P2=FL0

C  Convexity constraints

        ELSE IF (KIND.EQ.10) THEN
            P1=F2(B,N,K,DEL)
            P2=F3(A,N,K,DEL)

C  Concavity constraints

        ELSE IF (KIND.EQ.20) THEN
            P1=F3(B,N,K,DEL)
            P2=F2(A,N,K,DEL)

C  Increase and convexity

        ELSE IF (KIND.EQ.11) THEN
            P1=F2(B,N,K,DEL)
            P2=F3(A,N,K,DEL)

C  Increase and concavity

        ELSE IF (KIND.EQ.21) THEN
            P1=MAX(FL0,F3(B,N,K,DEL))
            P2=F2(A,N,K,DEL)

C  Decrease and convexity

        ELSE IF (KIND.EQ.12) THEN
            P1=F2(B,N,K,DEL)
            P2=MIN(FL0,F3(A,N,K,DEL))

C  Decrease and concavity

        ELSE IF (KIND.EQ.22) THEN
            P1=F3(B,N,K,DEL)
            P2=F2(A,N,K,DEL)

        END IF

        P1=MAX(P1,C)
        P2=MIN(P2,D)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SPRJ2 (A,B,C,D,I,INFO,COMM,PART,NP,P1,P2)

C  Given the set s=[a,b]x[c,d] and the integer i, SPRJ2 performs the
C  intersection of S with the domain D(i) and the projection of the
C  resulting set (subset of (i-1)i-plane) onto the (i-1)-axis .
C
C  The input parameters COMM,PART,NP are described in DBVSSC; the input
C  parameter INFO is described in SSTINF.
C
C  OUTPUT PARAMETERS:
C
C  P1  : floating variable, containing the left extreme of the
C        resulting interval.
C
C  P2  : floating variable, containing the right extreme of the
C        resulting interval.


        INTEGER I,COMM,PART,NP,KIND

        REAL A,B,C,D,INFO(1:COMM+PART*NP+NP+1),P1,P2,N,K,
     *       DEL,F1I,F2I,F3I,X

        SAVE FL0
        REAL FL0
        DATA FL0/0.0E00/

        F1I(X,N,K,DEL) = -X+N*DEL/K
        F2I(X,N,K,DEL) = -(N-K)*X/K+N*DEL/K
        F3I(X,N,K,DEL) = -K*X/(N-K)+N*DEL/(N-K)

        N=INFO(1)
        K=INFO(2)
        KIND=INFO(COMM+1+(I-1))
        DEL=INFO(COMM+NP+1+(I-1))

C  No constraints

        IF (KIND.EQ.0) THEN
            P1=A
            P2=B

C  Increase constraints

        ELSE IF (KIND.EQ.1) THEN
            P1=FL0
            P2=F1I(C,N,K,DEL)

C  Decrease constraints

        ELSE IF (KIND.EQ.2) THEN
            P1=F1I(D,N,K,DEL)
            P2=FL0

C  Convexity constraints

        ELSE IF (KIND.EQ.10) THEN
            P1=F2I(D,N,K,DEL)
            P2=F3I(C,N,K,DEL)

C  Concavity constraints

        ELSE IF (KIND.EQ.20) THEN
            P1=F3I(D,N,K,DEL)
            P2=F2I(C,N,K,DEL)

C  Increase and convexity

        ELSE IF (KIND.EQ.11) THEN
            P1=MAX(FL0,F2I(D,N,K,DEL))
            P2=F3I(C,N,K,DEL)

C  Increase and concavity

        ELSE IF (KIND.EQ.21) THEN
            P1=F3I(D,N,K,DEL)
            P2=F2I(C,N,K,DEL)

C  Decrease and convexity

        ELSE IF (KIND.EQ.12) THEN
            P1=F2I(D,N,K,DEL)
            P2=F3I(C,N,K,DEL)

C  Decrease and concavity

        ELSE IF (KIND.EQ.22) THEN
            P1=F3I(D,N,K,DEL)
            P2=MIN(FL0,F2I(C,N,K,DEL))

        END IF

        P1=MAX(P1,A)
        P2=MIN(P2,B)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SSCDRC (N,X,Y,D,OPT,NP,EPS,D20,D2NP,RHO,RHOI,
     *                     A1,A2,H,D2,ERRC)

C  SSCDRC computes the sequence  d2(i) , i=0,1,...,np , of second
C  derivatives at the knots. The basic idea is that the vanishing second
C  derivatives (which are admissible by virtue of the theory involved in
C  the routines called previously) can be locally changed to non-zero
C  values without modifying the monotonicity and/or convexity.
C  Let us consider the restriction to the i-th subinterval of the
C  Bernstein-Bezier net for the C(2) spline with zero derivatives given
C  by subroutine SALG3. Let A, B and C be the second, third and 
C  (int(n/2))-th point of the net, and let E, F, and G be given by a
C  symmetric construction.
C
C            B_______________C___G_______________F             
C           /           .             .           \                   
C          /      .                         .      \                  
C         /  .D                                 H.  \                 
C        A                                           E                
C       /                                             \               
C      /                                               \              
C     /                                                 \             
C
C  Then the 'intermediate net' obtained inserting the straight lines
C  trough A-C and E-F is shape-preserving and we can take as the 'final
C  net' the union of two convex combination of A-B-C , A-D-C and H-F-G ,
C  E-H-D respectively. Expressing the net in term of the second
C  derivatives, the points D, B and H,F lead to restriction like
C  d2(i).IN.[a1(1,i),a1(2,i)] , d2(i+1).IN.[a2(1,i),a2(2,i)]
C  This construction must be repeated for all the subintervals and so
C  d2(i) .IN. [a2(1,i-1),a2(2,i-1)].INT.[a1(1,i),a1(2,i)] .
C
C  The input parameters N,X,Y,D,OPT,NP,EPS,D20,D2NP,RHO,RHOI and the
C  input ones D2,ERRC are documented in subroutine DBVSSC.


        EXTERNAL RHO,RHOI

        INTEGER NP,N,OPT,ERRC,I,Q

        REAL X(0:NP),Y(0:NP),D20,D2NP,EPS,D(0:NP),D2(0:NP),
     *       H(0:NP-1),A1(1:2,0:NP),A2(1:2,0:NP),A,B,C,DD,E,
     *       F,G,HH,ALPHA,GAMMA,DIFF2,P1,P2,Q1,Q2,FL0,FL1,
     *       FL2,FL4,RHO,RHOI,SSL

        SAVE FL0,FL1,FL2,FL4
        DATA FL0,FL1,FL2,FL4/0.0E00,1.0E00,2.0E00,4.0E00/

        DO 10 I=0,NP-1
            H(I)=X(I+1)-X(I)
10      CONTINUE

        DO 20 I=0,NP-1

C  Compute the points of the 'original' and 'intermediate' net.

            A=Y(I)+H(I)*D(I)/N
            B=Y(I)+FL2*H(I)*D(I)/N
            E=Y(I+1)-H(I)*D(I+1)/N
            F=Y(I+1)-FL2*H(I)*D(I+1)/N

            ALPHA=((N-4)/2)/(N-FL4)
            C=ALPHA*F+(FL1-ALPHA)*B
            G=(FL1-ALPHA)*F+ALPHA*B

            GAMMA=FL1/(((N-4)/2)+FL1)
            DD=GAMMA*C+(FL1-GAMMA)*A
            HH=GAMMA*G+(FL1-GAMMA)*E

C  Define the left and the right restriction for the second finite
C  difference of the net.

            A1(1,I)=MIN(FL0,DD-FL2*A+Y(I))
            A1(2,I)=MAX(FL0,DD-FL2*A+Y(I))

            A2(1,I+1)=MIN(FL0,Y(I+1)-FL2*E+HH)
            A2(2,I+1)=MAX(FL0,Y(I+1)-FL2*E+HH)

20      CONTINUE

C  Take the intersection of the left and right restrictions for the
C  same second differences and translate it in terms of the second
C  derivatives.

        A1(1,0)=A1(1,0)*N*(N-1)/H(0)**2
        A1(2,0)=A1(2,0)*N*(N-1)/H(0)**2

        DO 30 I=1,NP-1

            CALL SINTRS (A1(1,I),A1(2,I),A2(1,I),A2(2,I),P1,P2)
            A1(1,I)=P1*N*(N-1)/H(I)**2
            A1(2,I)=P2*N*(N-1)/H(I)**2

30      CONTINUE
        
        A1(1,NP)=A2(1,NP)*N*(N-1)/H(NP-1)**2
        A1(2,NP)=A2(2,NP)*N*(N-1)/H(NP-1)**2

C  The internal derivatives are defined as the admissible value closest
C  to the central second divided difference of the data.

        DO 40 I=1,NP-1

            DIFF2=( (Y(I+1)-Y(I))/H(I) - (Y(I)-Y(I-1))/H(I-1) ) /
     *                    ( H(I)+H(I-1) )
            D2(I)=SSL(A1(1,I),A1(2,I),DIFF2)

40      CONTINUE

        Q=MOD(OPT,100)/10

        IF (Q.EQ.1) THEN

C  No boundary condition is required. Take the first and last
C  derivative as the middle of admissible values.

            D2(0)=(A1(1,0)+A1(2,0))/FL2
            D2(NP)=(A1(1,NP)+A1(2,NP))/FL2

        ELSE IF (Q.EQ.2) THEN

C  Non-separable boundary conditions are required. Check if these can be
C  satisfied by admissible derivatives.

            Q1=MIN(RHO(A1(1,0)),RHO(A1(2,0)))
                 Q2=MAX(RHO(A1(1,0)),RHO(A1(2,0)))
            CALL SINTRS(A1(1,NP),A1(2,NP),Q1,Q2,P1,P2)

            IF(P1.GT.P2+EPS) THEN

C  The boundary conditions cannot be satisfied. Set the error flag and
C  define the first and the last derivative as the nearest point to the
C  admissible and the boundary interval.

                IF (ERRC.EQ.0) THEN
                    ERRC=15
                ELSE IF (ERRC.EQ.6) THEN
                    ERRC=16
                ELSE IF (ERRC.EQ.9) THEN
                    ERRC=17
                ELSE IF (ERRC.EQ.10) THEN
                    ERRC=18
                END IF

                D2(NP)=SSL(A1(1,NP),A1(2,NP),(P1+P2)/FL2)
                D2(0)=SSL(A1(1,0),A1(2,0),RHOI(D2(NP)))

            ELSE

C  It is possible to satisfy the boundary conditions.

                D2(NP)=(P1+P2)/FL2
                D2(0)=RHOI(D2(NP))

            END IF

        ELSE IF (Q.EQ.3) THEN

C  Separable boundary conditions are required. Check if they are
C  compatible with the admissible intervals and, if not, set the
C  error flag and take the admissible points nearest to the boundary
C  values. Otherwise take simply the boundary values.

            IF ( D20.LT.A1(1,0)-EPS .OR. D20.GT.A1(2,0)+EPS .OR.
     *           D2NP.LT.A1(1,NP)-EPS .OR. D2NP.GT.A1(2,NP)+EPS ) THEN
                IF (ERRC.EQ.0) THEN
                    ERRC=11
                ELSE IF (ERRC.EQ.6) THEN
                    ERRC=12
                ELSE IF (ERRC.EQ.9) THEN
                    ERRC=13
                ELSE IF (ERRC.EQ.10) THEN
                    ERRC=14
                END IF
            END IF

            D2(0)=SSL(A1(1,0),A1(2,0),D20)
            D2(NP)=SSL(A1(1,NP),A1(2,NP),D2NP)

        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        REAL FUNCTION SSL(A,B,C)

C  Given the interval [a,b] and the number c, ssl is c if c belongs
C  to [a,b], otherwise, it is the nearest extreme to c.


        REAL A,B,C

        IF (C.LE.A) THEN
            SSL=A
        ELSE IF (C.GE.B) THEN
            SSL=B
        ELSE
            SSL=C
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SSQTAB (X,Y,NP,XTAB,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,D,D2,
     *                     TB,L,LAUX0,LAUX1,LAUX2,Y0TAB,Y1TAB,Y2TAB)

C  SSQTAB evaluates the spline and/or its first derivative and/or its
C  second derivative at the points  xtab(j) , j=0,1,...,ntab  using
C  a sequential search for finding the interval  [x(i),x(i+1)] in which
C  the tabulation point falls. The input (X,Y,NP,XTAB,NTAB,Y0OPT,
C  Y1OPT,Y2OPT,N,K,D,D2,TB) and the output (Y0TAB,Y1TAB,Y2TAB)
C  parameters have been explained in subroutine DBVSSE. For the others
C  see subroutines STRMB, SLSPIS.


        INTEGER NP,NTAB,Y0OPT,Y1OPT,Y2OPT,N,K,IND,IND1,J,I

        REAL X(0:NP),Y(0:NP),XTAB(0:NTAB),D(0:NP),D2(0:NP),
     *       TB(1:N*(N+1)/2+N),L(0:N),LAUX0(0:N),LAUX1(0:N),
     *       LAUX2(0:N),Y0TAB(0:NTAB),Y1TAB(0:NTAB),
     *       Y2TAB(0:NTAB),SBL,SBL1,SBL2

        IND=0
        IND1=1

        DO 30 J=0,NTAB

C  Compute the index  ind  such that  x(ind).LE.xtab(j).LT.x(ind+1) .

            IF(X(0).LE.XTAB(J)) THEN

               DO 20 I=IND1,NP-1
                    IF(X(I).LE.XTAB(J)) IND=I
20              CONTINUE

            END IF

C  Check if  ind  selects a new subinterval.

            IF(IND.NE.IND1.OR.J.EQ.0) THEN

C  Call subprogram  SLSPIS  to compute the linear shape-preserving
C  interpolating spline  l:=l(x)  at
C      x(ind)+p*(x(ind+1)-x(ind))/n , p=0,1,...,n .

                CALL SLSPIS(X,Y,D,D2,NP,N,K,IND,L)

                IF(Y0OPT.EQ.1) THEN

C  Evaluate the spline at  xtab(j)  using new values of  l .

                    Y0TAB(J)=SBL(XTAB(J),N,L,X(IND),X(IND+1),
     *                           TB(N*(N+1)/2),0,LAUX0)
                END IF

                IF(Y1OPT.EQ.1) THEN

C  Evaluate the first derivative of the spline at  xtab(j)  using new
C  values of  l .
C
                    Y1TAB(J)=SBL1(XTAB(J),N,L,X(IND),X(IND+1),
     *                            TB((N-1)*N/2),0,LAUX1)
                END IF

                IF(Y2OPT.EQ.1) THEN

C  Evaluate the second derivative of the spline at  xtab(j)  using new
C  values of  l .

                    Y2TAB(J)=SBL2(XTAB(J),N,L,X(IND),X(IND+1),
     *                            TB((N-2)*(N-1)/2),0,LAUX2)
                END IF

            ELSE

                IF(Y0OPT.EQ.1) THEN

C  Evaluate the spline at  xtab(j)  using old values of  l .

                    Y0TAB(J)=SBL(XTAB(J),N,L,X(IND),X(IND+1),
     *                           TB(N*(N+1)/2),1,LAUX0)
                END IF

                IF(Y1OPT.EQ.1) THEN

C  Evaluate the first derivative of the spline at  xtab(j)  using old
C  values of  l .

                    Y1TAB(J)=SBL1(XTAB(J),N,L,X(IND),X(IND+1),
     *                            TB((N-1)*N/2),1,LAUX1)
                END IF

                IF(Y2OPT.EQ.1) THEN

C  Evaluate the second derivative of the spline at  xtab(j)  using old
C  values of  l .

                    Y2TAB(J)=SBL2(XTAB(J),N,L,X(IND),X(IND+1),
     *                            TB((N-2)*(N-1)/2),1,LAUX2)
                END IF

            END IF

            IND1=IND
            IND=0

30      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE SSTINF (OPT,D0,DNP,CONSTR,N,K,X,Y,D,NP,COMM,PART,EPS,
     *                     BETA,BETAI,DAUX2,DAUX3,INFO,ERRC,DIAGN)

C  SSTINF computes the information needed in the other parts of the
C  program using the data-dependent input parameters and stores it in
C  the output array INFO.
C
C  The parameters OPT,N,K,X,Y,D,NP,COMM,PART,EPS,BETA,BETAI,ERRC,DIAGN
C  are described in subroutine DBVSSC .
C
C  Items of possible interest are:
C
C  INFO  : floating array, of bounds 1:COMM+PART*NP+NP+1. It is composed
C          of four parts: the first, of bounds 1:comm, contains the
C          global information n, k , the maximum of the first divided 
C          differences of initial points and the lower and upper bounds
C          for the derivatives, bounds which are used when no 
C          constraints are imposed (see the parameter OPT described in
C          DBVSSC) or when the constraints must be relaxed; the second,
C          of bounds  comm+1:comm+np, contains information about 
C          constraints in the interval (x(i),x(i+1)) , i=0,1,...,np-1 ;
C          if:
C          info((comm+1)+i)= 0 - no attribute;
C          info((comm+1)+i)= 1 - the data are increasing;
C          info((comm+1)+i)= 2 - the data are decreasing;
C          info((comm+1)+i)=10 - the data are convex;
C          info((comm+1)+i)=11 - the data are increasing and convex;
C          info((comm+1)+i)=12 - the data are decreasing and convex;
C          info((comm+1)+i)=20 - the data are concave;
C          info((comm+1)+i)=21 - the data are increasing and concave;
C          info((comm+1)+i)=22 - the data are decreasing and concave.
C          The third part, of bounds comm+np+1:comm+part*np, contains
C          the first divided differences of initial points
C              ( y(i+1)-y(i) ) / ( x(i+1)-x(i) ) ,  i=0,1,...,np-1 .
C          The fourth, of bounds comm+part*np+1:comm+part*np+np+1,
C          contains, eventually, the initial estimates of the first 
C          derivatives which are then used to compute the constrained
C          best approximation (see the description of the input
C          parameter OPT  and of the array D in subr. DBVSSC). More
C          precisely, having defined  p := opt/100 , if p=2 it contains
C          the Bessel estimates, if p=3 it contains a set of third order
C          accurate estimates giving a co-monotone cubic Hermite
C          interpolant (see subr. STDC described later), if p=4 it 
C          contains a set of values given by the user; if p=1 this part
C          of INFO is not referenced.


        EXTERNAL BETA,BETAI

        INTEGER NP
        INTEGER OPT,CONSTR(0:NP-1),N,K,COMM,PART,ERRC,DIAGN(0:NP-1),
     *          I,R,Q,P,IND1,IND2,IND3

        REAL D0,DNP,X(0:NP),Y(0:NP),D(0:NP),EPS,BETA,BETAI,
     *       DAUX2(1:NP-1),DAUX3(0:NP-1),
     *       INFO(1:COMM+PART*NP+NP+1),D2IM1,D2I,IAUX0,
     *       IAUXNP

        SAVE FL0,FLEM4,FL2
        REAL FL0,FLEM4,FL2
        DATA FL0,FLEM4,FL2/0.0E00,1.0E-04,2.0E00/

        R=MOD(OPT,10)
        Q=MOD(OPT,100)/10
        P=OPT/100

        IND1=COMM+1
        IND2=COMM+NP+1
        IND3=COMM+2*NP+1

C  Set the first and the second components of INFO to n and k
C  respectively.

        INFO(1)=N
        INFO(2)=K

C  Compute the first divided differences of the initial points and set
C  info(3) to their maximum.

        INFO(3)=FL0
        DO 10 I=0,NP-1
            INFO(IND2+I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
            INFO(3)=MAX(INFO(3),ABS(INFO(IND2+I)))
10      CONTINUE

C  Compute the lower and upper bounds for derivatives

        INFO(4)=-N*INFO(3)/K
        INFO(5)=N*INFO(3)/K

C  If eps.LE.0 it is reset to flem4.

        IF (EPS.LE.FL0) EPS=FLEM4

C  Compute the relative tollerance of the method.

        EPS=EPS*INFO(3)

C  Set the second part of INFO. Firstly, all the components are
C  initialized with 0.

        DO 20 I=0,NP-1
            INFO(IND1+I)=0
20      CONTINUE

C  Monotonicity is required: check if the initial points are increasing
C  or decreasing in each interval ( x(i), x(i+1) ) , i=0,1,...,np-1 .

        IF (R.EQ.1.OR.R.EQ.3.OR.R.EQ.4) THEN
            DO 30 I=0,NP-1
                IF (INFO(IND2+I).GE.FL0) THEN
                    INFO(IND1+I)=INFO(IND1+I)+1
                ELSE
                    INFO(IND1+I)=INFO(IND1+I)+2
                END IF
30          CONTINUE
        END IF

C  Convexity is required: check if the initial points are concave or
C  convex in each interval ( x(i), x(i+1) ) , i=1,...,np-2 .

        IF (R.EQ.2.OR.R.EQ.3.OR.R.EQ.4) THEN
            DO 40 I=1,NP-2
                D2IM1=INFO(IND2+I)-INFO(IND2+(I-1))
                D2I=INFO(IND2+(I+1))-INFO(IND2+I)
                IF (D2IM1.GE.EPS.AND.D2I.GE.-EPS.OR.D2IM1.GE.-EPS.AND.
     *              D2I.GE.EPS.OR.ABS(D2IM1).LE.EPS.AND.
     *              ABS(D2I).LE.EPS)                      THEN
                    INFO(IND1+I)=INFO(IND1+I)+10
                ELSE IF (D2IM1.LE.-EPS.AND.D2I.LE.EPS.OR.D2IM1.LE.EPS
     *                   .AND.D2I.LE.-EPS)                THEN
                    INFO(IND1+I)=INFO(IND1+I)+20
                END IF
40          CONTINUE

C  The convexity in the first and in the last interval is defined as the
C  second and the second to last, respectively.

            INFO(IND1)=INFO(IND1)+(INT(INFO(IND1+1))/10)*10
            INFO(IND1+(NP-1))=INFO(IND1+(NP-1))+
     *                        (INT(INFO(IND1+(NP-2)))/10)*10
        END IF

C  In the case  r=4 , that is when the constraint selection
C  is made on any interval, we compare the kind given by the data with
C  those given by the array CONSTR

        IF (R.EQ.4) THEN
            CALL SMSK1(INFO,CONSTR,COMM,PART,IND1,NP)
        END IF

C  In the case q=3, the kind in the first and last subinterval 
C  is compared with the boundary conditions

        IF (Q.EQ.3) THEN
            CALL SMSK2(INFO,COMM,PART,IND1,NP,D0,DNP,EPS,
     *                       ERRC,DIAGN)
        END IF

C  If p=2 the Bessel derivatives are stored in the fourth
C  part of INFO.

        IF (P.EQ.2) THEN
            DO 50 I=1,NP-1
                INFO(IND3+I)=((X(I+1)-X(I))*INFO(IND2+(I-1))+
     *                       (X(I)-X(I-1))*INFO(IND2+I))/(X(I+1)-X(I-1))
50          CONTINUE

            IF (Q.EQ.1) THEN

C  If no boundary condition is imposed, set the first and last
C  derivatives using the standard formulas for Bessel interpolation.

                INFO(IND3)=((X(1)-X(0))*(FL2*INFO(IND2)-INFO(IND2+1))+
     *                     (X(2)-X(1))*INFO(IND2))/(X(2)-X(0))
                INFO(IND3+NP)=((X(NP)-X(NP-1))*(FL2*INFO(IND2+(NP-1))-
     *                        INFO(IND2+(NP-2)))+(X(NP-1)-X(NP-2))*
     *                        INFO(IND2+(NP-1)))/(X(NP)-X(NP-2))
            ELSE

C  Compute the first and last derivatives taking into account both the
C  slopes of the data and the restriction imposed by the boundary
C  conditions

                IAUX0=((X(1)-X(0))*(FL2*INFO(IND2)-INFO(IND2+1))+
     *                (X(2)-X(1))*INFO(IND2))/(X(2)-X(0))
                IAUXNP=((X(NP)-X(NP-1))*(FL2*INFO(IND2+(NP-1))-
     *                 INFO(IND2+(NP-2)))+(X(NP-1)-X(NP-2))*
     *                 INFO(IND2+(NP-1)))/(X(NP)-X(NP-2))

                INFO(IND3)=((X(NP)-X(NP-1))*IAUX0+
     *                     (X(1)-X(0))*BETAI(IAUXNP))/
     *                     ((X(1)-X(0))+(X(NP)-X(NP-1)))
                INFO(IND3+NP)=((X(NP)-X(NP-1))*BETA(IAUX0)+
     *                        (X(1)-X(0))*IAUXNP)/
     *                        ((X(1)-X(0))+(X(NP)-X(NP-1)))
            END IF

C  If p=3 then the set of third order accurate estimates, computed by
C  subr. STDC, is stored in the fourth part of INFO.

        ELSE IF (P.EQ.3) THEN
            CALL STDC(NP,X,COMM,PART,INFO,DAUX2,DAUX3)

C  If p=4 then the set of values given by the user is stored in
C  the fourth part of INFO.

        ELSE IF (P.EQ.4) THEN
            DO 60 I=0,NP
                INFO(IND3+I)=D(I)
60          CONTINUE
        END IF

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE STDC(NP,X,COMM,PART,INFO,DD2,DD3)

C  Given the initial points ( x(i), y(i) ) , i=0,1,...,np , STDC
C  computes a sequence  ds(0),ds(1),...,ds(np)  of estimates of the
C  function's derivatives which are third order accurate and furnish a
C  cubic Hermite interpolant preserving monotonicity.
C  The method is composed by the two following fundamental steps:
C  1 - compute an initial estimate of ds(i), i=0,1,...,np , which is
C      third or higher order accurate;
C  2 - refine it imposing the monotonicity constraint.
C  The computation of ds(i) needs the points x(i+j), j = -3,-2,...,3 ,
C  consequently, the boundary values ds(0), ds(1), ds(2) and ds(np-2),
C  ds(np-1), ds(np) are computed in an approximate way. Although they 
C  are still third order accurate, may not preserve the monotonicity.
C  For more details see \3\ .
C
C  The input parameter NP,X,COMM,PART are described in subr. DBVSSC; the
C  parameter INFO is described in subr. SSTINF .
C
C  The computed values are stored in the last part of INFO.


        INTEGER NP,I,IND,COMM,PART,IND1

        REAL X(0:NP),INFO(1:COMM+PART*NP+NP+1),DD2(1:NP-1),
     *       DD3(0:NP-1),Q1,Q2,TIT,F,P1,P2,TI,TMAX,TMIN,SI,
     *       Q12,Q32,E1,E2,E3,D1,D2,SMNMOD,SMDIAN

        SAVE FL0,FL2,FL3
        REAL FL0,FL2,FL3
        DATA FL0,FL2,FL3/0.0E00,2.0E00,3.0E00/

        IND=COMM+NP+1
        IND1=COMM+PART*NP+1

C  Compute the second divided differences of the initial points.

        DO 10 I=1,NP-1
             DD2(I)=(INFO(IND+I)-INFO(IND+I-1))/(X(I+1)-X(I-1))
10      CONTINUE

C  Compute the third divided differences of the initial points

        DO 20 I=1,NP-2
             DD3(I)=(DD2(I+1)-DD2(I))/(X(I+2)-X(I-1))
20      CONTINUE

C  Compute approximate values for  f[x(-1),x(0),x(1),x(2)]  and
C  f[x(np-2),x(np-1),x(np),x(np+1)] ; they are needed for the
C  computation of ds(2) and ds(np-2).

        DD3(0)=DD3(1)+((DD3(2)-DD3(1))/(X(4)-X(0)))*
     *               (X(0)+X(1)-X(2)-X(3))
        DD3(NP-1)=DD3(NP-2)+((DD3(NP-3)-DD3(NP-2))/(X(NP-4)-X(NP)))*
     *            (X(NP)+X(NP-1)-X(NP-2)-X(NP-3))

        DO 50 I=2,NP-2

C  ds(i) : initialization

             E1=SMNMOD(DD3(I-2),DD3(I-1))
             E2=SMNMOD(DD3(I-1),DD3(I))
             E3=SMNMOD(DD3(I),DD3(I+1))
             Q1=INFO(IND+I-1)+(X(I)-X(I-1))*SMNMOD(DD2(I-1)+E1*
     *          (X(I)-X(I-2)),DD2(I)+E2*(X(I)-X(I+1)))
             Q2=INFO(IND+I)-(X(I+1)-X(I))*SMNMOD(DD2(I)+E2*
     *          (X(I)-X(I-1)), DD2(I+1)+E3*(X(I)-X(I+2)))
             F=(Q1+Q2)/FL2

C  Refinement

             TIT=SMNMOD(Q1,Q2)
             D1=SMNMOD(DD2(I-1),DD2(I))
             D2=SMNMOD(DD2(I),DD2(I+1))
             P1=INFO(IND+I-1)+D1*(X(I)-X(I-1))
             P2=INFO(IND+I)+D2*(X(I)-X(I+1))
             TI=SMNMOD(P1,P2)
             SI=SMNMOD(INFO(IND+I-1),INFO(IND+I))
             TMIN=MIN(FL0,FL3*SI,FL3*TI/FL2,TIT)
             TMAX=MAX(FL0,FL3*SI,FL3*TI/FL2,TIT)
             INFO(IND1+I)=F+SMNMOD(TMIN-F,TMAX-F)
50      CONTINUE

C  ds(1): initialization

        Q12=INFO(IND)+DD2(1)*(X(1)-X(0))+DD3(0)*(X(1)-X(0))*(X(1)-X(2))
        Q32=INFO(IND)+DD2(1)*(X(1)-X(0))+DD3(1)*(X(1)-X(0))*(X(1)-X(2))
        E1=SMNMOD(DD3(0),DD3(1))
        E2=SMNMOD(DD3(1),DD3(2))
        Q1=SMDIAN(INFO(IND),Q12,Q32)
        Q2=INFO(IND+1)-(X(2)-X(1))*SMNMOD(DD2(1)+E1*(X(1)-X(0)),
     *                                    DD2(2)+E2*(X(1)-X(3)))
        F=(Q1+Q2)/FL2

C  refinement

        TIT=SMNMOD(Q1,Q2)
        D2=SMNMOD(DD2(1),DD2(2))
        P1=INFO(IND)+DD2(1)*(X(1)-X(0))
        P2=INFO(IND+1)+D2*(X(1)-X(2))
        TI=SMNMOD(P1,P2)
        SI=SMNMOD(INFO(IND),INFO(IND+1))
        TMIN=MIN(FL0,FL3*SI,FL3*TI/FL2,TIT)
        TMAX=MAX(FL0,FL3*SI,FL3*TI/FL2,TIT)
        INFO(IND1+1)=F+SMNMOD(TMIN-F,TMAX-F)

C  ds(np-1): initialization

        E1=SMNMOD(DD3(NP-3),DD3(NP-2))
        E2=SMNMOD(DD3(NP-2),DD3(NP-1))
        Q1=INFO(IND+NP-2)+(X(NP-1)-X(NP-2))*SMNMOD(DD2(NP-2)+
     *           E1*(X(NP-1)-X(NP-3)),DD2(NP-1)+E2*(X(NP-1)-X(NP)))
        Q12=INFO(IND+NP-2)+DD2(NP-1)*(X(NP-1)-X(NP-2))+
     *      DD3(NP-1)*(X(NP-1)-X(NP-2))*(X(NP-1)-X(NP))
        Q32=INFO(IND+NP-2)+DD2(NP-1)*(X(NP-1)-X(NP-2))+
     *      DD3(NP-2)*(X(NP-1)-X(NP-2))*(X(NP-1)-X(NP))
        Q2=SMDIAN(INFO(IND+NP-1),Q12,Q32)
        F=(Q1+Q2)/FL2

C  Refinement

        TIT=SMNMOD(Q1,Q2)
        D1=SMNMOD(DD2(NP-2),DD2(NP-1))
        P1=INFO(IND+NP-2)+D1*(X(NP-1)-X(NP-2))
        P2=INFO(IND+NP-1)+DD2(NP-1)*(X(NP-1)-X(NP))
        TI=SMNMOD(P1,P2)
        SI=SMNMOD(INFO(IND+NP-2),INFO(IND+NP-1))
        TMIN=MIN(FL0,FL3*SI,FL3*TI/FL2,TIT)
        TMAX=MAX(FL0,FL3*SI,FL3*TI/FL2,TIT)
        INFO(IND1+NP-1)=F+SMNMOD(TMIN-F,TMAX-F)

C  ds(0):

        Q1=INFO(IND)+DD2(1)*(X(0)-X(1))+DD3(0)*(X(0)-X(1))*(X(0)-X(2))
        Q2=INFO(IND)+DD2(1)*(X(0)-X(1))+DD3(1)*(X(0)-X(1))*(X(0)-X(2))
        INFO(IND1)=SMDIAN(INFO(IND),Q1,Q2)

C  ds(np):

        Q1=INFO(IND+NP-1)+DD2(NP-1)*(X(NP)-X(NP-1))+
     *     DD3(NP-1)*(X(NP)-X(NP-2))*(X(NP)-X(NP-1))
        Q2=INFO(IND+NP-1)+DD2(NP-1)*(X(NP)-X(NP-1))+
     *     DD3(NP-2)*(X(NP)-X(NP-2))*(X(NP)-X(NP-1))
        INFO(IND1+NP)=SMDIAN(INFO(IND+NP-1),Q1,Q2)

        RETURN
        END

C  ---------------------------------------------------------------------

        SUBROUTINE STRMB(N,TB)

C  STRMB    computes the binomial terms
C      i!/(k!*(i-k)!) , i=1,2,...,n , k=0,1,...,i .
C
C  INPUT PARAMETERS
C
C  N     : integer variable, containing the largest binomial term
C          needed.
C
C  OUTPUT PARAMETERS
C
C  TB    : floating array, of bounds  1:N*(N+1)/2+N  , containing
C          the values   i!/(k!*(i-k)!)  , k=0,1,...,i , in the
C          elements   TB(i*(i+1)/2),...,TB((i*(i+1)/2)+i) .


        INTEGER N,I,K,IND,IND1

        REAL TB(1:N*(N+1)/2+N),FL1

        SAVE FL1
        DATA FL1/1.0E00/

        TB(1)=FL1
        TB(2)=FL1
        IND=1

        DO 20 I=2,N
            IND1=IND
            IND=I*(I+1)/2
            TB(IND)=FL1
            TB(IND+I)=FL1

            DO 10 K=1,I-1
                TB(IND+K)=TB(IND1+K)+TB((IND1+K)-1)
10          CONTINUE

20      CONTINUE

        RETURN
        END

C  ---------------------------------------------------------------------

        LOGICAL FUNCTION STST(A,B,C,D,EPS)

C  STST checks if two intervals [a,b] and [c,d] differ less than eps.
C  STST assumes  a.LE.b  and  c.LE.d .


        REAL A,B,C,D,EPS

        STST=ABS(A-C).LE.EPS.AND.ABS(B-D).LE.EPS

        RETURN
        END

