C......................................................................
C     Single Precision
C
C         L2WPMA: A FORTRAN 77 PACKAGE FOR WEIGHTED LEAST SQUARES
C                 PIECEWISE MONOTONIC DATA APPROXIMATION
C                              by
C                Ioannes C. DEMETRIOU, PhD (Cantab)
C
C     Home address: 65 Othonos Street, Kifissia 14561, Athens, Greece
C     Work address: University of Athens, Department of Economics
C                   8 Pesmazoglou Street, Athens 10559, Greece
C     Tel   : +30-210-8017732, 3689817
C     E-mail: demetri@econ.uoa.gr
C     URL: www.econ.uoa.gr
C......................................................................
C
C     D I R E C T O R Y  O F  S U B R O U T I N E S
C
C       L2WPMA  Interface to the user (main subroutine), calculates
C                a weighted least squares piecewise monotonic
C                approximation to a set of data.
C       L2WMON  Calculates the weighted least squares monotonic
C                approximation to a set of data.
C       MESSGW  Contains the messages of the calculation.
C       TRIVIA  Examines some trivial cases of data points, when
C                entering L2WPMA.
C       XTREMA  Calculates the local extrema of the data points.
C
C......................................................................
C     Given the noisy function measurements F(i), i=I1,...,N, with
C      corresponding positive weights WF(i), at strictly increasing
C      abscissae X(i), and a positive integer KSECTN, usually smaller
C      than N-I1+1, subroutine L2WPMA calculates integers t(j),j=0,1,
C      ...,KSECTN that satisfy the conditions
C                   I1=t(0)<=t(1)<=...<=t(KSECTN)=N
C      and an (N-I1+1)-component vector Y(.) that minimizes the
C      objective function
C                  Sum WF(i)*(Y(i)-F(i))^2, i=I1,...,N
C      such that the inequalities
C
C      Y(t(0)) <= Y(t(0)+1) <= ... <= Y(t(1))
C      Y(t(1)) >= Y(t(1)+1) >= ... >= Y(t(2))
C      ...
C      Y(t(KSECTN-1)) >(<)= Y(t(KSECTN-1)+1) >(<)= ... >(<)= Y(N)
C
C      are satisfied. It is important to note that the t(i)s are also
C      variables of the minimization calculation.
C
C     In this sense Y(.) is the weighted least squares piecewise
C      monotonic fit to F(.) with KSECTN monotonic sections.
C
C     INFORMATION: The user has some options in supplying the weights
C      and the abscissae together with the measurements (see MODEWF
C      in Input arguments below).
C
C.... REFERENCES
C     1. Demetriou, I.C. "L2WPMA: A Fortran 77 package for weighted
C          least squares piecewise monotonic data smoothing",
C          Report Department of Economics, University of Athens 2004.
C     2. Demetriou, I.C. "Discrete piecewise monotonic approximation
C          by a strictly convex distance function", Math. of
C          Computation, 64(1995), 157-180.
C     3. Demetriou, I.C. and M.J.D. Powell "Least squares smoothing
C          of univariate data to achieve piecewise monotonicity",
C          IMA J. Numer. Anal., 11(1991), 411-432.
C
C     Calls subroutines TRIVIA, XTREMA, L2WMON and MESSGW.
C......................................................................
      SUBROUTINE L2WPMA(I1,N,X,F,WF,MODEWF,KSECTN,IORDER,Y,WY,NK,IAKN,
     +                  NACT,IACT,PAR,ITAU,ITHETA,MODE,SS,G,RG,LOWER,
     +                  IUPPER,INDX,FT,WFT,FTNEG,WY1,Z,WZ,IW,IAKNW)
C
C
C.... I N P U T (They must be set by the user) ....
C     I1        Integer variable, lower data index, I1 = 1.
C     N         Integer variable, upper data index, N >= I1.
C     X(I1:N)   Real array of abscissae X(i), i=I1,...,N. The
C                abscissae are not directly needed in the required
C                calculation. There may be supposed, without loss of
C                generality, that they are the numbers 1,2,...,N taken
C                in their natural order. Further, the X(i)s may be
C                used in defining the weights (see input var MODEWF
C                below).
C     F(I1:N)   Real array of function measurements F(i), i=I1,...,N,
C                at strictly increasing abscissae.
C     WF(I1,N)  Real array of positive weights WF(i), i=I1,...,N,
C                associated with the function measurements
C                F(i), i=I1,...,N
C     MODEWF    Integer variable, specifying the weights WF(.) as
C			   follows:
C                MODEWF = 0, User defined weights, the default value.
C                MODEWF = 1, Reset weights to unity, ie. WF(i) = 1, for
C                             i=I1,...,N.
C                MODEWF = 2, Reset weights so as to reflect possible
C                             differences in abscissae spacing, ie.
C                             WF(i) = w(i)/(Sum of w(i), i=I1,...,N),
C                                                       for i=I1,...,N,
C                             where
C                             w(i) = 1/(X(i) - X(i-1)), i=I1+1,...,N
C                             and
C                             w(I1)=(Sum of w(i), i=I1+1,...,N)/(N-I1).
C                MODEWF = 3, Reset weights so as to reflect possible
C                             differences in abscissae spacing, ie.
C                             WF(i) = w(i)/(Sum of w(i), i=I1,...,N),
C                                                       for i=I1,...,N,
C                             where
C                             w(i) = 1/(X(i+1) - X(i)), i=I1,...,N-1
C                             and
C                             w(N)=(Sum of w(i), i=I1,...,N-1)/(N-I1).
C                             The option MODEWF=3 permutes cyclically
C                             the weights obtained by MODEWF=2 (ie.
C                             i --> i-1, for i= 2,...,n and 1 --> n).
C                MODEWF = 4, Reset weights so as to reflect possible
C                             differences in abscissae spacing, ie.
C                                     (X(I1+1)-X(I1))/2,       i=I1
C                             WF(i) = (X(i+1)-2X(i)+X(i-1))/2,
C                                                i=I1+1,I1+2,...,N-1
C                                     (X(N)-X(N-1))/2,         i=N
C                MODEWF = 5, Only F is given. The components of X are
C                             set to X(i)=i, i=I1,...,N, and the
C                             components of WF are set to WF(i)=1, i=
C                             I1,...,N.
C     KSECTN    Integer variable, required number of monotonic sections
C                in the optimal fit, such that 1 <= KSECTN. In certain
C                cases described in the paragraph under the heading
C                Output, KSECTN acts as an output variable.
C     IORDER    Integer variable, defining the order of first monotonic
C                section of optimal fit:
C                IORDER =  0, first section of the fit is increasing,
C                              the default value.
C                IORDER = -1, first section of the fit is decreasing.
C
C.... O U T P U T ....
C     KSECTN    If the data provides the required approximation with M
C                monotonic sections, where M < KSECTN (user defined),
C                then the package sets KSECTN = M. This setting is
C                necessary for presenting some output messages and
C                defining array ITHETA below. Note that the calculation
C                stops immediately after this setting.
C     Y(I1:N)   Real array containing the best piecewise monotonic fit
C                to F(.) at the end of the calculation.
C     WY(I1:N)  Real array containing in WY(i), i=I1,...,N, the weights
C                of Y(i), i=I1,...,N, respectively. They may be useful
C                to further processing of the smoothed values Y(.).
C     NK        Integer variable, the number of knots of spline
C                representation of Y(.) at the end of the calculation,
C                I1 <= NK <= N-I1.
C     IAKN(I1-1:N-I1+1) Integer array containing the knot indices of
C                spline representation of Y(.) at the end of the
C                calculation. Specifically, these indices are {IAKN(i),
C                i=I1,I1+1,...,NK}, where NK may never be zero. The
C                element IAKN(I1-1) acts as a switch in subroutine
C                L2WMON.
C     NACT      Integer variable, the number of constraints that are
C                satisfied as equations at the end of the calculation,
C                where 0 <= NACT <= N-I1.
C     IACT(1:NACT) Integer array that provides the indices of the
C                constraints satisfied as equations at the end of the
C                calculation.
C     PAR(I1+1:N) Real array containing the Lagrange parameters
C                associated with the constraints that are satisfied as
C                equations by Y(.) at the end of the calculation.
C                Specifically these parameters are PAR(IACT(k)), k=1,
C                2,...,NACT.
C     ITAU(0:KSECTN,I1:(N - I1 + 1)/2 + 1)
C               Integer array, such that ITAU(m,J) is the index of the
C                rightmost extremum of an optimal fit with m monotonic
C                sections to the first LOWER(J) data or IUPPER(J) data
C                where, 1 <= m <= KSECTN, I1 <= J< = (N -I1+1)/2+1, the
C                arrays LOWER and IUPPER being defined below.
C     ITHETA(0:KSECTN)
C               Integer array that holds in ITHETA(m) the greatest
C                value of ITAU(m,j) that has already been calculated.
C                At the end of the calculatn, ITHETA holds the optimal
C                values of the integer variables t(i), i=0,...,KSECTN.
C     MODE      Integer variable indicating the status of subroutine
C                termination:
C               0 = Unsuccessful return of L2WPMA, if KSECTN < 1.
C               1 = Successful return of L2WPMA.
C               2 = Successful return of L2WPMA due to that KSECTN is
C                    greater than the number of monotonic sectns in
C                    data. Hence the data itself is the required
C                    approximation.
C
C.... M E S S A G E S ....
C     Included in Subroutine MESSGW.
C
C.... W O R K I N G  S P A C E ....
C     SS(I1:N)  Real array, defined in subroutine L2WMON.
C     SS(N)     Contains the weighted sum of the squares of the
C                residuals (Y(i)-F(i)) at the end of the calculation,
C                when KSECTN >= 1. We call this sum, 'error'.
C     G(0:KSECTN,I1:N2)  Real array, keeps in G(m,j) the error of the
C                fit associated with ITAU(m,j), where
C                    N2 = (N - I1 + 1)/2 + 1
C                is the worst possible number of local minima (maxima)
C                in the data.
C     RG(I1:N2) Real array, providing in RG(j) temporary storage of
C                the error of a trial fit with m sections to the
C                first j data as m=0,1,..., KSECTN.
C     LOWER(I1:N2) Integer array, defined in subroutine XTREMA.
C     IUPPER(I1:N2) Integer array, defined in subroutine XTREMA.
C     INDX(I1:N) Integer array, defined in subroutine XTREMA.
C     FT(I1:N)  Real array, keeps the elements of F(.) in reverse
C                order.
C     WFT(I1:N) Real array, keeps the elements of WF(.) in reverse
C                order.
C     FTNEG(I1:N) Real array, keeps the negative of the elements
C                of FT(.).
C     WY1(I1:N) Real array, provides space for calculating the
C                weights of the optimal approximation.
C     Z(1:N + I1 - 1) Real array, defined in subroutine L2WMON.
C     WZ(1:N + I1 -1) Real array, defined in subroutine L2WMON.
C     IW(1:N + I1 - 1) Integer array, defined in subroutine L2WMON.
C     IAKNW(I1:N) Integer array that provides working space for IAKN.
C
C.... M E T H O D ....
C     It is described in the quoted references. It implements
C      Algorithm 1 of Demetriou (1995) except that L2WPMA allows the
C      data have positive weights and at the end of the calculation
C      it provides a spline representation of the solution and the
C      corresponding Lagrange multipliers.
C
C.... I N I T I A L I Z A T I O N .....................................
C
C     .. Scalar Arguments ..
      INTEGER I1,IORDER,KSECTN,MODE,MODEWF,N,NACT,NK
C     ..
C     .. Array Arguments ..
      REAL X(I1:N),F(I1:N),WF(I1:N),Y(I1:N),WY(I1:N),SS(I1:N),
     +     G(0:KSECTN,I1:(N - I1 + 1)/2 + 1),RG(I1:(N - I1 + 1)/2 + 1),
     +     FT(I1:N),WFT(I1:N),FTNEG(I1:N),WY1(I1:N),Z(1:N + I1 - 1),
     +     WZ(1:N + I1 - 1),PAR(I1:N)
C
      INTEGER ITAU(0:KSECTN,I1:(N - I1 + 1)/2 + 1),ITHETA(0:KSECTN),
     +     LOWER(I1:(N - I1 + 1)/2 + 1),IUPPER(I1:(N - I1 + 1)/2 + 1),
     +     INDX(I1:N),IW(1:N + I1 - 1),IAKN(I1-1:N - I1 + 1),
     +     IAKNW(1:N - I1 + 1),IACT(1:N-I1)
C     ..
C     .. Local Scalars ..
      REAL ANORM2,SUMWF
C
      INTEGER I,II,IL,ILU,IU,J,K,KSECT1,L1,LN,M,MINRG,N2,NL,NU
C     ..
C     .. External Subroutines ..
      EXTERNAL L2WMON,MESSGW,TRIVIA,XTREMA
C     ..
C
C                               Initialize NACT, NK and IAKN(I1-1)
      NACT = 0
      NK = I1-1
      IAKN(I1-1) = 1
C      
      IF (IORDER.EQ.-1) THEN
C                               MSG: First monotonic sectn is required
C                                     to be decreasing.
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,100)
C
C                               Reverse sign of F(.) if first monotonic
C                                section is required to be decreasing.
          DO 10 I = I1,N
             F(I) = -F(I)
   10     CONTINUE
C
      END IF
C
      NL = N - I1 + 1
      NU = 0
C                               The values of NL and NU will change by
C                                the call of XTREMA. Meanwhile, NL and
C                                NU may be used by MESSGW to indicate
C                                termination due to TRIVIAl cases.
C
C                               Reset the weights.
C
      IF (MODEWF.EQ.1) THEN
C                               The weights are set or reset to 1.0.
          DO 12 I = I1,N
             WF(I) = 1.0
   12     CONTINUE
      END IF
C                               If the data are not uniformly spaced
C                                it might be preferable to define (or
C                                reset) weights that are dependent on
C                                the data spacing (MODEWF=2,3,4).
      IF (MODEWF.EQ.2) THEN
          IF (N.EQ.I1) GOTO 18
C
C                               MSG: Automatic weight generation, due
C                                     to abscissae spacing (see also
C                                     alternative weightings below).
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,110)
C
          SUMWF = 0.0
          DO 15 I = I1 + 1,N
             WF(I) = 1/(X(I) - X(I - 1))
             SUMWF = SUMWF + WF(I)
   15     CONTINUE
          WF(I1) = SUMWF/(N - I1)
C                               Normalize the weights.
          SUMWF = SUMWF + WF(I1)
          DO 17 I = I1 ,N
             WF(I) = WF(I)/SUMWF
   17     CONTINUE
C
   18     CONTINUE
      END IF
C
      IF (MODEWF.EQ.3) THEN
          IF (N.EQ.I1) GOTO 25
C
C                               MSG: Automatic weight generation, due
C                                     to abscissae spacing.
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,110)
C
          SUMWF = 0.0
          DO 20 I = I1,N-1
             WF(I) = 1/(X(I+1) - X(I))
             SUMWF = SUMWF + WF(I)
   20     CONTINUE
          WF(N) = SUMWF/(N - I1)
C                               Normalize the weights.
          SUMWF = SUMWF + WF(N)
          DO 22 I = I1,N
             WF(I) = WF(I)/SUMWF
   22     CONTINUE
C
   25     CONTINUE
      END IF
C
C
      IF (MODEWF.EQ.4) THEN
          IF (N.EQ.I1) GOTO 27
C
C                               MSG: Automatic weight generation, due
C                                     to abscissae spacing.
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,110)
C
          WF(I1)=(X(I1+1) - X(I1))/2
          DO 26 I = I1+1,N-1
             WF(I) = (X(I+1) - X(I-1))/2
   26     CONTINUE
          WF(N)=(X(N) - X(N-1))/2
C
   27     CONTINUE
      END IF
C
C                               Set abscissae and weights, because
C                                only F is provided.
      IF (MODEWF.EQ.5) THEN
          DO 28 I = I1,N
             X(I) = I
             WF(I) = 1.0
   28     CONTINUE
      END IF
C
C                               Set some variables and arrays.
      MODE = 0
      SS(N) = 0.0
      N2 = (N - I1 + 1)/2 + 1
C                               We note that at the beginning of the
C                                calculation of the optimal integers
C                                (backwards tracing, below), N2 is set
C                                either to NL or to NU.
      DO 30 I = I1,N
         FT(I) = F(N + I1 - I)
         WFT(I) = WF(N + I1 - I)
         FTNEG(I) = -F(N + I1 - I)
         PAR(I) = 0.0
   30 CONTINUE
C                               Check on trivial cases.
      CALL TRIVIA(I1,N,F,KSECTN,Y,ITHETA,MODE,NL,NU)
C
      IF (MODE.EQ.0 .OR. MODE.EQ.2) THEN
C                               MSG: Terminate due to subrtn TRIVIA.
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,120)
          GO TO 410
      END IF
C
C.... M A I N  C A L C U L A T I O N ...................................
C
C
C                               Find the indices of lower and upper data
C                                values.
      CALL XTREMA(F,LOWER,IUPPER,INDX,KSECTN,I1,N,N2,NL,NU)
C
      IF (MODE.EQ.-1) THEN
C         The following assignment is explained inside L2WMON.
          IAKN(I1-1) = 1
C                               Calculate best weighted monotonic
C                                increasing approximation.
          CALL  L2WMON(SS,Y,WY,F,WF,Z,WZ,IW,IAKNW,I1,N,I1,N,1,NK,IAKN,
     +                 NACT,IACT,PAR)
          GO TO 320
      END IF
C
      IF (NL + NU - 2*I1 + 2 .LE. KSECTN + 1) THEN
C                               Terminate the calculation, because the
C                                data provide the reqired optimal fit.
          DO 35 I = I1,N
             Y(I) = F(I)
             WY(I) = WF(I)
   35     CONTINUE
C
C     IL  counter on LOWER.
C     IU  counter on IUPPER.
C
C                               Store the values of integer variables.
          IL = I1
          IU = I1
          ITHETA(0) = IL
          IL = IL + 1
          DO 40 J = 1,KSECTN
             IF (J.EQ.J/2*2) THEN
                 ITHETA(J) = LOWER(IL)
                 IL = IL + 1
             ELSE
                 ITHETA(J) = IUPPER(IU)
                 IU = IU + 1
             END IF
             IF (ITHETA(J).EQ.N) GO TO 50
   40     CONTINUE
C                               Set NK and store indices of knots.
   50     NK=N
          DO 55 I=I1,NK
             IAKN(I)=I
   55     CONTINUE
C                               Set number of active sonstraints.
          NACT=0
C                               Redefine KSECTN (and stop).
          KSECT1 = KSECTN
          KSECTN = J
C                               MSG: Data itself is the best piecewise
C                                     monotonic approximation, because
C                                     number of required monotonic
C                                     sections is either equal to, or
C                                     it exceeds number of monotonic
C                                     sections in data.
          MODE = 2
          CALL MESSGW(MODE,I1,N,KSECT1,NL,NU,130)
          GO TO 390
C
      END IF
C
      DO 60 M = 0,KSECTN
         ITHETA(M) = I1
   60 CONTINUE
C                               Calculate best weighted monotonic incrg
C                                appxtns sequentially to all data.
      CALL L2WMON(SS,Y,WY,F,WF,Z,WZ,IW,IAKNW,I1,N,I1,N,0,NK,IAKN,NACT,
     +            IACT,PAR)
C                               Initialize G(1,.), ITAU(1,.) for IUPPER
      DO 70 J = I1,NU
         G(1,J) = SS(IUPPER(J))
         ITAU(1,J) = I1
   70 CONTINUE
C                               Initialize counters.
      IL = I1
      IU = I1
C                               ILU data counter on LOWER and IUPPER.
      ILU = IUPPER(IU)
C                               Initialize G(.,I1), ITAU(.,I1) for M
C                                in [1,KSECTN-1].
      DO 80 M = 1,KSECTN - 1,2
         G(M,IU) = 0.
         ITAU(M,IU) = I1
   80 CONTINUE
C
      DO 85 M = 2,KSECTN - 1,2
         G(M,IL) = 0.
         ITAU(M,IL) = I1
   85 CONTINUE
C                               Initialize ITAU(0,.) for LOWER.
      DO 90 I = I1,NL
         ITAU(0,I) = I1
   90 CONTINUE
C
C     If LOWER has been exhausted or if KSECTN=2, branch to label 250.
C
  100 CONTINUE
      IF (ILU.GE.LOWER(NL) .OR. KSECTN.EQ.2) GO TO 250
      IL = IL + 1
      ILU = LOWER(IL)
C
C     Calculate best weighted monotonic decreasing approximations to
C      data on [ITHETA(2),ILU].
C
      CALL L2WMON(SS,Y,WY,FT,WFT,Z,WZ,IW,IAKNW,I1,N,N+I1-ILU,
     +            N+I1-ITHETA(2),0,NK,IAKN,NACT,IACT,PAR)
C
C     Calculate G(M,ILU), for even M in [1,KSECTN-1].
C
      DO 160 M = 2,KSECTN - 1,2
C                               Find lower bound on integer variable.
          IF (ITAU(M-2,IL).LT.ITHETA(M)) THEN
              K = INDX(ITHETA(M))
          ELSE
              K = INDX(ITAU(M - 2,IL))
          END IF
C
          DO 130 J = K,IU
             RG(J) = G(M - 1,J) + SS(N + I1 - IUPPER(J))
  130     CONTINUE
C                               Find position of least of RG(K:IU).
          MINRG = K
          IF (K.GE.IU) GO TO 150
          DO 140 II = K + 1,IU
             IF (RG(II).LT.RG(MINRG)) THEN
                 MINRG = II
             END IF
  140     CONTINUE
  150     G(M,IL) = RG(MINRG)
          ITAU(M,IL) = IUPPER(MINRG)
          ITHETA(M) = IUPPER(MINRG)
  160 CONTINUE
C
C     If IUPPER has been exhausted, branch to label 250.
C
      IF (ILU.GE.IUPPER(NU)) GO TO 250
      IU = IU + 1
      ILU = IUPPER(IU)
C
C     Calculate best weighted monotonic decreasing approximations to
C      data on [ITHETA(3),ILU].
C
      CALL L2WMON(SS,Y,WY,FTNEG,WFT,Z,WZ,IW,IAKNW,I1,N,N+I1-ILU,
     +            N+I1-ITHETA(3),0,NK,IAKN,NACT,IACT,PAR)
C
C     Calculate G(M,ILU) for odd M in [3,KSECTN-1].
C
      DO 230 M = 3,KSECTN - 1,2
C                               Find lower bound on integer variable.
          IF (ITAU(M - 2,IU).LT.ITHETA(M)) THEN
              K = INDX(ITHETA(M))
          ELSE
              K = INDX(ITAU(M - 2,IU))
          END IF
C
          DO 200 I = K,IL
              RG(I) = G(M - 1,I) + SS(N + I1 - LOWER(I))
  200     CONTINUE
C                               Find position of least of RG(K:IL).
          MINRG = K
          IF (K.GE.IL) GO TO 220
          DO 210 II = K + 1,IL
             IF (RG(II).LT.RG(MINRG)) THEN
                 MINRG = II
             END IF
  210     CONTINUE
  220     G(M,IU) = RG(MINRG)
          ITAU(M,IU) = LOWER(MINRG)
          ITHETA(M) = LOWER(MINRG)
C
  230 CONTINUE
C
C     Branch to label 100.
C
      GO TO 100
C
  250 CONTINUE
      ILU = N
      IF (KSECTN.EQ.KSECTN/2*2) THEN
          K = ITAU(KSECTN - 2,IL)
          CALL L2WMON(SS,Y,WY,FT,WFT,Z,WZ,IW,IAKNW,I1,N,I1,N + I1 - K,
     +                0,NK,IAKN,NACT,IACT,PAR)
C
C         Calculate G(KSECTN,NL), for even KSECTN.
C
          K = INDX(K)
          DO 260 J = K,NU
             RG(J) = G(KSECTN - 1,J) + SS(N + I1 - IUPPER(J))
  260     CONTINUE
C                               Find position of least of RG(K:NU).
          MINRG = K
          IF (K.GE.NU) GO TO 280
          DO 270 II = K + 1,NU
             IF (RG(II).LT.RG(MINRG)) THEN
                 MINRG = II
             END IF
  270     CONTINUE
  280     G(KSECTN,NL) = RG(MINRG)
          ITAU(KSECTN,NL) = IUPPER(MINRG)
          ITHETA(KSECTN) = IUPPER(MINRG)
      ELSE
          K = ITAU(KSECTN - 2,IU)
          CALL L2WMON(SS,Y,WY,FTNEG,WFT,Z,WZ,IW,IAKNW,I1,N,I1,
     +                N + I1 - K,0,NK,IAKN,NACT,IACT,PAR)
C
C         Calculate G(KSECTN,NU) for odd KSECTN.
C
          K = INDX(K)
          DO 290 I = K,NL
             RG(I) = G(KSECTN - 1,I) + SS(N + I1 - LOWER(I))
  290     CONTINUE
C                               Find position of least of RG(K:NL).
          MINRG = K
          IF (K.GE.NL) GO TO 310
          DO 300 II = K + 1,NL
             IF (RG(II).LT.RG(MINRG)) THEN
                 MINRG = II
             END IF
  300     CONTINUE
  310     G(KSECTN,NU) = RG(MINRG)
          ITAU(KSECTN,NU) = LOWER(MINRG)
          ITHETA(KSECTN) = LOWER(MINRG)
C
      END IF
C
C.... Calculate optimal sequence of integer variables..................
C
      MINRG = NU
      IF (KSECTN.EQ.KSECTN/2*2) MINRG = NL
      N2 = MINRG
C                               Use ITHETA(.) to store the values of
C                                the optimal integers.
  320 ITHETA(0) = I1
      ITHETA(KSECTN) = N
C
C                               Termination with one monotonic section
      IF (MODE.EQ.-1) THEN
C                               Set last index knot at N.
         IF (IAKN(NK).LT.N) THEN
C                               The condition "IAKN(NK) .LT. N" is
C                                imposed in order to avoid the case
C                                that a knot at N is counted twice
C                                whenever Y(N-1)<Y(N) occurs,because
C                                IAKN(NK) has already been set at N by
C                                subroutine L2WMON.
             NK = NK + 1
             IAKN(NK) = N
          END IF
C
          MODE = 1
          GO TO 390
      END IF
C                               Backwards tracing.
      ITHETA(KSECTN - 1) = ITAU(KSECTN,MINRG)
      MINRG = INDX(ITHETA(KSECTN - 1))
      DO 330 J = KSECTN - 2,0,-1
         ITHETA(J) = ITAU(J + 1,MINRG)
         MINRG = INDX(ITHETA(J))
  330 CONTINUE
C
C......................................................................
      ANORM2 = 0.
C                               Calculate the optimal approximation
C                                section by section.
C                               Set Y(i)=F(i) at i=ITHETA(.), i.e. Y(i)
C                                at positions of optimal turning points
      DO 350 I = 1,KSECTN - 1
         II = ITHETA(I)
         Y(II) = F(II)
         WY(II) = WF(II)
  350 CONTINUE
C                               Variable J counts the monotonic sectns.
      J = 1
C                               Initiate knot counter of spline
C                                representation of optimal fit.
      NK = I1-1
C	 
  360 L1 = ITHETA(J - 1) + 1
      IF (J.EQ.1) L1 = I1
      LN = ITHETA(J) - 1
      IF (J.EQ.KSECTN) LN = N
C                               Calculate monotonic section on [L1,LN].
      IF (L1.GT.LN) GO TO 380
C
      IF (J.NE.J/2*2) THEN
C
C         The following assignment is explained inside L2WMON.
          IAKN(I1-1) = 1
C
          CALL L2WMON(SS,Y,WY,F,WF,Z,WZ,IW,IAKNW,I1,N,L1,LN,1,NK,IAKN,
     +                NACT,IACT,PAR)
          ANORM2 = ANORM2 + SS(LN)
      ELSE
C                               Some explanations for the calling
C                                sequence of subr L2WMON below:
C                               FTNEG provides auxiliary space for
C                                calculating the components of
C                                Y(L1:LN).
C                               WY1 provides space for calculating
C                                the components of WY(L1:LN).
C
C         The following assignment is explained inside L2WMON.
          IAKN(I1-1) = -1
C
          CALL L2WMON(SS,FTNEG,WY1,FT,WFT,Z,WZ,IW,IAKNW,I1,N,
     +                N + I1 - LN,N + I1 - L1,1,NK,IAKN,NACT,IACT,PAR)
C
          ANORM2 = ANORM2 + SS(N + I1 - L1)
C                               Set optimal components.
          DO 370 I = L1,LN
             Y(I) = FTNEG(N + I1 - I)
             WY(I) = WY1(N + I1 - I)
  370     CONTINUE
C
      END IF
C                               Terminate if J=KSECTN.
  380 CONTINUE
C
      IF (J.LT.KSECTN) THEN
C                               Update IAKN. A knot is defined at a
C                                turning point index.
          NK = NK + 1
          IAKN(NK) = LN + 1
C
          J = J + 1
          GO TO 360
      END IF
C
      IF (J.EQ.KSECTN .AND. IAKN(NK).LT.N) THEN
C                               Set last index knot at LN=N. The
C                                condition ".AND. IAKN(NK).LT.N" is
C                                imposed in order to avoid the case
C                                that a knot at N is counted twice
C                                whenever a turning point occurs at N,
C                                because IAKN(NK) has already taken its
C                                value in the previous IF statement.
          NK = NK + 1
          IAKN(NK) = LN
      END IF
C                               Keep the weighted sum of the squares
C                                of the residuals of the best fit in
C                                SS(N).
      SS(N) = G(KSECTN,N2)
      MODE = 1
C
  390 CONTINUE
C                               Reverse sign of F(.), Y(.) and PAR if
C                                first monotonic section was required
C                                to be decreasing.
      IF (IORDER.EQ.-1) THEN
          DO 400 I = I1,N
             F(I) = -F(I)
             Y(I) = -Y(I)
             PAR(I) = -PAR(I)
  400     CONTINUE
      END IF
C
C
      CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,1000)
C                               MSG: RETURN FROM L2WPMA WITH MODE
C                                    NUMBER OF DATA (F)
C                                    NUMBER OF MONOTONIC SECTIONS
C                                    NUMBER OF LOCAL MINIMA IN F
C                                    NUMBER OF LOCAL MAXIMA IN F
C                                    DEGREES OF FREEDOM OF BEST FIT
C.... End of Subroutine L2WPMA ........................................
C
  410 RETURN
      END
C
C......................................................................
C     This subroutine calculates
C      1) the best weighted least squares monotonic increasing
C      approximation to univariate data indexed from L1 to LN,
C      where I1 <= L1 <= LN <=N,
C      2) the associated Lagrange multipliers, and
C      3) the spline order 1 representation of the best approxtn.
C......................................................................
      SUBROUTINE L2WMON(SS,Y,WY,F,WF,Z,WZ,IW,IAKNW,I1,N,L1,LN,MODE,NK,
     +                  IAKN,NACT,IACT,PAR)
C
C
C.... I N P U T ....
C     I1       Integer variable, lower data index.
C     N        Integer variable, upper data index, N >= I1.
C     L1       Integer variable, that has to be I1 <= L1.
C     LN       Integer variable, that has to be LN <= N.
C     F(I1:N)  Real array of function measurements, F(i), i=I1,...,N.
C     WF(I1:N) Real array of weights WF(i), i=I1,...,N, associated with
C               function measurements F(i), i=I1,...,N.
C     MODE     Integer variable that takes two values.
C              MODE = 0 neither the values {Y(i): i=L1,...,LN} nor
C                       their associated weights are calculated.
C              MODE = 1 three actions are taken:
C                       1) the values {Y(i): i=L1,...,LN} and their
C                       weights {WY(i): i=L1,...,LN} are calculated.
C                       2) array IAKN of knot indices of spline
C                       representation of fit is formed, and
C                       3) array IACT of indices of active constraints
C                       and, corresponding Lagrange multipliers
C                       (array PAR) are calculated.
C     IAKN(I1-1) Integer variable that takes two values: =-1, in order
C              to indicate that L2WMON calculates the decreasing fit,
C              or =1, in order to indicate that L2WMON calculates the
C              increasing fit to the data.
C
C.... O U T P U T ....
C     SS(I1:N) Real array containing in SS(i), i=L1,...,LN the weighted
C               sum of the squares of the residuals by the best monotnc
C               increasing approximation to the data F(L1),...,F(i),
C               L1 <= i <= LN.
C     Y(I1:N)  Real array that, if MODE=1, contains in Y(i), i=L1,...,
C               LN the best monotonic increasing approximation to F(i),
C               i=L1,...,LN, ie Y(.) satisfies the inequalities
C               Y(L1) <=...<= Y(LN).
C     WY(I1:N) Real array that, if MODE=1, contains in WY(i),
C               i=L1,...,LN the weights of Y(i), i=L1,...,LN,
C               respectively.
C     NK       Number of indices of knots.
C     IAKN(I1-1:N) Defined in subroutine L2WPMA.
C     NACT     Number of active constraints.
C     IACT(1:NACT) Defined in subroutine L2WPMA.
C     PAR(I1+1:N) Defined in subroutine L2WPMA.
C
C.... W O R K I N G  S P A C E ....
C     Z(1:N+I1-1) Real array, stores the distinct components of the
C               best monotonic approximation Y(L1:LN). The actual space
C               used by Z at the end of the calculation is Z(1:IS),
C               where IS is a positive integer that is explained below.
C     WZ(1:N+I1-1) Real array, stores the weights of the components
C               of Z(.).
C     IW(1:N+I1-1) Integer array such that Z(i) is repeated IW(i) times
C               in the best monotonic approximation Y(L1:LN). At the
C               end of the calculation IW(1)+IW(2)+...+IW(IS)=LN-L1+1.
C     IAKNW(1:N+I1-1) Integer array, stores temporarily the indices of
C               the knots of the spline representation of Y(.).
C
C.... M E T H O D ....
C     Implements a modification of Algorithm 1 of Demetriou & Powell,
C      IMA J. Numer. Anal., 11(1991), which takes account of positive
C      weights corresponding to the data. INFO: A generalization of
C      this algorithm for a strictly convex objective function is:
C      Algorithm 3, Demetriou, Math. of Computation, 64(1995).
C
C     The best monotonic increasing approximation Y(L1:LN) is
C      presented by the triple (IS,IW,Z), as follows: IW and Z have
C      each IS components. Y has IW(1) components equal to Z(1),
C      then IW(2) components equal to IW(2), and so on up to IW(IS)
C      components equal to IW(IS).
C     The spline representation of Y(L1:LN) is due to defining the
C      knot indices in IAKN(1:NK), where NK is the number of knots and
C      IAKN(NK)=N. We define a knot index to be at the leftist
C      index of equal components of Y(.). Since L2WMON is used for
C      calculating both the increasing and the decreasing case,
C      IAKN is defined as follows (we assume that I1=1):
C               If IAKN(0) .eq. 1, then IAKN(1)=L1,
C                                       IAKN(2)=IAKN(1)+IW(1)
C                                       IAKN(3)=IAKN(2)+IW(2), etc.
C      However, If IAKN(0) .eq. -1, then IAKN(1)=L1+IW(1)-1
C                                        IAKN(2)=IAKN(1)+IW(2)
C                                        IAKN(3)=IAKN(2)+IW(3), etc.
C      Therefore, the statements that include label 80 below arrange
C      for the knot indices.
C
C     IACT(1:NACT) keeps the indices of the constraints satisfied as
C      equalities by the best fit. Specifically, if j is an integer in
C      [I1+1,N] such that IACT(i)=j, for some i in [1,NACT], then
C      Y(j-1)=Y(j).
C
C.... W A R N I N G ....
C     This subrtne will fail if LN < L1 or L1 < I1 or N < LN. However,
C      whenever it is called by subroutine L2WPMA the conditions
C      I1 <= L1 <= LN <= N  are satisfied by default.
C.....................................................................
C     Initialization.
C
C     .. Scalar Arguments ..
      INTEGER I1,L1,LN,MODE,N,NACT,NK
C     ..
C     .. Array Arguments ..
      REAL SS(I1:N),Y(I1:N),WY(I1:N),F(I1:N),WF(I1:N),Z(1:N + I1 - 1),
     +     WZ(1:N + I1 - 1),PAR(I1:N)
      INTEGER IW(1:N + I1 - 1),IAKNW(1:N + I1 - 1),IAKN(I1-1:N-I1+1),
     +     IACT(1:N-I1)
C     ..
C     .. Local Scalars ..
      INTEGER I,IS,ITEMP,J,JJ,JNK,JSWI,NACT0
C     ..
      J = L1
      IS = 1
      Z(IS) = F(J)
      WZ(IS) = WF(J)
      IW(IS) = 1
      SS(J) = 0.
C
C     Termination.
C
   10 IF (J.EQ.LN) GO TO 50
C
      J = J + 1
      IS = IS + 1
      Z(IS) = F(J)
      WZ(IS) = WF(J)
      IW(IS) = 1
      SS(J) = SS(J-1)
C
C     If (IS.eq.1 .or. Z(IS - 1).le.Z(IS)) then go to 10
   30 IF (IS.EQ.1) GO TO 10
      IF (Z(IS - 1).LE.Z(IS)) GO TO 10
C
      SS(J) = SS(J) + (Z(IS) - Z(IS - 1))**2*WZ(IS)*WZ(IS - 1)/
     +        (WZ(IS) + WZ(IS - 1))
C
      IS = IS - 1
      Z(IS) = (WZ(IS)*Z(IS) + WZ(IS + 1)*Z(IS + 1))/(WZ(IS) + WZ(IS+1))
      WZ(IS) = WZ(IS) + WZ(IS + 1)
      IW(IS) = IW(IS) + IW(IS + 1)
      GO TO 30
C
C     Calculate the components and the weights of the best monotonic
C      increasing approximation. Also, define IAKN(.) and find NK.
C
   50 CONTINUE
      IF (MODE.EQ.0) GO TO 100
C
      NACT0 = NACT
C
      J = L1
      DO 70 I = 1,IS
C                               Switch for Lagrange parameters.
         JSWI = 0
C
         JJ = J + IW(I) - 1
C                               Calculate indices of knots.
         NK = NK + 1
         IF (IAKN(I1-1) .EQ. -1) THEN
C                               If IAKN(I1-1) = -1 at calling sequence
C                                of L2WMON, the approximation is
C                                decreasing. Hence knot indices are
C                                defined at righthand-end of equal
C                                components.
             IAKNW(I) = JJ
         ELSE
C                               IAKN(I1-1) = 1 at calling sequence
C                                of L2WMON, so the approximation is
C                                increasing. Hence knot indices are
C                                defined at lefthand-end of equal
C                                components.
             IAKN(NK) = J
         END IF
C
C        Calculate components, weights and Lagrange parameters
C         of the best approximation.
C
   60    Y(J) = Z(I)
         WY(J) = WZ(I)
         J = J + 1
         IF (J.LE.JJ) THEN
C	   		       			    	            
C            Calculate Lagrange parameter. Update IACT.
             NACT = NACT + 1
C		                               
             IF (IAKN(I1-1) .EQ. -1) THEN
C                Decreasing approximation (negative sign). Calculate
C                 backwards.
                 IF (JSWI.EQ.0) THEN
                     PAR(N+I1+1-J) = 2*WF(J-1)*(Z(I)-F(J-1))
C                               Switch for next Lagrange par.
                     JSWI = 1
                 ELSE
                     PAR(N+I1+1-J) = PAR(N+I1+1-J+1) +
     +                                        2*WF(J-1)*(Z(I)-F(J-1))
                 END IF
                 IACT(NACT) = N+I1+1-J
C
             ELSE
C                Increasing approximation.
                 IF (JSWI.EQ.0) THEN
                     PAR(J) = -2*WF(J-1)*(Z(I)-F(J-1))
C                               Switch for next Lagrange par.
                     JSWI = 1
                 ELSE
                     PAR(J) = PAR(J-1) - 2*WF(J-1)*(Z(I)-F(J-1))
                 END IF
                 IACT(NACT) = J
C
             END IF
C
             GO TO 60
         END IF
C
         J = JJ + 1
   70 CONTINUE
C
      IF (IAKN(I1-1) .EQ. -1) THEN
C                               Arrange knot indices in [L1,LN],
C                                where L1 and LN are defined by
C                                subroutine L2WPMA just before
C                                calling L2WMON.
          IAKN(NK) = N + I1 - IAKNW(1)
          JNK = NK
          IF (IS .EQ. 1) GOTO 85
          DO 80 I = 2,IS
             JNK = JNK -1
             IAKN(JNK) = N + I1 - IAKNW(I)
   80     CONTINUE
C
   85     CONTINUE
C                               Reverse IACT in [L1,LN].
          IF (NACT .GT. NACT0) THEN
              DO 90 I = 1,(NACT-NACT0)/2
                 ITEMP = IACT(NACT0+I)
                 IACT(NACT0+I) = IACT(NACT-I+1)
                 IACT(NACT-I+1) = ITEMP
   90         CONTINUE
          END IF
C
      END IF
C	         
C.... End of Subroutine L2WMON ........................................
C
  100 RETURN
      END
C
C......................................................................
C     This subroutine contains messages associated with the operations
C      of subroutines L2WPMA and TRIVIA.
C......................................................................
      SUBROUTINE MESSGW(MODE,I1,N,KSECTN,NL,NU,LABEL)
C
C.... I N P U T ....
C     MODE  Integer variable, indicating the termination status of
C            subroutine L2WPMA.
C     I1,N,KSECTN Defined in subroutine L2WPMA.
C     NL    Integer variable defined in subroutine L2WPMA. The value
C            NL = N - I1 + 1 indicates termination of L2WPMA due to
C            TRIVIA. Otherwise, I1 <= NL <= (N - I1 + 1)/2 (due to
C            subroutine XTREMA).
C     NU    Integer variable defined in subroutine L2WPMA. The value
C            NU = 0 indicates termination of L2WPMA due to TRIVIA.
C            Otherwise, I1 <= NU <= (N - I1 + 1)/2 (due to XTREMA).
C     LABEL Integer variable, indicating the level in increasing order
C            in subroutines TRIVIA and L2WPMA, where the message is
C            caused. The LABELs of TRIVIA are numbered up to 80,
C            while the LABELs of L2WPMA are numbered from 100.
C
C.... O U T P U T ....
C     Various messages that are printed when subroutine L2WPMA returns.
C......................................................................
C     Messages from subroutine TRIVIA.
C
C     .. Scalar Arguments ..
      INTEGER I1,KSECTN,LABEL,MODE,N,NL,NU
C     ..
      IF (LABEL.EQ.10) THEN
          WRITE (*,FMT=9000)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.15) THEN
          WRITE (*,FMT=9010) KSECTN,N - I1 + 1
          GO TO 10
      END IF
C
      IF (LABEL.EQ.20) THEN
          WRITE (*,FMT=9020)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.60) THEN
          WRITE (*,FMT=9030)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.70) THEN
          WRITE (*,FMT=9040)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.80) THEN
          WRITE (*,FMT=9050)
          GO TO 10
      END IF
C
C......................................................................
C     Messages from subroutine L2WPMA.
C
      IF (LABEL.EQ.100) THEN
          WRITE (*,FMT=9055)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.110) THEN
          WRITE (*,FMT=9058)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.120) THEN
          WRITE (*,FMT=9060)
          GO TO 10
      END IF
C
      IF (LABEL.EQ.130) THEN
          WRITE (*,FMT=9070) KSECTN
          GO TO 10
      END IF
C
   10 CONTINUE
C
      IF (LABEL.EQ.1000) THEN
          IF (NL.EQ.N - I1 + 1 .AND. NU.EQ.0) THEN
C                               NL = N - I1 + 1 and NU = 0 are used to
C                                indicate termination due to the cases
C                                checked by subroutine TRIVIA.
              WRITE (*,FMT=9080) MODE,N - I1 + 1,KSECTN
      ELSE
              WRITE (*,FMT=9090) MODE,N - I1 + 1,KSECTN,NL - I1 + 1,
     +                                                     NU - I1 + 1
          END IF
      END IF
C
      RETURN
C
 9000 FORMAT (//5X,'Calculation incomplete because required ',/5X,' nu',
     +       'mber of monotonic sections is less than one.',/5X,' It i',
     +       's recommended:',' 1<= Number Monotonic Sections < Number',
     +       ' of data')
 9010 FORMAT (//5X,'WARNING: ',/5X,' Required number of monotonic sect',
     +       'ions   =',I5,/5X,' is greater than/equal to number of da',
     +       'ta =',I5,/5X,' it is recommended:',' 1<= Number Monotoni',
     +       'c Sections < Number of data')
 9020 FORMAT (//5X,'There is only one data point, so data itself',/5X,
     +       ' is the required approximation.',/5X,' Now number of mon',
     +       'otonic sections is exactly 0.')
 9030 FORMAT (//5X,'Data itself is the required approximation.',/5X,
     +       ' (Best monotonic increasing)',/5X,' Now number of monoto',
     +       'nic sections is exactly 1.')
 9040 FORMAT (//5X,'Best monotonic approximation:',/5X,' Required ',
     +       'number of monotonic sections is',' exactly 1.')
 9050 FORMAT (//5X,'Data itself is the requires approximation.',/5X,
     +       ' (Best monotonic decreasing)',/5X,' Now number of monoto',
     +       'nic sections is exactly 2',/5X,' but the first section d',
     +       'egenerates to a single point.')
 9055 FORMAT (//5X,'(First monotonic section is required to be decreas',
     +       'ing.)')
 9058 FORMAT (//5X,'Automatic weight generation due to abscissae spaci',
     +       'ng')
 9060 FORMAT (//5X,'(Termination due to subroutine TRIVIA.)')
 9070 FORMAT (//5X,'Data itself is the best approximation, ',/5X,' bec',
     +       'ause number or required monotonic sectns=',I5,/5X,' is e',
     +       'ither equal to, or it exceeds number of,',/5X,
     +       ' monotonic',' sections in data.')
 9080 FORMAT (//5X,'----------------------------------------',/5X,'RET',
     +       'URN FROM L2WPMA (TRIVIA) WITH MODE = ',I5,/5X,'NUMBER OF',
     +       ' DATA                        = ',I5,/5X,'NUMBER OF MONOT',
     +       'ONIC SECTIONS          = ',I5,/5X,'---------------------',
     +       '-------------------')
 9090 FORMAT (//5X,'-----------------------------------------',/5X,'RE',
     +       'TURN FROM L2WPMA WITH MODE = ',I5,/5X,'NUMBER OF DATA (F',
     +       ')           = ',I5,/5X,'NUMBER OF MONOTONIC SECTIONS = ',
     +       I5,/5X,'NUMBER OF LOCAL MINIMA IN F  = ',I5,/5X,'NUMBER O',
     +       'F LOCAL MAXIMA IN F  = ',I5,/5X,'-----------------------',
     +       '------------------')
C
C.... End of Subroutine MESSGW ........................................
C
      END
C
C......................................................................
C     This subroutine investigates for trivial cases, when a piecewise
C      monotonic approximation with KSECTN sections to the data F(I1:N)
C      is required.
C......................................................................
      SUBROUTINE TRIVIA(I1,N,F,KSECTN,Y,ITHETA,MODE,NL,NU)
C
C
C.... I N P U T ....
C     I1,N,F,KSECTN     Defined in subroutine L2WPMA.
C
C.... O U T P U T ....
C     KSECTN    This is an input variable (user defined) that preserves
C                its value throughout the calculation (all subroutines)
C                except of some trivial cases defined at the
C                paragraph below under the heading "Method", where its
C                value changes to:
C                = 0
C                = 1
C                = 2
C                Note that immediately after this change, the
C                 calculation stops.
C     Y(I1:N)   Real array, that holds the required piecewise monotonic
C                approximation to F(I1:N).
C     ITHETA(0:KSECTN) Integer array that holds the indices of the
C                extrema of the optimal approximation derived due to
C                trivial cases. See subroutine L2WPMA for definition
C                of ITHETA.
C     MODE      Integer variable, to be used by subroutine L2WPMA. At
C                exit from TRIVIA, takes one of the following values:
C                = -1      then the best monotonic increasing apprxtn
C                           is calculated on return to subroutine
C                           L2WPMA.
C                =  0      indicating unsuccessful termination of
C                           L2WPMA.
C                =  2      indicating termination of L2WPMA because the
C                           data itself is the best approximation.
C                =  3      on return to subrtn L2WPMA, the calculation
C                           proceeds ahead.
C
C.... W O R K I N G  S P A C E ....
C     NL,NU     Defined in subroutine XTREMA.
C
C.... M E T H O D ....
C     Some trivial cases are examined sequentially:
C      1) If KSECTN < 1, TRIVIA sets MODE = 0 and exits (then L2WPMA
C          terminates unsuccessfully).
C      2) If KSECTN > N-I1, a warning message appears.
C      3) If N = I1 and KSECT >= 1, TRIVIA sets KSECTN = 0, MODE = 2
C          and exits (data is the required fit).
C      4) If F(i) <= F(i+1), i=I1,...,N-1 and KSECTN >= 1, TRIVIA sets
C          KSECTN = 1, MODE = 1 and exits (data is the required fit).
C      5) If KSECTN = 1, TRIVIA sets MODE = -1 and exits. Then subrtne
C          L2WMON calculates the best monotonic (increasing) fit to
C          the data. If IORDER = -1 at subroutine L2WPMA, then at the
C          end of the calculation the best monotonic decreasing fit to
C          the data is obtained.
C      6) If F(i) >= F(i+1), i=I1,...,N-1 and KSECTN > 1, TRIVIA sets
C          KSECTN = 2, MODE = 1 and exits (data is the required fit).
C......................................................................
C     .. Scalar Arguments ..
      INTEGER I1,KSECTN,MODE,N,NL,NU
C     ..
C     .. Array Arguments ..
      REAL F(I1:N),Y(I1:N)
      INTEGER ITHETA(0:KSECTN)
C     ..
C     .. Local Scalars ..
      INTEGER I,IA
C     ..
C     .. External Subroutines ..
      EXTERNAL MESSGW
C     ..
      IF (KSECTN.LT.1) THEN
C                             MSG: Calculation incomplete, because
C                                   required number of monotonic sectns
C                                   is less than one. It is recommended
C                                   1 <= KSECTN <= N-I1.
          MODE = 0
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,10)
          GO TO 70
      END IF
C
      IF (KSECTN.GT.N - I1) THEN
C                             MSG: WARNING: Required number of
C                                   monotonic sections is greater than
C                                   or equal to number of data. It is
C                                   recommended 1 <= KSECTN <= N-I1.
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,15)
      END IF
C
      IF (N.EQ.I1 .AND. KSECTN.GE.1) THEN
          KSECTN = 0
          Y(I1) = F(I1)
          ITHETA(0) = I1
          ITHETA(0) = N
C                             MSG: There is only one data point. So
C                                   data itself is the required apprxn,
C                                   and number of monotonic sections
C                                   is equal to 0.
          MODE = 2
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,20)
          GO TO 70
      END IF
C
C                             If F(I1) <=...<= F(N) then set Y(.)=F(.)
C                              At this stage it holds KSECTN >= 1.
      IA = I1
   30 IF (F(IA).LE.F(IA + 1)) THEN
          IA = IA + 1
          IF (IA.LT.N) GO TO 30
      END IF
C
      IF (IA.EQ.N) THEN
          KSECTN = 1
          DO 40 I = I1,N
             Y(I) = F(I)
   40     CONTINUE
          ITHETA(0) = I1
          ITHETA(KSECTN) = N
C                             MSG: Data itself is the required appxtn.
C                                   (Best monotonic increasing) and
C                                   number of monotonic sectns is
C                                   equal to 1.
          MODE = 2
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,60)
          GO TO 70
      END IF
C
      IF (KSECTN.EQ.1) THEN
C                             On return to subroutine L2WPMA, subrtn
C                              L2WMON calculates the best monotonic
C                              apprxmtion, increasing or decreasing,
C                              depending on the value of the argument
C                              IORDER.
C                             MSG: Best monotonic apprxtn.
C                                   Required number of monotonic sectns
C                                   is exactly 1.
          MODE = -1
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,70)
          GO TO 70
      END IF
C
C                             If F(I1) >=...>= F(N) then set Y(.)=F(.)
C                              At this stage it holds KSECTN > 1.
      IA = I1
   50 IF (F(IA).GE.F(IA + 1)) THEN
          IA = IA + 1
          IF (IA.LT.N) GO TO 50
      END IF
C
      IF (IA.EQ.N) THEN
          KSECTN = 2
          DO 60 I = I1,N
             Y(I) = F(I)
   60     CONTINUE
          ITHETA(0) = I1
          ITHETA(1) = I1
          ITHETA(KSECTN) = N
C                             MSG: Data itself is the required appxtn.
C                                   (Best monotonic decreasing)
C                                   Now number of monotonic sectns is
C                                   exactly to 2, but the first section
C                                   degenerates to a point.
C
          MODE = 2
          CALL MESSGW(MODE,I1,N,KSECTN,NL,NU,80)
          GO TO 70
C
      END IF
C......................................................................
C     At this stage of the calculation, the inequality  1 < KSECTN is
C      satisfied. Below we set MODE=3, in order to be able to activate
C      the main calculation of subroutine L2WPMA. Further, we note that
C      the case when KSECTN=2 is treated directly by subrtn L2WPMA.
C
      MODE = 3
C
C.... End of Subroutine TRIVIA ........................................
C
   70 RETURN
      END
C
C......................................................................
C     This subroutine forms the sets of indices of LOWER and UPPER
C      values of the data, ie the sets of the indices of the local
C      minima and the local maxima respectively of the data.
C......................................................................
C
      SUBROUTINE XTREMA(F,LOWER,IUPPER,INDX,KSECTN,I1,N,N2,NL,NU)
C
C
C.... I N P U T ....
C     F,KSECTN,I1,N,N2     As they are defined in Subroutine L2WPMA.
C
C.... O U T P U T ....
C     LOWER(I1:N2)  Integer array, indices of lower values of F(.).
C     IUPPER(I1:N2) Integer array, indices of upper values of F(.).
C     INDX(I1:N)    Integer array such that INDX(LOWER(I)) = I and
C                    INDX(IUPPER(I)) = I.
C     NL            Integer variable, such that (NL - I1 +1) is the
C                    number of elements of LOWER(.).
C     NU            Integer variable, such that (NU - I1 + 1) is the
C                    number of elements of IUPPER(.).
C
C.... M E T H O D ...
C     Described by Demetriou & Powell, IMA J. Numer. Anal., 11(1991).
C.....................................................................
C
C     .. Scalar Arguments ..
      INTEGER KSECTN,I1,N,N2,NL,NU
C     ..
C     .. Array Arguments ..
      REAL F(I1:N)
      INTEGER INDX(I1:N),IUPPER(I1:N2),LOWER(I1:N2)
C     ..
C     .. Local Scalars ..
      INTEGER I,IA,IB,II,IL,IU
C
C     IA   Such that F(I1) = F(I1 + 1) =...= F(IA), where I1 < IA < N.
C     IB   Such that F(IB) = F(IB + 1) =...= F(N), where IA < IB <N.
C     IL   Counter on LOWER(.).
C     IU   Counter on IUPPER(.).
C     ..
C                              Terminate, if KSECTN < 1.
      IF (KSECTN.LT.1) GO TO 60
C
C.... Form sets of indices of LOWER and UPPER values of data...........
C
      IL = I1
      LOWER(IL) = I1
      INDX(I1) = IL
      IU = I1 - 1
C
      IA = I1
   10 IF (F(IA).EQ.F(IA + 1)) THEN
          IA = IA + 1
          IF (IA.LT.N) GO TO 10
      END IF
C..............................Trivial case F(I1) = F(I1+1) =...= F(N).
C                               Thus IA = N (and, see below, IB = I1).
C
C                              Include I1 to IUPPER.
      IF (IA.EQ.N) THEN
          IU = I1
          IUPPER(IU) = I1
          INDX(I1) = IU
C                              Include N either to LOWER or to IUPPER.
          IF (KSECTN.EQ.KSECTN/2*2) THEN
              IL = IL + 1
              LOWER(IL) = N
              INDX(N) = IL
C
          ELSE
              IU = IU + 1
              IUPPER(IU) = N
              INDX(N) = IU
          END IF
C
          GO TO 50
C
      END IF
C.............................. Extrema inside the range of the data.
C
C                               Include I1 to IUPPER.
      IF (F(IA).GT.F(IA + 1)) THEN
          IU = I1
          IUPPER(IU) = I1
          INDX(I1) = IU
      END IF
C                               I1 <= IA < IB <= N.
      IB = N
   20 IF (F(IB).EQ.F(IB - 1)) THEN
          IB = IB - 1
          GO TO 20
      END IF
C
      DO 40 I = IA + 1,IB - 1
         II = I
   30    IF (F(II).EQ.F(II + 1) .AND. II.LT.IB - 1) THEN
             II = II + 1
             GO TO 30
         END IF
C                              Include I in LOWER.
         IF (F(I - 1).GT.F(I) .AND. F(II).LT.F(II + 1)) THEN
             IL = IL + 1
             LOWER(IL) = I
             INDX(I) = IL
         END IF
C                              Include I in IUPPER.
         IF (F(I - 1).LT.F(I) .AND. F(II).GT.F(II + 1)) THEN
             IU = IU + 1
             IUPPER(IU) = I
             INDX(I) = IU
         END IF
C
   40 CONTINUE
C                              Include N in LOWER.
      IF (KSECTN.EQ.KSECTN/2*2 .OR. F(IB - 1).GT.F(IB)) THEN
          IL = IL + 1
          LOWER(IL) = N
          INDX(N) = IL
      END IF
C                              Include N in IUPPER.
      IF (KSECTN.NE.KSECTN/2*2 .OR. F(IB - 1).LT.F(IB)) THEN
          IU = IU + 1
          IUPPER(IU) = N
          INDX(N) = IU
      END IF
C
   50 CONTINUE
      NL = IL
      NU = IU
   60 CONTINUE
C
C.... End of Subroutine XTREMA .......................................
C
      RETURN
      END
