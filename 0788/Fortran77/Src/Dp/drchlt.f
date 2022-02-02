
      SUBROUTINE DRCHLT(IOUT,IDBG,IBEG,IE,BDYFCN,CURVE,B,DPTS,NPTS,
     +                  EPS0,R_FORM,WORK,IWORK,NWORK,MWORK,UVEC,
     +                  ERROR,IER)
C     -----------------
C
C     This program will solve the interior and exterior Dirichlet
C     problems for Laplace's equation on a simply-connected region D
C     with smooth boundary C.  It is assumed that C is at least two
C     times continuously differentiable. To define a boundary, the user
C     must supply a parametric form of the external routine CURVE with 
C     parameter range [0,B]. The routine CURVE must also produce the 
C     first and second derivatives on the curve as defined below. The 
C     boundary condition must be given by the routine BDYFCN.
C
C     INPUT PARAMETERS:  
C
C     IOUT      Input.
C               The output unit number.
C     IDBG      Input.
C               The debugging parameter.
C               =0 produces a shortcut output.
C               =1 produces a debugging output. 
C     IBEG      Input.
C               =0 means this is a first call on DRCHLT for this
C               particular curve C, boundary function BDYFCN, and
C               error tolerance EPS.
C               =1 means that DRCHLT has been called previously for
C               this choice of parameters, and the solution U is
C               desired at a new set of points in DPTS. For such a
C               call on DRCHLT, change only DPTS, NPTS, AND IBEG. Do
C               not change any other input variable.
C     IE        Input.
C               =0 for the interior Dirichlet problem.
C               =1 for the exterior Dirichlet problem.
C     BDYFCN    External function.
C               Inputs : x,y
C               Outputs: BDYFCN
C               This is a function of two variables X,Y.  It gives
C               the value of the harmonic function U on the
C               boundary curve C. 
C     CURVE     External subroutine. 
C               Inputs : s
C               Outputs: X,Y,DX,DY,D2X,D2Y
C               This is a subroutine which defines the boundary
C               curve C.  The calling sequence for it is
C                   CALL CURVE(S,X,Y,DX,DY,D2X,D2Y)
C               The variable S is the parameterization variable for
C               the curve C. (X,Y) is the point on the curve
C               corresponding to the variable S. DX,DY,D2X,D2Y are
C               the first and second derivatives of X(S) and Y(S),
C               respectively. 
C     B         Input.
C               The limits of the variable S in defining the curve C
C               are  0 .LE. S .LE. B
C     DPTS      User supplied array.
C               This is a two dimensional array specifying the points
C               at which the harmonic function U is to be evaluated.
C               The dimension statement for DPTS is
C                      DIMENSION DPTS(2,NPTS)
C               The point #J is given by
C                   X(J)=DPTS(1,J),  Y(J)=DPTS(2,J)
C     NPTS      Input.
C               This is the number of points in DPTS.  It must be 
C               greater than zero.
C     EPS0      Input.
C               The desired absolute error tolerance in the solution U.
C     R_FORM    Input.
C               This specifies the way in which the variable RATE
C               is to be defined in the routines INTEQN and EVALU.
C               =0 means  we use the "normal" way to define RATE
C               based on estimating the rate of convergence in the
C               approximates calculated to date.  
C               =1 means we use a "conservative" error test
C               in which RATE=0.5 and the approximates are assumed to
C               have a very slow rate of convergence.  Use this for 
C               more ill-behaved problems and boundaries.
C     WORK      Real work array
C               Temporary work space. This vector should have a
C               dimension NWORK of at least 5,000. If some points of
C               DPTS are close to the boundary, then the dimension of
C               WORK should be increased further. A dimension of
C               NWORK=300,000 will allow a great many problems to be
C               treated very accurately.  For more details on
C               computing the needed dimension NWORK for WORK, see
C               the discussion in the driver program for DRCHLT or
C               the accompanying paper.
C     IWORK     Integer work array
C               Temporary integer work space. It will be used as
C               pivot and work array in LAPACK.
C     NWORK     Input.
C               Dimension of WORK. NWORK = 300,000 will solve  a lot of
C               problems very accurately.  See the preceding discussion
C               of WORK and the general discussion of NWORK given below.
C     MWORK     Input.
C               Dimension of IWORK. It is dependendent of NWORK. 
C                    MWORK .GE. (SQRT(NWORK+36)-6).
C     UVEC      Output vector. 
C               On exit, it will contain the approximate value of the 
C               solution function U(X(I),Y(I)) at the points (X(I),Y(I))
C               given in DPTS.
C     ERROR     Output.
C               This is an output vector containing predicted error
C               bounds for the solutions given in the vector UVEC.
C               The element UVEC(J) is the approximate solution at
C               the # J point of DPTS, and ERROR(J) is a predicted
C               error bound. If the desired error EPS was not
C               attained, then the sign of ERROR(J) is made negative
C               for such a solution value. Otherwise the sign of
C               ERROR(J) is positive.
C     IER       Output.
C               = 1 means some or all of the solution values in UVEC
C               may not satisfy the error tolerance EPS.
C               = 0 means the program was completed satisfactorily.
C               =-1 means EPSO < 0.
C               =-2 means NWORK < 0 or MWORK < SQRT(NWORK*36) - 6.
C               =-3 means B  .LE. 0.
C               =-4 means IE or IBEG are out of range. 
C               =-5 means NPTS<=0, an error.
C
C     *** Defining NWORK, the dimension of WORK ***
C     Introduce variables MAXSYS and MAXMES.  MAXSYS represents
C     the maximum order of the linear system to be solved in
C     subroutine INTEQN; and MAXMES represents the maximum number
C     of mesh points of C at which the double layer density
C     function is to be evaluated.  MAXFFT represents the maximum
C     degree of the Fourier expansion to be produced for the 
C     approximate density function.  As examples,
C            MAXSYS=512, MAXMES=8192
C     will solve many problems very accurately.  We define NWORK by
C            NWORK=MAX(MAXSYS**2+12*MAXSYS,7*MAXMES)
C     The program assumes 
C           MAXSYS .GE. 64,  MAXMES .GE. 128
C     and you should set NWORK accordingly.  These defaults can be 
C     changed by re-setting the respective parameters N_0 and M_0 
C     in the routines EVALU and INTEQN, respectively.  As can be
C     noted from the numbers given, these parameters should be 
C     chosen as powers of 2.
C
C     *** SOURCES OF INFORMATION ***
C
C     For additional information on this program, see the paper
C
C        K. Atkinson & Y. Jeon, "Automatic boundary integral equation 
C        programs for the planar Laplace equation", ACM Transactions on 
C        Mathematical Software. 
C
C     The authors can be contacted at the respective e_mail addresses:
C        Kendall-Atkinson@uiowa.edu
C        yjeon@madang.ajou.ac.kr
C
C     The web site for this program is located at the URL
C        http://www.math.uiowa.edu/~atkinson/laplace.html       
C
C     DATE OF LAST CHANGES TO THIS CODE: 7 April 1998
C
C     ================= END OF INTRODUCTORY REMARKS ===================
C
      INTEGER IBEG,IDBG,IE,IER,IOUT,MWORK,NPTS,NWORK,R_FORM
      DOUBLE PRECISION DPTS(2,NPTS),ERROR(NPTS),UVEC(NPTS),WORK(NWORK)
      INTEGER IWORK(MWORK)
      DOUBLE PRECISION B,BDYFCN,EPS0
      EXTERNAL BDYFCN,CURVE
C    
C     LOCAL VARIABLES
      DOUBLE PRECISION D1MACH,EP,EPS,FIW,FL,FN,RHOERR,U100
      DOUBLE PRECISION ZERO,TWO,SIX,SEVEN
      INTEGER I,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10,IEE,IER1,IER2,II,
     +        J,K,LB,LD1,LD2,LD3,LD4,LD5,LD6,LD7,LU,NB,NFINAL,NU
      INTEGER LBASE(7),NBASE(10)
      EXTERNAL D1MACH,EVALU,INTEQN
      INTRINSIC FLOAT,LOG,MAX0,SIGN,SQRT
      DATA ZERO/0.0D0/, TWO/2.0D0/, SIX/6.0D0/, SEVEN/7.0D0/
C      
C     TEST THE INPUT PARAMETERS.
C
      IF (EPS0 .LE. 0) THEN
          IER = -1
          RETURN
      END IF    
      IF (NWORK .LE. 0  .OR. 
     +       MWORK .LE. (SQRT(NWORK+SIX*SIX) - SIX)) THEN
          IER = -2
          RETURN
      END IF 
      IF (B .LE. ZERO) THEN
          IER = -3
          RETURN
      END IF
      IF ( (IE .LT. 0 .OR. IE .GT. 1) .OR. 
     +     (IBEG .LT. 0 .OR. IBEG .GT. 1) )  THEN
          IER = -4
          RETURN
      END IF
      IF( NPTS .LE. 0) THEN
          IER = -5
          RETURN
      END IF
C
C     SET MACHINE DEPENDENT CONSTANT. U100 IS 100 TIMES
C     MACHINE UNIT ROUND.
      U100 = 100*D1MACH(4)
      EPS = EPS0
      IF (IBEG .EQ. 1) GO TO 10
      EP = EPS/TWO
C     OBTAIN VALUE OF NUPPER FOR USE IN INTEQN
      FIW = NWORK
      FN = SQRT(FIW+SIX*SIX) - SIX
      IEE = LOG(FN)/LOG(TWO) 
      NU = 2**IEE
C     BREAK UP WORK INTO VECTORS AND MATRIX FOR USE IN INTEQ.
C     PRODUCE RELATIVE ADDRESSES WITHIN  WORK OF  RHO,X,...,D2Y,
C     OLDRHO,Z.
      NBASE(1) = 1
      NBASE(2) = 1 + NU*NU
      DO I = 3,10
          NBASE(I) = NBASE(I-1) + NU
      END DO 
      ID2 = NBASE(2)
      ID3 = NBASE(3)
      ID4 = NBASE(4)
      ID5 = NBASE(5)
      ID6 = NBASE(6)
      ID7 = NBASE(7)
      ID8 = NBASE(8)
      ID9 = NBASE(9)
      ID10 = NBASE(10)
C
      CALL INTEQN(IOUT,IDBG,IE,B,EP,R_FORM,BDYFCN,CURVE,NU,WORK(ID2),
     +            RHOERR,NFINAL,WORK(ID3),WORK(ID4),WORK(ID5),
     +            WORK(ID6),WORK(ID7),WORK(ID8),WORK(ID9),IWORK,
     +            WORK(ID10),WORK(1),IER1)
      IF (IDBG .EQ. 1)  THEN
          WRITE (IOUT,FMT=9000) NFINAL,RHOERR,IER1
          WRITE (IOUT,FMT=9010)
          DO I = 1,NFINAL
              WRITE (IOUT,FMT=9020) WORK(ID3+I-1),WORK(ID4+I-1),
     +               WORK(ID2+I-1)
          END DO
      END IF      
C
      IF (IER1 .EQ. 1) EP = RHOERR
C     OBTAIN LUPPER FOR USE IN EVALU.
      FL = FIW/SEVEN
      IEE = LOG(FL)/LOG(TWO) + U100
      LU = 2**IEE
C     OBTAIN RELATIVE ADDRESSES  FOR BREAKING UP WORK FOR USE
C     IN EVALU.
      LBASE(1) = 1
      DO  I = 2,7
          LBASE(I) = LBASE(I-1) + LU
      END DO
      LD1 = LBASE(1)
      LD2 = LBASE(2)
      LD3 = LBASE(3)
      LD4 = LBASE(4)
      LD5 = LBASE(5)
      LD6 = LBASE(6)
      LD7 = LBASE(7)
C     MOVE  RHO,X,Y,...,D2Y AROUND IN WORK, LENGTHEN EACH OF THEM.
      DO I = 1,7
          IF (LBASE(I)+NFINAL-1 .GE. NBASE(I+2)) THEN
              DO K = I,7
                  II = I + 7 - K
                  LB = LBASE(II) - 1
                  NB = NBASE(II+1) - 1
                  DO J = 1,NFINAL
                      WORK(LB+J) = WORK(NB+J)
                  END DO
              END DO
              GO TO 10
          END IF
          LB = LBASE(I) - 1
          NB = NBASE(I+1) - 1
          DO J = 1,NFINAL
              WORK(LB+J) = WORK(NB+J)
          END DO
      END DO 
C
   10 CALL EVALU(IOUT,IDBG,IBEG,IE,B,BDYFCN,CURVE,NFINAL,WORK(LD1),
     +           EP,R_FORM,WORK(LD2),WORK(LD3),WORK(LD4),WORK(LD5),
     +           WORK(LD6),WORK(LD7),LU,DPTS,NPTS,UVEC,ERROR,IER2)
C
      IER = MAX0(IER1,IER2)
      DO I = 1,NPTS
          ERROR(I) = ERROR(I) + SIGN(RHOERR, ERROR(I))
      END DO
      IF (IER1 .EQ. 0) RETURN
      DO I = 1,NPTS
          IF (ERROR(I) .GT. ZERO) ERROR(I) = -ERROR(I)
      END DO
      RETURN
 9000 FORMAT (/,' FROM SUBROUTINE DRCHLT: SUBROUTINE INTEQN RESULTS.',
     +   /,' NFINAL=',I3,5X, 'RHOERR=',1P,E8.2,5X,'IER1=',I1,/)
 9010 FORMAT (6X,'X',14X,'Y',19X,'RHO')
 9020 FORMAT (1P,D12.4,D15.4,D25.12)
      END

      SUBROUTINE EVALU(IOUT,IDBG,IBEG,IE,B,BDYFCN,CURVE,N,RHO,EPS,
     +                 R_FORM,X,Y,DX,DY,D2X,D2Y,LU,DPTS,NP,
     +                 U,ERROR,IER)
C     ----------------
C
C     This program evaluates the double layer potential U at the
C     given points in DPTS. The input is the density function RHO,
C     and DPTS at which U is evaluated. RHO is evaluated in the
C     subroutine INTEQN. These results are stored in U, along
C     with the predicted error bound in ERROR. The desired error
C     tolerance is EPS.  If the desired error bound is not attained,
C     then the corresponding entry in ERROR is made negative. Its
C     magnitude is still an estimated error bound.
C
C     IOUT      Input.
C               The output unit number, to the file DRCHLT.ANS
C     IDBG      Input.
C               The dubugging parameter. 
C               =0 produces a shortcut output file.
C               =1 produces a debugging output.
C     IBEG      Input.
C               =0 means this is a first call on DRCHLT for this
C               particular curve C, boundary function BDYFCN, and
C               error tolerance EPS.
C               =1 means that DRCHLT has been called previously for
C               this choice of parameters, and the solution U is
C               desired at a new set of points in DPTS. For such a
C               call on DRCHLT, change only DPTS, NPTS, and IBEG.
C               Do not change any other input variable.  
C     IE        Input.
C               =0 for the interior Dirichlet problem.
C               =1 for the exterior Dirichlet problem.
C     CURVE     External subroutine.
C               This program defines the curve of the boundary C.
C     B         Input.
C               For the parameterization of C defined in CURVE,
C               the parameterization interval is [0,B].
C     BDYFCN    External subroutine.
C               This program defines the Dirichlet data on the
C               boundary.
C     N         Input.
C               This is NFINAL as output from the subroutine INTEQ.
C     RHO       Input.
C               An array which contains the value of the double
C               layer density function defining U.
C     EPS       Input.
C               The user-supplied absolute error tolerance.
C     R_FORM    Input.
C               This specifies the way in which the variable RATE
C               is to be defined.
C               If R_FORM=0, we use the "normal" way to define RATE
C               based on estimating the rate of convergence in the
C               approximates calculated to date.  
C               If R_FORM=1, then we use a "conservative" error test
C               in which RATE=0.5 and the approximates are assumed to
C               have a very slow rate of convergence.  Use this for 
C               more ill-behaved problems and boundaries.
C     X,Y       Inputs.
C               Two arrays containing a sequence of points
c               (X(I),Y(I)) produced by calling the subroutine
C               NEWCUR.  They correspond to an even subdivision
C               of the parameterization interval [0,B].
C     DX,DY     Inputs.
C               The derivative values corresponding to the points
C               given in X,Y.
C     D2X,D2Y   Inputs.
C               The second derivative values corresponding to the
C               points given in X,Y.
C     LU        Input.
C               The upper bound of the size of the arrays X,Y,DX,DY,
C               D2X,D2Y,RHO.
C     DPTS      Input.
C               This is a two-dimensional array which supplies the
C               points at which U is to be evaluated.
C     NP        Input 
C               This is the number of points in DPTS.
C     U         Output.
C               An output array which contains U(P) for points P
C               given in DPTS.
C     ERROR     Output.
C               An output array which contains the predicted error
C               bound for the corresponding entries in U.
C     IER       Output.
C               =0 means the program was completed satisfactorily.
C               =1 means some or all of the solution values in U do
C               not satisfy the error tolerance EPS.
C
      INTEGER IBEG,IDBG,IE,IER,IOUT,LU,N,NP,R_FORM
      DOUBLE PRECISION DPTS(2,NP),X(LU),Y(LU),DX(LU),DY(LU),D2X(LU),
     +                D2Y(LU),ERROR(NP),RHO(LU),U(NP)
      DOUBLE PRECISION D1MACH,B,BDYFCN,EPS
      EXTERNAL BDYFCN,CURVE
C
C     LOCAL VARIABLES
      DOUBLE PRECISION AIT,BIT,CMIN,CUT,D2F,DF,DIFF,DISTPQ,DNORM,ERR,
     +                 FCNK,FCNM,FCNU,H,MESH,OLDIFF,OLDU,ONE,PASTRT,
     +                 PI,PROD,PX,PY,Q,QD2X,QD2Y,QDX,QDY,QSPD,QX,QY,R,
     +                 RATE,RTLOW,RTUP,S,SB,SLOPE,SONE,SUM,SZERO,T1,T2,
     +                 TQX,TQY,TR,TWO,TX,TY,VALUE,ZERO,PARM
      INTEGER I,ITR,J,JL,JM,JMIN,JSTEP,JU,K,KH,KSTEP,L,LD,LDM1,LOOP,
     +        LOOP1,LOOP2,M_0,M2
      PARAMETER (M_0=32)
C
C     M_0/2 DENOTES THE INITIAL NUMBER OF INTEGRATION NODES TO BE USED
C     IN THE EVALUATION OF THE DOUBLE LAYER POTENTIAL, AND IT WILL ALSO
C     BE PERFORMED WITH M_0 NODES, SO THAT THE VALUE OF NWORK IN THE 
C     CALLING PROGRAM NEEDS TO BE SET ACCORDINGLY.  ALWAYS SET M_0 TO
C     BE A POWER OF 2.
C
      DOUBLE PRECISION SQUARE(M_0)
      EXTERNAL NEWCUR
      INTRINSIC ABS,MAX,MIN,MIN0
      DATA ZERO/0.0D0/,ONE/1.0D0/,RTUP/0.5D0/,RTLOW/.1D0/,CUT/.01D0/,
     +     TWO/2.0D0/
C
C     DATA 'PI'
      PI=4.D0*ATAN(1.D0)     
C
C     INITIALIZE. LD IS THE RUNNING DIMENSION OF RHO,X,...,D2Y.
      IF (IBEG .EQ. 0) LD = N
      IER = 0
      IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9000)
C     BEGIN LOOP TO EVALUATE U AT POINTS P IN DPTS.
      DO I = 1,NP
          PX = DPTS(1,I)
          PY = DPTS(2,I)
C
          IF (IDBG .EQ. 1)  THEN
              WRITE (IOUT,FMT=9010) I,DPTS(1,I),DPTS(2,I),LD
          END IF 
C
C     IF IT IS THE EXTERIOR PROBLEM, CHANGE THE DPTS BY USING
C     THE KELVIN TRANSFORMATION.
          IF (IE .EQ. 1) THEN
              R = PX*PX + PY*PY
              PX = PX/R
              PY = PY/R
          END IF
C
C     BEGIN THE CALCULATION OF THE POINT Q OF C WHICH IS CLOSEST TO P.
C     SEARCH USING M2 EVENLY SPACED POINTS OF C, GIVEN IN X,Y, AT 
C     STEPS OF JSTEP. INITIALLY CALCULATE SQUARES OF APPROPRIATE 
C     ANGELS.
          M2 = MIN0(LD,M_0)
          JSTEP = LD/M2
          DO J = 1,M2
              JM = J*JSTEP
              T1 = X(JM) - PX
              T2 = Y(JM) - PY
              SQUARE(J) = T1*T1 + T2*T2
          END DO
          CMIN = 1.0D50
          DO J = 1,M2
              IF (SQUARE(J) .LT. CMIN) THEN
                  CMIN = SQUARE(J)
                  JMIN = J
              END IF
          END DO
C     THE POINT (X(K),Y(K)),K=JSTEP*JMIN, IS CLOSEST TO AMONG ALL
C     M2 POINTS SEARCHED.
          JM = JMIN*JSTEP
          T1 = X(JM) - PX
          T2 = Y(JM) - PY
          PROD = ((T1*DX(JM)+T2*DY(JM))**2)/
     +           (CMIN* (DX(JM)**2+DY(JM)**2))
          IF (PROD .GT. CUT) THEN
              IF (IDBG .EQ. 1) THEN
                  WRITE(IOUT,*) 'STEP 1 FOR FINDING Q FAILED. PROD=', 
     +                           PROD
              END IF
              GO TO 10
          END IF
C     THIS POINT IS ACCEPTABLY CLOSE, AND WILL BE CALLED Q.
C     EVALUATE RHO AT Q AND SAVE.
          VALUE = RHO(JSTEP*JMIN)
          Q = B/LD*(JSTEP*JMIN)         
C     FOR AN EXTERIOR PROBLEM, OBTAIN THE CORRESPONDING POINT ON THE
C     GIVEN CURVE USING THE KELVIN TRANSFORMATION.
          R = X(JM)*X(JM) + Y(JM)*Y(JM)
          TR = ONE
          IF (IE .EQ. 1) TR = ONE/R
          QX = TR*X(JM)
          QY = TR*Y(JM)
          IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9020) QX,QY,VALUE
          GO TO 60
C     LOOK MORE CAREFULLY AT THE POINTS IN X,Y FOR A CLOSEST POINT.
C     INITIALLY, SELECT AN INTERVAL IN WHICH TO SEARCH.
   10     IF (JMIN .EQ. M2) THEN
              JMIN = JMIN*JSTEP
              JL = JMIN - JSTEP + 1
              JU = JMIN
              DO J = JL,JU
                  T1 = X(J) - PX
                  T2 = Y(J) - PY
                  DNORM = T1*T1 + T2*T2
                  IF (DNORM .LT. CMIN) THEN
                      CMIN = DNORM
                      JMIN = J
                  END IF
              END DO
              JL = 1
              JU = JSTEP
          ELSE
              JMIN = JMIN*JSTEP
              JL = JMIN - JSTEP + 1
              JU = JMIN + JSTEP - 1
          END IF
          DO J = JL,JU
              T1 = X(J) - PX
              T2 = Y(J) - PY
              DNORM = T1*T1 + T2*T2
              IF (DNORM .LT. CMIN) THEN
                  CMIN = DNORM
                  JMIN = J
              END IF
          END DO
          T1 = X(JMIN) - PX
          T2 = Y(JMIN) - PY
          PROD = ((T1*DX(JMIN)+T2*DY(JMIN))**2)/
     +           (CMIN* (DX(JMIN)**2+DY(JMIN)**2))
          IF (PROD .GT. CUT) THEN
              IF (IDBG .EQ. 1) THEN
                  WRITE(IOUT,*) 'STEP 2 FOR FINDING Q FAILED.
     +                           PROD=', PROD
              END IF 
              GO TO 20
          END IF
C     THIS POINT Q=(X,Y) IS ACCEPTABLY CLOSE TO P. EVALUATE RHO AT Q.
          VALUE = RHO(JMIN)
          Q = B/LD*JMIN 
C     FOR AN EXTERIOR PROBLEM, OBTAIN A  CORRESPONDING POINTS ON THE
C     GIVEN CURVE USING THE KELVIN TRANSFORMATION.
          TR = ONE
          R = X(JMIN)*X(JMIN) + Y(JMIN)*Y(JMIN)
          IF (IE .EQ. 1) TR = ONE/R
          QX = TR*X(JMIN)
          QY = TR*Y(JMIN)
          IF(IDBG .EQ. 1) WRITE (IOUT,FMT=9030) QX,QY,VALUE
          GO TO 60
C     NO ACCEPTABLE POINT Q FOUND ON C USING THE VECTORS X,Y.
C     BEGIN ITERATION METHOD.
   20     H = B/LD
          SB = JMIN*H
          SZERO = SB
          LOOP = 0
          QX = X(JMIN)
          QY = Y(JMIN)
          QDX = DX(JMIN)
          QDY = DY(JMIN)
          QD2X = D2X(JMIN)
          QD2Y = D2Y(JMIN)
          T1 = QX - PX
          T2 = QY - PY
C     ITERATION LOOP.
   30     DF = TWO* (T1*QDX+T2*QDY)
          D2F = TWO* (QDX**2+QDY**2+T1*QD2X+T2*QD2Y)
          SONE = SZERO - DF/D2F
          LOOP = LOOP + 1
          CALL NEWCUR(IE,CURVE,SONE,QX,QY,QDX,QDY,QD2X,QD2Y)
          T1 = QX - PX
          T2 = QY - PY
          PROD = ((T1*QDX+T2*QDY)**2)/ ((T1*T1+T2*T2)*
     +           (QDX*QDX+QDY*QDY))
          IF (IDBG .EQ. 1) THEN
              WRITE(IOUT,*)'STEP 3 FOR FINDING Q. PROD=', PROD
          END IF
          IF (PROD .LE. CUT) THEN
             Q=SONE
             GO TO 50
          END IF
C     THE NEW POINT (QX,QY) IS NOT SUFFICIENTLY CLOSE. CHECK FOR
C     POSSIBLE DIVERGENCE OF ITERATION.
          IF ((SB-H .GT. SONE) .OR. (SB+H .LT. SONE)) GO TO 40
C     CONTINUE ITERATION.
          SZERO = SONE
          GO TO 30
C     PRIMARY ITERATION IS DIVERGING. GO TO A METHOD GUARANTEED
C     TO CONVERGE, THE BISECTION METHOD
   40     AIT = SB - H
          BIT = SB + H
C     BEGINNING OF ITERATION LOOP. WE LIMIT THE NUMBER OF ITERATION
C     TO A CERTAIN NUMBER( HERE 10)
          DO ITR = 1,10
              Q = (AIT+BIT)/2.0
              LOOP = LOOP + 1
              CALL NEWCUR(IE,CURVE,Q,QX,QY,QDX,QDY,QD2X,QD2Y)
              T1 = QX - PX
              T2 = QY - PY
              SLOPE = T1*QDX + T2*QDY
              PROD = (SLOPE*SLOPE)/ ((T1*T1+T2*T2)* (QDX*QDX+QDY*QDY)) 
              IF (IDBG .EQ. 1) WRITE(IOUT,*)'STEP 4 FOR FINDING Q. 
     +                         PROD=', PROD
              IF (PROD .LE. CUT) THEN
                 PARM=Q
                 GO TO 50
              END IF
              SLOPE = 2.0*SLOPE
              IF (SLOPE .LT. ZERO) THEN
                  AIT = Q
              ELSE
                  BIT = Q
              END IF
          END DO
C     THE ITERATION IS IN A TIGHT LOOP. SOMETHING IS WRONG ABOUT 
C     CURVE. FOR THE EXTERIOR PROBLEM, OBTAIN A CORRESPONDING POINT
C     ON GIVEN CURVE.
          IF (IE .EQ. 1) THEN
              R = QX*QX + QY*QY
              QX = QX/R
              QY = QY/R
          END IF
      IF(IDBG.EQ.1) WRITE (IOUT,FMT=9040) Q,QX,QY
C     A SUFFICIENTLY ACCURATE CLOSEST POINT TO P HAS BEEN FOUND, AND
C     NOW EVALUATE RHO AT THIS POINT, USING THE NYSTROM INTERPOLATION
C     FORMULA.
   50     SUM = ZERO
          KSTEP = LD/N
          DO K = KSTEP,LD,KSTEP
              T1 = X(K) - QX
              T2 = Y(K) - QY
              SUM = SUM + RHO(K)* (DY(K)*T1-DX(K)*T2)/ (T1*T1+T2*T2)
          END DO
          TQX = QX
          TQY = QY
C     FOR AN EXTERIOR PROBLEM, USE THE KELVIN TRANSFORMATION TO OBTAIN
C     THE CORRESPONDING POINT ON THE GIVEN CURVE.
          IF (IE .EQ. 1) THEN
              R = QX**2 + QY**2
              TQX = QX/R
              TQY = QY/R
          END IF
          VALUE = - (BDYFCN(TQX,TQY)+KSTEP*H*SUM)/PI
C         IF (IDBG .EQ. 1)  WRITE (IOUT,FMT=9050) TQX,TQY,VALUE,LOOP
C     CLOSEST POINT Q AND THE VALUE OF RHO AT Q HAVE BEEN CALCULATED.
C     NOW BEGIN EVALUATION OF U(PX,PY) USING NUMERICAL INTEGRATION.
C     INITIALIZE, AND BEGIN WITH M_0/2 NODES.
C
   60     CALL NEWCUR(IE,CURVE,Q,QX,QY,QDX,QDY,QD2X,QD2Y)
          DISTPQ = SQRT((PX-QX)*(PX-QX) + (PY-QY)*(PY-QY))
          QSPD = SQRT(QDX*QDX+QDY*QDY)
          RATE = RTUP
          PASTRT = RTUP
          L = M_0/2
          LOOP1 = 1
          LOOP2 = 1
C     CALCULATE NUMERICAL INTEGRAL WITH L SUBDIVISIONS OF (0,B).
   70     SUM = ZERO
          KSTEP = LD/L
          DO K = KSTEP,LD,KSTEP
              T1 = X(K) - PX
              T2 = Y(K) - PY
              FCNM = (DY(K)*T1-DX(K)*T2)/ (T1*T1+T2*T2)
              SUM = SUM + FCNM* (RHO(K)-VALUE)
          END DO
          FCNU = -TWO*PI*VALUE - (B/L)*SUM
          IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9060) FCNU,L,LD
          IF (LOOP1 .EQ. 1) GO TO 100
C         ESTIMATE ERROR IN FCNU.
          DIFF = ABS(FCNU-OLDU)
          IF (LOOP1 .EQ. 2) GO TO 80
C         UPDATE RATE OF CONVERGENCE OF NUMERICAL INTEGRATION.
          IF(R_FORM .EQ. 0) THEN
C             THE FOLLOWING IS A SOPHISTICATED ERROR ESTIMATOR,
C             USUALLY REASONABLY ACCURATE.
              RATE = MAX(PASTRT,RTLOW,MIN(RTUP,ABS(DIFF/OLDIFF)))
          ELSE
C             THE FOLLOWING IS A CONSERVATIVE ERROR ESTIMATOR.
              RATE = RTUP
          END IF
          PASTRT = MIN(RTUP,ABS(DIFF/OLDIFF))
   80     ERR = (RATE/ (ONE-RATE))*DIFF
          IF(IDBG .EQ. 1) WRITE (IOUT,FMT=9070) ERR,RATE
C
          IF (ERR .GT. EPS) GO TO 90
          IF (ERR .EQ. ZERO) THEN
              ERR = D1MACH(4)
              OLDIFF = ERR
          END IF
C     FOR A POINT CLOSE TO THE BOUNDARY, ITERATE TWO MORE TIMES.
          MESH = QSPD*B/L 
          IF ((MESH .GT. DISTPQ) .AND. (LOOP2 .LE. 2)) THEN
              LOOP2 = LOOP2 + 1
              GO TO 90
          END IF

C     THE VALUE OF FCNU IS SUFFICIENTLY ACCURATE.
          U(I) = FCNU
          ERROR(I) = ERR
          GO TO 120
C     FCNU IS NOT SUFFICIENTLY ACCURATE.
C     RE-INITIALIZE FOR ANOTHER NUMERICAL INTEGRATION.
   90     OLDIFF = DIFF
          IF(OLDIFF .EQ. ZERO) OLDIFF = D1MACH(4)
  100     OLDU = FCNU
          LOOP1 = LOOP1 + 1
          L = 2*L
          IF (L .LE. LD) GO TO 70
C     NOT SUFFICIENT VALUES IN RHO. THUS VALUES OF RHO ON A FINER
C     MESH MUST BE CREATED.
          LD = 2*LD
          IF (LD .GT. LU) GO TO 110
C     THERE IS SUFFICIENT SPACE IN RHO,X,Y,...,D2Y FOR AN INCREASED SUB-
C     DIVISION OF (0,B).
C     MOVE OLD VALUES OF RHO,...,D2Y TO MAKE ROOM FOR NEW VALUES.
          DO J = 2,LD,2
              K = LD + 2 - J
              KH = K/2
              RHO(K) = RHO(KH)
              X(K) = X(KH)
              Y(K) = Y(KH)
              DX(K) = DX(KH)
              DY(K) = DY(KH)
              D2X(K) = D2X(KH)
              D2Y(K) = D2Y(KH)
          END DO
C     PRODUCE NEW CURVE PARAMETERS FOR FINER SUBDIVISION.
          H = B/LD                                                     
          LDM1 = LD - 1
          DO J = 1,LDM1,2
              S = J*H
              CALL NEWCUR(IE,CURVE,S,X(J),Y(J),DX(J),DY(J),D2X(J),
     +                    D2Y(J))
          END DO
C     PRODUCE NEW VALUES OF RHO.
          H = B/N
          KSTEP = LD/N
          DO J = 1,LDM1,2
              SUM = ZERO
              DO K = KSTEP,LD,KSTEP
                  T1 = X(K) - X(J)
                  T2 = Y(K) - Y(J)
                  FCNK = (DY(K)*T1-DX(K)*T2)/ (T1*T1+T2*T2)
                  SUM = SUM + FCNK*RHO(K)
              END DO
C     FOR AN EXTERIOR PROBLEM, USE THE KELVIN TRANSFORMATION.
              TX = X(J)
              TY = Y(J)
              IF (IE .EQ. 1) THEN
                  R = X(J)*X(J) + Y(J)*Y(J)
                  TX = X(J)/R
                  TY = Y(J)/R
              END IF
              RHO(J) = - (BDYFCN(TX,TY)+H*SUM)/PI
          END DO
          GO TO 70
C     THE UPPER LIMITS FOR RHO,X,Y,...,D2Y HAVE BEEN REACHED.
C     MARK ERROR BOUND ACCORDINGLY AND CONTINUE ONTO NEXT POINT P.
  110     ERROR(I) = -ERR
          U(I) = FCNU
          IER = 1
          LD = LD/2
  120 END DO
      RETURN
 9000 FORMAT (/,' FROM SUBROUTINE EVALU.',/)
 9010 FORMAT (/,' I=',I3,4X,'PX=',1P,D11.4,3X,'PY=',D11.4,5X,'LD=',I6)
 9020 FORMAT (' QSTAGE1. QX=',1P,D11.4,3X,'QY=',D11.4,3X,'RHO=',
     +        D20.12)
 9030 FORMAT (' QSTAGE2. QX=',1P,D11.4,3X,'QY=',D11.4,3X,'RHO=',
     +        D20.12)
 9040 FORMAT (' PROD IS NOT CONVERGING TO ZERO IN LOOP BEGINNING AT 
     +        112',/,' Q=',1P,D11.4,5X,'QX=',D11.4,5X,'QY=',/)
C9050 FORMAT (' QSTAGE3. QX=',1P,D11.4,3X,'QY=',D11.4,3X,'RHO=',
C    +        D20.12,3X,'LOOPS=',I1)
 9060 FORMAT (' NUM INT =',1P,E20.12,5X,'L=',I6,5X,'LD=',I6)
 9070 FORMAT (' ERROR=',1P,D20.12,5X,'RATE=',D11.4)
      END

      SUBROUTINE INTEQN(IOUT,IDBG,IE,B,EPS,R_FORM,BDYFCN,CURVE,NUPPER,
     +                  RHO,ERROR,NFINAL,X,Y,DX,DY,D2X,D2Y,
     +                  OLDRHO,IWORK,WORK,KERMAT,IER)
C     -----------------
C
C     This program solves the second kind boundary integral equation
C     which arises from solving the interior Dirichlet problem as a
C     double layer potential.
C
C     The integral equation is solved using Nystrom's method with
C     the rectangular rule as the quadrature rule.   The resulting
C     linear system is solved directly using LAPACK routines.
C
C     The output is the double layer density function RHO.  This
C     is to be found with such accuracy that the resulting harmonic
C     function has an accuracy of EPS.
C
C     This routine assumes the boundary C is at least two times
C     continuously differentiable.  The boundary C is defined by
C     the subroutine CURVE.
C
C     The present routine calculates with the rectangular rule for
C     N=4,8,16,... until a sufficiently accurate value of RHO is
C     obtained.  This is subject to N .LE. NUPPER, with the latter
C     based on the size of the vector WORK supplied by the user in
C     calling DRCHLT.
C
C     IOUT     Input.
C              The output unit number, to the file DRCHLT.ANS
C     IDBG     Input.
C              The debugging parameter. 
C              =0 produces a shortcut output.
C              =1 produces a debugging output.
C     IE       Input.
C              =0 for the interior Dirichlet problem.
C              =1 for the exterior Dirichlet problem.
C     CURVE    Input.
C              This program defines the curve of the boundary C.
C     B        Input.
C              For the parameterization of C defined in CURVE,
C              the parameterization interval is [0,B].
C     EPS      Input.
C              The user-supplied error tolerance.
C     R_FORM   Input.
C              This specifies the way in which the variable RATE
C              is to be defined.
C              If R_FORM=0, we use the "normal" way to define RATE
C              based on estimating the rate of convergence in the
C              approximates calculated to date.  
C              If R_FORM=1, then we use a "conservative" error test
C              in which RATE=0.5 and the approximates are assumed to
C              have a very slow rate of convergence.  Use this for 
C              more ill-behaved problems and boundaries.
C     BDYFCN   Input.
C              This program defines the Dirichlet data on the
C              boundary.
C     NUPPER   Input.
C              This is the upper bound for the size of linear
C              system that can be constructed and solved.
C     RHO      Output.
C              An array which contains the value of the double
C              layer density function defining U.
C     ERROR    Output.
C              An output array which contains the predicted error
C              bound for the corresponding entries in U.
C     NFINAL   Output.
C              This is the dimension of the final linear system
C              constructed in solving for RHO.
C     X,Y      Output.
C              Two arrays containing a sequence of points
c              (X(I),Y(I)) produced by calling the subroutine
C              NEWCUR.  They correspond to an even subdivision
C              of the parameterization interval [0,B].
C     DX,DY    Output.
C              The derivative values corresponding to the points
C              given in X,Y.
C     D2X,D2Y  Output.
C              The second derivative values corresponding to the
C              points given in X,Y.
C     OLDRHO   Output.
C              An array containing the preceding value of RHO,
C              also produced in this program.
C     IWORK    Integer work space
C              This is an array for pivoting used for the subroutines
C              in LAPACK.
C     WORK     Real work space.
C              A work array for the subroutines in LAPACK.           
C     KERMAT   Output.
C              This is array contains the linear system associated
C              with the Nystrom method.
C     IER      Output.
C              =0 means the program was completed satisfactorily.
C              =1 means some or all of the solution values in U do
C              not satisfy the error tolerance EPS.
C
      INTEGER IDBG,IE,IER,INFO,IOUT,NFINAL,NUPPER,R_FORM
      DOUBLE PRECISION D2X(NUPPER),D2Y(NUPPER),DX(NUPPER),DY(NUPPER),
     +                 KERMAT(NUPPER,NUPPER),OLDRHO(NUPPER),
     +                 RHO(NUPPER),X(NUPPER),Y(NUPPER),WORK(4*NUPPER)
      INTEGER IWORK(NUPPER)
      DOUBLE PRECISION D1MACH,RTLOW,B,BDYFCN,EPS,ERROR
      EXTERNAL BDYFCN,CURVE
C
C     LOCAL VARIABLES.
      DOUBLE PRECISION DIFF,DIST,H,OLDIFF,ONE,PI,R,RATE,RCOND,
     +                 RTUP,SUM,SUMAX,T1,T2,TWO,TX,TY,ZERO
      INTEGER I,J,JH,N,NM1,N_0,NRHS
      EXTERNAL NEWCUR
      INTRINSIC ABS,MAX
      PARAMETER (N_0 = 32)
      DATA RTLOW/.1D0/,RTUP/.5D0/,ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,
     +     NRHS/1/
C
C     THE PARAMETER N_0 GIVES THE INITIAL NUMBER OF QUADRATURE POINTS
C     USED IN THE APPROXIMATION OF THE INTEGRAL EQUATION.  THE EQUATION
C     WILL ALSO BE SOLVED WITH 2*N_0 NODES.  IN THE PROGRAM CALLING
C     NEUMAN, THE PARAMETER NWORK SHOULD BE SET ACCORDING.  ALWAYS SET 
C     N_0 TO BE A POWER OF 2.
C
C     DATA 'PI'
      PI = 4.D0*ATAN(1.D0)
C     INITIAL CASE, N=N_0. INITIALIZE PARAMETERS.
      IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9000)
      N = N_0
      RATE = RTUP
C     DEFINE STEPSIZE AND POINTS ON CURVE.
      H = B/N
      DO I = 1,N
          CALL NEWCUR(IE,CURVE,I*H,X(I),Y(I),DX(I),DY(I),D2X(I),
     +                D2Y(I))
      END DO
      GO TO 20
C     DEFINE H AND POINTS ON CURVE FOR CONTINUING LOOP ON N.
   10 H = B/N
      DO I = 2,N,2
          J = N + 2 - I
          JH = J/2
          X(J) = X(JH)
          Y(J) = Y(JH)
          DX(J) = DX(JH)
          DY(J) = DY(JH)
          D2X(J) = D2X(JH)
          D2Y(J) = D2Y(JH)
      END DO
      NM1 = N - 1
      DO I = 1,NM1,2
          CALL NEWCUR(IE,CURVE,I*H,X(I),Y(I),DX(I),DY(I),D2X(I),
     +                D2Y(I))
      END DO
C     SET UP MATRIX EQUATION.
   20 DO I = 1,N
          RHO(I) = BDYFCN(X(I),Y(I))
          IF (IE .EQ. 1) THEN
C     FOR AN EXTERIOR DIRICHLET PROBLEM, EVALUATE THE BOUNDARY DATA
C     ON THE ORIGINALLY GIVEN  CURVE BY USING THE KELVIN 
C     TRANSFORMATION.
              R = X(I)*X(I) + Y(I)*Y(I)
              TX = X(I)/R
              TY = Y(I)/R
              RHO(I) = BDYFCN(TX,TY)
          END IF
          DO J = 1,N
              IF (I .EQ. J) THEN
C             DEFINE KERNEL FOR T(I) = T(J).
                  T1 = DX(I)
                  T2 = DY(I)
                  DIST = T1*T1 + T2*T2
                  KERMAT(I,I) = -PI - H* (T1*D2Y(I)-T2*D2X(I))/
     +                          (TWO*DIST)
              ELSE
C             DEFINE KERNEL FOR T(I) .NE. T(J)
                  T1 = X(J) - X(I)
                  T2 = Y(J) - Y(I)
                  DIST = T1*T1 + T2*T2
                  KERMAT(I,J) = -H* (DY(J)*T1-DX(J)*T2)/DIST
              END IF
          END DO
      END DO
      IF (N .EQ. N_0) THEN
          SUMAX = 1.D0 
          GO TO 30
      END IF
C     CALCULATE PI+NORM(INTEGRAL OPERATOR).
      SUMAX = ZERO
      DO I = 1,N
          SUM = ZERO
          DO J = 1,N
              IF (I .EQ. J) THEN
                  SUM = SUM + ABS(PI+KERMAT(I,I)) + PI
              ELSE
                  SUM = SUM + ABS(KERMAT(I,J))
              END IF
          END DO
          SUMAX = MAX(SUMAX,SUM)
      END DO
   30 CALL DGESV(N,NRHS,KERMAT,NUPPER,IWORK,RHO,NUPPER,INFO)
      CALL DGECON('I',N,KERMAT,NUPPER,SUMAX,RCOND,WORK,IWORK,INFO)
      IF(RCOND .EQ. ZERO) RCOND = D1MACH(4)
      IF (IDBG .EQ. 1) WRITE(IOUT,9030) ONE/RCOND
      IF (N .EQ. N_0) GO TO 60
C     CALCULATE NORM OF RHO-OLDRHO.
      DIFF = ZERO
      DO I = 2,N,2
          DIFF = MAX(DIFF,ABS(RHO(I)-OLDRHO(I/2)))
      END DO
      IF (N .EQ. 2*N_0) GO TO 40
C     MEASURE RATE OF CONVERGENCE.
      RATE = DIFF/OLDIFF
      IF(RATE .GE. ONE) THEN
          IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9020) N,ERROR,EPS
          IF (2*N .GT. NUPPER) THEN
C             EXIT FOR UNSUCCESSFUL RETURN.
              IER = 1
              NFINAL = N
              RETURN
          ELSE
              GO TO 50
          END IF
      END IF
      IF(R_FORM .EQ. 0) THEN
C         THE FOLLOWING IS A SOPHISTICATED ERROR ESTIMATOR,
C         USUALLY REASONABLY ACCURATE.
          RATE=MAX(RTLOW,MIN(RATE,RTUP))
      ELSE
C         THE FOLLOWING USES A RATE THAT ASSUMES THE RATE OF
C         CONVERGENCE IS PROPORTIONAL TO 1/N.  IT IS A CONSERVATIVE,
C         BUT SAFER, CHOICE.
          RATE = RTUP
      END IF
C     ESTIMATE ERROR IN RHO.
   40 ERROR = (ONE/RCOND)*(SUMAX)*DIFF*RATE/ (ONE-RATE)
      IF (IDBG .EQ. 1) WRITE (IOUT,FMT=9010) N,RATE,ERROR,DIFF,SUMAX
      IF (ERROR .LE. EPS) THEN
C         EXIT FOR SUCCESSFUL RETURN.
          NFINAL = N
          IER = 0
          RETURN
      ELSE IF (2*N .GT. NUPPER) THEN
C         EXIT FOR UNSUCCESSFUL RETURN.
          IER = 1
          NFINAL = N
          RETURN
      END IF
C     PREPARE FOR ANOTHER LOOP ON N.
   50 OLDIFF = DIFF
      IF(OLDIFF .EQ. ZERO) OLDIFF = D1MACH(4)
   60 DO I = 1,N
          OLDRHO(I) = RHO(I)
      END DO
      N = 2*N
      GO TO 10
 9000 FORMAT (/,' FROM SUBROUTINE INTEQN',/)
 9010 FORMAT (' N=',I3,3X,'RATE=',1P,D8.2,3X,'ERROR=',D8.2,3X,'DIFF=',
     +       D8.2,3X,'SUMAX=',D8.2)
 9020 FORMAT (' N=',I3,3X,'RATE > 1',25X,'ERROR=',1PD8.2,3X,
     +       'EPS=',D8.2)
 9030 FORMAT (/,' CONDITION NUMBER = ',1PD8.2)
      END

      SUBROUTINE KVTRNF(X,Y,DX,DY,D2X,D2Y,T,DT,D2T)
C     -----------------
C
C     Define the Kelvin transformation.
C
C     INPUTS: X, Y, DX, DY, D2X, D2Y
C     OUTPUTS: T, DT, D2T
C
      DOUBLE PRECISION D2T,D2X,D2Y,DIST,DT,DX,DY,T,X,Y
      DIST = X*X + Y*Y
      T = X/DIST
      DT = (DX* (Y*Y-X*X)-2*X*Y*DY)/ (DIST*DIST)
      D2T = D2X* (Y**4-X**4) - 2*X*DIST* (DY*DY+DX*DX+Y*D2Y)
      D2T = D2T - 4* (X*DX+Y*DY)* (DX* (Y*Y-X*X)-2*X*Y*DY)
      D2T = D2T/DIST**3
      RETURN
      END

      SUBROUTINE NEWCUR(IE,CURVE,S,X,Y,DX,DY,D2X,D2Y)
C     -----------------
C
C     If IE=0, the resulting curve C will be the same as that
C     defined in the subroutine CURVE.
C     If IE=1, the resulting curve will be that produced by
C     applying the Kelvin transformation to the original
C     boundary curve C.
C
C     INPUTS: IE, S
C     EXTERNAL SUROUTINE: CURVE
C     OUTPUTS: X, Y, DX, DY, D2X, D2Y
C
      DOUBLE PRECISION D2TX,D2TY,D2X,D2Y,DTX,DTY,DX,DY,S,TX,TY,X,Y
      INTEGER IE
      EXTERNAL CURVE,KVTRNF
      CALL CURVE(S,X,Y,DX,DY,D2X,D2Y)
      IF (IE .EQ. 1) THEN
          CALL KVTRNF(X,Y,DX,DY,D2X,D2Y,TX,DTX,D2TX)
          CALL KVTRNF(Y,X,DY,DX,D2Y,D2X,TY,DTY,D2TY)
          X = TX
          Y = TY
          DX = DTX
          DY = DTY
          D2X = D2TX
          D2Y = D2TY
      END IF
      RETURN
      END
