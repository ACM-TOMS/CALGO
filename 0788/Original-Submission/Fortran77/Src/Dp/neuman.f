
      SUBROUTINE NEUMAN(IOUT,IDBG,IE,BDYFCN,CURVE,B,DPTS,NPTS,IBEG,EPS0,
     +                  WORK,IWORK,NWORK,MWORK,MAXFFT,UVEC,ERROR,IER)
C     -----------------
C
C     This program will solve the interior and exterior Neumann
C     problems for Laplace's equation on a simply-connected region D
C     with smooth boundary C.  It is assumed that C is at least two
C     times continuously differentiable. To define a boundary, the user
C     must supply a parametric form of the external routine CURVE with 
C     parameter range [0,B]. The routine CURVE must also produce the 
C     first and second derivatives on the curve as defined below. The
C     boundary condition must be given by the routine BDYFCN.
C
C     SUBROUTINE PARAMETERS:
C     IOUT      Input.
C               The output unit number.
C     IDBG      Input.
C               The debugging parameter.
C               =0 produces a standard output.
C               =1 produces a debugging output.
C     IE        Input
C               =0 for the interior Neumann problem.
C               =1 for the exterior Neumann problem.
C     BDYFCN    External function.
C               Inputs : X,Y
C               Outputs: BDYFCN
C               This is a function of two variables X,Y.  It gives
C               the value of the normal derivative of the unknown
C               harmonic function U on the boundary curve C. This
C               is a user supplied function.
C     CURVE     External subroutine.
C               Inputs : S
C               Outputs: X,Y,DX,DY,D2X,D2Y
C               This is a subroutine which defines the boundary
C               curve C.  The calling sequence for it is
C                   CALL CURVE(S,X,Y,DX,DY,D2X,D2Y)
C               The variable S is the parameterization variable for
C               the curve C. (X,Y) is the point on the curve
C               corresponding to the variable S. DX,DY,D2X,D2Y are
C               the first and second derivatives of X(S) and Y(S),
C               respectively. The routine CURVE is supplied by the user.
C     B         Input.
C               The limits of the variable S in defining the curve C
C               are  0 .LE. S .LE. B
C     DPTS      Input array.
C               This is a two dimensional array specifying the points
C               at which the harmonic function U is to be evaluated,
C               with each column denoting a distinct point.  The
C               dimension statement for DPTS is
C                      DIMENSION DPTS(4,NPTS)
C               For point #J, its coordinates have the following
C               meaning:
C                         X=DPTS(1,J),    Y=DPTS(2,J)
C                         S=DPTS(3,J),    IBD=DPTS(4,J)
C               If (X,Y) is not on the boundary C, S need not be
C               set and the user must set IBD=0.
C               If (X,Y) is on the boundary, we need S to be the
C               parameter on [0,B] corresponding to (X,Y); and the
C               user must set IBD=1.
C     NPTS      Input.
C               This is the number of points in DPTS.  It must be 
C               greater than zero.
C     IBEG      Input. 
C               =0 means this is a first call on NEUMAN for this
C               particular curve C, boundary function BDYFCN, and
C               error tolerance EPS.
C               =1 means that NEUMAN has been called previously for
C               this choice of parameters, and the solution U is
C               desired at a new set of points in DPTS. For such a
C               call on NEUMAN, change only DPTS, NPTS, and IBEG.
C               DO NOT CHANGE any other input variable.
C     EPS0      Input.
C               The desired error tolerance in the solution U.
C     WORK      Input.
C               Temporary work space. This vector should have a
C               dimension NWORK of at least 10,000. If some points of
C               DPTS are close to the boundary, then the dimension of
C               WORK should be increased further. A dimension of
C               NWORK=300,000 will allow a wide variety of problems
C               to be treated very accurately.
C     NWORK     Input.
C               Dimension of WORK.  See the preceding discussion of WORK
C               and the general discussion of NWORK given below.
C     IWORK     Input.
C               Temporary work space for integer variables. It is used
C               as the pivot array in LU fcatorization by LAPACK 
C               subroutines in the subroutine INTEQN. It is also used 
C               in the subroutine EVALU for the FFT routines.
C     MWORK     Input.
C               Dimension of IWORK. 
C               MWORK .GE. (SQRT(NWORK+49) + 8 )
C     MAXFFT    Input.
C               Maximum dimension of FFT.  
C               MAXFFT .GE. (2*SQRT(NWORK + 49) - 14 )
C               is desirable to solve many problems accurately.
C     UVEC      Output.
C               This is an output vector. On exit, component #I will
C               contain the approximate value of the solution function
C               U(X(I),Y(I)) at the points (X(I),Y(I)) given in DPTS.
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
C               =-1 means EPS0 < 0.
C               =-2 means NWORK < 0 or MWORK < (SQRT(NWORK+49)+8).
C               =-3 means B  .LE. 0.
C               =-4 means IE or IBEG are out of range. 
C               =-5 means NPTS<=0, an error.
C               =-6 means MAXFFT < M_0, where M_0 is defined below.
C     
C     *** Defining NWORK, the dimension of WORK ***
C     Introduce variables MAXSYS and MAXMES.  MAXSYS represents
C     the maximum order of the linear system to be solved in
C     subroutine INTEQN; and MAXMES represents the maximum number
C     of mesh points of C at which the single layer density
C     function is to be evaluated.  MAXFFT represents the maximum
C     degree of the Fourier expansion to be produced for the 
C     approximate density function.  As examples,
C            MAXSYS=512, MAXMES=8192, MAXFFT = 1024
C     will solve many problems very accurately.  We define NWORK by
C            NWORK=MAX(MAXSYS**2+14*MAXSYS,8*MAXMES+4*MAXFFT)
C     The program assumes 
C           MAXSYS .GE. 64,  MAXMES .GE. 128,  MAXFFT .GE. 128
C     and you should set NWORK accordingly.  These defaults can be 
C     changed by re-setting the respective parameters N_0 and M_0,
C     as described below.  In particular,
C           MAXSYS .GE. 4*N_0,  MAXMES .GE. 4*M_0,  MAXFFT .GE. M_0
C
C     *** THE PROGRAM PARAMETERS K_0, M_0 and N_0 ***
C     N_0 Denotes the initial number of quadrature points used in the
C     approximation of the integral equation in the subroutine INTEQN.
C     K_0 and M_0 denote the initial degree of the  FFT and the initial 
C     number of integration nodes in the subroutine EVALU. Always set
C     M_0 and N_0 to be a power of 2, and K_0=N_0 is desirable.
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
      INTEGER IBEG,IDBG,IE,IER,IOUT,MAXFFT,MWORK,NPTS,NWORK
      DOUBLE PRECISION DPTS(4,NPTS),ERROR(NPTS),UVEC(NPTS),WORK(NWORK)
      INTEGER IWORK(MWORK)
      DOUBLE PRECISION B,BDYFCN,EPS0
      EXTERNAL BDYFCN,CURVE
C
C     LOCAL VARIABLES
      DOUBLE PRECISION LBASE(11),NBASE(11)
      DOUBLE PRECISION ZERO,TWO,SEVEN,EGHT
      DOUBLE PRECISION D1MACH,EP,EPS,FIW,FL,FN,RHOERR,U100
      INTEGER I,J,K,II,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10,
     +        ID11,IEE,IER1,IER2,K_0,LB,LD1,LD2,LD3,LD4,LD5,
     +        LD6,LD7,LD8,LD9,LD10,LD11,LU,M_0,N_0,NB,NFINAL,
     +        NFFT,NU
      EXTERNAL D1MACH,EVALU,INTEQN
      INTRINSIC FLOAT,LOG,MAX,SIGN,SQRT   
      COMMON /MACHCN/U100
      COMMON /BLKEVL/K_0,M_0
      COMMON /BLKINT/N_0
      DATA ZERO/0.0D0/,TWO/2.0D0/,SEVEN/7.D0/,EGHT/8.0D0/  
      K_0 = 16
      M_0 = 32
      N_0 = 16
C
C     TEST THE INPUT PARAMETERS.
C
      IF (EPS0 .LE. 0) THEN
          IER = -1
          RETURN
      END IF    
      IF (NWORK .LE. 0 .OR. MWORK .LT. 
     +           SQRT(FLOAT((NWORK+49)+8)) ) THEN
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
      IF (MAXFFT .LT. M_0) THEN
          IER = -6
          RETURN
      END IF
C     SET MACHINE DEPENDENT CONSTANT.  U100 IS 100 TIMES
C     MACHINE UNIT ROUND.
      U100 = 100*D1MACH(4)
      EPS = EPS0
C
      IF(IBEG .EQ. 1) GO TO 10
      EP = EPS/TWO
C     OBTAIN VALUE OF NUPPER FOR USE IN INTEQN
      FIW = NWORK
      FN = SQRT(FIW + SEVEN*SEVEN) - SEVEN
      IEE = LOG(FN)/LOG(TWO) + U100
      NU = 2**IEE
C     BREAK UP WORK INTO VECTORS AND MATRICES FOR USE IN INTEQN.
C     PRODUCE RELATIVE ADDRESSES WITHIN WORK OF RHO,X,...,D2Y,OLDRHO.
      NBASE(1) = 1
      NBASE(2) = 1 + NU*NU
      DO I = 3,11
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
      ID11 = NBASE(11)
C
      CALL INTEQN(IOUT,IDBG,IE,CURVE,B,EP,BDYFCN,NU,WORK(ID2),RHOERR,
     +            NFINAL,WORK(ID3),WORK(ID4),WORK(ID5),WORK(ID6),
     +            WORK(ID7),WORK(ID8),WORK(ID9),WORK(ID10),
     +            WORK(ID11),IWORK(16),WORK(1),IER1)
      IF (IDBG .EQ. 1) THEN
          WRITE (IOUT,9000) NFINAL,RHOERR,IER1
          WRITE (IOUT,9010)
          DO I = 1,NFINAL
            WRITE (IOUT,9020) WORK(ID3+I-1),WORK(ID4+I-1),WORK(ID2+I-1)
          END DO
      END IF
C
      IF (IER1 .EQ. 1) EP = RHOERR
C     OBTAIN LUPPER FOR USE IN EVALU.
      NFFT =MIN(4*NFINAL,MAXFFT)
      FL = (FIW-4.D0*NFFT)/EGHT
      IEE = LOG(FL)/LOG(TWO) + U100
      LU = 2**IEE
C     OBTAIN RELATIVE ADDRESSES FOR BREAKING UP WORK FOR USE
C     IN EVALU. 
      LBASE(1) = 1
      DO I = 2,9
        LBASE(I) = LBASE(I-1) + LU
      END DO 
          LBASE(10)=LBASE(9) + NFFT
          LBASE(11)=LBASE(10) + NFFT
      LD1 = LBASE(1)
      LD2 = LBASE(2)
      LD3 = LBASE(3)
      LD4 = LBASE(4)
      LD5 = LBASE(5)
      LD6 = LBASE(6)
      LD7 = LBASE(7)
      LD8 = LBASE(8)
      LD9 = LBASE(9)
      LD10 = LBASE(10)
      LD11 = LBASE(11)
C     MOVE RHO,X,Y,...,D2Y,SPD  AROUND IN WORK, LENGTHENING EACH OF
C     THEM.
      DO I = 1,8
        IF(LBASE(I)+NFINAL-1 .GE. NBASE(I+2)) THEN
            DO K = I,8
              II = I + 8 - K
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
   10 CALL EVALU(IOUT,IDBG,IE,CURVE,B,BDYFCN,NFINAL,WORK(LD1),EP,
     +           WORK(LD2),WORK(LD3),WORK(LD4),WORK(LD5),WORK(LD6),
     +           WORK(LD7),WORK(LD8),NFFT,WORK(LD9),WORK(LD10),
     +           WORK(LD11),IWORK(1),LU,DPTS,NPTS,IBEG,UVEC,ERROR,
     +           IER2)
C
      IER = MAX(IER1,IER2)
      DO I = 1,NPTS
        ERROR(I) = ERROR(I) + SIGN(RHOERR,ERROR(I))
      END DO
      IF(IER1 .EQ. 0) RETURN
      DO I = 1,NPTS
        IF (ERROR(I) .GT. ZERO) ERROR(I) = -ERROR(I)
      END DO
      RETURN

 9000 FORMAT (/,' FROM NEUMAN: RESULTS FROM SUBROUTINE INTEQN.',
     +     /,' NFINAL=',I3,5X,'RHOERR=',1P,E8.2,5X,'IER1=',I1,/)
 9010 FORMAT (6X,'X',14X,'Y',19X,'RHO')
 9020 FORMAT (1P,D12.4,D15.4,D25.12)
      END
C
      SUBROUTINE EVALU(IOUT,IDBG,IE,CURVE,B,BDYFCN,N,RHO,EPS,
     +                 X,Y,DX,DY,D2X,D2Y,SPD,NF,RFFT,OLDFFT,W,
     +                 IFAC,LU,DPTS,NP,IBEG,U,ERROR,IER)
C     ----------------
C
C     This program evaluates the single layer potential U at the
C     given points in DPTS. The input is the density function RHO,
C     and DPTS at which U is evaluated. RHO is evaluated in the
C     subroutine INTEQN. These results are stored in U, along
C     with the predicted error bound in ERROR. The desired error
C     tolerance is EPS.  If the desired error bound is not attained,
C     then the corresponding entry in ERROR is made negative. Its
C     magnitude is still an estimated error bound.
C     
C     IOUT      Input.
C               The output unit number.
C     IDBG      Input.
C               The debugging parameter.
C               = 0 produces a standard ouput.
C               = 1 produces a debugging parameter.
C     IE        Input.
C               =0 for the interior Neumann problem.
C               =1 for the exterior Neumann problem.
C     CURVE     External subroutine.
C               This program defines the curve of the boundary C.
C     B         Input.
C               For the parameterization of C defined in CURVE,
C               the parameterization interval is [0,B].
C     BDYFCN    Input.
C               This program defines the Neumann data on the
C               boundary.
C     N         Input.
C               This is NFINAL as output from the subroutine INTEQ.
C     RHO       Input.
C               An array which contains the value of the single
C               layer density function defining U.
C     EPS       Input.
C               The user-supplied absolute error tolerance.
C     X,Y       Input.
C               Two arrays containing a sequence of points
C               (X(I),Y(I)) produced by calling the subroutine
C               NEWCUR.  They correspond to an even subdivision
C               of the parameterization interval [0,B].
C     DX,DY     Input.
C               The derivative values corresponding to the points
C               given in X,Y.
C     D2X,D2Y   Input.
C               The second derivative values corresponding to the
C               points given in X,Y.
C     SPD       Input.
C               An array which contains SQRT(DX*DX+DY*DY).
C     NF        Input.
C               The size of the array RFFT.  NF .GE. 2*N.
C     RFFT      Output.
C               An array which contains the discrete Fourier
C               coefficients of RHO.
C     W         Input.
C               A work array which is needed in the fft subroutine.
C     IFAC      Input.
C               An integer array which is needed in the fft subroutine.
C     LU        Input.
C               The upper bound on the size of the arrays X,Y,DX,DY,
C               D2X,D2Y,SPD,RHO.
C     DPTS      Input.
C               This is a two dimensional array which supplies the
C               points at which U is to be evaluated.
C     NP        Input.
C               This is the number of points in DPTS.
C     IBEG      Input.
C               =0 means this is a first call on NEUMAN for this
C               particular curve C, boundary function BDYFCN, and
C               error tolerance EPS.
C               =1 means that NEUMAN has been called previously for
C               this choice of parameters, and the solution U is
C               desired at a new set of points in DPTS. For such a
C               call on NEUMAN, change only DPTS, NPTS, and IBEG.
C               Do not change any other input variable.
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
      INTEGER IBEG,IE,IER,IOUT,LU,N,NF,NP 
      DOUBLE PRECISION D2X(LU),D2Y(LU),DPTS(4,NP),DX(LU),DY(LU),
     +                 ERROR(NP),OLDFFT(NF),RFFT(NF),RHO(LU),SPD(LU),
     +                 U(NP),W(2*NF),X(LU),Y(LU)
      INTEGER IFAC(15)
      DOUBLE PRECISION D1MACH,B,BDYFCN,EPS
      EXTERNAL BDYFCN,CURVE
C
C     LOCAL VARIABLES
      DOUBLE PRECISION BDCON,DIFF,DIST,ERR,ERRFFT,FCNK,FCNM,FCNU,FUEVAL,
     +                 H,HH,OLDIFF,OLDU,ORHO,PASTRT,PI,PX,PY,
     +                 R,RATIO,RATE,RTLOW,RTUP,S,SINE,COEF,COSN,SING,
     +                 SUM,T1,T2,THETA,TT,U100
      DOUBLE PRECISION ZERO,ONE,TWO,FOUR
      INTEGER I,IBD,J,K,KH,KK,KSTEP,K_0,L,LB,LD,LDM1,LOOP,M_0,
     +        IDBG, JH, JLOOP
      EXTERNAL FUCOEF,NEWCUR,FUEVAL
      INTRINSIC ABS,LOG,MAX,MIN,SIN,SQRT
      COMMON /MACHCN/U100
      COMMON /BLKEVL/K_0,M_0
      DATA ZERO/0.0D0/,ONE/1.0D0/,FOUR/4.D0/,RTUP/0.5D0/,
     +     RTLOW/.1D0/,TWO/2.0D0/
C
C     K_0 IS THE INITIAL DEGREE OF THE FOURIER SERIES EXPANSION TO BE 
C     USED IN THE APPROXIMATION OF THE DENSITY FUNCTION. 
C     M_0 DENOTES THE INITIAL NUMBER OF INTEGRATION NODES TO BE USED
C     IN THE EVALUATION OF THE SINGLE LAYER POTENTIAL.   THESE OPERATIONS
C     WILL ALSO BE PERFORMED WITH THE PARAMETER 2*K_0 and 2*M_0, SO THAT 
C     THE VALUE OF NWORK IN THE CALLING PROGRAM NEEDS TO BE SET 
C     ACCORDINGLY.  ALWAYS SET K_0 AND M_0 TO BE A POWER OF 2.
C
C     STAGE 1: EVALUATE THE FOURIER COEFFICIENTS OF RHO. 
C     INITIALIZE. LD IS THE RUNNING DIMENSION OF RHO,X,...,D2Y.  LB 
C     IS THE RUNNING DIMENSION OF THE FOURIER COEFFICIENT VECTOR RFFT.
C     LB CAN BE AS BIG AS NF, WHERE NF IS INPUT FROM NEUMAN AND IS 
C     SET TO BE MIN(4*NFINAL,MAXFFT) WITH MAXFFT SET BY THE USER ON 
C     CALLING SUBROUTINE NEUMAN.
C
C     DATA 'PI'
      PI = FOUR*ATAN(ONE)
      IF (IDBG .EQ. 1) THEN
          WRITE(IOUT,9000) 
          WRITE(IOUT,9010) 
      END IF       
      IF (IBEG .GT. 0) GO TO 40
      RATE = RTUP
      LD = N
      LB = K_0
      LOOP = 1
      DO I = 1,LD
        RHO(I) = RHO(I)*SPD(I)
      END DO
      PASTRT=RTUP 
C     EVALUATE THE FOURIER COEFFICIENTS OF RHO*SPD FOR LATER USE.
C     FIRST ASSIGN THE INTERPOLATING POINTS TO RFFT. AFTER
C     CALLING FUCOEF, RFFT CONTAINS THE THE FOURIER COEFFICIENTS.
   10 KSTEP = LD/LB
      DO J = 1, LB-1
        JH = J*KSTEP            
        RFFT(J+1) = RHO(JH)
      END DO
      RFFT(1) = RHO(LD)
      CALL FUCOEF(LB,RFFT,W,IFAC) 
      DIFF = ZERO
      IF (LOOP .EQ. 1)  GO TO 20     
C     ESTIMATE THE ERROR IN THE FFT     
      DO I = 1, LB/2
        COEF = ONE/I
        IF (I .LT. LB/4)  THEN
            COSN = RFFT(2*I) -OLDFFT(2*I)
            SINE = RFFT(2*I+1) -OLDFFT(2*I+1)
        ELSE IF (I. EQ. LB/4) THEN
            COSN = RFFT(2*I) - OLDFFT(2*I)/2
            SINE = RFFT(2*I+1)
        ELSE IF (I .LT. LB/2) THEN
            COSN = RFFT(2*I)
            SINE = RFFT(2*I+1)
        ELSE 
            COSN = RFFT(LB)/2
            SINE = ZERO
        END IF
        DIFF = DIFF + TWO*COEF*SQRT(COSN*COSN + SINE*SINE)
      END DO 
      IF (LOOP .EQ. 2) GO TO 20
C     UPDATE THE RATE OF CONVERGENCE OF THE FFT
      RATE = MAX(PASTRT,RTLOW,MIN(RTUP,ABS(DIFF/OLDIFF)))
      PASTRT = MIN(RTUP,ABS(DIFF/OLDIFF))
      ERRFFT = RATE/(ONE-RATE)*DIFF
      IF (IDBG .EQ. 1) THEN
          WRITE(IOUT, 9020) LB, LD,DIFF, ERRFFT, RATE
      END IF 
      IF (ERRFFT .LT. EPS/20) GO TO 30
  20  OLDIFF = MAX(D1MACH(4), DIFF)
      DO I = 1, LB
        OLDFFT(I) = RFFT(I) 
      END DO 
      LB = 2*LB
      LOOP = LOOP + 1 
      IF ((LB .LE. LD) .AND. (LB .LE. NF))  GO TO 10
C     AN INSUFFICIENT NUMBER OF VALUES IN RHO. THUS VALUES OF RHO ON
C     A FINER MESH MUST BE CREATED.
      LD = 2*LD
      IF (LD .GT. NF) THEN
          LD = LD/2
          LB = LB/2
          LOOP = LOOP - 1 
          GO TO 30
      END IF
C     THERE IS SUFFICIENT SPACE IN RHO,X,Y,...,D2Y FOR AN
C     INCREASED SUB-DIVISION OF (0,B).  MOVE OLD VALUES OF
C     RHO,...,D2Y TO MAKE ROOM FOR NEW VALUES.
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
        SPD(K) = SPD(KH)
      END DO
C     PRODUCE NEW CURVE PARAMETERS FOR FINER SUBDIVISION.
      H = B/LD
      LDM1 = LD - 1
      DO J = 1, LDM1, 2 
        S = J*H
        CALL NEWCUR(IE,CURVE,S,X(J),Y(J),DX(J),DY(J),D2X(J),
     +                  D2Y(J))
        SPD(J) = SQRT(DX(J)*DX(J)+DY(J)*DY(J))
      END DO
C     PRODUCE NEW VALUES OF RHO.
      HH = B/N
      KSTEP = LD/N
      DO J = 1,LDM1,2
        SUM = ZERO
        DO K = KSTEP,LD,KSTEP
          T1 = X(J) - X(K)
          T2 = Y(J) - Y(K)
          FCNK = (DY(J)*T1-DX(J)*T2)/ ((T1*T1+T2*T2)*SPD(J))
          SUM = SUM + FCNK*RHO(K)
        END DO 
        IF(IE .EQ. 1) THEN
C           FOR THE EXTERIOR PROBLEM
            BDCON = BDYFCN(J*H)
        ELSE
C           FOR THE INTERIOR PROBLEM
            R = X(J)**2 + Y(J)**2
            BDCON = -ONE/R*BDYFCN(J*H)
        END IF
        ORHO = - (BDCON+HH*SUM)/PI
        RHO(J) = ORHO*SPD(J)
      END DO
      GO TO 10
   30 CONTINUE
      IF (IDBG .EQ. 1) THEN
          WRITE(IOUT,9030)
          WRITE(IOUT,9040) RFFT(1) 
          WRITE(IOUT,9050) (I, RFFT(2*I), RFFT(2*I +1),I=1,LB/2-1)
          WRITE(IOUT,9060) LB/2, RFFT(LB)
          WRITE(IOUT,9070)
      END IF 
C
C     STAGE2: BEGIN LOOP TO EVALUATE U AT POINTS P IN DPTS.
C       
   40 IER = 0
      DO I = 1,NP
        PX = DPTS(1,I)
        PY = DPTS(2,I)
        IBD = DPTS(4,I)
        IF (IBD.NE.0) THETA = DPTS(3,I)
C       IF IT IS AN INTERIOR PROBLEM.  CHANGE THE DPTS BY USING
C       KELVIN TRANSFORMATION. IF DPTS=(0,0), THEN ASSIGN SOME BIG
C       NUMBERS FOR (PX,PY).
        IF(IE .EQ. 0) THEN
            R = PX*PX + PY*PY
            IF(R .LT. U100) THEN
                R = U100
                PX = 1.D0/R
                PY = 1.D0/R
            ELSE
                PX = PX/R
                PY = PY/R
            END IF
        END IF
        IF (IDBG .EQ. 1) THEN
            WRITE (IOUT,9080) I,DPTS(1,I),DPTS(2,I),THETA,IBD,LD
        END IF
C
C       NOW BEGIN EVALUATION OF U(PX,PY) USING NUMERICAL INTEGRATION.
C       INITIALIZE, AND BEGIN WITH M_0 NODES.
        RATE = RTUP
        L = M_0
        LOOP = 1
        JLOOP = 1
        PASTRT = RTUP
C       CALCULATE NUMERICAL INTEGRAL WITH L SUBDIVISIONS OF (0,B).
   50   SUM = ZERO
        KSTEP = LD/L
        H = B/L
        DO K = KSTEP,LD,KSTEP
          KK = K/KSTEP
          S = KK*H
          T1 = X(K) - PX
          T2 = Y(K) - PY
          DIST = SQRT(T1*T1+T2*T2)
C         IF DPTS IS A BOUNDARY POINT(IBD=1), THE INTEGRAND HAS A
C         SINGULARITY AT S=THETA. DIVIDE THE INTEGRAND INTO TWO PARTS.
C         FCNM IS THE SMOOTH PART, AND SING IS THE SINGULAR PART.
          IF(IBD .EQ. 1) THEN
              IF(THETA .EQ. S) THEN
                  TT = SQRT(DX(K)*DX(K)+DY(K)*DY(K))
                  FCNM = LOG(TT)
              ELSE
                  TT = DIST/ABS(TWO*SIN((THETA-S)/TWO))
                  FCNM = LOG(TT)
              END IF
          ELSE
              FCNM = LOG(DIST)
          END IF
          FCNM = -FCNM*RHO(K)
          SUM = SUM + FCNM
        END DO
        IF(IBD .EQ. 1) THEN
            SING = PI*FUEVAL(B,THETA,RFFT,LB)
        ELSE
            SING = ZERO
        END IF
        FCNU = H*SUM + SING
        IF (IDBG. EQ. 1) THEN
            WRITE (IOUT,9090) FCNU,L,LD
        END IF
        IF(LOOP .EQ. 1) GO TO 80
C       ESTIMATE ERROR IN FCNU.
        DIFF = ABS(FCNU-OLDU)
        IF(LOOP .EQ. 2) GO TO 60
C       UPDATE RATE OF CONVERGENCE OF NUMERICAL INTEGRATION.
        RATIO = ABS(DIFF/OLDIFF)
C       RATE=MAX(PASTRT,RTLOW,MIN(RTUP,RATIO))
        PASTRT = MIN(RTUP,RATIO)
   60   ERR = (RATE/ (ONE-RATE))*DIFF
        IF (IDBG. EQ. 1) THEN
            WRITE (IOUT,9100) DIFF,ERR,RATE  
        END IF 
        IF (ERR .GT. EPS) GO TO 70        
C       THE VALUE OF FCNU IS SUFFICIENTLY ACCURATE.
        U(I) = FCNU
        ERROR(I) = ERR
        IF (IBD .EQ. 1) THEN 
            ERROR(I) = ERR + ERRFFT
            IF (ERROR(I) .GT. EPS) THEN
                ERROR(I) = -ERROR(I)
                IER = 1
            END IF
        END IF 
        GO TO 100
C       FCNU IS NOT SUFFICIENTLY ACCURATE.
C       RE-INITIALIZE FOR ANOTHER NUMERICAL INTEGRATION.
   70   OLDIFF = MAX(DIFF,D1MACH(4))
   80   OLDU = FCNU
        LOOP = LOOP + 1
        L = 2*L
        IF(L .LE. LD) GO TO 50
C       NOT SUFFICIENT VALUES IN RHO. THUS VALUES OF RHO ON A FINER
C       MESH MUST BE CREATED.
        LD = 2*LD
        IF(LD .GT. LU) GO TO 90
C       THERE IS SUFFICIENT SPACE IN RHO,X,Y,...,D2Y FOR AN
C       INCREASED SUB-DIVISION OF (0,B).  MOVE OLD VALUES OF
C       RHO,...,D2Y TO MAKE ROOM FOR NEW VALUES.
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
          SPD(K) = SPD(KH)
        END DO
C       PRODUCE NEW CURVE PARAMETERS FOR FINER SUBDIVISION.
        H = B/LD
        LDM1 = LD - 1
        DO J = 1,LDM1,2
          S = J*H
          CALL NEWCUR(IE,CURVE,S,X(J),Y(J),DX(J),DY(J),D2X(J),D2Y(J))
          SPD(J) = SQRT(DX(J)*DX(J)+DY(J)*DY(J))
        END DO
C       PRODUCE NEW VALUES OF RHO.
        HH = B/N
        KSTEP = LD/N
        DO J = 1,LDM1,2
          SUM = ZERO
          DO K = KSTEP,LD,KSTEP
            T1 = X(J) - X(K)
            T2 = Y(J) - Y(K)
            FCNK = (DY(J)*T1-DX(J)*T2)/ ((T1*T1+T2*T2)*SPD(J))
            SUM = SUM + FCNK*RHO(K)
          END DO  
          IF(IE .EQ. 1) THEN
C             FOR THE EXTERIOR PROBLEM
              BDCON = BDYFCN(J*H)
          ELSE
C             FOR THE INTERIOR PROBLEM
              R = X(J)**2 + Y(J)**2
              BDCON = -ONE/R*BDYFCN(J*H)
          END IF
          ORHO = - (BDCON+HH*SUM)/PI
          RHO(J) = ORHO*SPD(J)
        END DO
        GO TO 50
C       THE UPPER LIMITS FOR RHO,X,Y,...,D2Y HAVE BEEN REACHED.
C       MARK ERROR BOUND ACCORDINGLY AND CONTINUE ONTO NEXT POINT P. 
   90   ERROR(I) = -ERR
        U(I) = FCNU
        IER = 1
        LD = LD/2
        ERROR(I) = - ERR
        IF (IBD .EQ. 1) ERROR(I) = - (ERR + ERRFFT)
  100 END DO
      RETURN
 9000 FORMAT (/,' FROM SUBROUTINE EVALU.',/)
 9010 FORMAT (/,' STAGE 1: FOURIER COEFFICIENTS.',/)
 9020 FORMAT (' LB =', I4, 3X, 'LD = ', I4, 3X, 'DIFF =', 1PD9.2, 3X,
     *        'ERROR =', D9.2, 3X, 'RATE =', 1D9.2)
 9030 FORMAT (/,10X,'COEFF. OF COSINE',11X,'COEFF. OF SINE') 
 9040 FORMAT (/,2X,' 0',3X,E20.12)
 9050 FORMAT (1X, I3,3X,E20.12,5X,E20.12) 
 9060 FORMAT (1X, I3,3X,E20.12)  
 9070 FORMAT (/,'STAGE2: EVALUATE U(P)', /) 
 9080 FORMAT (/,' I=',I2,3X,'PX=',1PD11.4,3X,'PY=',D11.4,3X,
     +   'THETA=',D11.4,3X,'IBD=',I1,3X,'LD=',I5)
 9090 FORMAT (' NUM INT =',1P,E20.12,5X,'L=',I6,5X,'LD=',I6)
 9100 FORMAT (5X,'DIFF=',1PE9.2,5X,' ERROR=',D9.2,5X,'RATE=',D11.4)
      END
C
      SUBROUTINE FUCOEF(N,R,W,IFAC)
C     -----------------
C
C     This generates approximate Fourier coefficients for RHO,
C     with the calculations done with an FFT program.
C
C     The input is R(1),R(2),...R(N), which contain RHO at N evenly
C     space points in the interval [0,2*pi], with R(i) = i*2*pi/N.
C     The output is R(1)...R(2*N), which contains the Fourier
C     coefficients.  The array W is workspace for the FFT.
C
C     Inputs: N, R, W, IFAC
C     Outputs: R
C
      INTEGER J,N
      DOUBLE PRECISION R(N),W(2*N)
      INTEGER IFAC(15)
      EXTERNAL DRFFTF,DRFFTI
C
      CALL DRFFTI(N,W,IFAC)
      CALL DRFFTF(N,R,W,IFAC)
      DO 20 J = 1,N
        R(J) = R(J)/FLOAT(N)
   20   CONTINUE
      RETURN
      END

      DOUBLE PRECISION  FUNCTION FUEVAL(B,S,R,N)
C     ---------------------------------
C
C     Evaluate the singular part of the integral for the single
C     layer potential, evaluated at points on the boundary curve C.
C     The input are the Fourier coefficients R(1)...R(N) which is
C     evaluated in subroutine FUCOEF.
C
C     Inputs: B, S, R, N
C     Output: FUEVAL
C
      INTEGER I,M,N
      DOUBLE PRECISION R(N)
      DOUBLE PRECISION B,COEF,COSN,S,SINE,PI
      INTRINSIC COS,SIN,ATAN
C
      PI = 4.D0*ATAN(1.D0)
      M = N/2
      SINE = 0.D0
      COSN = 0.D0
      DO 10 I = 1,M-1
        COEF = 1.D0/I
        COSN = COSN + COEF*R(2*I)*COS(I*S*2*PI/B)
        SINE = SINE + COEF*R(2*I+1)*SIN(I*S*2*PI/B)
   10   CONTINUE
      FUEVAL= 2.D0*(-SINE + COSN) + R(N)*COS(M*S*2*PI/B)/M
      RETURN
      END
C
      SUBROUTINE INTEQN(IOUT,IDBG,IE,CURVE,B,EPS,BDYFCN,NUPPER,RHO,
     +                  ERROR,NFINAL,X,Y,DX,DY,D2X,D2Y,SPD,OLDRHO,
     +                  WORK,IWORK,KERMAT,IER)
C     -----------------
C
C     This program solves the second kind boundary integral equation
C     which arises from solving the exterior Neumann problem as a
C     single layer potential.
C
C     The integral equation is solved using Nystrom's method with
C     the rectangular rule as the quadrature rule.   The resulting
C     linear system is solved directly using LINPACK routines.
C
C     The output is the single layer density function RHO.  This
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
C     calling NEUMAN.
C
C     IOUT     Input.
C              The output unit number, to the file NEUMAN.ANS
C     IDBG     Input. 
C              The debugging parameter.
C               = 0 produces a standard ouput.
C               = 1 produces a debugging parameter.
C     IE       Input.
C              =0 for the interior Neumann problem.
C              =1 for the exterior Neumann problem.
C     CURVE    External subroutine.
C              This program defines the curve of the boundary C.
C     B        Input.
C              For the parameterization of C defined in CURVE,
C              the parameterization interval is [0,B].
C     EPS      Input.
C              The user-supplied absolute error tolerance.
C     BDYFCN   External function.
C              This program defines the Neumann data on the boundary.
C     NUPPER   Input.
C              This is the upper bound for the size of linear
C              system that can be constructed and solved.
C     RHO      Output.
C              An array which contains the value of the single
C              layer density function defining U.
C     ERROR    Output.
C              This is the predicted error estimate for RHO.
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
C     SPD      Output.
C              An array which contains SQRT(DX*DX+DY*DY).
C     OLDRHO   Output.
C              An array containing the preceding value of RHO,
C              also produced in this program.
C     WORK     Output.
C              This is a work array used in the LAPACK routine.
C     IWORK    Output.
C              This is an array for pivoting used in the LAPACK routine.
C     KERMAT   Output.
C              This is array contains the linear system associated
C              with the Nystrom method.
C     IER      Output.
C              =0 means the program was completed satisfactorily.
C              =1 means the desired error uniform error bound EPS
C              for the solution RHO was not attained.
C
      INTEGER IE,IER,INFO,IOUT,NFINAL,NUPPER
      DOUBLE PRECISION D2X(NUPPER),D2Y(NUPPER),DX(NUPPER),DY(NUPPER),
     +                 KERMAT(NUPPER,NUPPER),OLDRHO(NUPPER),
     +                 RHO(NUPPER),SPD(NUPPER),
     +                 X(NUPPER),Y(NUPPER)
      DOUBLE PRECISION D1MACH,B,BDYFCN,EPS,ERROR
      EXTERNAL BDYFCN,CURVE,NEWCUR
C
C     VARIABLES FOR THE LAPACK ROUTINES.
      DOUBLE PRECISION  WORK(4*NUPPER)
      INTEGER IWORK(NUPPER), NRHS, IDBG
C      
C     LOCAL VARIABLES
      DOUBLE PRECISION DIFF,DIST,H,OLDIFF,PI,R,RATE,RCOND,
     +                 RTLOW,RTUP,SUM,SUMAX,T1,T2
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      INTEGER I,J,JH,N,NM1,N_0
      INTRINSIC ABS,LOG,MAX,MIN,SQRT
      COMMON /BLKINT/N_0
      DATA RTLOW/.1D0/,RTUP/.5D0/,ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,
     +     FOUR/4.D0/, NRHS/1/
C
C     THE PARAMETER N_0 GIVES THE INITIAL NUMBER OF QUADRATURE POINTS
C     USED IN THE APPROXIMATION OF THE INTEGRAL EQUATION.  THE EQUATION
C     WILL ALSO BE SOLVED WITH 2*N_0 NODES.  IN THE PROGRAM CALLING
C     NEUMAN, THE PARAMETER NWORK SHOULD BE SET ACCORDING.  ALWAYS SET 
C     N_0 TO BE A POWER OF 2.
C
C     INITIAL CASE, N=N_0. INITIALIZE PARAMETERS.
      IF (IDBG .EQ. 1) THEN
          WRITE (IOUT,9000)
      END IF
C     DATA 'PI'
      PI = FOUR*ATAN(ONE) 
      N = N_0
      RATE = RTUP
C     DEFINE STEPSIZE AND POINTS ON CURVE.
      H = B/N
      DO I = 1,N
        CALL NEWCUR(IE,CURVE,I*H,X(I),Y(I),DX(I),DY(I),D2X(I),
     +                D2Y(I))
        SPD(I) = SQRT(DX(I)*DX(I)+DY(I)*DY(I))
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
        SPD(J) = SPD(JH)
      END DO
      NM1 = N - 1
      DO I = 1,NM1,2
        CALL NEWCUR(IE,CURVE,I*H,X(I),Y(I),DX(I),DY(I),D2X(I),
     +                D2Y(I))
        SPD(I) = SQRT(DX(I)*DX(I)+DY(I)*DY(I))
      END DO
C     SET UP MATRIX EQUATION, AND EVALUATE THE MAXIMUM NORM OF 
C     SINGLE LAYER POTENTIAL TO BE USED IN ERROR ESTIMATE. 
C     HERE 'SUM' REPRESENTS THE NORM.
   20 SUMAX = ZERO
      DO I = 1,N
        IF(IE .EQ. 1) THEN
C           BOUNDARY CONDITION FOR THE EXTERIOR PROBLEM
            RHO(I) = BDYFCN(I*H)
        ELSE
C           BOUNDARY CONDITION FOR THE INTERIOR PROBLEM
            R = X(I)**2 + Y(I)**2
            RHO(I) = -ONE/R*BDYFCN(I*H)
        END IF
        SUM = ZERO
        DO J = 1,N
          IF(I .EQ. J) THEN
C             DEFINE KERNEL FOR T(I) = T(J).
              T1 = DX(I)
              T2 = DY(I)
              DIST = T1*T1 + T2*T2
              KERMAT(I,I) = -PI - H* (T1*D2Y(I)-T2*D2X(I))/
     +                          (TWO*DIST)
          ELSE
C             DEFINE KERNEL FOR T(I) .NE. T(J)
              T1 = X(I) - X(J)
              T2 = Y(I) - Y(J)
              DIST = T1*T1 + T2*T2
              KERMAT(I,J) = -H* (DY(I)*T1-DX(I)*T2)*SPD(J)/
     +                          (SPD(I)*DIST)
              SUM = SUM + ABS(LOG(DIST))*SPD(J)
          END IF
          SUM = H*SUM
        END DO
        SUMAX = MAX(SUMAX,SUM)
      END DO
C     SOLVE LINEAR SYSTEM.
      CALL DGESV(N,NRHS,KERMAT,NUPPER,IWORK,RHO,NUPPER,INFO)
      CALL DGECON('I',N,KERMAT,NUPPER,SUMAX,RCOND,WORK,
     +            IWORK,INFO)
      IF(N .EQ. N_0) GO TO 40
C     CALCULATE NORM OF RHO-OLDRHO.
      DIFF = ZERO
      DO I = 2,N,2
        DIFF = MAX(DIFF,ABS(RHO(I)-OLDRHO(I/2)))
      END DO
      IF(N .EQ. 2*N_0) GO TO 30
C     MEASURE RATE OF CONVERGENCE.
      RATE = DIFF/OLDIFF
      RATE = MAX(RTLOW,MIN(RATE,RTUP))
C     ESTIMATE ERROR IN RHO.
   30 ERROR = (ONE/RCOND)*SUM*DIFF*RATE/ (ONE-RATE)
      IF (IDBG .EQ. 1) THEN
          WRITE(IOUT,9010) N,DIFF,ERROR
      END IF
      IF(ERROR .LE. EPS) THEN
C         EXIT FOR SUCCESSFUL RETURN.
          NFINAL = N
          IER = 0
          RETURN
      ELSE IF (2*N .GT. NUPPER) THEN
C     EXIT FOR UNSUCCESSFUL RETURN.
          IER = 1
          NFINAL = N
          RETURN
      END IF
C     PREPARE FOR ANOTHER LOOP ON N.
      OLDIFF = MAX(DIFF,D1MACH(4))
   40 DO I = 1,N
        OLDRHO(I) = RHO(I)
      END DO
      N = 2*N
      GO TO 10
 9000 FORMAT (/,' FROM SUBROUTINE INTEQN')
 9010 FORMAT (' N=',I3,5X,'DIFF=',1PD8.2,5X,'ERROR=',D8.2)
      END
C
      SUBROUTINE NEWCUR(IE,CURVE,S,X,Y,DX,DY,D2X,D2Y)
C     -----------------
C
C     Define a new curve if the given problem is an interior
C     Neumann problem(IE=0), using the Kelvin transformation.
C
C     Inputs: IE, S 
C     External subroutine: CURVE
C     Outputs: X, Y, X, Y, DX, DY, D2X, D2Y
C
      DOUBLE PRECISION D2X,D2Y,DX,DY,S,X,Y
      INTEGER IE
      EXTERNAL CURVE, KVTRNF
      DOUBLE PRECISION D2TX,D2TY,DTX,DTY,TX,TY
      CALL CURVE(S,X,Y,DX,DY,D2X,D2Y)
      IF (IE .EQ. 0) THEN
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
C
      SUBROUTINE KVTRNF(X,Y,DX,DY,D2X,D2Y,TX,TDX,TD2X)
C     -----------------
C
C     Define the Kelvin transformation.  If we disposition
C     X's and Y's, we have TY, TDY, TD2Y which are the
C     transformed values of Y, DY, D2Y.
C
C     Inputs: X, Y, DX, DY, D2X, D2Y
C     Outputs: TX, TDX, TD2X
C
      DOUBLE PRECISION D2X,D2Y,DX,DY,TD2X,TDX,TX,X,Y
      DOUBLE PRECISION DIST,T1,T2,XS,YS
      XS = X*X
      YS = Y*Y
      DIST = XS + YS
      TX = X/DIST
      TDX = (DX* (YS-XS)-2.D0*X*Y*DY)/ (DIST*DIST)
      T1 = D2X* (YS*YS-XS*XS) - 2.D0*X* (XS+YS)* (DY*DY+DX*DX+Y*D2Y)
      T2 = -4.D0* (X*DX+Y*DY)* (DX* (YS-XS)-2.D0*X*Y*DY)
      TD2X = (T1+T2)/DIST**3
      RETURN
      END
