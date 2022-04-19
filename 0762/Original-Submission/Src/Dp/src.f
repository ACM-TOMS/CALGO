      SUBROUTINE lldrlf(ww,dfn,dfd,case,l,dldw,d2ldw2)
C**********************************************************************
C
C     SUBROUTINE LLDRLF(WW,DFN,DFD,CASE,L,DLDW,D2LDW2)
C          Log-Likelihood and DeRivatives for Log-F Models
C
C
C                              Function
C
C
C      Let F(Z|DFN,DFD) be the cumulative F distribution with degrees
C      of freedom DFN and DFD. Let W have a Log-F distribtuion with
C      the same degrees of freedom. That is, exp(W) is distributed
C      as F. Let f(exp(W)|DFN,DFD) be the density corresponding to
C      F(exp(W)|DFN,DFD), f(exp(W)|DFN,DFD) = dF(exp(W)|DFN,DFD)/dw.
C
C      This subroutine calculates one of the following three values
C            log(f(exp(W)|DFN,DFD))
C            log(F(exp(W)|DFN,DFD))
C            log(1 - F(exp(W)|DFN,DFD))
C      and its first two derivatives with respect to W.
C
C
C
C                              Arguments
C
C
C      WW --> Argument at which log(f(exp(W)|DFN,DFD)),
C             log(F(exp(W)|DFN,DFD)) or log(1 - F(exp(W)|DFN,DFD))
c             is to be calculated
C                               DOUBLE PRECISION WW
C
C      DFN --> Numerator degrees of freedom
C                               DOUBLE PRECISION DFN
C
C      DFD --> Denominator degrees of freedom
C                               DOUBLE PRECISION DFD
C
C      CASE --> Indicates which function of W is to be calculated
C                 1 : log(f(exp(W)|DFN,DFD))
C                 2 : log(F(exp(W)|DFN,DFD))
C                 3 : log(1 - F(exp(W)|DFN,DFD))
C                               INTEGER CASE
C
C      L <-- log(f(exp(W)|DFN,DFD)), log(F(exp(W)|DFN,DFD)) or
C            log(1 - F(exp(W)|DFN,DFD)) depending on CASE.
C                               DOUBLE PRECISION L
C
C      DLDW <-- Derivative of L with respect to W
C                               DOUBLE PRECISION DLDW
C
C      D2LDW2 <-- Second derivative of L with respec to W (twice)
C                               DOUBLE PRECISION D2LDW2
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION center
      PARAMETER (center=0.1D0)
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0D0)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=1.0D-8)
      DOUBLE PRECISION thp
      PARAMETER (thp=0.03D0)
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
      DOUBLE PRECISION pt7
      PARAMETER (pt7=0.7D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
      DOUBLE PRECISION two
      PARAMETER (two=2.0D0)
      DOUBLE PRECISION fiften
      PARAMETER (fiften=1500.0D0)
      DOUBLE PRECISION forty
      PARAMETER (forty=40.0D0)
      DOUBLE PRECISION hundrd
      PARAMETER (hundrd=100.0D0)
      DOUBLE PRECISION maxx
      PARAMETER (maxx=0.74D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ww,dfn,dfd,l,dldw,d2ldw2
      INTEGER case
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,aswp,alopk,b,bswp,w,expwpk,k,logk,opk,opkemw,
     +                 prd,remain,rat,wswp,xld,xd1,tau,lbda,logz,logy,
     +                 logx,z,lcum,x,y,lambda,thpab,eps3,eps,epsold,
     +                 eps5,temp
      LOGICAL qcentr,qswtch
C     ..
C     .. External Functions ..
      DOUBLE PRECISION cenlf,dbetrm,dexpm1,dln1pe,dln1px,spmpar,dexp1
      EXTERNAL cenlf,dbetrm,dexpm1,dln1pe,dln1px,spmpar,dexp1
C
C     .. External Subroutines ..
      EXTERNAL ltlf,dlasym,dlfrac,dser,dextr
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,log,max,min,sqrt
C     ..
C     .. Save statement ..
      SAVE eps,epsold,eps3,eps5
C     ..
C     .. Data statements ..
      DATA eps/-1.0D0/,epsold/-2.0D0/
C     ..
C     .. Executable Statements ..
C
C     EPS, EPS3 and EPS5 are machine dependent and only need
C     to be calculated once per execution of a program calling llpder.
C
      IF (eps.NE.epsold) THEN
          eps = spmpar(1)
          epsold = eps
          eps3 = hundrd*eps
          eps5 = fiften*eps
      END IF
C
C     The following code calculates log(f(exp(W)|DFN,DFD)) or
C     log(1 - F(exp(W)|DFN,DFD)). To get log(F(exp(W)|DFN,DFD))
C     We take advantage of the fact that
C        log(F(exp(W)|DFN,DFD)) = log(1 - F(exp(-W)|DFD,DFN))
C
C       d log(F(exp(W)|DFN,DFD)) =    d log(1 - F(exp(-W)|DFD,DFN))
C       ________________________   - ______________________________
C                  dw                              dw
C     and
C       d^2 log(F(exp(W)|DFN,DFD)) =   d^2 log(1 - F(exp(-W)|DFD,DFN))
C       ________________________       ______________________________
C                  dw^2                            dw^2
C
C
      w = ww
C
C     Compute log(f(exp(W)|DFN,DFD)) and its first derivative
C     Note: These computations are needed to find the derivatives of
C           log(1 - F(exp(W)|DFN,DFD))
C
C
C     Calculate some functions of a and b
C
      a = half*dfn
      b = half*dfd
      IF (case.EQ.2) THEN
          temp = b
          b = a
          a = temp
          w = -w
      END IF

      aswp = max(a,b)
      bswp = min(a,b)
      qswtch = b .GT. a
C
C     Compute quantities for calculation of f(w) and df(w)/dw
C
      k = bswp/aswp
      logk = log(k)
      opk = one + k
      prd = bswp*opk
      alopk = dln1px(k)
      remain = dbetrm(a,b)
      tau = sqrt(one/ (one/a+one/b))
C
C     Reverse sign of W if B greater than A
C
      IF (qswtch) THEN
          wswp = -w

      ELSE
          wswp = w
      END IF

      qcentr = abs(wswp) .LE. center
C
C     Compute log(f(exp(W)|DFN,DFD))
C
      IF (qcentr) THEN
          xld = cenlf(wswp,k)

      ELSE
          xld = -k*wswp + opk* (alopk-dln1pe(logk-wswp))
      END IF

      xld = aswp*xld - remain
C
C     Compute first derivative of log(f(exp(W)|DFN,DFD))
C
      IF (wswp.GT.zero) THEN
          opkemw = one + k*exp(-wswp)

      ELSE
          expwpk = exp(wswp) + k
      END IF

      IF (wswp.GT.zero) THEN
C        XD1 = BSWP * (EXP(-WSWP) - ONE) / ( ONE + K*EXP(-WSWP) )
          xd1 = bswp*dexpm1(-wswp)/opkemw

      ELSE
C        XD1 = BSWP * ( ONE - EXP(WSWP) ) / ( EXP(WSWP) + K )
          xd1 = -bswp*dexpm1(wswp)/expwpk
      END IF

      IF (qswtch) xd1 = -xd1
C
C     Calculate log(f(exp(W)|DFN,DFD)) and its first and
C     its first and second derivatives if needed.
C
      IF (case.EQ.1) THEN
          l = xld
          dldw = xd1
C
C        Compute the second derivative of log(f(exp(W)|DFN,DFD))
C
          IF (wswp.GT.zero) THEN
              d2ldw2 = -prd*exp(-wswp)/opkemw**2

          ELSE
              d2ldw2 = -prd*exp(wswp)/expwpk**2
          END IF

      ELSE
C
C     Calculate log(1 - F(exp(W)|DFN,DFD)) and its first and
C     its first and second derivatives if needed.
C
C        Calculate log(1 - F(exp(W)|DFN,DFD))
C
          CALL ltlf(w,two*a,two*b,lcum,l)
C
C        Calculate first and second derivative
C
C        Calculate values based on A, B and W
C
          lbda = log(b/a)
          logz = lbda - w
          logx = -dln1pe(-logz)
          logy = -dln1pe(logz)
          x = exp(logx)
          y = exp(logy)
C
C       Try Extreme Value Calculations
C
          IF (a.EQ.one) THEN
              dldw = -b*y
              d2ldw2 = -b*x*y
              GO TO 10

          END IF

          IF ((abs(a-one)*x.LT.0.1D0) .AND. (x.LT.maxx)) THEN
              CALL dextr(a,b,x,y,dldw,d2ldw2)
              dldw = -b*y + dldw
              d2ldw2 = -b*x*y + d2ldw2
              GO TO 10

          END IF
C
C        Try asymptotic Calculation
C
          thpab = thp*bswp
          IF (qswtch) THEN
              lambda = (a+b)*y - a

          ELSE
              lambda = b - (a+b)*x
          END IF

          IF ((a.GT.hundrd) .AND. (b.GT.hundrd) .AND.
     +        (lambda.GT.tiny) .AND. (lambda.LE.thpab)) THEN
              CALL dlasym(b,a,x,y,lambda,eps3,dldw,d2ldw2)
              GO TO 10

          END IF
C
C        Try Continued Fraction Calculation
C
          IF ((a.GT.forty) .AND. (b.GT.forty) .AND.
     +        (lambda.GE.thpab) .AND. (a*x.GE.pt7)) THEN
              CALL dlfrac(b,a,x,y,lambda,eps5,dldw,d2ldw2)
              dldw = -lambda + dldw
              d2ldw2 = - (a+b)*x*y + d2ldw2
              GO TO 10

          END IF
C
C        Try series Calculation
C
          z = dexp1(logz)
          IF (z.LE.1.0D13) THEN
              rat = (a-two)*z/ (b+two)

          ELSE
              rat = one
          END IF

          IF (abs(rat).LE.0.1D0) THEN
              CALL dser(a,b,z,dldw,d2ldw2)
              IF (qswtch) THEN
                  dldw = (a-one) - (a+b-one)*y + dldw

              ELSE
                  dldw = (a+b-one)*x - b + dldw
              END IF

              d2ldw2 = - (b+a-one)*y*x + d2ldw2
              GO TO 10

          END IF
C
C        Use direct differentiation method
C
          dldw = -tau*exp(xld-l)
          d2ldw2 = dldw* (xd1-dldw)
C
C     If CASE = 2 set L = LCUM and negative DLDW
C
   10     IF (case.EQ.2) THEN
              dldw = -dldw
          END IF

      END IF

      RETURN

      END

      DOUBLE PRECISION FUNCTION algdiv(a,b)
C-----------------------------------------------------------------------
C
C     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8
C
C                         --------
C
C     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
C     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION c,c0,c1,c2,c3,c4,c5,d,h,s11,s3,s5,s7,s9,t,u,v,w,
     +                 x,x2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alnrel
      EXTERNAL alnrel
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog
C     ..
C     .. Data statements ..
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (a.LE.b) GO TO 10
      h = b/a
      c = 1.0D0/ (1.0D0+h)
      x = h/ (1.0D0+h)
      d = a + (b-0.5D0)
      GO TO 20

   10 h = a/b
      c = h/ (1.0D0+h)
      x = 1.0D0/ (1.0D0+h)
      d = b + (a-0.5D0)
C
C                SET SN = (1 - X**N)/(1 - X)
C
   20 x2 = x*x
      s3 = 1.0D0 + (x+x2)
      s5 = 1.0D0 + (x+x2*s3)
      s7 = 1.0D0 + (x+x2*s5)
      s9 = 1.0D0 + (x+x2*s7)
      s11 = 1.0D0 + (x+x2*s9)
C
C                SET W = DEL(B) - DEL(A + B)
C
      t = (1.0D0/b)**2
      w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t + c0
      w = w* (c/b)
C
C                    COMBINE THE RESULTS
C
      u = d*alnrel(a/b)
      v = a* (dlog(b)-1.0D0)
      IF (u.LE.v) GO TO 30
      algdiv = (w-v) - u
      RETURN

   30 algdiv = (w-u) - v
      RETURN

      END
      DOUBLE PRECISION FUNCTION alnrel(a)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION LN(1 + A)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,p3,q1,q2,q3,t,t2,w,x
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog
C     ..
C     .. Data statements ..
      DATA p1/-.129418923021993D+01/,p2/.405303492862024D+00/,
     +     p3/-.178874546012214D-01/
      DATA q1/-.162752256355323D+01/,q2/.747811014037616D+00/,
     +     q3/-.845104217945565D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      IF (abs(a).GT.0.375D0) GO TO 10
      t = a/ (a+2.0D0)
      t2 = t*t
      w = (((p3*t2+p2)*t2+p1)*t2+1.0D0)/ (((q3*t2+q2)*t2+q1)*t2+1.0D0)
      alnrel = 2.0D0*t*w
      RETURN
C
   10 x = 1.D0 + dble(a)
      alnrel = dlog(x)
      RETURN

      END
      DOUBLE PRECISION FUNCTION apser(a,b,x,eps)
C-----------------------------------------------------------------------
C     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
C     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
C     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aj,bx,c,g,j,s,t,tol
C     ..
C     .. External Functions ..
      DOUBLE PRECISION psi
      EXTERNAL psi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog
C     ..
C     .. Data statements ..
C--------------------
      DATA g/.577215664901533D0/
C     ..
C     .. Executable Statements ..
C--------------------
      bx = b*x
      t = x - bx
      IF (b*eps.GT.2.D-2) GO TO 10
      c = dlog(x) + psi(b) + g + t
      GO TO 20

   10 c = dlog(bx) + g + t
C
   20 tol = 5.0D0*eps*abs(c)
      j = 1.0D0
      s = 0.0D0
   30 j = j + 1.0D0
      t = t* (x-bx/j)
      aj = t/j
      s = s + aj
      IF (abs(aj).GT.tol) GO TO 30
C
      apser = -a* (c+s)
      RETURN

      END
      DOUBLE PRECISION FUNCTION basym(a,b,lambda,eps)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bsum,dsum,e0,e1,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,
     +                 t0,t1,u,w,w0,z,z0,z2,zn,znm1
      INTEGER i,im1,imj,j,m,mm1,mmj,n,np1,num
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a0(21),b0(21),c(21),d(21)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION bcorr,erfc1,rlog1
      EXTERNAL bcorr,erfc1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sqrt
C     ..
C     .. Data statements ..
C------------------------
C     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
C            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
C------------------------
C     E0 = 2/SQRT(PI)
C     E1 = 2**(-3/2)
C------------------------
      DATA num/20/
      DATA e0/1.12837916709551D0/,e1/.353553390593274D0/
C     ..
C     .. Executable Statements ..
C------------------------
      basym = 0.0D0
      IF (a.GE.b) GO TO 10
      h = a/b
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/b
      w0 = 1.0D0/sqrt(a* (1.0D0+h))
      GO TO 20

   10 h = b/a
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/a
      w0 = 1.0D0/sqrt(b* (1.0D0+h))
C
   20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      t = exp(-f)
      IF (t.EQ.0.0D0) RETURN
      z0 = sqrt(f)
      z = 0.5D0* (z0/e1)
      z2 = f + f
C
      a0(1) = (2.0D0/3.0D0)*r1
      c(1) = -0.5D0*a0(1)
      d(1) = -c(1)
      j0 = (0.5D0/e0)*erfc1(1,z0)
      j1 = e1
      sum = j0 + d(1)*w0*j1
C
      s = 1.0D0
      h2 = h*h
      hn = 1.0D0
      w = w0
      znm1 = z
      zn = z2
      DO 70 n = 2,num,2
          hn = h2*hn
          a0(n) = 2.0D0*r0* (1.0D0+h*hn)/ (n+2.0D0)
          np1 = n + 1
          s = s + hn
          a0(np1) = 2.0D0*r1*s/ (n+3.0D0)
C
          DO 60 i = n,np1
              r = -0.5D0* (i+1.0D0)
              b0(1) = r*a0(1)
              DO 40 m = 2,i
                  bsum = 0.0D0
                  mm1 = m - 1
                  DO 30 j = 1,mm1
                      mmj = m - j
                      bsum = bsum + (j*r-mmj)*a0(j)*b0(mmj)
   30             CONTINUE
                  b0(m) = r*a0(m) + bsum/m
   40         CONTINUE
              c(i) = b0(i)/ (i+1.0D0)
C
              dsum = 0.0D0
              im1 = i - 1
              DO 50 j = 1,im1
                  imj = i - j
                  dsum = dsum + d(imj)*c(j)
   50         CONTINUE
              d(i) = - (dsum+c(i))
   60     CONTINUE
C
          j0 = e1*znm1 + (n-1.0D0)*j0
          j1 = e1*zn + n*j1
          znm1 = z2*znm1
          zn = z2*zn
          w = w0*w
          t0 = d(n)*w*j0
          w = w0*w
          t1 = d(np1)*w*j1
          sum = sum + (t0+t1)
          IF ((abs(t0)+abs(t1)).LE.eps*sum) GO TO 80
   70 CONTINUE
C
   80 u = exp(-bcorr(a,b))
      basym = e0*t*u*sum
      RETURN

      END
      DOUBLE PRECISION FUNCTION bcorr(a0,b0)
C-----------------------------------------------------------------------
C
C     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
C     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
C     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a0,b0
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,c,c0,c1,c2,c3,c4,c5,h,s11,s3,s5,s7,s9,t,w,x,
     +                 x2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dmax1,dmin1
C     ..
C     .. Data statements ..
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C------------------------
      a = dmin1(a0,b0)
      b = dmax1(a0,b0)
C
      h = a/b
      c = h/ (1.0D0+h)
      x = 1.0D0/ (1.0D0+h)
      x2 = x*x
C
C                SET SN = (1 - X**N)/(1 - X)
C
      s3 = 1.0D0 + (x+x2)
      s5 = 1.0D0 + (x+x2*s3)
      s7 = 1.0D0 + (x+x2*s5)
      s9 = 1.0D0 + (x+x2*s7)
      s11 = 1.0D0 + (x+x2*s9)
C
C                SET W = DEL(B) - DEL(A + B)
C
      t = (1.0D0/b)**2
      w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t + c0
      w = w* (c/b)
C
C                   COMPUTE  DEL(A) + W
C
      t = (1.0D0/a)**2
      bcorr = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a + w
      RETURN

      END
      DOUBLE PRECISION FUNCTION betaln(a0,b0)
C-----------------------------------------------------------------------
C     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
C-----------------------------------------------------------------------
C     E = 0.5*LN(2*PI)
C--------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a0,b0
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,c,e,h,u,v,w,z
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,bcorr,gamln,gsumln
      EXTERNAL algdiv,alnrel,bcorr,gamln,gsumln
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog,dmax1,dmin1
C     ..
C     .. Data statements ..
      DATA e/.918938533204673D0/
C     ..
C     .. Executable Statements ..
C--------------------------
      a = dmin1(a0,b0)
      b = dmax1(a0,b0)
      IF (a.GE.8.0D0) GO TO 100
      IF (a.GE.1.0D0) GO TO 20
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .LT. 1
C-----------------------------------------------------------------------
      IF (b.GE.8.0D0) GO TO 10
      betaln = gamln(a) + (gamln(b)-gamln(a+b))
      RETURN

   10 betaln = gamln(a) + algdiv(a,b)
      RETURN
C-----------------------------------------------------------------------
C                PROCEDURE WHEN 1 .LE. A .LT. 8
C-----------------------------------------------------------------------
   20 IF (a.GT.2.0D0) GO TO 40
      IF (b.GT.2.0D0) GO TO 30
      betaln = gamln(a) + gamln(b) - gsumln(a,b)
      RETURN

   30 w = 0.0D0
      IF (b.LT.8.0D0) GO TO 60
      betaln = gamln(a) + algdiv(a,b)
      RETURN
C
C                REDUCTION OF A WHEN B .LE. 1000
C
   40 IF (b.GT.1000.0D0) GO TO 80
      n = a - 1.0D0
      w = 1.0D0
      DO 50 i = 1,n
          a = a - 1.0D0
          h = a/b
          w = w* (h/ (1.0D0+h))
   50 CONTINUE
      w = dlog(w)
      IF (b.LT.8.0D0) GO TO 60
      betaln = w + gamln(a) + algdiv(a,b)
      RETURN
C
C                 REDUCTION OF B WHEN B .LT. 8
C
   60 n = b - 1.0D0
      z = 1.0D0
      DO 70 i = 1,n
          b = b - 1.0D0
          z = z* (b/ (a+b))
   70 CONTINUE
      betaln = w + dlog(z) + (gamln(a)+ (gamln(b)-gsumln(a,b)))
      RETURN
C
C                REDUCTION OF A WHEN B .GT. 1000
C
   80 n = a - 1.0D0
      w = 1.0D0
      DO 90 i = 1,n
          a = a - 1.0D0
          w = w* (a/ (1.0D0+a/b))
   90 CONTINUE
      betaln = (dlog(w)-n*dlog(b)) + (gamln(a)+algdiv(a,b))
      RETURN
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .GE. 8
C-----------------------------------------------------------------------
  100 w = bcorr(a,b)
      h = a/b
      c = h/ (1.0D0+h)
      u = - (a-0.5D0)*dlog(c)
      v = b*alnrel(h)
      IF (u.LE.v) GO TO 110
      betaln = (((-0.5D0*dlog(b)+e)+w)-v) - u
      RETURN

  110 betaln = (((-0.5D0*dlog(b)+e)+w)-u) - v
      RETURN

      END
      DOUBLE PRECISION FUNCTION bfrac(a,b,x,y,lambda,eps)
C-----------------------------------------------------------------------
C     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda,x,y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,
     +                 t,w,yp1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION brcomp
      EXTERNAL brcomp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..
C--------------------
      bfrac = brcomp(a,b,x,y)
      IF (bfrac.EQ.0.0D0) RETURN
C
      c = 1.0D0 + lambda
      c0 = b/a
      c1 = 1.0D0 + 1.0D0/a
      yp1 = y + 1.0D0
C
      n = 0.0D0
      p = 1.0D0
      s = a + 1.0D0
      an = 0.0D0
      bn = 1.0D0
      anp1 = 1.0D0
      bnp1 = c/c1
      r = c1/c
C
C        CONTINUED FRACTION CALCULATION
C
   10 n = n + 1.0D0
      t = n/a
      w = n* (b-n)*x
      e = a/s
      alpha = (p* (p+c0)*e*e)* (w*x)
      e = (1.0D0+t)/ (c1+t+t)
      beta = n + w/s + e* (c+n*yp1)
      p = 1.0D0 + t
      s = s + 2.0D0
C
C        UPDATE AN, BN, ANP1, AND BNP1
C
      t = alpha*an + beta*anp1
      an = anp1
      anp1 = t
      t = alpha*bn + beta*bnp1
      bn = bnp1
      bnp1 = t
C
      r0 = r
      r = anp1/bnp1
      IF (abs(r-r0).LE.eps*r) GO TO 20
C
C        RESCALE AN, BN, ANP1, AND BNP1
C
      an = an/bnp1
      bn = bn/bnp1
      anp1 = r
      bnp1 = 1.0D0
      GO TO 10
C
C                 TERMINATION
C
   20 bfrac = bfrac*r
      RETURN

      END
      SUBROUTINE bgrat(a,b,x,y,w,eps,ierr)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
C     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
C     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,w,x,y
      INTEGER ierr
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,
     +                 t2,u,v,z
      INTEGER i,n,nm1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION c(30),d(30)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,gam1
      EXTERNAL algdiv,alnrel,gam1
C     ..
C     .. External Subroutines ..
      EXTERNAL grat1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp
C     ..
C     .. Executable Statements ..
C
      bm1 = (b-0.5D0) - 0.5D0
      nu = a + 0.5D0*bm1
      IF (y.GT.0.375D0) GO TO 10
      lnx = alnrel(-y)
      GO TO 20

   10 lnx = dlog(x)
   20 z = -nu*lnx
      IF (b*z.EQ.0.0D0) GO TO 70
C
C                 COMPUTATION OF THE EXPANSION
C                 SET R = EXP(-Z)*Z**B/GAMMA(B)
C
      r = b* (1.0D0+gam1(b))*exp(b*dlog(z))
      r = r*exp(a*lnx)*exp(0.5D0*bm1*lnx)
      u = algdiv(b,a) + b*dlog(nu)
      u = r*exp(-u)
      IF (u.EQ.0.0D0) GO TO 70
      CALL grat1(b,z,r,p,q,eps)
C
      v = 0.25D0* (1.0D0/nu)**2
      t2 = 0.25D0*lnx*lnx
      l = w/u
      j = q/r
      sum = j
      t = 1.0D0
      cn = 1.0D0
      n2 = 0.0D0
      DO 50 n = 1,30
          bp2n = b + n2
          j = (bp2n* (bp2n+1.0D0)*j+ (z+bp2n+1.0D0)*t)*v
          n2 = n2 + 2.0D0
          t = t*t2
          cn = cn/ (n2* (n2+1.0D0))
          c(n) = cn
          s = 0.0D0
          IF (n.EQ.1) GO TO 40
          nm1 = n - 1
          coef = b - n
          DO 30 i = 1,nm1
              s = s + coef*c(i)*d(n-i)
              coef = coef + b
   30     CONTINUE
   40     d(n) = bm1*cn + s/n
          dj = d(n)*j
          sum = sum + dj
          IF (sum.LE.0.0D0) GO TO 70
          IF (abs(dj).LE.eps* (sum+l)) GO TO 60
   50 CONTINUE
C
C                    ADD THE RESULTS TO W
C
   60 ierr = 0
      w = w + u*sum
      RETURN
C
C               THE EXPANSION CANNOT BE COMPUTED
C
   70 ierr = 1
      RETURN

      END
      DOUBLE PRECISION FUNCTION bpser(a,b,x,eps)
C-----------------------------------------------------------------------
C     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
C     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,apb,b0,c,n,sum,t,tol,u,w,z
      INTEGER i,m
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,betaln,gam1,gamln1
      EXTERNAL algdiv,betaln,gam1,gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp
C     ..
C     .. Executable Statements ..
C
      bpser = 0.0D0
      IF (x.EQ.0.0D0) RETURN
C-----------------------------------------------------------------------
C            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
C-----------------------------------------------------------------------
      a0 = dmin1(a,b)
      IF (a0.LT.1.0D0) GO TO 10
      z = a*dlog(x) - betaln(a,b)
      bpser = exp(z)/a
      GO TO 100

   10 b0 = dmax1(a,b)
      IF (b0.GE.8.0D0) GO TO 90
      IF (b0.GT.1.0D0) GO TO 40
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
C
      bpser = x**a
      IF (bpser.EQ.0.0D0) RETURN
C
      apb = a + b
      IF (apb.GT.1.0D0) GO TO 20
      z = 1.0D0 + gam1(apb)
      GO TO 30

   20 u = dble(a) + dble(b) - 1.D0
      z = (1.0D0+gam1(u))/apb
C
   30 c = (1.0D0+gam1(a))* (1.0D0+gam1(b))/z
      bpser = bpser*c* (b/apb)
      GO TO 100
C
C         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
C
   40 u = gamln1(a0)
      m = b0 - 1.0D0
      IF (m.LT.1) GO TO 60
      c = 1.0D0
      DO 50 i = 1,m
          b0 = b0 - 1.0D0
          c = c* (b0/ (a0+b0))
   50 CONTINUE
      u = dlog(c) + u
C
   60 z = a*dlog(x) - u
      b0 = b0 - 1.0D0
      apb = a0 + b0
      IF (apb.GT.1.0D0) GO TO 70
      t = 1.0D0 + gam1(apb)
      GO TO 80

   70 u = dble(a0) + dble(b0) - 1.D0
      t = (1.0D0+gam1(u))/apb
   80 bpser = exp(z)* (a0/a)* (1.0D0+gam1(b0))/t
      GO TO 100
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
C
   90 u = gamln1(a0) + algdiv(a0,b0)
      z = a*dlog(x) - u
      bpser = (a0/a)*exp(z)
  100 IF (bpser.EQ.0.0D0 .OR. a.LE.0.1D0*eps) RETURN
C-----------------------------------------------------------------------
C                     COMPUTE THE SERIES
C-----------------------------------------------------------------------
      sum = 0.0D0
      n = 0.0D0
      c = 1.0D0
      tol = eps/a
  110 n = n + 1.0D0
      c = c* (0.5D0+ (0.5D0-b/n))*x
      w = c/ (a+n)
      sum = sum + w
      IF (abs(w).GT.tol) GO TO 110
      bpser = bpser* (1.0D0+a*sum)
      RETURN

      END
      SUBROUTINE bratio(a,b,x,y,w,w1,ierr)
C-----------------------------------------------------------------------
C
C            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)
C
C                     --------------------
C
C     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
C     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
C
C                      W  = IX(A,B)
C                      W1 = 1 - IX(A,B)
C
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
C     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
C     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
C     ONE OF THE FOLLOWING VALUES ...
C
C        IERR = 1  IF A OR B IS NEGATIVE
C        IERR = 2  IF A = B = 0
C        IERR = 3  IF X .LT. 0 OR X .GT. 1
C        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
C        IERR = 5  IF X + Y .NE. 1
C        IERR = 6  IF X = A = 0
C        IERR = 7  IF Y = B = 0
C
C--------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN, VIRGINIA
C     REVISED ... NOV 1991
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,w,w1,x,y
      INTEGER ierr
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,b0,eps,lambda,t,x0,y0,z
      INTEGER ierr1,ind,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION apser,basym,bfrac,bpser,bup,fpser,spmpar
      EXTERNAL apser,basym,bfrac,bpser,bup,fpser,spmpar
C     ..
C     .. External Subroutines ..
      EXTERNAL bgrat
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dmax1,dmin1
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
C            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
C
      eps = spmpar(1)
C
C-----------------------------------------------------------------------
      w = 0.0D0
      w1 = 0.0D0
      IF (a.LT.0.0D0 .OR. b.LT.0.0D0) GO TO 270
      IF (a.EQ.0.0D0 .AND. b.EQ.0.0D0) GO TO 280
      IF (x.LT.0.0D0 .OR. x.GT.1.0D0) GO TO 290
      IF (y.LT.0.0D0 .OR. y.GT.1.0D0) GO TO 300
      z = ((x+y)-0.5D0) - 0.5D0
      IF (abs(z).GT.3.0D0*eps) GO TO 310
C
      ierr = 0
      IF (x.EQ.0.0D0) GO TO 210
      IF (y.EQ.0.0D0) GO TO 230
      IF (a.EQ.0.0D0) GO TO 240
      IF (b.EQ.0.0D0) GO TO 220
C
      eps = dmax1(eps,1.D-15)
      IF (dmax1(a,b).LT.1.D-3*eps) GO TO 260
C
      ind = 0
      a0 = a
      b0 = b
      x0 = x
      y0 = y
      IF (dmin1(a0,b0).GT.1.0D0) GO TO 40
C
C             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
C
      IF (x.LE.0.5D0) GO TO 10
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
C
   10 IF (b0.LT.dmin1(eps,eps*a0)) GO TO 90
      IF (a0.LT.dmin1(eps,eps*b0) .AND. b0*x0.LE.1.0D0) GO TO 100
      IF (dmax1(a0,b0).GT.1.0D0) GO TO 20
      IF (a0.GE.dmin1(0.2D0,b0)) GO TO 110
      IF (x0**a0.LE.0.9D0) GO TO 110
      IF (x0.GE.0.3D0) GO TO 120
      n = 20
      GO TO 140
C
   20 IF (b0.LE.1.0D0) GO TO 110
      IF (x0.GE.0.3D0) GO TO 120
      IF (x0.GE.0.1D0) GO TO 30
      IF ((x0*b0)**a0.LE.0.7D0) GO TO 110
   30 IF (b0.GT.15.0D0) GO TO 150
      n = 20
      GO TO 140
C
C             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
C
   40 IF (a.GT.b) GO TO 50
      lambda = a - (a+b)*x
      GO TO 60

   50 lambda = (a+b)*y - b
   60 IF (lambda.GE.0.0D0) GO TO 70
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
      lambda = abs(lambda)
C
   70 IF (b0.LT.40.0D0 .AND. b0*x0.LE.0.7D0) GO TO 110
      IF (b0.LT.40.0D0) GO TO 160
      IF (a0.GT.b0) GO TO 80
      IF (a0.LE.100.0D0) GO TO 130
      IF (lambda.GT.0.03D0*a0) GO TO 130
      GO TO 200

   80 IF (b0.LE.100.0D0) GO TO 130
      IF (lambda.GT.0.03D0*b0) GO TO 130
      GO TO 200
C
C            EVALUATION OF THE APPROPRIATE ALGORITHM
C
   90 w = fpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  100 w1 = apser(a0,b0,x0,eps)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  110 w = bpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  120 w1 = bpser(b0,a0,y0,eps)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  130 w = bfrac(a0,b0,x0,y0,lambda,15.0D0*eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  140 w1 = bup(b0,a0,y0,x0,n,eps)
      b0 = b0 + n
  150 CALL bgrat(b0,a0,y0,x0,w1,15.0D0*eps,ierr1)
      w = 0.5D0 + (0.5D0-w1)
      GO TO 250
C
  160 n = b0
      b0 = b0 - n
      IF (b0.NE.0.0D0) GO TO 170
      n = n - 1
      b0 = 1.0D0
  170 w = bup(b0,a0,y0,x0,n,eps)
      IF (x0.GT.0.7D0) GO TO 180
      w = w + bpser(a0,b0,x0,eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  180 IF (a0.GT.15.0D0) GO TO 190
      n = 20
      w = w + bup(a0,b0,x0,y0,n,eps)
      a0 = a0 + n
  190 CALL bgrat(a0,b0,x0,y0,w,15.0D0*eps,ierr1)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
  200 w = basym(a0,b0,lambda,100.0D0*eps)
      w1 = 0.5D0 + (0.5D0-w)
      GO TO 250
C
C               TERMINATION OF THE PROCEDURE
C
  210 IF (a.EQ.0.0D0) GO TO 320
  220 w = 0.0D0
      w1 = 1.0D0
      RETURN
C
  230 IF (b.EQ.0.0D0) GO TO 330
  240 w = 1.0D0
      w1 = 0.0D0
      RETURN
C
  250 IF (ind.EQ.0) RETURN
      t = w
      w = w1
      w1 = t
      RETURN
C
C           PROCEDURE FOR A AND B .LT. 1.E-3*EPS
C
  260 w = b/ (a+b)
      w1 = a/ (a+b)
      RETURN
C
C                       ERROR RETURN
C
  270 ierr = 1
      RETURN

  280 ierr = 2
      RETURN

  290 ierr = 3
      RETURN

  300 ierr = 4
      RETURN

  310 ierr = 5
      RETURN

  320 ierr = 6
      RETURN

  330 ierr = 7
      RETURN

      END
      DOUBLE PRECISION FUNCTION brcmp1(mu,a,b,x,y)
C-----------------------------------------------------------------------
C          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,x,y
      INTEGER mu
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,apb,b0,c,const,e,h,lambda,lnx,lny,t,u,v,x0,y0,
     +                 z
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,bcorr,betaln,esum,gam1,gamln1,rlog1
      EXTERNAL algdiv,alnrel,bcorr,betaln,esum,gam1,gamln1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp,sqrt
C     ..
C     .. Data statements ..
C-----------------
C     CONST = 1/SQRT(2*PI)
C-----------------
      DATA const/.398942280401433D0/
C     ..
C     .. Executable Statements ..
C
      a0 = dmin1(a,b)
      IF (a0.GE.8.0D0) GO TO 130
C
      IF (x.GT.0.375D0) GO TO 10
      lnx = dlog(x)
      lny = alnrel(-x)
      GO TO 30

   10 IF (y.GT.0.375D0) GO TO 20
      lnx = alnrel(-y)
      lny = dlog(y)
      GO TO 30

   20 lnx = dlog(x)
      lny = dlog(y)
C
   30 z = a*lnx + b*lny
      IF (a0.LT.1.0D0) GO TO 40
      z = z - betaln(a,b)
      brcmp1 = esum(mu,z)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .LT. 1 OR B .LT. 1
C-----------------------------------------------------------------------
   40 b0 = dmax1(a,b)
      IF (b0.GE.8.0D0) GO TO 120
      IF (b0.GT.1.0D0) GO TO 70
C
C                   ALGORITHM FOR B0 .LE. 1
C
      brcmp1 = esum(mu,z)
      IF (brcmp1.EQ.0.0D0) RETURN
C
      apb = a + b
      IF (apb.GT.1.0D0) GO TO 50
      z = 1.0D0 + gam1(apb)
      GO TO 60

   50 u = dble(a) + dble(b) - 1.D0
      z = (1.0D0+gam1(u))/apb
C
   60 c = (1.0D0+gam1(a))* (1.0D0+gam1(b))/z
      brcmp1 = brcmp1* (a0*c)/ (1.0D0+a0/b0)
      RETURN
C
C                ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   70 u = gamln1(a0)
      n = b0 - 1.0D0
      IF (n.LT.1) GO TO 90
      c = 1.0D0
      DO 80 i = 1,n
          b0 = b0 - 1.0D0
          c = c* (b0/ (a0+b0))
   80 CONTINUE
      u = dlog(c) + u
C
   90 z = z - u
      b0 = b0 - 1.0D0
      apb = a0 + b0
      IF (apb.GT.1.0D0) GO TO 100
      t = 1.0D0 + gam1(apb)
      GO TO 110

  100 u = dble(a0) + dble(b0) - 1.D0
      t = (1.0D0+gam1(u))/apb
  110 brcmp1 = a0*esum(mu,z)* (1.0D0+gam1(b0))/t
      RETURN
C
C                   ALGORITHM FOR B0 .GE. 8
C
  120 u = gamln1(a0) + algdiv(a0,b0)
      brcmp1 = a0*esum(mu,z-u)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .GE. 8 AND B .GE. 8
C-----------------------------------------------------------------------
  130 IF (a.GT.b) GO TO 140
      h = a/b
      x0 = h/ (1.0D0+h)
      y0 = 1.0D0/ (1.0D0+h)
      lambda = a - (a+b)*x
      GO TO 150

  140 h = b/a
      x0 = 1.0D0/ (1.0D0+h)
      y0 = h/ (1.0D0+h)
      lambda = (a+b)*y - b
C
  150 e = -lambda/a
      IF (abs(e).GT.0.6D0) GO TO 160
      u = rlog1(e)
      GO TO 170

  160 u = e - dlog(x/x0)
C
  170 e = lambda/b
      IF (abs(e).GT.0.6D0) GO TO 180
      v = rlog1(e)
      GO TO 190

  180 v = e - dlog(y/y0)
C
  190 z = esum(mu,- (a*u+b*v))
      brcmp1 = const*sqrt(b*x0)*z*exp(-bcorr(a,b))
      RETURN

      END
      DOUBLE PRECISION FUNCTION brcomp(a,b,x,y)
C-----------------------------------------------------------------------
C               EVALUATION OF X**A*Y**B/BETA(A,B)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,x,y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a0,apb,b0,c,const,e,h,lambda,lnx,lny,t,u,v,x0,y0,
     +                 z
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION algdiv,alnrel,bcorr,betaln,gam1,gamln1,rlog1
      EXTERNAL algdiv,alnrel,bcorr,betaln,gam1,gamln1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog,dmax1,dmin1,exp,sqrt
C     ..
C     .. Data statements ..
C-----------------
C     CONST = 1/SQRT(2*PI)
C-----------------
      DATA const/.398942280401433D0/
C     ..
C     .. Executable Statements ..
C
      brcomp = 0.0D0
      IF (x.EQ.0.0D0 .OR. y.EQ.0.0D0) RETURN
      a0 = dmin1(a,b)
      IF (a0.GE.8.0D0) GO TO 130
C
      IF (x.GT.0.375D0) GO TO 10
      lnx = dlog(x)
      lny = alnrel(-x)
      GO TO 30

   10 IF (y.GT.0.375D0) GO TO 20
      lnx = alnrel(-y)
      lny = dlog(y)
      GO TO 30

   20 lnx = dlog(x)
      lny = dlog(y)
C
   30 z = a*lnx + b*lny
      IF (a0.LT.1.0D0) GO TO 40
      z = z - betaln(a,b)
      brcomp = exp(z)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .LT. 1 OR B .LT. 1
C-----------------------------------------------------------------------
   40 b0 = dmax1(a,b)
      IF (b0.GE.8.0D0) GO TO 120
      IF (b0.GT.1.0D0) GO TO 70
C
C                   ALGORITHM FOR B0 .LE. 1
C
      brcomp = exp(z)
      IF (brcomp.EQ.0.0D0) RETURN
C
      apb = a + b
      IF (apb.GT.1.0D0) GO TO 50
      z = 1.0D0 + gam1(apb)
      GO TO 60

   50 u = dble(a) + dble(b) - 1.D0
      z = (1.0D0+gam1(u))/apb
C
   60 c = (1.0D0+gam1(a))* (1.0D0+gam1(b))/z
      brcomp = brcomp* (a0*c)/ (1.0D0+a0/b0)
      RETURN
C
C                ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   70 u = gamln1(a0)
      n = b0 - 1.0D0
      IF (n.LT.1) GO TO 90
      c = 1.0D0
      DO 80 i = 1,n
          b0 = b0 - 1.0D0
          c = c* (b0/ (a0+b0))
   80 CONTINUE
      u = dlog(c) + u
C
   90 z = z - u
      b0 = b0 - 1.0D0
      apb = a0 + b0
      IF (apb.GT.1.0D0) GO TO 100
      t = 1.0D0 + gam1(apb)
      GO TO 110

  100 u = dble(a0) + dble(b0) - 1.D0
      t = (1.0D0+gam1(u))/apb
  110 brcomp = a0*exp(z)* (1.0D0+gam1(b0))/t
      RETURN
C
C                   ALGORITHM FOR B0 .GE. 8
C
  120 u = gamln1(a0) + algdiv(a0,b0)
      brcomp = a0*exp(z-u)
      RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .GE. 8 AND B .GE. 8
C-----------------------------------------------------------------------
  130 IF (a.GT.b) GO TO 140
      h = a/b
      x0 = h/ (1.0D0+h)
      y0 = 1.0D0/ (1.0D0+h)
      lambda = a - (a+b)*x
      GO TO 150

  140 h = b/a
      x0 = 1.0D0/ (1.0D0+h)
      y0 = h/ (1.0D0+h)
      lambda = (a+b)*y - b
C
  150 e = -lambda/a
      IF (abs(e).GT.0.6D0) GO TO 160
      u = rlog1(e)
      GO TO 170

  160 u = e - dlog(x/x0)
C
  170 e = lambda/b
      IF (abs(e).GT.0.6D0) GO TO 180
      v = rlog1(e)
      GO TO 190

  180 v = e - dlog(y/y0)
C
  190 z = exp(- (a*u+b*v))
      brcomp = const*sqrt(b*x0)*z*exp(-bcorr(a,b))
      RETURN

      END
      DOUBLE PRECISION FUNCTION bup(a,b,x,y,n,eps)
C-----------------------------------------------------------------------
C     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
C     EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x,y
      INTEGER n
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ap1,apb,d,l,r,t,w
      INTEGER i,k,kp1,mu,nm1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION brcmp1,exparg
      EXTERNAL brcmp1,exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Executable Statements ..
C
C          OBTAIN THE SCALING FACTOR EXP(-MU) AND
C             EXP(MU)*(X**A*Y**B/BETA(A,B))/A
C
      apb = a + b
      ap1 = a + 1.0D0
      mu = 0
      d = 1.0D0
      IF (n.EQ.1 .OR. a.LT.1.0D0) GO TO 10
      IF (apb.LT.1.1D0*ap1) GO TO 10
      mu = abs(exparg(1))
      k = exparg(0)
      IF (k.LT.mu) mu = k
      t = mu
      d = exp(-t)
C
   10 bup = brcmp1(mu,a,b,x,y)/a
      IF (n.EQ.1 .OR. bup.EQ.0.0D0) RETURN
      nm1 = n - 1
      w = d
C
C          LET K BE THE INDEX OF THE MAXIMUM TERM
C
      k = 0
      IF (b.LE.1.0D0) GO TO 50
      IF (y.GT.1.D-4) GO TO 20
      k = nm1
      GO TO 30

   20 r = (b-1.0D0)*x/y - a
      IF (r.LT.1.0D0) GO TO 50
      k = nm1
      t = nm1
      IF (r.LT.t) k = r
C
C          ADD THE INCREASING TERMS OF THE SERIES
C
   30 DO 40 i = 1,k
          l = i - 1
          d = ((apb+l)/ (ap1+l))*x*d
          w = w + d
   40 CONTINUE
      IF (k.EQ.nm1) GO TO 70
C
C          ADD THE REMAINING TERMS OF THE SERIES
C
   50 kp1 = k + 1
      DO 60 i = kp1,nm1
          l = i - 1
          d = ((apb+l)/ (ap1+l))*x*d
          w = w + d
          IF (d.LE.eps*w) GO TO 70
   60 CONTINUE
C
C               TERMINATE THE PROCEDURE
C
   70 bup = bup*w
      RETURN

      END
      DOUBLE PRECISION FUNCTION cenlf(w,k)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION CENLF(W)
C          CENtral part of the Log F density
C
C     Returns Ninth Degree Taylor's series expansion of
C          W - (1+K) * LOG( 1 + (EXP(W)-1) / (1+K) )
C
C     Is called only for K < 1 and W <= 0.1
C
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (a-h,o-p,r-z),INTEGER (i-n),LOGICAL (q)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION k,w
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION denom,numer
      INTEGER i,icoef,idnum
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION coef(9),mult(3:9),numcof(7,3:9)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..
C     .. Save statement ..
      DOUBLE PRECISION kold
      SAVE coef,kold
C     ..
C     .. Data statements ..
      DATA kold/-1.0D0/
      DATA mult/0.5D0,0.16666666666666666667D0,
     +     0.41666666666666666667D-1,0.83333333333333333333D-2,
     +     0.13888888888888888889D-2,0.19841269841269841270D-3,
     +     0.24801587301587301587D-4/
      DATA (numcof(i,3),i=1,7)/1.0D0,6*0.0D0/
      DATA (numcof(i,4),i=1,7)/-1.0D0,1.0D0,5*0.0D0/
      DATA (numcof(i,5),i=1,7)/1.0D0,-4.0D0,1.0D0,4*0.0D0/
      DATA (numcof(i,6),i=1,7)/-1.0D0,11.0D0,-11.0D0,1.0D0,3*0.0D0/
      DATA (numcof(i,7),i=1,7)/1.0D0,-26.0D0,66.0D0,-26.0D0,1.0D0,
     +     2*0.0D0/
      DATA (numcof(i,8),i=1,7)/-1.0D0,57.0D0,-302.0D0,302.0D0,-57.0D0,
     +     1.0D0,0.0D0/
      DATA (numcof(i,9),i=1,7)/1.0D0,-120.0D0,1191.0D0,-2416.0D0,
     +     1191.0D0,-120.0D0,1.0D0/
C     ..
C     .. Executable Statements ..
C
C
C     If K changes, set up coefficients for new K
C
      IF (kold.NE.k) THEN
          kold = k
          DO 10,i = 1,2
              coef(i) = 0.0D0
   10     CONTINUE
          denom = 1.0D0
          DO 20,icoef = 3,9
              idnum = icoef - 2
              numer = -mult(icoef)*k*devlpl(numcof(1,icoef),idnum,k)
              denom = denom* (1.0D0+k)
              coef(icoef) = numer/denom
   20     CONTINUE
      END IF
C
C     Calculate CENLF
C
      cenlf = devlpl(coef,9,w)

      RETURN

      END
      DOUBLE PRECISION FUNCTION dbetrm(a,b)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DBETRM( A, B )
C          Double Precision Sterling Remainder for Complete
C                    Beta Function
C
C
C                              Function
C
C
C     Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
C     where Lgamma is the log of the (complete) gamma function
C
C     Let ZZ be approximation obtained if each log gamma is approximated
C     by Sterling's formula, i.e.,
C     Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
C
C     Returns Log(Beta(A,B)) - ZZ
C
C
C                              Arguments
C
C
C     A --> One argument of the Beta
C                    DOUBLE PRECISION A
C
C     B --> The other argument of the Beta
C                    DOUBLE PRECISION B
C
C**********************************************************************

C     .. Parameters ..
      DOUBLE PRECISION hln2pi
      PARAMETER (hln2pi=0.91893853320467274178D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dstrem
      EXTERNAL dstrem
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC max,min
C     ..
C     .. Executable Statements ..

C     Try to sum from smallest to largest
      dbetrm = -dstrem(a+b)
      dbetrm = dbetrm + dstrem(max(a,b))
      dbetrm = dbetrm + dstrem(min(a,b))
      dbetrm = dbetrm + hln2pi
      RETURN

      END
      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C              Double precision EVALuate a PoLynomial at X
C
C
C                              Function
C
C
C     returns
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)
C
C
C                              Arguments
C
C
C     A --> Array of coefficients of the polynomial.
C                                        A is DOUBLE PRECISION(N)
C
C     N --> Length of A, also degree of polynomial - 1.
C                                        N is INTEGER
C
C     X --> Point at which the polynomial is to be evaluated.
C                                        X is DOUBLE PRECISION
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION a(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
C     ..
C     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term
      RETURN

      END
      DOUBLE PRECISION FUNCTION dexp1(x)
C**********************************************************************
C
C      DOUBLE PRECISION FUNCTION dexp1(x)
C            Evaluation of the function EXP(X) with no overflow
C
C
C                              Arguments
C
C
C     X --> Argument at which exp(x) desired
C                    DOUBLE PRECISION X
C
C
C                              Method
C
C
C     If Exp(x) > largest double precision value for machine then dexp1
C     returns the largest double precision value for the machine;
C     otherwise it returns exp(x)
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION max,logmax
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,dlog
C     ..
C     .. Save statements ..
      SAVE max,logmax
C     ..
C     .. Data statements ..
      DATA max/0.0d0/,logmax/0.0d0/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (max.EQ.0.0d0)) GO TO 10
      max = spmpar(3)
      logmax = dlog(max)
   10 IF (.NOT. (x.GT.logmax)) GO TO 20
      dexp1 = max
      GO TO 30
C
   20 dexp1 = exp(x)
C
   30 RETURN

      END
      DOUBLE PRECISION FUNCTION dexpm1(x)
C**********************************************************************
C
C      DOUBLE PRECISION FUNCTION dexpm1(x)
C            Evaluation of the function EXP(X) - 1
C
C
C                              Arguments
C
C
C     X --> Argument at which exp(x)-1 desired
C                    DOUBLE PRECISION X
C
C
C                              Method
C
C
C     Renaming of function rexp from code of:
C
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,q1,q2,q3,q4,w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA p1/.914041914819518D-09/,p2/.238082361044469D-01/,
     +     q1/-.499999999085958D+00/,q2/.107141568980644D+00/,
     +     q3/-.119041179760821D-01/,q4/.595130811860248D-03/
C     ..
C     .. Executable Statements ..
C
      IF (abs(x).GT.0.15D0) GO TO 10
      dexpm1 = x* (((p2*x+p1)*x+1.0D0)/
     +         ((((q4*x+q3)*x+q2)*x+q1)*x+1.0D0))
      RETURN
C
   10 w = exp(x)
      IF (x.GT.0.0D0) GO TO 20
      dexpm1 = (w-0.5D0) - 0.5D0
      RETURN

   20 dexpm1 = w* (0.5D0+ (0.5D0-1.0D0/w))
      RETURN

      END
      SUBROUTINE dextr(a,b,x,y,d1,d2)
C--------------------------------------------------------------------
C
C    Calculates first and second derivatives of the log of
C          1 + A*[(1-A)(2-A)...(j-B)*X^j]/[j!(A+j)] j=1,2,...
C    This is the summation part of the extreme value formula
C          G(A,B)*[(X^B)/B]*SUM
C    which can be used to calculate IX(A,B).
C
C--------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
      DOUBLE PRECISION two
      PARAMETER (two=2.0D0)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=1.0D-13)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,x,y,d1,d2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION xd1,xd2,xnma,xnpb,term,termd1,termd2,termo,
     +                 trmd1o,trmd2o,sum,sumt,sumd1,sumd1t,sumd2,sumd2t,
     +                 temp,xn,rat,rat1,rat2
      LOGICAL qdone
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..
C
C     initialize variables
C
C
      qdone = .FALSE.
      xnma = one - a
      xnpb = one + b
      xd1 = -x*y
      xd2 = xd1* (two*x-one)
      term = xnma*x
      termd1 = xnma*xd1
      termd2 = xnma*xd2
      sum = zero
      sumd1 = zero
      sumd2 = zero
      xn = one
C
C     Compute sum
C
      GO TO 20

   10 IF (qdone .OR. (xn.GT.100D0)) GO TO 90
   20 sumt = term/xnpb
      sum = sum + sumt
      sumd1t = termd1/xnpb
      sumd1 = sumd1 + sumd1t
      sumd2t = termd2/xnpb
      sumd2 = sumd2 + sumd2t
      IF (.NOT. (sum.NE.zero)) GO TO 30
      rat = abs(sumt/sum)
      GO TO 40

   30 rat = one
   40 IF (.NOT. (sumd1.NE.zero)) GO TO 50
      rat1 = abs(sumd1t/sumd1)
      GO TO 60

   50 rat1 = one
   60 IF (.NOT. (sumd2.NE.zero)) GO TO 70
      rat2 = abs(sumd2t/sumd2)
      GO TO 80

   70 rat2 = one
   80 qdone = (((rat.LE.tiny).OR. (sumt.EQ.zero)) .AND.
     +        ((rat1.LE.tiny).OR. (sumd1t.EQ.zero)) .AND.
     +        ((rat2.LE.tiny).OR. (sumd2t.EQ.zero)))
      IF (.NOT.qdone) THEN
          xn = xn + one
          xnma = one + xnma
          xnpb = one + xnpb
          termo = term
          trmd1o = termd1
          trmd2o = termd2
          temp = xnma/xn
          term = temp*termo*x
          termd1 = temp* (termo*xd1+trmd1o*x)
          termd2 = temp* (termo*xd2+two*trmd1o*xd1+trmd2o*x)
      END IF

      GO TO 10

   90 temp = b/ (one+b*sum)

      d1 = sumd1*temp
      d2 = sumd2*temp - (d1*d1)

      RETURN

      END
      SUBROUTINE dlasym(a,b,x,y,lambda,eps,d1,d2)
C-----------------------------------------------------------------------
C
C     CALCULATES D1 AND D2 THE FIRST AND SECOND DERIATIVES
C     OF LOG(IX(A,B)) USING ASYMPTOTIC EXPANSION FOR LARGE A AND B.
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C
C     Derivatives were taken from LASYM which is a modification of
C     BASYM from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,x,y,eps,lambda,d1,d2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bsum,dsum,e0,e1,e2,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,
     +                 t0,t1,w,w0,z,z0,z2,zn,znm1,lambd1,lambd2,fd1,fd2,
     +                 j0d1,j0d2,znd1,znd2,znm1d1,znm1d2,j1d1,j1d2,t0d1,
     +                 t0d2,t1d1,t1d2,sumd1,sumd2,lambb,lamba,z12,z32,
     +                 zno,znd1o,znd2o,zn1o,zn1d1o,zn1d2o,temp
      INTEGER i,im1,imj,j,m,mm1,mmj,n,np1,num
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a0(21),b0(21),c(21),d(21)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION bcorr,erfc1,rlog1
      EXTERNAL bcorr,erfc1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt,log
C     ..
C     .. Data statements ..
C------------------------
C     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
C            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
C------------------------
C     E0 = 2/SQRT(PI)
C     E1 = 2**(-3/2)
C     E2 = 1/SQRT(PI)
C------------------------
      DATA num/20/
      DATA e0/1.12837916709551D0/,e1/.353553390593274D0/,
     +     e2/0.5641895835477563D0/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (a.GE.b) GO TO 10
      h = a/b
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/b
      w0 = 1.0D0/sqrt(a* (1.0D0+h))
      GO TO 20

   10 h = b/a
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/a
      w0 = 1.0D0/sqrt(b* (1.0D0+h))
C
   20 lambd1 = (a+b)*x*y
      lambd2 = lambd1* (2.0d0*x-1.0d0)
      lambb = 1.0d0/ (1.0d0+ (lambda/b))
      lamba = 1.0d0/ (0.5d0+ (0.5d0-lambda/a))
      temp = lamba - lambb
      f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      fd1 = lambd1*temp
      fd2 = lambd2*temp + (lambd1**2.0d0)*
     +      (((lamba**2.0d0)/a)+ ((lambb**2.0d0)/b))
      z0 = sqrt(f)
      z = 0.5D0* (z0/e1)
      z2 = f + f
      z12 = 1.0d0/z0
      z32 = z12/f
C
      a0(1) = (2.0D0/3.0D0)*r1
      c(1) = -0.5D0*a0(1)
      d(1) = -c(1)
      j0 = erfc1(1,z0)
      temp = (j0- (e2*z12))
      j0d1 = fd1*temp
      j0d2 = fd1* (j0d1+ (e2*0.5d0*fd1*z32)) + fd2*temp
      j0 = (0.5d0/e0)*j0
      j0d1 = (0.5d0/e0)*j0d1
      j0d2 = (0.5d0/e0)*j0d2
      j1 = e1
      j1d1 = 0.0d0
      j1d2 = 0.0d0
      sum = j0 + d(1)*w0*j1
      sumd1 = j0d1
      sumd2 = j0d2
C
      s = 1.0D0
      h2 = h*h
      hn = 1.0D0
      w = w0
      znm1 = z
      znm1d1 = fd1*0.25d0*z12/e1
      znm1d2 = (0.25d0/e1)* (z12*fd2-0.5d0*z32*fd1*fd1)
      zn = z2
      znd1 = 2.0d0*fd1
      znd2 = 2.0d0*fd2
      DO 70 n = 2,num,2
          hn = h2*hn
          a0(n) = 2.0D0*r0* (1.0D0+h*hn)/ (n+2.0D0)
          np1 = n + 1
          s = s + hn
          a0(np1) = 2.0D0*r1*s/ (n+3.0D0)
C
          DO 60 i = n,np1
              r = -0.5D0* (i+1.0D0)
              b0(1) = r*a0(1)
              DO 40 m = 2,i
                  bsum = 0.0D0
                  mm1 = m - 1
                  DO 30 j = 1,mm1
                      mmj = m - j
                      bsum = bsum + (j*r-mmj)*a0(j)*b0(mmj)
   30             CONTINUE
                  b0(m) = r*a0(m) + bsum/m
   40         CONTINUE
              c(i) = b0(i)/ (i+1.0D0)
C
              dsum = 0.0D0
              im1 = i - 1
              DO 50 j = 1,im1
                  imj = i - j
                  dsum = dsum + d(imj)*c(j)
   50         CONTINUE
              d(i) = - (dsum+c(i))
   60     CONTINUE
C
          j0 = e1*znm1 + (n-1.0D0)*j0
          j0d1 = e1*znm1d1 + (n-1.0d0)*j0d1
          j0d2 = e1*znm1d2 + (n-1.0d0)*j0d2
          j1 = e1*zn + n*j1
          j1d1 = e1*znd1 + n*j1d1
          j1d2 = e1*znd2 + n*j1d2
          zn1o = znm1
          zn1d1o = znm1d1
          zn1d2o = znm1d2
          znm1 = z2*zn1o
          znm1d1 = 2.0d0*f*zn1d1o + 2.0d0*zn1o*fd1
          znm1d2 = 2.0d0*f*zn1d2o + 4.0d0*fd1*zn1d1o + 2.0d0*zn1o*fd2
          zno = zn
          znd1o = znd1
          znd2o = znd2
          zn = z2*zno
          znd1 = 2.0d0*fd1*zno + 2.0d0*f*znd1o
          znd2 = 2.0d0*fd2*zno + 4.0d0*fd1*znd1o + 2.0d0*f*znd2o
          w = w0*w
          t0 = d(n)*w*j0
          t0d1 = d(n)*w*j0d1
          t0d2 = d(n)*w*j0d2
          w = w0*w
          t1 = d(np1)*w*j1
          t1d1 = d(np1)*w*j1d1
          t1d2 = d(np1)*w*j1d2
          sum = sum + (t0+t1)
          sumd1 = sumd1 + (t0d1+t1d1)
          sumd2 = sumd2 + (t0d2+t1d2)
          IF ((abs(t0+t1).LE.eps) .AND. (abs(t0d1+t1d1).LE.eps) .AND.
     +        (abs(t0d2+t1d2).LE.eps)) GO TO 80
   70 CONTINUE
      GO TO 90
C
   80 d1 = -fd1 + (sumd1/sum)
      d2 = -fd2 + (sumd2/sum) - ((sumd1/sum)**2.0d0)
      RETURN

   90 END
      SUBROUTINE dlfrac(a,b,x,y,lambda,eps,d1,d2)
C-----------------------------------------------------------------------
C
C     CALCULATES THE FIRST AND SECOND DERIATIVES OF LOG(IX(A,B)).
C     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
C
C     Derivatives of code in LBFRAC were taken. LBFRAC is a
C     modification of BFRAC from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION otnten
      PARAMETER (otnten=1.0D-10)
      DOUBLE PRECISION otten
      PARAMETER (otten=1.0D10)

C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda,x,y,d1,d2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,
     +                 t,w,yp1,xd1,xd2,and1,and2,bnd1,bnd2,anp1d1,
     +                 anp1d2,bnp1d1,bnp1d2,rd1,rd2,anp1o,an1d1o,an1d2o,
     +                 bnp1o,bn1d1o,bn1d2o,temp,betad1,betad2,alpad1,
     +                 alpad2,rd10,rd20
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..
C--------------------
C
      c = 1.0D0 + lambda
      c0 = b/a
      c1 = 1.0D0 + 1.0D0/a
      yp1 = y + 1.0D0
C
      xd1 = -x*y
      xd2 = xd1* (2.0d0*x-1.0d0)
C
      n = 0.0D0
      p = 1.0D0
      s = a + 1.0D0
      an = 0.0D0
      and1 = 0.0d0
      and2 = 0.0d0
      bn = 1.0D0
      bnd1 = 0.0d0
      bnd2 = 0.0d0
      anp1 = 1.0D0
      anp1d1 = 0.0d0
      anp1d2 = 0.0d0
      bnp1 = c/c1
      bnp1d1 = - (a+b)*xd1/c1
      bnp1d2 = - (a+b)*xd2/c1
      r = c1/c
      rd1 = (-r*bnp1d1+anp1d1)/bnp1
      rd2 = (-r*bnp1d2-bnp1d1*rd1+anp1d2)/bnp1 -
     +      bnp1d1* (-r*bnp1d1+anp1d1)/ (bnp1*bnp1)
      r0 = r
      rd10 = rd1
      rd20 = rd2
C
C        CONTINUED FRACTION CALCULATION
C
   10 n = n + 1.0D0
      t = n/a
      w = n* (b-n)*x
      e = a/s
      alpha = (p* (p+c0)*e*e)* (w*x)
      temp = 2.0d0* (p* (p+c0)*e*e)* (n* (b-n))
      alpad1 = temp*x*xd1
      alpad2 = temp* (x*xd2+xd1*xd1)
      e = (1.0D0+t)/ (c1+t+t)
      beta = n + w/s + e* (c+n*yp1)
      temp = n* (b-n)/s - e* (n+a+b)
      betad1 = temp*xd1
      betad2 = temp*xd2
      p = 1.0D0 + t
      s = s + 2.0D0
C
C        UPDATE AN, BN, ANP1, AND BNP1
C
      anp1o = anp1
      an1d1o = anp1d1
      an1d2o = anp1d2
      anp1 = alpha*an + beta*anp1o
      anp1d1 = alpha*and1 + alpad1*an + beta*an1d1o + betad1*anp1o
      anp1d2 = alpha*and2 + 2.0d0*alpad1*and1 + alpad2*an +
     +         beta*an1d2o + 2.0d0*betad1*an1d1o + betad2*anp1o
      an = anp1o
      and1 = an1d1o
      and2 = an1d2o
C
      bnp1o = bnp1
      bn1d1o = bnp1d1
      bn1d2o = bnp1d2
      bnp1 = alpha*bn + beta*bnp1o
      bnp1d1 = alpha*bnd1 + alpad1*bn + beta*bn1d1o + betad1*bnp1o
      bnp1d2 = alpha*bnd2 + 2.0d0*alpad1*bnd1 + alpad2*bn +
     +         beta*bn1d2o + 2.0d0*betad1*bn1d1o + betad2*bnp1o
      bn = bnp1o
      bnd1 = bn1d1o
      bnd2 = bn1d2o
C
      r = anp1/bnp1
      rd1 = (-r*bnp1d1+anp1d1)/bnp1
      rd2 = (-r*bnp1d2-bnp1d1*rd1+anp1d2)/bnp1 -
     +      bnp1d1* (-r*bnp1d1+anp1d1)/ (bnp1*bnp1)
      IF ((abs(r-r0).LE.eps*abs(r)) .AND.
     +    (abs(rd1-rd10).LE.eps*abs(rd1)) .AND.
     +    (abs(rd2-rd20).LE.eps*abs(rd2))) GO TO 70
C
C        RESCALE AN, BN, ANP1, AND BNP1 and DERIVATIVES
C
      r0 = r
      rd10 = rd1
      rd20 = rd2
C
   20 IF (.NOT. (abs(anp1).GT.otten)) GO TO 30
      anp1 = anp1*otnten
      anp1d1 = anp1d1*otnten
      anp1d2 = anp1d2*otnten
      an = an*otnten
      and1 = and1*otnten
      and2 = and2*otnten
      r0 = r0*otnten
      rd10 = rd10*otnten
      rd20 = rd20*otnten
      GO TO 20

   30 IF (.NOT. (abs(anp1).LT.otnten)) GO TO 40
      anp1 = anp1*otten
      anp1d1 = anp1d1*otten
      anp1d2 = anp1d2*otten
      an = an*otten
      and1 = and1*otten
      and2 = and2*otten
      r0 = r0*otten
      rd10 = rd10*otten
      rd20 = rd20*otten
      GO TO 30

   40 IF (.NOT. (abs(bnp1).GT.otten)) GO TO 50
      bnp1 = bnp1*otnten
      bnp1d1 = bnp1d1*otnten
      bnp1d2 = bnp1d2*otnten
      bn = bn*otnten
      bnd1 = bnd1*otnten
      bnd2 = bnd2*otnten
      r0 = r0*otten
      rd10 = rd10*otten
      rd20 = rd20*otten
      GO TO 40

   50 IF (.NOT. (abs(bnp1).LT.otnten)) GO TO 60
      bnp1 = bnp1*otten
      bnp1d1 = bnp1d1*otten
      bnp1d2 = bnp1d2*otten
      bn = bn*otten
      bnd1 = bnd1*otten
      bnd2 = bnd2*otten
      r0 = r0*otnten
      rd10 = rd10*otnten
      rd20 = rd20*otnten
      GO TO 50

   60 GO TO 10
C
C                 TERMINATION
C
   70 CONTINUE
      d1 = rd1/r
      d2 = - (d1*d1) + rd2/r
      RETURN

      END
      DOUBLE PRECISION FUNCTION dln1pe(x)
      IMPLICIT DOUBLE PRECISION (a-h,o-p,r-z),INTEGER (i-n),LOGICAL (q)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DLN1PE( X )
C          Double Precision LN(1 + exp(x))
C
C
C                              Arguments
C
C
C     X --> Argument
C                    DOUBLE PRECISION X
C
C**********************************************************************
C     LBREAK is log(0.375)
C     .. Parameters ..
      DOUBLE PRECISION lbreak
      PARAMETER (lbreak=-.98082925301172623686D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dln1px
      EXTERNAL dln1px
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,log
C     ..
C     .. Executable Statements ..

      IF (x.LE.lbreak) THEN
          dln1pe = dln1px(exp(x))

      ELSE IF (-x.LE.lbreak) THEN
          dln1pe = x + dln1px(exp(-x))

      ELSE
          dln1pe = log(one+exp(x))
      END IF

      RETURN

      END
      DOUBLE PRECISION FUNCTION dln1px(a)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DLN1PX(X)
C               Double precision LN(1+X)
C
C
C                              Function
C
C
C     Returns ln(1+x)
C     Note that the obvious code of
C               LOG(1.0+X)
C     won't work for small X because 1.0+X loses accuracy
C
C
C                              Arguments
C
C
C     X --> Value for which ln(1-x) is desired.
C                                        X is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Renames ALNREL from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C**********************************************************************
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION LN(1 + A)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,p3,q1,q2,q3,t,t2,w,x
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog
C     ..
C     .. Data statements ..
      DATA p1/-.129418923021993D+01/,p2/.405303492862024D+00/,
     +     p3/-.178874546012214D-01/
      DATA q1/-.162752256355323D+01/,q2/.747811014037616D+00/,
     +     q3/-.845104217945565D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      IF (abs(a).GT.0.375D0) GO TO 10
      t = a/ (a+2.0D0)
      t2 = t*t
      w = (((p3*t2+p2)*t2+p1)*t2+1.0D0)/ (((q3*t2+q2)*t2+q1)*t2+1.0D0)
      dln1px = 2.0D0*t*w
      RETURN
C
   10 x = 1.D0 + dble(a)
      dln1px = dlog(x)
      RETURN

      END
      DOUBLE PRECISION FUNCTION dlngam(a)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DLNGAM(X)
C                 Double precision LN of the GAMma function
C
C
C                              Function
C
C
C     Returns the natural logarithm of GAMMA(X).
C
C
C                              Arguments
C
C
C     X --> value at which scaled log gamma is to be returned
C                    X is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Renames GAMLN from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C**********************************************************************
C-----------------------------------------------------------------------
C            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS
C          NAVAL SURFACE WARFARE CENTER
C          DAHLGREN, VIRGINIA
C--------------------------
C     D = 0.5*(LN(2*PI) - 1)
C--------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION c0,c1,c2,c3,c4,c5,d,t,w
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION gamln1
      EXTERNAL gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog
C     ..
C     .. Data statements ..
C--------------------------
      DATA d/.418938533204673D0/
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
      IF (a.GT.0.8D0) GO TO 10
      dlngam = gamln1(a) - dlog(a)
      RETURN

   10 IF (a.GT.2.25D0) GO TO 20
      t = (a-0.5D0) - 0.5D0
      dlngam = gamln1(t)
      RETURN
C
   20 IF (a.GE.10.0D0) GO TO 40
      n = a - 1.25D0
      t = a
      w = 1.0D0
      DO 30 i = 1,n
          t = t - 1.0D0
          w = t*w
   30 CONTINUE
      dlngam = gamln1(t-1.0D0) + dlog(w)
      RETURN
C
   40 t = (1.0D0/a)**2
      w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a
      dlngam = (d+w) + (a-0.5D0)* (dlog(a)-1.0D0)
      END
      SUBROUTINE dser(a,b,z,d1,d2)
C-------------------------------------------------------------------
C
C    Calculates first and second derivatives of log of
C       1 + [Beta(A+1,n+1)/Beta(B-n-1,n+1)]*(x/(1-x))^(n+1)] n=0,1,...
C    This summation is part of Abramowitz and Stegun series 26.5.5
C       [(x^b)*((1-x)^(a-1))/(b*Beta(a,b))]*SUM
C    which can be used to calculate IX(A,B).
C
C-------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=1.0D-13)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,z,d1,d2
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION zd1,zd2,am2,bp2,term,termd1,termd2,termo,trmd1o,
     +                 trmd2o,sum,sumd1,sumd2,sump1,xn,ratio,ratd1,
     +                 ratd2,temp
      LOGICAL qdone
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dln1px
      EXTERNAL dln1px
C     ..
C     .. Executable Statements ..
C
C     initialize variables
C
C
      qdone = .FALSE.
      zd1 = -z
      zd2 = z
      am2 = a - one - one
      bp2 = b + one + one
      temp = (a-one)/ (b+one)
      term = temp*z
      termd1 = temp*zd1
      termd2 = temp*zd2
      sum = zero
      sumd1 = zero
      sumd2 = zero
      xn = zero
C
C     Compute sum
C
      GO TO 20

   10 IF (qdone .OR. (xn.GT.100D0)) GO TO 30
   20 sum = sum + term
      sumd1 = sumd1 + termd1
      sumd2 = sumd2 + termd2
      qdone = ((abs(term/sum).LE.tiny) .AND.
     +        (abs(termd1/sumd1).LE.tiny) .AND.
     +        (abs(termd2/sumd2).LE.tiny))
      IF (.NOT.qdone) THEN
          temp = (am2-xn)/ (bp2+xn)
          ratio = temp*z
          ratd1 = temp*zd1
          ratd2 = temp*zd2
          termo = term
          trmd1o = termd1
          trmd2o = termd2
          term = termo*ratio
          termd1 = trmd1o*ratio + termo*ratd1
          termd2 = trmd2o*ratio + 2.0D0*trmd1o*ratd1 + termo*ratd2
          xn = xn + one
      END IF

      GO TO 10

   30 sump1 = sum + one

      d1 = sumd1/sump1
      d2 = - (d1*d1) + sumd2/sump1

      RETURN

      END
      DOUBLE PRECISION FUNCTION dstrem(z)
      IMPLICIT DOUBLE PRECISION (a-h,o-p,r-z),INTEGER (i-n),LOGICAL (q)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DSTREM( Z )
C             Double precision Sterling Remainder
C
C
C                              Function
C
C
C     Returns   Log(Gamma(Z))  -  Sterling(Z)  where   Sterling(Z)  is
C     Sterling's Approximation to Log(Gamma(Z))
C
C     Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
C
C
C                              Arguments
C
C
C     Z --> Value at which Sterling remainder calculated
C           Must be positive.
C                  DOUBLE PRECISION Z
C
C
C                              Method
C
C
C
C     If Z >= 6 uses 9 terms of series in Bernoulli numbers
C     (Values calculated using Maple)
C     Otherwise computes difference explicitly
C
C**********************************************************************

C     .. Parameters ..
      DOUBLE PRECISION hln2pi
      PARAMETER (hln2pi=0.91893853320467274178D0)
      INTEGER ncoef
      PARAMETER (ncoef=10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sterl
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION coef(ncoef)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl,dlngam
      EXTERNAL devlpl,dlngam
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC log
C     ..
C     .. Data statements ..
      DATA coef/0.0D0,0.0833333333333333333333333333333D0,
     +     -0.00277777777777777777777777777778D0,
     +     0.000793650793650793650793650793651D0,
     +     -0.000595238095238095238095238095238D0,
     +     0.000841750841750841750841750841751D0,
     +     -0.00191752691752691752691752691753D0,
     +     0.00641025641025641025641025641026D0,
     +     -0.0295506535947712418300653594771D0,
     +     0.179644372368830573164938490016D0/
C     ..
C     .. Executable Statements ..

C    For information, here are the next 11 coefficients of the
C    remainder term in Sterling's formula
C            -1.39243221690590111642743221691
C            13.4028640441683919944789510007
C            -156.848284626002017306365132452
C            2193.10333333333333333333333333
C            -36108.7712537249893571732652192
C            691472.268851313067108395250776
C            -0.152382215394074161922833649589D8
C            0.382900751391414141414141414141D9
C            -0.108822660357843910890151491655D11
C            0.347320283765002252252252252252D12
C            -0.123696021422692744542517103493D14
C

      IF (z.LE.0.0D0) STOP 'Zero or negative argument in DSTREM'
      IF (.NOT. (z.GT.6.0D0)) GO TO 10
      dstrem = devlpl(coef,10,1.0D0/z**2)*z
      GO TO 20

   10 sterl = hln2pi + (z-0.5D0)*log(z) - z
      dstrem = dlngam(z) - sterl
   20 RETURN

      END
      DOUBLE PRECISION FUNCTION erf(x)
C-----------------------------------------------------------------------
C             EVALUATION OF THE REAL ERROR FUNCTION
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ax,bot,c,t,top,x2
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sign
C     ..
C     .. Data statements ..
C-------------------------
C-------------------------
C-------------------------
C-------------------------
      DATA c/.564189583547756D0/
      DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +     a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +     a(5)/.128379167095513D+00/
      DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +     b(3)/.375795757275549D+00/
      DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +     p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +     p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +     p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
      DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +     q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +     q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +     q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
      DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +     r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +     r(5)/2.82094791773523D-01/
      DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +     s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
C     ..
C     .. Executable Statements ..
C-------------------------
      ax = abs(x)
      IF (ax.GT.0.5D0) GO TO 10
      t = x*x
      top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
      bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
      erf = x* (top/bot)
      RETURN
C
   10 IF (ax.GT.4.0D0) GO TO 20
      top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+
     +      p(7))*ax + p(8)
      bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+
     +      q(7))*ax + q(8)
      erf = 0.5D0 + (0.5D0-exp(-x*x)*top/bot)
      IF (x.LT.0.0D0) erf = -erf
      RETURN
C
   20 IF (ax.GE.5.8D0) GO TO 30
      x2 = x*x
      t = 1.0D0/x2
      top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
      bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
      erf = (c-top/ (x2*bot))/ax
      erf = 0.5D0 + (0.5D0-exp(-x2)*erf)
      IF (x.LT.0.0D0) erf = -erf
      RETURN
C
   30 erf = sign(1.0D0,x)
      RETURN

      END
      DOUBLE PRECISION FUNCTION erfc1(ind,x)
C-----------------------------------------------------------------------
C         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
C
C          ERFC1(IND,X) = ERFC(X)            IF IND = 0
C          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER ind
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ax,bot,c,e,t,top,w
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg
      EXTERNAL exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,exp
C     ..
C     .. Data statements ..
C-------------------------
C-------------------------
C-------------------------
C-------------------------
      DATA c/.564189583547756D0/
      DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +     a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +     a(5)/.128379167095513D+00/
      DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +     b(3)/.375795757275549D+00/
      DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +     p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +     p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +     p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
      DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +     q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +     q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +     q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
      DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +     r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +     r(5)/2.82094791773523D-01/
      DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +     s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
C     ..
C     .. Executable Statements ..
C-------------------------
C
C                     ABS(X) .LE. 0.5
C
      ax = abs(x)
      IF (ax.GT.0.5D0) GO TO 10
      t = x*x
      top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
      bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
      erfc1 = 0.5D0 + (0.5D0-x* (top/bot))
      IF (ind.NE.0) erfc1 = exp(t)*erfc1
      RETURN
C
C                  0.5 .LT. ABS(X) .LE. 4
C
   10 IF (ax.GT.4.0D0) GO TO 20
      top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+p(6))*ax+
     +      p(7))*ax + p(8)
      bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+q(6))*ax+
     +      q(7))*ax + q(8)
      erfc1 = top/bot
      GO TO 40
C
C                      ABS(X) .GT. 4
C
   20 IF (x.LE.-5.6D0) GO TO 60
      IF (ind.NE.0) GO TO 30
      IF (x.GT.100.0D0) GO TO 70
      IF (x*x.GT.-exparg(1)) GO TO 70
C
   30 t = (1.0D0/x)**2
      top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
      bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
      erfc1 = (c-t*top/bot)/ax
C
C                      FINAL ASSEMBLY
C
   40 IF (ind.EQ.0) GO TO 50
      IF (x.LT.0.0D0) erfc1 = 2.0D0*exp(x*x) - erfc1
      RETURN

   50 w = dble(x)*dble(x)
      t = w
      e = w - dble(t)
      erfc1 = ((0.5D0+ (0.5D0-e))*exp(-t))*erfc1
      IF (x.LT.0.0D0) erfc1 = 2.0D0 - erfc1
      RETURN
C
C             LIMIT VALUE FOR LARGE NEGATIVE X
C
   60 erfc1 = 2.0D0
      IF (ind.NE.0) erfc1 = 2.0D0*exp(x*x)
      RETURN
C
C             LIMIT VALUE FOR LARGE POSITIVE X
C                       WHEN IND = 0
C
   70 erfc1 = 0.0D0
      RETURN

      END
      DOUBLE PRECISION FUNCTION esum(mu,x)
C-----------------------------------------------------------------------
C                    EVALUATION OF EXP(MU + X)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER mu
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp
C     ..
C     .. Executable Statements ..

      IF (x.GT.0.0D0) GO TO 10
C
      IF (mu.LT.0) GO TO 20
      w = mu + x
      IF (w.GT.0.0D0) GO TO 20
      esum = exp(w)
      RETURN
C
   10 IF (mu.GT.0) GO TO 20
      w = mu + x
      IF (w.LT.0.0D0) GO TO 20
      esum = exp(w)
      RETURN
C
   20 w = mu
      esum = exp(w)*exp(x)
      RETURN

      END
      DOUBLE PRECISION FUNCTION exparg(l)
C--------------------------------------------------------------------
C     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
C     EXP(W) CAN BE COMPUTED.
C
C     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
C     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
C
C     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
C--------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER l
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION lnb
      INTEGER b,m
C     ..
C     .. External Functions ..
      INTEGER ipmpar
      EXTERNAL ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Executable Statements ..
C
      b = ipmpar(4)
      IF (b.NE.2) GO TO 10
      lnb = .69314718055995D0
      GO TO 40

   10 IF (b.NE.8) GO TO 20
      lnb = 2.0794415416798D0
      GO TO 40

   20 IF (b.NE.16) GO TO 30
      lnb = 2.7725887222398D0
      GO TO 40

   30 lnb = dlog(dble(b))
C
   40 IF (l.EQ.0) GO TO 50
      m = ipmpar(9) - 1
      exparg = 0.99999D0* (m*lnb)
      RETURN

   50 m = ipmpar(10)
      exparg = 0.99999D0* (m*lnb)
      RETURN

      END
      DOUBLE PRECISION FUNCTION fpser(a,b,x,eps)
C-----------------------------------------------------------------------
C
C                 EVALUATION OF I (A,B)
C                                X
C
C          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
C
C-----------------------------------------------------------------------
C
C                  SET  FPSER = X**A
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION an,c,s,t,tol
C     ..
C     .. External Functions ..
      DOUBLE PRECISION exparg
      EXTERNAL exparg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp
C     ..
C     .. Executable Statements ..

      fpser = 1.0D0
      IF (a.LE.1.D-3*eps) GO TO 10
      fpser = 0.0D0
      t = a*dlog(x)
      IF (t.LT.exparg(1)) RETURN
      fpser = exp(t)
C
C                NOTE THAT 1/B(A,B) = B
C
   10 fpser = (b/a)*fpser
      tol = eps/a
      an = a + 1.0D0
      t = x
      s = t/an
   20 an = an + 1.0D0
      t = x*t
      c = t/an
      s = s + c
      IF (abs(c).GT.tol) GO TO 20
C
      fpser = fpser* (1.0D0+a*s)
      RETURN

      END
      DOUBLE PRECISION FUNCTION gam1(a)
C     ------------------------------------------------------------------
C     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
C     ------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bot,d,s1,s2,t,top,w
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p(7),q(5),r(9)
C     ..
C     .. Data statements ..
C     -------------------
C     -------------------
C     -------------------
C     -------------------
      DATA p(1)/.577215664901533D+00/,p(2)/-.409078193005776D+00/,
     +     p(3)/-.230975380857675D+00/,p(4)/.597275330452234D-01/,
     +     p(5)/.766968181649490D-02/,p(6)/-.514889771323592D-02/,
     +     p(7)/.589597428611429D-03/
      DATA q(1)/.100000000000000D+01/,q(2)/.427569613095214D+00/,
     +     q(3)/.158451672430138D+00/,q(4)/.261132021441447D-01/,
     +     q(5)/.423244297896961D-02/
      DATA r(1)/-.422784335098468D+00/,r(2)/-.771330383816272D+00/,
     +     r(3)/-.244757765222226D+00/,r(4)/.118378989872749D+00/,
     +     r(5)/.930357293360349D-03/,r(6)/-.118290993445146D-01/,
     +     r(7)/.223047661158249D-02/,r(8)/.266505979058923D-03/,
     +     r(9)/-.132674909766242D-03/
      DATA s1/.273076135303957D+00/,s2/.559398236957378D-01/
C     ..
C     .. Executable Statements ..
C     -------------------
      t = a
      d = a - 0.5D0
      IF (d.GT.0.0D0) t = d - 0.5D0
      IF (t) 40,10,20
C
   10 gam1 = 0.0D0
      RETURN
C
   20 top = (((((p(7)*t+p(6))*t+p(5))*t+p(4))*t+p(3))*t+p(2))*t + p(1)
      bot = (((q(5)*t+q(4))*t+q(3))*t+q(2))*t + 1.0D0
      w = top/bot
      IF (d.GT.0.0D0) GO TO 30
      gam1 = a*w
      RETURN

   30 gam1 = (t/a)* ((w-0.5D0)-0.5D0)
      RETURN
C
   40 top = (((((((r(9)*t+r(8))*t+r(7))*t+r(6))*t+r(5))*t+r(4))*t+r(3))*
     +      t+r(2))*t + r(1)
      bot = (s2*t+s1)*t + 1.0D0
      w = top/bot
      IF (d.GT.0.0D0) GO TO 50
      gam1 = a* ((w+0.5D0)+0.5D0)
      RETURN

   50 gam1 = t*w/a
      RETURN

      END
      DOUBLE PRECISION FUNCTION gamln(a)
C-----------------------------------------------------------------------
C            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS
C          NAVAL SURFACE WARFARE CENTER
C          DAHLGREN, VIRGINIA
C--------------------------
C     D = 0.5*(LN(2*PI) - 1)
C--------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION c0,c1,c2,c3,c4,c5,d,t,w
      INTEGER i,n
C     ..
C     .. External Functions ..
      DOUBLE PRECISION gamln1
      EXTERNAL gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dlog
C     ..
C     .. Data statements ..
C--------------------------
      DATA d/.418938533204673D0/
      DATA c0/.833333333333333D-01/,c1/-.277777777760991D-02/,
     +     c2/.793650666825390D-03/,c3/-.595202931351870D-03/,
     +     c4/.837308034031215D-03/,c5/-.165322962780713D-02/
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
      IF (a.GT.0.8D0) GO TO 10
      gamln = gamln1(a) - dlog(a)
      RETURN

   10 IF (a.GT.2.25D0) GO TO 20
      t = (a-0.5D0) - 0.5D0
      gamln = gamln1(t)
      RETURN
C
   20 IF (a.GE.10.0D0) GO TO 40
      n = a - 1.25D0
      t = a
      w = 1.0D0
      DO 30 i = 1,n
          t = t - 1.0D0
          w = t*w
   30 CONTINUE
      gamln = gamln1(t-1.0D0) + dlog(w)
      RETURN
C
   40 t = (1.0D0/a)**2
      w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a
      gamln = (d+w) + (a-0.5D0)* (dlog(a)-1.0D0)
      END
      DOUBLE PRECISION FUNCTION gamln1(a)
C-----------------------------------------------------------------------
C     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p0,p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6,r0,r1,r2,
     +                 r3,r4,r5,s1,s2,s3,s4,s5,w,x
C     ..
C     .. Data statements ..
C----------------------
      DATA p0/.577215664901533D+00/,p1/.844203922187225D+00/,
     +     p2/-.168860593646662D+00/,p3/-.780427615533591D+00/,
     +     p4/-.402055799310489D+00/,p5/-.673562214325671D-01/,
     +     p6/-.271935708322958D-02/
      DATA q1/.288743195473681D+01/,q2/.312755088914843D+01/,
     +     q3/.156875193295039D+01/,q4/.361951990101499D+00/,
     +     q5/.325038868253937D-01/,q6/.667465618796164D-03/
      DATA r0/.422784335098467D+00/,r1/.848044614534529D+00/,
     +     r2/.565221050691933D+00/,r3/.156513060486551D+00/,
     +     r4/.170502484022650D-01/,r5/.497958207639485D-03/
      DATA s1/.124313399877507D+01/,s2/.548042109832463D+00/,
     +     s3/.101552187439830D+00/,s4/.713309612391000D-02/,
     +     s5/.116165475989616D-03/
C     ..
C     .. Executable Statements ..
C----------------------
      IF (a.GE.0.6D0) GO TO 10
      w = ((((((p6*a+p5)*a+p4)*a+p3)*a+p2)*a+p1)*a+p0)/
     +    ((((((q6*a+q5)*a+q4)*a+q3)*a+q2)*a+q1)*a+1.0D0)
      gamln1 = -a*w
      RETURN
C
   10 x = (a-0.5D0) - 0.5D0
      w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/
     +    (((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x+1.0D0)
      gamln1 = x*w
      RETURN

      END
      SUBROUTINE grat1(a,x,r,p,q,eps)
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,eps,p,q,r,x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,
     +                 t,tol,w,z
C     ..
C     .. External Functions ..
      DOUBLE PRECISION erf,erfc1,gam1,rexp
      EXTERNAL erf,erfc1,gam1,rexp
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dlog,exp,sqrt
C     ..
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
C                      P(A,X) AND Q(A,X)
C
C     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
C     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
C-----------------------------------------------------------------------
      IF (a*x.EQ.0.0D0) GO TO 120
      IF (a.EQ.0.5D0) GO TO 100
      IF (x.LT.1.1D0) GO TO 10
      GO TO 60
C
C             TAYLOR SERIES FOR P(A,X)/X**A
C
   10 an = 3.0D0
      c = x
      sum = x/ (a+3.0D0)
      tol = 0.1D0*eps/ (a+1.0D0)
   20 an = an + 1.0D0
      c = -c* (x/an)
      t = c/ (a+an)
      sum = sum + t
      IF (abs(t).GT.tol) GO TO 20
      j = a*x* ((sum/6.0D0-0.5D0/ (a+2.0D0))*x+1.0D0/ (a+1.0D0))
C
      z = a*dlog(x)
      h = gam1(a)
      g = 1.0D0 + h
      IF (x.LT.0.25D0) GO TO 30
      IF (a.LT.x/2.59D0) GO TO 50
      GO TO 40

   30 IF (z.GT.-.13394D0) GO TO 50
C
   40 w = exp(z)
      p = w*g* (0.5D0+ (0.5D0-j))
      q = 0.5D0 + (0.5D0-p)
      RETURN
C
   50 l = rexp(z)
      w = 0.5D0 + (0.5D0+l)
      q = (w*j-l)*g - h
      IF (q.LT.0.0D0) GO TO 90
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
C              CONTINUED FRACTION EXPANSION
C
   60 a2nm1 = 1.0D0
      a2n = 1.0D0
      b2nm1 = x
      b2n = x + (1.0D0-a)
      c = 1.0D0
   70 a2nm1 = x*a2n + c*a2nm1
      b2nm1 = x*b2n + c*b2nm1
      am0 = a2nm1/b2nm1
      c = c + 1.0D0
      cma = c - a
      a2n = a2nm1 + cma*a2n
      b2n = b2nm1 + cma*b2n
      an0 = a2n/b2n
      IF (abs(an0-am0).GE.eps*an0) GO TO 70
      q = r*an0
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
C                SPECIAL CASES
C
   80 p = 0.0D0
      q = 1.0D0
      RETURN
C
   90 p = 1.0D0
      q = 0.0D0
      RETURN
C
  100 IF (x.GE.0.25D0) GO TO 110
      p = erf(sqrt(x))
      q = 0.5D0 + (0.5D0-p)
      RETURN

  110 q = erfc1(0,sqrt(x))
      p = 0.5D0 + (0.5D0-q)
      RETURN
C
  120 IF (x.LE.a) GO TO 80
      GO TO 90

      END
      DOUBLE PRECISION FUNCTION gsumln(a,b)
C-----------------------------------------------------------------------
C          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
C          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION x
C     ..
C     .. External Functions ..
      DOUBLE PRECISION alnrel,gamln1
      EXTERNAL alnrel,gamln1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Executable Statements ..
      x = dble(a) + dble(b) - 2.D0
      IF (x.GT.0.25D0) GO TO 10
      gsumln = gamln1(1.0D0+x)
      RETURN

   10 IF (x.GT.1.25D0) GO TO 20
      gsumln = gamln1(x) + alnrel(x)
      RETURN

   20 gsumln = gamln1(x-1.0D0) + dlog(x* (1.0D0+x))
      RETURN

      END
      INTEGER FUNCTION ipmpar(i)
C-----------------------------------------------------------------------
C
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C-----------------------------------------------------------------------
C
C     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
C     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
C     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
C     COLUMN 1.)
C
C-----------------------------------------------------------------------
C
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER i
C     ..
C     .. Local Arrays ..
      INTEGER imach(10)
C     ..
C     .. Data statements ..
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C     DATA IMACH( 1) /   2 /
C     DATA IMACH( 2) /  31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /  16 /
C     DATA IMACH( 5) /   6 /
C     DATA IMACH( 6) / -64 /
C     DATA IMACH( 7) /  63 /
C     DATA IMACH( 8) /  14 /
C     DATA IMACH( 9) / -64 /
C     DATA IMACH(10) /  63 /
C
C     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
C     PC 7300, AND AT&T 6300.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   33 /
C     DATA IMACH( 3) / 8589934591 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -256 /
C     DATA IMACH( 7) /  255 /
C     DATA IMACH( 8) /   60 /
C     DATA IMACH( 9) / -256 /
C     DATA IMACH(10) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /    8 /
C     DATA IMACH( 5) /   13 /
C     DATA IMACH( 6) /  -50 /
C     DATA IMACH( 7) /   76 /
C     DATA IMACH( 8) /   26 /
C     DATA IMACH( 9) /  -50 /
C     DATA IMACH(10) /   76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /      2 /
C     DATA IMACH( 2) /     39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /      8 /
C     DATA IMACH( 5) /     13 /
C     DATA IMACH( 6) /    -50 /
C     DATA IMACH( 7) /     76 /
C     DATA IMACH( 8) /     26 /
C     DATA IMACH( 9) / -32754 /
C     DATA IMACH(10) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS OPERATING SYSTEM).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   48 /
C     DATA IMACH( 3) / 281474976710655 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   48 /
C     DATA IMACH( 6) / -974 /
C     DATA IMACH( 7) / 1070 /
C     DATA IMACH( 8) /   95 /
C     DATA IMACH( 9) / -926 /
C     DATA IMACH(10) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS/VE OPERATING SYSTEM).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -4096 /
C     DATA IMACH( 7) /  4095 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -4096 /
C     DATA IMACH(10) /  4095 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    47 /
C     DATA IMACH( 6) / -8189 /
C     DATA IMACH( 7) /  8190 /
C     DATA IMACH( 8) /    94 /
C     DATA IMACH( 9) / -8099 /
C     DATA IMACH(10) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   23 /
C     DATA IMACH( 3) / 8388607 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   38 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
C     AND DPS 8/70 SERIES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   63 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   39 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   55 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -126 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
C     5/7/9 AND THE SEL SYSTEMS 85/86.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC.
C
C      DATA imach(1)/2/
C      DATA imach(2)/31/
C      DATA imach(3)/2147483647/
C      DATA imach(4)/2/
C      DATA imach(5)/24/
C      DATA imach(6)/-125/
C      DATA imach(7)/128/
C      DATA imach(8)/53/
C      DATA imach(9)/-1021/
C      DATA imach(10)/1024/
C
C     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
C     MACFORTRAN II.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   54 /
C     DATA IMACH( 9) / -101 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   62 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
C     SERIES (MIPS R3000 PROCESSOR).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   60 /
C     DATA IMACH( 9) /-1024 /
C     DATA IMACH(10) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
      ipmpar = imach(i)
      RETURN

      END
      DOUBLE PRECISION FUNCTION lasym(a,b,lambda,eps)
C-----------------------------------------------------------------------
C
C     ASYMPTOTIC EXPANSION FOR LOG(IX(A,B)) FOR LARGE A AND B.
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C
C     Modification of BASYM from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C-----------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION bsum,dsum,e0,e1,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t0,
     +                 t1,u,w,w0,z,z0,z2,zn,znm1
      INTEGER i,im1,imj,j,m,mm1,mmj,n,np1,num
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a0(21),b0(21),c(21),d(21)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION bcorr,erfc1,rlog1
      EXTERNAL bcorr,erfc1,rlog1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,sqrt,log
C     ..
C     .. Data statements ..
C------------------------
C     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
C            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
C------------------------
C     E0 = 2/SQRT(PI)
C     E1 = 2**(-3/2)
C------------------------
      DATA num/20/
      DATA e0/1.12837916709551D0/,e1/.353553390593274D0/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (a.GE.b) GO TO 10
      h = a/b
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/b
      w0 = 1.0D0/sqrt(a* (1.0D0+h))
      GO TO 20

   10 h = b/a
      r0 = 1.0D0/ (1.0D0+h)
      r1 = (b-a)/a
      w0 = 1.0D0/sqrt(b* (1.0D0+h))
C
   20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      z0 = sqrt(f)
      z = 0.5D0* (z0/e1)
      z2 = f + f
C
      a0(1) = (2.0D0/3.0D0)*r1
      c(1) = -0.5D0*a0(1)
      d(1) = -c(1)
      j0 = (0.5D0/e0)*erfc1(1,z0)
      j1 = e1
      sum = j0 + d(1)*w0*j1
C
      s = 1.0D0
      h2 = h*h
      hn = 1.0D0
      w = w0
      znm1 = z
      zn = z2
      DO 70 n = 2,num,2
          hn = h2*hn
          a0(n) = 2.0D0*r0* (1.0D0+h*hn)/ (n+2.0D0)
          np1 = n + 1
          s = s + hn
          a0(np1) = 2.0D0*r1*s/ (n+3.0D0)
C
          DO 60 i = n,np1
              r = -0.5D0* (i+1.0D0)
              b0(1) = r*a0(1)
              DO 40 m = 2,i
                  bsum = 0.0D0
                  mm1 = m - 1
                  DO 30 j = 1,mm1
                      mmj = m - j
                      bsum = bsum + (j*r-mmj)*a0(j)*b0(mmj)
   30             CONTINUE
                  b0(m) = r*a0(m) + bsum/m
   40         CONTINUE
              c(i) = b0(i)/ (i+1.0D0)
C
              dsum = 0.0D0
              im1 = i - 1
              DO 50 j = 1,im1
                  imj = i - j
                  dsum = dsum + d(imj)*c(j)
   50         CONTINUE
              d(i) = - (dsum+c(i))
   60     CONTINUE
C
          j0 = e1*znm1 + (n-1.0D0)*j0
          j1 = e1*zn + n*j1
          znm1 = z2*znm1
          zn = z2*zn
          w = w0*w
          t0 = d(n)*w*j0
          w = w0*w
          t1 = d(np1)*w*j1
          sum = sum + (t0+t1)
          IF (abs(t0+t1).LE.eps) GO TO 80
   70 CONTINUE
      GO TO 90
C
   80 u = -bcorr(a,b)
      lasym = log(e0) - f + u + log(sum)
      RETURN

   90 END
      DOUBLE PRECISION FUNCTION lbfrac(a,b,x,y,lambda,eps)
C-----------------------------------------------------------------------
C
C     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
C
C     Modification of BFRAC from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION otnten
      PARAMETER (otnten=1.0D-10)
      DOUBLE PRECISION otten
      PARAMETER (otten=1.0D10)
      DOUBLE PRECISION lotten
      PARAMETER (lotten=23.0258509299404568402D0)

C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,eps,lambda,x,y
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,
     +                 t,w,yp1
      INTEGER scala,scalb
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..
C--------------------
      scala = 0
      scalb = 0
C
      c = 1.0D0 + lambda
      c0 = b/a
      c1 = 1.0D0 + 1.0D0/a
      yp1 = y + 1.0D0
C
      n = 0.0D0
      p = 1.0D0
      s = a + 1.0D0
      an = 0.0D0
      bn = 1.0D0
      anp1 = 1.0D0
      bnp1 = c/c1
      r = c1/c
      r0 = r
C
C        CONTINUED FRACTION CALCULATION
C
   10 n = n + 1.0D0
      t = n/a
      w = n* (b-n)*x
      e = a/s
      alpha = (p* (p+c0)*e*e)* (w*x)
      e = (1.0D0+t)/ (c1+t+t)
      beta = n + w/s + e* (c+n*yp1)
      p = 1.0D0 + t
      s = s + 2.0D0
C
C        UPDATE AN, BN, ANP1, AND BNP1
C
      t = alpha*an + beta*anp1
      an = anp1
      anp1 = t
      t = alpha*bn + beta*bnp1
      bn = bnp1
      bnp1 = t
C
      r = anp1/bnp1
      IF (abs(r-r0).LE.eps*r) GO TO 70
C
C        RESCALE AN, BN, ANP1, AND BNP1
C
      r0 = r
   20 IF (.NOT. (abs(anp1).GT.otten)) GO TO 30
      anp1 = anp1*otnten
      an = an*otnten
      r0 = r0*otnten
      scala = scala + 1
      GO TO 20

   30 IF (.NOT. (abs(anp1).LT.otnten)) GO TO 40
      anp1 = anp1*otten
      an = an*otten
      r0 = r0*otten
      scala = scala - 1
      GO TO 30

   40 IF (.NOT. (abs(bnp1).GT.otten)) GO TO 50
      bnp1 = bnp1*otnten
      bn = bn*otnten
      r0 = r0*otten
      scalb = scalb + 1
      GO TO 40

   50 IF (.NOT. (abs(bnp1).LT.otnten)) GO TO 60
      bnp1 = bnp1*otten
      bn = bn*otten
      r0 = r0*otnten
      scalb = scalb - 1
      GO TO 50

   60 GO TO 10
C
C                 TERMINATION
C
   70 CONTINUE
      lbfrac = log(r) + (scala-scalb)*lotten
      RETURN

      END
      SUBROUTINE ltlf(w,dfn,dfd,lcum,ltail)
C**********************************************************************
C
C      SUBROUTINE LTLF( W, DFN, DFD, LCUM, LTAIL)
C           Logarithms of the Tails of the Log-F Distribution
C
C
C                              Function
C
C
C      Let F(Z|DFN,DFD) be the cumulative F distribution with degrees
C      of freedom DFN and DFD. Let W have a Log-F distribtuion with
C      the same degrees of freedom. That is, exp(W) is distributed
C      as F.
C
C      This routine calculates the cumulative and tail of the log-F
C            LCUM = log(F(exp(W)|DFN,DFD))
C            LTAIL = log(1 - F(exp(W)|DFN,DFD))
C
C
C                              Arguments
C
C
C      W --> Argument at which the tails of the Log-F are to
C            be evaluated
C                DOUBLE PRECISION W
C
C      DFN --> Numerator degrees of freedom of the F and Log-F
C              distribution
C                DOUBLE PRECISION DFN
C
C      DFD --> Denominator degrees of freedom of the F and Log-F
C              distribution
C                DOUBLE PRECISION DFD
C
C      LCUM <-- Cumulative of Log-F and -infinity to W, i.e.,
C               LOG(F(EXP(W)|DFN,DFD))
C                DOUBLE PRECISION LCUM
C
C      LTAIL <-- Tail of Log-F from WW to infinity, i.e.,
C                LOG(1-F(EXP(W)|DFN,DFD))
C                DOUBLE PRECISION LTAIL
C
C
C                              Method
C
C
C      The overall strategy is to compute accurately the smaller of
C      LCUM and or LTAIL and to derive the other from it. Five
C      steps are taken in this calculation.
C
C      (1) See if the argument to the F is small enough that the
C      associated incomplete beta density is well approximated
C      by a polynomial. This is done for both (W,DFN,DFD) and
C      (-W,DFD,DFN). If W cannot be exponentiated, this case
C      will return the answer.
C
C      (2) Use BRATIO to return both the cumulative and tail F.
C      Accept the answer unless the smaller of the two numbers is
C      zero.
C
C      (3) If BRATIO returns 0 try asymptotic calculation in
C      LAYMP. If a and b are large and x or y is approximately
C      the mean this method is used.
C
C      (4) If methods 1 to 3 fail try using the Continued Fraction
C      approximation to the incomplete beta.
C
C      (5) If methods 1 to 4 fail use series 26.6.5 of Abramowitz
C      and Stegun "Handbook of Mathematical Functions".
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0D0)
      DOUBLE PRECISION thp
      PARAMETER (thp=0.03D0)
      DOUBLE PRECISION half
      PARAMETER (half=0.5D0)
      DOUBLE PRECISION pt7
      PARAMETER (pt7=0.7D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
      DOUBLE PRECISION fiften
      PARAMETER (fiften=15.0D0)
      DOUBLE PRECISION forty
      PARAMETER (forty=40.0D0)
      DOUBLE PRECISION hundrd
      PARAMETER (hundrd=100.0D0)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=-1.0D-13)
      DOUBLE PRECISION loghaf
      PARAMETER (loghaf=-0.6931471805599453D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION w,dfd,dfn,lcum,ltail
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION eps,epsold,eps3,eps5,dfdold,dfnold,a,b,la,lb,bda,
     +                 lbda,l1padb,l1pbda,trm,thpab,logz,logx,logy,z,x,
     +                 y,lambda,frac,sum
      INTEGER ierr
      LOGICAL qagb,qtail
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,log,max,min
C     ..
C     .. Save statement ..
      SAVE dfdold,dfnold,a,b,la,lb,bda,lbda,l1padb,l1pbda,thpab,trm,
     +     qagb,eps,epsold,eps3,eps5
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dln1px,spmpar,dln1pe,dbetrm,lbfrac,lasym,series
      EXTERNAL dln1px,spmpar,dln1pe,dbetrm,lbfrac,lasym,series
C     ..
C     .. External Subroutines ..
      EXTERNAL bratio
C     ..
C     .. Statement Functions ..
C
C     The following statement functions are defined:
C       (1) CTC - Calculates LTAIL from LCUM or LCUM from LTAIL.
C
      DOUBLE PRECISION ctc

      ctc(x) = dln1px(-exp(x))
C     ..
C     .. Data statements ..
      DATA dfnold/-1.0D0/,dfdold/-1.0D0/
      DATA eps/-1.0D0/,epsold/-2.0D0/
C     ..
C     .. Executable Statements ..
C
C     EPS, EPS3 and EPS5 are machine dependent and only need to be
C     calculated once per execution of a program calling ltlf.
C
      IF (eps.NE.epsold) THEN
          eps = spmpar(1)
          epsold = eps
          eps3 = hundrd*eps
          eps5 = fiften*eps
      END IF
C
C     This subroutine is designed so that it can efficiently be
C     called multiple times with the same DFN and DFD but different
C     values of the other parameters.  The calculations
C     inside the next if statement are done the first time
C     the subroutine is called with a specific DFN and DFD
C     but not on subsequent calls.
C
      IF ((dfnold.NE.dfn) .OR. (dfdold.NE.dfd)) THEN
          dfnold = dfn
          dfdold = dfd
          a = half*dfn
          b = half*dfd
C
C     Compute values based on a and b only. The following values
C     are defined:
C       (1) LA = LOG(A)
C       (2) LB = LOG(B)
C       (3) BDA = B/A
C       (4) LBDA = LOG(B/A)
C       (5) L1PADB = LOG(1 + A/B)
C       (6) L1PBDA = LOG(1 + B/A)
C       (7) TRM = DBETRM(A,B)
C       (8) THPAB = THP*MIN(A,B)
C       (9) QAGB = A .GE. B
C
          la = log(a)
          lb = log(b)
          bda = b/a
          lbda = log(bda)
          trm = dbetrm(a,b)
          l1padb = dln1px(one/bda)
          l1pbda = dln1px(bda)
          thpab = thp*min(a,b)
          qagb = a .GE. b
      END IF
C
C     Calculate logx and logy.
C
      logz = lbda - w
      logx = -dln1pe(-logz)
      logy = -dln1pe(logz)
      x = exp(logx)
      y = exp(logy)
C
C     Try Extreme Value Calculations.
C
C     First try for ltail
C
      IF ((abs(a-one)*logy).GE.tiny) THEN
          IF (qagb) THEN
              ltail = -half*lb + (a+b-half)*l1pbda + b* (logy-w) - trm

          ELSE
              ltail = (-a+half)*la + (a-one)*lb + (a+b-half)*l1padb +
     +                b*logx - trm
          END IF

          IF (ltail.LE.loghaf) THEN
              lcum = ctc(ltail)
              RETURN

          END IF

      END IF
C
C     If this fails try lcum
C
      IF ((abs(b-one)*logx).GE.tiny) THEN
          IF (qagb) THEN
              lcum = (-b+half)*lb + (b-one)*la + (a+b-half)*l1pbda +
     +               a*logy - trm

          ELSE
              lcum = -half*la + (a+b-half)*l1padb + a* (logx+w) - trm
          END IF

          IF (lcum.LE.loghaf) THEN
              ltail = ctc(lcum)
              RETURN

          END IF

      END IF
C
C     If asymptotic don't work try using bratio to calculate cum and tai
C     We will accept the answer if both cum and tail are > 0
C
      CALL bratio(b,a,x,y,ltail,lcum,ierr)
      qtail = ltail .LE. lcum
      IF (min(ltail,lcum).GT.zero) THEN
          IF (qtail) THEN
              ltail = log(ltail)
              lcum = ctc(ltail)

          ELSE
              lcum = log(lcum)
              ltail = ctc(lcum)
          END IF

          RETURN

      END IF
C
C     Now asymptotics and bratio haven't worked. If qtail is true
C     we will calculate tail.
C
      IF (qtail) THEN
C
C     If a > 100 and b > 100 and x close to the mean we will use lasym
C     to calculate ltail.
C
          IF (qagb) THEN
              lambda = b - (a+b)*x

          ELSE
              lambda = (a+b)*y - a
          END IF
C
          IF ((a.GT.hundrd) .AND. (b.GT.hundrd) .AND.
     +        (lambda.GT.zero) .AND. (lambda.LE.thpab)) THEN
              ltail = lasym(b,a,lambda,eps3)
              lcum = ctc(ltail)
              RETURN

          END IF
C
C     Use LBFRAC
C
          IF ((a.GT.forty) .AND. (b.GT.forty) .AND.
     +        (lambda.GE.thpab) .AND. (a*x.GE.pt7)) THEN
              frac = lbfrac(b,a,x,y,lambda,eps5)
              IF (qagb) THEN
                  ltail = half*lb + (a+b-half)*l1pbda + (b+a)*logy -
     +                    b*w - trm + frac

              ELSE
                  ltail = half*la + (a+b-half)*l1padb + (b+a)*logx +
     +                    a*w - trm + frac
              END IF

              lcum = ctc(ltail)
              RETURN

          END IF
C
C     Since all other methods haven't work we must use a series
C     evaluation. Note, we are computing the incomplete beta with
C     arguement X and parameters B and A.
C
          z = exp(logz)
          sum = series(a,b,z)
          IF (qagb) THEN
              ltail = -half*lb + (a+b-half)*l1pbda + (b+a-one)*logy -
     +                b*w - trm + sum

          ELSE
              ltail = -half*la + (a+b-half)*l1padb + (b+a-one)*logx +
     +                (a-one)*w - trm + sum
          END IF

          lcum = ctc(ltail)
          RETURN
C
C     If qtail is false we will calculate tail.
C
      ELSE
C
C     LASYMP.
C
          IF (qagb) THEN
              lambda = (a+b)*x - b

          ELSE
              lambda = a - (a+b)*y
          END IF
C
          IF ((a.GT.hundrd) .AND. (b.GT.hundrd) .AND.
     +        (lambda.GT.zero) .AND. (lambda.LE.thpab)) THEN
              lcum = lasym(a,b,lambda,eps3)
              ltail = ctc(lcum)
              RETURN

          END IF
C
C     Use LBFRAC
C
          IF ((a.GT.forty) .AND. (b.GT.forty) .AND.
     +        (lambda.GE.thpab) .AND. (b*y.GE.0.7D0)) THEN
              frac = lbfrac(a,b,y,x,lambda,eps5)
              IF (qagb) THEN
                  lcum = half*lb + (a+b-half)*l1pbda + (b+a)*logy -
     +                   b*w - trm + frac

              ELSE
                  lcum = half*la + (a+b-half)*l1padb + (b+a)*logx +
     +                   a*w - trm + frac
              END IF

              ltail = ctc(lcum)
              RETURN

          END IF
C
C     Since all other methods haven't work we must use a series
C     evaluation.
C
          z = one/exp(logz)
          sum = series(b,a,z)
          IF (qagb) THEN
              lcum = -half*lb + (a+b-half)*l1pbda + (b+a-one)*logy -
     +               (b-one)*w - trm + sum

          ELSE
              lcum = -half*la + (a+b-half)*l1padb + (b+a-one)*logx +
     +               a*w - trm + sum
          END IF

          ltail = ctc(lcum)
          RETURN
C
      END IF

      END
      DOUBLE PRECISION FUNCTION psi(xx)
C---------------------------------------------------------------------
C
C                 EVALUATION OF THE DIGAMMA FUNCTION
C
C                           -----------
C
C     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
C     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
C     CODY, STRECOK AND THACHER.
C
C---------------------------------------------------------------------
C     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
C     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
C     A.H. MORRIS (NSWC).
C---------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION xx
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aug,den,dx0,piov4,sgn,upper,w,x,xmax1,xmx0,
     +                 xsmall,z
      INTEGER i,m,n,nq
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION p1(7),p2(4),q1(6),q2(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION spmpar
      INTEGER ipmpar
      EXTERNAL spmpar,ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,cos,dble,dlog,dmin1,int,sin
C     ..
C     .. Data statements ..
C---------------------------------------------------------------------
C
C     PIOV4 = PI/4
C     DX0 = ZERO OF PSI TO EXTENDED PRECISION
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
C
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
C
C---------------------------------------------------------------------
      DATA piov4/.785398163397448D0/
      DATA dx0/1.461632144968362341262659542325721325D0/
      DATA p1(1)/.895385022981970D-02/,p1(2)/.477762828042627D+01/,
     +     p1(3)/.142441585084029D+03/,p1(4)/.118645200713425D+04/,
     +     p1(5)/.363351846806499D+04/,p1(6)/.413810161269013D+04/,
     +     p1(7)/.130560269827897D+04/
      DATA q1(1)/.448452573429826D+02/,q1(2)/.520752771467162D+03/,
     +     q1(3)/.221000799247830D+04/,q1(4)/.364127349079381D+04/,
     +     q1(5)/.190831076596300D+04/,q1(6)/.691091682714533D-05/
      DATA p2(1)/-.212940445131011D+01/,p2(2)/-.701677227766759D+01/,
     +     p2(3)/-.448616543918019D+01/,p2(4)/-.648157123766197D+00/
      DATA q2(1)/.322703493791143D+02/,q2(2)/.892920700481861D+02/,
     +     q2(3)/.546117738103215D+02/,q2(4)/.777788548522962D+01/
C     ..
C     .. Executable Statements ..
C---------------------------------------------------------------------
C
C     MACHINE DEPENDENT CONSTANTS ...
C
C        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
C                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
C                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
C                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
C                 PSI MAY BE REPRESENTED AS ALOG(X).
C
C        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
C                 MAY BE REPRESENTED BY 1/X.
C
C---------------------------------------------------------------------
      xmax1 = ipmpar(3)
      xmax1 = dmin1(xmax1,1.0D0/spmpar(1))
      xsmall = 1.D-9
C---------------------------------------------------------------------
      x = xx
      aug = 0.0D0
      IF (x.GE.0.5D0) GO TO 50
C---------------------------------------------------------------------
C     X .LT. 0.5,  USE REFLECTION FORMULA
C     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
C---------------------------------------------------------------------
      IF (abs(x).GT.xsmall) GO TO 10
      IF (x.EQ.0.0D0) GO TO 100
C---------------------------------------------------------------------
C     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
C     FOR  PI*COTAN(PI*X)
C---------------------------------------------------------------------
      aug = -1.0D0/x
      GO TO 40
C---------------------------------------------------------------------
C     REDUCTION OF ARGUMENT FOR COTAN
C---------------------------------------------------------------------
   10 w = -x
      sgn = piov4
      IF (w.GT.0.0D0) GO TO 20
      w = -w
      sgn = -sgn
C---------------------------------------------------------------------
C     MAKE AN ERROR EXIT IF X .LE. -XMAX1
C---------------------------------------------------------------------
   20 IF (w.GE.xmax1) GO TO 100
      nq = int(w)
      w = w - dble(nq)
      nq = int(w*4.0D0)
      w = 4.0D0* (w-dble(nq)*.25D0)
C---------------------------------------------------------------------
C     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
C     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
C     QUADRANT AND DETERMINE SIGN
C---------------------------------------------------------------------
      n = nq/2
      IF ((n+n).NE.nq) w = 1.0D0 - w
      z = piov4*w
      m = n/2
      IF ((m+m).NE.n) sgn = -sgn
C---------------------------------------------------------------------
C     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
C---------------------------------------------------------------------
      n = (nq+1)/2
      m = n/2
      m = m + m
      IF (m.NE.n) GO TO 30
C---------------------------------------------------------------------
C     CHECK FOR SINGULARITY
C---------------------------------------------------------------------
      IF (z.EQ.0.0D0) GO TO 100
C---------------------------------------------------------------------
C     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
C     SIN/COS AS A SUBSTITUTE FOR TAN
C---------------------------------------------------------------------
      aug = sgn* ((cos(z)/sin(z))*4.0D0)
      GO TO 40

   30 aug = sgn* ((sin(z)/cos(z))*4.0D0)
   40 x = 1.0D0 - x
   50 IF (x.GT.3.0D0) GO TO 70
C---------------------------------------------------------------------
C     0.5 .LE. X .LE. 3.0
C---------------------------------------------------------------------
      den = x
      upper = p1(1)*x
C
      DO 60 i = 1,5
          den = (den+q1(i))*x
          upper = (upper+p1(i+1))*x
   60 CONTINUE
C
      den = (upper+p1(7))/ (den+q1(6))
      xmx0 = dble(x) - dx0
      psi = den*xmx0 + aug
      RETURN
C---------------------------------------------------------------------
C     IF X .GE. XMAX1, PSI = LN(X)
C---------------------------------------------------------------------
   70 IF (x.GE.xmax1) GO TO 90
C---------------------------------------------------------------------
C     3.0 .LT. X .LT. XMAX1
C---------------------------------------------------------------------
      w = 1.0D0/ (x*x)
      den = w
      upper = p2(1)*w
C
      DO 80 i = 1,3
          den = (den+q2(i))*w
          upper = (upper+p2(i+1))*w
   80 CONTINUE
C
      aug = upper/ (den+q2(4)) - 0.5D0/x + aug
   90 psi = aug + dlog(x)
      RETURN
C---------------------------------------------------------------------
C     ERROR RETURN
C---------------------------------------------------------------------
  100 psi = 0.0D0
      RETURN

      END
      DOUBLE PRECISION FUNCTION rexp(x)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION EXP(X) - 1
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,q1,q2,q3,q4,w
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA p1/.914041914819518D-09/,p2/.238082361044469D-01/,
     +     q1/-.499999999085958D+00/,q2/.107141568980644D+00/,
     +     q3/-.119041179760821D-01/,q4/.595130811860248D-03/
C     ..
C     .. Executable Statements ..
C-----------------------
      IF (abs(x).GT.0.15D0) GO TO 10
      rexp = x* (((p2*x+p1)*x+1.0D0)/ ((((q4*x+q3)*x+q2)*x+q1)*x+1.0D0))
      RETURN
C
   10 w = exp(x)
      IF (x.GT.0.0D0) GO TO 20
      rexp = (w-0.5D0) - 0.5D0
      RETURN

   20 rexp = w* (0.5D0+ (0.5D0-1.0D0/w))
      RETURN

      END
      DOUBLE PRECISION FUNCTION rlog1(x)
C-----------------------------------------------------------------------
C             EVALUATION OF THE FUNCTION X - LN(1 + X)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION a,b,h,p0,p1,p2,q1,q2,r,t,w,w1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,dlog
C     ..
C     .. Data statements ..
C------------------------
      DATA a/.566749439387324D-01/
      DATA b/.456512608815524D-01/
      DATA p0/.333333333333333D+00/,p1/-.224696413112536D+00/,
     +     p2/.620886815375787D-02/
      DATA q1/-.127408923933623D+01/,q2/.354508718369557D+00/
C     ..
C     .. Executable Statements ..
C------------------------
      IF (x.LT.-0.39D0 .OR. x.GT.0.57D0) GO TO 40
      IF (x.LT.-0.18D0) GO TO 10
      IF (x.GT.0.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
      h = x
      w1 = 0.0D0
      GO TO 30
C
   10 h = dble(x) + 0.3D0
      h = h/0.7D0
      w1 = a - h*0.3D0
      GO TO 30
C
   20 h = 0.75D0*dble(x) - 0.25D0
      w1 = b + h/3.0D0
C
C               SERIES EXPANSION
C
   30 r = h/ (h+2.0D0)
      t = r*r
      w = ((p2*t+p1)*t+p0)/ ((q2*t+q1)*t+1.0D0)
      rlog1 = 2.0D0*t* (1.0D0/ (1.0D0-r)-r*w) + w1
      RETURN
C
C
   40 w = (x+0.5D0) + 0.5D0
      rlog1 = x - dlog(w)
      RETURN

      END
      DOUBLE PRECISION FUNCTION series(a,b,z)
C----------------------------------------------------------------------
C
C    Calculates log of
C       1 + [Beta(A+1,n+1)/Beta(B-n-1,n+1)]*(x/(1-x))^(n+1)] n=0,1,...
C    This summation is part of Abramowitz and Stegun series 26.5.5
C       [(x^b)*((1-x)^(a-1))/(b*Beta(a,b))]*SUM
C    which can be used to calculate IX(A,B).
C
C----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0D0)
      DOUBLE PRECISION one
      PARAMETER (one=1.0D0)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=1.0D-13)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION a,b,z
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION am2,bp2,term,sum,xn,ratio
      LOGICAL qdone
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. External Functions ..
      DOUBLE PRECISION dln1px
      EXTERNAL dln1px
C     ..
C     .. Executable Statements ..
C
C     initialize variables
C
C
      qdone = .FALSE.

      am2 = a - one - one
      bp2 = b + one + one
      term = (a-one)*z/ (b+one)
      sum = zero
      xn = zero
C
C     Compute sum
C
      GO TO 20

   10 IF (qdone .OR. (xn.GT.100D0)) GO TO 30
   20 sum = sum + term
      qdone = abs(term/sum) .LE. tiny
      IF (.NOT.qdone) THEN
          ratio = (am2-xn)*z/ (bp2+xn)
          term = term*ratio
          xn = xn + one
      END IF

      GO TO 10

   30 series = dln1px(sum)


      RETURN

      END
      DOUBLE PRECISION FUNCTION spmpar(i)
C-----------------------------------------------------------------------
C
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C-----------------------------------------------------------------------
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C-----------------------------------------------------------------------
C     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
C     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
C     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER i
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION b,binv,bm1,one,w,z
      INTEGER emax,emin,ibeta,m
C     ..
C     .. External Functions ..
      INTEGER ipmpar
      EXTERNAL ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble
C     ..
C     .. Executable Statements ..
C
      IF (i.GT.1) GO TO 10
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
      RETURN
C
   10 IF (i.GT.2) GO TO 20
      b = ipmpar(4)
      emin = ipmpar(9)
      one = dble(1)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
      RETURN
C
   20 ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
C
      b = ibeta
      bm1 = ibeta - 1
      one = dble(1)
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
C
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b
      RETURN

      END
