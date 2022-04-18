  MODULE GammaError
  USE set_precision
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: errorfunction, inverfc, gamma, &
             loggam, quotgamm, lnec, hypfun, auxgam
  CONTAINS 
    RECURSIVE FUNCTION errorfunction (x, erfcc, expo) RESULT(errfu)
    !-----------------------------------------------------------------
    ! Computation of: the error function, the complementary error 
    !                 function or the scaled complementary error
    !                 function. 
    !-----------------------------------------------------------------
    ! Inputs:
    !       x, real number.
    !       erfcc, logical variable.  
    !       expo,  logical variable.
    !
    ! The meaning of the logical variables erfcc and expo 
    !   is the following:
    !    When erfcc=.true. and expo=.false., the function computes 
    !         the complementary error function erfc(x). 
    !    When erfcc=.true., expo=.true. and x>0, the function computes
    !         the scaled complementary error function exp(x*x}erfc(x). 
    !    When erfcc=.false. and expo=.false., the function computes 
    !         the error function erf(x).
    !------------------------------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, y, z, r(0:8), s(0:8), errfu
    LOGICAL erfcc, expo
    IF (erfcc) THEN
      IF (x < -6.5_r8) THEN
        y= 2.0_r8 
      ELSEIF (x < 0.0_r8) THEN
        y= 2.0_r8 - errorfunction(-x, .true., .false.) 
      ELSEIF (x == 0.0_r8) THEN
        y= 1.0_r8 
      ELSEIF (x < 0.5_r8) THEN
        IF (expo) THEN
          y=exp(x*x)
        ELSE
          y=1.0_r8
        ENDIF
        y=y*(1.0_r8-errorfunction(x, .false., .false.))
      ELSEIF (x < 4.0_r8) THEN
        IF (expo) THEN
          y= 1.0_r8 
        ELSE
          y= exp(-x*x)
        ENDIF
        r(0)= 1.230339354797997253e3_r8
        r(1)= 2.051078377826071465e3_r8
        r(2)= 1.712047612634070583e3_r8
        r(3)= 8.819522212417690904e2_r8
        r(4)= 2.986351381974001311e2_r8
        r(5)= 6.611919063714162948e1_r8
        r(6)= 8.883149794388375941_r8
        r(7)= 5.641884969886700892e-1_r8
        r(8)= 2.153115354744038463e-8_r8
        s(0)= 1.230339354803749420e3_r8
        s(1)= 3.439367674143721637e3_r8
        s(2)= 4.362619090143247158e3_r8
        s(3)= 3.290799235733459627e3_r8
        s(4)= 1.621389574566690189e3_r8
        s(5)= 5.371811018620098575e2_r8
        s(6)= 1.176939508913124993e2_r8
        s(7)= 1.574492611070983473e1_r8
        y=y*fractio(x,8,r,s)
      ELSE
        z=x*x
        IF (expo) THEN
          y=1.0_r8 
        ELSE
          y= exp(-z)
        ENDIF
        z=1.0_r8/z
        r(0)=6.587491615298378032e-4_r8
        r(1)=1.608378514874227663e-2_r8
        r(2)=1.257817261112292462e-1_r8
        r(3)=3.603448999498044394e-1_r8
        r(4)=3.053266349612323440e-1_r8
        r(5)=1.631538713730209785e-2_r8
        s(0)=2.335204976268691854e-3_r8
        s(1)=6.051834131244131912e-2_r8
        s(2)=5.279051029514284122e-1_r8
        s(3)=1.872952849923460472_r8
        s(4)=2.568520192289822421_r8
        y=y*(oneoversqrtpi-z*fractio(z,5,r,s))/x       
      ENDIF
      errfu=y
    ELSE
      IF (x == 0.0_r8) THEN 
        y=0.0_r8
      ELSEIF (abs(x) > 6.5_r8) THEN 
        y=x/abs(x)
      ELSEIF (x > 0.5_r8) THEN
        y=1.0_r8 - errorfunction(x, .true., .false.) 
      ELSEIF (x < -0.5_r8) THEN
        y=errorfunction(-x, .true., .false.)-1.0_r8
      ELSE
        r(0)=3.209377589138469473e3_r8
        r(1)=3.774852376853020208e2_r8
        r(2)=1.138641541510501556e2_r8
        r(3)=3.161123743870565597e0_r8
        r(4)=1.857777061846031527e-1_r8
        s(0)=2.844236833439170622e3_r8
        s(1)=1.282616526077372276e3_r8
        s(2)=2.440246379344441733e2_r8
        s(3)=2.360129095234412093e1_r8
        z=x*x
        y=x*fractio(z,4,r,s)
      ENDIF  
      errfu= y
    ENDIF        
    END FUNCTION errorfunction

    RECURSIVE FUNCTION inverfc(x) RESULT(y)
    !---------------------------------------------
    ! Inverse of the complementary error function
    !---------------------------------------------
    ! Input:
    !       x, real positive number.
    !---------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) ::  x, y, y0, y02, h, r, f, fp, c1, c2, c3, c4, c5;
    IF (x==0.0_r8) THEN
      y=giant
    ELSE
      IF (x > 1) THEN
        y=-inverfc(2.0_r8-x)
      ELSE
        y0=0.70710678_r8*invq(x*0.5_r8);
        f=errorfunction(y0,.true.,.false.)-x;
        y02= y0*y0;
        fp=-2.0_r8/sqrtpi*exp(-y02);
        c1=-1.0_r8/fp;
        c2= y0;
        c3=(4.0_r8*y02+1.0_r8)*0.333333333333333333333333333333_r8;
        c4=y0*(12*y02+7.0_r8)*.166666666666666666666666666667_r8;
        c5=(8*y02+7)*(12*y02+1.0_r8)*0.333333333333333333333333333333_r8;
        r= f*c1;
        h=r*(1+r*(c2+r*(c3+r*(c4+r*c5))));
        y=y0+h
      ENDIF
    ENDIF
    END FUNCTION inverfc
   
    RECURSIVE FUNCTION gamma(x) RESULT(gam)
    !----------------------------------------
    ! Computation of the Euler gamma function 
    ! Gamma(x).
    !----------------------------------------
    !  Input:
    !       x, real number.
    !----------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, dw, gam, z
    INTEGER :: k, k1, n
    k=nint(x)
    k1=k-1
    IF (k==0) THEN
      dw=dwarf
    ELSE
      dw=machtol
    ENDIF
    IF ((k <= 0).AND.(abs(k - x)<= dw)) THEN
      IF (mod(k,2)>0) THEN
        ! k is odd
        gam=sign(1.0_r8,k-x)*giant
      ELSE
        ! k is even
        gam=sign(1.0_r8,x-k)*giant
      ENDIF
    ELSEIF (x<0.45_r8) THEN
      gam=pi/(sin(pi*x)*gamma(1.0_r8-x))
    ELSEIF ((abs(k-x)<dw).AND.(x<21.0_r8)) THEN
      gam=1;
      DO n=2,k1 
        gam=gam*n
      ENDDO
    ELSEIF ((abs(k-x-0.5_r8)<dw).AND.(x<21.0_r8)) THEN
      gam=sqrt(pi);
      DO n=1,k1 
        gam=gam*(n-0.5)
      ENDDO
    ELSEIF (x<3.0_r8) THEN
      IF (k>x) THEN
        k=k1
      ENDIF
      k1=3-k;
      z=k1+x;
      gam=gamma(z);
      DO n=1,k1 
        gam=gam/(z-n)
      ENDDO
    ELSE
      gam=sqrttwopi*exp(-x+(x-0.5_r8)*log(x)+stirling(x))
    ENDIF
    END FUNCTION gamma

    FUNCTION loggam(x)
    !------------------------------------
    ! Computation of the logarithm of
    ! the Gamma Function, log(gamma(x))
    !------------------------------------
    ! Input:
    !       x, real positive number.
    !------------------------------------
    USE Someconstants
    IMPLICIT NONE  
    REAL(r8) :: loggam, x
    IF (x>=3) THEN
      loggam=(x-0.5_r8)*log(x)- x+lnsqrttwopi+stirling(x)
    ELSEIF (x >= 2) THEN
      loggam=(x-2.0_r8)*(3.0_r8-x)*auxloggam(x-2.0_r8)+logoneplusx(x-2.0_r8)
    ELSEIF (x>=1) THEN
      loggam=(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8)
    ELSEIF (x>0.5_r8) THEN
      loggam=x*(1.0_r8-x)*auxloggam(x)-logoneplusx(x-1.0_r8)
    ELSEIF (x>0.0_r8) THEN
      loggam=x*(1.0_r8-x)*auxloggam(x)-log(x)
    ELSE
      loggam=giant
    ENDIF
    END FUNCTION loggam
  
    FUNCTION gamstar(x)
    !-----------------------------------------------
    ! Computation of the regulated gamma function
    !-----------------------------------------------
    ! Input: x, real positive number.
    !-----------------------------------------------
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: gamstar, x
    IF (x>=3.0_r8) THEN
      gamstar=exp(stirling(x))
    ELSEIF (x>0.0_r8) THEN
      gamstar=gamma(x)/(exp(-x+(x-0.5_r8)*log(x))*sqrttwopi)
    ELSE
      gamstar=giant
    ENDIF
    END FUNCTION gamstar

    RECURSIVE FUNCTION quotgamm(x,y) RESULT(quotgam)
    !----------------------------------------------------
    ! Computation of the quotient of two gamma functions
    !  gamma(x)/gamma(y)
    !----------------------------------------------------
    USE Someconstants
    IMPLICIT NONE  
    REAL(r8) :: quotgam,x,y
    INTEGER :: n
    IF ((x <= 0.0_r8).OR.(y <= 0.0_r8)) THEN
      quotgam=qratio(x,y)
    ELSEIF (x>y) THEN
      quotgam=1.0_r8/quotgamm(y,x)
    ELSE
      n=int(y-x);
      IF (n==(y - x)) THEN
        quotgam=1.0_r8/shiftfact(x, n)
      ELSEIF (n>15) THEN
        quotgam=qratio(x,y)
      ELSEIF (n>0) THEN
        quotgam=quotgamm(x+n,y)/shiftfact(x,n)
      ELSEIF (x<26.0_r8) THEN
        quotgam=qratio(x,y)
      ELSE
        quotgam=qasym(x,y)
      ENDIF
    ENDIF
    END FUNCTION quotgamm

    FUNCTION qratio(x,y)
    USE Someconstants
    IMPLICIT NONE                                        
    REAL(r8) :: qratio,x,y
    REAL(r8) :: b, c, q
    IF (((x<=0.0_r8).OR.(y<=0.0_r8)).OR.(y>1.5_r8*x)) THEN
      qratio=gamma(x)/gamma(y)
    ELSE
      b=y-x;
      c=b/x;
      q=lnec(c);
      q=exp(-b*log(x)-(b-0.5_r8)*log(1.0_r8+c)-x*q);
      qratio=q*gamstar(x)/gamstar(y)
    ENDIF
    END FUNCTION qratio

    FUNCTION qasym(x,y)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: qasym, x, y   
    REAL(r8) :: s, t, u, w, w2, r, r2, r3, r4, r5, r6, r7, cc(0:7) 
    INTEGER :: k                          
    w=(x+y-1)/2.0_r8;
    w2=1.0_r8/(w*w);
    r=(x-y+1)/2.0_r8;
    r2=r*r;
    r3=r2*r;
    r4=r2*r2;
    r5=r4*r;
    r6=r5*r;
    r7=r6*r;
    cc(0)=1.0_r8;
    cc(1)=r/12.0_r8;
    cc(2)=r/1440.0_r8+r2/288.0_r8;
    cc(3)=r/90720.0_r8+r2/17280.0_r8+r3/10368.0_r8;
    cc(4)=r/4838400.0_r8+101.0_r8*r2/87091200.0_r8+r3/414720.0_r8+r4/497664.0_r8;
    cc(5)=r/239500800.0_r8+13*r2/522547200.0_r8+61*r3/1045094400.0_r8+&
          r4/14929920.0_r8+r5/29859840.0_r8;
    cc(6)=691.0_r8*r/7846046208000.0_r8+7999.0_r8*r2/14485008384000.0_r8+&
          59.0_r8*r3/41803776000.0_r8+143.0_r8*r4/75246796800.0_r8&
          +r5/716636160.0_r8+r6/2149908480.0_r8;
    cc(7)=r/523069747200.0_r8 + 2357.0_r8*r2/188305108992000.0_r8+&
          5941.0_r8*r3/173820100608000.0_r8+11.0_r8*r4/214990848000.0_r8+&
          41.0_r8*r5/902961561600.0_r8+r6/42998169600.0_r8+r7/180592312320.0_r8;
    s=1.0_r8;
    k=1;
    t=1.0_r8;
    u=1.0_r8
    DO WHILE ((abs(u)>machtol).AND.(k<8))                                   
      t=-4.0_r8*w2*(k-r)*(k-r-0.5)*t;
      u=t*cc(k);
      s=s+u;
      k=k+1;
      qasym=s*exp((x-y)*log(w))
    ENDDO
    END FUNCTION qasym

    RECURSIVE FUNCTION shiftfact(x,n) RESULT(shift)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: shift, x, s
    INTEGER :: k,n 
    IF (n==0) THEN
      shift=1.0_r8
    ELSEIF (n<0) THEN
      shift=1.0_r8/shiftfact(x-n,n)
    ELSE
      s=1.0_r8;
      k=0;
      DO WHILE (k<n)
        s=s*(x+k);
        k=k+1
      ENDDO                                                                                                     
      shift=s
    ENDIF
    END FUNCTION shiftfact

    FUNCTION invq(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: invq, x, t
    !Abramowitx & Stegun 26.2.23; 
    t=sqrt(-2*log(x));
    t=t-(2.515517_r8+t*(0.802853_r8+t*0.010328_r8))/&
      (1.0_r8+t*(1.432788+t*(0.189269_r8+t*0.001308_r8)))
    invq=t
    END FUNCTION invq

    FUNCTION factor(x,n)
    IMPLICIT NONE
    REAL(r8) :: x, facto, factor
    INTEGER :: n, i
    facto=1
    DO i=1,n
      facto=facto*(x/i)
    ENDDO
    factor=facto
    END FUNCTION factor

    FUNCTION pol(fjkm,d,v)
    IMPLICIT NONE
    REAL(r8) :: pol,v,s,fjkm(0:32)
    INTEGER :: d,m
    m=d; s= fjkm(d);
    DO WHILE (m > 0)
      m=m-1; s=s*v + fjkm(m)
    ENDDO
    pol=s
    END FUNCTION pol

    FUNCTION ps(mu,mulnmu,lnx,y,lny,n, a, b)
    USE Someconstants
    REAL(r8) :: mu,mulnmu,lneps,lnx,y,lny,ps, n, f
    INTEGER :: a, b 
    lneps=log(epss) 
    IF ((a==0).AND.(b==0)) THEN
      f=(n-lneps)/(log(n)-lnx) 
    ELSEIF ((a==0).AND.(b==1)) THEN
      f=(2*n-lneps+mulnmu-mu*log(mu+n))/(log(n)-lnx-lny+log(mu+n)) 
    ELSEIF ((a==1).AND.(b==0)) THEN
      f=(2*n-lneps-y+mu*lny-mu*log(mu+n)+mu)/(log(n)-lnx-lny+log(mu+n));
    ENDIF
    ps= f
    END FUNCTION ps
   
    SUBROUTINE hypfun(x,sinh,cosh);
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,sinh,cosh,ss
    REAL(r8) :: ax, y, f, f2
    ax=abs(x);
    IF (ax<0.21_r8) THEN
      IF (ax<0.07_r8) THEN
        y=x*x
      ELSE
        y=x*x/9.0_r8;                                                          
      ENDIF
      f=2.0_r8+y*(y*28+2520.0_r8)/(y*(y+420)+15120.0_r8);
      f2=f*f;
      sinh=2*x*f/(f2-y);
      cosh=(f2+ y)/(f2-y);
      IF (ax>=0.07_r8) THEN
        ss=2.0_r8*sinh/3.0_r8
        f=ss*ss                                                      
        sinh=sinh*(1.0_r8+f/3.0_r8);                                        
        cosh=cosh*(1.0_r8+f)                                           
      ENDIF
    ELSE
      y=exp(x);
      f=1.0_r8/y;
      cosh=(y+f)/2.0_r8;
      sinh=(y-f)/2.0_r8
    END IF
    END SUBROUTINE hypfun        

    FUNCTION startkbes(x,eps)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,eps
    REAL(r8) :: y, p, q, r, s, c, del, r2
    INTEGER :: startkbes
    IF (eps<machtol) THEN
      del=-log(machtol/2.0_r8)
    ELSE
      del=-log(eps/2.0_r8);
    ENDIF
    p=x+del;
    q=p/x;
    IF (q<2.0_r8) THEN
      r=log((q+1.0_r8)/(q-1.0_r8))/2.0_r8;
      y=r*(1.0_r8+2.0_r8/(1.0_r8+r*(q+1.0_r8/q)))
    ELSE
      r=2.0_r8*x/p;
      r2=r*r
      y=r*(1.0_r8+r2*r2/45.0_r8)
    ENDIF
    CALL hypfun(y,s,c);
    startkbes=1+int(x/(2.0_r8*s*s))
    END FUNCTION  startkbes

    FUNCTION startijbes(x,n,t,eps)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x,eps
    REAL(r8) :: p, q, r, y, del
    INTEGER :: startijbes, n,t, s
    IF (x<=0.0_r8) THEN
      startijbes=0
    ELSE
      s=2*t-1;
    ENDIF
    IF (eps<machtol) THEN
      del=-log(machtol/2.0_r8)
    ELSE
      del=-log(eps/2.0_r8);
    ENDIF
    p=del/x-t;
    r=n/x;
    IF ((r>1.0_r8).OR.(t==1)) THEN
      q=sqrt(r*r+s);
      r=r*log(q+r)-q
    ELSE
      r=0;
    ENDIF
    q=del/(2.0_r8*x)+r;
    IF (p>q) THEN
      r=p
    ELSE
      r=q;
    ENDIF
    y=alfinv(t,r,p,q)
    CALL hypfun(y,p,q)
    IF (t==0) THEN
      s=int(x*q)+1
    ELSE
      s=int(x*p)+1;
    ENDIF
    IF (MOD(s,2)>0) THEN
      s=s+1;
    ENDIF
    startijbes=s
    END FUNCTION startijbes
        
    FUNCTION alfinv(t,r,p,q)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: alfinv, r, p, q, a, b, a2, lna
    INTEGER :: t
    p=0.0_r8
    IF ((t+r)<2.7_r8) THEN
      IF (t==0) THEN
        a=exp(log(3.0_r8*r)/3.0_r8);
        a2=a*a;                                 
        b=a*(1.0_r8+a2*(-1.0_r8/30.0_r8+0.004312_r8*a2))                       
      ELSE
        a=sqrt(2.0_r8*(1.0_r8+r));
        a2=a*a;
        b=a/(1.0_r8+a2/8.0_r8)                     
      ENDIF
    ELSE
      a=log(0.7357589_r8*(r+t));
      lna=log(a)/a
      b=1.0_r8+a+log(a)*(1.0_r8/a-1.0_r8)+0.5_r8*lna*lna
    ENDIF
    DO WHILE (abs(a/b-1.0_r8)>1.0e-2_r8)                                        
      a=b;
      b=fi(a,r,t,q)
    ENDDO
    alfinv=b
    END FUNCTION alfinv

    FUNCTION falfa(al,r,t,df)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: falfa, al, r, df, sh, ch
    INTEGER :: t
    CALL hypfun(al,sh,ch)
    IF (t==1) THEN
      falfa=al*sh/ch-1.0_r8-r/ch;
      df=(sh+(al+r*sh)/ch)/ch
    ELSE
      falfa=al-(sh+r)/ch;
      df=sh*(r+sh)/(ch*ch)
    ENDIF
    END FUNCTION falfa

    FUNCTION fi(al,r,t,q)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: fi, al, r, q, p
    INTEGER :: t
    p=falfa(al,r,t,q);
    fi=al-p/q
    END FUNCTION fi

    FUNCTION recipgam(x,q,r)
    USE Someconstants
    IMPLICIT NONE                                                                                                       
    !recipgam(x)=1/gamma(x+1)=1 + x * (q + x * r); -0.5<=x<=0.5}
    REAL(r8) :: recipgam,x,q,r
    REAL(r8) :: t, tx, c(0:8)
    IF (x==0) THEN
      q=0.5772156649015328606e-0_r8;
      r=-0.6558780715202538811e-0_r8
    ELSE
      c(0)=+1.142022680371167841_r8;
      c(1)=-6.5165112670736881e-3_r8;
      c(2)=-3.087090173085368e-4_r8;
      c(3)=+3.4706269649043e-6_r8;
      c(4)=-6.9437664487e-9_r8;
      c(5)=-3.67795399e-11_r8;
      c(6)=+1.356395e-13_r8;
      c(7)=+3.68e-17_r8;
      c(8)=-5.5e-19_r8;
      tx=2.0_r8*x
      t=2*tx*tx-1;
      q=chepolsum(8,t,c);
      c(0)=-1.270583625778727532_r8;
      c(1)=+2.05083241859700357e-2_r8;
      c(2)=-7.84761097993185e-5_r8;
      c(3)=-5.377798984020e-7_r8;
      c(4)=+3.8823289907e-9_r8;
      c(5)=-2.6758703e-12_r8;
      c(6)=-2.39860e-14_r8;
      c(7)=+3.80e-17_r8;
      c(8)=+4e-20_r8;
      r=chepolsum(8,t,c)
    ENDIF               
    recipgam=1+x*(q+x*r)
    END FUNCTION recipgam   
    
    FUNCTION xpowy(x,y)
    IMPLICIT NONE
    REAL(r8) :: x,y,xpowy
    xpowy=x**y
    END FUNCTION xpowy

    FUNCTION fractio(x,n,r,s)
    IMPLICIT NONE
    INTEGER n,k
    REAL(r8) :: x, fractio, r(0:n), s(0:n), a, b
    a=r(n); b=1
    DO k=n-1,0,-1 
      a=a*x+r(k); b=b*x+s(k) 
    ENDDO
    fractio=a/b
    END FUNCTION fractio
   
    FUNCTION chepolsum(n,t,ak)
    USE Someconstants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(r8), INTENT(IN) :: t, ak(0:n)
    REAL(r8) :: chepolsum, u0, u1, u2, s, tt;
    INTEGER :: k
    u0=0; u1=0; k=n; tt=t+t;
    DO WHILE (k>=0)
      u2=u1; 
      u1=u0; 
      u0=tt*u1-u2+ ak(k); 
      k= k-1 
    ENDDO
    s=(u0-u2)/2.0_r8
    chepolsum=s
    END FUNCTION chepolsum

    FUNCTION oddchepolsum(n, x, ak)
    USE Someconstants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(r8), INTENT(IN) :: x, ak(0:n)
    REAL(r8) :: oddchepolsum, h, r, s, y;
    INTEGER :: k
    IF (n==0) THEN 
      s=ak(0)*x 
    ELSEIF  (n == 1) THEN
      s=x*(ak(0)+ak(1)*(4*x*x - 3)) 
    ELSE
      y=2*(2*x*x - 1); 
      r= ak(n); h= ak(n-1)+ r*y;  
      k=n-2;
      DO WHILE (k >= 0)
        s=r; r= h; h= ak(k)+r*y-s; 
        k= k-1 
      ENDDO
      s=x*(h-r)
    ENDIF
    oddchepolsum=s
    END FUNCTION oddchepolsum

    FUNCTION logoneplusx(t)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: t
    REAL(r8) ::  logoneplusx, ck(0:100), c, p, p2, pj, x, y
    INTEGER :: j
    IF ((-0.2928_r8<t).AND.(t<0.4142_r8)) THEN
      p= twoexp14
      p=(p-1.0_r8)/(p+1.0_r8);
      pj=p; ck(0)= pj; p2= p*p; j=1; c=1.0_r8;
      DO WHILE ((abs(c)> 1.0e-20_r8).AND.(j<1000))   
        pj=pj*p2; c=pj/(2.0_r8*j+1.0_r8); 
        ck(j)=c; j= j+1 
      ENDDO
      x=t/(2.0_r8+t)*(1.0_r8+p2)/(2.0_r8*p);
      y=4*oddchepolsum(j-1, x, ck)
    ELSE 
      y=log(1.0_r8+t) 
    ENDIF
    logoneplusx=y
    END FUNCTION logoneplusx

    FUNCTION xminsinx(x) 
    !  {(x-sin(x))/(x^3/6)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: xminsinx, f, fk(0:8), t
    IF (abs(x)>1.0_r8) THEN
      f=6*(x-sin(x))/(x*x*x) 
    ELSE
      fk(0)=1.95088260487819821294e-0_r8;
      fk(1)=-0.244124470324439564863e-1_r8;
      fk(2)=0.14574198156365500e-3_r8;
      fk(3)=-0.5073893903402518e-6_r8;
      fk(4)=0.11556455068443e-8_r8;
      fk(5)=-0.185522118416e-11_r8;
      fk(6)=0.22117315e-14_r8;
      fk(7)=-0.2035e-17_r8;
      fk(8)=0.15e-20_r8;
      t=2*x*x-1.0_r8;
      f=chepolsum(8,t,fk)
    ENDIF
    xminsinx=f
    END FUNCTION xminsinx
        
  
    RECURSIVE FUNCTION sinh(x,eps) RESULT(sinhh)
    USE Someconstants  
    IMPLICIT NONE
    !to compute hyperbolic function sinh (x)}
    REAL(r8) :: sinhh, x, eps
    REAL(r8) :: ax, e, t, x2, y
    INTEGER  :: u, k
    ax=abs(x);
    IF (x==0.0_r8) THEN
      y=0.0_r8
    ELSEIF (ax<0.12) THEN
      e=eps*0.1_r8;
      x2=x*x;
      y=1;
      t=1;
      u=0;
      k=1;
      DO WHILE(t>e)
        u=u+8*k-2;
        k=k+1;
        t=t*x2/u;
        y=y+t
      END DO
      y=x*y
    ELSEIF (ax<0.36_r8) THEN
      t=sinh(x*0.333333333333333333333333333333_r8,eps);
      y=t*(3.0_r8+4.0_r8*t*t);
    ELSE
      t=exp(x);
      y=(t-1.0_r8/t)*0.5_r8
    ENDIF
    sinhh=y
    END FUNCTION sinh

    FUNCTION exmin1(x,eps)
    USE Someconstants  
    IMPLICIT NONE
    !computes (exp(x)-1)/x 
    REAL(r8) :: exmin1, x, eps
    REAL(r8) :: t, y
    IF (x==0) THEN
      y=1.0_r8
    ELSEIF ((x<-0.69_r8).OR.(x>0.4_r8)) THEN
      y=(exp(x)-1.0_r8)/x
    ELSE
      t=x*0.5_r8;
      y=exp(t)*sinh(t,eps)/t
    ENDIF
    exmin1=y
    END FUNCTION exmin1

    FUNCTION exmin1minx(x,eps)
    USE Someconstants    
    IMPLICIT NONE
    !computes (exp(x)-1-x)/(0.5*x*x) 
    REAL(r8) :: exmin1minx, x, eps
    REAL(r8) :: t, t2, y
    IF (x==0) THEN
      y=1.0_r8
    ELSEIF (abs(x)>0.9_r8) THEN
      y=(exp(x)-1.0_r8-x)/(x*x*0.5_r8)
    ELSE
      t=sinh(x*0.5_r8,eps);
      t2=t*t;
      y=(2*t2+(2.0_r8*t*sqrt(1.0_r8+t2)-x))/(x*x*0.5_r8)
    ENDIF
    exmin1minx=y
    END FUNCTION exmin1minx

    FUNCTION lnec(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: lnec, ln1, x,  y0, z, e2, r, s
    !x>-1; lnec:=ln1:=ln(1+x)-x
    z=logoneplusx(x);
    y0=z-x;
    e2=exmin1minx(z,machtol);
    s=e2*z*z*0.5_r8;
    r=(s+y0)/(s+1+z);
    ln1=y0-r*(6.0_r8-r)/(6.0_r8-4.0_r8*r);
    lnec=ln1
    END FUNCTION lnec

    FUNCTION alfa(x)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, alfa, lnx
    lnx=log(x)
    IF (x>0.25_r8) THEN
      alfa=x+0.25_r8
    ELSEIF (x>=dwarf) THEN
      alfa=-0.6931_r8/lnx
    ELSE
      alfa=-0.6931_r8/log(dwarf)
    ENDIF
    END FUNCTION alfa

    FUNCTION  dompart(a,x,qt)
    !dompart is approx. of  x^a*exp(-x)/gamma(a+1)   
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: dompart, a, x
    REAL(r8) :: lnx, c, dp, la, lambda, mu, r 
    LOGICAL :: qt
    lnx=log(x)
    IF (a<=1.0_r8) THEN                     
      r=-x+a*lnx
    ELSE
      IF (x==a) THEN
        r=0
      ELSE
        la=x/a
        r=a*(1.0_r8-la+log(la))
      ENDIF
      r=r-0.5_r8*log(6.2832_r8*a)
    ENDIF
    IF (r<explow) THEN
      dp=0.0_r8
    ELSE
      dp=exp(r)
    ENDIF
    IF (qt) THEN
      dompart=dp
    ELSE
      IF (a<8.0_r8) THEN
        dompart=exp(a*lnx-x)/gamma(a+1.0_r8)
      ELSE
        dompart=1.0_r8/(sqrt(a*twopi)*gamstar(a))
        lambda=x/a
        IF ((lambda>0.3_r8).AND.(lambda<2.36_r8)) THEN
          mu=lambda-1.0_r8
          c=lnec(mu);
          dompart=dompart*exp(a*c)          
        ELSE
          dompart=dompart*exp(a*log(lambda)+a-x)
        ENDIF
      ENDIF
    ENDIF
    END FUNCTION dompart

    FUNCTION eps1(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eps1, eta, la, ak(0:4), bk(0:4) 
    IF (abs(eta)<1.0) THEN
      ak(0)=-3.333333333438e-1_r8;  bk(0)= 1.000000000000e+0_r8;     
      ak(1)=-2.070740359969e-1_r8;  bk(1)= 7.045554412463e-1_r8;     
      ak(2)=-5.041806657154e-2_r8;  bk(2)= 2.118190062224e-1_r8;     
      ak(3)=-4.923635739372e-3_r8;  bk(3)= 3.048648397436e-2_r8;     
      ak(4)=-4.293658292782e-5_r8;  bk(4)= 1.605037988091e-3_r8;     
      eps1=ratfun(eta,ak,bk)
    ELSE
      la=lambdaeta(eta);
      eps1=log(eta/(la-1.0_r8))/eta
    ENDIF
    END FUNCTION eps1
     FUNCTION eps2(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, eps2, x, lnmeta, ak(0:4), bk(0:4)
    IF (eta < -5.0) THEN
      x=eta*eta;
      lnmeta=log(-eta)
      eps2=(12.0_r8-x-6.0*(lnmeta*lnmeta))/(12.0*x*eta)
    ELSEIF (eta<-2.0) THEN
      ak(0)=-1.72847633523e-2_r8;  bk(0)=1.00000000000e+0_r8;     
      ak(1)= -1.59372646475e-2_r8;  bk(1)= 7.64050615669e-1_r8;     
      ak(2)= -4.64910887221e-3_r8;  bk(2)= 2.97143406325e-1_r8;     
      ak(3)= -6.06834887760e-4_r8;  bk(3)= 5.79490176079e-2_r8;     
      ak(4)= -6.14830384279e-6_r8;  bk(4)= 5.74558524851e-3_r8;     
      eps2= ratfun(eta,ak,bk)
    ELSEIF (eta < 2.0) THEN
      ak(0)=-1.72839517431e-2_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)=-1.46362417966e-2_r8;  bk(1)= 6.90560400696e-1_r8;     
      ak(2)=-3.57406772616e-3_r8;  bk(2)= 2.49962384741e-1_r8;     
      ak(3)=-3.91032032692e-4_r8;  bk(3)= 4.43843438769e-2_r8;     
      ak(4)=2.49634036069e-6_r8;   bk(4)= 4.24073217211e-3_r8;     
      eps2= ratfun(eta,ak,bk)
   ELSEIF (eta < 1000.0) THEN
      ak(0)= 9.99944669480e-1_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 1.04649839762e+2_r8;  bk(1)= 1.04526456943e+2_r8;     
      ak(2)= 8.57204033806e+2_r8;  bk(2)= 8.23313447808e+2_r8;     
      ak(3)= 7.31901559577e+2_r8;  bk(3)= 3.11993802124e+3_r8;     
      ak(4)= 4.55174411671e+1_r8;  bk(4)= 3.97003311219e+3_r8;     
      x=1.0_r8/eta;
      eps2=ratfun(x,ak,bk)/(-12.0*eta)
    ELSE
      eps2=-1.0/(12.0*eta)
    ENDIF
    END FUNCTION

    FUNCTION eps3(eta)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, eps3, eta3, x, y, ak(0:4), bk(0:4)
    IF (eta <-8.0) THEN
      x=eta*eta
      y=log(-eta)/eta;
      eps3=(-30.0+eta*y*(6.0_r8*x*y*y-12.0+x))/(12.0*eta*x*x)
    ELSEIF (eta <-4.0) THEN
      ak(0)= 4.95346498136e-2_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 2.99521337141e-2_r8;  bk(1)= 7.59803615283e-1_r8;     
      ak(2)= 6.88296911516e-3_r8;  bk(2)= 2.61547111595e-1_r8;     
      ak(3)= 5.12634846317e-4_r8;  bk(3)= 4.64854522477e-2_r8;     
      ak(4)= -2.01411722031e-5_r8; bk(4)= 4.03751193496e-3_r8;     
      eps3=ratfun(eta,ak,bk)/(eta*eta)
    ELSEIF (eta <-2.0) THEN
      ak(0)=4.52313583942e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)=1.20744920113e-3_r8;  bk(1)= 9.12203410349e-1_r8;     
      ak(2)=-7.89724156582e-5_r8; bk(2)= 4.05368773071e-1_r8;     
      ak(3)=-5.04476066942e-5_r8; bk(3)= 9.01638932349e-2_r8;     
      ak(4)=-5.35770949796e-6_r8; bk(4)= 9.48935714996e-3_r8;     
      eps3=ratfun(eta,ak,bk)
    ELSEIF  (eta < 2.0) THEN
      ak(0)= 4.39937562904e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= 4.87225670639e-4_r8;  bk(1)= 7.94435257415e-1_r8;     
      ak(2)= -1.28470657374e-4_r8; bk(2)= 3.33094721709e-1_r8;     
      ak(3)= 5.29110969589e-6_r8;  bk(3)= 7.03527806143e-2_r8;     
      ak(4)= 1.57166771750e-7_r8;  bk(4)= 8.06110846078e-3_r8;     
      eps3= ratfun(eta,ak,bk)
    ELSEIF (eta < 10.0) THEN
      ak(0)= -1.14811912320e-3_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= -1.12850923276e-1_r8;  bk(1)= 1.42482206905e+1_r8;     
      ak(2)= 1.51623048511e+0_r8;   bk(2)= 6.97360396285e+1_r8;     
      ak(3)= -2.18472031183e-1_r8;  bk(3)= 2.18938950816e+2_r8;     
      ak(4)= 7.30002451555e-2_r8;   bk(4)= 2.77067027185e+2_r8;     
      x= 1.0_r8/eta;
      eps3= ratfun(x,ak,bk)/(eta*eta)
    ELSEIF (eta < 100.0) THEN
      ak(0)= -1.45727889667e-4_r8;  bk(0)= 1.00000000000e+0_r8;     
      ak(1)= -2.90806748131e-1_r8;  bk(1)= 1.39612587808e+2_r8;     
      ak(2)= -1.33085045450e+1_r8;  bk(2)= 2.18901116348e+3_r8;     
      ak(3)= 1.99722374056e+2_r8;   bk(3)= 7.11524019009e+3_r8;     
      ak(4)= -1.14311378756e+1_r8;  bk(4)= 4.55746081453e+4_r8;     
      x= 1.0_r8/eta;
      eps3= ratfun(x,ak,bk)/(eta*eta)
    ELSE
     eta3=eta*eta*eta
     eps3=-log(eta)/(12.0*eta3)
    ENDIF
    END FUNCTION eps3

    FUNCTION lambdaeta(eta)
! lambdaeta is the positive number satisfying 
! eta^2/2=lambda-1-ln(lambda)
! with sign(lambda-1)=sign(eta); 
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: eta, lambdaeta, ak(6), q, r, s, L, la
    REAL(r8) :: L2, L3, L4, L5
    s=eta*eta*0.5_r8
    IF (eta==0) THEN
      la=1 
    ELSEIF (eta < -1) THEN 
      r=exp(-1-s);
      ak(1)=1.0_r8;
      ak(2)=1.0_r8;
      ak(3)=1.5_r8;
      ak(4)=2.66666666666666666666666666667_r8;
      ak(5)=5.20833333333333333333333333333_r8;
      ak(6)=10.8_r8;
      la=r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ELSEIF (eta<1) THEN 
      ak(1)=1.0_r8;
      ak(2)=0.333333333333333333333333333333_r8;
      ak(3)=0.0277777777777777777777777777778_r8;
      ak(4)=-0.00370370370370370370370370370370_r8;
      ak(5)=0.000231481481481481481481481481481_r8;
      ak(6)=0.0000587889476778365667254556143445_r8;
      r=eta;
      la=1.0_r8+r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ELSE
      r=11+s; L=log(r); la=r+L; r=1.0_r8/r;
      L2=L*L
      L3=L2*L
      L4=L3*L 
      L5=L4*L
      ak(1)= 1;
      ak(2)=(2-L)*0.5_r8;
      ak(3)=(-9*L+6+2*L2)*0.166666666666666666666666666667_r8;
      ak(4)= -(3*L3+36*L-22*L2-12)*0.0833333333333333333333333333333_r8;
      ak(5)=(60+350*L2-300*L-125*L3+12*L4)*0.0166666666666666666666666666667_r8;
      ak(6)=-(-120-274*L4+900*L-1700*L2+1125*L3+20*L5)*&
            0.00833333333333333333333333333333_r8;
      la=la+L*r*(ak(1)+r*(ak(2)+r*(ak(3)+r*(ak(4)+r*(ak(5)+r*ak(6))))))
    ENDIF
    r= 1;
    IF (((eta>-3.5_r8).AND.(eta<-0.03_r8)).OR.((eta>0.03_r8).AND.(eta<40.0_r8))) THEN
      r=1.0_r8; 
      q=la;
      DO WHILE (r>1.0e-8_r8)
        la=q*(s+log(q))/(q-1.0_r8);
        r=abs(q/la-1.0_r8); 
        q=la
      ENDDO
    ENDIF
    lambdaeta=la
    END FUNCTION lambdaeta

    RECURSIVE FUNCTION auxloggam(x) RESULT(auxloggamm)
    !function g in ln(Gamma(1+x))=x*(1-x)*g(x), 0<=x<=1}
    USE Someconstants
    IMPLICIT NONE    
    REAL(r8) :: auxloggamm, x
    REAL(r8) :: ak(0:25)
    REAL(r8) :: g, t
    IF (x<-1.0_r8) THEN 
      g=giant
    ELSEIF (abs(x)<=dwarf) THEN 
      g=-eulmasc
    ELSEIF (abs(x - 1.0_r8)<=machtol) THEN
      g=eulmasc-1.0_r8
    ELSEIF (x<0) THEN
      g=-(x*(1+x)*auxloggam(x+1.0_r8)+logoneplusx(x))/(x*(1.0_r8-x))
    ELSEIF (x<1) THEN
      ak(0)=-0.98283078605877425496_r8;
      ak(1)=0.7611416167043584304e-1_r8;
      ak(2)=-0.843232496593277796e-2_r8;
      ak(3)=0.107949372632860815e-2_r8;
      ak(4)=-0.14900748003692965e-3_r8;
      ak(5)=0.2151239988855679e-4_r8;
      ak(6)=-0.319793298608622e-5_r8;
      ak(7)=0.48516930121399e-6_r8;
      ak(8)=-0.7471487821163e-7_r8;
      ak(9)=0.1163829670017e-7_r8;
      ak(10)=-0.182940043712e-8_r8;
      ak(11)= 0.28969180607e-9_r8;
      ak(12)=-0.4615701406e-10_r8;
      ak(13)= 0.739281023e-11_r8;
      ak(14)= -0.118942800e-11_r8;
      ak(15)= 0.19212069e-12_r8;
      ak(16)= -0.3113976e-13_r8;
      ak(17)= 0.506284e-14_r8;
      ak(18)= -0.82542e-15_r8;
      ak(19)= 0.13491e-15_r8;
      ak(20)= -0.2210e-16_r8;
      ak(21)= 0.363e-17_r8;
      ak(22)= -0.60e-18_r8;
      ak(23)= 0.98e-19_r8;
      ak(24)= -0.2e-19_r8;
      ak(25)= 0.3e-20_r8;
      t=2*x-1;
      g=chepolsum(25, t, ak)
    ELSEIF (x<1.5) THEN
      g=(logoneplusx(x-1.0_r8)+(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8))/(x*(1.0_r8-x));
    ELSE
      g=(log(x)+(x-1.0_r8)*(2.0_r8-x)*auxloggam(x-1.0_r8))/(x*(1.0_r8-x));
    ENDIF
    auxloggamm= g
    END FUNCTION auxloggam

    RECURSIVE FUNCTION auxgam(x) RESULT(auxgamm)
    USE Someconstants
    IMPLICIT NONE    
    !function g in 1/gamma(x+1)=1+x*(x-1)*g(x), -1<=x<=1}
    REAL(r8) :: auxgamm, x
    REAL(r8) :: t, dr(0:17)
    IF (x<0.0_r8) THEN
      auxgamm=-(1.0_r8+(1.0_r8+x)*(1.0_r8+x)*auxgam(1.0_r8+x))/(1.0_r8-x)
    ELSE
      dr(0)= -1.013609258009865776949_r8;
      dr(1)= 0.784903531024782283535e-1_r8;
      dr(2)= 0.67588668743258315530e-2_r8;
      dr(3)= -0.12790434869623468120e-2_r8;
      dr(4)= 0.462939838642739585e-4_r8;
      dr(5)= 0.43381681744740352e-5_r8;
      dr(6)= -0.5326872422618006e-6_r8;
      dr(7)= 0.172233457410539e-7_r8;
      dr(8)= 0.8300542107118e-9_r8;
      dr(9)= -0.10553994239968e-9_r8;
      dr(10)= 0.39415842851e-11_r8;
      dr(11)= 0.362068537e-13_r8;
      dr(12)= -0.107440229e-13_r8;
      dr(13)= 0.5000413e-15_r8;
      dr(14)= -0.62452e-17_r8;
      dr(15)= -0.5185e-18_r8;
      dr(16)= 0.347e-19_r8;
      dr(17)= -0.9e-21_r8;
      t=2*x-1.0_r8;
      auxgamm=chepolsum(17,t,dr);
     ENDIF
    END FUNCTION auxgam

    FUNCTION lngam1(x)
    !ln(gamma(1+x)), -1<=x<=1}
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: x, lngam1   
    lngam1=-logoneplusx(x*(x-1.0_r8)*auxgam(x))
    END FUNCTION lngam1

    FUNCTION stirling(x)            
    !Stirling series, function corresponding with}
    !asymptotic series for log(gamma(x))}
    !that is:  1/(12x)-1/(360x**3)...; x>= 3}
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: stirling, x, a(0:17), c(0:6), z
    IF (x<dwarf) THEN
      stirling=giant
    ELSEIF (x<1.0_r8) THEN
      stirling= lngam1(x)-(x+0.5_r8)*log(x)+x-lnsqrttwopi
    ELSEIF (x<2.0_r8) THEN
      stirling=lngam1(x-1.0_r8)-(x-0.5_r8)*log(x)+x-lnsqrttwopi
    ELSEIF (x<3.0_r8) THEN
      stirling=lngam1(x-2.0_r8)-(x-0.5_r8)*log(x)+x&
               -lnsqrttwopi+log(x-1)
    ELSEIF (x<12.0_r8) THEN
      a(0)=1.996379051590076518221_r8;
      a(1)=-0.17971032528832887213e-2_r8;
      a(2)=0.131292857963846713e-4_r8;
      a(3)=-0.2340875228178749e-6_r8;
      a(4)=0.72291210671127e-8_r8;
      a(5)=-0.3280997607821e-9_r8;
      a(6)=0.198750709010e-10_r8;
      a(7)=-0.15092141830e-11_r8;
      a(8)=0.1375340084e-12_r8;
      a(9)=-0.145728923e-13_r8;
      a(10)=0.17532367e-14_r8;
      a(11)=-0.2351465e-15_r8;
      a(12)=0.346551e-16_r8;
      a(13)=-0.55471e-17_r8;
      a(14)=0.9548e-18_r8;
      a(15)=-0.1748e-18_r8;
      a(16)=0.332e-19_r8;
      a(17)=-0.58e-20_r8;
      z=18.0_r8/(x*x)-1.0_r8;
      stirling=chepolsum(17,z,a)/(12.0_r8*x);
    ELSE
      z=1.0_r8/(x*x);
      IF (x<1000.0_r8) THEN
        c(0)=0.25721014990011306473e-1_r8;
        c(1)=0.82475966166999631057e-1_r8;
        c(2)=-0.25328157302663562668e-2_r8;
        c(3)=0.60992926669463371e-3_r8;
        c(4)=-0.33543297638406e-3_r8;
        c(5)=0.250505279903e-3_r8;
        c(6)=0.30865217988013567769_r8;
        stirling=((((((c(5)*z+c(4))*z+c(3))*z+c(2))*z+c(1))*z+c(0))/(c(6)+z)/x)
      ELSE
        stirling=(((-z*0.000595238095238095238095238095238_r8+&
                 0.000793650793650793650793650793651_r8)*z&
                -0.00277777777777777777777777777778_r8)*z+&
                 0.0833333333333333333333333333333_r8)/x
      ENDIF
    ENDIF 
    END FUNCTION stirling

    FUNCTION ratfun(x,ak,bk)
    USE Someconstants
    IMPLICIT NONE
    REAL(r8) :: ratfun, x, ak(0:4), bk(0:4), p, q
    p= ak(0)+x*(ak(1)+x*(ak(2)+x*(ak(3)+x*ak(4))));
    q= bk(0)+x*(bk(1)+x*(bk(2)+x*(bk(3)+x*bk(4))));
    ratfun=p/q
    END FUNCTION ratfun
 END MODULE GammaError


 
