  MODULE IncgamNEG
  USE set_precision
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: incgamstar
  CONTAINS     
     SUBROUTINE incgamstar(a,z,igam,ierr)
    ! ------------------------------------------------------------------
    ! Calculation of the incomplete gamma star function for negative
    ! real values of the argument z and real values (positive
    ! and negative) of the a parameter.
    !
    ! The incomplete gamma star function is defined as
    ! gamma*(a,z)=(z**(-a)/GAMMA(a))*gamma(a,z),
    ! where gamma(a,z) is the lower incomplete
    ! gamma funcion.
    !
    ! -------------------------------------------------------------------
    ! Inputs:
    !   a ,    parameter of the function
    !   z ,    argument of the function 
    ! Outputs:
    !   igam,  incomplete gamma star function. 
    !   ierr , error flag
    !          ierr=0, computation succesful.
    !          ierr=1, overflow/underflow problems. The function
    !                  value is set to 0. 
    !          ierr=2, argument z out of range.
    !                  The function value is set to 0. 
    ! -------------------------------------------------------------------
    USE Someconstants    
    USE GammaError
    IMPLICIT NONE
    REAL(r16), INTENT(IN) :: a
    REAL(r8),  INTENT(IN) :: z
    REAL(r8),  INTENT(OUT) :: igam 
    INTEGER,   INTENT(OUT) :: ierr
    REAL(r8) :: aa, zz, epsi, p, q, sinpia, cospia, igamnew  
    REAL(r16) :: epsic
    INTEGER :: i, j
    ierr=0
    IF (z>0) THEN
      ierr=2
      igam=0.0_r8
    ENDIF
    IF (ierr==0) THEN
      IF (a>0) THEN
        aa=a
        zz=-z
        IF (zz>50.0_r8) THEN
          CALL gexpan(aa,zz,igam,ierr)
        ELSE
          CALL gseries(aa,z,igam,ierr)
        ENDIF 
      ELSE
        aa=-a
        zz=-z
        IF (abs(nint(a)-a)<dwarf) THEN
          IF (aa*log(zz)>loggiant) THEN
            igam=0.0_r8
            ierr=1
          ELSE
            igam=(-1)**int(aa)*zz**aa
          ENDIF
        ELSE 
          epsic=a-nint(a);
          epsi=epsic
          IF ((aa<5.0_r8).OR.(zz<1.5_r8)) THEN
            IF (zz<100.0_r8) THEN
            ! Series expansion
              IF (abs(epsi)>0.1_r8) THEN
                CALL gseries(-aa,z,igam,ierr) 
              ELSE
                CALL gsereps(-aa,-epsi,z,igam,ierr)
              ENDIF
            ELSE 
              ! scaled UAE+Recursion
              j=20.0_r8-aa+1
              CALL incgamnegascal(aa+j,zz,igamnew,ierr)
              q=igamnew
              DO i=j-1,0,-1
                p=((aa+i)/zz)*(q-1.0_r8/pi);
                q=p;
              ENDDO
              p=q
              IF (abs(epsi)<0.1_r8) THEN
                sinpia=(-1)**nint(aa)*sin(-pi*epsi)
                cospia=cos(pi*aa)
              ELSE
                sinpia=sin(pi*aa)
                cospia=cos(pi*aa)
              ENDIF
              IF ((aa*log(zz)<loggiant).AND.(zz<loggiant)) THEN
                igam=exp(aa*log(zz))*cospia+(sinpia*exp(zz))*gamma(aa)*p
              ELSE
                igam=0.0_r8
                ierr=1
              ENDIF
            ENDIF
          ELSE
          ! Uniform Asymptotic Expansion
            CALL incgamnegaeps(aa,zz,-epsi,igam,ierr) 
          ENDIF
        ENDIF 
      ENDIF
    ENDIF
    END SUBROUTINE incgamstar

    SUBROUTINE incgamnegaeps(a,x,epsil,igam,ierr)
    ! ---------------------------------------------------------------
    ! Calculation of the incomplete gamma function for negative
    ! input values using UAE
    ! ---------------------------------------------------------------
    ! Inputs:
    !   a ,    argument of the functions
    !   x ,    argument of the functions
    !   epsil, parameter in a=-n+epsil. 
    ! Outputs:
    !   igam,  incomplete gamma function 
    !   ierr , error flag
    !          ierr=0, computation succesful.
    !          ierr=1, overflow/underflow problems. The function
    !                  value is set to 0. 
    !          ierr=2, possible loss of accuracy in the computation
    !                  of the function value.
    ! -------------------------------------------------------------------
    USE Someconstants    
    USE GammaError
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(IN) :: epsil
    REAL(r8), INTENT(OUT) :: igam 
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: lambda, mu, eta, p, q, sinpia, cospia, alogx
    ierr=0
    eta=0.0_r8
    lambda=x/a; mu= lambda-1.0_r8;
    IF (mu==0.0_r8) THEN
      eta=0.0_r8 
    ELSEIF (mu > 0.0_r8) THEN
      IF (mu>0.1_r8) THEN
        eta=sqrt(2.0_r8*(mu-log(1.0_r8+mu)))
      ELSE 
        eta=sqrt(2.0_r8*(-lnec(mu)))
      ENDIF 
    ELSE
      eta=-sqrt(2.0_r8*(-lnec(mu))) 
    ENDIF   
    sinpia=(-1)**nint(a)*sin(pi*epsil)
    cospia=cos(pi*a)
    IF (a*(mu-log(1.0_r8+mu))>loggiant) THEN
      igam=0.0_r8
      ierr=1
    ELSE
      p=sqrt(2.0_r8*a/pi)*exp(a*(mu-log(1.0_r8+mu)))*sinpia;
      q=sqrt(2.0_r8/a)*fz(eta*sqrt(a/2.0_r8))+Taeta(a,mu,eta)/a;
      alogx=a*log(x)
      IF ((alogx<log(dwarf)).OR.(alogx>loggiant)) THEN
        igam=0.0_r8
        ierr=1
      ELSE
        IF ((alogx+log(abs(cospia-p*q)))>loggiant) THEN
          igam=0.0_r8
          ierr=1
        ELSE
          igam=exp(alogx)*(cospia-p*q)
        ENDIF
      ENDIF

    ENDIF
    END SUBROUTINE incgamnegaeps


    SUBROUTINE incgamnegascal(a,x,igam,ierr)
    ! ---------------------------------------------------------------
    ! Calculation of the scaled incomplete gamma function for negative
    ! input values using UAE
    ! ---------------------------------------------------------------
    ! Inputs:
    !   a ,    argument of the functions
    !   x ,    argument of the functions
    ! Outputs:
    !   igam,  incomplete gamma function (plain or scaled) 
    !   ierr , error flag
    !          ierr=0, computation succesful
    !          ierr=1, overflow/underflow problems. 
    ! -------------------------------------------------------------
    USE Someconstants    
    USE GammaError
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: a
    REAL(r8), INTENT(IN) :: x
    REAL(r8), INTENT(OUT) :: igam 
    INTEGER,  INTENT(OUT) :: ierr
    REAL(r8) :: lambda, mu, eta, p, q
    ierr=0
    eta=0.0_r8
    lambda=x/a; mu= lambda-1.0_r8;
    IF (mu==0.0_r8) THEN
      eta=0.0_r8 
    ELSEIF (mu > 0.0_r8) THEN
      IF (mu>0.1_r8) THEN
        eta=sqrt(2.0_r8*(mu-log(1.0_r8+mu)))
      ELSE 
        eta=sqrt(2.0_r8*(-lnec(mu)))
      ENDIF 
    ELSE
      eta=-sqrt(2.0_r8*(-lnec(mu))) 
    ENDIF 
    ! scaled functions
    p=-a/pi/gamstar(a);
    q=sqrt(2.0_r8/a)*fz(eta*sqrt(a/2.0_r8))+Taeta(a,mu,eta)/a;
    igam=p*q     
    END SUBROUTINE incgamnegascal

    FUNCTION fz(z)
    ! Dawson integral 
    ! F(z) continued fraction. Cuyt et al. (13.1.13b)
    USE Someconstants    
    IMPLICIT NONE
    REAL(r8) :: fz, z
    REAL(r8) :: a, b, Ak(0:1000), Bk(0:1000), r0, r1, z2
    INTEGER :: k
    z2=z*z
    a=1.0_r8; b= 1.0_r8+2*z2; 
    Ak(0)= 0.0_r8; Bk(0)= 1.0_r8;
    Ak(1)= a; Bk(1)= b;
    r0= 1.0_r8; r1= 2.0_r8; k= 2;
    DO WHILE (abs((r1/r0)-1.0_r8)>1.0e-17_r8) 
      r0= r1;
      a=-4.0_r8*(k-1.0_r8)*z2/((2*k-3.0_r8)*(2.0_r8*k-1.0_r8));
      b=1.0_r8+2*z2/(2*k-1.0_r8);
      Ak(k)=b*Ak(k-1)+a*Ak(k-2);
      Bk(k)= b*Bk(k-1)+a*Bk(k-2);
      r1= Ak(k)/Bk(k);
      k= k+1;
    ENDDO
    fz=z*r1
    END FUNCTION fz

    FUNCTION gamstar(a) 
    ! gk[k]=(-1)^k\gamma_k
    REAL(r8) :: gamstar, a
    REAL(r8) :: s, b, t, gk(0:20)
    INTEGER :: k
    gk(0)= 1.0_r8;
    gk(1)=0.00347222222222222222222222222222_r8;
    gk(2)=-0.00268132716049382716049382716049_r8;
    gk(3)=-0.000229472093621399176954732510288_r8;
    gk(4)=0.000784039221720066627474034881442_r8;
    gk(5)=0.0000697281375836585777429398828576_r8;
    gk(6)=-0.000592166437353693882864836225604_r8;
    gk(7)=-0.0000517179090826059219337057843002_r8;
    gk(8)=0.000839498720672087279993357516765_r8;
    gk(9)=0.0000720489541602001055908571930225_r8;
    gk(10)=-0.00191443849856547752650089885833_r8;
    gk(11)=-0.000162516262783915816898635123980_r8;
    gk(12)=0.00640336283380806979482363809027_r8;
    gk(13)=0.000540164767892604515180467508570_r8;
    gk(14)=-0.0295278809456991205054406510547_r8;
    gk(15)=-0.00248174360026499773091565836874_r8;
    gk(16)=0.179540117061234856107699407722_r8;
    gk(17)=0.0150561130400264244123842218771_r8;
    gk(18)=-1.39180109326533748139914776354_r8;
    gk(19)=-0.116546276599463200850734036907_r8;
    s= gk(0)
    b= 1.0_r8;  k=0; t=s;
    DO WHILE ((abs(t/s) > 1.0e-20_r8).AND.(k<20))
      b=b/a; k=k+1; t=gk(k)*b; s=s + t; 
    ENDDO   
    gamstar=s 
    END FUNCTION gamstar

    FUNCTION Taetarec(a, eta) 
    ! T_a(eta) for abs(eta)<1
    REAL(r8) :: Taetarec, a, eta
    REAL(r8) :: y, s, t, fm(0:32), bm(0:32);
    INTEGER ::  m               
    fm(0)= 1.0_r8;
    fm(1)=-0.333333333333333333333333333333_r8;
    fm(2)=0.0833333333333333333333333333333_r8;
    fm(3)=-0.0148148148148148148148148148148_r8;
    fm(4)=0.00115740740740740740740740740741_r8;
    fm(5)=0.000352733686067019400352733686067_r8;
    fm(6)=-0.000178755144032921810699588477366_r8;
    fm(7)=0.0000391926317852243778169704095630_r8
    fm(8)=-0.00000218544851067999216147364295512_r8;
    fm(9)=-0.00000185406221071515996070179883623_r8                                      
    fm(10)=8.29671134095308600501624213166e-7_r8;                                     
    fm(11)=-1.76659527368260793043600542457e-7_r8;                                     
    fm(12)=6.70785354340149858036939710030e-9_r8;                                     
    fm(13)=1.02618097842403080425739573227e-8;                                     
    fm(14)=-4.38203601845335318655297462245e-9_r8;                                    
    fm(15)=9.14769958223679023418248817633e-10_r8                                     
    fm(16)=-2.55141939949462497668779537994e-11_r8;                                     
    fm(17)=-5.83077213255042506746408945040e-11_r8                                    
    fm(18)=2.43619480206674162436940696708e-11_r8                                     
    fm(19)=-5.02766928011417558909054985926e-12_r8                                    
    fm(20)=1.10043920319561347708374174497e-13_r8                                    
    fm(21)=3.37176326240098537882769884169e-13_r8               
    fm(22)=-1.39238872241816206591936618490e-13_r8                                    
    fm(23)=2.85348938070474432039669099053e-14_r8                                     
    fm(24)=-5.13911183424257261899064580300e-16_r8                                     
    fm(25)=-1.97522882943494428353962401581e-15_r8                                    
    fm(26)=8.09952115670456133407115668703e-16_r8
    fm(27)=-.1652253122e-15_r8;
    fm(28)=.2530543010e-17_r8;
    fm(29)=.1168693974e-16_r8;
    fm(30)=-.4770037050e-17_r8;
    fm(31)=.9699126059e-18_r8;
    fm(32)=-.1293256554e-19_r8;
    bm(31)= fm(32);
    bm(30)= fm(31);
    DO m=30,1,-1 
       bm(m-1)= fm(m)-(m + 1)*bm(m + 1)/a
    ENDDO
    s=bm(0);
    y=1.0_r8;  m=0; t=s;
    DO WHILE ((abs(t/s) > 1.0e-20_r8).AND.(m < 31))
      y=y*eta; m= m+1;
      t= bm(m)*y;
      s= s + t;
    ENDDO 
    Taetarec=s/(1.0_r8-bm(1)/a)
    END FUNCTION Taetarec

    FUNCTION Taeta(a, mu, eta) 
    ! T_a(eta) for abs(eta)>1 
    USE Someconstants  
    IMPLICIT NONE
    REAL(r8) :: Taeta, a, mu, eta 
    REAL(r8) :: b, s, t, Ck(0:15)
    REAL(r8) :: mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10, mu11,&
                mu12, mu13, mu14, mu15, mu16, mu17, mu18, mu19, mu20,&
                mu21, mu22, mu23, mu24, mu25, mu26, mu27, mu28, mu29, mu30,&
                mu31
    REAL(r8) :: eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, eta10, eta11, &
                eta12, eta13, eta14, eta15, eta16, eta17, eta18, eta19, eta20, &
                eta21, eta22, eta23, eta24, eta25, eta26, eta27, eta28, eta29, eta30, &
                eta31
    INTEGER :: k, irec
    IF (a<4.0_r8) THEN
      IF ((mu>-1.0_r8).AND.(mu < 1.8_r8)) THEN
        IF (abs(eta)<1.5_r8) THEN
          irec=1
        ELSE
          irec=0
        ENDIF
      ELSE
        irec=0
      ENDIF
    ELSEIF (a<11.0_r8) THEN
      IF ((a*(mu+1.0_r8))>0.0_r8) THEN
        IF ((mu>-1.0_r8).AND.(mu <1.8_r8)) THEN
          IF (abs(eta)<1.5_r8) THEN
            irec=1
           ELSE
            irec=0
          ENDIF
        ELSE
          irec=0
        ENDIF
      ELSE
        irec=0
      ENDIF 
    ELSEIF (a<13.0_r8) THEN
      IF ((a*(mu+1.0_r8))>1.0_r8) THEN
        IF ((mu>-0.9_r8).AND.(mu <1.8_r8)) THEN
          IF (abs(eta)<1.5_r8) THEN
            irec=1
           ELSE
            irec=0
          ENDIF
        ELSE
          irec=0
        ENDIF
      ELSE
        irec=0
      ENDIF 
    ELSE
      IF ((mu>-0.7_r8).AND.(mu < 1.35_r8)) THEN
        IF (abs(eta)<1.5_r8) THEN
          irec=1
        ELSE
          irec=0
        ENDIF
      ELSE
        irec=0
      ENDIF
    ENDIF
    IF (irec==1) THEN
      t=Taetarec(a,eta) 
    ELSE
      mu2=mu*mu
      mu3=mu2*mu
      mu4=mu3*mu
      mu5=mu4*mu
      mu6=mu5*mu
      mu7=mu6*mu
      mu8=mu7*mu
      mu9=mu8*mu
      mu10=mu9*mu
      mu11=mu10*mu
      mu12=mu11*mu
      mu13=mu12*mu
      mu14=mu13*mu
      mu15=mu14*mu
      mu16=mu15*mu; 
      mu17=mu16*mu; 
      mu18=mu17*mu; 
      mu19=mu18*mu; 
      mu20=mu19*mu
      mu21=mu20*mu; 
      mu22=mu21*mu; 
      mu23=mu22*mu; 
      mu24=mu23*mu; 
      mu25=mu24*mu;
      mu26=mu25*mu; 
      mu27=mu26*mu; 
      mu28=mu27*mu; 
      mu29=mu28*mu; 
      mu30=mu29*mu
      mu31=mu30*mu
      eta2=eta*eta
      eta3=eta2*eta
      eta4=eta3*eta
      eta5=eta4*eta
      eta6=eta5*eta
      eta7=eta6*eta
      eta8=eta7*eta
      eta9=eta8*eta
      eta10=eta9*eta
      eta11=eta10*eta
      eta12=eta11*eta
      eta13=eta12*eta
      eta14=eta13*eta
      eta15=eta14*eta
      eta16=eta15*eta
      eta17=eta16*eta
      eta18=eta17*eta
      eta19=eta18*eta
      eta20=eta19*eta
      eta21=eta20*eta
      eta22=eta21*eta
      eta23=eta22*eta
      eta24=eta23*eta
      eta25=eta24*eta
      eta26=eta25*eta
      eta27=eta26*eta
      eta28=eta27*eta
      eta29=eta28*eta
      eta30=eta29*eta
      eta31=eta30*eta  
      Ck(0)=1.0_r8/mu-1.0_r8/eta;
      Ck(1)=-0.833333333333333333333333333333e-1_r8*(12.0_r8*mu+12.0_r8+mu2)/mu3+1.0_r8/eta3;
      Ck(2)=0.347222222222222222222222222222e-2_r8*(600.0_r8*mu2+24.0_r8*mu3+&
            1440.0_r8*mu+864.0_r8+mu4)/mu5-3.0_r8/eta5;
      Ck(3)=0.192901234567901234567901234568e-4_r8*(-332640.0_r8*mu3-8820.0_r8*mu4-1360800.0_r8*mu2&
             -180.0_r8*mu5-1814400.0_r8*mu-777600.0_r8+139*mu6)/mu7+15.0_r8/eta7;
      Ck(4)=-0.401877572016460905349794238683e-6_r8*(-65136960.0_r8*mu4-1287360.0_r8*mu5&
             -390458880.0_r8*mu3-10608.0_r8*mu6-849139200.0_r8*mu2+6672.0_r8*mu7&
             -783820800.0_r8*mu-261273600.0_r8+571.0_r8*mu8)/mu9-105.0_r8/eta9;
      Ck(5)=-0.478425680971977268273564569861e-8_r8*(27790076160.0_r8*mu5+&
             435226176.0_r8*mu6+224148798720.0_r8*mu4+1552320.0_r8*mu7+&
             696085125120.0_r8*mu3-1168860.0_r8*mu8+1026021427200.0_r8*mu2&
             -47964.0_r8*mu9+724250419200.0_r8*mu+197522841600.0_r8+163879.0_r8*mu10)&
              /mu11+945.0_r8/eta11;
      Ck(6)=0.132896022492215907853767936072e-10_r8*(60809971622400.0_r8*mu6+785642457600.0_r8*mu7&
             +624881537280000.0_r8*mu5+972972000.0_r8*mu8+2569580133120000.0_r8*mu4&
             -1296902880.0_r8*mu9+5329034584473600.0_r8*mu3+24462360.0_r8*mu10&
             +5931610933248000.0_r8*mu2+58996440.0_r8*mu11+3389491961856000.0_r8*mu&
             +782190452736000.0_r8+5246819.0_r8*mu12)/mu13-10395.0_r8/eta13;
      Ck(7)=0.110746685410179923211473280060e-11_r8*(-610108553134080000.0_r8*mu&
             -1422456793325568000.0_r8*mu3&
            -1271059485696000000.0_r8*mu2-337503281955840000.0_r8*mu5&
            -916998804513792000.0_r8*mu4-5164603873228800.0_r8*mu7-65096665195161600.0_r8*mu6&
            +3873018240.0_r8*mu9-56624635267200.0_r8*mu8-2296559520.0_r8*mu11+61370693280.0_r8*mu10&
            -62961828.0_r8*mu13-1478876388.0_r8*mu12-122021710626816000.0_r8+&
            534703531.0_r8*mu14)/mu15+135135.0_r8/eta15;
      Ck(8)=-0.115361130635604086678618000063e-13_r8*(-175711263302615040000.0_r8&
             -995697158714818560000.0_r8*mu&
             -3224952464059662336000.0_r8*mu3-2406268133560811520000.0_r8*mu2&
             -1292353888244170752000.0_r8*mu5-2607020963477618688000.0_r8*mu4&
             -60209934503259340800.0_r8*mu7-380246669406226022400.0_r8*mu6&
             -38049524041052160.0_r8*mu9-4004467529539276800.0_r8*mu8+28576053918720.0_r8*mu11&
             +31688791280640.0_r8*mu10-438005070720.0_r8*mu13-1307795255424.0_r8*mu12&
             +51331538976.0_r8*mu15+39242868000.0_r8*mu14+4483131259.0_r8*mu16)/mu17&
             -2027025.0_r8/eta17;
      Ck(9)=-0.194210657635697115620569023675e-17_r8*(112374381332554422681600000.0_r8*mu&
             +482585537611469826293760000.0_r8*mu3+309029548664524662374400000.0_r8*mu2&
             +293433143953786206289920000.0_r8*mu5+469501205710943234949120000.0_r8*mu4&
             +28421787488496413147136000.0_r8*mu7+116964302533047299506176000.0_r8*mu6&
             +215886947511600536371200.0_r8*mu9+3790548943622774581248000.0_r8*mu8&
             -2336070503112192000.0_r8*mu11+1806795762489349632000.0_r8*mu10&
             +49248519566400000.0_r8*mu13-979609042577088000.0_r8*mu12&
             -1309126590794880.0_r8*mu15+9707692572547200.0_r8*mu14-26629799678460.0_r8*mu17&
             -636448482713340.0_r8*mu16+17743323368298066739200000.0_r8+432261921612371.0_r8*mu18)&
             /mu19+34459425.0_r8/eta19;
      Ck(10)=0.115601581926010187869386323616e-19_r8*(56636688191607429031526400000.0_r8&
             +396456817341252003220684800000.0_r8*mu+2179778316085513328818913280000.0_r8*mu3&
             +1222408520135527009930444800000.0_r8*mu2+1873297792970882109265674240000.0_r8*mu5&
             +2480332963491207845149409280000.0_r8*mu4+312748360308984071026114560000.0_r8*mu7&
             +945604791311480459315380224000.0_r8*mu6+7367624519634376327249920000.0_r8*mu9&
             +64303258024623656330772480000.0_r8*mu8+2728735514127713857536000.0_r8*mu11&
             +365421947012372797747200000.0_r8*mu10-1102377726347724288000.0_r8*mu13&
             -4291698990253441536000.0_r8*mu12+7274728691925488640.0_r8*mu15+57796969483870848000.0_r8*mu14&
             -329717647979485920.0_r8*mu17-1200503104301682720.0_r8*mu16+72620002830878328.0_r8*mu19&
             +63672390138915768.0_r8*mu18+6232523202521089.0_r8*mu20)/mu21-654729075.0_r8/eta21;
      Ck(11)=0.741035781576988383778117459077e-22_r8*(-185541790515705937507280486400000.0_r8&
             -1422487060620412187555817062400000.0_r8*mu-9744036365249823484757346877440000.0_r8*mu3&
             -4860164123786408307482374963200000.0_r8*mu2-11253594310434004950063356313600000.0_r8*mu5&
             -12698660530746804632659741900800000.0_r8*mu4-2895757630583685485919023923200000.0_r8*mu7&
             -6888466502924186019525112627200000.0_r8*mu6-144199200374692327523418439680000.0_r8*mu9&
             -813449426188757986159840788480000.0_r8*mu8-631320888475270954551951360000.0_r8*mu11&
             -14419257161828784205546045440000.0_r8*mu10+7401312784797791827968000.0_r8*mu13&
             -4250801856656917785839616000.0_r8*mu12-69923436732029223383040.0_r8*mu15&
             +1312653111805572945408000.0_r8*mu14+1142136233694511735680.0_r8*mu17&
             -5872753634286944845440.0_r8*mu16-52456119468246617760.0_r8*mu19+&
              175945133754186634656.0_r8*mu18-972273619593289884.0_r8*mu21&
              -23629714502827328220.0_r8*mu20&
              +25834629665134204969.0_r8*mu22)/mu23+13749310575.0_r8/eta23;
       Ck(12)=-0.102921636330137275524738535983e-24_r8*(-3072572050940090325120564854784000000.0_r8&
              -25604767091167419376004707123200000000.0_r8*mu&
              -213799805211247951789639304478720000000.0_r8*mu3&
              -96017876591877822660017651712000000000.0_r8*mu2&
              -319564258323840991527606367027200000000.0_r8*mu5&
              -314031799720213745555291064238080000000.0_r8*mu4&
              -117673957900116093676774534152192000000.0_r8*mu7&
              -230161412259016741031808476381184000000.0_r8*mu6&
              -10238781742615484911962338387558400000.0_r8*mu9&
              -42144381707162643048313436110848000000.0_r8*mu8&
              -140418859511143761211239845068800000.0_r8*mu11&
              -1588492186811616081599968857292800000.0_r8*mu10&
              -33613061252672244762688389120000.0_r8*mu13&
              -5488278827131063836192709509120000.0_r8*mu12+8103233168923624359557529600.0_r8*mu15&
              +61795444215044213847613440000.0_r8*mu14-24664649787245911322880000.0_r8*mu17&
              -432357673893294528707328000.0_r8*mu16+482328857446521625612800.0_r8*mu19&
              +5567431011075362582899200.0_r8*mu18-52440257338321366388160.0_r8*mu21&
              -202113807394657288104000.0_r8*mu20+18600933358896627577680.0_r8*mu23&
              +17200859346682290144720.0_r8*mu22+1579029138854919086429.0_r8*mu24)/mu25&
              -316234143225.0_r8/eta25;
       Ck(13)=-0.857680302751143962706154466524e-26_r8&
            *(921771615282027097536169456435200000000.0_r8+&
            8295944537538243877825525107916800000000.0_r8*mu+&
            82944082515127738326629648254894080000000.0_r8*mu3&
            +33875106861614495834454227523993600000000.0_r8*mu2&
            +155831435527215701846558876274524160000000.0_r8*mu5&
            +135579162105263323152398124570378240000000.0_r8*mu4&
            +77894376901480893189435631991783424000000.0_r8*mu7&
            +129172223992777654921877860701437952000000.0_r8*mu6&
            +10563299962843352284952709937024204800000.0_r8*mu9&
            +34015028774686255416039240343486464000000.0_r8*mu8&
            +309518962023963046571482688284262400000.0_r8*mu11&
            +2251774688208263997784763364723916800000.0_r8*mu10&
            +861011777852830761691889811456000000.0_r8*mu13&
            +24446539894904597841934356654981120000.0_r8*mu12&
            -9129386616656671151031877632000.0_r8*mu15+&
            4832123821748417409599243059200000.0_r8*mu14+&
            49062435160051416587387904000.0_r8*mu17&
            -925693351490359114046512128000.0_r8*mu16&
            -502391882666480014007654400.0_r8*mu19+&
             1900142174645277030032947200.0_r8*mu18&
            +14643960795918862872871680.0_r8*mu21&
            -22600849292470119757881600.0_r8*mu20&
            -1065653337094081507074240.0_r8*mu23&
            +1897901415758863141421760.0_r8*mu22&
            -18948349666259029037148.0_r8*mu25&
            -465370750279778090901468.0_r8*mu24+&
             746590869962651602203151.0_r8*mu26)/mu27+7905853580625.0_r8/eta27;
       Ck(14)=0.357366792812976651127564361052e-27_r8*&
             (597308006702753559203437807770009600000000.0_r8+&
             5773977398126617738966565475110092800000000.0_r8*mu+&
             68100855645682274776811213909215150080000000.0_r8*mu3&
             +25501733508392561680435664181736243200000000.0_r8*mu2&
             +157118695440475244955106851435797544960000000.0_r8*mu5&
             +122615489010818931656262442157667778560000000.0_r8*mu4&
             +102492101805070766811555545149580771328000000.0_r8*mu7&
             +147381798850729828655609528466473484288000000.0_r8*mu6&
             +20074198705205260656813464289424284057600000.0_r8*mu9&
             +52900154033967761200642996952685871104000000.0_r8*mu8&
             +1037579354206173520979632805108514816000000.0_r8*mu11&
             +5482069656737299898195754145601814528000000.0_r8*mu10&
             +9090054319524206359024843372442419200000.0_r8*mu13&
             +127656035779367465106545720696949964800000.0_r8*mu12&
             +1504993369039909110503466653712384000.0_r8*mu15&
             +290807579990936642160269940483686400000.0_r8*mu14&
             -232608060355042466127306104832000.0_r8*mu17&
             -2873646390390576097605459959808000.0_r8*mu16+&
             313971468251415675797646950400.0_r8*mu19+&
             12185415148135719819460213555200.0_r8*mu18&
             -1688212326522663865630586880.0_r8*mu21&
             -100256184153099142808793753600.0_r8*mu20&
             +125445449530031752291484160.0_r8*mu23&
             +2336478524503379830664133120.0_r8*mu22&
             -34416214804124455938688800.0_r8*mu25&
             -135809414381175847224032736.0_r8*mu24&
             +17918180879103638452875624.0_r8*mu27&
             +17008660095123205059092520.0_r8*mu26&
             +1511513601028097903631961.0_r8*mu28)/mu29-213458046676875.0_r8/eta29;
       Ck(15)=0.333675810282891364264765976706e-32_r8*&
             (-1855178938018082279529957487152872816640000000000.0_r8&
             -19170182359520183555142894033913019105280000000000.0_r8*mu&
             -263377005417185632954824316365927201374208000000000.0_r8*mu3&
             -91058366207720871886928746661086840750080000000000.0_r8*mu2&
             -732160866586665259142249339807729486659584000000000.0_r8*mu5&
             -517936614457434551853151195128755984007168000000000.0_r8*mu4&
             -604537475382110527688750126952906347249664000000000.0_r8*mu7&
             -766903479869200356559969602402038143411814400000000.0_r8*mu6&
             -161976670064346560043260604954232609151385600000000.0_r8*mu9&
             -360469626710904632477753339927728006103040000000000.0_r8*mu8&
             -13155709023533392876583076232573327245312000000000.0_r8*mu11&
             -54154397171044399889732980618831948834406400000000.0_r8*mu10&
             -248000061425431023736312765275954130452480000000.0_r8*mu13&
             -2232668823382155657464023915361131663196160000000.0_r8*mu12&
             -469438964312978155350762557287618550169600000.0_r8*mu15&
             -16043899459192119732455445262404371742720000000.0_r8*mu14&
             +4299925748509099678162895707970764800000.0_r8*mu17&
             -2252586079669098921188201818674307891200000.0_r8*mu16&
             -14691901028515957708085357473259520000.0_r8*mu19&
             +284592241582281885947507728340705280000.0_r8*mu18&
             +98083396226536882153468585933824000.0_r8*mu21&
             -239626506592093925025883584783360000.0_r8*mu20&
             -1832269195688182263071408321280000.0_r8*mu23&
             -305197628455565458976613761280000.0_r8*mu22&
             +87469847823206583112603812048000.0_r8*mu25&
             -7885304466878737834038191088000.0_r8*mu24&
             -9302956832867085142092385336800.0_r8*mu27&
             +9279023933523831138647855244000.0_r8*mu26&
             -161883106670109285478983023100.0_r8*mu29&
             -3999957450974108642084941683900.0_r8*mu28+&
             8849272268392873147705987190261.0_r8*mu30)/mu31+6190283353629375.0_r8/eta31;  
       s=Ck(0);
       b=1.0_r8;  k=0; t=s;
       DO WHILE ((abs(t/s) > 1.0e-20_r8).AND.(k < 15)) 
         b=-b/a; k=k+1; t=Ck(k)*b; s=s+t  
       ENDDO
       t=s
      ENDIF
      Taeta=t;       
      END FUNCTION Taeta

      SUBROUTINE gseries(a,x,igam,ierr)
      USE Someconstants
      USE GammaError
      IMPLICIT NONE
      REAL(r8) :: a, x, igam
      REAL(r8) :: eps, p, q, t, v
      REAL(r8) :: logamma
      INTEGER ::  ierr, k, m
      eps=epss            
      t=1.0_r8/a;
      v=t;
      k=0
      m=0
      DO WHILE ((abs(t/v)>eps).AND.(m==0))
        p=(a+k)/(a+k+1);
        q=k+1;
        t=-x*t*p/q;
        v=v+t
        k=k+1
        IF (t>giant) m=1
      ENDDO
      IF (m==0) THEN
        IF (a>0.0_r8) THEN
          logamma=loggam(a)
          IF (logamma<loggiant) THEN
            igam=v/gamma(a)
          ELSE
            igam=0.0_r8
            ierr=1
          ENDIF
        ELSE
          IF (1-a<170) THEN
            igam=v/gamma(a)
          ELSE
            igam=0.0_r8
            ierr=1
          ENDIF
        ENDIF
      ELSE
        igam=0.0_r8
        ierr=1
      ENDIF
      END SUBROUTINE gseries

      SUBROUTINE gsereps(a,epsil,x,igam,ierr)
      USE Someconstants
      USE GammaError
      IMPLICIT NONE
      REAL(r8) :: a, x, epsil, igam
      REAL(r8) :: eps, p, q, t, v, piepsi
      REAL(r8) :: frac, piep2, piep4, piep6,&
                  piep8, piep10
      REAL(r8) :: piep3, piep5, piep7,&
                  piep9, piep11, reffor
      INTEGER :: k, l, m, ninta, ierr
      ierr=0
      eps=epss           
      ninta=abs(nint(a)) 
      piepsi=pi*epsil
      IF (ninta>0) THEN
        t=1.0_r8/a;
        v=t;
        k=0
        l=0
        m=0
        DO WHILE ((abs(t/v)>eps).AND.(m==0))
          IF (k+1<ninta) THEN
            p=(a+k)/(a+k+1);
            q=k+1;
            t=-x*t*p/q;
            v=v+t
          ELSEIF (k+1>ninta) THEN
            IF (l==0) THEN
              p=(a+k-1)/(a+k+1);
              q=(k+1)*k;
              t=x*x*t*p/q;
              v=v+t
              l=1
            ELSE
              p=(a+k)/(a+k+1);
              q=k+1;
              t=-x*t*p/q;
              v=v+t
            ENDIF 
          ENDIF
          IF (t>giant) m=1
          k=k+1
        ENDDO
      ELSE
        t=-x/(a+1.0_r8);
        v=t;
        k=1
        m=0
        DO WHILE ((abs(t/v)>eps).AND.(m==0))
          p=(a+k)/(a+k+1);
          q=k+1;
          t=-x*t*p/q;
          v=v+t
          k=k+1
           IF (t>giant) m=1
        ENDDO
      ENDIF
      IF (m==0) THEN
        IF (abs(piepsi)>0.1_r8) THEN
          frac=sin(piepsi)/piepsi
          IF (abs(gamma(a))>dwarf) THEN 
            reffor=1.0_r8/gamma(a)
          ELSE
            m=1
          ENDIF
        ELSE
          IF ((1.0_r8-a)<170) THEN 
            piep2=piepsi*piepsi
            piep4=piep2*piep2
            piep6=piep4*piep2
            piep8=piep6*piep2
            piep10=piep8*piep2
            frac=1.0_r8-piep2/6.0_r8+piep4/120.0_r8-piep6/5040.0_r8+&
                 piep8/(362880.0_r8)-piep10/39916800.0_r8
            piep3=piep2*piepsi
            piep5=piep4*piepsi
            piep7=piep6*piepsi
            piep9=piep8*piepsi
            piep11=piep10*piepsi
            reffor=gamma(1.0_r8-a)/pi*(piepsi-piep3/6.0_r8+piep5/120.0_r8-piep7/5040.0_r8+&
              piep9/(362880.0_r8)-piep11/39916800.0_r8)
            IF (MOD(ninta,2)==0) reffor=-reffor 
          ELSE
            m=1
          ENDIF
        ENDIF
        IF (m==0) THEN
          igam=v*reffor+(x**ninta*quotgamm(1.0_r8-a,1.0_r8+ninta)*frac)
        ENDIF
      ENDIF
      IF (m==1) THEN
        igam=0.0_r8
        ierr=1
      ENDIF
      END SUBROUTINE gsereps

      SUBROUTINE gexpan(a,x,igam,ierr)
      USE Someconstants
      USE GammaError
      IMPLICIT NONE
      REAL(r8) :: a, x, igam
      REAL(r8) :: eps, p, t, v
      INTEGER ::  ierr, k, m
      eps=epss            
      ierr=0
      t=1.0_r8;
      v=t;
      k=0
      m=0
      DO WHILE ((abs(t/v)>eps).AND.(m==0))
        p=(1.0_r8-a+k);
        t=t*p/x;
        v=v+t
        k=k+1
        IF (t>giant) m=1
      ENDDO
      IF (m==0) THEN
        IF ((x>loggiant).OR.(loggam(a)>loggiant)) THEN
          m=1
        ELSE 
          igam=v*exp(x)/(x*gamma(a));
        ENDIF
      ENDIF
      IF (m==1) THEN
        igam=0.0_r8
        ierr=1
      ENDIF
      END SUBROUTINE gexpan

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
   END MODULE IncgamNEG


 