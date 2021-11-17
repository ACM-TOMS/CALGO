  PROGRAM test
! --------------------------------------------------
! This is a test program for the Fortran 90 
! module Parabolic.
! --------------------------------------------------
! 1)First, we compare the values obtained from 
!   the routine parab with 17(non-scaled)+17(scaled)
!   pre-computed values of the functions U(a,x), 
!   U'(a,x), V(a,x), V'(a,x).
!   Relative errors are shown.
! 2)Second, we test the Wronskian relation:
!     W[U(a,x),V(a,x)]=(2/pi)**(1/2)
!   for several values of (a,x).
! -------------------------------------------------- 
  USE Parabolic
  USE Someconstants
  IMPLICIT NONE
  REAL(r8) a,x,wmax,wuvax
  REAL(r8) x1(17),a1(17),x2(17),a2(17)
  REAL(r8) u(17),up(17),v(17),vp(17)
  REAL(r8) us(17),ups(17),vs(17),vps(17)
  REAL(r8) err1,err2,err3,err4
  REAL(r8)uaxx(2),vaxx(2)
  INTEGER i,j,ierr1,mode
  DATA x1/1.0_r8,2.0_r8,3.0_r8,4.0_r8,5.0_r8,6.0_r8,&
          7.0_r8,8.0_r8,9.0_r8,10.0_r8,11.0_r8,13.0_r8,&
          15.0_r8,17.0_r8,19.0_r8,21.0_r8,23.0_r8/
  DATA a1/-90.0_r8,-60.0_r8,-45.0_r8,-30.0_r8,-21.0_r8,&
          -17.0_r8,-11.0_r8,-7.0_r8,-3.0_r8,-1.0_r8,&
          1.0_r8,3.0_r8,7.0_r8,11.0_r8,17.0_r8,21.0_r8,&
          30.0_r8/
  DATA x2/10.0_r8,30.0_r8,50.0_r8,70.0_r8,90.0_r8,110.0_r8,&
          130.0_r8,150.0_r8,170.0_r8,190.0_r8,210.0_r8,230.0_r8,&
          300.0_r8,350.0_r8,400.0_r8,450.0_r8,500.0_r8/
  DATA a2/-1000.0_r8,-600.0_r8,-450.0_r8,-300.0_r8,-210.0_r8,&
          -170.0_r8,-110.0_r8,-70.0_r8,-30.0_r8,-10.0_r8,&
          10.0_r8,30.0_r8,70.0_r8,110.0_r8,170.0_r8,210.0_r8,&
          300.0_r8/
  DATA u/0.763343158893198e+68_r8,-0.912326832009841e+40_r8,&
         0.139768190422770e+28_r8,0.268775889888618e+16_r8,&
        -0.752051240774644e+09_r8,-0.447103531084652e+07_r8,&
         0.981861674474585e+03_r8,0.615550805644269e-01_r8,&
         0.381034094970342e-06_r8,0.439719298311775e-10_r8,&
         0.196758064761047e-14_r8,0.540241715864902e-22_r8,&
         0.490941122431228e-33_r8,0.233650929214398e-45_r8,&
         0.174348672418738e-61_r8,0.291444218947009e-76_r8,&
         0.456223518736535e-99_r8/     
  DATA v/0.547984047326091e-69_r8,-0.494654818446871e-41_r8,&
         0.261630582278972e-28_r8,0.117561853536393e-16_r8,&
         0.118868795045892e-09_r8,0.235964703545892e-07_r8,&
         0.329557754211340e-03_r8,0.216653346048654e+01_r8,&
         0.252191000884935e+06_r8,0.185227194475932e+10_r8,&
         0.362737153165276e+14_r8,0.109781802370711e+22_r8,&
         0.102178304531530e+33_r8,0.187134996691091e+45_r8,&
         0.220950696718772e+61_r8,0.119483159450560e+76_r8,&
         0.686500903792830e+98_r8/ 
  DATA up/-0.811788733386297e+69_r8,0.406928885132567e+41_r8,&
        -0.302983767150835e+28_r8,-0.283970202966543e+16_r8,&
        -0.509425621238514e+10_r8,-0.614513678183800e+07_r8,&
        -0.148339523088500e+04_r8,-0.191087078708465e+00_r8,&
        -0.160680102637161e-05_r8,-0.217671833842493e-09_r8,&
        -0.110847133757808e-13_r8,-0.365336028683281e-21_r8,&
        -0.391892823358278e-32_r8,-0.213780646480184e-44_r8,&
        -0.180943244177870e-60_r8,-0.334472861998361e-75_r8,&
        -0.581932602765032e-98_r8/
  DATA vp/0.462488306857959e-68_r8,-0.653927086297065e-40_r8,&
         0.514147560477428e-27_r8,0.284437900181384e-15_r8,&
        -0.255749146804786e-09_r8,-0.146024573172453e-06_r8,&
         0.314728813577488e-03_r8,0.623649595442493e+01_r8,&
         0.103052143344241e+07_r8,0.897611570216431e+10_r8,&
         0.201161169136914e+15_r8,0.734508410128058e+22_r8,&
         0.809577972092385e+33_r8,0.170265171392036e+46_r8,&
         0.228329356687162e+62_r8,0.136645639848812e+77_r8,&
         0.873228333693647e+99_r8/  
  DATA us/0.998854680728957e-02_r8,-0.814402524588766e-01_r8,&
         0.194204060401539e+00_r8,0.128202906588794e+00_r8,&
         0.108328305270713e+00_r8,0.967316794277160e-01_r8,&
         0.882841004029830e-01_r8,0.819042813146982e-01_r8,&
         0.767753191479941e-01_r8,0.725669801362307e-01_r8,&
         0.689903312667889e-01_r8,0.659002408200904e-01_r8,&
         0.576899698194329e-01_r8,0.534041953518731e-01_r8,&
         0.499468993867898e-01_r8,0.470916053977231e-01_r8,&
         0.446677878430093e-01_r8/
  DATA ups/-0.789619549439823e+01_r8,-0.602078021535948e+01_r8,&
        -0.257597314338556e+01_r8,-0.390035125536867e+01_r8,&
        -0.461576224446812e+01_r8,-0.516905224966607e+01_r8,&
        -0.566362053329541e+01_r8,-0.610475587772802e+01_r8,&
        -0.651256618771261e+01_r8,-0.689023390923585e+01_r8,&
        -0.724743340053157e+01_r8,-0.758726144185510e+01_r8,&
        -0.866704186262370e+01_r8,-0.936257933549912e+01_r8,&
        -0.100106469720780e+02_r8,-0.106176166616135e+02_r8,&
        -0.111937600037488e+02_r8/
  DATA vs/0.100889139898120e+00_r8,0.124027541224908e+00_r8,&
         0.155288668216602e+00_r8,0.102315493005571e+00_r8,&
         0.864429972742099e-01_r8,0.771859559162881e-01_r8,&
         0.704437800772827e-01_r8,0.653523835119920e-01_r8,&
         0.612594420632175e-01_r8,0.579012781107862e-01_r8,&
         0.550472549893944e-01_r8,0.525815275132513e-01_r8,&
         0.460303179344434e-01_r8,0.426106423571517e-01_r8,&
         0.398520454364344e-01_r8,0.375738031430389e-01_r8,&
         0.356398443968751e-01_r8/ 
  DATA vps/0.124561553309290e+00_r8,-0.627969497787122e+00_r8,&
         0.204869620739215e+01_r8,0.311083586032971e+01_r8,&
         0.368217924843460e+01_r8,0.412384364761534e+01_r8,&
         0.451856811921189e+01_r8,0.487061490597449e+01_r8,&
         0.519604990472043e+01_r8,0.549742059304264e+01_r8,&
         0.578245152704741e+01_r8,0.605361618687910e+01_r8,&
         0.691520327647638e+01_r8,0.747018166960289e+01_r8,&
         0.798727864434808e+01_r8,0.847158042589134e+01_r8,&
         0.893128393676945e+01_r8/
  open(unit=1,file='Res',status='unknown') 
! ------------------------------------------------
! We check the values of the unscaled functions
! ------------------------------------------------
  write(1,*)'ERR means relative error'
  write(1,*)'********************************************'
  write(1,*)'Test of the values of the unscaled functions'
  write(1,*)'********************************************'
  write(1,30)'x','a','ERR(Uax)','ERR(Uaxp)','ERR(Vax)','ERR(Vaxp)'
  DO i=1,17
    x=x1(i)
    a=a1(i) 
    ierr1=0
    mode=0
    CALL parab(a,x,mode,uaxx,vaxx,ierr1)
    err1=abs(1.0_r8-uaxx(1)/u(i))
    err2=abs(1.0_r8-uaxx(2)/up(i))
    err3=abs(1.0_r8-vaxx(1)/v(i))
    err4=abs(1.0_r8-vaxx(2)/vp(i))
    write(1,31)int(x),int(a),err1,err2,err3,err4
  ENDDO
! ------------------------------------------------
! We check the values of the scaled functions
! ------------------------------------------------
  write(1,*)'********************************************'
  write(1,*)'Test of the values of the scaled functions'
  write(1,*)'********************************************'
  write(1,30)'x','a','ERR(Uax)','ERR(Uaxp)','ERR(Vax)','ERR(Vaxp)'
  DO i=1,17
    x=x2(i)
    a=a2(i) 
    ierr1=0
    mode=1
    CALL parab(a,x,mode,uaxx,vaxx,ierr1)
    err1=abs(1.0_r8-uaxx(1)/us(i))
    err2=abs(1.0_r8-uaxx(2)/ups(i))
    err3=abs(1.0_r8-vaxx(1)/vs(i))
    err4=abs(1.0_r8-vaxx(2)/vps(i))
    write(1,31)int(x),int(a),err1,err2,err3,err4
  ENDDO 
 ! ---------------------------------------------------------
 !  Test of the Wronskian relation:  W[Uax,Vax]=(2/pi)**(1/2)
 ! ---------------------------------------------------------
 ! We test the relation for 0<x<=1000,  -1000<a<=1000        
 ! ---------------------------------------------------------
 ! This is a Wronskian test for the scaled functions       
  mode=1
 ! ---------------------------------------------------------
  write(1,*)' '
  write(1,*)'--------------------------------------------------'
  write(1,*)'Test of the Wronskian relation for U(a,x), V(a,x)' 
  write(1,*)'--------------------------------------------------'
  wmax=0.0_r8
  DO i=1,100,1
    DO j=-100,100,20
      x=i*0.95_r8
      a=j*0.9_r8
      ierr1=0
      CALL parab(a,x,mode,uaxx,vaxx,ierr1)    
      IF (ierr1==0) THEN
        wuvax=abs((uaxx(1)*vaxx(2)-uaxx(2)*vaxx(1))/sqrt2opi-1.0_r8);   
        IF (wuvax>wmax) THEN
          wmax=wuvax
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  write(1,*)' '
  write(1,*)'Maximum value of the Wronskian check =',wmax
 30 format(3x,a1,3x,a1,8x,a8,8x,a9,8x,a8,8x,a9)
 31 format(i4,1x,i5,d16.9,1x,d16.9,1x,d16.9,1x,d16.9)    
  END PROGRAM test

     
