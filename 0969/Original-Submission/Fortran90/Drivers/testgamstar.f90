   PROGRAM test
! -----------------------------------------------
! This is a test program for the Fortran 90 
! module incgamNEG. This module includes the
! public routine incgamstar for the computation
! of the gamma star function, defined as 
! gamma*(a,z)=(z**(-a)/GAMMA(a))*gamma(a,z),
! where gamma(a,z) is the lower incomplete
! gamma funcion.
!
! The gamma star function is a single-valued 
! entire function of both a and z and is
! real for positive and negative values
! of a and z.
!------------------------------------------------   
   USE IncgamNEG
   USE Someconstants
   IMPLICIT NONE
   REAL(r16) :: aa(10), a
   REAL(r8)  :: zz1(10), zz2(10), igams(40), igam, errr, z, ar
   INTEGER   :: i, ierr
   DATA aa/0.03_r16,2.0_r16,8.3_r16,25.8_r16,35.9_r16,47.1_r16,&
        67.0_r16,85.5_r16,100.0_r16,130.0_r16/ 
   DATA zz1/0.03_r8,2.0_r8,8.5_r8,15.3_r8,23.0_r8,47.0_r8,&
        55.0_r8,85.0_r8,90.3_r8,100.0_r8/
   DATA zz2/310.0_r8,300.0_r8,250.0_r8,210.0_r8,&
        195.0_r8,180.0_r8,150.0_r8,135.0_r8,120.0_r8,110.0_r8/
   DATA igams/1.0176203871469862_r8,2.0972640247326617_r8,&     
        3.2497009683534236E-002_r8,1.3338766060778538E-020_r8,&
        2.3028505609784507E-032_r8,3.4134982653835186E-040_r8,&
        1.1629710903946985E-071_r8,1.5853438781338887E-093_r8,&
        9.2989354551562908E-120_r8,2.3539206797983229E-177_r8,&
        0.98117651308252274_r8,4.0_r8,31982809.488827940_r8,&     
        -1.7098511726589705E+031_r8,-1.6218042956823193E+049_r8,&
        -5.5089763097616628E+078_r8,-4.0206863555398567E+116_r8,&
        -6.6562356422802371E+163_r8,3.7048875538708741E+195_r8,&
        9.9999999999999993E+259_r8,4.2229702515153779E+130_r8,&
        6.4531721353015019E+127_r8,1.5686255996358435E+102_r8,&
        8.3420489514260343E+063_r8,2.9263521794594825E+042_r8,&
        8.1460946436242266E+017_r8,1.1836585792928070E-030_r8,&
        6.3547792434104463E-072_r8,6.3678090276855394E-107_r8,&
        4.9689181373288279E-173_r8,-4.0799264281243891E+130_r8,&
        90000.0_r8,-3.0882590201738765E+110_r8,3.4090169825076780E+114_r8,&
        7.8766316802761328E+122_r8,4.2384404807027509E+134_r8,&
        -6.2822375846822284E+145_r8,7.6261599326043147E+185_r8,&
        8.2817974522014552E+207_r8,2.4046344822913874E+265_r8/
   open(unit=1,file='testGammStar.dat',status='unknown') 
   write(1,*)'****************************************************************'
   write(1,*)'                  TEST                        '
   write(1,*)'Test of the computed values of the gamma star function   '
   write(1,*)'defined as  gamma*(a,z)=(z**(-a)/GAMMA(a))*gamma(a,z),'
   write(1,*)'where gamma(a,z) is the lower incomplete gamma funcion:'
   write(1,*)'comparison against precomputed values. ERR means relative error'
   write(1,*)'****************************************************************'
   write(1,30)'a','z','gamma*(a,z)','ERR(gamma*(a,z))'
   DO i=1,10
     a=aa(i)
     z=-zz1(i)
     ar=a 
     CALL incgamstar(a,z,igam,ierr)
     errr=abs(1.0_r8-igam/igams(i))
     write(1,31)ar,z,igam,errr
   ENDDO
   DO i=1,10
     a=-aa(i)
     z=-zz1(i)
     ar=a 
     CALL incgamstar(a,z,igam,ierr)
     errr=abs(1.0_r8-igam/igams(i+10))
     write(1,31)ar,z,igam,errr
   ENDDO
   DO i=1,10
     a=aa(i)
     z=-zz2(i)
     ar=a 
     CALL incgamstar(a,z,igam,ierr)
     errr=abs(1.0_r8-igam/igams(i+20))
     write(1,31)ar,z,igam,errr
   ENDDO
   DO i=1,10
    a=-aa(i)
    z=-zz2(i)
    ar=a 
    CALL incgamstar(a,z,igam,ierr)
    errr=abs(1.0_r8-igam/igams(i+30))
    write(1,31)ar,z,igam,errr
  ENDDO
 30 format(2x,a2,8x,a2,8x,a15,6x,a17)
 31 format(d10.4,1x,d10.4,1x,d20.13,1x,d16.9)
   END PROGRAM test
