! File examination.f90
!
! This file includes PROGRAM examination, which is a sample program for 
! execution of MODULE mod_zbes in file mod_zbes.f90.
!
! The symbols znu, zz and so forth are explained in the comment sentences in
! the beginning of File mod_zbes.f90.
! The reference below is also listed there.
!
! The method of the implementation of PROGRAM examination depends on the
! compiler that the user uses. The general method is as follows.
! Compile files examination.f90 and mod_zbes.f90, which are included in archive
! XXX.zip. Link these object programs, and execute the executable program.
! Here we must know that the the program is written in Fortran 90.
! The kind type parameter kp can be selected from the three statements in 
! the beginning of MODULE mod_zbes of file mod_zbes.f90.
! The default value of kp is KIND(1D0).
! When kp=KIND(1.), the output file of PROGRAM examination is 
! examination_04.out. 
! When kp=KIND(1D0), the output file is examination_08.out.
! When kp=SELECTED_REAL_KIND(20,550), the output file is examination_16.out.
!
! PROGRAM examination includes the following tests.
! Tests by Wronskians:
! zbessel1(znu,zz)*zbessel2(znu+1,zz)
!     -zbessel1(znu+1,zz)*zbessel12(znu,zz)=-2/(pi*zz),      (1)
!
! zhankel1(znu,zz)*zhankel2(znu+1,zz)
!   -zhankel1(znu+1,zz)*zhankel2(znu,zz)=(0,4)/(pi*zz).      (2)
!
! error(1) for Eq. (1) is given by
! error(1)=ABS(zbessel1(znu,zz)*zbessel2(znu+1,zz)
!       -zbessel1(znu+1,zz)*zbessel2(znu,zz)+2/(pi*zz))
!       /MAX(ABS(zbessel1(znu,zz)*zbessel2(znu+1,zz)),
!        ABS(zbessel1(znu+1,zz)*zbessel2(znu,zz)),ABS(2/(pi*zz))).
! Here, zbessel1(znu,zz), zbessel2(znu+1,zz) and so forth designate the 
! numerical values calculated by the present algorithm. 
! error(2) for Eq. (2) is also defined similarly.
! error(2)=ABS(zhankel1(znu,zz)*zhankel2(znu+1,zz)
!       -zhankel1(znu+1,zz)*zhankel2(znu,zz)-(0,4)/(pi*zz))
!       /MAX(ABS(zhankel1(znu,zz)*zhankel2(znu+1,zz)),
!        ABS(zhankel1(znu+1,zz)*zhankel2(znu,zz)),ABS(4/(pi*zz))).
!
! error(1) is equal to ABS(eb), where eb is defined in Eq. (42a) of Ref. (1).
! error(2) is equal to ABS(eh), where eh is defined in Eq. (42b) of Ref. (1).
! The test points are put at regular intervals in the region -60<re_znu<60,
! -60<ai_znu<60, -60<re_zz<60 and -60<ai_zz<60.
!
! Explanation of the output files examination_04.out, examination_08.out,
! and examination_16.out.
! Column Equation indicates equation (1) or (2).
! We attempt to calculate error(1) and error(2) at all the test points.
! Column num_t shows numbers of test points.
! Column num_s shows numbers of successful test points. A successful test point
! means that all the cylindrical functions included in equation (1) or (2)  
! were calculated with info=0 and Eqs. (46), (47) are satisfied at this test
! point.
! Column err_max indicates the maximum values of error(1) or error(2).
! Columns znu and zz indicate the values of znu and zz at the test point
! where error(i)=err_max(i) (i=1,2).
! The successful test points in file examination_04.out are very small because
! of overflows and underflows.
!
PROGRAM examination
USE mod_zbes,ONLY: kp,bessel1,bessel2,hankel1,hankel2
IMPLICIT NONE
COMPLEX(kp):: znu_max(4),zz_max(4),znu,zz
REAL(kp):: err_max(4)
REAL(kp):: tiny2,c1,epsilon1
REAL(kp),PARAMETER:: tiny1=TINY(1._kp)
REAL(kp),PARAMETER:: huge1=HUGE(1._kp)/3
INTEGER:: nscf(4)
INTEGER:: i,info,j,k,l,ntry,ig
CHARACTER(LEN=2):: char
WRITE(char,8) kp;  8 FORMAT(I2.2)
OPEN(7,FILE='examination_'//char//'.out')
c1=huge1**0.487     !  =1E150_kp
epsilon1=MAX(EPSILON(1._kp),1E-20_kp)
tiny2=tiny1/SQRT(epsilon1)
WRITE(7,*) 'The output from PROGRAM examination'
WRITE(7,*)
WRITE(7,*) 'Main processor-dependent constants of the processor used ', &
           'here are as follows.'
WRITE(7,*) 'RADIX(1._kp)=',RADIX(1._kp)
WRITE(7,*) 'EPSILON(1._kp)=',EPSILON(1._kp)
WRITE(7,*) 'HUGE(1._kp)=',HUGE(1._kp)
WRITE(7,*) 'TINY(1._kp)=',TINY(1._kp)
WRITE(7,*) 'HUGE(1)=',HUGE(1)
WRITE(7,*)
!
WRITE(7,*) 'Wronskian tests in the region -60<re_znu<60, -60<ai_znu<60,'
WRITE(7,*) '-60<re_zz<60 and -60<ai_zz<60.'
WRITE(7,*)
err_max=0;  nscf=0;  ntry=0;  znu_max=0;  zz_max=0;  ig=3
DO i=-ig,ig
  DO j=-ig,ig
    znu=CMPLX(i*18.3,j*17.8,kp)
    DO k=-ig,ig
      DO l=-ig,ig
        zz=CMPLX(k*18.6+1.1,l*17.4+0.21,kp)
        ntry=ntry+1
        CALL sub_wronski(znu,zz,info)
      ENDDO
    ENDDO
  ENDDO
ENDDO
WRITE(7,*) 'Equation  num_t   num_s    err_max         znu    ',  &
  '            zz'
WRITE(7,26) '(1)   ',ntry,nscf(1),err_max(1),znu_max(1),zz_max(1)
WRITE(7,26) '(2)   ',ntry,nscf(2),err_max(2),znu_max(2),zz_max(2)
26 FORMAT(A8,2I8,ES11.2,2('  (',F7.2,',',F7.2,')'),I5)
CONTAINS

SUBROUTINE sub_wronski(znu,zz,info)
! Tests by the use of Wronskians.
COMPLEX(kp),INTENT(IN):: znu,zz
INTEGER,INTENT(OUT):: info
COMPLEX(kp):: z0,z1,z2,z3,za0,za1,zb0,zb1
REAL(kp):: b1,b2,error
REAL(kp),PARAMETER:: pi=3.14159265358979323846_kp  ! the circle ratio
INTEGER:: info1,info2,info3,info4
CALL bessel1(znu,  zz,za0,info1)
CALL bessel1(znu+1,zz,za1,info2)
CALL bessel2(znu,  zz,zb0,info3)
CALL bessel2(znu+1,zz,zb1,info4)
info=MAX(info1,info2,info3,info4)
IF(info > 21)  WRITE(*,1) 'b12',znu,zz,info
1 FORMAT(A,4ES11.3,2I5)
b1=MIN(abs1(za0),abs1(zb1),abs1(zb0),abs1(za1))
IF(info<2 .AND. b1>tiny2) THEN
  nscf(1)=nscf(1)+1
  IF(abs1(za0) > abs1(zb1)) THEN
    z0=za0;  za0=zb1;  zb1=z0
  ENDIF
  z3=-2/(pi*zz)
  IF(abs1(zb1)>c1 .AND. abs1(za0)>1) THEN
    IF(abs1(zb0) > abs1(za1)) THEN
      z0=zb0;  zb0=za1;  za1=z0
    ENDIF
    z1=za0;  z2=zb0*(za1/zb1);  z3=z3/zb1
  ELSE
    z1=za0*zb1;  z2=zb0*za1
  ENDIF
  b2=MAX(abs1(z1),abs1(z2),abs1(z3))
  error=abs1(z1-z2-z3)/b2
  IF(error > err_max(1)) THEN
    znu_max(1)=znu;  zz_max(1)=zz;  err_max(1)=error
  ENDIF
ENDIF
CALL hankel1(znu,  zz,za0,info1)
CALL hankel1(znu+1,zz,za1,info2)
CALL hankel2(znu,  zz,zb0,info3)
CALL hankel2(znu+1,zz,zb1,info4)
info=MAX(info1,info2,info3,info4)
IF(info > 21)  WRITE(*,1) 'h12',znu,zz,info
b1=MIN(abs1(za0),abs1(zb1),abs1(zb0),abs1(za1))
IF(info<2 .AND. b1>tiny2) THEN
  nscf(2)=nscf(2)+1
  IF(abs1(za0) > abs1(zb1)) THEN
    z0=za0;  za0=zb1;  zb1=z0
  ENDIF
  z3=(0,4)/(pi*zz)
  IF(abs1(zb1)>c1 .AND. abs1(za0)>1) THEN
    IF(abs1(zb0) > abs1(za1)) THEN
      z0=zb0;  zb0=za1;  za1=z0
    ENDIF
    z1=za0;  z2=zb0*(za1/zb1);  z3=z3/zb1
  ELSE
    z1=za0*zb1;  z2=zb0*za1
  ENDIF
  b2=MAX(abs1(z1),abs1(z2),abs1(z3))
  error=abs1(z1-z2-z3)/b2
  IF(error > err_max(2)) THEN
    znu_max(2)=znu;  zz_max(2)=zz;  err_max(2)=error
  ENDIF
ENDIF
END SUBROUTINE sub_wronski

FUNCTION abs1(za) RESULT(ab)
COMPLEX(kp),INTENT(IN):: za
REAL(kp):: ab
ab=ABS(REAL(za))+ABS(AIMAG(za))
END FUNCTION abs1

END PROGRAM examination
