! Program test_bes examines effects of modification of Module mod_bes.
! Module mod_bes_old is the old module, which is renamed from mod_bes
! temporarily for this test.
! Module mod_bes is the modified module.
! Module mod_bes_old is stored in the file pack_bes_old.f95.
! Module mod_bes is stored in the file pack_bes.f95.
! Subroutine hankel2_old is the old subroutine in Module mod_bes_old.
! Subroutine hankel2 is the modified subroutine in Module mod_bes.
! zhan2_old is the Hankel function calculated by Subroutine hankel2_old.
! zhan2 is the Hankel function calculated by Subroutine hankel2.
! For implementing test_bes, connect files pack_bes_old.f95, pack_bes.f95 and
! test_bes.f95 in this order, and make them one file. Compile the connected
! file with a Fortran 95 compiler, and link the object programs, and execute
! the generated executable program.
! We know from the file test_bes.out that the values of zhan2 and info were
! improved.

PROGRAM test_bes
USE mod_bes, ONLY: kp,hankel2,num_region
USE mod_bes_old, ONLY: hankel2_old=>hankel2
IMPLICIT NONE
COMPLEX(kp):: znu,znu0,zhan2,zhan2_old
REAL(kp):: xx,xx0,dd
INTEGER:: i,info,info_old,j
OPEN(7,file='test_bes.out')
znu0=(39.60473566960330_kp,-17.33456476536629_kp)
xx0=30
WRITE(7,7) znu0,xx0
7 FORMAT('znu0=(',F20.14,',',F18.14,'),   xx0=',,F6.3)
dd=5E-9
WRITE(7,6) dd
6 FORMAT('dd=',ES8.1)
WRITE(7,*) 'znu=znu0+dd*j,  xx=xx0+dd*i'
WRITE(7,*)
WRITE(7,*) '   i    j          zhan2_old       info_old ', &
  '           zhan2        info'
DO i=-2,2
  DO j=-3,3
    znu=znu0+dd*i;  xx=xx0+dd*j
    CALL hankel2_old(znu,xx,zhan2_old,info_old)
    CALL hankel2(znu,xx,zhan2,info)
    WRITE(7,8) i,j,zhan2_old,info_old,zhan2,info
  ENDDO
  WRITE(7,8)
ENDDO
8 format(2I5,2(ES14.3,ES12.3,I5))
END PROGRAM test_bes

