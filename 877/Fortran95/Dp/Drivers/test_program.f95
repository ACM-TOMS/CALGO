! File test_program.f95
!
! The file test_program.f95 contains PROGRAM test_program, which tests the 
! present algorithm. This program contains tests 1-9.
!
! The notations zbessel(znu,xx), zneumann(znu,xx), zhankel1(znu,xx),
! zhankel2(znu,xx), znu, xx, info, nregion, kp and so forth are explained in
! the comment lines in the beginning of file pack_bes.f95.
! The references are also described there.
!
! The default value of kp is defined in MODULE mod_bes in file pack_bes.f95.
! The default value of kp is kp=KIND(1D0).
! We can redefine kp, and kp can have an arbitrary value if the processor 
! permits this kp. 
! The value of kp in this file is equal to kp defined in MODULE mod_bes.
!
! A method of the implementation of PROGRAM test_program is as follows. Connect
! files pack_bes.f95, test_program.f95 in this order, and make them one file.
! Compile the connected file with a Fortran 95 compiler, and link the object
! programs, and execute the executable program.
!
! The output from PROGRAM test_program will be stored in file 
! test_program_ff.out, where ff denotes double figures expressing kp.
! The tables stated below are stored in file test_program_ff.out.
!
! The test program examines the value of info and sometimes outputs a 
! functional name, znu, xx, info to the standard output unit in order to give 
! us attention. When kp=KIND(1.0), messages for attention issue from 
! PRGRAM test_program frequently, but we do not need to worry about the
! messages.
!
! If info=35, this info informs us that the present algorithm has some defects
! that we should modify.
!
! The output from PROGRAM test_program is influenced not only by the value
! of kp but also the processor. The influence of the processor is written in
! Section 3.2 of Ref. (1).
!
! PROGRAM test_program examines the present algorithm under extreme conditions.
! Hence, depending on the processor, the implementation of the test programs may 
! discontinue because of the reasons stated in Section 3.2 of Ref. (1), but the 
! results computed by the present algorithm is still reliable if info=0.
!
 MODULE mod_test
 USE mod_bes, ONLY: kp,abs2,D1MACH,I1MACH,bes_series,neu_series,bes_han_dby, &
   bes_olver,bes_recur,han2_temme,han2_olver,num_region, &
   sqrt_huge1,hugelog,epsilon1,err_range1,err_range2,err_range3, &
   theta_lim,zunit,huge1,tiny1,epsilon0,pi,ihuge,bessel,neumann,hankel1,hankel2
 COMPLEX(kp):: znu_max(4),zans_max(4)
 REAL(kp):: err_max(4),xx_max(4)
 INTEGER:: nreg_max(4),id(7)
 END MODULE mod_test

 PROGRAM test_program
 USE mod_test,ONLY: bessel,neumann,hankel1,hankel2,err_max,id,nreg_max,xx_max, &
   znu_max,tiny1,hugelog,huge1,zunit,pi, kp,epsilon1,epsilon0,err_range1, &
   err_range2,err_range3,sqrt_huge1,theta_lim,abs2,num_region
 IMPLICIT NONE
 COMPLEX(kp):: zans,zbes,zhan1,zhan2,zneu,znu,znus,zans0,zans3,za,znua,znu0
 REAL(kp):: bb,error(4),xx,xxs,xx0,phi,theta,rr,anu,huge2,tiny2,a1,a2,rp,eps0
 REAL(kp),EXTERNAL:: amp
 INTEGER:: info,info0,info2,info3,i5,i6,nregion
 INTEGER:: i,j,k,nregionm1,nregionm2,nregionp1,nregionp2,nregionp3,nrem, &
          nrep,itotal,itotalm1,itotalm2,itotalp1,itotalp2,itotalp3
 CHARACTER(LEN=2) kp_num
1 FORMAT(A8,ES11.2,'  (',ES9.2,',',ES9.2,')',ES11.2E3,I5)
2 FORMAT(A,I1,A,4ES9.2)
5 FORMAT(A9,' (',ES10.2E3,',',ES10.2E3,')',ES11.2E3,' (',ES10.2E3,',',&
     ES10.2E3,')',I4)
6 FORMAT(A8,ES11.2,'  (',ES9.2,',',ES9.2,')',ES10.2,I5)
8 FORMAT(' za = (0,',i2,')')
 WRITE(kp_num,'(I2.2)') kp
 OPEN(UNIT=7,FILE='test_program_'//kp_num//'.out')
 WRITE(7,*) 'File test_program_'//kp_num//'.out'
 WRITE(7,*) 'File test_program_'//kp_num//'.out is an output ', &
            'file from PROGRAM test_program.'
 WRITE(7,*)
 WRITE(7,*) 'Main processor-dependent constants of the processor used ', &
             'here are'
 WRITE(7,*) 'HUGE(1)=',HUGE(1)
 WRITE(7,*) 'RADIX(1._kp)=',RADIX(1._kp)
 WRITE(7,*) 'EPSILON(1._kp)=',EPSILON(1._kp)
 WRITE(7,*) 'HUGE(1._kp)=',HUGE(1._kp)
 WRITE(7,*) 'TINY(1._kp)=',TINY(1._kp)

! Test 1
! Tests by Wronskians:
! zbessel(znu,xx)*zneumann(znu+1,xx)
!     -zbessel(znu+1,xx)*zneumann(znu,xx)=-2/(pi*xx),      (1)
!
! zhankel1(znu,xx)*zhankel2(znu+1,xx)
!   -zhankel1(znu+1,xx)*zhankel2(znu,xx)=(0,4)/(pi*xx).    (2)
!
! For example, error(1) for Equation (1) is given by
! error(1)=ABS(zbessel(znu,xx)*zneumann(znu+1,xx)
!       -zbessel(znu+1,xx)*zneumann(znu,xx)+2/(pi*xx))
!       /MAX(ABS(zbessel(znu,xx)*zneumann(znu+1,xx)),
!        ABS(zbessel(znu+1,xx)*zneumann(znu,xx)),ABS(2/(pi*xx))).
! Here, zbessel(znu,xx), zneumann(znu+1,xx) and so forth designate the numerical
! values are calculated by the present algorithm. error(2) for equation (2) is
! also defined similarly.
!
! error(i) (i=1,2) are computed at the test points shown in Table 1-1.
! Line (a) of Table 1-1 shows nregion and Total.
! Line (b) shows the number of the test points in each nregion and the total of
! the test points.
!
! In Table 1-2, column Equation indicates the equation number (i), and column
! err_max indicates the maximum among the errors error(i) evaluated at the 
! test points.
! Table 1-2 also indicates znu, xx, nregion when error(i) reaches the maximum.
!
 err_max=0;  id=0
 DO i=0,100
   DO j=0,100
     znu=CMPLX(i*0.643,j*0.63)
     DO k=0,100
       xx=k*0.97+0.024
       CALL sub_wronsky(znu,xx,info)
       IF(info > 2) EXIT
     ENDDO
   ENDDO
 ENDDO
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 1'
 WRITE(7,*) 'Tests by Wronskians in the region 0<=re_znu<=64,'
 WRITE(7,*) '0<=aim_znu<=63,  0.024<=xx<=97.'
 WRITE(7,*)
 WRITE(7,*) '    Table 1-1'
 WRITE(7,'(A,7I8,A)') '(a)',(i,i=1,7),'   Total'
 itotal=0
 DO i=1,7
   itotal=itotal+id(i)
 ENDDO
 WRITE(7,'(A,8I8)') '(b)',id,itotal
 WRITE(7,*)
 WRITE(7,*) '    Table 1-2'
 WRITE(7,*) 'Equation   err_max           znu              xx    nregion'
 WRITE(7,6) '(1)   ',err_max(1),znu_max(1),xx_max(1),nreg_max(1)
 WRITE(7,6) '(2)   ',err_max(2),znu_max(2),xx_max(2),nreg_max(2)


! Test 2
! Tests by the Wronskian equations (i) (i=1,2) of test 1 when xx<<1.
! error(i) are defined in test 1.
! error(i) are computed at the test points shown in Table 2-1.
! Line (a) of Table 2-1 shows nregion and Total.
! Line (b) shows the number of the test points for each nregion and the total of
! the test points.
!
! In Table 2-2, column Equation indicates the equation number (i), and column
! err_max indicates the maximum among the errors error(i) evaluated at the
! test points.
! Table 2-2 also indicates znu, xx, nregion when error(i) reaches the maximum.
!
 err_max=0;  id=0
 DO i=0,20
   DO j=0,20
     znu=CMPLX(i*0.99,j*1.21)
     xx=1E-30_kp
     DO k=1,15
       CALL sub_wronsky(znu,xx,info)
       IF(info > 2) EXIT
       xx=xx*1E-10
     ENDDO
   ENDDO
 ENDDO
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 2'
 WRITE(7,*) 'Tests by Wronskians in the region 0<=re_znu<=20, &
            &0<=aim_znu<=24,'
 WRITE(7,*) '1E-180<=xx<=1E-30 under the condition of info=0.'
 WRITE(7,*)
 WRITE(7,*) '    Table 2-1'
 WRITE(7,'(A,7I7,A)') '(a)',(i,i=1,7),'  Total'
 itotal=0
 DO i=1,7
   itotal=itotal+id(i)
 ENDDO
 WRITE(7,'(A,8I7)') '(b)',id,itotal
 WRITE(7,*)
 WRITE(7,*) '    Table 2-2'
 WRITE(7,*) 'Equation   err_max           znu               xx    nregion'
 WRITE(7,1) '(1) ',err_max(1),znu_max(1),xx_max(1),nreg_max(1)
 WRITE(7,1) '(2) ',err_max(2),znu_max(2),xx_max(2),nreg_max(2)


! Test 3
! Tests for ABS(znu)+xx>23 by the Wronsky equation (i) (i=1,2) defined in
! test 1.
! error(i) are computed at the test points shown in Table 3-1.
! Line (a) shows nregion and Total.
! Line (b) shows the number of the test points for each nregion and the total of
! the test points.
!
! In Table 3-2, column Equation indicates the equation number (i), and column
! err_max indicates the maximum among the errors error(i) evaluated at the
! test points.
! Table 3-2 also indicates znu, xx, nregion when error(i) reaches the maximum.
!
 id=0;   err_max=0;  i5=30
 DO j=0,i5
   theta=(pi*j)/(2*i5)
   i6=i5*SIN(theta)+1
   DO i=0,i6
     phi=(pi*i)/(2*i6)
     rr=20
     DO
       rr=rr*1.17
       xx=rr*COS(theta)+0.03
       znu=rr*CMPLX(COS(phi),SIN(phi))*SIN(theta)
       CALL sub_wronsky(znu,xx,info)
       IF(rr>3000 .OR. info>2) EXIT
     ENDDO
   ENDDO
 ENDDO

 WRITE(7,'(//)')
 WRITE(7,*) 'Test 3'
 WRITE(7,*) 'Tests by Wronskians in the region ABS(znu)+xx>20'
 WRITE(7,*) 'under the condition of info=0.'
 WRITE(7,*)
 WRITE(7,*) '    Table 3-1'
 WRITE(7,'(A,7I7,A)') '(a)',(i,i=1,7),'   Total'
 itotal=0
 DO i=1,7
   itotal=itotal+id(i)
 ENDDO
 WRITE(7,'(A,8I7)') '(b)',id,itotal
 WRITE(7,*)
 WRITE(7,*) '    Table 3-2'
 WRITE(7,*) 'Equation   err_max           znu              xx    nregion'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1),nreg_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2),nreg_max(2)


! Test 4
! Here, we define a function fun_p(zcylind(znu,xx)),
! where zcylind(znu,xx) is a cylindrical function, that is, one of
! zbessel(znu,xx), zneumann(znu,xx), zhankel1(znu,xx), zhankel2(znu,xx).
! The value of fun_p(zcylind(znu,xx)) is given by
! fun_p(zcylind(znu,xx))=amp('cyl',znu,xx,zans),
! where FUNCTION amp is a subprogram of PROGRAM test_program, and is stored
! in this file.  This fun_p is used in tests 4,3,6.
! The function fp in Section 4.1.3 of Ref. (1) corresponds to 
! fun_p(zcylind(znu,xx)).
!
! Tests on Eqs. (8a) of Ref. (1):
! CONJG(zbessel(znu,xx)) -  zbessel(CONJG(znu),xx)=0,       (1)
! CONJG(zneumann(znu,xx))- zneumann(CONJG(znu),xx)=0.       (2)
!
! error(1) of equation (1) is defined by
! error(1)=ABS(CONJG(zbessel(znu,xx))-zbessel(CONJG(znu),xx)) &
!            /fun_p(zbessel(znu,xx))
! error(2) is also defined similarly.
! error(i) (i=1,2) are computed at the test points shown in Table 4-1.
! Line (a) of Table 4-1 shows nregion and Total.
! Line (b) indicates the number of the test points in each nregion and the total
! of the test points.
!
! In Table 4-2, column Equation indicates the equation number (i), and column
! err_max indicates the maximum among the errors error(i) evaluated at the
! test points.
! Table 4-2 also shows znu, xx, nregion when error(i) reaches the maximum.
!
 err_max=0;   znu_max=0
 xx_max=0;    nreg_max=0
 id=0
 DO i=-45,45
   DO j=0,50
     znu=CMPLX(i*1.2,j*1.1)
     DO k=0,50
       xx=k*1.1+0.03
       CALL sub_conjg(znu,xx,info2)
       IF(info2 > 2) EXIT
     ENDDO
   ENDDO
 ENDDO
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 4'
 WRITE(7,*) 'Tests on Eqs. (8a) of Ref. (1) in the region'
 WRITE(7,*) 'ABS(re_znu)<=54,  0<=aim_znu<=55,  0.03<=xx<=55.'
 WRITE(7,*)
 WRITE(7,*) '    Table 4-1'
 WRITE(7,'(A,7I7,A)') '(a)',(i,i=1,7),'  Total'
 itotal=0
 DO i=1,7
   itotal=itotal+id(i)
 ENDDO
 WRITE(7,'(A,8I7)') '(b)',id,itotal
 WRITE(7,*)
 WRITE(7,*) '    Table 4-2'
 WRITE(7,*) 'Equation   err_max           znu              xx    nregion'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1),nreg_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2),nreg_max(2)


! Test 5
! Tests on the first of Eqs. (8b) of Ref. (1).
! Tests on the equations:
! EXP(-pi*znu*zunit)*zhankel1(-znu,xx)-zhankel1(znu,xx)=0,     (1)
!
! error(1) of equation (1) is defined by
! error(1)=ABS(EXP(-pi*znu*zunit)*zhankel1(-znu,xx)-zhankel1(znu,xx))
!           /fun_p(zhankel1(znu,xx))
!
! error(1) are computed at the test points shown in Table 5-1.
! Line (a) of Table 5-1 indicates nregion and Total.
! Line (b) shows the number of the test points in each nregion and the total of
! the test points.
!
! In Table 5-2, column Equation indicates the equation number (1), and
! column err_max indicates the maximum among the errors error(1) evaluated at
! the test points.
! Table 5-2 also indicates znu, xx, nregion when error(1) reaches the maximum.
!
 err_max=0;   id=0
 DO i=0,50
   DO j=-45,45
     znu=CMPLX(i*1.2+0.25,j*1.1+0.1)
     DO k=0,50
       xx=k*1.1+0.18
       CALL hankel1(znu,xx,zans0,info0)
       CALL test_info('zhankel1',znu,xx,info0)
       CALL hankel1(-znu,xx,zans3,info3)
       CALL test_info('zhankel1',-znu,xx,info3)
       info=MAX(info0,info3)
       nregion=num_region(znu,xx)
       IF(info<2 .AND. ABS(REAL(zunit*pi*znu))<hugelog) THEN
         id(nregion)=id(nregion)+1
         IF(REAL(-zunit*pi*znu) > hugelog) THEN
           error(1)=0
         ELSE
           error(1)=abs2(EXP(-zunit*pi*znu)*zans3-zans0)/amp('hn1',znu,xx,zans0)
         ENDIF
         IF(error(1) > err_max(1)) THEN
           znu_max(1)=znu;  xx_max(1)=xx
           err_max(1)=error(1)
           nreg_max(1)=nregion
         ENDIF
       ENDIF
     ENDDO
   ENDDO
 ENDDO

 WRITE(7,'(//)')
 WRITE(7,*) 'Test 5'
 WRITE(7,*) 'Tests on the first of Eqs. (8b) of Ref. (1) in the region'
 WRITE(7,*) '0<=re_znu<=60,  ABS(aim_znu)<=50,  0.18<=xx<=55.'
 WRITE(7,*)
 WRITE(7,*) '    Table 5-1'
 WRITE(7,'(A,7I7,A)') '(a)',(i,i=1,7),'  Total'
 itotal=0
 DO i=1,7
   itotal=itotal+id(i)
 ENDDO
 WRITE(7,'(A,8I7)') '(b)',id,itotal
 WRITE(7,*)
 WRITE(7,*) '    Table 5-2'
 WRITE(7,*) 'Equation   err_max           znu              xx    nregion'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1),nreg_max(1)


! Test 6
! The methods for evaluating a cylindrical function on both sides of
! the boundary of any two regions Rm, Rn (m,n=1,...,7) are different.
! Test 6 examines the relative errors between the two values calculated on both
! the sides of each boundary.
! The following equations should hold on the boundary of regions Rm and Rn
! zbessel_bdry(m,znu,xx) -zbessel_bdry(n,znu,xx) =0,          (1)
! zneumann_bdry(m,znu,xx)-zneumann_bdry(n,znu,xx)=0,          (2)
! zhankel2_bdry(m,znu,xx)-zhankel2_bdry(n,znu,xx)=0,          (3)
! where the point (znu,xx) must exist on the boundary of regions Rm and Rn, and
! zbessel_bdry(m,znu,xx) is the value of the bessel function calculated 
! by the method at nregion=m. The method is described in Section 2.2 of Ref.
! (1). zneumann_bdry(m,znu,xx), zhankel1_bdry(m,znu,xx), zhankel2_bdry(m,znu,xx)
! are also defined similarly.
!
! error(1) of equation (1) is defined by
! error(1)=ABS(zbessel_bdry(m,znu,xx)-zbessel_bdry(n,znu,xx))
!                      /fun_p(zbessel(znu,xx))
! error(2), error(3) are also defined similarly.
! Test 6 examines error(1), (2), (3) under the condition of info=0.
!
! In Tables 6-1,....,6-9, columns Equation indicate the equation number (i),
! where i=1,2,3. 
! Columns err_max indicate the maximums among the errors error(i) calculated at
! the test points.
! Each of tables 6-1,....,6-9 also indicates znu, xx when error(i) reaches the
! maximum.  Each Table also indicates the total numbers of the test points on
! the boundaries.
!
! The statements marked with !* in the program below ascertain that the point 
! (znu,xx) exists on the boundary between two regions Rm, Rn.
!
! SUBROUTINE boundary_test(nregion1,nregion2,znu,xx,itotal) computes the maximum
! value of error(i).
! nregion1= region number m
! nregion2= region number n
! (znu,xx)= a point on the boundary of regions Rm, Rn.
! The maximum value of the errors error(i) is stored in err_max(i).
!
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 6'
 WRITE(7,*)
 eps0=5*EPSILON(1.0)
! Test 6-1-1
! Tests on the boundary between regions R1 and R2.  0=<INT(re_znu)=<32.
 nregionm1=2
 nregionp1=1
 err_max=0;   znu_max=0;   xx_max=0;   itotalp1=0
 DO i=0,32
   znua=0.3
   DO j=1,4
     znu=i+znua
     DO k=1,3
       xx=k*0.3
       nrem=num_region(znu-znua*eps0,xx)      !*
       IF(nrem /= nregionm1) CYCLE            !*
       nrep=num_region(znu+znua*eps0,xx)      !*
       IF(nrep /= nregionp1) CYCLE            !*
       itotalp1=itotalp1+1
       CALL boundary_test(nrem,nrep,znu,xx)
     ENDDO
     znua=znua*(0,1)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-1-1'
 WRITE(7,*) 'Tests on the boundary between regions R1 and R2.  &
             &0=<INT(re_znu)=<32.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R2=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R1=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-1-2
! Tests on the boundary between regions R1 and R2.  30=<INT(re_znu)=<120.
 nregionm1=2
 nregionp1=1
 err_max=0;   znu_max=0;   xx_max=0;   itotalp1=0
 DO i=6, 24
   znua=0.3
   DO j=1,4
     znu=i*5+znua
     DO k=1,3
       xx=k*0.3
       nrem=num_region(znu-znua*eps0*10,xx)    !*
       IF(nrem /= nregionm1) CYCLE             !*
       nrep=num_region(znu+znua*eps0*10,xx)    !*
       IF(nrep /= nregionp1) CYCLE             !*
       itotalp1=itotalp1+1
       CALL boundary_test(nrem,nrep,znu,xx)
     ENDDO
     znua=znua*zunit
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-1-2'
 WRITE(7,*) 'Tests on the boundary between regions R1 and R2.  ',&
             '30=<INT(re_znu)=<120.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R2=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R1=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-2-1
! Tests on the boundary plane xx=1 in the region 0<re_znu<32, 0<aim_znu<32.
 nregionm1=1;    nregionm2=2
 nregionp1=6;    nregionp2=7
 err_max=0;   znu_max=0;   xx_max=0
 itotalm1=0;  itotalm2=0;  itotalp1=0;  itotalp2=0
 xx=1
 DO i=0,40
   DO j=0,37
     znu=CMPLX(i*0.8, j*0.87)
     nrem=num_region(znu,xx*(1-eps0))                      !*
     IF(nrem/=nregionm1 .AND. nrem/=nregionm2) CYCLE       !*
     nrep=num_region(znu,xx*(1+eps0))                      !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2) CYCLE       !*
     IF(nrem == nregionm1) THEN
       itotalm1=itotalm1+1
     ELSE
       itotalm2=itotalm2+1
     ENDIF
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE
       itotalp2=itotalp2+1
     ENDIF
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-2-1'
 WRITE(7,*) 'Tests on the boundary plane xx=1 in the region 0<re_znu<32, &
              &0<aim_znu<32.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R1=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R2=',&
               itotalm2
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R6=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R7=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1+itotalp2
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-2-2
! Tests on the boundary plane xx=1 in the region ABS(znu)>32.
 nregionm1=1;    nregionm2=2
 nregionp1=4
 err_max=0;   znu_max=0;   xx_max=0
 itotalm1=0;  itotalm2=0;  itotalp1=0;  itotalp2=0
 xx=1;  rr=32
 DO i=0,15
   DO j=0,8
     znu=rr*EXP(pi*(zunit/16)*j)
     nrem=num_region(znu,xx*(1-eps0))                        !*
     IF(nrem/=nregionm1 .AND. nrem/=nregionm2) CYCLE         !*
     nrep=num_region(znu,xx*(1+eps0))                        !*
     IF(nrep /= nregionp1) CYCLE                             !*
     IF(nrem == nregionm1) THEN
       itotalm1=itotalm1+1
     ELSE
       itotalm2=itotalm2+1
     ENDIF
     itotalp1=itotalp1+1
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
   rr=rr*1.1
 ENDDO
 WRITE(7,*) 'Table 6-2-2'
 WRITE(7,*) 'Test on the boundary plane xx=1 in the region ABS(znu)>32.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R1=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R2=',&
               itotalm2
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R4=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-3-1
! Tests on the boundary between regions R3,R4 and region R5.  1<xx<32
 nregionm1=5
 nregionp1=3
 nregionp2=4
 err_max=0;   znu_max=0;   xx_max=0
 itotalm1=0;  itotalp1=0;  itotalp2=0
 DO i=1,30
   xx=i*1.1
   IF(xx < 20) THEN
     rp=0.512900114*xx+12
   ELSE
     rp=8.2*xx**0.33333
   ENDIF
   DO j=0,7
     theta=(pi/4)*j
     znu=CMPLX(SIN(theta),COS(theta),kp)*rp+xx
     IF(REAL(znu) < 0) CYCLE                                          !*
     nrem=num_region(znu-eps0*rp*CMPLX(SIN(theta),COS(theta),kp),xx)  !*
     IF(nrem /= nregionm1) CYCLE                                      !*
     nrep=num_region(znu+eps0*rp*CMPLX(SIN(theta),COS(theta),kp),xx)  !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2) CYCLE                  !*
     itotalm1=itotalm1+1
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE
       itotalp2=itotalp2+1
     ENDIF
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-3-1'
 WRITE(7,*) 'Tests on the boundary between regions R3,R4 and region R5.&
            &  1<xx<32.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R5=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R3=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R4=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
                itotalm1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-3-2
! Tests on the boundary between regions R3,R4 and region R5.   36<xx<4100.
 nregionm1=5;  nregionp1=3;  nregionp2=4
 err_max=0;   znu_max=0;   xx_max=0
 itotalm1=0;  itotalp1=0;  itotalp2=0
 xx=30
 DO i=1,27
   xx=xx*1.2
   IF(xx < 20.00) THEN
     rp=0.512900114*xx+12
   ELSE
     rp=8.2*xx**0.33333
   ENDIF
   DO j=0,7
     theta=(pi/4)*j
     znu=CMPLX(SIN(theta),COS(theta),kp)*rp+xx
     nrem=num_region(znu-2*eps0*rp*CMPLX(SIN(theta),COS(theta),kp),xx)  !*
     IF(nrem /= nregionm1) CYCLE                                        !*
     nrep=num_region(znu+2*eps0*rp*CMPLX(SIN(theta),COS(theta),kp),xx)  !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2) CYCLE                    !*
     itotalm1=itotalm1+1
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE
       itotalp2=itotalp2+1
     ENDIF
       CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-3-2'
 WRITE(7,*) 'Tests on the boundary between regions R3,R4 and region R5.  &
            &36<xx<4100.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R5=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R3=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R4=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalm1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-4
! Tests on the boundary plane re_znu=20.
 nregionm1=6;   nregionm2=7
 nregionp1=4;   nregionp2=5
 err_max=0;    znu_max=0;   xx_max=0
 itotalm1=0;   itotalm2=0;  itotalp1=0;  itotalp2=0
 DO i=0,5
   znu=CMPLX(20, i*0.97)
   DO j=1,22
     xx=j*0.96
     nrem=num_region(CMPLX(REAL(znu)*(1-eps0),AIMAG(znu),kp),xx) !*
     IF(nrem/=nregionm1 .AND. nrem/=nregionm2) CYCLE             !*
     nrep=num_region(CMPLX(REAL(znu)*(1+eps0),AIMAG(znu),kp),xx) !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2) CYCLE             !*
     IF(nrem == nregionm1) THEN
       itotalm1=itotalm1+1
     ELSE
       itotalm2=itotalm2+1
     ENDIF
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE
       itotalp2=itotalp2+1
     ENDIF
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-4'
 WRITE(7,*) 'Tests on the boundary plane re_znu=20.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R6=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R7=',&
               itotalm2
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R4=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R5=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1+itotalp2
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-5
! Tests on the boundary plane aim_znu=0.3.
! Errors on the boundary between regions R6 and R7
 nregionm1=7;  nregionp1=6
 err_max=0;   znu_max=0;   xx_max=0;   itotalp1=0
 DO i=0,21
   znu=CMPLX(i*0.97,0.3)
   DO j=1,26
     xx=j*1.1
     nrem=num_region(CMPLX(REAL(znu),AIMAG(znu)*(1-eps0),kp),xx)  !*
     IF(nrem /= nregionm1) CYCLE                                  !*
     nrep=num_region(CMPLX(REAL(znu),AIMAG(znu)*(1+eps0),kp),xx)  !*
     IF(nrep /= nregionp1) CYCLE                                  !*
     itotalp1=itotalp1+1
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-5'
 WRITE(7,*) 'Tests on the boundary plane aim_znu=0.3.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R7=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R6=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalp1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-6
! Tests on the boundary plane aim_znu=20.
 nregionm1=6;  nregionp1=3;  nregionp2=4;  nregionp3=5
 err_max=0;    znu_max=0;    xx_max=0
 itotalm1=0;   itotalp1=0;   itotalp2=0;  itotalp3=0
 DO i=0,22
   znu=CMPLX(i*0.97,20)
   DO j=1,26
     xx=j*1.08
     nrem=num_region(CMPLX(REAL(znu),AIMAG(znu)*(1-eps0),kp),xx)  !*
     IF(nrem /= nregionm1) CYCLE                                  !*
     nrep=num_region(CMPLX(REAL(znu),AIMAG(znu)*(1+eps0),kp),xx)  !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2 .AND. &             !*
            nrep/=nregionp3) CYCLE                                !*
     itotalm1=itotalm1+1
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE IF(nrep == nregionp2) THEN
       itotalp2=itotalp2+1
     ELSE
       itotalp3=itotalp3+1
     ENDIF
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-6'
 WRITE(7,*) 'Tests on the boundary plane aim_znu=20.'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R6=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R3=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R4=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R5=',&
               itotalp3
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalm1
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)
 WRITE(7,*)

! Test 6-7
! Tests on the boundary plane re_znu=20-4*(xx-23).
 nregionm1=6;  nregionm2=7
 nregionp1=3;  nregionp2=5
 err_max=0;    znu_max=0;   xx_max=0
 itotalm1=0;   itotalm2=0;  itotalp1=0;  itotalp2=0;  itotalp3=0
 DO i=0,21
   a1=i*0.97_kp
    xx=28-a1/4    ! xx=(20-a1)/4+23
   DO j=0,21
     znu=CMPLX(a1,j*0.967,kp)
     nrem=num_region(znu,xx*(1-eps0))                  !*
     IF(nrem/=nregionm1 .AND. nrem/=nregionm2) CYCLE   !*
     nrep=num_region(znu,xx*(1+eps0))                  !*
     IF(nrep/=nregionp1 .AND. nrep/=nregionp2) CYCLE   !*
     IF(nrem == nregionm1) THEN
       itotalm1=itotalm1+1
     ELSE
       itotalm2=itotalm2+1
     ENDIF
     IF(nrep == nregionp1) THEN
       itotalp1=itotalp1+1
     ELSE
       itotalp2=itotalp2+1
     ENDIF
     CALL boundary_test(nrem,nrep,znu,xx)
   ENDDO
 ENDDO
 WRITE(7,*) 'Table 6-7'
 WRITE(7,*) 'Tests on the boundary plane re_znu=20-4*(xx-23).'
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R6=',&
               itotalm1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R7=',&
               itotalm2
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R3=',&
               itotalp1
 WRITE(7,'(A,I4)') ' The number of test points on the surface of region R5=',&
               itotalp2
 WRITE(7,'(A,I4)') ' The total number of test points                      =',&
               itotalm1+itotalm2
 WRITE(7,*) 'Equation   err_max           znu              xx'
 WRITE(7,6) '(1) ',err_max(1),znu_max(1),xx_max(1)
 WRITE(7,6) '(2) ',err_max(2),znu_max(2),xx_max(2)
 WRITE(7,6) '(3) ',err_max(3),znu_max(3),xx_max(3)


! Test 7
! Tests for xx>>1, when znu=xx+za*xx**0.2,
! where za are complex constants shown in Table 7.
! Test 7 ascertains that zans*znu**(1./3) converges to a value as xx tends to
! an infinity, and also ascertains that info gives correct values at the limit.
! Table 7 shows zans*znu**(1./3) and info for each test point (znu,xx), za,
! and each function.
!
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 7'
 WRITE(7,*) 'Tests for xx >> 1, when znu = xx + za*xx**0.2,'
 WRITE(7,*) 'where za are complex constants shown in Table 7.'
 WRITE(7,*) 'Test 7 ascertains that zans*znu**(1./3) converges to a value '
 WRITE(7,*) 'as xx tends to an infinity.'
 WRITE(7,*)
 WRITE(7,*) 'Table 7'
 WRITE(7,*)
 DO i=0,1
 WRITE(7,8) i
 WRITE(7,*) 'Function            znu                xx     ', &
    '    zans*znu**(1./3)   info'
   za = CMPLX(0,i)
   xx=10
   DO k=1,10
     znu=xx+za*xx**0.2
     CALL bessel(znu,xx,zbes,info)
     CALL test_info('zbessel ',znu,xx,info)
     WRITE(7,5) 'zbessel ',znu,xx,zbes*znu**(1/3._kp),info
     IF(info>25 .OR. xx>huge1/200) EXIT
     xx=xx*100
   ENDDO
   WRITE(7,*)
 ENDDO

 DO i=0,1
 WRITE(7,8) i
 WRITE(7,*) 'Function            znu                xx     ', &
    '    zans*znu**(1./3)   info'
   za=CMPLX(0,i)
   xx=10
   DO k=1,10
     znu=xx+za*xx**0.2
     CALL neumann(znu,xx,zneu,info)
     CALL test_info('zneumann',znu,xx,info)
     WRITE(7,5) 'zneumann',znu,xx,zneu*znu**(1/3._kp),info
     IF(info>25 .OR. xx>huge1/200) EXIT
     xx=xx*100
   ENDDO
   WRITE(7,*)
 ENDDO


! Test 8
! Test 8 ascertains that each subroutine offers correct zans, info in various
! cases.  Tests 8-6, 8-7 also show examples of underflows of zans.
! Test 8-1: REAL(znu) increases with xx=0.2.
! Test 8-2: REAL(znu) increases with xx=1.2.
! Test 8-3: ABS(znu) increases with xx=0.05.
!           Test 8-3 shows cases of info=0,5,10,30.
! Test 8-4: ABS(znu) increases with xx=2.
!           Test 8-4 shows cases of info=0,5,10,30.
! Test 8-5: xx increases with znu=(2, 0).
!           Test 8-5 shows cases of info=0,5,10,30.
! Test 8-6: xx tends to zero with znu=(1.5, 0).
! Test 8-7: xx tends to zero with znu=(1.5, 2).
! Test 8-8: Tests on zbes, info of SUBROUTINE bessel at xx=0.
! Test 8-9: Negative arguments: xx<0.
!           Test 8-9 shows cases of info=30.
!
! Tables 8 show zans, info for each test point (znu,xx) and each function.
!
! huge2, tiny2 are defined as
! huge2=MIN(HUGE(1._kp)/3,1E500),   tiny2=MAX(TINY(1._kp),1E-500).
 a2=HUGE(1._kp)/3;  a1=1
 DO i=1,55
   IF(a1 > a2*1E-20) THEN;  huge2=a2;   EXIT;   ENDIF
   a1=a1*1E20
   IF(i >= 25) THEN;   huge2=a1;   EXIT;   ENDIF
 ENDDO
 a2=TINY(1._kp)
 a1=1
 DO i=1,55
   IF(a1 < a2*1E+20) THEN;     tiny2=a2;   EXIT;   ENDIF
   a1=a1*1E-20
   IF(i >= 25) THEN;  tiny2=a1;   EXIT;   ENDIF
 ENDDO
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 8'
 WRITE(7,*) 'Tests 8 ascertains that each subroutine offers correct zans,' 
 WRITE(7,*) 'info in various cases.'
 WRITE(7,*)
! Tests 8-1,2
 DO i=0,1
   WRITE(7,2) ' Table 8-', i+1
   WRITE(7,*) 'Function            znu                xx       ', &
      '       zans         info'
   xx=0.2+i
   DO k=8,13
     znu=10*k+i*40
     CALL bessel(znu,xx,zans,info)
     CALL test_info('zbessel ',znu,xx,info)
     WRITE(7,5) 'zbessel ',znu,xx,zans,info
   ENDDO
   WRITE(7,*)
   DO k=8,13
     znu=10*k+i*40
     CALL neumann(znu,xx,zans,info)
     CALL test_info('zneumann',znu,xx,info)
     WRITE(7,5) 'zneumann',znu,xx,zans,info
   ENDDO
   WRITE(7,*)
 ENDDO

! Test 8-3
 WRITE(7,2) ' Table 8-3'
 WRITE(7,*) 'Function            znu                xx       ', &
    '       zans         info'
 xx=0.05_kp
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL bessel(znu,xx,zans,info)
   CALL test_info('zbessel ',znu,xx,info)
   WRITE(7,5) 'zbessel ',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL neumann(znu,xx,zans,info)
   CALL test_info('zneumann',znu,xx,info)
   WRITE(7,5) 'zneumann',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL hankel1(znu,xx,zans,info)
   CALL test_info('zhankel1',znu,xx,info)
   WRITE(7,5) 'zhankel1',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL hankel2(znu,xx,zans,info)
   CALL test_info('zhankel2',znu,xx,info)
   WRITE(7,5) 'zhankel2',znu,xx,zans,info
 ENDDO
 WRITE(7,*)

! Test 8-4
 WRITE(7,2) ' Table 8-4'
 WRITE(7,*) 'Function            znu                xx       ', &
    '       zans         info'
 xx=2
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL bessel(znu,xx,zans,info)
   CALL test_info('zbessel ',znu,xx,info)
   WRITE(7,5) 'zbessel ',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL neumann(znu,xx,zans,info)
   CALL test_info('zneumann',znu,xx,info)
   WRITE(7,5) 'zneumann',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL hankel1(znu,xx,zans,info)
   CALL test_info('zhankel1',znu,xx,info)
   WRITE(7,5) 'zhankel1',znu,xx,zans,info
 ENDDO
 WRITE(7,*)
 DO i=1,10
   bb=0.5_kp/i
   znu=zunit*EXP((pi/2-bb)/TAN(bb)+1-zunit*bb)*xx/2
   CALL hankel2(znu,xx,zans,info)
   CALL test_info('zhankel2',znu,xx,info)
   WRITE(7,5) 'zhankel2',znu,xx,zans,info
 ENDDO
 WRITE(7,*)

! Test 8-5
 WRITE(7,2) ' Table 8-5'
 WRITE(7,*) 'Function            znu                xx       ', &
    '       zans         info'
 znu=2;  xx=100
 DO k=1,9
   CALL hankel1(znu,xx,zans,info)
   CALL test_info('zhankel1',znu,xx,info)
   WRITE(7,5) 'zhankel1',znu,xx,zans,info
   xx=xx*100
 ENDDO
 WRITE(7,*)
 xx=100
 DO k=1,9
   CALL hankel2(znu,xx,zans,info)
   CALL test_info('zhankel2',znu,xx,info)
   WRITE(7,5) 'zhankel2',znu,xx,zans,info
   xx=xx*100
 ENDDO
 WRITE(7,*)

! Tests 8-6,7
 DO i=0,1
   WRITE(7,2) ' Table 8-', i+6
   WRITE(7,*) 'Function            znu                xx       ', &
      '       zans         info'
   znu=CMPLX(1.5, i*2)
   xx=tiny2**0.520    !   1E-160_kp
   DO k=1,7
     CALL bessel(znu,xx,zans,info)
     CALL test_info('zbessel ',znu,xx,info)
     WRITE(7,5) 'zbessel ',znu,xx,zans,info
     xx=xx*1E-20
   ENDDO
   WRITE(7,*)
   xx=tiny2**0.520    !   1E-160_kp
   DO k=1,7
     CALL neumann(znu,xx,zans,info)
     CALL test_info('zneumann',znu,xx,info)
     WRITE(7,5) 'zneumann',znu,xx,zans,info
     xx=xx*1E-20
   ENDDO
   WRITE(7,*)
 ENDDO

! Test 8-8.     Test for x=0.
 WRITE(7,*) 'Table 8-8' 
 xx=0
 WRITE(7,*) 'Function            znu                xx       ', &
    '       zans          info'
 DO i=-2,2
   znu=i*1.5
   CALL bessel(znu,xx,zbes,info)
   CALL test_info('zbessel ',znu,xx,info)
   WRITE(7,7) 'zbessel ',znu,xx,zbes,info
 ENDDO
 DO i=-1,1
   DO j=-2,2
     znu=CMPLX(i*2,j*1.5)
     CALL bessel(znu,xx,zbes,info)
     CALL test_info('zbessel ',znu,xx,info)
     WRITE(7,7) 'zbessel ',znu,xx,zbes,info
   ENDDO
 ENDDO
 WRITE(7,*)
7 FORMAT(A9, ' (', ES10.2E3, ',', ES10.2E3, ')', ES11.2E3, ' (', ES10.1E4, ',',&
     ES10.1E4, ')', I4)

! Test 8-9.     Test for x<0.
 WRITE(7,*) 'Table 8-9'
 WRITE(7,*) 'Function            znu                xx       ', &
    '       zans         info'
 znu=1;  xx=-1
 CALL bessel(znu,xx,zbes,info)
 CALL test_info('zbessel ',znu,xx,info)
 WRITE(7,5) 'zbessel ',znu,xx,zbes,info
 CALL neumann(znu, xx, zneu, info)
 CALL test_info('zneumann',znu,xx,info)
 WRITE(7,5) 'zneumann',znu,xx,zneu,info
 CALL hankel1(znu,xx,zhan1,info)
 CALL test_info('zhankel1',znu,xx,info)
 WRITE(7,5) 'zhankel1',znu,xx,zhan1,info
 CALL hankel2(znu,xx,zhan2,info)
 CALL test_info('zhankel2',znu,xx,info)
 WRITE(7,5) 'zhankel2',znu,xx,zhan2,info


! Test 9
! Tests for large ABS(znu) or large xx
! Test 9 ascertains that any overflows do not occur, and that the processor 
! does not stop the execution of the test program halfway.
!
 WRITE(7,'(//)')
 WRITE(7,*) 'Test 9'
 WRITE(7,*) 'Test 9 ascertains that the processor does not stop the execution&
          & of the test'
 WRITE(7,*) 'program halfway.'
 WRITE(7,*) 
 i5=8;   itotal=0
 DO j=0,i5
   theta=(pi/(2*i5))*j
   i6=NINT(i5*4*COS(theta))+1
   DO i=1,i6
     phi=(2*pi/i6)*i;  rr=5
     DO
       itotal=itotal+4
       IF(rr > huge2/1.5) EXIT
       rr=rr*1.38
       xx=rr*SIN(theta)+0.01
       znu=rr*CMPLX(COS(phi),SIN(phi))*COS(theta)

       CALL bessel(znu,xx,zans,info)
       CALL test_info('zbessel ',znu,xx,info)

       CALL neumann(znu,xx,zans,info)
       CALL test_info('zneumann',znu,xx,info)

       CALL hankel1(znu,xx,zans,info)
       CALL test_info('zhankel1',znu,xx,info)

       CALL hankel2(znu,xx,zans,info)
       CALL test_info('zhankel2',znu,xx,info)

     ENDDO
   ENDDO
 ENDDO

 xx=0.05
 DO j=0,200
   xx=xx*1E-2
   IF(xx < tiny2) EXIT
   rr=3
   DO
     rr=rr*1.57
     i6=rr
     DO i=1,i6
       itotal=itotal+4
       phi=(2*pi*i)/i6
       znu=rr*CMPLX(COS(phi),SIN(phi))

       CALL bessel(znu,xx,zans,info)
       CALL test_info('zbessel ',znu,xx,info)

       CALL neumann(znu,xx,zans,info)
       CALL test_info('zneumann',znu,xx,info)

       CALL hankel1(znu,xx,zans,info)
       CALL test_info('zhankel1',znu,xx,info)

       CALL hankel2(znu,xx,zans,info)
       CALL test_info('zhankel2',znu,xx,info)

     ENDDO
     IF(rr > 200) EXIT
   ENDDO
 ENDDO

 znu=0
 DO i=90,105
   itotal=itotal+2
   xx=i*0.01
   CALL bessel(znu,xx,zans,info)
   CALL test_info('zbessel ',znu,xx,info)
   CALL neumann(znu,xx,zans,info)
   CALL test_info('zneumann',znu,xx,info)
 ENDDO

 znus=(35.177_kp,20.535_kp)
 xxs=24
 DO i=-10,10
   DO j=-10,10
     itotal=itotal+2
     xx=xxs+i*1e-3
     znu=znus+j*1e-3
     CALL hankel1(znu,xx,zhan1,info)
     CALL test_info('zhankel1',znu,xx,info)
     CALL bessel(znu,xx,zbes,info)
     CALL test_info('zbessel ',znu,xx,info)
   ENDDO
 ENDDO

!  Testing overflows at SUBROUTINE neumann.
 DO i=-1,1
   DO j=-1,1
     znu=CMPLX(i*0.1-9.8,j*0.1-3.8)
     DO k=-1,1
       itotal=itotal+1
       xx=8.00E-31*(1+k*0.05)
       CALL neumann(znu,xx,zans,info)
       CALL test_info('zneumann',znu,xx,info)
     ENDDO
   ENDDO
 ENDDO

! Testing SUBROUTINE neumann at xx<<1 and nregion=2.
 xx=tiny1**0.292    !  1E-90_kp
 DO j=-4,6
   anu=j*0.03+0.15
   DO i=-2,2
     itotal=itotal+1
     theta=(pi/100)*i+pi/2
     znu=anu*CMPLX(COS(theta),SIN(theta))
     CALL neumann(znu,xx,zans,info)
     CALL test_info('zneumann',znu,xx,info)
   ENDDO
 ENDDO

! Test when xx is nearly equal to 0.
 xx=tiny1*3
 DO i=-6,6
   znu=-110.75_kp+0.1*i
   CALL bessel(znu,xx,zans,info)
 ENDDO

! Tests when zhankel2(znu,xx) is nearly equal to 0.
 znu0=(39.60473566960330_kp,-17.33456476536629_kp)
 xx0=30;  a1=3E-9
 DO i=-3,3
   DO j=-3,3
     DO k=-3,3
       znu=znu0+CMPLX(i,j)*a1
       xx=xx0+k*a1
       CALL hankel2(znu,xx,zhan2,info)
       IF(info > 27) PRINT *,'info(Test 9-5)=',info,i,j,k
     ENDDO
   ENDDO
 ENDDO

 WRITE(7,'(A,I7)') ' The number of the examined cylindrical functions=',itotal
 WRITE(7,*) 'The termination of test 9 without trouble'

 WRITE(7,*)
 WRITE(7,*) 'The end of file test_program_'//kp_num//'.out'
 END PROGRAM test_program

 SUBROUTINE sub_wronsky(znu,xx,info)
! Tests by the use of Wronskians.
! This subroutine is used for tests 1, 2, 3.
 USE mod_test, ONLY: bessel,neumann,hankel1,hankel2,err_max,id, &
   nreg_max,xx_max,znu_max,zunit,pi,kp,huge1,abs2,num_region
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 INTEGER,INTENT(OUT):: info
 COMPLEX(kp):: z0,z1,z2,z3,za1,za2,zb1,zb2
 REAL(kp):: a1,c1,error(4)
 INTEGER:: info0,info1,info2,info3,nregion
 c1=huge1**0.487     !  1E150_kp
 CALL bessel(znu,xx,za1,info0)
 CALL test_info('zbessel ',znu,xx,info0)
 CALL bessel(znu+1,xx,zb2,info1)
 CALL test_info('zbessel ',znu+1,xx,info1)
 CALL neumann(znu,xx,zb1,info2)
 CALL test_info('zneumann',znu,xx,info2)
 CALL neumann(znu+1,xx,za2,info3)
 CALL test_info('zneumann',znu+1,xx,info3)
 info=MAX(info0,info1,info2,info3)
 IF(info > 32) STOP 'info > 32'
 nregion=num_region(znu,xx)
 IF(info < 2) THEN
    id(nregion)=id(nregion)+1
    IF(abs2(za1) > abs2(za2)) THEN
       z0=za1;  za1=za2;  za2=z0
    ENDIF
    z3=-2/(pi*xx)
    IF(abs2(za2)>c1 .AND. abs2(za1)>1) THEN
      IF(abs2(zb1) > abs2(zb2)) THEN
        z0=zb1;  zb1=zb2;  zb2=z0
      ENDIF
      z1=za1;  z2=zb1*(zb2/za2);  z3=z3/za2
    ELSE
       z1=za1*za2;  z2=zb1*zb2
    ENDIF
    a1=MAX(abs2(z1),abs2(z2),abs2(z3))
    error(1)=abs2(z1-z2-z3)/a1
    IF(error(1) > err_max(1)) THEN
       znu_max(1)=znu;       xx_max(1)=xx
       err_max(1)=error(1);  nreg_max(1)=nregion
    ENDIF
 ENDIF
 CALL hankel1(znu,xx,za1,info0)
 CALL test_info('zhankel1',znu,xx,info0)
 CALL hankel1(znu+1,xx,zb2,info1)
 CALL test_info('zhankel1',znu+1,xx,info1)
 CALL hankel2(znu,xx,zb1,info2)
 CALL test_info('zhankel2',znu,xx,info2)
 CALL hankel2(znu+1,xx,za2,info3)
 CALL test_info('zhankel2',znu+1,xx,info3)
 info=MAX(info0,info1,info2,info3)
 IF(info > 32) STOP 'info > 32'
 IF(info < 2) THEN
   IF(abs2(za1) > abs2(za2)) THEN
     z0=za1;  za1=za2;  za2=z0
   ENDIF
   z3=4*zunit/(pi*xx)
   IF(abs2(za2)>c1 .AND. abs2(za1)>1) THEN
     IF(abs2(zb1) > abs2(zb2)) THEN
       z0=zb1;  zb1=zb2;  zb2=z0
     ENDIF
     z1=za1;  z2=zb1*(zb2/za2);  z3=z3/za2
   ELSE
     z1=za1*za2;  z2=zb1*zb2
   ENDIF
   a1=MAX(abs2(z1),abs2(z2),abs2(z3))
   error(2)=abs2(z1-z2-z3)/a1
   IF(error(2) > err_max(2)) THEN
     znu_max(2)=znu;       xx_max(2)=xx
     err_max(2)=error(2);  nreg_max(2)=nregion
   ENDIF
 ENDIF
 END SUBROUTINE sub_wronsky

 SUBROUTINE sub_conjg(znu, xx, info)
! This subroutine is used for test 4.
 USE mod_test, ONLY: kp,bessel,neumann,id,err_max,xx_max,znu_max,nreg_max, &
   tiny1,abs2,num_region
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 INTEGER,INTENT(OUT):: info
 COMPLEX(kp):: zans,zconans,zconnu
 REAL(kp):: error(4),a1
 REAL(kp),EXTERNAL:: amp
 INTEGER:: info0,info1,info2,nregion
 zconnu=CONJG(znu)
 nregion=num_region(znu,xx)
 CALL bessel(znu,xx,zans,info0)
 CALL test_info('zbessel ',znu,xx,info0)
 CALL bessel(zconnu,xx,zconans,info1)
 CALL test_info('zbessel ',zconnu,xx,info1)
 info2=MAX(info0,info1)
 IF(info2 <= 2) THEN
   id(nregion)=id(nregion)+1
   a1 = amp('bes',znu,xx,zans)
   IF(a1 > tiny1) THEN
     error(1)=abs2(zans-CONJG(zconans))/a1
     IF(error(1) > err_max(1)) THEN
       znu_max(1)=znu;        xx_max(1)=xx
       err_max(1)=error(1);   nreg_max(1)=nregion
     ENDIF
   ENDIF
 ENDIF
 info=info2
 CALL neumann(znu,xx,zans,info0)
 CALL test_info('zneumann',znu,xx,info0)
 CALL neumann(zconnu,xx,zconans,info1)
 CALL test_info('zneumann',zconnu,xx,info1)
 info2=MAX(info0, info1)
 IF(info2 <= 2) THEN
   a1=amp('neu',znu,xx,zans)
   IF(a1 > tiny1) THEN
     error(2)=abs2(zans-CONJG(zconans))/a1
     IF(error(2) > err_max(2)) THEN
       znu_max(2)=znu;       xx_max(2)=xx
       err_max(2)=error(2);  nreg_max(2)=nregion
     ENDIF
   ENDIF
 ENDIF
 info=MAX(info,info2)
5 FORMAT(A8,' (',F8.2,',',F8.2,')',F8.2,I5)
 END SUBROUTINE sub_conjg

 SUBROUTINE boundary_test(nregion1,nregion2,znu,xx)
! This SUBROUTINE is used for test 6.
! This SUBROUTINE invokes SUBROUTINEs bessel_boundary, neumann_boundary
! hankel2_boundary.
! The variable nregion in SUBROUTINE bessel_boundary is given by one of the 
! arguments, while the variable nregion in SUBROUTINE bessel is determined by
! FUNCTION num_region that is invoked by SUBROUTINE bessel. 
! SUBROUTINE bessel_boundary is the same as SUBROUTINE bessel except this.
! SUBROUTINEs neumann_boundary, hankel2_boundary are explained similarly.
 USE mod_test, ONLY: kp,err_max,xx_max,znu_max,tiny1, abs2
 IMPLICIT NONE
 INTEGER,INTENT(IN):: nregion1,nregion2
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp):: zans0,zans1
 REAL(kp):: a1,err
 REAL(kp),EXTERNAL:: amp
 INTEGER:: info,info0,info1
!
 CALL bessel_boundary(znu,xx,nregion1,zans0,info0)
 IF(info0 > 2) PRINT 1,'zbessel ',znu,xx,info0
 CALL bessel_boundary(znu,xx,nregion2,zans1,info1)
 IF(info1 > 2) PRINT 1,'zbessel ',znu,xx,info1
 info=MAX(info0,info1)
 a1=amp('bes',znu,xx,zans0)
 IF(a1 > tiny1) THEN
   IF(info > 2) THEN
     WRITE(*,1) 'zbessel ',znu,xx,info
   ELSE
     err=abs2(zans1-zans0)/a1
     IF(err > err_max(1)) THEN
       err_max(1)=err;   znu_max(1)=znu;  xx_max(1)=xx
     ENDIF
   ENDIF
 ENDIF

 CALL neumann_boundary(znu,xx,nregion1,zans0,info0)
 IF(info0 > 2) PRINT 1,'zneumann',znu,xx,info0
 CALL neumann_boundary(znu,xx,nregion2,zans1,info1)
 IF(info1 > 2) PRINT 1,'zneumann',znu,xx,info1
 info=MAX(info0,info1)
 a1=amp('neu',znu,xx,zans0)
 IF(a1 > tiny1) THEN
   IF(info > 2) THEN
     WRITE(*,1) 'zneumann',znu,xx,info
   ELSE
     err=abs2(zans1-zans0)/a1
     IF(err > err_max(2)) THEN
       err_max(2)=err;   znu_max(2)=znu;  xx_max(2)=xx
     ENDIF
   ENDIF
 ENDIF

 CALL hankel2_boundary(znu,xx,nregion1,zans0,info0)
 IF(info0 > 2) PRINT 1,'zhankel2',znu,xx,info0
 CALL hankel2_boundary(znu,xx,nregion2,zans1,info1)
 IF(info1 > 2) PRINT 1,'zhankel2',znu,xx,info1
 info=MAX(info0,info1)
 a1=amp('hn2',znu,xx,zans0)
 IF(a1 > tiny1) THEN
   IF(info > 2) THEN
     WRITE(*,1) 'zhankel2',znu,xx,info
   ELSE
     err=abs2(zans1-zans0)/a1
     IF(err > err_max(3)) THEN
       err_max(3)=err;   znu_max(3)=znu;  xx_max(3)=xx
     ENDIF
   ENDIF
 ENDIF
1 FORMAT(A8,' (',F8.2,',',F8.2,')',F8.2,I5)
 END SUBROUTINE boundary_test

 SUBROUTINE bessel_boundary(znu,xx,nregion,zans,info)
! This subroutine is used for test 6.
! This SUBROUTINE is invoked by SUBROUTINE boundary_test.
! SUBROUTINE bessel_boundary calculates zans.
! zans=zbessel_bdry(nregion,znu,xx)
! zbessel_bdry(nregion,znu,xx) is defined at test 6 of PROGRAM test_program.
! The variable nregion in SUBROUTINE bessel_boundary is given by one of the
! arguments, while the variable nregion in SUBROUTINE bessel is determined by 
! FUNCTION num_region. SUBROUTINE bessel_boundary is the same as SUBROUTINE 
! bessel except this.
 USE mod_test, ONLY: err_range1,err_range2,err_range3,hugelog,huge1,ihuge,kp, &
   pi,zunit,abs2,bes_series,neu_series,bes_han_dby,han2_olver,bes_olver, &
   bes_recur,han2_temme
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
 COMPLEX(kp):: zbesa,zepspi,zhan2a,zlogbes,zlogmu,zneua,znua,zconnu,zconnupi, &
   znupi,zss1,zss2,zsum,zconepspi
 REAL(kp):: aim_znu,re_znu,error,error1,error2
 INTEGER:: info1,info2,nregion,nn
! This SUBROUTINE calculates zans.  zans=zbessel(znu,xx).
! This is invoked only by SUBROUTINE boundary_test.
! The dummy arguments znu, xx, zans, info, and the outline of the calculation
! in SUBROUTINE neumann are explained in the beginning of this file.
! Determination of nregion
 zans=0
 IF(nregion == 0) THEN;   info=30;   RETURN;  ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 zconnupi=zconnu*pi;   znua=-znu
 IF(re_znu < -ihuge) THEN;   info=30;     RETURN;  ENDIF
 IF(re_znu < 0) THEN
   nn=NINT(re_znu)
   zepspi=(znu-nn)*pi
   zconepspi=CONJG(zepspi)
 ENDIF
 SELECT CASE(nregion)
 CASE(1)
   CALL bes_series(znu,xx,zsum,zlogbes,error,info)
   IF(info > 17) THEN;    zans=huge1;    RETURN;   ENDIF
   IF(REAL(zlogbes) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans = zsum*EXP(zlogbes)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)
   IF(re_znu >= 0) THEN
     CALL bes_series(znu,xx,zsum,zlogbes,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(REAL(zlogbes) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=zsum*EXP(zlogbes)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0, zbessel(znu,xx) is calculated using Eq. (40a) of Ref. (1).
   CALL bes_series(znua,xx,zsum,zlogbes,error,info1)
   IF(info1 > 17) THEN;   zans=huge1;   info=info1;   RETURN;   ENDIF
   IF(REAL(zlogbes) > hugelog) THEN;   zans=huge1;   info=20;    RETURN;   ENDIF
   IF(error > err_range3) THEN;   info=30;     RETURN;   ENDIF
   zbesa=zsum*EXP(zlogbes)
   CALL neu_series(znua,xx,zneua,info2)
   IF(info2 > 17) THEN;   zans=huge1;   info=info2;   RETURN;   ENDIF
   info=MAX(info1,info2)
   error=abs2(znupi)+error
   IF(error > err_range3) THEN;  info=30;  RETURN;   ENDIF
   zans=zbesa*COS(znupi)+zneua*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
   IF(error > err_range1) info=5
 CASE(3)  !  by Debye's method       |znu/xx|<=1
! zbessel(znu,xx) is calculated using Eq. (40b) of Ref. (1).
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(ABS(REAL(zlogmu)) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=(zss1*EXP(zlogmu)+zss2*EXP(-zlogmu))/2
   IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
   IF(error > err_range1) info=5
 CASE(4)        !  by Debye's method       |znu/xx|>1
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;     zans=huge1;     RETURN;       ENDIF
       IF(REAL(zlogmu) > hugelog) THEN;  zans=huge1;  info=20;   RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans=zss1*EXP(zlogmu)/2
       IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zbessel(znu,xx) is calculated 
! using Eq. (40c) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;    RETURN;    ENDIF
     IF(REAL(zlogmu) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;       info=30;     RETURN;     ENDIF
     zans=CONJG(zss1*EXP(zlogmu))/2
     IF(error > err_range2) THEN;  info=10;  RETURN;   ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zbessel(znu,xx) is calculated 
! using Eq. (40d) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;    zans=huge1;     RETURN;     ENDIF
     IF(AIMAG(zconnupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(MAX(REAL(zlogmu-zunit*zconnupi), &
              REAL(-zlogmu)+ABS(AIMAG(zconepspi))) > hugelog) THEN
       zans=huge1;    info=20;    RETURN
     ENDIF
     error=error+abs2(znupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=EXP(zlogmu-zunit*zconnupi)*zss1/2 &
          +zunit*zss2*EXP(-zlogmu)*SIN(zconnupi)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zbessel(znu,xx) is calculated by 
! using Eq. (40e) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(MAX(REAL(zlogmu+zunit*zepspi), &
           REAL(-zlogmu)+ABS(AIMAG(zepspi))) > hugelog) THEN
     zans=huge1;     info=20;     RETURN
   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=EXP(zlogmu-zunit*znupi)*zss1/2 &
        +zunit*zss2*EXP(-zlogmu)*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)        !  Method by the Olver's asymptotic series
   IF(re_znu >= 0) THEN
     CALL bes_olver(znu, xx, zans, error, info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zbessel(znu,xx) is calculated by 
! using Eq. (40d) of Ref. (1).
     CALL bes_olver(-zconnu,xx,zbesa,error1,info1)
     CALL han2_olver(-zconnu,xx,zhan2a,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     error=ABS(REAL(znupi,kp))+error2
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(ABS(AIMAG(zconepspi)) > hugelog) THEN;  info=30;  RETURN;  ENDIF
     zans=zbesa*EXP(-zunit*zconnupi)+zunit*zhan2a*SIN(zconnupi)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zbessel(znu,xx) is calculated by 
! using Eq. (40e) of Ref. (1).
   CALL bes_olver(znua,xx,zbesa,error1,info1)
   CALL han2_olver(znua,xx,zhan2a,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   error=ABS(REAL(znupi,kp))+error2
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=EXP(-zunit*znupi)*zbesa+zunit*zhan2a*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)
   CALL bes_recur(znu,xx,0,zans,info)
 CASE(7)
   IF(re_znu >= 0) THEN
     CALL bes_recur(znu,xx,0,zans,info)
     RETURN
   ENDIF
! When re_znu<0, zbessel(znu,xx) is calculated using Eq. (40e) of Ref. (1).
   CALL bes_recur(znua,xx,1,zbesa,info1)
   CALL han2_temme(znua,xx,zhan2a,info2)
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zbesa+zunit*zhan2a*SIN(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 END SELECT
 END SUBROUTINE bessel_boundary

 SUBROUTINE neumann_boundary(znu,xx,nregion,zans,info)
! This subroutine is used for test 6.
! This SUBROUTINE is invoked by SUBROUTINE boundary_test.
! SUBROUTINE neumann_boundary calculates zans.
! zans=zneumann_bdry(nregion,znu,xx)
! zneumann_bdry(nregion,znu,xx) is defined at test 6 of PROGRAM test_program.
! The variable nregion in SUBROUTINE neumann_boundary is given by one of the
! arguments, while the variable nregion in SUBROUTINE neumann is determined by 
! FUNCTION num_region. SUBROUTINE neumann_boundary is the same as SUBROUTINE 
! bessel except this.
 USE mod_bes, ONLY: err_range1,err_range2,err_range3,hugelog,theta_lim, &
   huge1,kp,pi,zunit,han2_temme, &
   abs2,bes_series,neu_series,bes_han_dby,bes_olver,han2_olver,bes_recur
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
 COMPLEX(kp):: zarg1,zarg2,zbes,zbesa,zhan2,zhan2a,zlogbes1, &
  zlogbes2,zlogmu,zneua,znua,zconnu,zconnupi,znupi,zpart1,zpart2, &
  zss1,zss2,zsum1,zsum2,z1
 REAL(kp):: aim_znu,re_znu,abszpart1,abszpart2
 REAL(kp):: error,error1,error2,abes,ahan2
 INTEGER:: info1,info2,nregion
! This SUBROUTINE calculates zans.  zans=zneumann(znu,xx).
! This is invoked only by SUBROUTINE boundary_test.
! The dummy arguments znu, xx, zans, info, and the outline of the calculation
! in SUBROUTINE neumann are explained in the beginning of this file.
! Determination of nregion
 zans=0
 IF(nregion == 0) THEN;  info=30;   RETURN;   ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 zconnupi=zconnu*pi;   znua=-znu
 SELECT CASE(nregion)
 CASE(1)        !  With Taylor's expansion.
   CALL bes_series(znu,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;     zans=huge1;     info=info1;     RETURN;   ENDIF
   IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   CALL bes_series(znua,xx,zsum2,zlogbes2,error2,info2)
   IF(info2 > 17) THEN;  zans=huge1;   info=info2;     RETURN;   ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>=0, zneumann(znu,xx) is calculated using Eq. (41b) of Ref. (1).
     zarg1=zlogbes1+2*zunit*znupi
     zarg2=zlogbes2+zunit*znupi
     IF(REAL(zarg1) > hugelog .OR. REAL(zarg2) > hugelog) THEN
       zans=huge1;     info=20;     RETURN;     ENDIF
     IF((ABS(AIMAG(zarg1))>theta_lim) .OR. &
          (ABS(AIMAG(zlogbes1))>theta_lim)) THEN;  info=30;  RETURN;  ENDIF
     zpart1=zsum1*(EXP(zarg1)+EXP(zlogbes1))
     zpart2=-2*zsum2*EXP(zarg2)
     zans=zunit*(zpart1+zpart2)/(EXP(2*zunit*znupi)-1)
   ELSE
! When aim_znu<0, zneumann(znu,xx) is calculated using Eq. (41c) of Ref. (1).
     zarg1=zlogbes1-2*zunit*znupi
     zarg2=zlogbes2-zunit*znupi
     IF(REAL(zarg1) > hugelog .OR. REAL(zarg2) > hugelog) THEN
       zans=huge1;  info=20;   RETURN;   ENDIF
     zpart1=zsum1*(EXP(zlogbes1)+EXP(zarg1))
     zpart2=-2*zsum2*EXP(zarg2)
     zans=zunit*(zpart1+zpart2)/(1-EXP(-2*zunit*znupi))
   ENDIF
   info=MAX(info1,info2)
   error2=error2+abs2(znupi)
   abszpart1=abs2(zpart1)
   abszpart2=abs2(zpart2)
   error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)        !  With Taylor's expansion.   ABS(zeps) < 0.3
   IF(re_znu >= 0) THEN
      CALL neu_series(znu,xx,zans,info)
      RETURN
   ENDIF
! When re_znu<0, zneumann(znu,xx) is calculated using Eq. (41a) of Ref. (1).
   CALL bes_series(znua,xx,zsum1,zlogbes1,error,info1)
   IF(info1 > 17) THEN;   zans=huge1;    info=info1;     RETURN;   ENDIF
   IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zbesa=zsum1*EXP(zlogbes1)
   CALL neu_series(znua,xx,zneua,info2)
   IF(info2 > 17) THEN;  zans=huge1;   info=info2;   RETURN;   ENDIF
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=-zbesa*SIN(znupi)+zneua*COS(znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(3)        !  With Debye's method         |znu/xx|<=1
! zneumann(znu,xx) is calculated using Eq. (41d) of Ref. (1).
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(ABS(REAL(zlogmu)) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=-zunit*(zss1*EXP(zlogmu)-zss2*EXP(-zlogmu))/2
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info = 5
 CASE(4)        !  With Debye's method         |znu/xx|>1
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
! When re_znu>=0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41e) of Ref. (1).
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;  zans=huge1;   RETURN;  ENDIF
       IF(ABS(REAL(zlogmu)) > hugelog) THEN
         zans=huge1;  info=20;  RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans = zunit*(zss2*EXP(-zlogmu)-zss1*EXP(zlogmu)/2)
       IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41f) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(ABS(REAL(zlogmu,kp)) > hugelog) THEN
       zans=huge1;  info=20;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=-zunit*CONJG(zss2*EXP(-zlogmu)-zss1*EXP(zlogmu)/2)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41g) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;     RETURN;     ENDIF
     IF(AIMAG(zconnupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
     IF(MAX(REAL(zlogmu-zconnupi*zunit), &
              REAL(-zlogmu)+ABS(AIMAG(znupi))) > hugelog) THEN
       zans=huge1;     info=20;     RETURN;    ENDIF
     error=error+abs2(zconnupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=zunit*(-EXP(zlogmu-zunit*zconnupi)*zss1/2 &
           +zss2*EXP(-zlogmu)*COS(zconnupi))
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41h) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;    zans=huge1;     RETURN;   ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   IF(MAX(REAL(zlogmu-znupi*zunit),REAL(-zlogmu)+ABS(AIMAG(znupi))) &
       > hugelog) THEN;  zans=huge1;   info=20;   RETURN;   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zunit*(-EXP(zlogmu-zunit*znupi)*zss1/2 &
         +zss2*EXP(-zlogmu)*COS(znupi))
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)        !  Method by the Olver's asymptotic series
   IF(re_znu >= 0) THEN
! When re_znu>=0, zneumann(znu,xx) is calculated using Eq. (41e) of Ref. (1).
     CALL bes_olver(znu,xx,zbes,error1,info1)
     CALL han2_olver(znu,xx,zhan2,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     abes=abs2(zbes);  ahan2=abs2(zhan2)
     error=(abes*error1+ahan2*error2)/(abes+ahan2)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     zans = zunit*(zhan2-zbes)
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zneumann(znu,xx) is calculated 
! using Eq. (41g) of Ref. (1).
     CALL bes_olver(-zconnu,xx,zbesa,error1,info1)
     CALL han2_olver(-zconnu,xx,zhan2a,error2,info2)
     info=MAX(info1,info2)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     error=ABS(REAL(znupi,kp))+error2
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(ABS(AIMAG(zconnupi)) > hugelog) THEN;  info=30;  RETURN;  ENDIF
     zans=zunit*(-EXP(-zunit*zconnupi)*zbesa+COS(zconnupi)*zhan2a)
     zans=CONJG(zans)
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zneumann(znu,xx) is calculated 
! using Eq. (41h) of Ref. (1).
   CALL bes_olver(znua,xx,zbesa,error1,info1)
   CALL han2_olver(znua,xx,zhan2a,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(AIMAG(znupi) < -hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
   error = ABS(REAL(znupi,kp))+error2
   IF(error > err_range3) THEN;  info=30;  RETURN;   ENDIF
   zans = zunit*(-EXP(-zunit*znupi)*zbesa+COS(znupi)*zhan2a)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)        !  By recurrence formula.
   CALL bes_recur( znu,xx,0,zbes,info1)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>=0, zneumann(znu,xx) is calculated using Eq. (41b) of Ref. (1).
     CALL bes_recur(-znu,xx,-1,zbesa,info2)
     info=MAX(info1,info2)
     z1=EXP(2*zunit*znupi)
     zans=zunit*(zbes*(z1+1)-2*zbesa)/(z1-1)
   ELSE
! When aim_znu<0, zneumann(znu,xx) is calculated using Eq. (41c) of Ref. (1).
     CALL bes_recur(-znu,xx,1,zbesa,info2)
     info=MAX(info1,info2)
     z1=EXP(-2*zunit*znupi)
     zans=zunit*(zbes*(z1+1)-2*zbesa)/(1-z1)
   ENDIF
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(7)        !  By recurrence formula.   ABS(AIMAG(znu)) < 0.5
   IF(re_znu >= 0) THEN
! When re_znu>=0, zneumann(znu,xx) is calculated using Eq. (41e) of Ref. (1).
     CALL bes_recur(znu,xx,0,zbes,info1)
     CALL han2_temme(znu,xx,zhan2,info2)
     info=MAX(info1,info2)
     zans=zunit*(zhan2-zbes)
     RETURN
   ENDIF
! When re_znu<0, zneumann(znu,xx) is calculated using Eq. (41h) of Ref. (1).
   CALL bes_recur(-zconnu,xx,1,zbesa,info1)
   CALL han2_temme(-zconnu,xx,zhan2a,info2)
   info=MAX(info1,info2)
   error=ABS(REAL(znupi,kp))
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zunit*(-zbesa+COS(zconnupi)*zhan2a)
   zans=CONJG(zans)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 END SELECT
 END SUBROUTINE neumann_boundary

 SUBROUTINE hankel2_boundary(znu,xx,nregion,zans,info)
! This subroutine is used for test 6.
! This SUBROUTINE is invoked by SUBROUTINE boundary_test.
! SUBROUTINE hankel2_boundary calculates zans.
! zans=zhankel2_bdry(nregion,znu,xx)
! zhankel2_bdry(nregion,znu,xx) is defined at test 6 of PROGRAM test_program.
! The variable nregion in SUBROUTINE hankel2_boundary is given by one of the
! arguments, while the variable nregion in SUBROUTINE hankel2 is determined by
! FUNCTION num_region. SUBROUTINE hankel2_boundary is the same as SUBROUTINE 
! hankel2 except this.
 USE mod_bes, ONLY: err_range1,err_range2,err_range3,hugelog,huge1,kp,pi, &
   tiny1,zunit,han2_temme,abs2,bes_series,neu_series,bes_han_dby,bes_olver, &
   han2_olver,bes_recur
 IMPLICIT NONE
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp),INTENT(OUT):: zans
 INTEGER,INTENT(OUT):: info
 COMPLEX(kp):: zarg1,zarg2,zbes,zbes1,zbes2,zlogbes1,zlogbes2, &
  zlogmu,zneu,znua,zconnu,znupi,zpart1,zpart2,zss1,zss2,zsum1,zsum2,z1,z2
 REAL(kp):: aim_znu,re_znu,abszpart1,abszpart2
 REAL(kp):: error,error1,error2
 INTEGER:: info1,info2,nregion
! This SUBROUTINE calculates zans.  zans=zhankel2(znu,xx).
! This is invoked by  SUBROUTINE boundary_test and SUBROUTINE hankel1.
! The dummy arguments znu, xx, zans, info, and the outline of the calculation
! in SUBROUTINE hankel2 are explained in the beginning of this file.
! Determination of nregion
 zans=0
 IF(nregion == 0) THEN;   info=30;   RETURN;  ENDIF
 re_znu=REAL(znu,kp);  aim_znu=AIMAG(znu)
 znupi=znu*pi;         zconnu=CONJG(znu)
 znua=-znu
 SELECT CASE(nregion)
 CASE(1)
   CALL bes_series(znu,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;   zans=huge1;   info=info1;    RETURN;   ENDIF
   IF(error1 > err_range3) THEN;   info=30;     RETURN;    ENDIF
   CALL bes_series(-znu,xx,zsum2,zlogbes2,error2,info2)
   info=MAX(info1,info2)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   IF(error2 > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(aim_znu >= 0) THEN
! When aim_znu>0, zhankel2(znu,xx) is calculated using Eq. (42a) of Ref. (1).
     zarg1=zlogbes1+2*zunit*znupi
     zarg2=zlogbes2+zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;     info=20;     RETURN
     ENDIF
     zpart1=zsum1*EXP(zarg1);  zpart2=-zsum2*EXP(zarg2)
     zans=2*(zpart1+zpart2)/(EXP(2*zunit*znupi)-1)
     error1=error1+2*abs2(znupi);  error2=error2+abs2(znupi)
     abszpart1=abs2(zpart1);  abszpart2=abs2(zpart2)
     IF(abszpart1+abszpart2 < tiny1) THEN
       info=20;  RETURN
     ELSE
       error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
     ENDIF
   ELSE
! When aim_znu<=0, zhankel2(znu,xx) is calculated using Eq. (42b) of Ref. (1).
     zarg1=zlogbes1
     zarg2=zlogbes2-zunit*znupi
     IF(REAL(zarg1)>hugelog .OR. REAL(zarg2)>hugelog) THEN
       zans=huge1;     info=20;     RETURN
     ENDIF
     zpart1=zsum1*EXP(zarg1);    zpart2=-zsum2*EXP(zarg2)
     zans=2*(zpart1+zpart2)/(1-EXP(-2*zunit*znupi))
     error2=abs2(znupi)+error2
     abszpart1=abs2(zpart1);     abszpart2=abs2(zpart2)
     error=(abszpart1*error1+abszpart2*error2)/(abszpart1+abszpart2)
   ENDIF
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(2)
   IF(re_znu >= 0) THEN
! When re_znu>=0, zhankel2(znu,xx) is calculated using Eq. (42c) of Ref. (1).
     CALL bes_series(znu,xx,zsum1,zlogbes1,error,info1)
     IF(info1 > 17) THEN;     zans=huge1;     info=info1;     RETURN;    ENDIF
     IF(REAL(zlogbes1) > hugelog) THEN;  zans=huge1;  info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;  info=30;   info=info1;   RETURN;   ENDIF
     zbes=zsum1*EXP(zlogbes1)
     IF(error > err_range1) info1=5
     IF(error > err_range2) info1=10
     CALL neu_series(znu,xx,zneu,info2)
     info=MAX(info1,info2)
     zans=zbes-zunit*zneu
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42d) of Ref. (1).
   CALL bes_series(znua,xx,zsum1,zlogbes1,error1,info1)
   IF(info1 > 17) THEN;  zans=huge1;   info=info1;   RETURN;   ENDIF
   zarg1=zlogbes1+zunit*znupi
   IF(REAL(zarg1) > hugelog) THEN;   zans=huge1;   info=20;    RETURN;   ENDIF
   error=abs2(znupi)+error1
   IF(error > err_range3) THEN;    info=30;     RETURN;   ENDIF
   CALL neu_series(znua,xx,zneu,info2)
   info=MAX(info1,info2)
   zans=EXP(zarg1)*zsum1-EXP(zunit*znupi)*zunit*zneu
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info = 5
 CASE(3)
   CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   IF(REAL(-zlogmu) > hugelog) THEN;   zans=huge1;   info=20;   RETURN;   ENDIF
   IF(error > err_range3) THEN;   info=30;   RETURN;   ENDIF
   zans=zss2*EXP(-zlogmu)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(4)
   IF(re_znu >= 0) THEN
     IF(aim_znu >= 0) THEN
       CALL bes_han_dby(znu,xx,zss1,zss2,zlogmu,error,info)
       IF(info > 17) THEN;   zans=huge1;     RETURN;    ENDIF
       IF(REAL(-zlogmu) > hugelog) THEN;  zans=huge1;  info=20;  RETURN;  ENDIF
       IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
       zans=zss2*EXP(-zlogmu)
       IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
       IF(error > err_range1) info=5
       RETURN
     ENDIF
! When re_znu>=0 and aim_znu<0, zhankel2(znu,xx) is calculated
! using Eq. (42f) of Ref. (1).
     CALL bes_han_dby(zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;   zans=huge1;     RETURN;     ENDIF
     IF(ABS(REAL(zlogmu)) > hugelog) THEN
       zans=huge1;  info=20;   RETURN;   ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=CONJG(zss1*EXP(zlogmu)-zss2*EXP(-zlogmu))
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
   IF(aim_znu >= 0) THEN
! When re_znu<0 and aim_znu>=0, zhankel2(znu,xx) is calculated 
! using Eq. (42g) of Ref. (1).
     CALL bes_han_dby(-zconnu,xx,zss1,zss2,zlogmu,error,info)
     IF(info > 17) THEN;    zans=huge1;     RETURN;     ENDIF
     z1=zlogmu+CONJG(zunit*znupi)
     z2=-zlogmu+CONJG(zunit*znupi)
     IF(MAX(REAL(z1),REAL(z2)) > hugelog) THEN
       zans=huge1;  info=20;  RETURN;  ENDIF
     error=error+abs2(znupi)
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     zans=CONJG(zss1*EXP(z1)-zss2*EXP(z2))
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0 and aim_znu<0, zhankel2(znu,xx) is calculated
! using Eq. (42e) of Ref. (1).
   CALL bes_han_dby(znua,xx,zss1,zss2,zlogmu,error,info)
   IF(info > 17) THEN;     zans=huge1;     RETURN;   ENDIF
   z1=zunit*znupi-zlogmu
   IF(REAL(z1) > hugelog) THEN;  zans=huge1;  info=20;    RETURN;   ENDIF
   error=error+abs2(znupi)
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   zans=zss2*EXP(z1)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(5)
   IF(re_znu >= 0) THEN
     CALL han2_olver(znu,xx,zans,error,info)
     IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
     IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
     IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
     IF(error > err_range1) info=5
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42e) of Ref. (1).
   CALL han2_olver(znua,xx,zans,error,info)
   IF(info > 17) THEN;  zans=huge1;  RETURN;  ENDIF
   error=ABS(REAL(znupi,kp))+error
   IF(error > err_range3) THEN;  info=30;  RETURN;  ENDIF
   IF(REAL(zunit*znupi) > hugelog) THEN;  info=30;  RETURN;  ENDIF
   zans=zans*EXP(zunit*znupi)
   IF(error > err_range2) THEN;  info=10;  RETURN;  ENDIF
   IF(error > err_range1) info=5
 CASE(6)
   IF(aim_znu > 0) THEN
! zhankel2(znu,xx) is calculated using Eq. (42a) of Ref. (1).
     CALL bes_recur( znu,xx, 1,zbes1,info1)
     CALL bes_recur(-znu,xx,-1,zbes2,info2)
     info=MAX(info1,info2)
     z1=EXP(zunit*znupi)
     zans=2*(zbes1*z1-zbes2)/(z1*z1-1)
   ELSE
! zhankel2(znu,xx) is calculated using Eq. (42b) of Ref. (1).
     CALL bes_recur( znu,xx,0,zbes1,info1)
     CALL bes_recur(-znu,xx,1,zbes2,info2)
     info=MAX(info1,info2)
     zans=2*(zbes1-zbes2)/(1-EXP(-2*zunit*znupi))
   ENDIF
 CASE(7)
   IF(re_znu >= 0) THEN
     CALL han2_temme(znu,xx,zans,info)
     RETURN
   ENDIF
! When re_znu<0, zhankel2(znu,xx) is calculated using Eq. (42e) of Ref. (1).
   CALL han2_temme(znua,xx,zans,info)
   zans=zans*EXP(zunit*znupi)
 END SELECT
 END SUBROUTINE hankel2_boundary

 REAL(kp) FUNCTION amp(cyl,znu,xx, zans0)
! This FUNCTION amp is explained in the beginning of test 4.
 USE mod_test,ONLY: kp,abs2,bessel,neumann,hankel1,hankel2
 IMPLICIT NONE
 CHARACTER(LEN=3),INTENT(IN):: cyl
 COMPLEX(kp),INTENT(IN):: zans0,znu
 REAL(kp),INTENT(IN):: xx
 COMPLEX(kp):: zansm,zansp
 REAL(kp):: ansm,ans0,ansp
 INTEGER:: info
 IF(cyl == 'bes') THEN
   CALL bessel(znu-0.5,xx,zansm,info)
   CALL test_info('zbessel ',znu-0.5,xx,info)
   CALL bessel(znu+0.5,xx,zansp,info)
   CALL test_info('zbessel ',znu+0.5,xx,info)
 ELSE IF(cyl == 'neu') THEN
   CALL neumann(znu-0.5,xx,zansm,info)
   CALL test_info('zneumann',znu-0.5,xx,info)
   CALL neumann(znu+0.5,xx,zansp,info)
   CALL test_info('zneumann',znu+0.5,xx,info)
 ELSE IF(cyl == 'hn1') THEN
   CALL hankel1(znu-0.5,xx,zansm,info)
   CALL test_info('zhankel1',znu-0.5,xx,info)
   CALL hankel1(znu+0.5,xx,zansp,info)
   CALL test_info('zhankel1',znu+0.5,xx,info)
 ELSE IF(cyl == 'hn2') THEN
   CALL hankel2(znu-0.5,xx,zansm,info)
   CALL test_info('zhankel2',znu-0.5,xx,info)
   CALL hankel2(znu+0.5,xx,zansp,info)
   CALL test_info('zhankel2',znu+0.5,xx,info)
 ELSE
   STOP 'Input error(amp)'
 ENDIF
 ansm=abs2(zansm);   ans0=abs2(zans0);   ansp=abs2(zansp)
 IF(ans0<ansm .AND. ans0<ansp) THEN
   amp=MIN(ansm,ansp)
 ELSE
    amp=ans0
 ENDIF
 END FUNCTION amp

 SUBROUTINE test_info(zfunction,znu,xx,info)
! SUBROUTINE test_info examines the correctness of the value of info.
! If the value of info is not normal, this FUNCTION issues a warning to the 
! standard output unit in order to give us attention.
 USE mod_test,ONLY: kp
 IMPLICIT NONE
 CHARACTER(LEN=8),INTENT(IN):: zfunction
 COMPLEX(kp),INTENT(IN):: znu
 REAL(kp),INTENT(IN):: xx
 INTEGER,INTENT(IN):: info
 IF(info == 0) RETURN
 IF(-50<REAL(znu) .AND. REAL(znu)<80 .AND. ABS(AIMAG(znu))<65 .AND. &
      1E-19_kp<xx .AND. xx<1E3 .AND. (info==5 .OR. info==10 .OR. 22<info)) THEN
   WRITE(*,5) zfunction,znu,xx,info
 ENDIF
 IF(info ==  5) RETURN;   IF(info == 10) RETURN
 IF(info == 20) RETURN;   IF(info == 30) RETURN
 PRINT *,'info=',info;    STOP
5 FORMAT(A8,' (',F8.2,',',F8.2,')',F8.2,I5)
 END SUBROUTINE test_info
