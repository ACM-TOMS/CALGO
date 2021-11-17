  PROGRAM testmarq
! ------------------------------------------------------
! This is a test program for the Fortran 90 
! module MarcumQ. This module includes the 
! public routine marcum for the computation of the 
! Marcum-Q function Q(mu,a,x) and its complementary, 
! the function P(mu,x,y)  
! ------------------------------------------------------ 
  USE MarcumQ
  USE Someconstants
  IMPLICIT NONE
  REAL(r8) mu, x, y, pm1, p0, p1, p2, qm1, q0, q1, q2
  REAL(r8) err1, err2, d0, delta
  REAL(r8) mu1(18),x1(18),y1(18),p(18),q(18)
  INTEGER  ierr1, ierr2, ierr3, ierr4, i, j, k
  DATA mu1/1.0_r8,3.0_r8,4.0_r8,6.0_r8,8.0_r8,&
          10.0_r8,20.0_r8,22.0_r8,25.0_r8,27.0_r8,&
          30.0_r8,32.0_r8,40.0_r8,50.0_r8,200.0_r8,&
          350.0_r8,570.0_r8,1000.0_r8/                           
  DATA y1/0.01_r8,0.1_r8,50.0_r8,10.0_r8,15.0_r8,&
          25.0_r8,30.0_r8,150.0_r8,60.0_r8,205.0_r8,&
          90.0_r8,100.0_r8,120.0_r8,150.0_r8,190.0_r8,&
          320.0_r8,480.0_r8,799.0_r8/                             
  DATA x1/0.3_r8,2.0_r8,8.0_r8,25.0_r8,13.0_r8,45.0_r8,&
          47.0_r8,100.0_r8,85.0_r8,120.0_r8,130.0_r8,&
          140.0_r8,30.0_r8,40.0_r8,0.01_r8,100.0_r8,&
          1.0_r8,0.08_r8/                                
  DATA p/.7382308441994e-02_r8,.2199222796783e-04_r8,&
         .9999999768807_r8,.1746869995977e-03_r8,&
         .1483130042637_r8,.1748328323235e-03_r8,&
         .1340769184710e-04_r8,.9646591215441_r8,&
         .1783991673043e-04_r8,.9994542406431_r8,&
         .1220231636641e-05_r8,.1757487653307e-05_r8,&
         .9999894753719_r8,.9999968347378_r8,&
         .2431297758708_r8,.3851423018735e-09_r8,&
         .2984493152360e-04_r8,.4191472999694e-11_r8/
  DATA q/.9926176915580_r8,.9999780077720_r8,&
         .2311934913546e-07_r8,.9998253130004_r8,&
         .8516869957363_r8,.9998251671677_r8,&
         .9999865923082_r8,.3534087845586e-01_r8,&
         .9999821600833_r8,.5457593568564e-03_r8,&
         .9999987797684_r8,.9999982425123_r8,&
         .1052462813144e-04_r8,.3165262228904e-05_r8,&
         .7568702241292_r8,.9999999996149_r8,&
         .9999701550685_r8,.9999999999958_r8/
  open(unit=1,file='testQP.dat',status='unknown') 
! ----------------------------------------------------------------
! 1) We check some values of the functions Q(mu,x,y), P(mu,x,y)
! ----------------------------------------------------------------
  write(1,*)'ERR means relative error'
  write(1,*)'********************************************'
  write(1,*)'Test of the values of the functions'
  write(1,*)'********************************************'
  write(1,30)'mu','x','y','ERR(Q(mu,x,y))','ERR(P(mu,x,y))'
  DO i=1,18
    x=x1(i)
    y=y1(i)
    mu=mu1(i) 
    ierr1=0
    CALL marcum(mu,x,y, p0,q0,ierr1)
    err1=abs(1.0_r8-p0/p(i))
    err2=abs(1.0_r8-q0/q(i))
    write(1,31)int(mu),x,y,err1,err2
  ENDDO
! ----------------------------------------------------------------------
! 2) Test of the recurrence relation
!      a) If y > x+mu  
!          xQ(mu+2,x,y)=(x-mu)Q(mu+1,x,y)+(y+mu)Q(mu,x,y)-yQ(mu-1,x,y)
!      B) If y <= x+mu
!          xP(mu+2,x,y)=(x-mu)P(mu+1,x,y)+(y+mu)P(mu,x,y)-yP(mu-1,x,y)
! ----------------------------------------------------------------------
  d0=-1; 
  DO i=0,10
    mu=i*50.0_r8+10.0_r8;
    DO j=1,10 
      x=j*50.18_r8+5.0_r8;
      DO k=0,10
        y=k*19.15_r8+2.0_r8;
        CALL marcum(mu,x,y, p0,q0,ierr1)
        CALL marcum(mu-1,x,y,pm1,qm1,ierr2)
        CALL marcum(mu+1,x,y,p1,q1,ierr3)
        CALL marcum(mu+2,x,y,p2,q2,ierr4)
        IF (((ierr1==0).AND.(ierr2==0)).AND.((ierr3==0).AND.(ierr4==0))) THEN
          IF (y > x + mu) THEN
            delta=abs(((x-mu)*q1+(y+mu)*q0)/(x*q2+y*qm1)-1.0_r8)
          ELSE
            delta=abs(((x-mu)*p1+(y+mu)*p0)/(x*p2+y*pm1)-1.0_r8)
          ENDIF
          IF (delta>d0) THEN
            d0=delta
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO 
  write(1,*)' '
  write(1,*)'Maximum value of the recurrence check =',d0 
 30 format(3x,a2,3x,a1,6x,a1,8x,a15,2x,a15)
 31 format(i4,1x,d8.2,1x,d8.2,1x,d16.9,1x,d16.9)
  END PROGRAM testmarq

     


