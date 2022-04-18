    PROGRAM testm

!                             David M. Smith              3-23-97

!    Test program using the FM Fortran-90 module for doing
!    arithmetic using the FM, IM, and ZM derived types.

!    Any errors will be noted in file Test90.LOG.
!    After a successful run of this program, there should be
!    one line in Test90.LOG:
!    603 cases tested.  No errors were found.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Local Structures ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Local Scalars ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL test1, test10, test11, test12, test13, test14, test15, test16, &
        test17, test18, test19, test2, test3, test4, test5, test6, test7, &
        test8, test9, zmset
! ..
      CALL zmset(50)
      kdebug = 1
      kwarn = 2

!             Write output to the standard FM output (unit KW, defined
!             in subroutine FMSET), and also to the file Test90.LOG.

      klog = 11
      OPEN (klog,file='Test90.LOG')

!             NERROR is the number of errors found.
!             NCASE is the number of cases tested.

      nerror = 0
      ncase = 0

      i1 = 131
      r1 = 241.21
      d1 = 391.61D0
      z1 = (411.11D0,421.21D0)
      c1 = (431.11D0,441.21D0)
      CALL fm_st2m('581.21',mfm1)
      CALL fm_st2m('-572.42',mfm2)
      CALL im_st2m('661',mim1)
      CALL im_st2m('-602',mim2)
      CALL zm_st2m('731.51 + 711.41 i',mzm1)
      CALL zm_st2m('-762.12 - 792.42 i',mzm2)

!             Test the '=' assignment operator.

      CALL test1(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3,mzm4, &
        nerror,ncase,klog)

!             Test the '.EQ.' logical operator.

      CALL test2(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,mzm1,mzm2,nerror,ncase, &
        klog)

!             Test the '.NE.' logical operator.

      CALL test3(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,mzm1,mzm2,nerror,ncase, &
        klog)

!             Test the '.GT.' logical operator.

      CALL test4(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.GE.' logical operator.

      CALL test5(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.LT.' logical operator.

      CALL test6(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.LE.' logical operator.

      CALL test7(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '+' arithmetic operator.

      CALL test8(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '-' arithmetic operator.

      CALL test9(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '*' arithmetic operator.

      CALL test10(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '/' arithmetic operator.

      CALL test11(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '**' arithmetic operator.

      CALL test12(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test functions ABS, ..., CEILING.

      CALL test13(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3, &
        mzm4,nerror,ncase,klog)

!             Test functions CMPLX, ..., EXPONENT.

      CALL test14(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm3,mzm4,nerror,ncase,klog)

!             Test functions FLOOR, ..., MIN.

      CALL test15(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm3,mzm4,nerror,ncase,klog)

!             Test functions MINEXPONENT, ..., RRSPACING.

      CALL test16(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3,mzm4, &
        nerror,ncase,klog)

!             Test functions SCALE, ..., TINY.

      CALL test17(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4,mzm1, &
        mzm3,mzm4,nerror,ncase,klog)

!             Test functions TO_FM, ..., TO_ZM.

      CALL test18(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3,mzm4, &
        nerror,ncase,klog)

!             Test derived-type interface routines.

      CALL test19(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3,mzm4, &
        nerror,ncase,klog)

      IF (nerror==0) THEN
        WRITE (kw,*) ncase, ' cases tested.  No errors were found. '
        WRITE (klog,*) ncase, ' cases tested.  No errors were found. '
      ELSE
        WRITE (kw,*) ncase, ' cases tested.  ', nerror, ' error(s) found. '
        WRITE (klog,*) ncase, ' cases tested.  ', nerror, ' error(s) found. '
      END IF

    END PROGRAM testm
    SUBROUTINE test1(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3, &
        mzm4,nerror,ncase,klog)

!             Test the '=' assignment operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm3, mfm4
      TYPE (im) :: mim1, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (fm) :: msmall
! ..
! .. Local Scalars ..
      COMPLEX (kind(0.0D0)) :: c3
      COMPLEX :: z3
      REAL (KIND(0.0D0)) :: d3, dsmall
      REAL :: r3, rsmall
      INTEGER :: i3
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0
      msmall = epsilon(to_fm(1))*10000.0
      ncase = 1
      i3 = mfm1
      IF (i3/=581) CALL prterr(kw,klog,ncase,nerror)

      ncase = 2
      i3 = mim1
      IF (i3/=661) CALL prterr(kw,klog,ncase,nerror)

      ncase = 3
      i3 = mzm1
      IF (i3/=731) CALL prterr(kw,klog,ncase,nerror)

      ncase = 4
      r3 = mfm1
      IF (abs((r3-581.21)/581.21)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 5
      r3 = mim1
      IF (abs((r3-661.0)/661.0)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 6
      r3 = mzm1
      IF (abs((r3-731.51)/731.51)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 7
      d3 = mfm1
      IF (abs((d3-581.21D0)/581.21D0)>dsmall) CALL prterr(kw,klog,ncase,nerror &
        )

      ncase = 8
      d3 = mim1
      IF (abs((d3-661.0D0)/661.0D0)>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 9
      d3 = mzm1
      IF (abs((d3-731.51D0)/731.51D0)>dsmall) CALL prterr(kw,klog,ncase,nerror &
        )

      ncase = 10
      z3 = mfm1
      IF (abs((z3-581.21)/581.21)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 11
      z3 = mim1
      IF (abs((z3-661.0)/661.0)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 12
      z3 = mzm1
      IF (abs((z3-(731.51,711.41))/(731.51,711.41))>rsmall) CALL prterr(kw, &
        klog,ncase,nerror)

      ncase = 13
      c3 = mfm1
      IF (abs((c3-581.21D0)/581.21D0)>dsmall) CALL prterr(kw,klog,ncase,nerror &
        )

      ncase = 14
      c3 = mim1
      IF (abs((c3-661.0D0)/661.0D0)>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 15
      c3 = mzm1
      IF (abs((c3-(731.51D0,711.41D0))/(731.51D0,711.41D0))>dsmall) CALL &
        prterr(kw,klog,ncase,nerror)

      ncase = 16
      mfm3 = i1
      CALL fm_i2m(131,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_abs(mfm4,mfm4)
      CALL fm_st2m('0',mfm3)
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 17
      mfm3 = r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 18
      mfm3 = d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 19
      mfm3 = z1
      CALL fm_st2m('411.11',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 20
      mfm3 = c1
      CALL fm_st2m('431.11',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 21
      mfm3 = mfm1
      CALL fm_st2m('581.21',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      CALL fm_eq(msmall,mfm3)
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 22
      mfm3 = mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_abs(mfm4,mfm4)
      CALL fm_st2m('0',mfm3)
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 23
      mfm3 = mzm1
      CALL fm_st2m('731.51',mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = msmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 24
      mim3 = i1
      CALL im_i2m(131,mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 25
      mim3 = r1
      CALL im_st2m('241',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 26
      mim3 = d1
      CALL im_st2m('391',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 27
      mim3 = z1
      CALL im_st2m('411',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 28
      mim3 = c1
      CALL im_st2m('431',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 29
      mim3 = mfm1
      CALL im_st2m('581',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 30
      mim3 = mim1
      CALL im_st2m('661',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 31
      mim3 = mzm1
      CALL im_st2m('731',mim4)
      CALL im_sub(mim3,mim4,mim4)
      CALL im_st2m('0',mim3)
      IF (im_comp(mim4,'GT',mim3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 32
      mzm3 = i1
      CALL zm_i2m(131,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_abs(mzm4,mfm4)
      CALL fm_st2m('0',mfm3)
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 33
      mzm3 = r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 34
      mzm3 = d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 35
      mzm3 = z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 36
      mzm3 = c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 37
      mzm3 = mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = msmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 38
      mzm3 = mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_abs(mzm4,mfm4)
      CALL fm_st2m('0',mfm3)
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 39
      mzm3 = mzm1
      CALL zm_st2m('731.51 + 711.41 i',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = msmall
      IF (fm_comp(mfm4,'GT',mfm3)) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test1

    SUBROUTINE test2(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,mzm1,mzm2,nerror, &
        ncase,klog)

!             Test the '.EQ.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
      TYPE (zm) :: mzm1, mzm2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 40
      IF (i1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 41
      IF (i1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 42
      IF (i1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 43
      IF (r1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 44
      IF (r1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 45
      IF (r1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 46
      IF (d1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 47
      IF (d1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 48
      IF (d1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 49
      IF (z1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 50
      IF (z1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 51
      IF (z1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 52
      IF (c1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 53
      IF (c1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 54
      IF (c1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 55
      IF (mfm1==i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 56
      IF (mfm1==r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 57
      IF (mfm1==d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 58
      IF (mfm1==z1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 59
      IF (mfm1==c1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 60
      IF (mfm1==mfm2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 61
      IF (mfm1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 62
      IF (mfm1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 63
      IF (mim1==i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 64
      IF (mim1==r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 65
      IF (mim1==d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 66
      IF (mim1==z1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 67
      IF (mim1==c1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 68
      IF (mim1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 69
      IF (mim1==mim2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 70
      IF (mim1==mzm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 71
      IF (mzm1==i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 72
      IF (mzm1==r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 73
      IF (mzm1==d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 74
      IF (mzm1==z1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 75
      IF (mzm1==c1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 76
      IF (mzm1==mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 77
      IF (mzm1==mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 78
      IF (mzm1==mzm2) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test2

    SUBROUTINE test3(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,mzm1,mzm2,nerror, &
        ncase,klog)

!             Test the '.NE.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
      TYPE (zm) :: mzm1, mzm2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 79
      IF ( .NOT. (i1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 80
      IF ( .NOT. (i1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 81
      IF ( .NOT. (i1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 82
      IF ( .NOT. (r1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 83
      IF ( .NOT. (r1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 84
      IF ( .NOT. (r1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 85
      IF ( .NOT. (d1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 86
      IF ( .NOT. (d1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 87
      IF ( .NOT. (d1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 88
      IF ( .NOT. (z1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 89
      IF ( .NOT. (z1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 90
      IF ( .NOT. (z1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 91
      IF ( .NOT. (c1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 92
      IF ( .NOT. (c1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 93
      IF ( .NOT. (c1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 94
      IF ( .NOT. (mfm1/=i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 95
      IF ( .NOT. (mfm1/=r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 96
      IF ( .NOT. (mfm1/=d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 97
      IF ( .NOT. (mfm1/=z1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 98
      IF ( .NOT. (mfm1/=c1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 99
      IF ( .NOT. (mfm1/=mfm2)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 100
      IF ( .NOT. (mfm1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 101
      IF ( .NOT. (mfm1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 102
      IF ( .NOT. (mim1/=i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 103
      IF ( .NOT. (mim1/=r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 104
      IF ( .NOT. (mim1/=d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 105
      IF ( .NOT. (mim1/=z1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 106
      IF ( .NOT. (mim1/=c1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 107
      IF ( .NOT. (mim1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 108
      IF ( .NOT. (mim1/=mim2)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 109
      IF ( .NOT. (mim1/=mzm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 110
      IF ( .NOT. (mzm1/=i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 111
      IF ( .NOT. (mzm1/=r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 112
      IF ( .NOT. (mzm1/=d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 113
      IF ( .NOT. (mzm1/=z1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 114
      IF ( .NOT. (mzm1/=c1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 115
      IF ( .NOT. (mzm1/=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 116
      IF ( .NOT. (mzm1/=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 117
      IF ( .NOT. (mzm1/=mzm2)) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test3

    SUBROUTINE test4(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.GT.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 118
      IF (i1>mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 119
      IF (i1>mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 120
      IF (r1>mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 121
      IF (r1>mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 122
      IF (d1>mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 123
      IF (d1>mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 124
      IF ( .NOT. (mfm1>i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 125
      IF ( .NOT. (mfm1>r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 126
      IF ( .NOT. (mfm1>d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 127
      IF ( .NOT. (mfm1>mfm2)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 128
      IF (mfm1>mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 129
      IF ( .NOT. (mim1>i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 130
      IF ( .NOT. (mim1>r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 131
      IF ( .NOT. (mim1>d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 132
      IF ( .NOT. (mim1>mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 133
      IF ( .NOT. (mim1>mim2)) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test4

    SUBROUTINE test5(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.GE.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 134
      IF (i1>=mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 135
      IF (i1>=mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 136
      IF (r1>=mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 137
      IF (r1>=mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 138
      IF (d1>=mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 139
      IF (d1>=mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 140
      IF ( .NOT. (mfm1>=i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 141
      IF ( .NOT. (mfm1>=r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 142
      IF ( .NOT. (mfm1>=d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 143
      IF ( .NOT. (mfm1>=mfm2)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 144
      IF (mfm1>=mim1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 145
      IF ( .NOT. (mim1>=i1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 146
      IF ( .NOT. (mim1>=r1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 147
      IF ( .NOT. (mim1>=d1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 148
      IF ( .NOT. (mim1>=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 149
      IF ( .NOT. (mim1>=mim2)) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test5

    SUBROUTINE test6(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.LT.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 150
      IF ( .NOT. (i1<mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 151
      IF ( .NOT. (i1<mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 152
      IF ( .NOT. (r1<mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 153
      IF ( .NOT. (r1<mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 154
      IF ( .NOT. (d1<mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 155
      IF ( .NOT. (d1<mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 156
      IF (mfm1<i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 157
      IF (mfm1<r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 158
      IF (mfm1<d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 159
      IF (mfm1<mfm2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 160
      IF ( .NOT. (mfm1<mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 161
      IF (mim1<i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 162
      IF (mim1<r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 163
      IF (mim1<d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 164
      IF (mim1<mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 165
      IF (mim1<mim2) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test6

    SUBROUTINE test7(i1,r1,d1,z1,c1,mfm1,mfm2,mim1,mim2,nerror,ncase,klog)

!             Test the '.LE.' logical operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2
      TYPE (im) :: mim1, mim2
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 166
      IF ( .NOT. (i1<=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 167
      IF ( .NOT. (i1<=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 168
      IF ( .NOT. (r1<=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 169
      IF ( .NOT. (r1<=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 170
      IF ( .NOT. (d1<=mfm1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 171
      IF ( .NOT. (d1<=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 172
      IF (mfm1<=i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 173
      IF (mfm1<=r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 174
      IF (mfm1<=d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 175
      IF (mfm1<=mfm2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 176
      IF ( .NOT. (mfm1<=mim1)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 177
      IF (mim1<=i1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 178
      IF (mim1<=r1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 179
      IF (mim1<=d1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 180
      IF (mim1<=mfm1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 181
      IF (mim1<=mim2) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test7

    SUBROUTINE test8(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '+' arithmetic operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dsmall
      REAL :: rsmall
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0

      ncase = 182
      mfm3 = i1 + mfm1
      CALL fm_st2m('131',mfm4)
      CALL fm_add(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 183
      mim3 = i1 + mim1
      CALL im_st2m('131',mim4)
      CALL im_add(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 184
      mzm3 = i1 + mzm1
      CALL zm_st2m('131',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 185
      mfm3 = r1 + mfm1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_add(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 186
      CALL fm_st2m('241.21',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_add(mfm4,mfm3,mfm4)
      mfm3 = r1 + mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 187
      mzm3 = r1 + mzm1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 188
      mfm3 = d1 + mfm1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_add(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 189
      CALL fm_st2m('391.61',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_add(mfm4,mfm3,mfm4)
      mfm3 = d1 + mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 190
      mzm3 = d1 + mzm1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 191
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = z1 + mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 192
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = z1 + mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 193
      mzm3 = z1 + mzm1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 194
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = c1 + mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 195
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = c1 + mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 196
      mzm3 = c1 + mzm1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 197
      mfm3 = mfm1 + i1
      CALL fm_st2m('131',mfm4)
      CALL fm_add(mfm1,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 198
      mfm3 = mfm1 + r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_add(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 199
      mfm3 = mfm1 + d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_add(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 200
      CALL zm_st2m('581.21',mzm3)
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_add(mzm3,mzm4,mzm4)
      mzm3 = mfm1 + z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 201
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('581.21',mzm4)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = mfm1 + c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 202
      mfm3 = mfm1 + mfm2
      CALL fm_add(mfm1,mfm2,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 203
      mfm3 = mfm1 + mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_add(mfm1,mfm4,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 204
      mzm3 = mfm1 + mzm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 205
      mim3 = mim1 + i1
      CALL im_st2m('131',mim4)
      CALL im_add(mim1,mim4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 206
      CALL fm_st2m('241.21',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_add(mfm4,mfm3,mfm4)
      mfm3 = mim1 + r1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 207
      CALL fm_st2m('391.61',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_add(mfm4,mfm3,mfm4)
      mfm3 = mim1 + d1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 208
      CALL zm_st2m('411.11 + 421.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = mim1 + z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 209
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_add(mzm4,mzm3,mzm4)
      mzm3 = mim1 + c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 210
      mfm3 = mim1 + mfm1
      CALL fm_st2m('661',mfm4)
      CALL fm_add(mfm4,mfm1,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 211
      mim3 = mim1 + mim2
      CALL im_add(mim1,mim2,mim4)
      IF (mim4/=mim3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 212
      mzm3 = mim1 + mzm1
      CALL zm_st2m('661',mzm4)
      CALL zm_add(mzm4,mzm1,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 213
      mzm3 = mzm1 + i1
      CALL zm_st2m('131',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 214
      mzm3 = mzm1 + r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 215
      mzm3 = mzm1 + d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 216
      mzm3 = mzm1 + z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 217
      mzm3 = mzm1 + c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 218
      mzm3 = mzm1 + mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 219
      mzm3 = mzm1 + mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_add(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 220
      mzm3 = mzm1 + mzm2
      CALL zm_add(mzm1,mzm2,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 221
      mfm3 = + mfm1
      CALL fm_eq(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 222
      mim3 = + mim1
      CALL im_eq(mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 223
      mzm3 = + mzm1
      CALL zm_eq(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test8

    SUBROUTINE test9(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '-' arithmetic operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dsmall
      REAL :: rsmall
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0

      ncase = 224
      mfm3 = i1 - mfm1
      CALL fm_st2m('131',mfm4)
      CALL fm_sub(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 225
      mim3 = i1 - mim1
      CALL im_st2m('131',mim4)
      CALL im_sub(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 226
      mzm3 = i1 - mzm1
      CALL zm_st2m('131',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 227
      mfm3 = r1 - mfm1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_sub(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 228
      CALL fm_st2m('241.21',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_sub(mfm4,mfm3,mfm4)
      mfm3 = r1 - mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 229
      mzm3 = r1 - mzm1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 230
      mfm3 = d1 - mfm1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_sub(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 231
      CALL fm_st2m('391.61',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_sub(mfm4,mfm3,mfm4)
      mfm3 = d1 - mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 232
      mzm3 = d1 - mzm1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 233
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = z1 - mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 234
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = z1 - mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 235
      mzm3 = z1 - mzm1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 236
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = c1 - mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 237
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = c1 - mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 238
      mzm3 = c1 - mzm1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 239
      mfm3 = mfm1 - i1
      CALL fm_st2m('131',mfm4)
      CALL fm_sub(mfm1,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 240
      mfm3 = mfm1 - r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_sub(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 241
      mfm3 = mfm1 - d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_sub(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 242
      CALL zm_st2m('581.21',mzm3)
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      mzm3 = mfm1 - z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 243
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('581.21',mzm4)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = mfm1 - c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 244
      mfm3 = mfm1 - mfm2
      CALL fm_sub(mfm1,mfm2,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 245
      mfm3 = mfm1 - mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_sub(mfm1,mfm4,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 246
      mzm3 = mfm1 - mzm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 247
      mim3 = mim1 - i1
      CALL im_st2m('131',mim4)
      CALL im_sub(mim1,mim4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 248
      CALL fm_st2m('241.21',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_sub(mfm4,mfm3,mfm4)
      mfm3 = mim1 - r1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 249
      CALL fm_st2m('391.61',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_sub(mfm4,mfm3,mfm4)
      mfm3 = mim1 - d1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 250
      CALL zm_st2m('411.11 + 421.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = mim1 - z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 251
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_sub(mzm4,mzm3,mzm4)
      mzm3 = mim1 - c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 252
      mfm3 = mim1 - mfm1
      CALL fm_st2m('661',mfm4)
      CALL fm_sub(mfm4,mfm1,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 253
      mim3 = mim1 - mim2
      CALL im_sub(mim1,mim2,mim4)
      IF (mim4/=mim3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 254
      mzm3 = mim1 - mzm1
      CALL zm_st2m('661',mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 255
      mzm3 = mzm1 - i1
      CALL zm_st2m('131',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 256
      mzm3 = mzm1 - r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 257
      mzm3 = mzm1 - d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 258
      mzm3 = mzm1 - z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 259
      mzm3 = mzm1 - c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 260
      mzm3 = mzm1 - mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 261
      mzm3 = mzm1 - mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_sub(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 262
      mzm3 = mzm1 - mzm2
      CALL zm_sub(mzm1,mzm2,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 263
      mfm3 = -mfm1
      CALL fm_i2m(0,mfm4)
      CALL fm_sub(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 264
      mim3 = -mim1
      CALL im_i2m(0,mim4)
      CALL im_sub(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 265
      mzm3 = -mzm1
      CALL zm_i2m(0,mzm4)
      CALL zm_sub(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test9

    SUBROUTINE test10(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '*' arithmetic operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dsmall
      REAL :: rsmall
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0

      ncase = 266
      mfm3 = i1*mfm1
      CALL fm_st2m('131',mfm4)
      CALL fm_mpy(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 267
      mim3 = i1*mim1
      CALL im_st2m('131',mim4)
      CALL im_mpy(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 268
      mzm3 = i1*mzm1
      CALL zm_st2m('131',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 269
      mfm3 = r1*mfm1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_mpy(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 270
      CALL fm_st2m('241.21',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_mpy(mfm4,mfm3,mfm4)
      mfm3 = r1*mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 271
      mzm3 = r1*mzm1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 272
      mfm3 = d1*mfm1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_mpy(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 273
      CALL fm_st2m('391.61',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_mpy(mfm4,mfm3,mfm4)
      mfm3 = d1*mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 274
      mzm3 = d1*mzm1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 275
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = z1*mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 276
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = z1*mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 277
      mzm3 = z1*mzm1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 278
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = c1*mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 279
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = c1*mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 280
      mzm3 = c1*mzm1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 281
      mfm3 = mfm1*i1
      CALL fm_st2m('131',mfm4)
      CALL fm_mpy(mfm1,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 282
      mfm3 = mfm1*r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_mpy(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 283
      mfm3 = mfm1*d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_mpy(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 284
      CALL zm_st2m('581.21',mzm3)
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_mpy(mzm3,mzm4,mzm4)
      mzm3 = mfm1*z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 285
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('581.21',mzm4)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = mfm1*c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 286
      mfm3 = mfm1*mfm2
      CALL fm_mpy(mfm1,mfm2,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 287
      mfm3 = mfm1*mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_mpy(mfm1,mfm4,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 288
      mzm3 = mfm1*mzm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 289
      mim3 = mim1*i1
      CALL im_st2m('131',mim4)
      CALL im_mpy(mim1,mim4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 290
      CALL fm_st2m('241.21',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_mpy(mfm4,mfm3,mfm4)
      mfm3 = mim1*r1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 291
      CALL fm_st2m('391.61',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_mpy(mfm4,mfm3,mfm4)
      mfm3 = mim1*d1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 292
      CALL zm_st2m('411.11 + 421.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = mim1*z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 293
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_mpy(mzm4,mzm3,mzm4)
      mzm3 = mim1*c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 294
      mfm3 = mim1*mfm1
      CALL fm_st2m('661',mfm4)
      CALL fm_mpy(mfm4,mfm1,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 295
      mim3 = mim1*mim2
      CALL im_mpy(mim1,mim2,mim4)
      IF (mim4/=mim3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 296
      mzm3 = mim1*mzm1
      CALL zm_st2m('661',mzm4)
      CALL zm_mpy(mzm4,mzm1,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 297
      mzm3 = mzm1*i1
      CALL zm_st2m('131',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 298
      mzm3 = mzm1*r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 299
      mzm3 = mzm1*d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 300
      mzm3 = mzm1*z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 301
      mzm3 = mzm1*c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 302
      mzm3 = mzm1*mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 303
      mzm3 = mzm1*mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_mpy(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 304
      mzm3 = mzm1*mzm2
      CALL zm_mpy(mzm1,mzm2,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test10

    SUBROUTINE test11(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '/' arithmetic operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dsmall
      REAL :: rsmall
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0

      ncase = 305
      mfm3 = i1/mfm1
      CALL fm_st2m('131',mfm4)
      CALL fm_div(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 306
      mim3 = i1/mim1
      CALL im_st2m('131',mim4)
      CALL im_div(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 307
      mzm3 = i1/mzm1
      CALL zm_st2m('131',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 308
      mfm3 = r1/mfm1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_div(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 309
      CALL fm_st2m('241.21',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_div(mfm4,mfm3,mfm4)
      mfm3 = r1/mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 310
      mzm3 = r1/mzm1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 311
      mfm3 = d1/mfm1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_div(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 312
      CALL fm_st2m('391.61',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_div(mfm4,mfm3,mfm4)
      mfm3 = d1/mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 313
      mzm3 = d1/mzm1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 314
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = z1/mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 315
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = z1/mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 316
      mzm3 = z1/mzm1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 317
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = c1/mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 318
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = c1/mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 319
      mzm3 = c1/mzm1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 320
      mfm3 = mfm1/i1
      CALL fm_st2m('131',mfm4)
      CALL fm_div(mfm1,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 321
      mfm3 = mfm1/r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_div(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 322
      mfm3 = mfm1/d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_div(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 323
      CALL zm_st2m('581.21',mzm3)
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_div(mzm3,mzm4,mzm4)
      mzm3 = mfm1/z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 324
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('581.21',mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = mfm1/c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 325
      mfm3 = mfm1/mfm2
      CALL fm_div(mfm1,mfm2,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 326
      mfm3 = mfm1/mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_div(mfm1,mfm4,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 327
      mzm3 = mfm1/mzm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 328
      mim3 = mim1/i1
      CALL im_st2m('131',mim4)
      CALL im_div(mim1,mim4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 329
      CALL fm_st2m('241.21',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      mfm3 = mim1/r1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 330
      CALL fm_st2m('391.61',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      mfm3 = mim1/d1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 331
      CALL zm_st2m('411.11 + 421.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = mim1/z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 332
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      mzm3 = mim1/c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 333
      mfm3 = mim1/mfm1
      CALL fm_st2m('661',mfm4)
      CALL fm_div(mfm4,mfm1,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 334
      mim3 = mim1/mim2
      CALL im_div(mim1,mim2,mim4)
      IF (mim4/=mim3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 335
      mzm3 = mim1/mzm1
      CALL zm_st2m('661',mzm4)
      CALL zm_div(mzm4,mzm1,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 336
      mzm3 = mzm1/i1
      CALL zm_st2m('131',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 337
      mzm3 = mzm1/r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 338
      mzm3 = mzm1/d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 339
      mzm3 = mzm1/z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 340
      mzm3 = mzm1/c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 341
      mzm3 = mzm1/mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 342
      mzm3 = mzm1/mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_div(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 343
      mzm3 = mzm1/mzm2
      CALL zm_div(mzm1,mzm2,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test11

    SUBROUTINE test12(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm2,mzm3,mzm4,nerror,ncase,klog)

!             Test the '**' arithmetic operator.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm2, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dsmall
      REAL :: rsmall
      INTEGER :: i3
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
!             Use a larger error tolerance for large exponents.

      rsmall = epsilon(1.0)*10000.0
      dsmall = epsilon(1.0D0)*10000.0

      ncase = 344
      mfm3 = i1**mfm1
      CALL fm_st2m('131',mfm4)
      CALL fm_pwr(mfm4,mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 345
      i3 = 13
      mim3 = i3**mim1
      CALL im_st2m('13',mim4)
      CALL im_pwr(mim4,mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 346
      mzm3 = i1**mzm1
      CALL zm_st2m('131',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 347
      mfm3 = r1**mfm1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_pwr(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)


      ncase = 348
      CALL fm_st2m('241.21',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_pwr(mfm4,mfm3,mfm4)
      mfm3 = r1**mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 349
      mzm3 = r1**mzm1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 350
      mfm3 = d1**mfm1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_pwr(mfm4,mfm1,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 351
      CALL fm_st2m('391.61',mfm4)
      CALL fm_st2m('661',mfm3)
      CALL fm_pwr(mfm4,mfm3,mfm4)
      mfm3 = d1**mim1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 352
      mzm3 = d1**mzm1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 353
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = z1**mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 354
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = z1**mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 355
      mzm3 = z1**mzm1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 356
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('581.21',mzm3)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = c1**mfm1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 357
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_st2m('661',mzm3)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = c1**mim1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 358
      mzm3 = c1**mzm1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 359
      mfm3 = mfm1**i1
      CALL fm_st2m('131',mfm4)
      CALL fm_pwr(mfm1,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 360
      mfm3 = mfm1**r1
      CALL fm_st2m('241.21',mfm4)
      CALL fm_pwr(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 361
      mfm3 = mfm1**d1
      CALL fm_st2m('391.61',mfm4)
      CALL fm_pwr(mfm1,mfm4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 362
      CALL zm_st2m('581.21',mzm3)
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_pwr(mzm3,mzm4,mzm4)
      mzm3 = mfm1**z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 363
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('581.21',mzm4)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = mfm1**c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 364
      mfm3 = mfm1**mfm2
      CALL fm_pwr(mfm1,mfm2,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 365
      mfm3 = mfm1**mim1
      CALL fm_st2m('661',mfm4)
      CALL fm_pwr(mfm1,mfm4,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 366
      mzm3 = mfm1**mzm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 367
      i3 = 17
      mim3 = mim1**i3
      CALL im_st2m('17',mim4)
      CALL im_pwr(mim1,mim4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 368
      CALL fm_st2m('241.21',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_pwr(mfm4,mfm3,mfm4)
      mfm3 = mim1**r1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 369
      CALL fm_st2m('391.61',mfm3)
      CALL fm_st2m('661',mfm4)
      CALL fm_pwr(mfm4,mfm3,mfm4)
      mfm3 = mim1**d1
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 370
      CALL zm_st2m('411.11 + 421.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = mim1**z1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 371
      CALL zm_st2m('431.11 + 441.21 i',mzm3)
      CALL zm_st2m('661',mzm4)
      CALL zm_pwr(mzm4,mzm3,mzm4)
      mzm3 = mim1**c1
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 372
      mfm3 = mim1**mfm1
      CALL fm_st2m('661',mfm4)
      CALL fm_pwr(mfm4,mfm1,mfm4)
      IF (mfm4/=mfm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 373
      mim4 = 19
      mim3 = mim1**mim4
      CALL im_pwr(mim1,mim4,mim4)
      IF (mim4/=mim3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 374
      mzm3 = mim1**mzm1
      CALL zm_st2m('661',mzm4)
      CALL zm_pwr(mzm4,mzm1,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 375
      mzm3 = mzm1**i1
      CALL zm_st2m('131',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 376
      mzm3 = mzm1**r1
      CALL zm_st2m('241.21',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 377
      mzm3 = mzm1**d1
      CALL zm_st2m('391.61',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 378
      mzm3 = mzm1**z1
      CALL zm_st2m('411.11 + 421.21 i',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 379
      mzm3 = mzm1**c1
      CALL zm_st2m('431.11 + 441.21 i',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      IF (mfm4>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 380
      mzm3 = mzm1**mfm1
      CALL zm_st2m('581.21',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 381
      mzm3 = mzm1**mim1
      CALL zm_st2m('661',mzm4)
      CALL zm_pwr(mzm1,mzm4,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

      ncase = 382
      mzm3 = mzm1**mzm2
      CALL zm_pwr(mzm1,mzm2,mzm4)
      IF (mzm4/=mzm3) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test12

    SUBROUTINE test13(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim3,mim4,mzm1, &
        mzm3,mzm4,nerror,ncase,klog)

!             Test functions ABS, ..., CEILING.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Scalars ..
      INTEGER :: j, jerr
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 383
      mfm3 = abs(mfm1)
      CALL fm_abs(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 384
      mim3 = abs(mim1)
      CALL im_abs(mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 385
      mfm3 = abs(mzm1)
      CALL zm_abs(mzm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 386
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = acos(mfm4)
      CALL fm_acos(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 387
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = acos(mzm4)
      CALL zm_acos(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 388
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mfm3 = aimag(mzm4)
      CALL zm_imag(mzm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 389
      mfm3 = aint(mfm1)
      CALL fm_int(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 390
      mzm3 = aint(mzm1)
      CALL zm_int(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 391
      mfm3 = anint(mfm1)
      CALL fm_nint(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 392
      mzm3 = anint(mzm1)
      CALL zm_nint(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 393
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = asin(mfm4)
      CALL fm_asin(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 394
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = asin(mzm4)
      CALL zm_asin(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 395
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = atan(mfm4)
      CALL fm_atan(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 396
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = atan(mzm4)
      CALL zm_atan(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 397
      mfm3 = atan2(mfm1,mfm2)
      CALL fm_atn2(mfm1,mfm2,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 398
      jerr = -1
      DO j = 0, 10
        IF (btest(661,j)) THEN
          IF ( .NOT. btest(mim1,j)) jerr = j
        ELSE
          IF (btest(mim1,j)) jerr = j
        END IF
      END DO
      IF (jerr>=0) CALL prterr(kw,klog,ncase,nerror)

      ncase = 399
      CALL fm_st2m('12.37654',mfm4)
      mfm3 = ceiling(mfm4)
      CALL fm_st2m('13',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 400
      CALL fm_st2m('-12.7654',mfm4)
      mfm3 = ceiling(mfm4)
      CALL fm_st2m('-12',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 401
      CALL zm_st2m('12.37654 - 22.54 i',mzm4)
      mzm3 = ceiling(mzm4)
      CALL zm_st2m('13 - 22 i',mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 402
      CALL zm_st2m('-12.7654 + 22.31 i',mzm4)
      mzm3 = ceiling(mzm4)
      CALL zm_st2m('-12 + 23 i',mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test13

    SUBROUTINE test14(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm3,mzm4,nerror,ncase,klog)

!             Test functions CMPLX, ..., EXPONENT.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (fm) :: mfmv1(3), mfmv2(3)
      TYPE (im) :: mimv1(3), mimv2(3)
      TYPE (zm) :: mzmv1(3), mzmv2(3)
! ..
! .. Local Scalars ..
      INTEGER :: j
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 403
      mzm3 = cmplx(mfm1,mfm2)
      CALL zm_cmpx(mfm1,mfm2,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 404
      mzm3 = cmplx(mim1,mim2)
      CALL im_i2fm(mim1,mfm3)
      CALL im_i2fm(mim2,mfm4)
      CALL zm_cmpx(mfm3,mfm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 405
      mzm3 = cmplx(mfm1)
      CALL fm_i2m(0,mfm4)
      CALL zm_cmpx(mfm1,mfm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 406
      mzm3 = cmplx(mim1)
      CALL im_i2fm(mim1,mfm3)
      CALL fm_i2m(0,mfm4)
      CALL zm_cmpx(mfm3,mfm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 407
      mzm3 = conjg(mzm1)
      CALL zm_conj(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 408
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = cos(mfm4)
      CALL fm_cos(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 409
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = cos(mzm4)
      CALL zm_cos(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 410
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = cosh(mfm4)
      CALL fm_cosh(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 411
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = cosh(mzm4)
      CALL zm_cosh(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 412
      mfm3 = dble(mfm1)
      CALL fm_eq(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 413
      mfm3 = dble(mim1)
      CALL im_i2fm(mim1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 414
      mfm3 = dble(mzm1)
      CALL zm_real(mzm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 415
      j = digits(mfm1)
      IF (j/=ndig) CALL prterr(kw,klog,ncase,nerror)

      ncase = 416
      j = digits(mim1)
      IF (j/=ndigmx) CALL prterr(kw,klog,ncase,nerror)

      ncase = 417
      j = digits(mzm1)
      IF (j/=ndig) CALL prterr(kw,klog,ncase,nerror)

      ncase = 418
      mfm3 = dim(mfm1,mfm2)
      CALL fm_dim(mfm1,mfm2,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 419
      mim3 = dim(mim1,mim2)
      CALL im_dim(mim1,mim2,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 420
      mfm3 = dint (mfm1)
      CALL fm_int(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 421
      mzm3 = dint (mzm1)
      CALL zm_int(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 422
      CALL fm_st2m('1.23',mfmv1(1))
      CALL fm_st2m('2.23',mfmv1(2))
      CALL fm_st2m('3.23',mfmv1(3))
      CALL fm_st2m('4.23',mfmv2(1))
      CALL fm_st2m('5.23',mfmv2(2))
      CALL fm_st2m('6.23',mfmv2(3))
      mfm3 = dotproduct(mfmv1,mfmv2)
      mfm4 = 0
      DO j = 1, 3
        mfm4 = mfm4 + mfmv1(j)*mfmv2(j)
      END DO
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 423
      CALL im_st2m('12',mimv1(1))
      CALL im_st2m('23',mimv1(2))
      CALL im_st2m('34',mimv1(3))
      CALL im_st2m('-14',mimv2(1))
      CALL im_st2m('-5',mimv2(2))
      CALL im_st2m('16',mimv2(3))
      mim3 = dotproduct(mimv1,mimv2)
      mim4 = 0
      DO j = 1, 3
        mim4 = mim4 + mimv1(j)*mimv2(j)
      END DO
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 424
      CALL zm_st2m('1.23 + 1.67 i',mzmv1(1))
      CALL zm_st2m('2.23 - 2.56 i',mzmv1(2))
      CALL zm_st2m('3.23 + 3.45 i',mzmv1(3))
      CALL zm_st2m('4.23 - 4.34 i',mzmv2(1))
      CALL zm_st2m('5.23 + 5.23 i',mzmv2(2))
      CALL zm_st2m('6.23 - 6.12 i',mzmv2(3))
      mzm3 = dotproduct(mzmv1,mzmv2)
      mzm4 = 0
      DO j = 1, 3
        mzm4 = mzm4 + mzmv1(j)*mzmv2(j)
      END DO
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 425
      mfm3 = epsilon(mfm1)
      CALL fm_i2m(1,mfm4)
      CALL fm_ulp(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 426
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = exp(mfm4)
      CALL fm_exp(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 427
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = exp(mzm4)
      CALL zm_exp(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 428
      j = exponent(mfm1)
      IF (j/=int(mfm1%mfm(1))) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test14

    SUBROUTINE test15(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm3,mzm4,nerror,ncase,klog)

!             Test functions FLOOR, ..., MIN.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (fm) :: mfma(3,3), mfmb(3,3), mfmc(3,3)
      TYPE (im) :: mima(2,2), mimb(2,2), mimc(2,2)
      TYPE (zm) :: mzma(2,3), mzmb(3,4), mzmc(2,4)
! ..
! .. Local Scalars ..
      INTEGER :: i, j
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 429
      CALL fm_st2m('12.37654',mfm4)
      mfm3 = floor(mfm4)
      CALL fm_st2m('12',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 430
      CALL fm_st2m('-12.7654',mfm4)
      mfm3 = floor(mfm4)
      CALL fm_st2m('-13',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 431
      CALL im_st2m('12',mim4)
      mim3 = floor(mim4)
      CALL im_st2m('12',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 432
      CALL im_st2m('-123',mim4)
      mim3 = floor(mim4)
      CALL im_st2m('-123',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 433
      CALL zm_st2m('12.37654 - 22.54 i',mzm4)
      mzm3 = floor(mzm4)
      CALL zm_st2m('12 - 23 i',mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 434
      CALL zm_st2m('-12.7654 + 22.31 i',mzm4)
      mzm3 = floor(mzm4)
      CALL zm_st2m('-13 + 22 i',mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 435
      CALL fm_st2m('12.37654',mfm4)
      mfm3 = fraction(mfm4)
      mfm4%mfm(1) = 0
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 436
      CALL zm_st2m('12.37654 - 22.54',mzm4)
      mzm3 = fraction(mzm4)
      mzm4%mzm(1) = 0
      mzm4%mzm(kptimu+1) = 0
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 437
      mfm3 = huge(mfm1)
      CALL fm_big(mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 438
      mim3 = huge(mim1)
      CALL im_big(mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 439
      mzm3 = huge(mzm1)
      CALL fm_big(mfm4)
      CALL zm_cmpx(mfm4,mfm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 440
      mim3 = int(mfm1)
      CALL fm_int(mfm1,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 441
      mim3 = int(mim1)
      CALL im_eq(mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 442
      mim3 = int(mzm1)
      CALL zm_int(mzm1,mzm4)
      CALL zm_real(mzm4,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 443
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = log(mfm4)
      CALL fm_ln(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 444
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = log(mzm4)
      CALL zm_ln(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 445
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = log10(mfm4)
      CALL fm_lg10(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 446
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = log10(mzm4)
      CALL zm_lg10(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 447
      DO i = 1, 3
        DO j = 1, 3
          mfma(i,j) = 3*(j-1) + i
          mfmb(i,j) = 3*(i-1) + j + 10
        END DO
      END DO
      mfmc = matmul(mfma,mfmb)
      mfm3 = abs(mfmc(1,1)-186) + abs(mfmc(1,2)-198) + abs(mfmc(1,3)-210) + &
        abs(mfmc(2,1)-228) + abs(mfmc(2,2)-243) + abs(mfmc(2,3)-258) + &
        abs(mfmc(3,1)-270) + abs(mfmc(3,2)-288) + abs(mfmc(3,3)-306)
      IF (mfm3/=0) CALL prterr(kw,klog,ncase,nerror)

      ncase = 448
      DO i = 1, 2
        DO j = 1, 2
          mima(i,j) = 2*(j-1) + i + 20
          mimb(i,j) = 2*(i-1) + j + 30
        END DO
      END DO
      mimc = matmul(mima,mimb)
      mim3 = abs(mimc(1,1)-1410) + abs(mimc(1,2)-1454) + abs(mimc(2,1)-1474) + &
        abs(mimc(2,2)-1520)
      IF (mim3/=0) CALL prterr(kw,klog,ncase,nerror)

      ncase = 449
      DO i = 1, 2
        DO j = 1, 3
          mzma(i,j) = cmplx(to_fm(2*(j-1)+i+10),to_fm(2*(j-1)+i+20))
        END DO
      END DO
      DO i = 1, 3
        DO j = 1, 4
          mzmb(i,j) = cmplx(to_fm(4*(i-1)+j+50),to_fm(4*(i-1)+j+30))
        END DO
      END DO
      mzmc = matmul(mzma,mzmb)
      mfm3 = abs(mzmc(1,1)-to_zm('-270 + 5192 i')) + &
        abs(mzmc(1,2)-to_zm('-300 + 5300 i')) + abs(mzmc(1,3)-to_zm( &
        '-330 + 5408 i')) + abs(mzmc(1,4)-to_zm('-360 + 5516 i')) + &
        abs(mzmc(2,1)-to_zm('-210 + 5462 i')) + abs(mzmc(2,2)-to_zm( &
        '-240 + 5576 i')) + abs(mzmc(2,3)-to_zm('-270 + 5690 i')) + &
        abs(mzmc(2,4)-to_zm('-300 + 5804 i'))
      IF (mfm3/=0) CALL prterr(kw,klog,ncase,nerror)

      ncase = 450
      mfm3 = max(mfm1,mfm2)
      CALL fm_max(mfm1,mfm2,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 451
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = max(mfm2,mfm1,mfm4)
      CALL fm_max(mfm1,mfm4,mfm4)
      CALL fm_max(mfm2,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 452
      mim3 = max(mim1,mim2)
      CALL im_max(mim1,mim2,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 453
      CALL im_st2m('7654',mim4)
      CALL im_st2m('-1654',mim3)
      mim3 = max(mim2,mim1,mim3,mim4)
      CALL im_st2m('7654',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 454
      j = maxexponent(mfm1)
      IF (j/=int(mxexp)+1) CALL prterr(kw,klog,ncase,nerror)

      ncase = 455
      mfm3 = min(mfm1,mfm2)
      CALL fm_min(mfm1,mfm2,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 456
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = min(mfm2,mfm1,mfm4)
      CALL fm_min(mfm1,mfm4,mfm4)
      CALL fm_min(mfm2,mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 457
      mim3 = min(mim1,mim2)
      CALL im_min(mim1,mim2,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 458
      CALL im_st2m('7654',mim4)
      CALL im_st2m('-1654',mim3)
      mim3 = min(mim2,mim1,mim3,mim4)
      CALL im_st2m('-1654',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test15

    SUBROUTINE test16(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3, &
        mzm4,nerror,ncase,klog)

!             Test functions MINEXPONENT, ..., RRSPACING.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm3, mfm4
      TYPE (im) :: mim1, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (fm) :: mfm5
! ..
! .. Local Scalars ..
      INTEGER :: j
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 459
      j = minexponent(mfm1)
      IF (j/=-int(mxexp)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 460
      CALL fm_st2m('8',mfm3)
      CALL fm_st2m('5',mfm4)
      mfm3 = mod(mfm3,mfm4)
      CALL fm_st2m('3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 461
      CALL fm_st2m('-8',mfm3)
      CALL fm_st2m('5',mfm4)
      mfm3 = mod(mfm3,mfm4)
      CALL fm_st2m('-3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 462
      CALL fm_st2m('8',mfm3)
      CALL fm_st2m('-5',mfm4)
      mfm3 = mod(mfm3,mfm4)
      CALL fm_st2m('3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 463
      CALL fm_st2m('-8',mfm3)
      CALL fm_st2m('-5',mfm4)
      mfm3 = mod(mfm3,mfm4)
      CALL fm_st2m('-3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 464
      CALL im_st2m('8',mim3)
      CALL im_st2m('5',mim4)
      mim3 = mod(mim3,mim4)
      CALL im_st2m('3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 465
      CALL im_st2m('-8',mim3)
      CALL im_st2m('5',mim4)
      mim3 = mod(mim3,mim4)
      CALL im_st2m('-3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 466
      CALL im_st2m('8',mim3)
      CALL im_st2m('-5',mim4)
      mim3 = mod(mim3,mim4)
      CALL im_st2m('3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 467
      CALL im_st2m('-8',mim3)
      CALL im_st2m('-5',mim4)
      mim3 = mod(mim3,mim4)
      CALL im_st2m('-3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 468
      CALL fm_st2m('8',mfm3)
      CALL fm_st2m('5',mfm4)
      mfm3 = modulo(mfm3,mfm4)
      CALL fm_st2m('3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 469
      CALL fm_st2m('-8',mfm3)
      CALL fm_st2m('5',mfm4)
      mfm3 = modulo(mfm3,mfm4)
      CALL fm_st2m('2',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 470
      CALL fm_st2m('8',mfm3)
      CALL fm_st2m('-5',mfm4)
      mfm3 = modulo(mfm3,mfm4)
      CALL fm_st2m('-2',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 471
      CALL fm_st2m('-8',mfm3)
      CALL fm_st2m('-5',mfm4)
      mfm3 = modulo(mfm3,mfm4)
      CALL fm_st2m('-3',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 472
      CALL im_st2m('8',mim3)
      CALL im_st2m('5',mim4)
      mim3 = modulo(mim3,mim4)
      CALL im_st2m('3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 473
      CALL im_st2m('-8',mim3)
      CALL im_st2m('5',mim4)
      mim3 = modulo(mim3,mim4)
      CALL im_st2m('2',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 474
      CALL im_st2m('8',mim3)
      CALL im_st2m('-5',mim4)
      mim3 = modulo(mim3,mim4)
      CALL im_st2m('-2',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 475
      CALL im_st2m('-8',mim3)
      CALL im_st2m('-5',mim4)
      mim3 = modulo(mim3,mim4)
      CALL im_st2m('-3',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 476
      CALL fm_st2m('0',mfm4)
      CALL fm_st2m('1',mfm3)
      CALL fm_big(mfm5)
      CALL fm_div(mfm3,mfm5,mfm5)
      mfm3 = nearest(mfm4,mfm3)
      IF (mfm3/=mfm5) CALL prterr(kw,klog,ncase,nerror)

      ncase = 477
      CALL fm_st2m('0',mfm4)
      CALL fm_st2m('-1',mfm3)
      CALL fm_big(mfm5)
      CALL fm_div(mfm3,mfm5,mfm5)
      mfm3 = nearest(mfm4,mfm3)
      IF (mfm3/=mfm5) CALL prterr(kw,klog,ncase,nerror)

      ncase = 478
      CALL fm_st2m('2.345',mfm4)
      CALL fm_st2m('1',mfm3)
      mfm3 = nearest(mfm4,mfm3)
      CALL fm_ulp(mfm4,mfm5)
      CALL fm_add(mfm4,mfm5,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 479
      CALL fm_st2m('2.345',mfm4)
      CALL fm_st2m('-1',mfm3)
      mfm3 = nearest(mfm4,mfm3)
      CALL fm_ulp(mfm4,mfm5)
      CALL fm_sub(mfm4,mfm5,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 480
      CALL fm_st2m('1',mfm4)
      CALL fm_st2m('-1',mfm3)
      mfm3 = nearest(mfm4,mfm3)
      CALL fm_st2m('0.99',mfm5)
      CALL fm_ulp(mfm5,mfm5)
      CALL fm_sub(mfm4,mfm5,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 481
      CALL fm_st2m('-1',mfm4)
      CALL fm_st2m('12',mfm3)
      mfm3 = nearest(mfm4,mfm3)
      CALL fm_st2m('-0.99',mfm5)
      CALL fm_ulp(mfm5,mfm5)
      CALL fm_sub(mfm4,mfm5,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 482
      mim3 = nint(mfm1)
      CALL fm_nint(mfm1,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 483
      mim3 = nint(mim1)
      CALL im_eq(mim1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 484
      mim3 = nint(mzm1)
      CALL zm_nint(mzm1,mzm4)
      CALL zm_real(mzm4,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 485
      j = precision(mfm1)
      IF (j/=int(log10(real(mbase))*(ndig-1)+1)) CALL prterr(kw,klog,ncase, &
        nerror)

      ncase = 486
      j = precision(mzm1)
      IF (j/=int(log10(real(mbase))*(ndig-1)+1)) CALL prterr(kw,klog,ncase, &
        nerror)

      ncase = 487
      j = radix(mfm1)
      IF (j/=int(mbase)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 488
      j = radix(mim1)
      IF (j/=int(mbase)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 489
      j = radix(mzm1)
      IF (j/=int(mbase)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 490
      j = range(mfm1)
      IF (j/=int(mxexp*log10(real(mbase)))) CALL prterr(kw,klog,ncase,nerror)

      ncase = 491
      j = range(mim1)
      IF (j/=int(ndigmx*log10(real(mbase)))) CALL prterr(kw,klog,ncase,nerror)

      ncase = 492
      j = range(mzm1)
      IF (j/=int(mxexp*log10(real(mbase)))) CALL prterr(kw,klog,ncase,nerror)

      ncase = 493
      mfm3 = real(mfm1)
      CALL fm_eq(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 494
      mfm3 = real(mim1)
      CALL im_i2fm(mim1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 495
      mfm3 = real(mzm1)
      CALL zm_real(mzm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 496
      mfm3 = rrspacing(mfm1)
      CALL fm_abs(mfm1,mfm4)
      mfm4%mfm(1) = ndig
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test16

    SUBROUTINE test17(i1,r1,d1,z1,c1,mfm1,mfm2,mfm3,mfm4,mim1,mim2,mim3,mim4, &
        mzm1,mzm3,mzm4,nerror,ncase,klog)

!             Test functions SCALE, ..., TINY.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm2, mfm3, mfm4
      TYPE (im) :: mim1, mim2, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      ncase = 497
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = scale(mfm4,1)
      CALL fm_mpyi(mfm4,int(mbase),mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 498
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = scale(mzm4,-2)
      CALL zm_divi(mzm4,int(mbase),mzm4)
      CALL zm_divi(mzm4,int(mbase),mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 499
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = setexponent(mfm4,1)
      CALL fm_mpyi(mfm4,int(mbase),mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 500
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = sign(mfm4,mfm2)
      CALL fm_sign(mfm4,mfm2,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 501
      CALL im_st2m('231',mim4)
      mim3 = sign(mim4,mim2)
      CALL im_sign(mim4,mim2,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 502
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = sin(mfm4)
      CALL fm_sin(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 503
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = sin(mzm4)
      CALL zm_sin(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 504
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = sinh(mfm4)
      CALL fm_sinh(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 505
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = sinh(mzm4)
      CALL zm_sinh(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 506
      CALL fm_st2m('-0.7654',mfm4)
      mfm3 = spacing(mfm4)
      CALL fm_ulp(mfm4,mfm4)
      CALL fm_abs(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 507
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = sqrt(mfm4)
      CALL fm_sqrt(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 508
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = sqrt(mzm4)
      CALL zm_sqrt(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 509
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = tan(mfm4)
      CALL fm_tan(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 510
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = tan(mzm4)
      CALL zm_tan(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 511
      CALL fm_st2m('0.7654',mfm4)
      mfm3 = tanh(mfm4)
      CALL fm_tanh(mfm4,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 512
      CALL zm_st2m('0.7654 - 0.3456 i',mzm4)
      mzm3 = tanh(mzm4)
      CALL zm_tanh(mzm4,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 513
      CALL fm_big(mfm4)
      CALL fm_i2m(1,mfm3)
      CALL fm_div(mfm3,mfm4,mfm4)
      mfm3 = tiny(mfm1)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 514
      mim3 = tiny(mim1)
      CALL im_i2m(1,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 515
      CALL fm_big(mfm4)
      CALL fm_i2m(1,mfm3)
      CALL fm_div(mfm3,mfm4,mfm4)
      CALL zm_cmpx(mfm4,mfm4,mzm4)
      mzm3 = tiny(mzm1)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test17

    SUBROUTINE test18(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3, &
        mzm4,nerror,ncase,klog)

!             Test functions TO_FM, TO_IM, TO_ZM, ..., TO_DPZ.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm3, mfm4
      TYPE (im) :: mim1, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (fm) :: mfm5
! ..
! .. Local Scalars ..
      COMPLEX (kind(0.0D0)) :: c2
      COMPLEX :: z2
      REAL (KIND(0.0D0)) :: d2, d3, dsmall
      REAL :: r2, rsmall
      INTEGER :: i2
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..

      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0

      ncase = 516
      mfm3 = to_fm(123)
      CALL fm_i2m(123,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 517
      mfm3 = to_fm(123.4)
      CALL fm_sp2m(123.4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 518
      mfm3 = to_fm(123.45D0)
      CALL fm_dp2m(123.45D0,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 519
      mfm3 = to_fm(cmplx(123.4,567.8))
      CALL fm_sp2m(123.4,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 520
      mfm3 = to_fm(cmplx(123.4D0,567.8D0,kind(1.0D0)))
      CALL fm_dp2m(123.4D0,mfm4)
      CALL fm_sub(mfm3,mfm4,mfm4)
      CALL fm_div(mfm4,mfm3,mfm4)
      CALL fm_abs(mfm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 521
      mfm3 = to_fm(mfm1)
      CALL fm_eq(mfm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 522
      mfm3 = to_fm(mim1)
      CALL im_i2fm(mim1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 523
      mfm3 = to_fm(mzm1)
      CALL zm_real(mzm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 524
      mfm3 = to_fm('-123.654')
      CALL fm_st2m('-123.654',mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 525
      mim3 = to_im(123)
      CALL im_i2m(123,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 526
      mim3 = to_im(123.4)
      CALL fm_sp2m(123.4,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 527
      mim3 = to_im(123.45D0)
      CALL fm_dp2m(123.45D0,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 528
      mim3 = to_im(cmplx(123.4,567.8))
      CALL fm_sp2m(123.4,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 529
      mim3 = to_im(cmplx(123.4D0,567.8D0,kind(1.0D0)))
      CALL fm_dp2m(123.4D0,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 530
      mim3 = to_im(mfm1)
      CALL fm_eq(mfm1,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 531
      mim3 = to_im(mim1)
      CALL im_i2fm(mim1,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 532
      mim3 = to_im(mzm1)
      CALL zm_real(mzm1,mfm4)
      CALL im_fm2i(mfm4,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 533
      mim3 = to_im('-123654')
      CALL im_st2m('-123654',mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 534
      mzm3 = to_zm(123)
      CALL zm_i2m(123,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 535
      mzm3 = to_zm(123.4)
      CALL fm_sp2m(123.4,mfm4)
      CALL fm_i2m(0,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 536
      mzm3 = to_zm(123.45D0)
      CALL fm_dp2m(123.45D0,mfm4)
      CALL fm_i2m(0,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 537
      mzm3 = to_zm(cmplx(123.4,567.8))
      CALL fm_sp2m(123.4,mfm4)
      CALL fm_sp2m(567.8,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = rsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 538
      mzm3 = to_zm(cmplx(123.4D0,567.8D0,kind(1.0D0)))
      CALL fm_dp2m(123.4D0,mfm4)
      CALL fm_dp2m(567.8D0,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      CALL zm_sub(mzm3,mzm4,mzm4)
      CALL zm_div(mzm4,mzm3,mzm4)
      CALL zm_abs(mzm4,mfm4)
      mfm3 = dsmall
      IF (fm_comp(mfm4,'gt',mfm3)) CALL prterr(kw,klog,ncase,nerror)

      ncase = 539
      mzm3 = to_zm(mfm1)
      CALL fm_eq(mfm1,mfm4)
      CALL fm_i2m(0,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 540
      mzm3 = to_zm(mim1)
      CALL im_i2fm(mim1,mfm4)
      CALL fm_i2m(0,mfm5)
      CALL zm_cmpx(mfm4,mfm5,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 541
      mzm3 = to_zm(mzm1)
      CALL zm_eq(mzm1,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 542
      mzm3 = to_zm('-123.654 + 98.7 i')
      CALL zm_st2m('-123.654 + 98.7 i',mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 543
      CALL fm_m2i(mfm1,i2)
      IF (to_int(mfm1)/=i2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 544
      CALL im_m2i(mim1,i2)
      IF (to_int(mim1)/=i2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 545
      CALL zm_m2i(mzm1,i2)
      IF (to_int(mzm1)/=i2) CALL prterr(kw,klog,ncase,nerror)

      ncase = 546
      CALL fm_m2sp(mfm1,r2)
      IF (abs((to_sp(mfm1)-r2)/r2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 547
      CALL im_m2dp(mim1,d2)
      r2 = d2
      IF (abs((to_sp(mim1)-r2)/r2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 548
      CALL zm_real(mzm1,mfm4)
      CALL fm_m2sp(mfm4,r2)
      IF (abs((to_sp(mzm1)-r2)/r2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 549
      CALL fm_m2dp(mfm1,d2)
      IF (abs((to_dp(mfm1)-d2)/d2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 550
      CALL im_m2dp(mim1,d2)
      IF (abs((to_dp(mim1)-d2)/d2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 551
      CALL zm_real(mzm1,mfm4)
      CALL fm_m2dp(mfm4,d2)
      IF (abs((to_dp(mzm1)-d2)/d2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 552
      CALL fm_m2sp(mfm1,r2)
      z2 = r2
      IF (abs((to_spz(mfm1)-z2)/z2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 553
      CALL im_m2dp(mim1,d2)
      z2 = d2
      IF (abs((to_spz(mim1)-z2)/z2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 554
      CALL zm_m2z(mzm1,z2)
      IF (abs((to_spz(mzm1)-z2)/z2)>rsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 555
      CALL fm_m2dp(mfm1,d2)
      c2 = d2
      IF (abs((to_dpz(mfm1)-c2)/c2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 556
      CALL im_m2dp(mim1,d2)
      c2 = d2
      IF (abs((to_dpz(mim1)-c2)/c2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

      ncase = 557
      CALL zm_real(mzm1,mfm4)
      CALL fm_m2dp(mfm4,d2)
      CALL zm_imag(mzm1,mfm4)
      CALL fm_m2dp(mfm4,d3)
      c2 = cmplx(d2,d3,kind(0.0D0))
      IF (abs((to_dpz(mzm1)-c2)/c2)>dsmall) THEN
          CALL prterr(kw,klog,ncase,nerror)
      ENDIF

    END SUBROUTINE test18

    SUBROUTINE test19(i1,r1,d1,z1,c1,mfm1,mfm3,mfm4,mim1,mim3,mim4,mzm1,mzm3, &
        mzm4,nerror,ncase,klog)

!             Test the derived-type interface routines that are not
!             used elsewhere in this program.

! .. Use Statements ..
      USE fmzm
! ..
! .. Intrinsic Functions ..
!      INTRINSIC kind
! ..
! .. Structure Arguments ..
      TYPE (fm) :: mfm1, mfm3, mfm4
      TYPE (im) :: mim1, mim3, mim4
      TYPE (zm) :: mzm1, mzm3, mzm4
! ..
! .. Scalar Arguments ..
      COMPLEX (kind(0.0D0)) :: c1
      COMPLEX :: z1
      REAL (KIND(0.0D0)) :: d1
      REAL :: r1
      INTEGER :: i1, klog, ncase, nerror
! ..
! .. Local Structures ..
      TYPE (im) :: mim2
      TYPE (fm) :: msmall
! ..
! .. Local Scalars ..
      COMPLEX :: z3, z4
      REAL (KIND(0.0D0)) :: d3, d4, dsmall
      REAL :: r3, r4, rsmall
      INTEGER :: i3, i4
      CHARACTER (80) :: string
! ..
! .. External Subroutines ..
      EXTERNAL prterr
! ..
      rsmall = epsilon(1.0)*100.0
      dsmall = epsilon(1.0D0)*100.0
      msmall = epsilon(to_fm(1))*10000.0

      ncase = 558
      mfm3 = mfm1 + 123
      mfm4 = mfm1
      CALL fm_addi(mfm4,123)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 559
      CALL fm_chsh(mfm1,mfm4,mfm3)
      mfm3 = cosh(mfm1)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 560
      CALL fm_chsh(mfm1,mfm3,mfm4)
      mfm3 = sinh(mfm1)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 561
      CALL fm_cssn(mfm1,mfm4,mfm3)
      mfm3 = cos(mfm1)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 562
      CALL fm_cssn(mfm1,mfm3,mfm4)
      mfm3 = sin(mfm1)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 563
      mfm3 = mfm1/123
      CALL fm_divi(mfm1,123,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 564
      mfm3 = 123.45D0
      CALL fm_dpm(123.45D0,mfm4)
      IF (abs((mfm3-mfm4)/mfm4)>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 565
      CALL fm_form('F70.56',mfm1,string)
      CALL fm_st2m(string(1:70),mfm4)
      IF (abs((mfm1-mfm4)/mfm4)>msmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 566
      mfm3 = mfm1**123
      CALL fm_ipwr(mfm1,123,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 567
      mfm3 = log(to_fm(123))
      CALL fm_lni(123,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 568
      d3 = mfm1
      CALL fm_m2dp(mfm1,d4)
      IF (abs((d3-d4)/d3)>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 569
      i3 = mfm1
      CALL fm_m2i(mfm1,i4)
      IF (i3/=i4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 570
      r3 = mfm1
      CALL fm_m2sp(mfm1,r4)
      IF (abs((r3-r4)/r3)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 571
      mfm3 = 2.67
      CALL fm_mod(mfm1,mfm3,mfm4)
      mfm3 = mod(mfm1,mfm3)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 572
      CALL fm_pi(mfm4)
      mfm3 = 4*atan(to_fm(1))
      IF (abs((mfm3-mfm4)/mfm4)>msmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 573
      mfm3 = mfm1**(to_fm(1)/to_fm(3))
      CALL fm_rpwr(mfm1,1,3,mfm4)
      IF (abs((mfm3-mfm4)/mfm4)>msmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 574
      CALL fm_sqr(mfm1,mfm4)
      mfm3 = mfm1*mfm1
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 575
      mim3 = mim1/13
      CALL im_divi(mim1,13,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 576
      mim3 = 13
      CALL im_divr(mim1,mim3,mim3,mim4)
      mim3 = mod(mim1,mim3)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 577
      mim3 = 13
      CALL im_divr(mim1,mim3,mim3,mim4)
      mim4 = mim1/13
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 578
      mim3 = mim1/13
      CALL im_dvir(mim1,13,mim4,i4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 579
      i3 = mod(mim1,to_im(13))
      CALL im_dvir(mim1,13,mim4,i4)
      IF (i3/=i4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 580
      CALL im_form('I70',mim1,string)
      CALL im_st2m(string(1:70),mim4)
      IF (mim1/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 581
      mim3 = 40833
      mim4 = 16042
      CALL im_gcd(mim3,mim4,mim4)
      IF (mim4/=13) CALL prterr(kw,klog,ncase,nerror)

      ncase = 582
      d3 = mim1
      CALL im_m2dp(mim1,d4)
      IF (abs((d3-d4)/d3)>dsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 583
      i3 = mim1
      CALL im_m2i(mim1,i4)
      IF (i3/=i4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 584
      mim3 = 6
      CALL im_mod(mim1,mim3,mim4)
      mim3 = mod(mim1,mim3)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 585
      mim3 = mim1*123
      CALL im_mpyi(mim1,123,mim4)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 586
      mim2 = 3141
      mim3 = 133
      CALL im_mpym(mim1,mim2,mim3,mim4)
      mim3 = mod(mim1*mim2,mim3)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 587
      mim2 = 31
      mim3 = 147
      CALL im_pmod(mim1,mim2,mim3,mim4)
      mim3 = mod(mim1**mim2,mim3)
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 588
      CALL im_sqr(mim1,mim4)
      mim3 = mim1*mim1
      IF (mim3/=mim4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 589
      mzm3 = mzm1 + 123
      mzm4 = mzm1
      CALL zm_addi(mzm4,123)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 590
      mfm3 = atan2(aimag(mzm1),real(mzm1))
      CALL zm_arg(mzm1,mfm4)
      IF (mfm3/=mfm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 591
      CALL zm_chsh(mzm1,mzm4,mzm3)
      mzm3 = cosh(mzm1)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 592
      CALL zm_chsh(mzm1,mzm3,mzm4)
      mzm3 = sinh(mzm1)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 593
      CALL zm_cssn(mzm1,mzm4,mzm3)
      mzm3 = cos(mzm1)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 594
      CALL zm_cssn(mzm1,mzm3,mzm4)
      mzm3 = sin(mzm1)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 595
      CALL zm_form('F35.26','F35.26',mzm1,string)
      CALL zm_st2m(string(1:75),mzm4)
      IF (abs((mzm1-mzm4)/mzm4)>msmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 596
      mzm3 = to_zm('123-456i')
      CALL zm_2i2m(123,-456,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 597
      mzm3 = mzm1**123
      CALL zm_ipwr(mzm1,123,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 598
      i3 = mzm1
      CALL zm_m2i(mzm1,i4)
      IF (i3/=i4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 599
      z3 = mzm1
      CALL zm_m2z(mzm1,z4)
      IF (abs((z3-z4)/z3)>rsmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 600
      mzm3 = mzm1*123
      CALL zm_mpyi(mzm1,123,mzm4)
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 601
      mzm3 = mzm1**(to_zm(1)/to_zm(3))
      CALL zm_rpwr(mzm1,1,3,mzm4)
      IF (abs((mzm3-mzm4)/mzm4)>msmall) CALL prterr(kw,klog,ncase,nerror)

      ncase = 602
      CALL zm_sqr(mzm1,mzm4)
      mzm3 = mzm1*mzm1
      IF (mzm3/=mzm4) CALL prterr(kw,klog,ncase,nerror)

      ncase = 603
      mzm3 = z1
      CALL zm_z2m(z1,mzm4)
      IF (abs((mzm3-mzm4)/mzm3)>rsmall) CALL prterr(kw,klog,ncase,nerror)

    END SUBROUTINE test19

    SUBROUTINE prterr(kw,klog,ncase,nerror)
! .. Scalar Arguments ..
      INTEGER :: klog, kw, ncase, nerror
! ..

      WRITE (kw,*) ' Error in case ', ncase
      WRITE (klog,*) ' Error in case ', ncase
      nerror = nerror + 1

    END SUBROUTINE prterr
