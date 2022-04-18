! Last Change:  Mofreh Zaghloul   October 6, 2018

!.. For compilation
! >ifort set_rk.f90 rk_erfcx_Cody.f90 wofz.f90 wofz_v2.f90 Faddeyeva_v2_mod_rk.f90
!  wofz_680_driver.f90 -o wofz_680_driver
!!  /fpe:0 /traceback
!! /warn:all
!! /warn:declarations
!! /CU

PROGRAM wofz_680_driver

  ! .. Use Statements ..
  USE set_rk, ONLY :rk
  USE Faddeyeva_v2_mod_rk, ONLY : Faddeyeva_v2_rk
  USE wofz_v2_mod, ONLY : wofz_v2
  USE wofzmod, ONLY : wofz
  ! ..
  ! .. Parameters ..
  INTEGER, PARAMETER :: m_max = 40001, n_max = 71, nrepeats=10

  ! ..
  ! .. Local Scalars ..
  REAL(rk) :: time_begin, time_end, execution_time, xtmp, ytmp
  REAL(rk), DIMENSION(m_max):: harvest
  REAL(rk), DIMENSION(m_max*n_max):: tmp_err

  INTEGER :: IM, jn, icase
  LOGICAL :: flag(m_max*n_max), flag1(40)
  ! ..
  ! .. Local Arrays ..
  COMPLEX(rk), DIMENSION(m_max*n_max)::z,w,w_ref
  COMPLEX(rk), DIMENSION(m_max,n_max)::z2
  COMPLEX(rk), DIMENSION(40)::z1,w1,w2,w1_ref

  ! .. Intrinsic Functions ..
  INTRINSIC CMPLX, CPU_TIME, REAL

  ! ..
  ! .. Equivalences ..
  EQUIVALENCE (z,z2)

  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(harvest)

  z1 = (/ (1.83e0_rk,1.0e-20_rk),(1.84e0_rk,1.0e-20_rk),(2.3e0_rk,1.0e-20_rk),(2.5e0_rk,1.0e-20_rk), &
    (3.9e0_rk,1.0e-20_rk),(4.3e0_rk,1.0e-20_rk),(4.5e0_rk,1.0e-20_rk), &
    (1.0e-20_rk,6.3e0_rk),(0.0_rk,4.0e0_rk),(0.0e0_rk,4.3e0_rk),(0.0e0_rk,4.0e0_rk), &
    (6.3e+0_rk,2.0e-1_rk),(6.3e+0_rk,3.0e-4_rk),(6.3e+0_rk,4.5e-6_rk),(6.3e+0_rk,5.0e-11_rk), &
    (6.3e+0_rk,-1.0e-20_rk),(-6.3e+0_rk,1.0e-18_rk),(5.3e+0_rk,1.0e-16_rk), &
    (5.3e+0_rk,1.e-14_rk),(5.3e+0_rk,1.e-12_rk),(5.3e+0_rk,1.0e-10_rk),(5.3e+0_rk,1.0e-8_rk), &
    (5.3e+0_rk,1.0e-6_rk),(5.3e+0_rk,1.0e-4_rk),(5.3e0_rk,1.0e-2_rk),(5.3e+0_rk,1.0e-1_rk), &
    (7.3e+0_rk,1.0e-20_rk),(7.3e+0_rk,1.0e-18_rk),(7.3e+0_rk,1.0e-16_rk), &
    (7.3e+0_rk,1.e-14_rk),(7.3e+0_rk,1.e-12_rk),(7.3e+0_rk,1.0e-10_rk),(17.3e+0_rk,1.0e-20_rk), &
    (21.3e+0_rk,1.0e-20_rk),(22.3e+0_rk,0.0e0_rk),(27.0e+0_rk,0.0e0_rk),(7.3e+0_rk,1.0e-1_rk),&
    (7.3_rk,0.0_rk),(5.3_rk,0.0_rk),(6.4_rk,1.0e-20_rk)/)

  write(*,*)'*********************************************************************'
  write(*,*)'Sample results showing the relative error in the real and imaginary parts'
  write(*,*)'from the "Original 680 code" and the "Revised 680 code" for points scattered'
  write(*,*)'over the region in which the Original 680 code loses its accuracy'
  write(*,*)'*********************************************************************'
  WRITE(*,'((50x,A22,8x,A22))')'Original 680','Revised 680'
  WRITE(*,'((10x,A2,22x,A2,15x,A12,2x,A12,4x,A12,2x,A12))')'x', 'y ','Error V ','Error L','Error V ','Error L'

  call wofz(z1,w1,flag1)
  call wofz_v2(z1,w2,flag1)

  ! Generate best accuracy using Algorithm 916 with Ndgts=13
  call Faddeyeva_v2_rk(z1,w1_ref,13)

  do jn=1,40
    if (REAL(w1_ref(jn)) .ne. 0.0_rk .and. Aimag(w1_ref(jn)) &
        .ne. 0.0_rk)then
      WRITE (*,'((2es23.15E3,4x),4es15.5)') REAL(z1(jn)), AIMAG(z1(jn)),&
        (REAL(w1(jn))-REAL(w1_ref(jn)))/REAL(w1_ref(jn)),&
        (AIMAG(w1(jn))-AIMAG(w1_ref(jn)))/AIMAG(w1_ref(jn)),&
        (REAL(w2(jn))-REAL(w1_ref(jn)))/REAL(w1_ref(jn)),&
        (AIMAG(w2(jn))-AIMAG(w1_ref(jn)))/AIMAG(w1_ref(jn))
    endif

  enddo
  write(*,*)'   '
  write(*,*)'   '


  WRITE (*,'(/3a)') 'Summary of execution time and accuracy',&
    'in computing 2,840,071 input values',&
    'for z with 71 values for y and 40001 values for x.'
  write(*,*)'___________________________'

  ! ..
  ! .. define x and y arrays
  DO icase=1,4     ! Loop on the four tested cases
    DO IM = 1, m_max
      DO jn = 1, n_max
        !.. Case 1
        IF (icase==1)THEN
          xtmp = -5.0e2_rk+REAL((IM-1),KIND=rk)*(10.0e2_rk/&
            REAL((m_max-1),KIND=rk))
          ytmp = 1.0e1_rk**(-5.0e0_rk+REAL((jn-1),KIND=rk)*1.0e1_rk/&
            REAL((n_max-1),KIND=rk))

          !.. Case 2
        ELSEIF (icase==2)THEN
          xtmp = -2.0e2_rk+REAL((IM-1),KIND=rk)*(4.0e2_rk/&
            REAL((m_max-1),KIND=rk))
          ytmp = 1.0e1_rk**(-2.0e1_rk+REAL((jn-1),KIND=rk)*2.4e1_rk/&
            REAL((n_max-1),KIND=rk))

          !.. Case 3
        ELSEIF (icase==3)THEN
          xtmp = -1.0e1_rk+REAL((IM-1),KIND=rk)*(2.0e1_rk/&
            REAL((m_max-1),KIND=rk))
          ytmp = 1.0e1_rk**(-5.0e0_rk+REAL((jn-1),KIND=rk)*1.0e1_rk/&
            REAL((n_max-1),KIND=rk))

          !.. Case 4 (|z|^2<=36)
        ELSEIF (icase==4)THEN

          ytmp = 1.0e1_rk**(-2.0e1_rk+REAL((jn-1),KIND=rk)*&
            2.077815e+001_rk/REAL((n_max-1),KIND=rk))
          xtmp = -6.0_rk+harvest(IM)*2.0_rk*&
            (36.0_rk-ytmp*ytmp)**0.5_rk

        END IF
        z2(IM,jn) = CMPLX(xtmp,ytmp,KIND=rk)
      END DO
    END DO

    write(*,*)'________'
    write(*,'(/a,I5)')'case',icase
    if (icase==1) then
      WRITE (*,'(/2a)')'y uniformly spaced on the logarithmic scale between 1e-5 and 1e5',&
        ' and x uniformly spaced on the linear scale between -500 and 500'
    elseif(icase==2) then
      WRITE (*,'(/2a)')'y uniformly spaced on the logarithmic scale between 1e-20 and 1e4',&
        ' and x uniformly spaced on the linear scale between -200 and 200'
    elseif(icase==3) then
      WRITE (*,'(/2a)')'y uniformly spaced on the logarithmic scale between 1e-5 and 1e5',&
        ' and x uniformly spaced on the linear scale between -10 and 10'
    elseif(icase==4) then
      WRITE (*,'(/2a)')'y uniformly spaced on the logarithmic scale between 1e-20 and 1e6',&
        ' and x randomly generated to satisfy |z|<=6'
    endif


    ! Generate best accuracy using Zaghloul with Ndgts=13
    CALL Faddeyeva_v2_rk(z,w_ref,13)

    CALL CPU_TIME (time_begin)
    DO IM = 1,nrepeats
      CALL wofz(z,w,flag)
    ENDDO
    CALL CPU_TIME (time_end)
    execution_time=time_end-time_begin
    WRITE (*,'(/a)') 'Original 680:'
    WRITE (*,'(/A,g13.6,a)')'Total execution time  =', execution_time, &
      ' processor dependent units'
    WRITE (*,'(a,g13.6,a)') 'Average time per call =', &
      execution_time/(m_max*n_max*real(nrepeats, kind=rk)), &
      ' processor dependent units'

    WRITE (*,*)'  '

    WHERE (ABS(REAL(w_ref)) .ne. 0.0_rk)
      tmp_err=((ABS(REAL(w)-REAL(w_ref))/ABS(REAL(w_ref))))
    ENDWHERE
    PRINT '(a,g13.6,a)','max_rel_err_Re    =', MAXVAL(tmp_err)

    WHERE (ABS(AIMAG(w_ref)) .ne. 0.0_rk)
      tmp_err=((ABS(AIMAG(w)-AIMAG(w_ref))/ABS(AIMAG(w_ref))))
    ENDWHERE
    PRINT '(a,g13.6,a)','max_rel_err_Im    =', MAXVAL(tmp_err)
    WRITE (*,*)'  '
    WRITE (*,*)'  '


    CALL CPU_TIME (time_begin)
    DO IM = 1,nrepeats
      CALL wofz_v2(z,w,flag)
    ENDDO
    CALL CPU_TIME (time_end)
    execution_time=time_end-time_begin
     WRITE (*,'(/a)') 'Revised 680:'
    WRITE (*,'(/A,g13.6,a)')'Total execution time  =', execution_time, &
      ' processor dependent units'
    WRITE (*,'(a,g13.6,a)') 'Average time per call =', &
      execution_time/(m_max*n_max*real(nrepeats, kind=rk)), &
      ' processor dependent units'

    WRITE (*,*)'  '

    WHERE (ABS(REAL(w_ref)) .ne. 0.0_rk)
      tmp_err=((ABS(REAL(w)-REAL(w_ref))/ABS(REAL(w_ref))))
    ENDWHERE
    PRINT '(a,g13.6,a)','max_rel_err_Re    =', MAXVAL(tmp_err)

    WHERE (ABS(AIMAG(w_ref)) .ne. 0.0_rk)
      tmp_err=((ABS(AIMAG(w)-AIMAG(w_ref))/ABS(AIMAG(w_ref))))
    ENDWHERE
    PRINT '(a,g13.6,a)','max_rel_err_Im    =', MAXVAL(tmp_err)
    WRITE (*,*)'  '
    WRITE (*,*)'  '

  ENDDO
  WRITE (*,*)'  '

END PROGRAM wofz_680_driver
