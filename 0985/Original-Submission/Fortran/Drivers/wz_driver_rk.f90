! Last Change:  Mofreh Zaghloul   June 05 2017
! For compilation
! Using ifort
! >ifort set_rk.f90 rk_erfcx_Cody.f90 humlicek0mod.f90 Faddeyeva_v2_mod_rk.f90
!   wz_mod_rk.f90  wz_driver_rk.f90 -o wz_driver_rk
!!!   -warn /fpe:0 /traceback

program wz_driver_rk
    use set_rk, only:rk
    use wz_mod_rk, only: wz_rk
    use Faddeyeva_v2_mod_rk, only: Faddeyeva_v2_rk
    USE humlicekmod, ONLY : humlicek0

    implicit none
    integer :: im, jn, icase
    integer, parameter :: m_max=40001, n_max=71

    complex(rk), dimension(m_max*n_max)::z,w,w_ref
    complex(rk), dimension(m_max,n_max)::z2
    complex(rk), dimension(40)::z1,w1,w1_ref
    equivalence ( z, z2 )
    real(rk), dimension(m_max*n_max)::dvdx,dvdy,tmp_err
    real(rk), dimension(m_max):: harvest
    real(rk) xtmp, ytmp, time_begin, time_end
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(harvest)

    z1 = (/ (6.3e-2_rk,1.0e-20_rk),(6.3e-2_rk,1.0e-14_rk),(6.3e-2_rk,1.0e-12_rk), (6.3e-2_rk,1.0e-10_rk), &
        (6.3e-2_rk,1.0e-06_rk),(6.3e-2_rk,1.0e-2_rk),(6.3e-2_rk,1.0e+1_rk), &
        (6.3e-2_rk,1.2e+01_rk),(6.3e-2_rk,1.5e+01_rk),(6.3e-2_rk,2.0e+02_rk),(6.3e-2_rk,1.0e+05_rk), &
        (6.3e+0_rk,1.0e-20_rk),(6.3e+0_rk,1.0e-14_rk),(6.3e+0_rk,1.0e-12_rk),(6.3e+0_rk,1.0e-10_rk), &
        (6.3e+0_rk,1.0e-06_rk),(6.3e+0_rk,1.0e-2_rk),(6.3e+0_rk,1.0e+1_rk), &
        (6.3e+0_rk,1.2e+01_rk),(6.3e+0_rk,1.5e+01_rk),(6.3e+0_rk,2.0e+02_rk),(6.3e+0_rk,1.0e+05_rk), &
        (6.3e+2_rk,1.0e-20_rk),(6.3e+2_rk,1.0e-14_rk),(6.3e+2_rk,1.0e-12_rk),(6.3e+2_rk,1.0e-10_rk), &
        (6.3e+2_rk,1.0e-06_rk),(6.3e+2_rk,1.0e-2_rk),(6.3e+2_rk,1.0e+1_rk), &
        (6.3e+2_rk,1.2e+01_rk),(6.3e+2_rk,1.5e+01_rk),(6.3e+2_rk,2.0e+02_rk),(6.3e+2_rk,1.0e+05_rk), &
        (1.0e+0_rk,1.0e-20_rk),(5.5e+0_rk,1.0e-14_rk),(3.9e+4_rk,1.0e+0_rk),(1.0e+0_rk,2.8e+4_rk), &
        (3.2927225E-10_rk,4.2949670E+09_rk),(4.2949670E+09_rk,3.2927225E-10_rk),(4.2949670E+09_rk,4.2949670E+09_rk)/)


    write(*,*)'Sample Results from wz'
    call wz_rk(z1,w1)
    ! Generate best accuracy using Algorithm 916 with Ndgts=13
    call Faddeyeva_v2_rk(z1,w1_ref,13)
    do jn=1,40
        write (*,'(3(2es23.15E3,5x),2es15.5)') real(z1(jn)), aimag(z1(jn)), &
            real( w1(jn)), aimag(w1(jn)), real( w1_ref(jn)),aimag(w1_ref(jn)), &
            (real(w1(jn))-real(w1_ref(jn)))/real(w1_ref(jn)),&
            (aimag(w1(jn))-aimag(w1_ref(jn)))/aimag(w1_ref(jn))
    enddo
    write(*,*)'   '

    ! ..
    ! .. define x and y arrays
    DO icase=1,4     ! Loop on the four tested cases
        DO im = 1, m_max
            DO jn = 1, n_max
                !.. Case 1
                IF (icase==1)THEN
                    xtmp = -5.0e2_rk+REAL((IM-1),KIND=rk)*(10.0e2_rk/REAL((m_max-1),KIND=rk))
                    ytmp = 1.0e1_rk**(-5.0e0_rk+REAL((jn-1),KIND=rk)*1.0e1_rk/REAL((n_max-1),KIND=rk))

                    !.. Case 2
                ELSEIF (icase==2)THEN
                    xtmp = -2.0e2_rk+REAL((IM-1),KIND=rk)*(4.0e2_rk/REAL((m_max-1),KIND=rk))
                    ytmp = 1.0e1_rk**(-2.0e1_rk+REAL((jn-1),KIND=rk)*2.4e1_rk/REAL((n_max-1),KIND=rk))

                    !.. Case 3
                ELSEIF (icase==3)THEN
                    xtmp = -1.0e1_rk+REAL((IM-1),KIND=rk)*(2.0e1_rk/REAL((m_max-1),KIND=rk))
                    ytmp = 1.0e1_rk**(-5.0e0_rk+REAL((jn-1),KIND=rk)*1.e1_rk/REAL((n_max-1),KIND=rk))

                    !.. Case 4 (|z|^2<=36)
                ELSEIF (icase==4)THEN
                    ytmp = 1.0e1_rk**(-2.0e1_rk+REAL((jn-1),KIND=rk)*2.077815e+001_rk/REAL((n_max-1),KIND=rk))
                    xtmp = -6.0_rk+harvest(IM)*2.0_rk*(36.0_rk-ytmp*ytmp)**0.5_rk
                END IF
                z2(IM,jn) = CMPLX(xtmp,ytmp,KIND=rk)
            END DO
        END DO

        write(*,*)'!!------------------------------------------'
        WRITE (*,'(/a,I5)') 'Results for Case', icase

        ! Generate best accuracy using Algorithm 916 with Ndgts=13
        call Faddeyeva_v2_rk(z,w_ref,13,dvdx,dvdy)

        call cpu_time (time_begin)
        do im = 1,100
            call wz_rk(z,w)
        enddo
        call cpu_time (time_end)

        print '(a,g13.6,a)','elapsed time Wz      =',(time_end-time_begin),' processor dependent units'
        print '(a,g13.6,a)','Average time per cal =',(time_end-time_begin)/(m_max*n_max*100),&
            ' processor dependent units'
        write (*,*)'  '

        where (abs(real(w_ref)) .ne. 0.0_rk)
            tmp_err=((abs(real(w)-real(w_ref))/abs(real(w_ref))))
        endwhere
        print '(a,g13.6,a)','max_rel_err_Re_wz    =', maxval(tmp_err)
        where (abs(aimag(w_ref)) .ne. 0.0_rk)
            tmp_err=((abs(aimag(w)-aimag(w_ref))/abs(aimag(w_ref))))
        endwhere
        print '(a,g13.6,a)','max_rel_err_Im_wz    =', maxval(tmp_err)
        write (*,*)'  '
        write (*,*)'  '

        call cpu_time (time_begin)
        do im = 1,100
            call humlicek0(z,w)
        enddo
        call cpu_time (time_end)

        print '(a,g13.6,a)','elapsed time w4      =',(time_end-time_begin),' processor dependent units'
        print '(a,g13.6,a)','Average time per cal =',(time_end-time_begin)/(m_max*n_max*100),&
            ' processor dependent units'
        write (*,*)'  '

        where (abs(real(w_ref)) .ne. 0.0_rk)
            tmp_err=((abs(real(w)-real(w_ref))/abs(real(w_ref))))
        endwhere
        print '(a,g13.6,a)','max_rel_err_Re_ w4   =', maxval(tmp_err)
        where (abs(aimag(w_ref)) .ne. 0.0_rk)
            tmp_err=((abs(aimag(w)-aimag(w_ref))/abs(aimag(w_ref))))
        endwhere
        print '(a,g13.6,a)','max_rel_err_Im_ w4   =', maxval(tmp_err)
        write (*,*)'  '
        write (*,*)'  '
    end do
    write (*,*)'  '
end program wz_driver_rk
