!     Last change:  CPB  22 JAN 2002   11:10 am
!     Last change:  RH   25 Apr 2001    2:18 pm
!     Last change:  RH    6 Dec 2000    4:23 pm
! Testing elementary use of Pthreads.  Print summary of tests.
! Any failures print 'failed' and cause an abort.
! This code is part of the package "A Fortran Interface to Posix Threads,"  to be
! published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
! This is test4.f90.

     MODULE global_test4
        USE fpthrd
        IMPLICIT NONE
        INTEGER test_number
        Type(fpthrd_t) base_thread, test_thread
        Type(fpthrd_attr_t) saa, sab
        Type(fpthrd_condattr_t) sca, scb
        Type(fpthrd_mutexattr_t) sma
        Type(fsched_param) param
     END MODULE

     program fmain4
        USE global_test4
        implicit none
        integer inheritsched, policy, prioceiling, result, scope, status

        test_number=24

! Test 8. Exercise various attributes fetches and settings.
        call fpthrd_condattr_init(sca, status)
        test_number=test_number+1
! Initizalize, destroy and initialize once again.
        call ferr_abort(test_number, status, " initializing condition attribute structure")
        call fpthrd_condattr_destroy(sca, status)
        call ferr_abort(test_number, status, " destroying condition attribute structure")
        call fpthrd_condattr_init(sca, status)
        call ferr_abort(test_number, status, " re-initializing condition attribute structure")
        call ftest_comment(test_number, ". Initialized, destroyed and re-initialized condition attribute structure.")

! Set condition attribute to process-private.  Fetched and checked for consistency.
        call fpthrd_condattr_setpshared(sca, FPTHRD_PROCESS_PRIVATE, status)
        call ferr_abort(test_number, status, " setting condition attribute structure")
        scb=sca
        call fpthrd_condattr_getpshared(scb, result, status)
        test_number=test_number+1
        if(result /= FPTHRD_PROCESS_PRIVATE)call ferr_abort(test_number, 1, " wrong condition attribute result")
        call ftest_comment(test_number, ". Set attribute, fetched it and checked for correctness.")

! Get the main thread's ID.  Use it to fetch various settings.
        Call fpthrd_self(base_thread)
        call fpthrd_getschedparam(base_thread, policy, param, status)
        call ferr_abort(test_number, status, " getting default schedule parameters")
        call fpthrd_setschedparam(base_thread, policy, param, status)
        call ferr_abort(test_number, status, " setting default schedule parameters")
        test_number=test_number+1
        call ftest_comment(test_number, ". Retrieved and then set default schedule parameters.")

! Attempt to reset the scheduling scope of a thread to have the highest system-wide
!   priority.  Then return it to the initial setting.
        test_number=test_number+1
        call fpthrd_attr_init(saa, status)
        call ferr_abort(test_number, status, " initializing thread attribute")
        call fpthrd_attr_getscope(saa, scope, status)
        call ferr_abort(test_number, status, " getting thread scope")
! This gives this thread highest priority, system wide.
        call fpthrd_attr_setscope(saa, FPTHRD_SCOPE_SYSTEM, status)

! See if there is permission to set the thread scope.  If not, just quit.
        if(status /= EPERM)THEN

           if(status /= ENOSYS) call ferr_abort(test_number, status, " setting thread scope, system wide")


! Immediately return scope to its default.
           call fpthrd_attr_setscope(saa, scope, status)
           if(status /= ENOSYS)call ferr_abort(test_number, status, " returning thread scope to default setting")

        END IF

        call ftest_comment(test_number, ". Attempted to retrieve, set and reset thread scope.")

! Exercise various schedule and policy routines.
        test_number=test_number+1
        call fpthrd_attr_getinheritsched(saa, inheritsched, status)
        call ferr_abort(test_number, status, " getting inherited schedule")
        call fpthrd_attr_setinheritsched(saa, inheritsched, status)
        call ferr_abort(test_number, status, " setting inherited schedule")
        call ftest_comment(test_number,". Retrieved and then set inherited schedule.")

        test_number=test_number+1
        call fpthrd_attr_getschedparam(saa, param, status)
        call ferr_abort(test_number, status, " getting schedule parameters")
        call fpthrd_attr_getschedparam(saa, param, status)
        call ferr_abort(test_number, status, " setting schedule parameters")
        call ftest_comment(test_number,". Retrieved and then set schedule parameters.")

        test_number=test_number+1
        call fpthrd_attr_getschedpolicy(saa, policy, status)
        call ferr_abort(test_number, status, " getting schedule policy")
        call fpthrd_attr_setschedpolicy(saa, policy, status)
        call ferr_abort(test_number, status, " setting schedule policy")
        call ftest_comment(test_number,". Retrieved and then set schedule policy.")

        test_number=test_number+1
        call fpthrd_mutexattr_getprioceiling(sma, prioceiling, status)
        if(status /= ENOSYS .and. status /= ENOTSUP)&
          call ferr_abort(test_number, status, " getting mutex priority ceiling")

        call fpthrd_mutexattr_setprioceiling(sma, prioceiling, status)
        if(status /= ENOSYS .and. status /= ENOTSUP)&
          call ferr_abort(test_number, status, " setting mutex priority ceiling")
        call ftest_comment(test_number,". Attempted to retrieve and then set mutex priority ceiling.")
        call skip
     END PROGRAM

