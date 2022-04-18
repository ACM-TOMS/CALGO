!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Authors:
!!!  Richard J. Hanson (koolhans@rice.edu)
!!!  Rice University,  Rice Center for High Performance Software Research
!!!
!!!  Clay P. Breshears (clay.breshears@intel.com)
!!!  KAI Software Labs, a division of Intel Americas, Inc.
!!!
!!!  Henry A. Gabb (henry.gabb@intel.com)
!!!  KAI Software Labs, a division of Intel Americas, Inc.
!!!
!!!  Last change:  CPB   Tue Jan 22 11:40:46 CST 2002
!!!  Last change:  RJH  21 Dec 2001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     module fpthrdA
        implicit none

        include "constants.inc"

        type C_NULL
           INTEGER POINTER_VALUE
        end type

        TYPE fpthrd_once_t
          private
          TYPE(FPTHRD_MUTEX_T) MUTEX
          LOGICAL :: FLAG
        END TYPE fpthrd_once_t


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These routines define overloaded assignment for a derived type.
!     I=TYPE(C_NULL)(= value of)

        INTERFACE ASSIGNMENT (=)
           MODULE PROCEDURE INTEQCPTR, CPTREQINT
        END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTERFACE fpthrd_join
        MODULE PROCEDURE fpthrd1_join, fpthrd2_join
      END INTERFACE

      INTERFACE fpthrd_mutex_init ! Second argument can be of two types:
        MODULE PROCEDURE fpthrd1_mutex_init, fpthrd2_mutex_init
      END INTERFACE


      INTERFACE fpthrd_cond_init ! Second argument can be of two types:
        MODULE PROCEDURE fpthrd1_cond_init, fpthrd2_cond_init
      END INTERFACE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This gets the initialization for the condition variable.
        TYPE(fpthrd_cond_t),PARAMETER :: FPTHRD_COND_INITIALIZER=fpthrd_cond_t(cond_I)
! This gets the initialization for the mutex variable.
        TYPE(fpthrd_mutex_t),PARAMETER :: FPTHRD_MUTEX_INITIALIZER=fpthrd_mutex_t(mutex_I)
! This gets the initialization structure for the 'once' functionality.
        TYPE(fpthrd_once_t), PARAMETER :: &
          FPTHRD_ONCE_INIT=fpthrd_once_t(fpthrd_mutex_t(mutex_I), .FALSE.)

! This is the reserved value that denotes a NULL on the C side.
        TYPE(C_NULL) :: NULL = C_NULL(huge(1))
!           NULL%POINTER_VALUE=huge(1)



! This is a table of descriptors for Pthreads error codes.
  character(LEN=12), private :: desc_pterrors(13) = &
(/  'ESRCH       ',  'EINVAL      ',  'EFAULT      ',   'ENOTSUP     ', &
    'EAGAIN      ',  'EDEADLK     ',  'ENOSYS      ',   'EPERM       ', &
    'EBUSY       ',  'ENOMEM      ',  'ETIMEDOUT   ',   'EINTR       ', &
    'ENOSPC      ' /)

! This is a brief meaning clause for each error code.
  character(LEN=32), private :: desc_meanings(13) = &
(/  'No such thread exists.          ',&
    'Invalid argument                ',&
    'Illegal address                 ',&
    'Unsupported option - ignore?    ',&
    'Resource temporarily unavailable',&
    'Program would otherwise deadlock',&
    'Unsupported function            ',&
    'No permission for the operation ',&
    'A "try" function failed         ',&
    'Not enough memory               ',&
    'A time limit was reached.       ',&
    'Interrupted by a signal         ',&
    'No space on a device            '/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CONTAINS


        SUBROUTINE ftest_comment (test_number, text)
           IMPLICIT NONE
           INTEGER test_number, status, I, J
           INTEGER(kind=CINT), PARAMETER :: MESSAGE_SIZE=256  !force 4-bit
           INTEGER(kind=CINT) MESSAGE(MESSAGE_SIZE)           !integer for
                                                              !C code
           CHARACTER(LEN=*)text
           CHARACTER(LEN=MESSAGE_SIZE+20) message_text


           WRITE(*,"('Testing',I3,A)") test_number,trim(text)
           RETURN

        ENTRY ferr_abort (test_number, status, text)
          IF(status == 0) RETURN
! This routine gives integer values of the separate characters in the message.
! The ensuing message may be in a language other than English.
          message=32
          call fpthrd_strerror(status, message, message_size)
          message_text='Unix error summary: '
          If(message(1) == 0) message_text=trim(message_text)//'   None available.'
          J=23
          do I=1,MESSAGE_SIZE
            if(message(i) == 0) exit ! C character strings end with 0.
            message_text(J:J)=char(message(i))
            J=J+1
          end do

          WRITE(*,"('Failed ',I3,' with value',I3,2x,A)") test_number,status,trim(text)
          DO I=1,13
! Write out a matched error message corresponding to the value of status.
            if(fpthrd_errors(I) /= status) CYCLE
            WRITE(*,'(I6,2x,A,1x,A))')&
              fpthrd_errors(I), desc_pterrors(I),trim(desc_meanings(I))
            EXIT
          END DO
          WRITE(*,'(2x,A)') trim(message_text)
          STOP "Abort"

        ENTRY SKIP
           WRITE(*,"(/)")
           RETURN
        END SUBROUTINE

! These subroutines support overloaded assignment used with some of the derived types.

        SUBROUTINE INTEQCPTR (I, S)
           IMPLICIT NONE
           TYPE(C_NULL), INTENT(IN) :: S
           INTEGER, INTENT(INOUT) :: I
           I = s % POINTER_VALUE
        END SUBROUTINE


        SUBROUTINE CPTREQINT (S, I)
           IMPLICIT NONE
           TYPE(C_NULL), INTENT(INOUT) :: S
           INTEGER, INTENT(IN) :: I
           s % POINTER_VALUE=I
        END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some routines have arguments that can be of type C_NULL or another choice.
! These are implemented as GENERIC interfaces.


           SUBROUTINE fpthrd1_join(THREAD, EXITCODE, status)
             TYPE(FPTHRD_T) THREAD
             INTEGER EXITCODE
             INTEGER, OPTIONAL :: status
             INTEGER Lstatus


             call fpthr_join(THREAD, EXITCODE, Lstatus)
             if(present(status)) status=Lstatus
          END SUBROUTINE


          SUBROUTINE fpthrd2_join(THREAD, EXITCODE, status)
             TYPE(FPTHRD_T) THREAD
             TYPE(C_NULL) EXITCODE
             INTEGER, OPTIONAL :: status
             INTEGER Lstatus


             call fpthr_join(THREAD, EXITCODE, Lstatus)
             if(present(status)) status=Lstatus
          END SUBROUTINE



           SUBROUTINE FPTHRD1_MUTEX_INIT(MUTEX, ATTR, STATUS)
             TYPE(FPTHRD_MUTEX_T) MUTEX
             TYPE(FPTHRD_MUTEXATTR_T) ATTR
             INTEGER, OPTIONAL :: STATUS
             INTEGER Lstatus


             CALL FPTHR_MUTEX_INIT(MUTEX, ATTR, Lstatus)
             IF(PRESENT(status)) STATUS=Lstatus
           END SUBROUTINE


           SUBROUTINE FPTHRD2_MUTEX_INIT(MUTEX, ATTR, STATUS)
             TYPE(FPTHRD_MUTEX_T) MUTEX
             TYPE(C_NULL) ATTR
             INTEGER, OPTIONAL :: STATUS
             INTEGER Lstatus


             CALL FPTHR_MUTEX_INIT(MUTEX, ATTR, Lstatus)
             IF(PRESENT(status)) STATUS=Lstatus


           END SUBROUTINE


           SUBROUTINE FPTHRD1_COND_INIT(COND, ATTR, STATUS)
             TYPE(FPTHRD_COND_T) COND
             TYPE(FPTHRD_CONDATTR_T) ATTR
             INTEGER, OPTIONAL :: STATUS
             INTEGER Lstatus


             CALL FPTHR_COND_INIT(COND, ATTR, Lstatus)
             IF(PRESENT(status)) STATUS=Lstatus
           END SUBROUTINE


           SUBROUTINE FPTHRD2_COND_INIT(COND, ATTR, STATUS)
             TYPE(FPTHRD_COND_T) COND
             TYPE(C_NULL) ATTR
             INTEGER, OPTIONAL :: STATUS
             INTEGER Lstatus


             CALL FPTHR_COND_INIT(COND, ATTR, Lstatus)
             IF(PRESENT(status)) STATUS=Lstatus


           END SUBROUTINE


           SUBROUTINE fpthrd_attr_destroy(ATTR, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus


              Call fpthr_attr_destroy(ATTR, Lstatus)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getdetachstate(ATTR, CREATESTATE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              INTEGER, INTENT(OUT) :: CREATESTATE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getdetachstate(ATTR, CREATESTATE, Lstatus)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getinheritsched(ATTR, INHERITSCHED, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              INTEGER, INTENT(OUT) :: INHERITSCHED
              INTEGER,OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getinheritsched(ATTR, INHERITSCHED, Lstatus)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getschedparam(ATTR, PARAM, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              TYPE(FSCHED_PARAM), INTENT(INOUT) :: PARAM
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getschedparam(ATTR, PARAM, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getschedpolicy(ATTR, POLICY, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              INTEGER, INTENT(OUT) :: POLICY
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getschedpolicy(ATTR, POLICY, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getscope(ATTR, SCOPE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              INTEGER, INTENT(OUT) :: SCOPE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getscope(ATTR, SCOPE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_getstacksize(ATTR, STACKSIZE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(IN) :: ATTR
              TYPE(FSIZE_T), INTENT(OUT) :: STACKSIZE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_getstacksize(ATTR, STACKSIZE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_init(ATTR, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_init(ATTR, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setdetachstate(ATTR, DETACHSTATE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: DETACHSTATE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_setdetachstate(ATTR, DETACHSTATE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setinheritsched(ATTR, INHERIT, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: INHERIT
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_setinheritsched(ATTR, INHERIT, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setschedparam(ATTR, PARAM, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              TYPE(FSCHED_PARAM), INTENT(IN) :: PARAM
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              call fpthr_attr_setschedparam(ATTR, PARAM, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setschedpolicy(ATTR, POLICY, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: POLICY
              INTEGER, OPTIONAL  :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_setschedpolicy(ATTR, POLICY, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setscope(ATTR, SCOPE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: SCOPE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_setscope(ATTR, SCOPE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_attr_setstacksize(ATTR, STACKSIZE, STATUS)
              TYPE(FPTHRD_ATTR_T), INTENT(INOUT) :: ATTR
              TYPE(FSIZE_T), INTENT(IN) :: STACKSIZE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_attr_setstacksize(ATTR, STACKSIZE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cancel(THREAD, STATUS)
              TYPE(FPTHRD_T) THREAD
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cancel(THREAD, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_condattr_destroy(ATTR, STATUS)
              TYPE(FPTHRD_CONDATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_condattr_destroy(ATTR, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_condattr_getpshared(ATTR, PSHARED, STATUS)
              TYPE(FPTHRD_CONDATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(OUT) :: PSHARED
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_condattr_getpshared(ATTR, PSHARED, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_condattr_init(ATTR, STATUS)
              TYPE(FPTHRD_CONDATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_condattr_init(ATTR, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_condattr_setpshared(ATTR, PSHARED, STATUS)
              TYPE(FPTHRD_CONDATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: PSHARED
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_condattr_setpshared(ATTR, PSHARED, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cond_broadcast(COND, STATUS)
              TYPE(FPTHRD_COND_T), INTENT(INOUT) :: COND
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cond_broadcast(COND, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cond_destroy(COND, STATUS)
              TYPE(FPTHRD_COND_T), INTENT(INOUT) :: COND
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cond_destroy(COND, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cond_signal(COND, STATUS)
              TYPE(FPTHRD_COND_T), INTENT(INOUT) :: COND
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cond_signal(COND, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cond_timedwait(COND, MUTEX, TIMESPEC, STATUS)
              TYPE(FPTHRD_COND_T), INTENT(INOUT) :: COND
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              TYPE(FTIMESPEC), INTENT(IN) :: TIMESPEC
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cond_timedwait(COND, MUTEX, TIMESPEC, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_cond_wait(COND, MUTEX, STATUS)
              TYPE(FPTHRD_COND_T), INTENT(INOUT) :: COND
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_cond_wait(COND, MUTEX, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_detach(THREAD, STATUS)
              TYPE(FPTHRD_T), INTENT(IN) :: THREAD
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_detach(THREAD, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_getschedparam(THREAD, POLICY, PARAM, STATUS)
              TYPE(FPTHRD_T), INTENT(IN) :: THREAD
              TYPE(FSCHED_PARAM), INTENT(INOUT) :: PARAM
              INTEGER, INTENT(OUT) :: POLICY
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_getschedparam(THREAD, POLICY, PARAM, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_destroy(ATTR, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_destroy(ATTR, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_getprioceiling(ATTR, PRIOCEILING, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(OUT) :: PRIOCEILING
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              call fpthr_mutexattr_getprioceiling(ATTR, PRIOCEILING, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_getprotocol(ATTR, PROTOCOL, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(OUT) :: PROTOCOL
               INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_getprotocol(ATTR, PROTOCOL, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_getpshared(ATTR, PSHARED, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(OUT) :: PSHARED
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_getpshared(ATTR, PSHARED, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_init(ATTR, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_init(ATTR, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_setprioceiling(ATTR, PRIOCEILING, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: PRIOCEILING
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_setprioceiling(ATTR, PRIOCEILING, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_setprotocol(ATTR, PROTOCOL, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: PROTOCOL
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_setprotocol(ATTR, PROTOCOL, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutexattr_setpshared(ATTR, PSHARED, STATUS)
              TYPE(FPTHRD_MUTEXATTR_T), INTENT(INOUT) :: ATTR
              INTEGER, INTENT(IN) :: PSHARED
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutexattr_setpshared(ATTR, PSHARED, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_destroy(MUTEX, STATUS)
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_destroy(MUTEX, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_getprioceiling(mutex, priceiling, status)
              TYPE(FPTHRD_MUTEX_T), INTENT(IN) :: mutex
              INTEGER, intent(OUT) :: priceiling
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_getprioceiling(mutex, priceiling, Lstatus)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_lock(MUTEX, STATUS)
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_lock(MUTEX, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_setprioceiling &
                (mutex, priceiling, oldceiling, status)
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: mutex
              INTEGER, intent(IN) :: priceiling
              INTEGER, intent(INOUT) :: oldceiling
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_setprioceiling &
               (mutex, priceiling, oldceiling, Lstatus)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_trylock(MUTEX, STATUS)
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_trylock(MUTEX, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_mutex_unlock(MUTEX, STATUS)
              TYPE(FPTHRD_MUTEX_T), INTENT(INOUT) :: MUTEX
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_mutex_unlock(MUTEX, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


        SUBROUTINE FPTHRD_ONCE(once_block, once_routine, status)

! Alternative code for the pthread_once() function.  (This version bypasses
! the copy of the compile-time initialization structure to module FPTHRDA.)

        TYPE(FPTHRD_ONCE_T), INTENT(INOUT) :: once_block
! The variable once_block must be set with an assignment statement-
!       once_block = FPTHRD_ONCE_INIT
! before FPTHRD_ONCE() is entered.  Derived type variables FPTHRD_ONCE_T
! consist of (FPTHRD_MUTEX_T, LOGICAL).  This is defined in module FPTHRDA.

! This use of FPTHRD_ONCE_INIT must be made after the call to
! fpthrd_data_exchange(), so that the mutex and LOGICAL of
! FPTHRD_ONCE_INIT are initialized.

! This is the routine that a single thread will call.
          INTERFACE
            subroutine once_routine()
            end subroutine
          END INTERFACE

        INTEGER, OPTIONAL, INTENT(OUT) :: STATUS
        INTEGER local_status

        local_status=0
LOOP:   DO
! An alternate thread already called once_routine().  So exit immediately.
! No need to acquire the mutex.
        IF(VOLATILE_L(once_block % flag)) THEN
          IF(PRESENT(STATUS)) STATUS=local_status
          RETURN
        END IF
! Acquire the mutex and check again.
        call fpthrd_mutex_lock(once_block % mutex, local_status)
! If there is an exception in locking the mutex, exit.
          IF(local_status /= 0) EXIT LOOP
! A thread may have missed the first check and then acquire the mutex.
! But the call to once_routine was already made.
! The mutex will then be unlocked after the loop is exited.
          IF(VOLATILE_L(once_block % flag)) EXIT LOOP
! Make the call to once_routine(), set the flag.
          call once_routine()
          once_block % flag=.TRUE.
! Exit the loop, which unlocks the mutex.
          EXIT LOOP
        END DO LOOP

! Release the mutex.  Only threads that had it locked will release.
! The value of local_status (exception flag) will be optionally returned to
! the calling routine as status.  An exception in locking or unlocking
! the mutex will cause local_status to have a non-zero value.
        IF(local_status == 0) &
          call fpthrd_mutex_unlock(once_block % mutex, local_status)
        IF(PRESENT(STATUS)) STATUS=local_status
        RETURN

        END SUBROUTINE

        FUNCTION VOLATILE_L(FLAG)
          IMPLICIT NONE
! This function avoides code optimization that may result in errors
! with multiple threads.  The routine is called from FPTHRD_ONCE().
          LOGICAL VOLATILE_L
          LOGICAL, INTENT(INOUT) :: FLAG
          FLAG=.NOT. (.NOT. FLAG)
          VOLATILE_L=FLAG
        END FUNCTION

           SUBROUTINE fpthrd_setcancelstate(STATE, OLDSTATE, STATUS)
              INTEGER, INTENT(IN) :: STATE
              INTEGER, INTENT(OUT) :: OLDSTATE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_setcancelstate(STATE, OLDSTATE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_setcanceltype(TYPE, OLDTYPE, STATUS)
              INTEGER, INTENT(IN) :: TYPE
              INTEGER, INTENT(OUT) :: OLDTYPE
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_setcanceltype(TYPE, OLDTYPE, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE fpthrd_setschedparam(THREAD, POLICY, PARAM, STATUS)
              TYPE(FPTHRD_T), INTENT(INOUT) :: THREAD
              TYPE(FSCHED_PARAM), INTENT(IN):: PARAM
              INTEGER, INTENT(IN) :: POLICY
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_setschedparam(THREAD, POLICY, PARAM, LSTATUS)
              If(present(status))status=Lstatus
           END SUBROUTINE


           SUBROUTINE FPTHRD_setconcurrency(NTHREADS, STATUS)
              INTEGER, INTENT(IN) :: NTHREADS
              INTEGER, OPTIONAL :: STATUS
              INTEGER Lstatus
              Call fpthr_setconcurrency(NTHREADS, Lstatus)
              If(present(status)) status=Lstatus
           END SUBROUTINE

     end module


     MODULE FPTHRD
        USE fpthrdA
        INTERFACE

           SUBROUTINE fpthrd_equal(THREAD1, THREAD2, FLAG)
              USE fpthrdA, ONLY : FPTHRD_T
              TYPE(FPTHRD_T), INTENT(IN) :: THREAD1, THREAD2
              INTEGER, INTENT(OUT) :: FLAG
           END SUBROUTINE


           SUBROUTINE fpthrd_self(THREAD)
              USE fpthrdA, ONLY : FPTHRD_T
              TYPE(FPTHRD_T), INTENT(OUT) :: THREAD
           END SUBROUTINE


           SUBROUTINE fpthrd_testcancel()
           END SUBROUTINE

           SUBROUTINE fpthrd_getconcurrency(NTHREADS)
             INTEGER, INTENT(OUT) :: NTHREADS
           END SUBROUTINE

           SUBROUTINE fpthrd_set_ftimespec(change_sec, &
             Change_nanosec, waittime)
             USE fpthrdA, ONLY : FTIMESPEC
              INTEGER, INTENT(IN) :: change_sec, change_nanosec
              TYPE(FTIMESPEC), INTENT(INOUT) :: waittime
           END SUBROUTINE

           SUBROUTINE fpthrd_set_fsched_param(schedule_value, param)
             USE fpthrdA, ONLY: FSCHED_PARAM
             INTEGER, INTENT(IN) :: schedule_value
             TYPE(FSCHED_PARAM), INTENT(INOUT) :: param
           END SUBROUTINE

           SUBROUTINE fpthrd_get_fsched_param(schedule_value, param)
             USE fpthrdA, ONLY: FSCHED_PARAM
             INTEGER, INTENT(OUT) :: schedule_value
             TYPE(FSCHED_PARAM), INTENT(IN) :: param
           END SUBROUTINE

           SUBROUTINE fpthrd_set_fsize(size_value, size)
             USE fpthrdA, ONLY: FSIZE_T
             INTEGER, INTENT(IN) :: size_value
             TYPE(FSIZE_T), INTENT(INOUT) :: size
           END SUBROUTINE

           SUBROUTINE fpthrd_get_fsize(size_value, size)
             USE fpthrdA, ONLY: FSIZE_T
             INTEGER, INTENT(OUT) :: size_value
             TYPE(FSIZE_T), INTENT(IN) :: size
           END SUBROUTINE

        END INTERFACE

     end module
