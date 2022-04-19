!     Last change:  CPB  22 JAN 2002   11:12 am
!     Last change:  RH   25 JAN 2002    2:22 pm
! Benchmarking elementary use of Pthreads.  Print summary of tests.
! This code is part of the package "A Fortran Interface to Posix Threads,"  to be
! published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
! This is bench1.f90.

     MODULE global_test2
        USE fpthrd
        IMPLICIT NONE
! These are control values used in the tests.
        INTEGER test_number, thread_counter
! This is the number of threads launched at one time.
        INTEGER, PARAMETER :: NTHREADS=1023
        INTEGER, PARAMETER :: NTRIES = 16
        INTEGER id_no(NTHREADS)

        TYPE(FPTHRD_t) test_id(NTHREADS), first_id
        TYPE(FPTHRD_attr_t) attribute, start_attr
        TYPE(FPTHRD_mutex_t) cancel_mutex, recursion_mutex, sum_mutex
        TYPE(FPTHRD_once_t) :: testing_once=FPTHRD_ONCE_INIT
        TYPE(FPTHRD_cond_t) launch_done
        TYPE(FSCHED_PARAM) param
        TYPE(fsize_t) stacksize


! Declare a derived type that will carry the problem dope.
! These are the arguments of a dot product function.

        TYPE FUNCTION_ARGUMENTS
           INTEGER N
           REAL, DIMENSION(:,:), POINTER :: SMATRIX
           REAL, DIMENSION(:), POINTER :: SX
           INTEGER INCY
           REAL, dimension(:), POINTER :: SY
           INTEGER ROW_INDEX
        END TYPE

        TYPE(FUNCTION_ARGUMENTS) inb(NTHREADS), INA, INC

        interface

           SUBROUTINE FSGEMV (N, SMATRIX, SX, SY, ROW_START, ROW_END)
              IMPLICIT NONE
              INTEGER N, ROW_START, ROW_END
              REAL, POINTER, DIMENSION(:) :: SMATRIX(:,:), SX, SY
           END SUBROUTINE
           SUBROUTINE TEST6(NUMBER)
             INTEGER NUMBER
           END SUBROUTINE
           SUBROUTINE ONCE6()
           END SUBROUTINE
           RECURSIVE SUBROUTINE test7(LIMITS)
              INTEGER LIMITS(2)
           END SUBROUTINE
           RECURSIVE SUBROUTINE test8(LIMITS)
              INTEGER LIMITS(2)
           END SUBROUTINE

        end interface

     END MODULE

     SUBROUTINE FSGEMV (N, SMATRIX, SX, SY, COL_START, COL_END)
        USE global_test2, only : sum_mutex, fpthrd_mutex_unlock, fpthrd_mutex_lock
        IMPLICIT NONE
        INTEGER N, I, J, L, COL_START, COL_END, status
        REAL, pointer :: SMATRIX(:,:), SX(:)
        REAL, pointer :: SY(:)
        REAL T(N)
        IF(COL_START <= 0) RETURN
        IF(COL_END > N) RETURN
        IF(COL_END < COL_START) RETURN
        T=0E0
        L=IAND(COL_END-COL_START+1,7)
         DO J=0,L-1
           DO I=1,N
              T(I)=T(I)+smatrix(I,COL_START+J)*SX(COL_START+J)
           END DO
        END DO

        DO J=COL_START+L, COL_END, 8
            DO I=1,N
              T(I)=T(I)+smatrix(I,J)*SX(J) +smatrix(I,J+1)*SX(J+1) +smatrix(I,J+2)*SX(J+2) &
              +smatrix(I,J+3)*SX(J+3) +smatrix(I,J+4)*SX(J+4)+smatrix(I,J+5)*SX(J+5) &
              +smatrix(I,J+6)*SX(J+6)+smatrix(I,J+7)*SX(J+7)
           END DO
        END DO

        call fpthrd_mutex_lock(sum_mutex, status)
        DO I=1,N
           SY(I)=SY(I)+T(I)
        END DO
        call fpthrd_mutex_unlock(sum_mutex, status)
     END SUBROUTINE

     SUBROUTINE FSDSDOT (N, SMATRIX, SX, INCY, SY, ROW_INDEX)

        IMPLICIT NONE

        INTEGER N, INCY, J, Q, ROW, ROW_INDEX
        REAL, pointer :: SMATRIX(:,:), SX(:)
        REAL, pointer :: SY(:)
        DOUBLE PRECISION T
        ROW=ROW_INDEX
        T=0D0
        Q=1

        DO J=1,N
           T=T+SMATRIX(ROW, J)*SX(Q)
           Q=Q+INCY
        END DO
        SY(ROW)=T
     END SUBROUTINE

     subroutine once6()
        use global_test2, dummy => once6
        implicit none

! Use with FPTHRD_once() to initialize any data element.
        integer status
        recursion_mutex=FPTHRD_MUTEX_INITIALIZER
        sum_mutex =FPTHRD_MUTEX_INITIALIZER
        test_number=test_number+1
        call ftest_comment(test_number, ". Once initializing a mutex that is used later.")
     end subroutine


     subroutine test6(arg_in)
        use global_test2, dummy => test6
        implicit none

        INTERFACE
           SUBROUTINE FSDSDOT (N, SMATRIX, SX, INCY, SY, ROW_INDEX)
              INTEGER N, INCY, ROW_INDEX
              REAL, POINTER, DIMENSION(:) :: SMATRIX(:,:), SX, SY
           END SUBROUTINE
        END INTERFACE

        integer status, value, arg_in, I
        call FPTHRD_once (testing_once, once6, status)
        call ferr_abort(test_number, status, "starting once function for thread data initialization")
        I=arg_in

! Use the structure holding the arguments for this function.
        call fsdsdot (inb(I) % n, inb(I) % SMATRIX, inb(I) % SX, inb(I) % incy, inb(I) % SY, inb(I)% row_INDEX)
     end subroutine

     recursive subroutine test7(LIMITS)
        use global_test2, my_test=>test7 ! Symbol is reset to a dummy. Avoids conflict.
        implicit none

        INTEGER :: IDEAL, K, J, LIMITS(2), LIMITS_L(2), LIMITS_R(2), status
        TYPE(FPTHRD_T) THREAD_L, THREAD_R

        call FPTHRD_mutex_lock(recursion_mutex, status)
        call ferr_abort(test_number, status, " locking mutex in test7")
        K=LIMITS(2)-LIMITS(1)+1;J=(LIMITS(1)+LIMITS(2))/2

! The problem limits are split into two equally sized groups.
        LIMITS_L=LIMITS;LIMITS_R=LIMITS
        LIMITS_L(2)=J;LIMITS_R(1)=J+1

        call FPTHRD_mutex_unlock(recursion_mutex, status)
        call ferr_abort(test_number, status, " unlocking mutex in test7")
        IDEAL=(INA%N+7)/8
        IF(K <= IDEAL) THEN
! This is where the work actually gets done.  The above value of IDEAL is
! problem dependent.
           CALL FSGEMV (INA % N, INA % SMATRIX, INA % SX, INA % SY, LIMITS(1), LIMITS(2))
        ELSE
           call FPTHRD_create(THREAD_L, attribute, my_test, limits_l, status)

           call ferr_abort(test_number, status, " recursive create-L in test7")
           call FPTHRD_create(THREAD_R, attribute, my_test, limits_r, status)
           call ferr_abort(test_number, status, " recursive create-R in test7")
           call FPTHRD_join(THREAD_L, NULL, status)
           call ferr_abort(test_number, status, " recursive join-L in test7")
           call FPTHRD_join(THREAD_R, NULL, status)
           call ferr_abort(test_number, status, " recursive join-R in test7")
        END IF

     END SUBROUTINE
     recursive subroutine test8(LIMITS)
        use global_test2, my_test=>test8 ! Symbol is reset to a dummy. Avoids conflict.
        implicit none

        INTEGER :: IDEAL, I, J, K, LIMITS(2), LIMITS_L(2), LIMITS_R(2), status
        TYPE(FPTHRD_T) THREAD_L, THREAD_R

        call FPTHRD_mutex_lock(recursion_mutex, status)
        call ferr_abort(test_number, status, " locking mutex in test7")
        K=LIMITS(2)-LIMITS(1)+1;J=(LIMITS(1)+LIMITS(2))/2

! The problem limits are split into two equally sized groups.
        LIMITS_L=LIMITS;LIMITS_R=LIMITS
        LIMITS_L(2)=J;LIMITS_R(1)=J+1
        call FPTHRD_mutex_unlock(recursion_mutex, status)
        call ferr_abort(test_number, status, " unlocking mutex in test7")
        IDEAL=(INA%N+7)/8

        IF(K <= IDEAL) THEN
! This is where the work actually gets done.  The above value of IDEAL is
! problem dependent.
           DO I=1,INA%N
              DO J=max(1,limits(1)),min(ina%n,limits(2))
                 INC % SMATRIX(I,J)=INA % smatrix(J,I)
              END DO
           END DO

        ELSE
           call FPTHRD_create(THREAD_L, attribute, my_test, limits_l, status)
           call ferr_abort(test_number, status, " recursive create-L in test7")
           call FPTHRD_create(THREAD_R, attribute, my_test, limits_r, status)
           call ferr_abort(test_number, status, " recursive create-R in test7")
           call FPTHRD_join(THREAD_L, NULL, status)
           call ferr_abort(test_number, status, " recursive join-L in test7")
           call FPTHRD_join(THREAD_R, NULL, status)
           call ferr_abort(test_number, status, " recursive join-R in test7")
        END IF

     END SUBROUTINE


     program benchmain1
        USE global_test2
        implicit none

        REAL, POINTER :: matrix_a(:,:), matrix_b(:,:), vector(:),y_serial(:),y_thread(:)
        REAL errnorm, norm, temp

        integer create_state, i, j, k, L, n, status, value, LIMITS(2), schedule_value
        INTEGER TIMEE(3), TIMES(3)

        call FPTHRD_setconcurrency(NTHREADS+1, STATUS)
! Initialize mutex and condition variable.
        cancel_mutex=FPTHRD_MUTEX_INITIALIZER
        launch_done =FPTHRD_COND_INITIALIZER

! Test 6. Start a peer function that calls a library function.  Arguments
!           are packed into a derived type.
        N=NTHREADS

! No claims of randomness are made for this sequence.  It is used to generate
! a non-repeatable sequence of matrix and vector values.
        ALLOCATE(MATRIX_A(N,N), VECTOR(N), Y_SERIAL(N), Y_THREAD(N))
        call random_number(matrix_a)
        call random_number(vector)
        DO J=1,N
           vector(J)=2e0*vector(J)-1e0
        END DO

! This value for stacksize is used only for a test.
        call fpthrd_set_fsize(64000, stacksize)
        call FPTHRD_attr_init(attribute, status)
        call FPTHRD_attr_setstacksize(attribute, stacksize, status)

! Use threads to compute each row (times) vector concurrently.
        DO I=1,N
! This is a typical situation:  Use a structure to pack all arguments into one
! object.  Then pass the peer code that object.  It unpacks the structure to get
! the arguments.

           inb(I) % n=n;inb(I) % SMATRIX=> matrix_a;
           inb(I) % SX=>vector;inb(I) % incy=1;inb(I) % SY=>y_thread;
           inb(I) % row_index=I
           ID_NO(I)=I

! Create the set of threads for the product.
           call FPTHRD_create (test_id(i), attribute, test6, ID_NO(I), status)
           call ferr_abort (test_number, status, " creating threads")
        END DO

        DO I=1,NTHREADS
           call FPTHRD_join(test_id(i), NULL, status)
           call ferr_abort(test_number, status, " joining a single thread")
        END DO

! Compute the matrix-vector product for comparison.
        CALL SYSTEM_CLOCK(TIMES(1))
        DO L=1,NTRIES
           y_serial=0e0
         do j=1,n
          do i=1,n
            y_serial(i)=y_serial(i)+matrix_a(I,j)*vector(j)
          end do
         end do
        END DO
        CALL SYSTEM_CLOCK(TIMEE(1))

        CALL SYSTEM_CLOCK(TIMES(2))
        DO J=1,NTRIES
! Time the intrinsic matrix-vector multiply, provided with Fortran 90.
           y_serial=matmul(matrix_a, vector)
        END DO
        CALL SYSTEM_CLOCK(TIMEE(2))

        errnorm=sum((y_serial-y_thread)**2)
        norm=sum(y_serial**2)

! The results are correct even if they do not completely agree.
! This test will fail only with a blunder.  Small relative errors will be allowed.
        test_number=test_number+1

        if(errnorm <= EPSILON(NORM)*norm) THEN
           call ftest_comment (test_number, ". Matrix-vector product with each entry using a separate thread.")
        else
           status=1
           call ferr_abort(test_number,status," Serial and threaded matrix-vector product gave different results")
        END IF
        call skip
! Test 6 is complete. */

! Use threads to compute each row (times) vector concurrently.
! This exericse uses one thread to call a routine.  The routine uses divide and
! conquer  (recursion) to reduce the problem size to one with 'good' properties.

        ina % N=n;ina % SMATRIX=> matrix_a;ina % SX=>vector;ina % SY=>y_thread

        call FPTHRD_attr_setinheritsched(attribute, FPTHRD_EXPLICIT_SCHED, status)
        call ferr_abort(test_number, status, " setting inherit schedule")

        call FPTHRD_attr_setschedpolicy(attribute, FSCHED_FIFO, status)
! This may not be supported:
        IF(status /= ENOTSUP) call ferr_abort(test_number, status, " setting schedule policy")
        call FPTHRD_attr_setinheritsched(attribute, FPTHRD_INHERIT_SCHED, status)
! This may not be supported:
        IF(status /= ENOTSUP) call ferr_abort(test_number, status, " setting inherit schedule")
        schedule_value=50
        call fpthrd_set_fsched_param(schedule_value, param)
        call fpthrd_attr_setschedparam(attribute, param, status)
        call ferr_abort(test_number, status, " setting schedule parameter")

! Create the thread for the product.
        CALL SYSTEM_CLOCK(TIMES(3))
        DO J=1,NTRIES
           LIMITS=(/1,N/)
           y_thread=0e0
           call FPTHRD_create (first_id, attribute, test7, LIMITS, status)
           call ferr_abort (test_number, status, " creating a recursive thread")

           call FPTHRD_join(first_id, NULL, status)
           call ferr_abort(test_number, status, " joining a single recursive thread")
        END DO
        CALL SYSTEM_CLOCK(TIMEE(3))

! Check results for correctness:
        errnorm=sum((y_serial-y_thread)**2)
        norm=sum(y_serial**2)
! The results are correct even if they do not completely agree.
! This test will fail only with a blunder.  Small relative errors will be allowed.
        test_number=test_number+1

        if(errnorm <= EPSILON(NORM)*norm) THEN
           call ftest_comment (test_number, ". Matrix-vector product with recursive threads.")
        else
           status=1
           call ferr_abort(test_number,status," Serial and threaded matrix-vector product gave different results")
        END IF
        call skip
! Test 7 is complete.
        call system_clock(count_rate=I)
        write(*,*) "    Matrix size, repeats     ", NTHREADS, NTRIES
        write(*,*) " 0. Clock Rate, Ticks per S. ", I
        write(*,*) " 1. Ordinary product time =  ", TIMEE(1)-TIMES(1)
        write(*,*) " 2. Matmul product time =    ", TIMEE(2)-TIMES(2)
        write(*,*) " 3. THREADED product time =  ", TIMEE(3)-TIMES(3)

        write(*,*) "Matrix-vector Time ratios: [1.]/[3.]: ", &
          REAL(TIMEE(1)-TIMES(1))/REAL(TIMEE(3)-TIMES(3))," [2.]/[3.]: ", REAL(TIMEE(2)-TIMES(2))&
          /REAL(TIMEE(3)-TIMES(3))

        ALLOCATE(matrix_b(N,N))

        CALL SYSTEM_CLOCK(TIMES(1))
        DO K=1,NTRIES
           DO J=1,N
              DO I=1,N
                 matrix_b(j,i)=matrix_a(i,j)
              end do
           end do
        END DO
        CALL SYSTEM_CLOCK(TIMEE(1))


        CALL SYSTEM_CLOCK(TIMES(2))
        DO J=1,NTRIES
           matrix_b=transpose(matrix_a)
        END DO
        CALL SYSTEM_CLOCK(TIMEE(2))

        inc % SMATRIX=> matrix_b
        CALL SYSTEM_CLOCK(TIMES(3))
        DO J=1,NTRIES
           LIMITS=(/1,N/)
           y_thread=0e0
           call FPTHRD_create (first_id, attribute, test8, LIMITS, status)
           call ferr_abort (test_number, status, " creating a recursive thread")

           call FPTHRD_join(first_id, NULL)
           call ferr_abort(test_number, status, " joining a single recursive thread")
        END DO
        CALL SYSTEM_CLOCK(TIMEE(3))
        call skip
        write(*,*) " 1. Ordinary transpose time = ",                    TIMEE(1)-TIMES(1)
        write(*,*) " 2. Intrinsic function TRANSPOSE() time = ",        TIMEE(2)-TIMES(2)
        write(*,*) " 3. THREADED panel transpose time = ",              TIMEE(3)-TIMES(3)
        write(*,*) "Transpose Time ratios: [1.]/[3.]: ", REAL(TIMEE(1)-TIMES(1))&
         /REAL(TIMEE(3)-TIMES(3))," [2.]/[3.]: ", REAL(TIMEE(2)-TIMES(2))/REAL(TIMEE(3)-TIMES(3))

     end program


