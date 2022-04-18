 program thrdtmessy
! Time-stamp:  <2015-12-10 19:30:15 m>
! This code almost entirely written by Richard J. Hanson.
! Small test of threaded use of messy.  Output is to ?_err#, ?_mes#, and ?_odd#,
! where ? is s, d, or q depending on the precision, # is 0, 1, ... numt_-1, with
! only odd numbers for ?_odd#, and only even numbers for ?_err#, and ?_mes#.
  use messy_, only: messy_ty, messy, rk
  use omp_lib
  implicit none
  integer, parameter :: nlen=1+int(log10(real(numt_))) ! Needed for numbers
  character (len=*), parameter :: root="Drivers/NewResults/"
  type(messy_ty), allocatable :: et(:)
  character(len=nlen) :: buf
  integer :: i, j, nt, tn
  real(rk) :: fdat(0:numt_+3)
  ! buf  Used for creating character representation of integers
  !      and lines of working character strings
  ! nt   Max number of threads in use in a parallel section
  ! tn   The thread index used in a parallel do section
  nt = numt_
  call omp_set_num_threads(nt)
  allocate(et(0:nt-1))

  ! Assign units, one for each thread
#ifdef MSWin
! On Windows, remove results directory and then create it again.
! Contents of results directory are deleted and then written.
   call system("rd /s /q ")  !remove, via dos commands
   call system("mkdir ")    !create
#endif

  do i = 0, nt - 1
    buf = ' '
    write(buf,'(I0)') i
    et(i)%ename = "thrdtmessy_" // trim(buf)
    if (mod(i,2) == 1) then ! To show two ways of arranging outputs
      open(newunit=et(i)%munit,&
        & file= root // "odd" // trim(buf) // "." // plet_, status='REPLACE')
      et(i)%eunit = et(i)%munit
    else
      open(newunit=et(i)%munit, &
        & file= root // "mes" // trim(buf)//"." // plet_, status='REPLACE')
      open(newunit=et(i)%eunit, &
        & file= root // "err" // trim(buf)//"." // plet_, status='REPLACE')
    end if
    et(i)%fpprec = min(32,et(i)%fpprec) ! So results with NAG compiler match
  end do
  fdat(0:nt+3)=[(-1.9375+real(j,rk), j = 0, nt+3)]
  !Execute a parallel do, with OpenMP.  Each thread writes a simple message.
  !Every third value of the loop index writes an exception or "error" message.
  ! Note that OpenMP types the loop index i as THREADPRIVATE, by default.
! Below can use: static, dynamic, guided, runtime, or auto (only static checked)
!$omp parallel do schedule(static),private(tn)
  do i=0,3*nt-1 ! Upper limit here is arbitrary > 0.
    tn=omp_get_thread_num()
    call messy(et(tn),'$L80$D6messy says Hi from Thread=$I,&
      & with loop index $I, and reals:$N$V',idat=[tn,i],&
      rdat=fdat(0:tn+3))
    if (mod(i,3) == 1) then
      call messy(et(tn),"$E14Same for errors from Thread=$I,&
      & with loop index $I,&
        & more $N$V", idat=[99,tn,i], rdat=fdat(0:tn+3))
    end if
  end do
!$omp end parallel do
   do i = 0, nt-1 ! Close all the output files
     if (et(i)%munit /= et(i)%eunit) close(et(i)%eunit)
     close(et(i)%munit)
   end do
#ifdef MSWin
! On Windows, run the batch file checkthrdtmessy.bat that compares
! output in directory results with the distribution files.
! The results of the comparison are in d_summary.
   call system("echo off")
   call system("checkthrdtmessy > d_summary")
   call system("type d_summary")
#endif

end program thrdtmessy
