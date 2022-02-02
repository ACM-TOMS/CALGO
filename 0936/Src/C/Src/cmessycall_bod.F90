! This code almost entirely written by Richard J. Hanson.
! This is a wrapper routine for C users to use messy() with their
! C-style optional argument interface, cmessy().
  USE ISO_C_BINDING, only: C_LOC, C_NULL_CHAR, C_CHAR, C_INT, C_DOUBLE,&
    & C_FLOAT, C_LONG_DOUBLE, C_PTR, C_F_POINTER

  IMPLICIT NONE
! These define the maximum allowed sizes for C routines to use as
! string arguments for the cmessy text, cmessy insert or ptext, and the
! name of the file to open for output of messages and error messages.
  integer, parameter, private :: max_slen=4096, max_plen=80, max_flen=128

! This derived type functions to transform C addresses of assumed size arrays
! into Fortran assumed shape arrays, required by messy().  An array of this type
! is allocated by the C routine calling Fortran routine
! "allocate_cmessy_interface".  Its size is the maximum number of threads
! created within OpenMP parallel sections in C code that call "cmessy()" within
! that section.  If not using OpenMP this size is =1.
  type, private :: threadargs
    character(len=max_slen) C
    integer, pointer :: ip(:)
    integer, pointer :: iap(:,:)
    real(rk), pointer :: dp(:)
    real(rk), pointer :: dap(:,:)
    complex(rk), pointer :: zp(:)
    complex(rk), pointer :: zap(:,:)
    integer, pointer :: ixp(:)
    character(len=max_plen) P
  end type threadargs

! This C structure essentially matches the public components of messy_ty.
! Differences of Fortran character data and C string data are accounted for by
! conversion of a string to a character variable.  The components "cinit" and
! "cstruct" are used to note initialization and the index of a copy of messy_ty
! to use in calls to messy from C or C++.
  integer, parameter :: numdig = ceiling(-log10(epsilon(1.0_rk)))
  type, bind(c) :: CMESSY_TY_
    character(C_CHAR) :: ename(32)! Name printed in error messages - defined at initialization
    integer(c_int) :: fpprec = numdig ! Default for floating point precision
    integer(c_int) :: kdf = numdig ! Current default real precision.
    integer(c_int) :: line_len = 128 ! Default for line length
    integer(c_int) :: munit = OUTPUT_UNIT ! Message unit number
    integer(c_int) :: eunit = OUTPUT_UNIT ! ERROR_UNIT mixes up output with piping
    integer(c_int) :: maxerr = 0 ! Max value of 1000 * (10*stop + print) + |index|
    integer(c_int) :: lstop = 3 ! Stop indexes <= this do not stop
    integer(c_int) :: lprint = 3 ! Print indexes <= this do not print
    integer(c_int) :: errcnt = 0 ! Count of the number of error messages, incremented
!                       by 1000000 for internal errors inside messy.
    integer(c_int) :: dblev = 3 ! If 0, an immediate return is made (unless text
!                   starts with "$E"), else a $K<integer(c_int)> will behave
!                   as if reaching the end of text if <integer(c_int)> is > dblev.
    integer(c_int) :: cinit = 1234565 ! Used as signature within cmessy() that the struct is initialized
    integer(c_int) :: cthread=0 ! Used for thread index of an OpenMP thread when calling cmessy()
    integer(c_int) :: cstruct   ! Used to pick out copy CE of type messy_ty, array SE(:), used to call messy
  end type CMESSY_TY_ 

  
! Stubs for arrays passed as arguments to messy().
! The C code may not have defined them.
  integer, target :: idat(1), imat(1,1)
  real(rk), target :: rdat(1), rmat(1,1)
  complex(rk), target :: zdat(1), zmat(1,1)
  
! The contents hold pointers to assumed shape arrays
! created for passing to messy(). 
! The allocated size here must be no smaller than the
! number of parallel OpenMP threads used in a C program.
! It is set by the C routine allocate_cmessy_interface.
  type(threadargs),allocatable,private:: TA(:)

! The allocated limit here must be no smaller than the
! number of different C structures cmessy_ty used in
! C or C++ applications.
  type(messy_ty),allocatable,private,save :: SE(:)

! For identifying what local copy of messy_ty type to associate with the
! input argument cmessy_ty.
  integer, private, save :: maxstructs=0, maxthreads=1
! A temp used for error messages, if C interfacing has one.
  type(messy_ty), save :: U
  integer, private, save :: cindex=0
  pdeft_
  
CONTAINS
  SUBROUTINE CALLMESSY_(E, STRING, slen, ND, D,&
    & MDA, NDA, DA,&
    & NZ, Z,&
    & MZA, NZA, ZA,&
    & NI, I,&
    & MIA, NIA, IA,&
    & NX, IX,&
    & PT, plen) bind(C)
! Full use of messy() from a C code, using cmessy().  All supported
! rank-1,-2 double, complex, integer and character arrays are printed.
! This is also suitable for error messages.


    type(CMESSY_TY_), intent(inout) :: E
    character(c_char), intent(in) :: string(*)
    integer, intent(in), value :: nd, mda, nda, nz, mza, nza,&
      & ni, mia, nia, nx, plen, slen
    real(CTYP_), target, intent(in) :: d(*), da(*)
    complex(CTYP_), target, intent(in) :: z(*), za(*)
    integer(c_int), target, intent(in)    :: i(*), ia(*), ix(*)
    character(c_char), target, intent(in) :: pt(*)
    integer, save :: j, tn, ce

#if numt_ > 1
!$OMP THREADPRIVATE(j,tn,ce)
!$OMP CRITICAL
#endif
!----------------------------------------------------------------------

! It determines the entry of ta(:) used for calls to messy.
! With a single thread in C code.
    if(.not. (ALLOCATED(TA).and.ALLOCATED(SE))) then
      U%ENAME=QCALLMESSY_
      call messy(U,"$E36Before any use of $Pmessy(), you first must execute &
        &'allocate_$Pmessy_interface(maxcstructs, maxcthreads).'",&
        & idat=[01], ptext=T_)
    end if
    tn=E % cthread ! The C OpenMP thread number is an argument from cmessy.
    if(tn > maxthreads) then
        call messy(U,"$E36You must first increase maxctreads to &
          &account for maximum number of OpenMP threads in 'allocate_$P&
          & messy_interface (maxcstructs, maxcthreads).'",idat=[02], ptext=T_)
    end if
    CE = E % CSTRUCT ! Index used to pick the local copy of messy_ty.
    SE(CE)%ENAME=' '
    DO J=1,32 ! TRANSFORM C NAME (as a string) FROM INCOMING STRUCT TO FORTRAN NAME (as character).
      IF(E%ENAME(J)==C_NULL_CHAR) EXIT
      SE(CE)%ENAME(J:J)=E%ENAME(J)
    END DO
! Transfer contents of incoming C struct to Fortran messy_ty type.
    SE(CE)%FPPREC = E%FPPREC; SE(CE)%KDF = E%KDF;       SE(CE)%LINE_LEN = E%LINE_LEN
    SE(CE)%MUNIT  = E%MUNIT;  SE(CE)%EUNIT = E%EUNIT;   SE(CE)%MAXERR= E%MAXERR
    SE(CE)%LSTOP =  E%LSTOP;  SE(CE)%LPRINT = E%LPRINT; SE(CE)%ERRCNT=E%ERRCNT
    SE(CE)%DBLEV =  E%DBLEV

    if (slen > max_slen) then
      call messy(se(ce), "$E99In a call to $Pmessycall from $Pmessy, slen=$I&
        & > max_slen=$I, either decrease slen, or increase max_slen in&
        & $Pmessycall_m.F90.", idat=[03, slen, max_slen], ptext=T_)
    end if
! Transform C assumed size arrays to Fortran assumed shape arrays. 
    DO J=1,slen ! Transform C argument string to Fortran characters
      ta(tn)%C(J:J) = STRING(J)
    END DO


! Move data from C assumed size arrays to Fortran assumed shape arrays.
    ta(tn)%DP=>rdat ! Setting this deals with a bug on an unnamed compiler
    if (nd > 0) THEN
      CALL C_F_POINTER(C_LOC(D), ta(tn)%DP, [ND]) ! Creates rank1 real assumed shape array
    END IF

    ta(tn)%DAP=>rmat ! Setting this deals with a bug on an unnamed compiler
    if (mda > 0) THEN
      CALL C_F_POINTER(C_LOC(DA), ta(tn)%DAP, [mda, nda]) ! Creates rank2 real assumed shape array
    END IF


    ta(tn)%ZP=>zdat ! Setting this deals with a bug on an unnamed compiler
    if (nz > 0) THEN
        CALL C_F_POINTER(C_LOC(Z), ta(tn)%ZP, [NZ]) ! Creates rank1 complex assumed shape array
    END IF
    ta(tn)%ZAP=>zmat ! Setting this deals with a bug on an unnamed compiler
    if (mza > 0) THEN
      CALL C_F_POINTER(C_LOC(ZA), ta(tn)%ZAP, [mza, nza]) ! Creates rank2 complex assumed shape array
    END IF
    ta(tn)%IP=>idat ! Setting this deals with a bug on an unnamed compiler
    if (ni > 0) THEN
      CALL C_F_POINTER(C_LOC(I), ta(tn)%IP, [ni]) ! Creates rank1 integer
! assumed shape array
    END IF

    ta(tn)%IAP=>imat ! Setting this deals with a bug on an unnamed compiler
    if (mia > 0) THEN
      CALL C_F_POINTER(C_LOC(IA), ta(tn)%IAP, [mia, nia]) ! Creates rank2 integer assumed shape array
    end if

    if (plen > max_plen) then
      call messy(se(ce), "$E99In a call to callmessy from cmessy, plen=$I >&
        & max_plen=$I, either decrease slen, or increase max_slen in&
        & c2fmessy_m.f90.", idat=[04, plen, max_plen])
    end if
    DO J=1,plen ! Convert inserted text string to Fortran character
      ta(tn)%P(J:J)=PT(J)
    END DO

    ta(tn)%IXP=>idat ! Setting this deals with a bug on an unnamed compiler
    if(nx > 0) then ! messy will get confused on errors if ix is present,
      CALL C_F_POINTER(C_LOC(IX), ta(tn)%IXP, [NX]) ! but not used.
      call messy(se(ce),ta(tn)%C(1:slen),ta(tn)%ip,ta(tn)%DP,ta(tn)%iap,ta(tn)%dap,&
        & ta(tn)%zp,ta(tn)%zap,ta(tn)%ixp,ta(tn)%p(1:plen))
    ELSE
      call messy(se(ce),ta(tn)%C(1:slen),ta(tn)%ip,ta(tn)%DP,ta(tn)%iap,ta(tn)%dap,&
        & ta(tn)%zp,ta(tn)%zap,ptext=ta(tn)%p(1:plen))
    END IF
    flush(unit=se(ce)%munit)
    flush(unit=se(ce)%eunit)
! End of Transformations C assumed size arrays to Fortran assumed shape arrays. 
    
! Transfer contents back to C struct from Fortran messy_ty type.
    E%FPPREC = SE(CE)%FPPREC; E%KDF = SE(CE)%KDF;       E%LINE_LEN = SE(CE)%LINE_LEN
    E%MUNIT =  SE(CE)%MUNIT;  E%EUNIT = SE(CE)%EUNIT;   E%MAXERR = SE(CE)%MAXERR
    E%LSTOP =  SE(CE)%LSTOP;  E%LPRINT = SE(CE)%LPRINT; E%ERRCNT = SE(CE)%ERRCNT
    E%DBLEV =  SE(CE)%DBLEV
#if numt_ > 1
!$OMP END CRITICAL
#endif
!----------------------------------------------------------------------
  END SUBROUTINE CALLMESSY_
  subroutine GET_CMESSY_DEFAULTS_ (S) bind(C)
    type(CMESSY_TY_), intent(inout) :: S
    type(CMESSY_TY_) :: T

! Here space for S is passed (undefined) from C -- i.e. it is not initialized.
! But T (local) is initialized by Fortran.
! So this local object is copied into S.

! The C code needs default settings for its copy of struct S.
! So S (with no defaults set) is written by the Fortran copy, T,
! which has defaults. This makes S acquire
! the Fortran defaults.
    character(len=9) :: undef='Undefined'
    integer j
    if(.not. (ALLOCATED(TA).and.ALLOCATED(SE))) then
      U%ENAME=QGET_CMESSY_DEFAULTS_
      call messy(U,"$E36Before any call you first must execute 'allocate_$P&
        &messy_interface(maxcstructs, maxcthreads).'",idat=[05], ptext=T_)
    end if
    S=T
    DO j=1,9
      S % ename(J)=undef(j:j)
    END DO
    S % ename(10)=C_NULL_CHAR ! Just the C string "Undefined" for the default name.
    cindex=cindex+1
    if(cindex > maxstructs) THEN
      U%ename=QGET_CMESSY_DEFAULTS_
      call messy(U, "$E36The number $I of different copies of structure&
        & $Pmessy_ty available to use is restricted to $I.  Need to increase&
        & value maxcstructs to the number used in C applications.$N",&
        & idat=[06,cindex,maxstructs], ptext=T_)
    end if
    S % cstruct = cindex
  end subroutine GET_CMESSY_DEFAULTS_

  subroutine OPEN_CMESSY_FILES_(S, task, file_name) bind(C)
! Open or close a file with name of the string in FILE_NAME.
! Place a unit number on or two components in S if the file is opened.
! (This is done with unit_number=NEWUNIT in the open statement.)
! Named files are opened with STATUS='REPLACE'.

! Task: Default units are OUTPUT_UNIT for all files.  Possible changes are
!     1: Use file_name for plain message files.
!     2: Use file_name for error message files.
!     3: Use file_name for both types.
!     0: Use ERROR_UNIT for error message files.
! FILE_NAME: C char string with the name for the file.
    type(CMESSY_TY_), intent(inout) :: S
    integer(c_int), value :: task
    character(c_char), intent(in) :: file_name(*)

    character(max_flen) :: Fortran_Name
    integer:: j, stat

    U%ename=QOPEN_CMESSY_FILES_
    if (task /= 0) then
      Fortran_name=' '
      do j = 1, max_flen
        if (file_name(j) == c_null_char) exit
        fortran_name(j:j) = file_name(j)
      end do
    end if
    select case(task)
    case(0) ! Set ERROR_UNIT for output files
      s%eunit = ERROR_UNIT
    case(1) ! Message files
      open(file=trim(fortran_name),status='REPLACE',&
        & NEWUNIT=S%munit,iostat=stat)
    case(2) ! Error files
      open(file=trim(Fortran_Name),status='REPLACE',&
        & NEWUNIT=S%eunit,iostat=stat)
    case(3) ! Both Message and Error files
      open(file=trim(Fortran_Name),status='REPLACE',&
        & NEWUNIT=S%munit,iostat=stat)
      S%eunit=S%munit
    case default
      call messy(U, '$E99In a call to open_$Pmessy_files, the TASK argument&
        & is $I, but it must be 0, 1, 2 or 3.  Use 1 for messages,&
        & 2 for errors and 3 for both writing or closing the same file.',&
        & idat=[07,task], ptext=T_)
    end select
    if(stat /= 0) then
      call messy(U, "$E99In a call to open_$Pmessy_files, the OPEN&
        & statement with task=$I, stat=$I, failed.$C",&
        & idat=[08,task,stat], ptext=T_)
      call messy(U, "Possible cause is that the name '$P' is a non-existent&
        & or invalid file name or subdirectory.", ptext=trim(fortran_name))
    end if
    return
  end subroutine OPEN_CMESSY_FILES_

  subroutine  CLOSE_CMESSY_FILES_(s) bind(C)
    type(CMESSY_TY_), intent(inout) :: S
    integer :: stat
    U%ename= QLOSE_CMESSY_FILES_
    if (s%eunit /= s%munit) then
      if (s%eunit /= ERROR_UNIT) close(s%eunit, iostat=stat)
      if (stat /= 0) call messy(u, "$E36Close of file failed for&
        & error unit. Stat=$I.",idat=[09,stat])
    end if

    if (s%munit /= OUTPUT_UNIT) close(s%munit, iostat=stat)
    if (stat /= 0) call messy(u, "$E36Close of file failed for&
      & message unit. Stat=$I.",idat=[10,stat])
    return
  end subroutine CLOSE_CMESSY_FILES_

  subroutine ALLOCATE_CMESSY_INTERFACE_(maxcstructs, maxcthreads) bind(C)
    integer, intent(in), value :: maxcstructs, maxcthreads
    integer stat

! The value of maxcstructs must be no smaller than the number
! of different C structures cmessy_ty used in a C program.

! The value of maxcthreads should be no smaller than the
! number of parallel OpenMP threads used in a C program.
    U%ename=QALLOCATE_CMESSY_INTERFACE_
    ALLOCATE(TA(0:maxcthreads-1),STAT=stat)
    if(stat /= 0) THEN
      call messy(u,"$E36Allocation of array TA(:) for threaded arguments failed.&
        &  Possible cause: array already allocated.  STAT=$I",idat=[11,stat])
    end if
    ALLOCATE(SE(maxcstructs),stat=stat)
    if(stat /= 0) THEN
      call messy(u,"$E36Allocation of array SE(:) for storing interface copies&
        & of $Pmessy_ty failed. Possible cause: array already allocated.&
        &  STAT=$I", idat=[12,stat], ptext=T_)
    end if
    maxstructs=maxcstructs
    maxthreads=maxcthreads
  end subroutine ALLOCATE_CMESSY_INTERFACE_

  subroutine DEALLOCATE_CMESSY_INTERFACE_() bind(C)
    integer stat
    U%ename=QALLOCATE_CMESSY_INTERFACE_
    DEALLOCATE(TA,STAT=stat)
    if(stat /= 0) THEN
      call messy(u,"$E36Deallocation of array TA(:) for threaded arguments&
        & failed. Possible cause: array was not allocated.  STAT=$I",&
        & idat=[13,stat])
    end if
    DEALLOCATE(SE,STAT=stat)
    if(stat /= 0) THEN
      call messy(u,"$E36Deallocation of array SE(:) for storing interface&
        & copies of $Pmessy_ty failed. Possible cause: array was not&
        & allocated.  STAT=$I", idat=[14,stat], ptext=T_)
    end if
    maxstructs=0
    maxthreads=1
  end subroutine DEALLOCATE_CMESSY_INTERFACE_
