! Fortran_Virtual_Memory
! Integer version

! Authors: J.K. Reid and J.A. Scott
!          Computational Science and Engineering Department
!          Rutherford Appleton Laboratory
!          Chilton, Didcot
!          Oxfordshire OX10 0QX
!          U.K.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Summary
! =======
! This package provides read/write facilities for 
! one or more direct-access files through a single in-core buffer, so
! that actual input-output operations are often avoided. The
! buffer is divided into fixed-length pages and all 
! input-output is performed by transferring a single page to or from a
! single record of a file (the length of a record is equal to the
! length of a page). 

! Each set of data is addressed as a virtual array, that is, as if it
! were a very large array. The lower bound of the virtual array is 1.
! Each element of the virtual array has initial value zero. 
! Any contiguous section of the virtual array
! may be read or written, without regard to page boundaries. 

! The virtual array is permitted to be too large to be accommodated 
! in a single file, in which case FVM opens secondary files 
! with names that it constructs from the name of the primary file by 
! appending 1, 2, ... . We refer to the set of files as a superfile.
! Each superfile is identified by the name 
! of its primary file or the index that it is given when it is opened. 
! To allow the secondary files to reside on different devices, the user is
! required to supply an array of path names; the full name of a file is the
! concatenation of a path name with the file name. 

! To facilitate finite-element assembly and the multifrontal method, 
! there is an option to add data from the virtual array to a given array
! under the control of a map. 

! Full details of the calling sequences and argument lists are given
! in the accompanying user documentation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Fortran_Virtual_Memory_integer
  implicit none
  private

  public FVM_data, FVM_data_private
  public FVM_initialize, FVM_open, FVM_close, FVM_read, FVM_write, &
    FVM_end


  !!! Parameters !!!
  integer, parameter  :: default_lpage = 2**12 ! Default page length.

  integer, parameter  :: default_maxfiles = 10 ! Maximum number of open files.

  integer, parameter  :: default_npage = 1600 ! Default number of pages
!          in the buffer.

  integer, parameter  :: long = selected_int_kind(18) ! Long integer.

  integer, parameter  :: default_file_size = 2**21 ! Default
!          target length of each file.

  integer, parameter  :: nsup = 2 ! Initial size of the array
!          data%private%filename.

  integer, parameter  :: wp = kind(0) ! Defines data type (single/double).

  integer, parameter  :: ihash = 3 ! Constant used in the hashing function.

  integer, parameter  :: maxpath = 400 ! Max. length of path name.

  integer, parameter  :: maxname = 400 ! Max. length of superfile name.

  integer, parameter  :: maxlen = maxpath+maxname+10 ! Max. length of file name
!          There has to be a limit on the lengths of path and file names because
!          Fortran 95 requires the file name in an OPEN statement to be a scalar
!          character variable.

  integer(wp), parameter :: zero = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type FVM_data_private

! Note: Age is measured since last reference.
    integer(wp), allocatable :: buffer(:,:) ! of shape (lpage,npage)

    logical,  allocatable :: differs(:) ! of size npage.
!             differs(I)=.true. if the file version of page I
!             of the buffer is different from the buffer version.

    character(maxname), allocatable :: filename(:)
!             Set to hold a copy of the superfile names given to FVM_open.

    integer(long), allocatable :: first(:) ! of size npage. first(K) is the
!            first page with hash code K or zero if there are none.

    integer :: free !  Start of linked list of free file indices.

    integer(long), allocatable :: left(:),right(:) ! of size maxfiles.
!            A sequence of discarded entries on a superfile is recorded as
!            left(superfile):right(superfile).

    integer(long), allocatable :: highest(:) ! of size maxfiles.
!            highest(I) holds:
!            -1 after call of FVM_initialize,
!            -2 after call of FVM_close for file I, or
!            >0 highest page number read or written on file I.

    integer, allocatable :: index(:) ! of size npage. index(I) contains the
!            index of the (super)file associated with page I of the buffer.

    integer :: iolength ! iolength of a page

    integer :: maxfiles = 0 ! Maximum number of open files.

    integer, allocatable :: name(:) ! of size maxfiles. For a primary file,
!            name(I) holds the position of the superfile name in filename.

    integer(long), allocatable :: next(:) ! of size npage.
!            next(I) is the next page with same hash code
!            or zero for last in list. Otherwise undefined.

    integer :: nfiles !  Number of file indices in use.

    integer(long) :: nrec ! Number of records in each file.

    integer, allocatable :: nextfile(:) ! of size maxfiles. Used for the
!            linked list of free file indices and linked lists of indices of
!            files in a superfile.  nextfile(I) is the next file index in
!            its list or zero if there are none.

    integer(long), allocatable :: older(:) ! of size npage. older(I) is the
!            next older page to page I, or the youngest page if I is oldest.

    integer(long), allocatable :: page(:) ! of size npage. page(I) is the
!            file address (page number) of page I in the buffer.

    integer(long), allocatable :: page_list(:) ! of size npage.
!            list of pages for this read or write
!            that were in the buffer at the time of call.

    character(maxpath), allocatable :: path(:)
!            The paths given to FVM_initialize

    integer(long), allocatable :: prev(:) ! of size npage. prev(I) is the
!            previous page with same hash code or -(hash code)
!            for first in list. Otherwise undefined.

    integer, allocatable :: unit(:) ! of size maxfiles. unit(I) holds the
!            unit for file I or zero if not in use.

    integer(long), allocatable :: younger(:) ! of size npage.
!            younger(I) is the next younger page to page I, or
!            the oldest page if I is the youngest.

    integer :: youngest ! most recently referenced page in the buffer

  end type FVM_data_private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type FVM_data

    integer :: entry = 0 ! Index of the entry to the package.

    integer :: iostat    ! Fortran iostat parameter.
    integer :: lpage     ! length of each page of the buffer, that is, the
!            number of scalar variables in each page. The default is 1024.

    integer(long) :: ncall_read  ! number of calls to FVM_read.

    integer(long) :: ncall_write ! number of calls to FVM_write.

    integer(long) :: nio_read ! number of records read by
!            FVM_read and FVM_write.

    integer(long) :: nio_write ! number of records written by
!            FVM_read and FVM_write.

    integer(long) :: npage ! number of pages in the in-core buffer.
!            The default value is 20. It is a long integer, but
!            we do not anticipate it being very large.

    integer(long) :: file_size ! target length of each file.
!             The default value is default_file_size.

    integer(long) :: nwd_read ! number of scalars read by FVM_read.

    integer(long) :: nwd_write ! number of scalars written by FVM_write.

    type (FVM_data_private) :: private
! In Fortran 2003, this should be
!   type (FVM_data_private), private :: private

    integer :: stat ! Fortran stat parameter.

  end type FVM_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface FVM_initialize
    module procedure FVM_initialize_integer
  end interface
  interface FVM_open
    module procedure FVM_open_integer
  end interface
  interface FVM_close
    module procedure FVM_close_integer
  end interface
  interface FVM_read
    module procedure FVM_read_integer
  end interface
  interface FVM_write
    module procedure FVM_write_integer
  end interface
  interface FVM_end
    module procedure FVM_end_integer
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_initialize_integer(iflag,data,path,file_size,lpage,npage,lp)

! This subroutine must be called once to initialize
! a structure of derived type FVM_data
! and to allocate and initialize its array components.

    integer, intent (out) :: iflag ! Value 0 on successful return.
!            The only possible negative values are:
!            -1 Allocation error. The STAT parameter returned in data%stat.
!            -8 Deallocation error. The STAT parameter returned in data%stat.
!           -16 Character length of path too great.

    type (FVM_data) , intent (inout) :: data

    character(*), optional, intent (in) :: path(:) ! Paths for the files.
!            If absent, (\ '' \) is used.

    integer(long), optional, intent (in) :: file_size ! If present and positive,
!            the target length of a file. Default value
!            used if absent or present and not positive.

    integer, optional, intent (in) :: lpage ! If present and positive, the
!            size of each page in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: npage ! If present and positive, the
!            number of pages in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: lp ! If present and not negative, the
!            unit number for diagnostic messages. 6 is
!            used if absent or present and negative.

! local variables
    integer :: i !  do loop variable
    integer :: size_path ! size of array path

!!!!!!!!!!!!!!!!!!!!!!!
! Initialise
    iflag = 0
    data%entry       = 1
    data%ncall_read  = 0
    data%ncall_write = 0
    data%nio_read    = 0
    data%nio_write   = 0
    data%nwd_read    = 0
    data%nwd_write   = 0

! If napge/lpage supplied, use instead of the default values
    data%npage = default_npage
    if (present(npage)) then
      if (npage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%npage = npage
    end if

    data%lpage = default_lpage
    if (present(lpage)) then
      if (lpage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%lpage = lpage
    end if

! Find record length
    if (allocated(data%private%buffer)) then
       deallocate (data%private%buffer,stat=data%stat)
    end if
    allocate (data%private%buffer(data%lpage,1),stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -1; call print_iflag(data,iflag,lp); return
    end if
    inquire (iolength=data%private%iolength) data%private%buffer(:,1)
    deallocate (data%private%buffer,stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -8; call print_iflag(data,iflag,lp); return
    end if

    size_path = 1
    if (present(path)) then
       size_path = size(path)
       if (len(path)>maxpath) then
         iflag = -16; call print_iflag(data,iflag,lp); return
       end if
    end if

    data%file_size = default_file_size
    if (present(file_size)) then
      if (file_size < data%lpage) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%file_size = file_size
    end if
    data%private%nrec = data%file_size/data%lpage
    data%file_size = data%private%nrec*data%lpage

! Allocate private arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,   data%private%nextfile, &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%path,      data%private%name,     &
                  data%private%left,      data%private%right,    &
                  data%private%filename, stat=data%stat)
       if (data%stat /= 0) then
         iflag = -8; call print_iflag(data,iflag,lp); return
       end if
    end if
    data%private%maxfiles = default_maxfiles
    allocate (data%private%highest(data%private%maxfiles),  &
              data%private%nextfile(data%private%maxfiles), &
              data%private%index(1:data%npage),             &
              data%private%older(1:data%npage),             &
              data%private%younger(1:data%npage),           &
              data%private%differs(1:data%npage),           &
              data%private%prev(1:data%npage),              &
              data%private%next(1:data%npage),              &
              data%private%first(1:data%npage),             &
              data%private%page(1:data%npage),              &
              data%private%page_list(1:data%npage),         &
              data%private%unit(data%private%maxfiles),     &
              data%private%path(size_path),                 &
              data%private%name(data%private%maxfiles),     &
              data%private%left(data%private%maxfiles),     &
              data%private%right(data%private%maxfiles),    &
              data%private%filename(nsup),                  &
              data%private%buffer(data%lpage,data%npage), stat=data%stat)
    if (data%stat /= 0) then
      data%private%maxfiles = 0
      iflag = -1; call print_iflag(data,iflag,lp); return
    end if

! Initialise arrays
    data%private%free = 0
    data%private%nfiles = 0
    data%private%filename(:) = ''
    do i = 1, data%npage
      data%private%index(i) = -1
      data%private%older(i) = i + 1
      data%private%younger(i) = i - 1
      data%private%page(i) = 0
      data%private%first(i) = 0
      data%private%next(i) = 0
      data%private%prev(i) = 0
      data%private%differs(i) = .false.
    end do
    data%private%youngest = 1
    data%private%younger(1) = data%npage
    data%private%older(data%npage) = 1
    if (present(path)) then
       data%private%path(:) = path(:)
    else
       data%private%path(:) = (/ '' /)
    end if

  end subroutine FVM_initialize_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_open_integer(filename,ifile,iflag,data,lenw,lp)

! This subroutine must be called for each superfile that is to be
! accessed through Fortran_Virtual_Memory. It gives the superfile an index and opens
! its files.

    character(len=*), intent (in) :: filename ! Name of the superfile.

    integer, intent (out) :: ifile ! Index of the superfile.

    integer, intent (out) ::  iflag ! Value 0 on successful return.
!             A negative value is associated with an error message
!               on unit lp. Possible negative values are:
!            -1  Allocation error. The STAT parameter returned in data%stat.
!            -5  Error in inquire statement.
!            -7  Error in open statement.
!            -8  Deallocation error. The STAT parameterreturned in data%stat.
!            -11 lenw > 0, but not enough files exist.
!            -12 File filename already exists, but lenw is not present
!                or lenw <=  0.
!            -13 filename is too long.
!            -17 open unsuccessful for all elements of path.

    type (FVM_data), intent (inout) :: data

    integer(long), optional, intent (in) :: lenw ! length in pages of the
!            part of the file that has been written and is not
!            regarded as having been overwritten by zeros.
!            Pages beyond this are regarded as containing zeros.
!            If not present or lenw <=  0, a new file is opened.

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.
!            If lp < 0, messages are suppressed.

! Local variables
    character(10) :: ci ! Filename extension
    integer :: i ! do loop variable
    integer :: k ! number of secondary files
    integer :: l ! file index
    integer(long) :: lenw_copy ! no. of pages in the virtual array
    integer :: nout ! unit for error messages
    integer :: m ! Previous file index in linked list
    character(maxname), allocatable :: temp(:)

    iflag = 0
    nout = 6
    if (present(lp)) nout = lp
    data%entry = 2

    if (len(filename) > maxname) then
       iflag = -13; go to 100
    end if

! Find the length of the virtual array
    lenw_copy = 0
    if (present(lenw)) lenw_copy = max(0_long,lenw)

! Find number of secondary files
    k = (lenw_copy-1)/data%private%nrec

! Find suitable indices and units for all the files and open them
    m = 0
    ci=''
    do i = 0,k
      if (lenw_copy > 0) then
        call open_old (data,l,trim(filename)//adjustl(ci),iflag)
      else
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
      end if
      if (iflag /= 0) go to 100
      if (i==0) ifile = l
      data%private%highest(l) = max(0_long,lenw_copy)
      data%private%left(l) = 1
      data%private%right(l) = 0
      if (m > 0) data%private%nextfile(m) = l
      m = l
      data%private%highest(l) = lenw_copy
      write(ci,'(i5)') i + 1
    end do

! Look for a place and store filename there
    k = size(data%private%filename)
    do i = 1, k
      if (data%private%filename(i)== '') then
          data%private%filename(i) = filename
          data%private%name(ifile) = i
          return
      end if
    end do
! Increase size of data%private%filename
    allocate(temp(k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    temp(:) = data%private%filename(:)
    deallocate(data%private%filename,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
    allocate(data%private%filename(2*k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    data%private%filename(1:k) = temp(:)
    data%private%filename(k+1:) = ''
    deallocate(temp,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
! Store filename
    data%private%filename(k+1) = filename
    data%private%name(ifile) = k+1
    return

100 call print_iflag(data,iflag,lp)

  end subroutine FVM_open_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_read_integer &
       (ifile,loc,n,read_array,iflag,data,lp,map,discard)

! This subroutine performs the read from a superfile.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be read.

    integer(wp), intent (inout) :: read_array(*)
!            Array into which data from the file is read.
!            If map is present, file data is added into it thus:
!               read_array(map) = read_array(map) + file_data(1:n)

    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not open.
!            -5  Error in Fortran INQUIRE. The IOSTAT parameter is returned in
!                data%iostat.
!            -6  Error in Fortran READ. The IOSTAT parameter is returned in
!                data%iostat.
!            -7  Error in Fortran OPEN statement. The IOSTAT parameter
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.

    type (FVM_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

    integer, optional, intent (in) :: map(n) ! map array.

    logical, optional, intent (in) :: discard ! Whether data is to be discarded.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos  ! position within first page of first value
!           required.
    integer(long) :: left,right ! Range of discarded entries
    integer(long) :: ip  ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer
!           pages that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer       :: jfile ! file associated with a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!           required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer       :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: my_discard ! Whether the data is to be discarded.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (n <= 0) return
    iflag = 0
    data%entry = 3

! Check for errors in incoming data
    if (ifile < 1) then
       iflag = -9; go to 200
    else if (ifile > data%private%nfiles) then
       iflag = -4; go to 200
    else if (data%private%highest(ifile) < 0) then
       iflag = -4; go to 200
    else if (loc <= 0) then
       iflag = -3; go to 200
    end if

    my_discard = .false.
    if (present(discard)) then
      if (discard) then
        my_discard = .true.
        left = loc
	  right = loc+n-1
	  if (loc == data%private%right(ifile)+1) then
	    left = data%private%left(ifile)
        else if (loc+n == data%private%left(ifile)) then
	    right = data%private%right(ifile)
	  end if
        data%private%left(ifile) = left
	  data%private%right(ifile) = right
      end if
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_read = data%ncall_read + 1
    data%nwd_read = data%nwd_read + n

    lenl = 0
! Look for required pages that are in the buffer
    look1: do ip = first_page, last_page

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp) == ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
10    call read_from_buffer
      if (my_discard) then
        if ((ip-1)*lpage+1 >= left .and. ip*lpage+1 <= right ) &
          data%private%differs(jp) = .false.
      end if
   end do look1

! look for pages not yet treated and treat each of them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw.le.lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if
! the oldest page is overwritten and becomes the youngest
      jp = data%private%younger(data%private%youngest)
      data%private%youngest = jp

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile>=0) then
! write page from buffer to record page(jp) in the file
! if it has been altered (ie if differs(,jp) is set to 1).
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 200

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip > data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 200
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m > 0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .false.

! perform actual transfer to read_array from buf
      call read_from_buffer

    end do look2
!    call buffer
    return

200 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_from_buffer

      integer :: l1,l2 ! Temporaries
! Read page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip == last_page) l2 = last_pos
      if (.not. present(map)) then
        call FVM_icopy(l2-l1+1,data%private%buffer(l1,jp),read_array(m))
      else
        call FVM_mcopy(l2-l1+1,data%private%buffer(l1,jp),read_array,map(m))
      end if
    end subroutine read_from_buffer

  end subroutine FVM_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_write_integer &
             (ifile,loc,n,write_array,iflag,data,lp,inactive)

! This subroutine performs the writing to a superfile.

    integer, intent (in) :: ifile ! Index of the (super)file.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be written.

    integer(wp), intent (in) :: write_array(*) ! array from which data
!            is written to the file.
    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not FVM_open.
!            -5  Error in INQUIRE.  The IOSTAT value
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (FVM_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.

    integer(long), optional, intent (in) :: inactive ! If present, it is
!            assumed that the data transferred are unlikely to be needed before
!            other data in the buffer. If inactive<loc, it is assumed that this
!            is true also for the data in the interval inactive:loc. If
!            inactive loc+n, it is assumed that this  is true also for the data
!            in the interval loc+n:inactive.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos !  position within first page of first value
!            required.
    integer(long) :: ip ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer pages
!            that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer :: jfile ! file index of a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!            required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: inactive_page ! Whether this page is fully inactive.
    logical       :: inactive_first ! Whether first page is fully inactive.
    logical       :: inactive_middle ! Whether middle pages are fully inactive.
    logical       :: inactive_last ! Whether last page is fully inactive.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (n <= 0) return
    iflag = 0
    data%entry = 4

! Check for errors in incoming data
    if (ifile < 1) then
      iflag = -9; go to 300
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 300
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 300
    else if (loc <= 0) then
      iflag = -3; go to 300
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

    if (present(inactive)) then
       if (first_page==last_page) then
          inactive_first = .false.
       else
          inactive_first = inactive<=1+(first_page-1)*lpage
          inactive_middle = .true.
          inactive_last = inactive>=last_page*lpage
       end if
    else
       inactive_middle = .false.
       inactive_first = .false.
       inactive_last = .false.
    end if

! If there is any overlap with the range of discarded entries,
! reset the range to be null
    if (loc>data%private%right(ifile)) then
    else if (loc+n-1<data%private%left(ifile)) then
    else
       data%private%left(ifile) = 1
       data%private%right(ifile) = 0
    end if

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_write = data%ncall_write + 1
    data%nwd_write = data%nwd_write + n

    lenl = 0
! Look for required pages in the buffer
    look1: do ip = first_page, last_page
      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if
!     if (inactive_page) write(10,'(i3,i9,3i5)') ifile,loc,n,ip,data%lpage

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp)==ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          call write_to_buffer
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
      call write_to_buffer
! If the page is fully inactive, make it the oldest
 10   if (inactive_page) then
         data%private%youngest = data%private%older(jp)
      end if
    end do look1

! look for pages not yet treated and treat them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw <= lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if

      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if

! the oldest page is overwritten
      jp = data%private%younger(data%private%youngest)
! if transfer is active or the page is partially active, page becomes
! the youngest
     if (.not.inactive_page) then
        data%private%youngest = jp
     end if

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile >= 0) then
! write page from buffer to record page(jp) in the file
! if it has been altered.
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 300

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip>data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 300
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m>0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .true.

! perform actual transfer to write_array from buf
      call write_to_buffer

    end do look2
!    call buffer
    return

300 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_to_buffer

      integer :: l1,l2 ! Temporaries

! Write page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip==last_page) l2 = last_pos
      call FVM_icopy(l2-l1+1,write_array(m),data%private%buffer(l1,jp))
      data%private%differs(jp) = .true.

    end subroutine write_to_buffer
  end subroutine FVM_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_close_integer(ifile,lenw,num_file,iflag,data,lp,lkeep)

! This subroutine should be called when a superfile is finished with.
! Its data in the buffer are unloaded.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (out) :: lenw ! length in pages of the part of
!            the file that has been written.
    integer, intent (out) :: num_file ! total number of additional files
!            that have been opened. They are named
!            filenamej, j = 1, ..., num_file.
    integer, intent (out) :: iflag ! iflag = 0 for a successful return.
!            A negative value is associated with an error
!            message on unit lp. Possible negative values:
!            -4  Attempt to access a file that is not FVM_open.
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -9  ifile out of range.
!            -14 Error in CLOSE statement. The IOSTAT value
!                is returned in data%iostat.
!            -15 Error in WRITE.  The IOSTAT value
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (FVM_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0 , messages are suppressed.

    logical, optional, intent (in) :: lkeep ! If present and set to .FALSE.,
!            the files are deleted on being closed.
!            Otherwise, the files are kept.

! Local variables
    integer(long) :: i !  do loop variable
    integer :: j       !  do loop variable
    integer :: jfile   !  index of a file
    integer(long) :: l ! page number
    integer :: lpage   ! length of each page (= data%lpage)
    integer(long) :: m ! page number
    integer(long) :: oldest ! oldest page
    integer(long) :: oldest2 ! second oldest page
    character (maxlen) name  ! File name
    character(len=6) :: status ! File status on closing

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 5
! Check for errors
    if (ifile < 1) then
      iflag = -9; go to 400
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 400
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 400
    end if

    status = 'keep'
    if (present(lkeep)) then
      if (.not. lkeep) status = 'delete'
    end if

! search the pages of the buffer dealing with those associated
! with ifile.
    search:do i = 1, data%npage
! data%private%index(i) is index of the (super)file associated
! with page i of buffer
      if (ifile /= data%private%index(i)) cycle
      data%private%index(i) = -1
      if (data%private%differs(i)) then
! File version of page i is not the same as the buffer version so
! write page i to file
        data%private%differs(i) = .false.
        lpage = data%lpage
        if (status == 'keep') then
          call page_write(data, &
          data%private%buffer(:,i),ifile,data%private%page(i),iflag)
          if (iflag /= 0) then
            call print_iflag(data,iflag,lp)
            return
          end if
        end if
      end if

! remove page from hash list
      l = data%private%prev(i) ! Previous
      m = data%private%next(i) ! Next
      if (m>0) data%private%prev(m) = l
      if (l>0) then
        data%private%next(l) = m
      else if (l < 0) then
        l = -l  ! i is first in the list with code l	
        data%private%first(l) = m
      end if

! youngest is the most recently referenced page in the buffer.
! younger(youngest) is the oldest page in the buffer.
      oldest = data%private%younger(data%private%youngest)
! if page i is not already the oldest, make it so.
      if (oldest==i) cycle
! remove page from its present position
      l = data%private%younger(i)
      m = data%private%older(i)
! l is next youngest to page I and m is the next oldest
      data%private%older(l) = m
      data%private%younger(m) = l

! insert page i as the oldest page
! older(oldest) is the youngest page
      data%private%youngest = data%private%older(oldest)

! Page i becomes the oldest page. The page that was the oldest
! becomes the second (ie next) oldest page
      oldest2 = oldest
      oldest = i
      data%private%older(oldest2) = oldest
      data%private%younger(oldest) = oldest2
      data%private%older(oldest) = data%private%youngest
      data%private%younger(data%private%youngest) = oldest
    end do search

    lenw = data%private%highest(ifile)
    if (lenw <= 0 ) status = 'delete'
    num_file = (lenw-1)/data%private%nrec

! remove file from list of declared files.
! close files
    inquire(unit=data%private%unit(ifile),name=name,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -5; go to 400
    end if
    close (unit=data%private%unit(ifile),status=status,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -14; go to 400
    end if
    data%private%unit(ifile) = 0
    data%private%highest(ifile) = -1
    jfile = ifile
    if (num_file > 0) then
      do j = 1, num_file
        m = data%private%nextfile(jfile)
        inquire(unit=data%private%unit(m),name=name,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -5; go to 400
        end if
        close (unit=data%private%unit(m),status=status,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -14; go to 400
        end if
        data%private%unit(m) = 0
        data%private%highest(m) = -1
        jfile = m
      end do
    end if
! Put the freed indices at the front of the list list of free indices
    data%private%nextfile(jfile) = data%private%free
    data%private%free = ifile

! Clear the name for reuse
    data%private%filename(data%private%name(ifile))=''

    return

400 call print_iflag(data,iflag,lp)

  end subroutine FVM_close_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_end_integer(iflag,data,lp)

! This subroutine deallocates the private array components of data.

    integer, intent (out) :: iflag
!             0  successful return.
!            -8  Error in Fortran DEALLOCATE statement. The STAT
!                value is returned in data%stat.
!            -10 FVM_close has not been called for a superfile.

    type (FVM_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0, messages are suppressed.

! Local variables
    integer :: i

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 6

    do i = 1, data%private%nfiles
      if (data%private%highest(i) >= 0) then
        iflag = -10; call print_iflag(data,iflag,lp)
        return
      end if
    end do

! Deallocate arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,  data%private%nextfile,  &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%buffer,    data%private%path,     &
                  data%private%left,      data%private%right,    &
                  data%private%name,      data%private%filename, &
                  stat=data%stat)
      data%private%maxfiles = 0
      if (data%stat /= 0) then
        iflag = -8; call print_iflag(data,iflag,lp)
      end if
    end if

  end subroutine FVM_end_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_write(data,buffer,ifile,locp,iflag)

! This subroutine transfers a page of the buffer to a file. If the
! page is within the existing virtual array then a simple write
! is performed. To write a page beyond the present end may require
! writing pages full of zeros in front of its position.

    type (FVM_data), intent (inout) :: data

    integer(wp), intent (inout) :: buffer(data%lpage) ! array from which
!             the transfer is made.
    integer, intent (in) :: ifile ! the index of the primary file.

    integer(long), intent (in) :: locp ! page number in the file.

    integer, intent (out) :: iflag ! error flag. Possible negative values:
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -7  Error in OPEN. IOSTAT valueis returned in data%iostat.
!            -15 Error in WRITE. IOSTAT value is returned in data%iostat.

! Local variables
    integer(long) :: hpage ! highest page read or written for ifile
    integer(long) :: l ! do loop variable

    hpage = data%private%highest(ifile)
    if (locp <= hpage+1) then
! Write straight to file
      call actual_page_write(locp)
      if (iflag /= 0) return
    else
! Fill the gap with zero page
! first write out buffer to page hpage+1
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
      buffer(:) = zero
! write pages of zeros
      do l = hpage+2, locp-1
        call actual_page_write(l)
        if (iflag /= 0) return
      end do
! read buffer back in from page hpage+1
      call page_read(data,buffer,ifile,hpage+1,iflag)
      if (iflag /= 0) return
! Now write buffer to the required location
      call actual_page_write(locp)
      if (iflag /= 0) return
! Finally write zeros to hpage+1
      buffer(:) = zero
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
    end if

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine actual_page_write(locp)

! Perform an actual i/o write of a page to a superfile. If necessary, open
! intermediate secondary files.

    integer(long), intent (in) :: locp ! page index in the superfile.

! Local variables
    character(10) :: ci ! Filename extension.
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or >0 for k-th secondary file.
    integer l ! index of the file that is written.
    integer m ! Previous file in linked list.
    character(len=maxlen) :: filename ! name of the primary file.

    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec
    l = ifile
    do i = 1, k
      m = l
      l = data%private%nextfile(l)
      if (l == 0) then
! open new file (filename1, filename2, etc)
        filename = data%private%filename(data%private%name(ifile))
        write(ci,'(i5)')i
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
        if (iflag /= 0)  return

! Add to linked list for ifile
        if (m > 0) data%private%nextfile(m) = l
        m = l
      end if
    end do

    data%nio_write = data%nio_write + 1
    write (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
!    write(*,*) 'write record',new_locp,data%iostat
    if (data%iostat /= 0) then
      iflag = -15; return
    end if
    data%private%highest(ifile) = max(data%private%highest(ifile),locp)

  end subroutine actual_page_write

  end subroutine page_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_read(data,buffer,ifile,locp,iflag)

! This subroutine performs an actual i/o read of a page from a superfile

    type (FVM_data), intent (inout) :: data

    integer(wp), intent (out) :: buffer(data%lpage) ! array to receive page.

    integer, intent (in) :: ifile ! index of the superfile.

    integer(long), intent (in) :: locp ! page index in the superfile.

    integer, intent (out) :: iflag ! Error flag. Possible negative value:
!            -6 error in read statement

! Local variables
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or > 0 for k-th secondary file.
    integer l ! index of the file that is read.

!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec

! Find index of the kth linked file
    l = ifile
    do i = 1, k
      l = data%private%nextfile(l)
    end do

    data%nio_read = data%nio_read + 1
    read (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
    if (data%iostat /= 0) iflag = -6

  end subroutine page_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_new (data,ifile,filename,iflag)

! Open a new file

    type (FVM_data) :: data

    integer, intent(out) :: ifile ! File index.

    character (*) :: filename!    ! File name.

    integer, intent (inout) :: iflag ! Error flag. Possible nonzero value:
!            -7 Error in OPEN.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: j ! Do index.
    integer :: unit ! Unit.
    logical :: open ! open variable for inquire.

! Check file does not already exist
    do i = 1, size(data%private%path)
       inquire (file=trim(data%private%path(i))//filename,exist=ex,&
          iostat=data%iostat)
       if (data%iostat /= 0) then
          iflag = -5; return
       else if (ex) then
          iflag = -12; data%iostat=i; return
       end if
    end do

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

! Open the file on the first available path
    outer: do i = 1, size(data%private%path)
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       if (data%iostat /= 0) then
         close(unit,iostat=data%iostat,status='delete')
         cycle
       end if

       if (size(data%private%path) == 1) return
! Check that there is room for the file
       do j = 1, data%private%nrec
         write(unit,rec=j,iostat=data%iostat)data%private%buffer(:,1)
         if (data%iostat /= 0) then
           close(unit,iostat=data%iostat,status='delete')
           cycle outer
         end if
       end do
       close(unit,iostat=data%iostat)
       if (data%iostat /= 0) then
         inquire(unit,opened=open)
         if (.not.open) &
         open (unit,file=trim(data%private%path(i))//filename, &
             access='direct',iostat=data%iostat,recl=data%private%iolength)
         close(unit,iostat=data%iostat,status='delete')
         cycle outer
       end if
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       return
    end do outer
    iflag = -17

  end subroutine open_new

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_old (data,ifile,filename,iflag)

! Open an old file.

    type (FVM_data) :: data

    integer, intent(out) :: ifile ! File index.
    character (*) :: filename     ! File name.

    integer, intent (inout) :: iflag ! Possible nonzero values:
!            -5  Error in inquire statement
!            -7  Error in OPEN.  The IOSTAT value is returned in data%iostat.
!            -11 File does not exist

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: unit ! Unit.

! Find the file.
    do i = 1, size(data%private%path)
      inquire (file=trim(data%private%path(i))//filename,exist=ex,&
             iostat=data%iostat)
      if (data%stat /= 0) then
        iflag = -5; return
      end if
      if (ex) exit
    end do
    if (.not. ex) then
! File not found
      iflag = -11; return
    end if

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

    open (unit,file=trim(data%private%path(i))//filename, access='direct', &
          iostat=data%iostat, recl=data%private%iolength)
    if (data%iostat /= 0) then
      iflag = -7; return
    end if

  end subroutine open_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_unit(data,unit,ifile,iflag)

! Find a free unit and a free file index.

    type (FVM_data) :: data
    integer, intent (out) :: unit ! Free unit.

    integer, intent (out) :: ifile ! Free file.

    integer, intent (inout) :: iflag ! Possible nonzero value:
!           -5 Error in INQUIRE.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: open

! Find free unit
    do unit = 12, huge(0)
      inquire (unit=unit,iostat=data%iostat,opened=open)
      if (data%iostat /= 0) then
! Error in inquire
        ifile = -1
        iflag = -5; return
      end if
      if (.not. open) exit
    end do

! Find free file index
    ifile = data%private%free
    if (ifile > 0) then
      data%private%free = data%private%nextfile(ifile)
    else
      ifile =  data%private%nfiles+1
      if (ifile>data%private%maxfiles) then
! Reallocate files
        data%private%maxfiles = 0
        call reallocate_long(data%private%highest)
        if (iflag < 0) return
        call reallocate_long(data%private%left)
        if (iflag < 0) return
        call reallocate_long(data%private%right)
        if (iflag < 0) return
        call reallocate(data%private%nextfile)
        if (iflag < 0) return
        call reallocate(data%private%name)
        if (iflag < 0) return
        call reallocate(data%private%unit)
        if (iflag < 0) return
        data%private%maxfiles = size(data%private%unit)
      end if
      data%private%nfiles =  ifile
    end if
    data%private%nextfile(ifile) = 0
    data%private%highest(ifile) = -1
    data%private%left(ifile) = 1
    data%private%right(ifile) = 0
    data%private%unit(ifile) = unit

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    subroutine reallocate(array)

      integer, allocatable :: array(:)

      integer, allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)
      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reallocate_long(array)

      integer(long), allocatable :: array(:)

      integer(long), allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)

      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate_long

  end subroutine find_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_iflag(data,iflag,lp)

    type (FVM_data) :: data

    integer, intent (in) :: iflag ! Error flag.

    integer, optional :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

! Local variables
    integer :: i, nout
    character(10), parameter :: names(6) = (/'initialize','open      ',&
                    'read      ','write     ','close     ','end       '/)
!!!!!!!!!!!!!!!!!!!!!!!!!!

    nout = 6
    if (present(lp)) then
      if (lp < 0) return
      nout = lp
! check lp not equal to the unit number of a file that has been declared
      do i = 1, data%private%nfiles
        if (data%private%unit(i) == lp) then
          nout = 6; exit
        end if
      end do
    end if

    write (nout,'(3a,i3)') ' Error return from FVM_', &
        trim(names(data%entry)), '. Error flag = ', iflag

    if (iflag == -1) then
      write (nout,'(a,i3)') ' Allocation error. stat parameter = ', data%stat
    else if (iflag == -2) then
      write (nout,'(a,i3)') ' Violation of restriction on optional argument '
    else if (iflag == -3) then
      write (nout,'(a,i3)') ' loc out of range'
    else if (iflag == -4) then
      write (nout,'(a,i3)') &
        ' the superfile is not open under FVM'
    else if (iflag == -5) then
      write (nout,'(2a,i3)') ' INQUIRE statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -6) then
      write (nout,'(2a,i3)') ' READ statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -7) then
      write (nout,'(2a,i3)') ' OPEN statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -8) then
      write (nout,'(a,i3)') ' Deallocation error. stat parameter = ', &
        data%stat
    else if (iflag == -9) then
      write (nout,'(a)') &
        ' ifile is out of its range'
    else if (iflag == -10) then
      write (nout,'(2a)') ' one or more superfiles are open through FVM'
    else if (iflag == -11) then
      write (nout,'(2a)') ' lenw is positive but one or more of the ', &
             'required files does not exist'
    else if (iflag == -12) then
      write (nout,'(2a)') ' filename already exists in path ', &
             trim(data%private%path(data%iostat))
    else if (iflag == -13) then
      write (nout,'(a,i5)') ' file name is longer than ', maxname
    else if (iflag == -14) then
      write (nout,'(2a,i3)') ' CLOSE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -15) then
      write (nout,'(2a,i3)') ' WRITE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -16) then
      write (nout,'(a,i5)') ' path name is longer than ', maxpath
    else if (iflag == -17) then
      write (nout,'(a,i5)') ' unable to open file of given length'
    end if

  end subroutine print_iflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_icopy(len,buffer,array)
! Copy buffer into array, thus:
!    do l = 1,len
!      array(i) =  buffer(i)
!    end do

    integer, intent(in) :: len
    integer(wp), intent(in) ::  buffer(len)
    integer(wp), intent(inout) :: array(*)

    integer :: i,k
    k = mod(len,7)
    do i = 1, k
      array(i) =  buffer(i)
    end do

    do i = k + 1, len, 7
       array(i)   = buffer(i)
       array(i+1) = buffer(i+1)
       array(i+2) = buffer(i+2)
       array(i+3) = buffer(i+3)
       array(i+4) = buffer(i+4)
       array(i+5) = buffer(i+5)
       array(i+6) = buffer(i+6)
    end do

  end subroutine FVM_icopy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_mcopy(len,buffer,array,map)

! Add buffer into array under the control of map, thus:
!    do l = 1,len
!      i = map(l)
!      array(i) = array(i) + buffer(l)
!    end do

    integer, intent(in) :: len
    integer(wp), intent(in) :: buffer(len)
    integer(wp), intent(inout) :: array(*)
    integer, intent(in) :: map(len)

! Local variables
    integer :: i,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    k = mod(len,7)
    do l = 1, k
      i = map(l); array(i) = array(i) + buffer(l)
    end do

    do l = k + 1, len, 7
       i = map(l);    array(i) = array(i) + buffer(l)
       i = map(l+1);  array(i) = array(i) + buffer(l+1)
       i = map(l+2);  array(i) = array(i) + buffer(l+2)
       i = map(l+3);  array(i) = array(i) + buffer(l+3)
       i = map(l+4);  array(i) = array(i) + buffer(l+4)
       i = map(l+5);  array(i) = array(i) + buffer(l+5)
       i = map(l+6);  array(i) = array(i) + buffer(l+6)
    end do

  end subroutine FVM_mcopy


end module Fortran_Virtual_Memory_integer

