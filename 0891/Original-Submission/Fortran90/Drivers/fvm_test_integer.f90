
! Test Fortran_Virtual_Memory

      MODULE RANDOM
!  Portable random number generator
         IMPLICIT NONE
         PRIVATE
         PUBLIC :: RANDOM_integer, RANDOM_GET_SEED, RANDOM_INITIALIZE

         TYPE, PUBLIC :: RANDOM_seed
            PRIVATE
            INTEGER :: ix = 65535
         END TYPE
         integer, parameter :: a = 16807, b15 = 32768
         integer, parameter :: b16 = 65536, p = 2147483647
         integer, parameter :: b30 = 1073741824, q = 1073741823

      CONTAINS
!!!!!!!!!!!!!!!!!!!
         subroutine RANDOM_INTEGER ( seed, n, random )
!  Integer random number in the range [1,n] if n > 1.
!  Otherwise, the value n is returned

         TYPE (RANDOM_seed), intent( INOUT ) :: seed
         integer, intent( IN ) :: n
         integer, intent( OUT ) :: random

         INTEGER :: be1, be2, c, d, f, fhi, g, k, leftlo
         INTEGER :: mhi, mlo, mu, nu, xalo, xhi, xlo

         IF ( n > 1 ) THEN

            xhi = seed%ix / b16

!  Get 16 lo bits of seed%ix and form lo product

            xalo = ( seed%ix - xhi * b16 ) * a

!  Get 15 hi order bits of lo product

            leftlo = xalo / b16

!  Form the 31 highest bits of full product

            fhi = xhi * a + leftlo

!  Get overflopast 31st bit of full product

            k = fhi / b15

!  Assemble all the parts and presubtract P. The parentheses are essential

            seed%ix = (((xalo-leftlo*b16) - p) + (fhi-k*b15) * b16) +  k

!  Add p back in if neccessary

            IF ( seed%ix < 0 ) seed%ix = seed%ix + p

!  Multiply by n and divide by 2**31-1 in integer arithmetic.
!  Split seed%ix and n into hi and lo parts

            xhi = seed%ix / b15 ; xlo = seed%ix - b15 * xhi
            mhi = n / b15 ; mlo = n - b15 * mhi

!  Calculate intermediate product and split into hi and lo parts.
!  Presubtract p

            f = ( xhi * mlo - p ) + xlo * mhi

!  f is > 0 if intermediate product would have overflowed

            IF ( f <= 0 ) THEN
               f = f + p ; be1 = f / b15 ; be2 = f - be1 * b15
            ELSE
               f = f - 1 ; be1 = f / b15 ; be2 = f - be1 * b15; be1 = be1 + b16
            END IF

!  Form product of lo parts and add in lo part of intermediate product
!  to get lo part of complete product

            g = b15 * be2 + xlo * mlo

!  Represent lo part of full product in base 2**30

            d = g / b30 ; c = xhi / 2

!  Calculate full product divided by 2**30

            f = (( 2 * ( c * mhi - q ) - 1) + mhi * ( xhi - 2 * c )) + d + be1

!  Get full product divided in base 2**31

            IF ( f <= 0 ) THEN
               f = f + p ; nu = f / 2 ; mu = f - nu * 2
            ELSE
               f = f - 1 ; nu = f / 2 ; mu = f - 2 * nu ; nu = nu + b30
            ENDIF

!  Calculate remainder of product divided by 2**31

            f = ( b30 * mu - p ) + nu + ( g - b30 * d )
            random = nu + 1

!  Add one if remainder is not < 2**31-1
            IF ( f >= 0 ) random = random + 1
         ELSE

!  If n is less than or equal to 1, set random to n.
            random = n
         END IF
         END subroutine RANDOM_INTEGER
!!!!!!!!!!!!!!!!!!!
         subroutine RANDOM_GET_SEED ( seed, value )
!  Determine the current word generator.
         TYPE (RANDOM_seed), intent( IN ) :: seed
         integer, intent( OUT ) :: value

         value = seed%ix
         END subroutine RANDOM_GET_SEED
!!!!!!!!!!!!!!!!!!!
         subroutine RANDOM_INITIALIZE ( seed )
!  Set the word generator to its default value.
         TYPE (RANDOM_seed), intent( OUT ) :: seed

         seed%ix = 65535
         END subroutine RANDOM_INITIALIZE

      END MODULE RANDOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program test_FVM_integer

      use random
      use Fortran_Virtual_Memory_integer
      implicit none
      integer, parameter :: wp = kind(0)
      integer, parameter :: long = selected_int_kind(18)
      type (FVM_data) :: data
      type (random_seed) :: seed
      logical :: active
      integer(long) :: inactive
      integer :: i,icase,io,iflag,j,k,lpage,npage,value, &
                 l,len,loc_short,lvarray,m,nsec,part,st
      integer :: ifile(10),num_file(10)
      integer :: iopen(10) ! iopen(j)=0 if filename(j) not open, 1 if open
      integer(long):: lenw(10),loc,nrec,file_size
      character(len=12) :: filename(10)
      character(len=512) :: long_name = 'long', long_path(1)=(/'long'/)
      character(len=400) :: path(3) 
      integer(wp),allocatable:: jarray(:) ! Array for reading or writing
      integer(wp),allocatable:: karray(:) ! Old value of jarray
      integer(wp),allocatable :: array(:,:) ! Copy of all the virtual arrays
      integer,allocatable :: map(:) ! Map array for reading
      logical :: discard

      call random_initialize(seed)
      call random_get_seed(seed,value)
      npage = 8
      lpage = 10
      nrec = 50

! Set path
      path(1) = 'xtmp'
      path(2) = 'wtmp'
      path(3) = '     '
!!! The following were used to write to different directories.
!   The directories /xtmp and /wtmp were set up to have very small amounts
!   of space available and so could be used to test writing
!   files to different devices, but their inclusion makes
!   the test deck non portable. 
!      path(1) = '/xtmp/'
!      path(2) = '/wtmp/'
!      path(3) = '   '

      file_size = lpage*nrec

! Choose names for primary files
      do i = 1,10
        write (filename(i),'(a,i1,a)') 'test',i-1,'name'
        iopen(i) = 0
        lenw(i) = -1
      end do

! Let lvarray be max. length of virtual array
! we are going to allow 3 secondary files (so 4 files in total)
      nsec = 4
      lvarray = nsec*nrec*lpage

      deallocate(array,jarray,karray,map,stat=st)
      allocate(array(lvarray,12),jarray(lvarray),karray(lvarray), &
               map(lvarray),stat=st)
      jarray(1:10) = 34

! Simple test of a write followed by a read with default sizes
      call FVM_initialize(iflag,data,path)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_initialize.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      call FVM_open('test2',ifile(1),iflag,data)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_open.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      loc = 3*data%lpage+7
      call FVM_write(ifile(1),loc,10,jarray,iflag,data)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_write.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lkeep=.true.)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_close.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      call FVM_open('test2',ifile(1),iflag,data,lenw(1))
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_open.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      call FVM_read(ifile(1),loc-1,12,jarray,iflag,data)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_read.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      if(any(jarray(1:12)/=(/0,(34,i=1,10),0/))) then
         write (6,'(a)') 'Initial read error'; stop
      end if
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lkeep=.false.)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_close.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      call FVM_end(iflag,data)
      if (iflag < 0) then
         write(6,'(a,i1)') ' Failure for call to FVM_end.'
         write(6,'(a,i4,a)') ' iflag= ',iflag; stop
      end if
      lenw(1) = -1

      array(:,:) = 0.0d0
      jarray(:) = 0.0d0

! Initialize FVM
      write(6,'(3(a,i4))') ' FVM_initialize with npage=',npage, &
         ' lpage=',lpage, ' file_size=',file_size
      call FVM_initialize(iflag,data,path,npage=npage,lpage=lpage,&
           file_size=file_size)
      if (iflag /= 0) then
          write(6,'(a,i1)')' Failure for call to FVM_initialize'
          stop
      end if

! Loop making random read/writes
cases:  do icase = 1,500
         active = mod(icase,7)/=1
         if (.not.active) then
             inactive = 1 + mod(icase,5)*file_size
         end if
! if io = 1, we read; if io = 2 we write; if io = 3, we read with map
        call random_integer(seed,3,io)
        call random_integer(seed,3,l)
        if (icase==1 .or. l/=2) then
! Choose random start position and length
           call random_integer(seed,3,l)
           if (l==2) then
             call random_integer(seed,lpage,k)
           else
             call random_integer(seed,lvarray,l)
             call random_integer(seed,l,k)
           end if
           len = k
           k = lvarray - k + 1
           call random_integer(seed,k,loc_short)
           len = len - 2
           loc = loc_short
! Choose random superfile
           call random_integer(seed,10,j)
        end if

        if (io == 2) then
          do i = 1,len
            k = loc + i - 1
            jarray(i) = real(i)
            array(k,j) = real(i)
          end do
        else if (io == 3) then
          do i = 1,lvarray
            map(i) = i
            karray(i) = jarray(i)
          end do
          do i = 1,len
            call random_integer(seed,lvarray,l)
            m = map(i)
            map(i) = map(l)
            map(l) = m
          end do
        end if

        if (iopen(j) == 0) then
          call FVM_open(filename(j),ifile(j),iflag,data,lenw=lenw(j))
          write(6,'(i4,a,i3,a,i3)') icase,': FVM_open for superfile ', &
             ifile(j),','//filename(j)//' with lenw=',lenw(j)
          if (iflag < 0) then
            write(6,'(a,i1)') ' Failure for call to FVM_open. j = ',j
            write(6,'(a,i4,a)') ' iflag= ',iflag
            stop
          end if
          iopen(j) = 1
        end if

        discard = mod(icase,11)==1
        if (io == 1) then
          write(6,'(i4,a,i5,a,i5,a,i4)') icase,': FVM_read  for locations ', &
             loc,' to', loc+len-1,' on superfile',ifile(j)
          if (.not.discard) then
             call FVM_read(ifile(j),loc,len,jarray,iflag,data)
          else
             write(6,'(a)') '****** without retention **************'
	     part = len/3
             call FVM_read(ifile(j),loc,part,jarray,iflag,data,discard=discard)
             call FVM_read(ifile(j),loc+2*part,len-2*part,jarray(1+2*part:),&
                            iflag,data,discard=discard)
             call FVM_read(ifile(j),loc+part,part,jarray(1+part:),&
                            iflag,data,discard=discard)
         end if
          if (iflag < 0) then
            write(6,'(a)') ' Failure for call to FVM_read'
            stop
          end if
        else if (io == 3) then
          write(6,'(i4,a,i5,a,i5,a,i4)') icase, &
             ': FVM_read with map for locations ', &
             loc,' to', loc+len-1,' on superfile',ifile(j)
          if (.not.discard) then
             write(6,'(a,i4)') 'Not active, inactive=',inactive
             call FVM_read(ifile(j),loc,len,jarray,iflag,data,map=map)
          else
            write(6,'(a)') '****** without retention **************'
 	      part = len/3
            call FVM_read(ifile(j),loc+2*part,len-2*part,jarray,iflag,data, &
                         discard=discard,map=map(1+2*part:))
            call FVM_read(ifile(j),loc,part,jarray,iflag,data, &
                         discard=discard,map=map)
            call FVM_read(ifile(j),loc+part,part,jarray,iflag,data, &
                         discard=discard,map=map(1+part:))
          end if
          if (iflag < 0) then
            write(6,'(a)') ' Failure for call to FVM_read with map'
            stop
          end if
        else
          discard = .false.
          write(6,'(i4,a,i5,a,i5,a,i4)')icase,': FVM_write for locations ', &
             loc,' to', loc+len-1,' on superfile',ifile(j)
          call FVM_write(ifile(j),loc,len,jarray,iflag,data, &
                      inactive=inactive)
          if (iflag < 0) then
            write(6,'(a)') ' Failure for call to FVM_write'
            stop
          end if
        end if

        if (io == 2 .or. len <= 0) cycle

! Check for correct read
        if (io == 1) then
          do i = 1,len
            k = loc + i - 1
            if (array(k,j) == jarray(i)) cycle
            if (array(k,j) == -9) cycle
            write(6,'(a,i7,a,i3,a,es10.3,a,es10.3)') &
           ' Read error. Element',k,' of superfile ',j,' is ',array(k,j), &
                  ' but was read as ',jarray(I)
             stop
          end do
        else
          do i = 1,len
            k = loc + i - 1
            if ( jarray(map(i)) == karray(map(i)) + array(k,j) ) cycle
            if (array(k,j) == -9) cycle
            write(6,'(a,i7,a,es10.3,/,a,es10.3,a,es10.3)') &
           ' Read error for element',k,' jarray(map(i)) =',jarray(map(i)), &
             ' karray(map(i)) =',karray(map(i)),' array(k,j) =',array(k,j)
            stop
          end do
        end if

! Flag data read without retention
        if (discard) then
          do i = 1,len
            k = loc + i - 1
            array(k,j) = -9
          end do
        end if

        call random_integer(seed,2,i)
        if (i==1) then
! Close the superfile
          call FVM_close(ifile(j),lenw(j),num_file(j),iflag,data)
          write(6,'(i4,a,i3,a,i5,a,i4)')icase,': FVM_close for superfile', &
              ifile(j),' lenw=', lenw(j), ' num_file=', num_file(j)
          if (iflag < 0) then
            write(6,'(a,i1,a)') ' Failure for call ',j,' to FVM_close'
            stop
          end if
          iopen(j) = 0
        end if

      end do cases

! Write out contents of data
      write(6,'(a,i2)') 'data%entry =',data%entry
      write(6,'(a,i5)') 'data%file_size =',data%file_size
      write(6,'(a,i2)') 'data%iostat =',data%iostat
      write(6,'(a,i3)') 'data%lpage =',data%lpage
      write(6,'(a,i5)') 'data%ncall_read  =',data%ncall_read
      write(6,'(a,i5)') 'data%ncall_write =',data%ncall_write
      write(6,'(a,i6)') 'data%nio_read  =',data%nio_read
      write(6,'(a,i5)') 'data%nio_write =',data%nio_write
      write(6,'(a,i5)') 'data%npage = ',data%npage
      write(6,'(a,i7)') 'data%nwd_read  =',data%nwd_read
      write(6,'(a,i6)') 'data%nwd_write =',data%nwd_write
      write(6,'(a,i5)') 'data%stat =',data%stat

! FVM has finished with superfile
      do j = 1,10
        if (iopen(j) == 0) then
          call FVM_open(filename(j),ifile(j),iflag,data,lenw=lenw(j))
          if (iflag < 0) then
            write(6,'(a,i1)') ' Failure for call to FVM_open. j = ',j
            stop
          end if
        end if
        call FVM_close(ifile(j),lenw(j),num_file(j),iflag,data,lkeep=.false.)
        if (iflag < 0) then
          write(6,'(a,i1)') ' Failure for call to FVM_close. j = ',j
          stop
        end if
      end do

! Deallocate arrays
      call FVM_end(iflag,data)
      if (iflag /= 0) then
        write(6,'(a)') ' Failure for call to FVM_end.'
        stop
      end if

! Test error conditions
      write(6,'(/,a)')' Test error conditions'
      npage = 3
      lpage = 3
      nrec = 5
      file_size = lpage*nrec

! Invalid optional arguments for FVM_initialize
      call FVM_initialize(iflag,data,path,npage=0)
      call check(1,-2)
      call FVM_initialize(iflag,data,path,lpage=0)
      call check(1,-2)
      call FVM_initialize(iflag,data,path,0_long)
      call check(1,-2)

      call FVM_initialize(iflag,data,path,npage=npage,lpage=lpage, &
              file_size=file_size)

! Reopen for a superfile that does not exist
      call FVM_open(filename(1),ifile(1),iflag,data,lenw=9_long)
      call check(4,-11)

! Open for a new superfile with the name of a file that exists
      call FVM_open(filename(1),ifile(1),iflag,data)
      call check(5,0)
      call FVM_open(filename(1),ifile(2),iflag,data)
      call check(6,-12)
      write(6,'(2a)') ' Path is ', trim(path(data%iostat))

! Open with too long a name
      call FVM_open(long_name,ifile(1),iflag,data,lp=6)
      call check(7,-13)

! Write with negative loc, first with message suppressed, then try
! to an open superfile.
      loc = -1000
      call FVM_write(ifile(1),loc,len,jarray,iflag,data,lp=-6)
      call check(8,-3)
      call FVM_write(ifile(1),loc,len,jarray,iflag,data,lp=12)
      call check(9,-3)
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lp=6,&
                      lkeep=.false.)
      call check(10,0)
      loc = 1

! Try to write to a superfile that has been closed
      call FVM_write(ifile(1),loc,len,jarray,iflag,data)
      call check(11,-4)

! Write to superfile with invalid index
      call FVM_write(1111,loc,len,jarray,iflag,data)
      call check(12,-4)
      call FVM_write(0,loc,len,jarray,iflag,data)
      call check(13,-9)
      call FVM_open(filename(1),ifile(1),iflag,data)
      call check(14,0)

      loc = 100000
      len = lpage*npage*nrec
      jarray(1:len) = -1

! Read with negative loc
      loc = -1000
      call FVM_read(ifile(1),loc,len,jarray,iflag,data,lp=6)
      call check(18,-3)

! Try to read from a superfile that has been closed
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lkeep=.false.)
      call check(19,0)
      loc = 1
      call FVM_read(ifile(1),loc,len,jarray,iflag,data)
      call check(20,-4)

! Read from a superfile with invalid index
      call FVM_read(1111,loc,len,jarray,iflag,data)
      call check(21,-4)
      call FVM_read(0,loc,len,jarray,iflag,data)
      call check(21,-9)

! Try to close a superfile that has been closed
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lkeep=.false.)
      call check(25,-4)

! Close a superfile with invalid index
      call FVM_close(111,lenw(1),num_file(1),iflag,data,lkeep=.false.)
      call check(26,-4)
      call FVM_close(0,lenw(1),num_file(1),iflag,data,lkeep=.false.)
      call check(26,-9)

      if(loc>=0)then
! Read from a superfile that has not been opened
      call FVM_initialize(iflag,data,path,npage=npage,lpage=lpage, &
                           file_size=file_size)
      call check(28,0)
      call FVM_read(0,loc,len,jarray,iflag,data)
      call check(28,-9)

! Read from a superfile that has been closed
      call FVM_read(ifile(1),loc,len,jarray,iflag,data)
      call check(28,-4)
      end if

! End with an open superfile
      call FVM_open(filename(1),ifile(1),iflag,data)
      call check(29,0)
      call FVM_write(ifile(1),1_long,1,jarray,iflag,data,lp=6)
      call check(29,0)
      call FVM_end(iflag,data)
      call check(30,-10)
      call FVM_close(ifile(1),lenw(1),num_file(1),iflag,data,lkeep=.false.)
      call check(31,0)
      call FVM_end(iflag,data)
      call check(32,0)

! Initialize with too long a path name
      call FVM_initialize(iflag,data,long_path)
      call check(33,-16)

      write(6,'(a)') ' Finished'
      deallocate(array,jarray,stat=st)

contains

      subroutine check(index,value)
! Check value of iflag
        integer index,value
        if (iflag/=value) then
          write(6,'(a,i3)') ' check call',index
          write(6,'(a,i4,a,i4)') ' iflag should be',value,' but is',iflag
!         stop
      end if

      end subroutine check

      end program test_FVM_integer


