! Simple driver program for alg 421
! Just computes a few values and prints them along with
! values published in Abramowitz and Stegun
      program test
      double precision carg(2), cans(2), ctrue(2)
      real error, err1
      integer lf0

10    read(*,*,end=20)carg(1), carg(2), error, lf0
      read(*,*)ctrue(1), ctrue(2)
      err1 = error
      call cdlgam(carg, cans, err1, lf0)
      IF (lf0 == 0) THEN
	write(*,'(" loggamma(",e14.4,",",e14.4,") = (",e16.8,'//&
		  &'",",e16.8,")")')carg(1:2),cans(1:2)
	write(*,'("    error requested: ",e14.4,'//&
		  &'" error estimate returned: ", e14.4)') error, err1
	write(*,'("    A and S value: = (",e16.8,",",e16.8,")")')ctrue(1:2)
	write(*,'()')
      ELSE
	write(*,'(" gamma(",e14.4,",",e14.4,") = (",e16.8,'//&
		  &'",",e16.8,")")')carg(1:2),cans(1:2)
	write(*,'("    error requested: ",e14.4,'//&
		  &'" error estimate returned: ", e14.4)') error, err1
	write(*,'("    A and S value: = (",e16.8,",",e16.8,")")')ctrue(1:2)
	write(*,'()')
      ENDIF
      goto 10
20    end
