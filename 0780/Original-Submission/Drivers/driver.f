      PROGRAM DRIVER
!
!     The return variable, x, is in COMMON to make certain that no
!     overly-enthuiastic optimizer eliminates the function calls.
!
      common x
!
      n = 1000000
      call msec(mtime0)
      do i=1,n
        x = rexp_lg()
      enddo
      call msec(mtime1)
      dt_lg = 0.001*float(mtime1-mtime0)
      dta_lg = 1.E6*dt_lg/float(n)
      print *, 'Old logarithmic method:'
      print *, '  LG Run time =',dt_lg,' seconds'
      print *, '  LG Average time =',dta_lg
!
      call msec(mtime0)
      do i=1,n
        x = rexpu()
      enddo
      call msec(mtime1)
      dt_eau = 0.001*float(mtime1-mtime0)
      dta_eau = 1.E6*dt_eau/float(n)
      print *, 'Unstructured:'
      print *, '  EA Run time =',dt_eau,' seconds'
      print *, '  EA Average time =',dta_eau
!
      call msec(mtime0)
      do i=1,n
        x = rexps()
      enddo
      call msec(mtime1)
      dt_eas = 0.001*float(mtime1-mtime0)
      dta_eas = 1.E6*dt_eas/float(n)
      print *, 'Structured:'
      print *, '  EA Run time =',dt_eas,' seconds'
      print *, '  EA Average time =',dta_eas
!
      stop
      end
      SUBROUTINE MSEC(M)
      integer m, iv(8)
      call date_and_time(values=iv)     ! F90 intrinsic
      m = iv(8) + 1000*(iv(7) + 60*(iv(6) + 60*iv(5)))
      return
      end
      FUNCTION REXP_LG()                ! Old method for comparison
      call random_number(u)
      rexp_lg = -alog(u)
      return
      end
