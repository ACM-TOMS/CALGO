      PROGRAM tepsln
      EXTERNAL depsln, deps
      DOUBLE PRECISION depsln, deps
*
      INCLUDE 'stdio.inc'
*
      WRITE (stdout,'('' EISPACK depsln:'' ,1p,e10.2)') depsln(1.0d0)
      WRITE (stdout,'('' Our deps:      '' ,1p,e10.2)') deps(1.0d0)
      END
