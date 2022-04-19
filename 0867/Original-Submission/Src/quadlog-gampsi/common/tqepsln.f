      PROGRAM tqepsl
      EXTERNAL qepsln, qeps
      REAL*16          qepsln, qeps
*
      INCLUDE 'stdio.inc'
*
      WRITE (stdout,'('' EISPACK qepsln:'' ,1p,e10.2)') qepsln(1.0q0)
      WRITE (stdout,'('' Our qeps:      '' ,1p,e10.2)') qeps(1.0q0)
      END
