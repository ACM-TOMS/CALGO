      PROGRAM tisdin
      DOUBLE PRECISION inf
      LOGICAL isdinf, result
      INTEGER iinf(2)
      EQUIVALENCE (inf, iinf(1))
*
      INCLUDE 'stdio.inc'
*
      inf = 0.0d+00
      inf = 1.0d+00 / inf
*
      WRITE (stdout,*) 'Test of isdinf()'
*
      result =  isdinf(inf)
      WRITE (stdout,*) 'isdinf(', inf, ') = ', result
      WRITE (stdout,'(1x, 2z9.8)') iinf
*
      END


      LOGICAL FUNCTION isdinf(x)
************************************************************************
*     (Is x infinite?)
*     Return .TRUE. if x is infinite, and .FALSE. otherwise.
*     [12-Jun-2000]
************************************************************************
      DOUBLE PRECISION    x
*
*     Parameter variables
*
      DOUBLE PRECISION    zero
      PARAMETER           (zero = 0.0d+00)
*
      isdinf = (x .NE. zero) .AND. ((x + x) .EQ. x)
*
      END
