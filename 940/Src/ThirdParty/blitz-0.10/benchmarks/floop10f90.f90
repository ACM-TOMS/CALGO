
      SUBROUTINE floop10_F90(N, x, a, b, c, u)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), u

      x = u+a+b+c
      RETURN
      END


      SUBROUTINE floop10_F90Overhead(N, x, a, b, c, u)
      INTEGER i, N
      REAL*4 x(N), a(N), b(N), c(N), u

      RETURN
      END
