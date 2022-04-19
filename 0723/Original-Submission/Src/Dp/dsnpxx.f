c--D replaces "?": ?SNPXX, ?COSPX, ?SINPX
c>> 1996-01-08 DSNPXX WV Snyder Original code
      double precision function DSNPXX (X)

c SIN(PI * X * X / 2) carefully to avoid loss of precision for large X

      double precision X

c DSINPX is used to compute SIN(PI * X)

      double precision D1MACH, DCOSPX, DSINPX
      external D1MACH, DCOSPX, DSINPX

c BIGX = 1 / round-off = biggest integer exactly representable by F.P.
c    If X > BIGX then to the working precision x**2 is an integer (which
c    we assume to be a multiple of four), so sin(pi/2 * x**2) = 0.
c N = [X], and later [K F]
c F = X - N = fractional part of X
c K = [ N / 2 ]
c J = N mod 2
c G = K F - M = fractional part of K F

      integer J, K, N
      double precision BIGX, F, G
      save BIGX
      data BIGX /-1.0/

      if (bigx .lt. 0.0d0) bigx = 1.0d0 / d1mach(4)
      f = abs(x)
      if (f .gt. bigx) then
c       Assume x is an even integer.
        dsnpxx = 0.0d0
        return
      endif
      n = f
      f = f - n
      k = n / 2
      j = mod(n,2)
      g = k * f
      n = g
      g = g - n
      if (j .eq. 0) then
        dsnpxx = dsinpx(0.5d0*f*f + g + g)
      else
        dsnpxx = dcospx(0.5d0*f*f + f + g + g)
      endif
      return
      end
