C     THIS FILE CONTAINS THE ROUTINES THAT CONSTRUCT AND 
C     APPLY GIVENS ROTATIONS IN THE FOUR ARITHMETICS
C     SROT, SROTG, DROT, DROTG ARE EXTRACTED FROM THE BLAS LIBRARY 
C     CROT, ZROT ARE EXTRACTED FROM THE LAPACK LIBRARY
C
      subroutine srot (n,sx,incx,sy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sx(*),sy(*),stemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c     code for unequal increments or equal increments not equal
c     to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp  = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1
c
   20 do 30 i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
   30 continue
      return
      end
C
C
      subroutine srotg(sa,sb,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      real sa,sb,c,s,roe,scale,r,z
c
      roe = sb
      if( abs(sa) .gt. abs(sb) ) roe = sa
      scale = abs(sa) + abs(sb)
C
      if( scale .ne. 0.0 ) go to 10
C
      c = 1.0
      s = 0.0
      r = 0.0
      z = 0.0
      go to 20
C
   10 r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
      r = sign(1.0,roe)*r
      c = sa/r
      s = sb/r
      z = 1.0
      if( abs(sa) .gt. abs(sb) ) z = s
      if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0 ) z = 1.0/c
C
   20 sa = r
      sb = z
      return
      end
C
C
      subroutine  drot (n,dx,incx,dy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c     code for unequal increments or equal increments not equal
c     to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
C
C
      subroutine drotg(da,db,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,db,c,s,roe,scale,r,z
c
      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
C
      if( scale .ne. 0.0d0 ) go to 10
C
      c = 1.0d0
      s = 0.0d0
      r = 0.0d0
      z = 0.0d0
      go to 20
C
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
      z = 1.0d0
      if( dabs(da) .gt. dabs(db) ) z = s
      if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
C
   20 da = r
      db = z
      return
      end
C
C
      subroutine crot( n, cx, incx, cy, incy, c, s )
C
C  -- lapack auxiliary routine (version 3.0) --
C     univ. of tennessee, univ. of california berkeley, nag ltd.,
C     courant institute, argonne national lab, and rice university
C     october 31, 1992
C
C     .. scalar arguments ..
      integer            incx, incy, n
      real               c
      complex            s
C     ..
C     .. array arguments ..
      complex            cx( * ), cy( * )
C     ..
C
C  purpose
C  =======
C
C  crot   applies a plane rotation, where the cos (c) is real and the
C  sin (s) is complex, and the vectors cx and cy are complex.
C
C  arguments
C  =========
C
C  n       (input) integer
C          the number of elements in the vectors cx and cy.
C
C  cx      (input/output) complex array, dimension (n)
C          on input, the vector x.
C          on output, cx is overwritten with c*x + s*y.
C
C  incx    (input) integer
C          the increment between successive values of cy.  incx <> 0.
C
C  cy      (input/output) complex array, dimension (n)
C          on input, the vector y.
C          on output, cy is overwritten with -conjg(s)*x + c*y.
C
C  incy    (input) integer
C          the increment between successive values of cy.  incx <> 0.
C
C  c       (input) real
C  s       (input) complex
C          c and s define a rotation
C             [  c          s  ]
C             [ -conjg(s)   c  ]
C          where c*c + s*conjg(s) = 1.0.
C
C =====================================================================
C
C     .. local scalars ..
      integer            i, ix, iy
      complex            stemp
C     ..
C     .. intrinsic functions ..
      intrinsic          conjg
C     ..
C     .. executable statements ..
C
      if( n.le.0 )return
      if( incx.eq.1 .and. incy.eq.1 )go to 20
C
C     code for unequal increments or equal increments not equal 
C     to 1
C
      ix = 1
      iy = 1
      if( incx.lt.0 )ix = ( -n+1 )*incx + 1
      if( incy.lt.0 )iy = ( -n+1 )*incy + 1
      do 10 i = 1, n
         stemp  = c*cx(ix) + s*cy(iy)
         cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
         cx(ix) = stemp
         ix = ix + incx
         iy = iy + incy
   10 continue
      return
C
C     code for both increments equal to 1
C
   20 do 30 i = 1, n
         stemp = c*cx(i) + s*cy(i)
         cy(i) = c*cy(i) - conjg(s)*cx(i)
         cx(i) = stemp
   30 continue
      return
      end
C
C
      subroutine crotg(ca,cb,c,s)
      complex ca,cb,s
      real c
      real norm,scale
      complex alpha
C
      if (cabs(ca) .ne. 0.) go to 10
C
      c = 0.
      s = (1.,0.)
      ca = cb
      go to 20
C
   10 continue
      scale = cabs(ca) + cabs(cb)
      norm = scale * sqrt((cabs(ca/scale))**2 + (cabs(cb/scale))**2)
      alpha = ca /cabs(ca)
      c = cabs(ca) / norm
      s = alpha * conjg(cb) / norm
      ca = alpha * norm
C
   20 continue
      return
      end
C
C
      subroutine zrot( n, cx, incx, cy, incy, c, s )
C
C  -- lapack auxiliary routine (version 3.0) --
C     univ. of tennessee, univ. of california berkeley, nag ltd.,
C     courant institute, argonne national lab, and rice university
C     october 31, 1992
C
C     .. scalar arguments ..
      integer            incx, incy, n
      double precision   c
      complex*16         s
C     ..
C     .. array arguments ..
      complex*16         cx( * ), cy( * )
C     ..
C
C  purpose
C  =======
C
C  zrot   applies a plane rotation, where the cos (c) is real and the
C  sin (s) is complex, and the vectors cx and cy are complex.
C
C  arguments
C  =========
C
C  n       (input) integer
C          the number of elements in the vectors cx and cy.
C
C  cx      (input/output) complex*16 array, dimension (n)
C          on input, the vector x.
C          on output, cx is overwritten with c*x + s*y.
C
C  incx    (input) integer
C          the increment between successive values of cy.  incx <> 0.
C
C  cy      (input/output) complex*16 array, dimension (n)
C          on input, the vector y.
C          on output, cy is overwritten with -conjg(s)*x + c*y.
C
C  incy    (input) integer
C          the increment between successive values of cy.  incx <> 0.
C
C  c       (input) double precision
C  s       (input) complex*16
C          c and s define a rotation
C             [  c          s  ]
C             [ -conjg(s)   c  ]
C          where c*c + s*conjg(s) = 1.0.
C
C =====================================================================
C
C     .. local scalars ..
      integer            i, ix, iy
      complex*16         stemp
C     ..
C     .. intrinsic functions ..
      intrinsic          dconjg
C     ..
C     .. executable statements ..
C
      if( n.le.0 )return
      if( incx.eq.1 .and. incy.eq.1 )go to 20
C
C     code for unequal increments or equal increments not equal to 1
C
      ix = 1
      iy = 1
      if( incx.lt.0 )ix = ( -n+1 )*incx + 1
      if( incy.lt.0 )iy = ( -n+1 )*incy + 1
      do 10 i = 1, n
         stemp  = c*cx(ix) + s*cy(iy)
         cy(iy) = c*cy(iy) - dconjg(s)*cx(ix)
         cx(ix) = stemp
         ix = ix + incx
         iy = iy + incy
   10 continue
      return
C
C     code for both increments equal to 1
C
   20 continue
      do 30 i = 1, n
         stemp = c*cx(i) + s*cy(i)
         cy(i) = c*cy(i) - dconjg(s)*cx(i)
         cx(i) = stemp
   30 continue
      return
      end
C
C
      subroutine zrotg(ca,cb,c,s)
      double complex ca,cb,s
      double precision c
      double precision norm,scale
      double complex alpha
C
      if (cdabs(ca) .ne. 0.0d0) go to 10
C
      c = 0.0d0
      s = (1.0d0,0.0d0)
      ca = cb
      go to 20
C
   10 continue
      scale = cdabs(ca) + cdabs(cb)
      norm = scale*dsqrt((cdabs(ca/dcmplx(scale,0.0d0)))**2 +
     &                      (cdabs(cb/dcmplx(scale,0.0d0)))**2)
      alpha = ca /cdabs(ca)
      c = cdabs(ca) / norm
      s = alpha * dconjg(cb) / norm
      ca = alpha * norm
C
   20 continue
      return
      end
