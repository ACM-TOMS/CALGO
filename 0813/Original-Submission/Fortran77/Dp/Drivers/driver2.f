      program spgma2

C     SPG Test driver. This master file uses SPG to solve
C     the location problem described in Test Problems
C     section. It generates the problems, calls the optimizer
C     to solve them and writes the reports (tables).
C
C     This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
C     Reformatted 03 OCT 2000 by Tim Hopkins.
C     Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

C     PARAMETERS
      integer npmax
      parameter (npmax=50000)
      integer nvsmax
      parameter (nvsmax=30)
      integer nmax
      parameter (nmax=npmax*2)

C     COMMON SCALARS
      integer np,totnvs

C     COMMON ARRAYS
      double precision edges(npmax*nvsmax*3),vert(npmax*nvsmax*2)
      integer nvs(npmax)

C     LOCAL SCALARS
      double precision pginfn,pgtwon,eps,eps2,f,prob,rmax,rmin,xstep,
     +       ystep
      real time
      integer fcnt,flag,gcnt,i,iter,m,maxfc,maxit,n,nvsvma,nvsvmi,nx,ny,
     +        pnum
      logical output

C     LOCAL ARRAYS
      double precision x(nmax)
      real dum(2)

C     EXTERNAL SUBROUTINES
      external dtime,genpro,spg

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

C     COMMON BLOCKS
      common /polyg/nvs,vert,edges,np,totnvs

C     DEFINE THE PROBLEM

1     continue

      write(*,fmt=9000)
      read (*,fmt=*,end=2) pnum,nx,ny,prob,nvsvmi,nvsvma

      if (nvsvmi.lt.3) then
          write (*,fmt=*) 'NVSVMI MUST BE GREATER THAN OR EQUAL TO 3'
          read (*,fmt=*)
      end if

      if (nvsvma.gt.nvsmax) then
          write (*,fmt=*) 'NVSVMA MUST BE LESS THAN OR EQUAL TO',nvsmax
          write (*,fmt=*) 'REDUCE NVSVMA OR INCREASE NVSMAX'
          read (*,fmt=*)
      end if

      if (nx*ny*prob.gt.npmax) then
          write (*,fmt=*) 'NX*NY*PROB MUST BE LESS THAN ',npmax
          write (*,fmt=*) 'REDUCE NX, NY OR PROB OR INCREASE NPMAX'
          read (*,fmt=*)
      end if

      xstep = 5.0d0
      ystep = 5.0d0

      rmin = 1.0d0
      rmax = 2.0d0

C     GENERATE THE PROBLEM

      call genpro(nx,ny,xstep,ystep,prob,nvsvmi,nvsvma,rmin,rmax)

      write (*,fmt=*) 'NUMBER OF POLYGONS: ',np
      write (*,fmt=*) 'NUMBER OF VERTICES (AND EDGES): ',totnvs

C     DEFINE INITIAL POINT

      n = 2*np
      do i = 1,n
          x(i) = 0.0d0
      end do

C     SET UP THE INPUT DATA OF THE OPTIMIZATION ALGORITHM

      output= .false.
      maxit = 1000
      maxfc = 2000
      eps = 0.0d0
      eps2 = 1.0d-06
      m = 10

C     CALL THE OPTIMIZER

      call dtime(dum)
      call spg(n,x,m,eps,eps2,maxit,maxfc,output,f,pginfn,pgtwon,iter,
     +         fcnt,gcnt,flag)
      call dtime(dum)

      time = dum(1) + dum(2)

C     WRITE STATISTICS

      write (*,fmt=9010) f,pginfn,sqrt(pgtwon),flag
      write (*,fmt=9020) iter,fcnt,gcnt
      write (*,fmt=9030) time

C     WRITE SOLUTION

      goto 1
2     continue

      stop

 9000 format (//'---------------------------------------'//)
 9010 format (/' F= ',8X,1P,D17.10,/' PGINFNORM= ',1X,1P,D16.10,
     +       /' PGTWONORM= ',1X,1P,D16.10,/' FLAG= ',6X,I1)
 9020 format (/' ITER= ',1X,I10,/' FCNT= ',1X,I10,/' GCNT= ',1X,I10)
 9030 format (/' TIME= ',F12.2,' SECONDS')
      end


      subroutine evalf(n,x,f,inform)

C     This subroutine computes the objective function.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point at which the function will be evaluated.
C
C     On Return
C
C     f     double precision,
C           function value at x,
C
C     inform integer,
C           termination parameter:
C           0= the function was successfully evaluated,
C           1= some error occurs in the function evaluation.

C     SCALAR ARGUMENTS
      double precision f
      integer n,inform

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      double precision diff1,diff2,dist
      integer i,ndist

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

      inform = 0

      f = 0.0d0
      ndist= n/2-1
      do i = 1,ndist
          diff1 = x(1) - x(2*i+1)
          diff2 = x(2) - x(2*i+2)
          dist = sqrt(diff1**2+diff2**2)
          if (dist.le.1.0d-4) then
              write (*,fmt=*)
     +          'ERROR IN PROBLEM DEFINITION (DIST TOO SMALL)'
              inform = 1
          end if
          f = f + dist
      end do

      f = f / ndist

      end


      subroutine evalg(n,x,g,inform)

C     This subroutine computes the gradient of the
C     objective function.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point at which the gradient will be evaluated.
C
C     On Return
C
C     g     double precision g(n),
C           gradient vector at x,
C
C     inform integer,
C           termination parameter:
C           0= the gradient was successfully evaluated,
C           1= some error occurs in the gradient evaluation.

C     SCALAR ARGUMENTS
      integer n,inform

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      double precision diff1,diff2,dist
      integer i,ndist

C     INTRINSIC FUNCTIONS
      intrinsic sqrt

      inform = 0

      g(1) = 0.0d0
      g(2) = 0.0d0
      ndist= n/2-1
      do i = 1,ndist
          diff1 = x(1) - x(2*i+1)
          diff2 = x(2) - x(2*i+2)
          dist = sqrt(diff1**2+diff2**2)
          g(2*i+1) = - (diff1/dist) / ndist
          g(2*i+2) = - (diff2/dist) / ndist
          g(1) = g(1) - g(2*i+1)
          g(2) = g(2) - g(2*i+2)
      end do

      end


      subroutine proj(n,x,inform)

C     This subroutine computes the projection of an arbitrary 
C     point onto the feasible set. Since the feasible set can be 
C     described in many ways, its information is in a common 
C     block.
C
C     In this particular implementation of the projection subroutine
C     for the location problem, each point z_i in R^2 (stored at 
C     positions 2i-1 and 2i of x) must be projected onto polygon P_i.
C     See the Test Problems section for details.
C
C     On Entry:
C
C     n     integer,
C           size of the problem,
C
C     x     double precision x(n),
C           point that will be projected.
C
C     On Return
C
C     x     double precision x(n),
C           projected point,
C
C     inform integer,
C           termination parameter:
C           0= the projection was successfully done,
C           1= some error occurs in the projection.

C     PARAMETERS
      integer npmax
      parameter (npmax=50000)
      integer nvsmax
      parameter (nvsmax=30)

C     SCALAR ARGUMENTS
      integer n,inform

C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON SCALARS
      integer np,totnvs

C     COMMON ARRAYS 
      double precision edges(npmax*nvsmax*3),vert(npmax*nvsmax*2)
      integer nvs(npmax)

C     LOCAL SCALARS
      double precision xproj,yproj
      integer base,i

C     EXTERNAL SUBROUTINES
      external projp

C     COMMON BLOCKS
      common /polyg/nvs,vert,edges,np,totnvs

      inform = 0

      base = 0

      do i = 1,np

        ! PROJECT z_i ONTO P_i

          call projp(i,base,x(2*i-1),x(2*i),xproj,yproj,inform)

          x(2*i-1) = xproj
          x(2*i) = yproj

          base = base + nvs(i)

      end do

      end


      subroutine projp(p,base,x,y,xproj,yproj,inform)

C     This subroutine complements proj subroutine projecting
C     a point z_i in R^2 onto its corresponding polygon P_i.

C     On Entry:
C
C     p     integer,
C           index of the point z_p to be projected,
C
C     base  integer,
C           index which indicates the positions where the edges 
C           of the polygon  P_p are stored inside the "edges" 
C           vector which is part of the polygons description.
C
C     x     double precision,
C           first coordinate of z_i, i.e., z_i^1,
C
C     y     double precision,
C           second coordinate of z_i, i.e., z_i^2,
C
C     On Return
C
C     xproj double precision,
C           the projection of x,
C
C     yproj double precision,
C           the projection of y.
C
C     inform integer,
C           termination parameter:
C           0= the projection was successfully done,
C           1= some error occurs in the projection.

C     PARAMETERS
      integer npmax
      parameter (npmax=50000)
      integer nvsmax
      parameter (nvsmax=30)

C     SCALAR ARGUMENTS
      double precision x,xproj,y,yproj
      integer base,inform,p

C     COMMON SCALARS
      integer np,totnvs

C     COMMON ARRAYS
      double precision edges(npmax*nvsmax*3),vert(npmax*nvsmax*2)
      integer nvs(npmax)

C     LOCAL SCALARS
      double precision a,b,c,dist,mindis,nearvx,nearvy,px,py,v1x,v1y,
     +                 v2x,v2y,vx,vy
      integer i
      logical fact

C     LOCAL ARRAYS
      logical satisf(nvsmax)

C     INTRINSIC FUNCTIONS
      intrinsic max,min,sqrt

C     COMMON BLOCKS
      common /polyg/nvs,vert,edges,np,totnvs

C     TEST ALL EDGES OF THE POLYGON FOR SATISFIABILITY

      fact = .true.
      do i = base + 1,base + nvs(p)
          a = edges(3*i-2)
          b = edges(3*i-1)
          c = edges(3*i)
          if (a*x+b*y+c.le.0.0d0) then
              satisf(i-base) = .true.

          else
              fact = .false.
              satisf(i-base) = .false.
          end if
      end do

      if (fact) then
          xproj = x
          yproj = y
          return

      end if

C     SEARCH FOR THE CLOSEST VERTEX

      mindis = 1.0d+99
      do i = base + 1,base + nvs(p)
          vx = vert(2*i-1)
          vy = vert(2*i)
          dist = sqrt((vx-x)**2+ (vy-y)**2)
          if (dist.le.mindis) then
              mindis = dist
              nearvx = vx
              nearvy = vy
          end if

      end do

      xproj = nearvx
      yproj = nearvy

C     PROJECT ONTO THE VIOLATED CONSTRAINTS

      do i = base + 1,base + nvs(p)
          if (.not.satisf(i-base)) then
              a = edges(3*i-2)
              b = edges(3*i-1)
              c = edges(3*i)
              if (a.ne.0.0d0) then
                  py = (y-b* (c/a+x)/a)/ (b**2/a**2+1)
                  px = - (c+b*py)/a

              else
                  if (b.ne.0.0d0) then
                      px = (x-a* (c/b+y)/b)/ (a**2/b**2+1)
                      py = - (c+a*px)/b

                  else
                      write (*,fmt=*)
     +                  'ERROR IN PROBLEM DEFINITION (a=b=0 for a ',
     +                  'constraint of the type a x + b y + c = 0)'
                      inform= 1
                  end if

              end if

              v1x = vert(2*i-1)
              v1y = vert(2*i)
              if (i.ne.base+nvs(p)) then
                  v2x = vert(2*i+1)
                  v2y = vert(2*i+2)

              else
                  v2x = vert(2*base+1)
                  v2y = vert(2*base+2)
              end if

              if (min(v1x,v2x).le.px .and. px.le.max(v1x,v2x) .and.
     +            min(v1y,v2y).le.py .and. py.le.max(v1y,v2y)) then
                  dist = sqrt((px-x)**2+ (py-y)**2)
                  if (dist.le.mindis) then
                      mindis = dist
                      xproj = px
                      yproj = py
                  end if

              end if

          end if

      end do

      end


      subroutine genpro(nx,ny,xstep,ystep,prob,nvsvmi,nvsvma,rmin,rmax)

C     This subroutine generates a location problem (see the Figure).
C     First, a regular grid with nx horizontal points and ny vertical 
C     points in the positive orthant is considered. The points of the 
C     grid start at the origin with an horizontal distance of xstep and 
C     a vertical distance of ystept. This grid will be the space at which
C     the cities (represented by polygons) will be distributed. Before
C     building the cities, an area (rectangle) of preservation where 
C     almost nothing can be done, is defined. This area of preservation 
C     will receive, after the construction of the cities, an hydraulic 
C     plant of energy generation (to supply the energy to the cities). Then, 
C     in the rest of the space, the cities are built. At each point of the 
C     grid (out of the central region) a city (represented by a polygon)
C     will be built with probability prob. The definition of the polygon
C     uses variables nvsvmi, nvsvma, rmin and rmax in a way described in
C     the genpol (generate polygon) subroutine. To transmit the energy
C     from the plant to the cities, a tower inside each city and a tower
C     inside the central region must be built. The objective of the
C     problem is to determine the location of this towers in order 
C     to minimize the sum of the distances from each city tower to the
C     central one. 
C
C     On Entry:
C
C     nx    integer,
C           number of horizontal points in the grid,
C
C     ny    integer,
C           number of vertical points in the grid,
C
C     xstep double precision,
C           horizontal distance between points of the grid,
C
C     ystep double precision,
C           vertical distance between points of the grid,
C
C     prob  double precision, 
C           probability of defining a city at point of the grid
C           (0 <= prob <= 1),
C
C     nvsvmi integer,
C           parameter for the polygon generation (described
C           in genpol subroutine),
C
C     nvsvma integer,
C           parameter for the polygon generation (described
C           in genpol subroutine),
C
C     rmin  double precision,
C           parameter for the polygon generation (described
C           in genpol subroutine),
C
C     rmax  double precision,
C           parameter for the polygon generation (described
C           in genpol subroutine).
C
C     On output:
C
C     As described in the genpol subroutine, the output is
C     saved in the polyg common block.


C     PARAMETERS
      integer npmax
      parameter (npmax=50000)
      integer nvsmax
      parameter (nvsmax=30)

C     SCALAR ARGUMENTS
      double precision prob,rmax,rmin,xstep,ystep
      integer nvsvma,nvsvmi,nx,ny

C     COMMON SCALARS
      integer np,totnvs

C     COMMON ARRAYS
      double precision edges(npmax*nvsmax*3),vert(npmax*nvsmax*2)
      integer nvs(npmax)

C     LOCAL SCALARS
      double precision cx,cy,lx,ly,ux,uy
      integer i,j,seed

C     EXTERNAL FUNCTIONS
      real ran
      external ran

C     EXTERNAL SUBROUTINES
      external genpol

C     COMMON BLOCKS
      common /polyg/nvs,vert,edges,np,totnvs

C     SEED FOR THE RANDOM GENERATION

      seed = 760013

C     DEFINE BOX CONSTRAINTS FOR THE CENTRAL POINT

      lx = 0.40d0* (nx-1)*xstep
      ux = 0.60d0* (nx-1)*xstep
      ly = 0.40d0* (ny-1)*ystep
      uy = 0.60d0* (ny-1)*ystep

C     DEFINE CENTRAL POINT POLYGON

      nvs(1) = 4

      vert(1) = lx
      vert(2) = ly
      vert(3) = lx
      vert(4) = uy
      vert(5) = ux
      vert(6) = uy
      vert(7) = ux
      vert(8) = ly

      edges(1) = -1.0d0
      edges(2) = 0.0d0
      edges(3) = lx
      edges(4) = 0.0d0
      edges(5) = 1.0d0
      edges(6) = -uy
      edges(7) = 1.0d0
      edges(8) = 0.0d0
      edges(9) = -ux
      edges(10) = 0.0d0
      edges(11) = -1.0d0
      edges(12) = ly

C     DEFINE CITY-POLYGONS CENTERED AT THE GRID POINTS AND OUTSIDE
C     THE CENTRAL REGION

      np = 1
      totnvs = 4

      do i = 0,nx - 1
          do j = 0,ny - 1

              cx = i*xstep
              cy = j*ystep

              if ((cx.lt.lx-xstep.or.cx.gt.ux+xstep.or.
     +            cy.lt.ly-ystep.or.cy.gt.uy+ystep) .and.
     +            ran(seed).le.prob) then

                ! GENERATE A NEW POLYGON

                  np = np + 1
                  call genpol(cx,cy,nvsvmi,nvsvma,rmin,rmax,seed)
                  totnvs = totnvs + nvs(np)

              end if

          end do
      end do

      end


      subroutine genpol(cx,cy,nvsvmi,nvsvma,rmin,rmax,seed)

C     This subroutine generates a polygon in R^2 with its 
C     vertices in a sphere centered at point (cx,cy). The
C     number of vertices is randomly generated 
C     satisfying nvsvmi <= number of vertices  <= 
C     nvsvma. The ratio of the sphere is also randomly 
C     generated satisfying rmin <= ratio < rmax. The 
C     generated polygon is stored in the common block "polyg".
C 
C     On Entry:
C
C     cx    double precision,
C           first coordinate of the center of the sphere,
C
C     cy    double precision,
C           second coordinate of the center of the sphere,
C
C     nvsvmi integer,
C           minimum number of vertices,
C
C     nvsvma integer,
C           maximum number of vertices,
C
C     rmin double,
C           minimum ratio of the sphere,
C
C     rmax double,
C           maximum ratio of the sphere,
C
C     seed double,
C           seed for the random generation.
C
C     On Output:
C     
C     The generated polygon is stored in the polyg common block
C     described below.
C
C     Common block polyg:
C
C     common /polyg/nvs,vert,edges,np,totnvs
C 
C     This structure represents, at any time, np polygons.
C     Position i of array nvs indicates the number of vertices
C     of polygon i. Arrays vert and edges store the 
C     vertices and edges of the polygons.
C
C     For example, if nvs(1) = 3 it indicates that the first
C     polygon has 3 vertices (edges). Then, if the vertices
C     are (x1,y1), (x2,y2) and (x3,y3), we have that vert(1) 
C     = x1, vert(2) = y1, vert(3) = x2, vert(4) = y2, vert(5) 
C     = x3, and vert(6) = y3. And, if the edges (written as
C     ax + by + c = 0) are a1 x + b1 y + c1 = 0, a2 x + b2 y 
C     + c2 = 0, and a3 x + b3 y + c3 = 0 then edges(1) = a1,
C     edges(2) = b1, edges(3) = c1, edges(4) = a2, edges(5) = 
C     b2, edges(6) = c2, edges(7) = a3, edges(8) = b3 and
C     edges(9) = c3.
C
C     totnvs indicates the total number of vertices  
C     of the set of polygons. This information is used when
C     a new polygon is created to know the first free position
C     of arrays vert and edges at which the vertices and edges 
C     of the new polygon will be saved.
C
C     Two additional details: 
C
C     1) For each polygon, the vertices are ordered clockwise 
C     and edge i corresponds to the edge between vertices 
C     i and i+1 (0 if i=n).
C
C     2) For each edge of the form ax + bx + c = 0, constants
C     a, b and c are chosen in such a way that 
C     (|a| = 1 or |b| = 1) and (a cx + b cy + c <= 0).

C     PARAMETERS
      integer npmax
      parameter (npmax=50000)
      integer nvsmax
      parameter (nvsmax=30)
      double precision pi
      parameter (pi=3.1415926535898d0)

C     SCALAR ARGUMENTS
      double precision cx,cy,rmax,rmin
      integer nvsvma,nvsvmi,seed

C     COMMON SCALARS
      integer np,totnvs

C     COMMON  ARRAYS
      double precision edges(npmax*nvsmax*3),vert(npmax*nvsmax*2)
      integer nvs(npmax)

C     LOCAL SCALARS
      double precision r
      integer i

C     LOCAL ARRAYS
      double precision angl(nvsmax)

C     EXTERNAL FUNCTIONS
      real ran
      external ran

C     EXTERNAL SUBROUTINES
      external class,constr

C     INTRINSIC FUNCTIONS
      intrinsic dcos,dsin,int

C     COMMON BLOCKS
      common /polyg/nvs,vert,edges,np,totnvs

C     GENERATE THE NUMBER OF VERTICES 

      nvs(np) = nvsvmi + int((nvsvma-nvsvmi+1)*ran(seed))

C     GENERATE THE RATIO OF THE SPHERE

      r = rmin + (rmax-rmin)*ran(seed)

C     GENERATE ALL ANGLES SATISFYING 0 <= ANGLE_I < 2*PI

      do i = 1,nvs(np)
          angl(i) = 2*pi*ran(seed)
      end do

C     CLASSIFY THE ANGLES IN DECREASING ORDER

      call class(nvs(np),angl)

C     CONSTRUCT THE VERTICES

      do i = 1,nvs(np)
          vert(2* (totnvs+i)-1) = cx + r*dcos(angl(i))
          vert(2* (totnvs+i)) = cy + r*dsin(angl(i))
      end do

C     CONSTRUCT THE EDGES

      do i = totnvs + 1,totnvs + nvs(np) - 2
          call constr(vert(2*i-1),vert(2*i),vert(2*i+1),vert(2*i+2),
     +                vert(2*i+3),vert(2*i+4),edges(3*i-2),edges(3*i-1),
     +                edges(3*i))
      end do

      i = totnvs + nvs(np) - 1

      call constr(vert(2*i-1),vert(2*i),vert(2*i+1),vert(2*i+2),
     +            vert(2*totnvs+1),vert(2*totnvs+2),edges(3*i-2),
     +            edges(3*i-1),edges(3*i))

      i = totnvs + nvs(np)

      call constr(vert(2*i-1),vert(2*i),vert(2*totnvs+1),
     +            vert(2*totnvs+2),vert(2*totnvs+3),vert(2*totnvs+4),
     +            edges(3*i-2),edges(3*i-1),edges(3*i))

      end


      subroutine constr(x1,y1,x2,y2,x3,y3,a,b,c)

C     This subroutine computes the real constants a, b and c of 
C     the straight line ax + by + c = 0 in R^2 defined by the 
C     points (x1,y1) and (x2,y2); such that the point (x3,y3) 
C     satisfies the constraint ax + by + c <= 0.
C
C     On Entry:
C
C     x1    double precision,
C           first coordinate of point (x1,y1),
C
C     y1    double precision,
C           second coordinate of point (x1,y1),
C
C     x2    double precision,
C           first coordinate of point (x2,y2),
C
C     y2    double precision,
C           second coordinate of point (x2,y2),
C
C     x3    double precision,
C           first coordinate of point (x3,y3),
C
C     y3    double precision,
C           second coordinate of point (x3,y3).
C
C     On Return
C
C     a,b,c double precision
C           the desired constants.

C     SCALAR ARGUMENTS
      double precision a,b,c,x1,x2,x3,y1,y2,y3

      if (x1.eq.x2 .and. y1.eq.y2) then
          write (*,fmt=*)
     +      'ERROR IN FUNCTION CONSTRAINT: X1=X2 AND Y1=Y2'
      end if

      if (y1.ne.y2) then
          a = 1.0d0
          b = - (x2-x1)/ (y2-y1)
          c = - (x1+b*y1)

      else
          a = 0.0d0
          b = 1.0d0
          c = -y1
      end if

      if (a*x3+b*y3+c.gt.0.0d0) then
          a = -a
          b = -b
          c = -c
      end if

      end


      subroutine class(n,x)

C     This subroutine classifies the elements of a vector in
C     decreasing order, i.e., on output: x(1) >= x(2) >= 
C     ... >= x(n).
C
C     On Entry:
C
C     n     integer,
C           number of elements of the vector to be classified,
C
C     x     double precision x(n),
C           vector to be classified.
C
C     On Return
C
C     x     double precision x(n),
C           classified vector.

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     LOCAL SCALARS
      double precision aux,xmax
      integer i,j,pos

      do i = 1,n

          xmax = x(i)
          pos = i

          do j = i + 1,n
              if (x(j).gt.xmax) then
                  xmax = x(j)
                  pos = j
              end if
          end do

          if (pos.ne.i) then
              aux = x(i)
              x(i) = x(pos)
              x(pos) = aux
          end if

      end do

      end













