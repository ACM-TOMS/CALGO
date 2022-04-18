      subroutine bacolr(t0, tout, atol, rtol, npde, kcol, nintmx, nint, 
     &                  x, mflag, rpar, lrp, ipar, lip, cpar, lcp, y, 
     &                  idid)

c-----------------------------------------------------------------------
c Purpose:
c       The purpose of BACOLR is to solve NPDE dimensional systems of 
c       second order parabolic partial differential equations (PDEs) 
c       in one space variable of the form:
c
c            dU    
c            -- (t,x) = f ( t, x, U(t,x), U_x (t,x), U_{xx} (t,x) ) ,
c            dt  
c
c       where xa < x < xb and t > t0, with initial conditions at
c       time t = t0 given by:
c
c                          u(t0,x) = u_0(x),
c
c       for x_a <= x <= x_b, subject to separated boundary conditions 
c       given by:
c
c                   b_{xa} ( t, U(t,xa), U_x (t,xa) ) = 0,
c
c                   b_{xb} ( t, U(t,xb), U_x (t,xb) ) = 0,
c
c       for t > t0 and x = x_a, x = x_b, respectively.
c         
c       Guide to the above notation:
c          dU
c          -- (t,x)     denotes the first partial derivative of U(t,x)
c          dt           with respect to the time variable t.
c
c          U_x (t,x)    denotes the first partial derivative of U(t,x)
c                       with respect to space variable x.
c
c          U_{xx} (t,x) denotes the second partial derivative of U(t,x)
c                       with respect to space variable x.
c
c       Also, the above functions are NPDE dimensional vector functions.
c
c       BACOLR is a method of lines algorithm which uses bspline 
c       collocation to discretize the spatial domain [x_a,x_b].
c       The output is a vector of bspline coefficients which
c       can be used to calculate the approximate solution U(t,x) and
c       its spatial derivatives at (tout,x) where x_a <= x <= x_b
c       and t0 < tout.
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Setup of BACOLR:
c       BACOLR requires that the user specifies the system of PDEs and
c       the related initial and boundary conditions as also sets
c       input parameters (which define the bspline space and the
c       requested error tolerances) and allocates work storage.
c
c       The calling sequence of BACOLR is:
c
c       call bacolr(t0, tout, atol, rtol, npde, kcol, nint, nintmx, x, 
c    &             mflag, rpar, lrp, ipar, lip, cpar, lcp, y, idid)
c
c       which will generate the vector y of bspline coefficients for the
c       approximation of U at t = tout, upon successful completion.
c       Generally, the call to BACOLR will be followed by a 
c       call to VALUES to calculate the solution at a set of points:
c
c       call values(kcol, xsol, nint, x, npde, npts, nderiv, 
c     &             usol, y, work)
c
c       The details of the parameters to VALUES are documented within 
c       the source code for that routine. The input parameters for
c       BACOLR are dealt with in detail below, but a quick summary is:
c
c       [t0, tout] is the time domain of the problem.
c       atol is the absolute error tolerance.
c       rtol is the relative error tolerance.
c       npde is the number of components in the PDE system.
c       kcol, nint define the bspline space.
c       nintmx is the maximum number of subintervals allowed.
c       mflag(1:6) is used to control the operation of BACOLR. 
c       x is the spatial mesh.
c       rpar(lrp) is a floating point work array.
c       ipar(lip) is an integer work array.
c       cpar(lcp) is a complex work array.
c       y is the set of coefficients in the bspline approximation to U at
c       t = tout.
c       idid is an exit status flag provided by BACOLR.
c       The user must check idid to determine what further action needs
c       to be taken.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
c
        double precision        t0
c       On input, t0 < tout is the initial time. On output, t0 is the
c       current time, t0 <= tout. 
c
        double precision        tout
c       tout is the desired final output time. After a successful 
c       return from BACOLR, the time stepping may be resumed by
c       changing tout so that t0 < tout and setting mflag(1) = 1
c       to indicate a continuation of the previous problem.
c
        integer                 npde
c       npde is the number of components in the system of PDEs. npde > 0.
c
        double precision        atol(npde)
c       atol is the absolute error tolerance requested and 
c       is a scalar quantity if mflag(2) = 0.
c
c       If the PDE components vary in importance, then vector error
c       tolerances may be used by setting mflag(2) = 1. In this 
c       case, the dimension of atol must be npde. The user will 
c       define atol(1), atol(2), ..., atol(npde) appropriately. 
c       Note that a change from scalar to vector tolerances (or vice 
c       versa) constitutes a new problem, and BACOLR will have to
c       be reinitialized.
c
        double precision        rtol(npde)
c       rtol is the relative error tolerance request and is a scalar 
c       quantity if mflag(2) = 0.
c
c       If the PDE components vary in importance, then vector error 
c       tolerances may be used by setting mflag(2) = 1. In this 
c       case, the dimension of rtol must be npde. The user will define
c       rtol(1), rtol(2), ..., rtol(npde) appropriately. 
c       Note that a change from scalar to vector tolerances (or vice 
c       versa) constitutes a new problem, and BACOLR will have to
c       be reinitialized.
c
        integer                 kcol
c       kcol is the number of collocation points to be used in each 
c       subinterval. 1 < kcol <= mxkcol.
c
c       The order of the bsplines used will be (kcol+2). 
c  
        integer                 nint
c       at input, nint is the number of subintervals defined by the 
c       spatial mesh x at the initial time t0. at output, nint is 
c       the number of subintervals at tout. nint >= 1.
c
        integer                 nintmx
c       the maximum number of subintervals that the user requires.
c
        double precision        x(nintmx+1)
c       x is the spatial mesh which divides the interval [x_a,x_b] 
c       as: x_a = x(1) < x(2) < x(3) < ... < x(nint+1) = x_b.
c       At input, x(1:nint+1) stores the mesh points at the initial
c       time t0.  at output, x(1:nint+1) stores the mesh points at tout.
c
        integer                 mflag(6)
c       This vector determines the interaction of BACOLR with RADAU5.
c
c       How to set mflag(1):
c
c       On the initial call to BACOLR with a new problem, set 
c       mflag(1) = 0, which indicates that BACOLR and RADAU5 should
c       perform the initialization steps that are required by each code,
c       respectively.
c
c       In order to continue time stepping in the current problem after 
c       a successful return from BACOLR, set mflag(1) = 1,
c       idid = 1, and ensure that t0 < tout.
c
c       How to set mflag(2):
c
c       If scalar absolute and relative error tolerances (atol and rtol)
c       are desired, then set mflag(2) = 0.
c
c       For vector absolute and relative error tolerances, set 
c       mflag(2) = 1, define atol(1), ..., atol(npde), and
c       rtol(1), ..., rtol(npde), as described above, ensuring that 
c       the dimension of each of atol and rtol is at least npde.
c
c       How to set mflag(3):
c
c       It is not used in BACOLR.
c
c       How to set mflag(4):

c       It is not used in BACOLR.
c
c       How to set mflag(5):
c
c       If both boundary conditions are dirichlet, set mflag(5) = 1;
c       else, set mflag(5) = 0.
c
c       How to set mflag(6):
c
c       If the user wants to specify an initial stepsize, set 
c       mflag(6) = 1, and define rpar(2) = the initial stepsize;
c       else, set mflag(6) = 0;
c                                                       
        integer                 lrp
c       lrp is the size of the rpar storage array and must satisfy:
c       lrp >= 74 + 24*npde*npde*nintmx*kcol 
c            + 8*npde*npde*nintmx*kcol*kcol
c            + 29*npde*nintmx*kcol+61*npde+14*kcol
c            + 35*nintmx*kcol
c            + 35*nintmx+21*nintmx*npde+8*nintmx*kcol*kcol
c            + 37*npde*npde+12*npde*npde*nintmx
c
        integer                 lip
c       lip is the size of the ipar integer work array and must satisfy:
c       lip>=100+3*npde*(nintmx*(2*kcol+1)+4)
c
        integer                 lcp
c       lcp is the size of the cpar complex work array and must satisfy:
c       lcp>= npde*(4+(2*kcol+1)*nintmx) 
c           + npde*npde*(8+nintmx*(2*kcol*(kcol+3)+3))
c
c       Work Storage:
        double precision        rpar(lrp)
c       rpar is a floating point work array of size lrp.
c
        integer                 ipar(lip)
c       ipar is an integer work array of size lip.
c
        double complex          cpar(lcp)
c       cpar is a complex work array of size lcp.
c
c       Output:
        double precision        y(npde*(kcol*nintmx+2))
c       On successful return from BACOLR, y(1:npde*(kcol*nint+2)) is 
c       the vector of bspline coefficients at the current time.
c
        integer                 idid
c       idid is the BACOLR exit status flag which is based on the exit
c       status from RADAU5 plus some additional status codes based on 
c       error checking performed by BACOLR on initialization. Positive 
c       values of idid indicate a successful return. Negative values of
c       idid indicate an error which may or may not be fatal.
c       The exact descriptions of idid return values will be discussed below.
c
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
c-----------------------------------------------------------------------

c Constants:
        integer                 nconti 
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points. 
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
        integer                 maxrsh
        parameter              (maxrsh = 20)
c                               maxrsh is the maximum number of
c                               remesh times at one time step,
c                               i.e., ipar(icount) must be less than or
c                               equal to maxrsh
c
        double precision        point1
        parameter              (point1 = 0.1D0)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Local variables:  
c
        integer                 neq1
c                               neq1=npde*ncpts1 is the number of  
c                               bspline coefficients (or DAEs) when  
c                               using radau_{kcol}.
c
        integer                 neq2
c                               neq2=neq1+npde*nint is the number of
c                               bspline coefficients (or DAEs) when
c                               using radau_{kcol+1}.
c
        integer                 neq
c                               neq = neq1 + neq2.
c
        integer                 leniw
c                               leniw = 20 + neq is the length of the
c                               integer work array required by RADAU5.
c                               
        integer                 lenpd1
c                               lenpd1 is the size of the Almost Block
c                               Diagonal (ABD) Jacobian required by
c                               radau_{kcol}.
c                               lenpd1=npde*npde*(2*nconti 
c                                      +kcol*(kcol+nconti)*nint) 
c
        integer                 lenpd2
c                               lenpd2 is the size of the Almost Block
c                               Diagonal (ABD) Jacobian required by
c                               radau_{kcol+1}.
c                               lenpd2=lenpd1+npde*npde*nint
c                                      *(2*kcol+nconti+1)
c
        integer                 lenpd
c                               lenpd = lenpd1 + lenpd2 .
c
        integer                 lenrw
c                               lenrw = 20+12*neq+3*lenpd
c                               is the total size of the floating point
c                               work array required by radau_{kcol}.
c
        integer                 lenin1
c                               lenin1 is the size of the floating
c                               point work array used by INIY and INIYP 
c                               when using radau_{kcol}.
c                               lenin1>=lenpd1+2*neq1+npde*2+2*npde*npde
c
        integer                 lenin2
c                               lenin2 is the size of the floating
c                               point work array used by INIY and INIYP 
c                               when using radau_{kcol+1}.
c                               lenin2>=lenpd2+2*neq2+npde*2+2*npde*npde
c
        integer                 lenri1
c                               lenri1 is the size of the floating
c                               point work array used by REINIT when
c                               using radau_{kcol}.
c 
        integer                 lenri2
c                               lenri2 is the size of the floating
c                               point work array used by REINIT when
c                               using radau_{kcol+1}.
c 
        integer                 lenrj
c                               lenrj is the size of the floating
c                               point work array used by RES and JAC.
c                               lenrj>=4*npde+5*npde*npde.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c                               lenerr>=2*npde*necpts+npde*nint.
c
        integer                 ncpts1
c                               ncpts1=(kcol*nint+nconti) is the total
c                               number of collocation points when using
c                               radau_{kcol}.
c
        integer                 ncpts2
c                               ncpts2=ncpts1+nint is the total number
c                               of collocation points when using
c                               radau_{kcol+1}.
c
        integer                 necpts
c                               necpts=(kcol+3)*nint is the total number
c                               of collocation points used for
c                               error estimate.
c
        integer                 icflag
c                               This is the status flag from the almost
c                               block diagnonal factorization routine,
c                               CRDCMP.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c
        double precision        torign
c                               torign is the initial time, i.e. = t0
c                               at the beginning.
c
c       integer                 isstep
c                               isstep is the number of accepted time
c                               steps since we restart BACOLR in the
c                               case that mflag(4) = 1.
c
        integer                 ninold
c                               ninold is the number of subintervals 
c                               before the current remeshing.
c 
        integer                 ninpre
c                               ninpre is the number of subintervals 
c                               when ipar(icount) = 0 before remeshing.
c 
        integer                 neqpre
c                               neqpre is the number of bspline 
c                               coefficients when ipar(icount)=0 before
c                               remeshing when using radau_{kcol+1}.
c
        integer                 irold
c                               irold is the value of ipar(ixold) before
c                               remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 iiniss
        parameter              (iiniss =  2)
c                               rpar(iiniss) = the initial stepsize when
c                               mflag(6) = 1.
c
        integer                 ierrat
        parameter              (ierrat =  3)
c                               rpar(ierrat) = the value of the largest
c                               component of rpar(ipar(iercom)).
c
        integer                 it0
        parameter              (it0    =  4)
c                               rpar(it0)    = t0 at the last accepted
c                               time step.
c 
        integer                 irpstr
        parameter              (irpstr = 11)
c                               rpar(1:irpstr-1) are reserved to store
c                               floating point scalar quantities.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt1
        parameter              (incpt1 =  4)
c                               ipar(incpt1) = ncpts1.
c
        integer                 ineq1
        parameter              (ineq1  =  5)
c                               ipar(ineq1) = neq1.
c
        integer                 iipstp
        parameter              (iipstp =  6)
c                               ipar(iipstp) = the minimum size of ipar.
c
        integer                 irpstp
        parameter              (irpstp =  7)
c                               ipar(irpstp) = the minimum size of rpar.
c
        integer                 irshin
        parameter              (irshin =  9)
c                               ipar(irshin) is the number of remeshing
c                               times at the initial step.
c
        integer                 isteps
        parameter              (isteps = 10)
c                               ipar(isteps) is the number of time steps
c                               on the current problem.
c
        integer                 irmesh
        parameter              (irmesh = 11)
c                               ipar(irmesh) is the number of remeshing
c                               times after BACOLR starts the initial
c                               step.
c
        integer                 istblc
        parameter              (istblc = 13)
c                               ipar(istblc) is the number of steps
c                               BACOLR has taken before the latest cold
c                               start.
c
        integer                 irshfg
        parameter              (irshfg = 14)
c                               ipar(irshfg) is a flag for redefining
c                               all the pointers.
c                               ipar(irshfg) = 0, the initial step or
c                                                 any step not needing
c                                                 remesh;
c                                            = 1, a step needing remesh.
c 
        integer                 icount
        parameter              (icount = 15)
c                               ipar(icount) is the number of remeshing 
c                               times at the current step.
c
        integer                 istart
        parameter              (istart = 16)
c                               ipar(istart) is a flag to begin the
c                               code.
c                               ipar(istart) = 0, the initial step;
c                                            = 1, not the initial step.
c
        integer                 iradi
        parameter              (iradi  = 61)
c                               ipar(iradi) stores, before remeshing,
c                               the first 20 elements of the integer
c                               point work array in RADAU5.
c
        integer                 iiwork
        parameter              (iiwork = 81)
c                               ipar(iiwork) is the integer work array  
c                               for RADAU5.
c
        integer                 ipivot
        parameter              (ipivot = 101)
c                               ipar(ipivot-1+i), i = 1, neq, contains
c                               the pivoting information from the
c                               factorization of the temporary matrix
c                               for RADAU5.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 ih
        parameter              (ih     = 21) 
c                               rpar(ipar(ih)) stores the mesh step 
c                               size sequence.
c
        integer                 ixcol1
        parameter              (ixcol1 = 22)
c                               rpar(ipar(ixcol1)) stores the 
c                               collocation points when using
c                               radau_{kcol}. 
c
        integer                 ixbs1
        parameter              (ixbs1  = 23)
c                               rpar(ipar(ixbs1)) stores the breakpoint
c                               sequence when using radau_{kcol}.
c
        integer                 iy1
        parameter              (iy1    = 24)
c                               rpar(ipar(iy1)) stores the vector of
c                               solution components to the DAE system
c                               when using radau_{kcol}.
c
        integer                 iabtp1
        parameter              (iabtp1 = 26)
c                               rpar(ipar(iabtp1)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_{kcol}.
c
        integer                 iabbk1
        parameter              (iabbk1 = 27)
c                               rpar(ipar(iabbk1)) stores the nint 
c                               blocks in the middle of the ABD 
c                               collocation matrices when using 
c                               radau_{kcol}.
c
        integer                 iabbt1
        parameter              (iabbt1 = 28)
c                               rpar(ipar(iabbt1)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_{kcol}.
c
        integer                 irwork
        parameter              (irwork = 29)
c                               rpar(ipar(irwork)) stores the floating
c                               point work array for RADAU5. And it
c                               is also used to be a work storage for
c                               the subroutine INIY and INIYP to get
c                               the initial condition.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by RES and JAC.
c
        integer                 ibasi1
        parameter              (ibasi1 = 31)
c                               rpar(ipar(ibasi1)) stores the basis
c                               function values at the collocation
c                               points when using radau_{kcol}.
c                               rpar(ipar(ibasi1)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts1). A(k,j,i) stores
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 iatol
        parameter              (iatol  = 32)
c                               rpar(ipar(iatol)) = atol. 
c
        integer                 irtol
        parameter              (irtol  = 33)
c                               rpar(ipar(irtol)) = rtol. 
c
        integer                 iexcol
        parameter              (iexcol = 34)
c                               rpar(ipar(iexcol)) stores the 
c                               collocation points which are used for
c                               error estimate.
c
        integer                 iewts
        parameter              (iewts  = 35)
c                               rpar(ipar(iewts)) stores the gaussian
c                               weights which are used for error 
c                               estimate.
c
        integer                 iebas1
        parameter              (iebas1 = 36)
c                               rpar(ipar(iebas1)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               radau_{kcol}.
c
        integer                 iebas2
        parameter              (iebas2 = 37)
c                               rpar(ipar(iebas2)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               radau_{kcol+1}.
c
        integer                 iercom
        parameter              (iercom = 38)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 39)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 40)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 ixcol2
        parameter              (ixcol2 = 41)
c                               rpar(ipar(ixcol2)) stores the 
c                               collocation points when using
c                               radau_{kcol+1}.
c
        integer                 ixbs2
        parameter              (ixbs2  = 42)
c                               rpar(ipar(ixbs2)) stores the breakpoint
c                               sequence when using radau_{kcol+1}.
c
        integer                 iy2
        parameter              (iy2    = 43)
c                               rpar(ipar(iy2)) stores the vector of
c                               solution components to the DAE system
c                               when using radau_{kcol+1}.
c
        integer                 iabtp2
        parameter              (iabtp2 = 45)
c                               rpar(ipar(iabtp2)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_{kcol+1}.
c
        integer                 iabbk2
        parameter              (iabbk2 = 46)
c                               rpar(ipar(iabbk2)) stores the nint 
c                               blocks in the middle of the ABD 
c                               collocation matrices when using
c                               radau_{kcol+1}.
c
        integer                 iabbt2
        parameter              (iabbt2 = 47)
c                               rpar(ipar(iabbt2)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_{kcol+1}.
c
        integer                 ibasi2
        parameter              (ibasi2 = 48)
c                               rpar(ipar(ibasi2)) stores the basis
c                               function values at the collocation
c                               points when using radau_{kcol+1}.
c                               rpar(ipar(ibasi2)) contains
c                               a three dimensional array A of size
c                               (kcol+1+nconti,3,ncpts2). A(k,j,i) 
c                               stores the values of the (j-1)st 
c                               derivative (j=1,2,3) of the k-th 
c                               non-zero basis function (k=1,...,
c                               kcol+1+nconti) at the i-th collocation 
c                               point.
c
        integer                 ixold
        parameter              (ixold  = 51)
c                               rpar(ipar(ixold)) stores the mesh point
c                               sequence when ipar(icount) = 0 before
c                               remeshing.
c
        integer                 iypre
        parameter              (iypre  = 53)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy2)) at the previous 6 steps.
c                               It is required for a hot restart after
c                               remeshing.
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               colpnt
c                               RADAU5
c                               iniy
c                               meshsq
c                               reinit
c                               remesh
        external                radfcn
        external                radjac
        external                radmas
        external                solout
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Oct 18, 2006.
c
c-----------------------------------------------------------------------

c     Check validity of the mflag vector.
      do 1 i = 1, 6
         if ((mflag(i) .lt. 0) .or. (mflag(i) .gt. 1)) goto 710
   1  continue

      ipar(irshfg) = 0

c     Check for continuation of a previous problem.
      if (mflag(1) .eq. 1) then
         ipar(istart) = 1

         neq1   = ipar(ineq1)
         neq2   = neq1 + npde*nint
         neq    = neq1 + neq2
         lenpd1 = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint)
         lenpd  = lenpd1+lenpd1+npde*npde*nint*(kcol+kcol+nconti+1)
         leniw  = 20 + neq*2
         lenrw  = 20 + 12*neq + 3*lenpd
         
         goto 200
      else

c        Check if the user specifies an initial stepsize
         if (mflag(6) .eq. 1) then
            if (((tout-t0)*rpar(iiniss)) .lt. zero) goto 720
            if (rpar(iiniss) .eq. zero) goto 725
         endif      

         do 10 i = 1, iradi
            ipar(i) = 0
   10    continue
         do 20 i = irpstr, lrp
            rpar(i) = zero
   20    continue
         torign = t0
      endif

c-----------------------------------------------------------------------
c     On the initial call or after remeshing, check for valid input and
c     initialize the workspace.
c-----------------------------------------------------------------------

  100 continue

c     Check validity of npde, kcol, and nint.
      if (npde .le. 0) goto 730
      if ((kcol .le. 1) .or. (kcol .gt. mxkcol)) goto 740
      if ((nint .le. 0) .or. (nint .gt. nintmx)) goto 750

c     Check for a monotone mesh.
      do 110 i = 1, nint
         if (x(i) .ge. x(i+1)) goto 760
  110 continue

c-----------------------------------------------------------------------
c     Calculate the extra storage requirements of res and jac.
      lenrj = (4 + 5 * npde) * npde

c     Calculate the number of collocation points when using
c     radau_{kcol}.
      ncpts1 = nint * kcol + nconti

c     Calculate the number of DAEs when using radau_{kcol}.
      neq1 = npde * ncpts1

c     Size of the ABD iteration matrix when using radau_{kcol}.
      lenpd1 = npde*npde*(nconti+nconti+kcol*(kcol+nconti)*nint) 

c     Calculate the extra storage requirements of iniy when
c     using radau_{kcol}.
      lenin1 = lenpd1 + 2 * neq1 + 2 * npde * (1 + npde)

c-----------------------------------------------------------------------
c     Calculate the number of collocation points when using
c     radau_{kcol+1}.
      ncpts2 = ncpts1 + nint

c     Calculate the number of DAEs when using radau_{kcol+1}.
      neq2 = neq1 + npde * nint

c     Size of the ABD iteration matrix when using radau_{kcol+1}.
      lenpd2 = lenpd1 + npde * npde * nint * (kcol + kcol + nconti + 1) 

c     Calculate the extra storage requirements of iniy when
c     using radau_{kcol+1}.
      lenin2 = lenpd2 + 2 * neq2 + 2 * npde * (1 + npde)

c-----------------------------------------------------------------------
c     Calculate the total number of variables given to RADAU5.
      neq = neq1 + neq2

c     Size of the total iteration matrix in RADAU5.
      lenpd = lenpd1 + lenpd2

c     Total size of the RADAU5 floating point work array.
      lenrw = 20 + 12 * neq + 3 * lenpd

c     Total size of the RADAU5 integer work array.
      leniw = 20 + neq * 2

c-----------------------------------------------------------------------
c     Calculate the number of collocation point used for error estimate.
      necpts = (kcol + 3) * nint

c     Calculate the extra storage requirements of errest.
      lenerr = (2 * necpts + nint) * npde

c-----------------------------------------------------------------------
c     Save the input parameters in ipar, the integer communication
c     storage array.
      ipar(inpde)  = npde
      ipar(ikcol)  = kcol
      ipar(inint)  = nint
      ipar(incpt1) = ncpts1
      ipar(ineq1)  = neq1

c-----------------------------------------------------------------------
c     Calculate the offsets into rpar, the floating point storage array.
c-----------------------------------------------------------------------
      ipar(iatol)  = irpstr
      ipar(irtol)  = ipar(iatol)  + npde

      ipar(ih)     = ipar(irtol)  + npde

      ipar(iy1)    = ipar(ih)     + nint
      ipar(iy2)    = ipar(iy1)    + neq1

      ipar(ixcol1) = ipar(iy2)    + neq2
      ipar(ixbs1)  = ipar(ixcol1) + ncpts1
      ipar(iabtp1) = ipar(ixbs1)  + ncpts1 + kcol + nconti
      ipar(iabbk1) = ipar(iabtp1) + npde * npde * nconti
      ipar(iabbt1) = ipar(iabbk1) + npde * npde * nint * kcol
     &                              * (kcol + nconti)
      ipar(ibasi1) = ipar(iabbt1) + npde * npde * nconti

      ipar(irwork) = ipar(ibasi1) + (kcol + nconti) * 3 * ncpts1
      ipar(iwkrj)  = ipar(irwork) + lenrw
 
      ipar(iexcol) = ipar(iwkrj)  + lenrj
      ipar(iewts)  = ipar(iexcol) + necpts
      ipar(ierint) = ipar(iewts)  + necpts 
      ipar(iercom) = ipar(ierint) + nint
      ipar(iebas1) = ipar(iercom) + npde
      ipar(iebas2) = ipar(iebas1) + (kcol + nconti) * necpts
      ipar(iework) = ipar(iebas2) + (kcol + 1 + nconti) * necpts
 
      ipar(ixcol2) = ipar(iework) + lenerr
      ipar(ixbs2)  = ipar(ixcol2) + ncpts2
      ipar(iabtp2) = ipar(ixbs2)  + ncpts2 + kcol + 1 + nconti 
      ipar(iabbk2) = ipar(iabtp2) + npde * npde * nconti
      ipar(iabbt2) = ipar(iabbk2) + npde * npde * nint * (kcol+1)
     &                              * (kcol + 1 + nconti)
      ipar(ibasi2) = ipar(iabbt2) + npde * npde * nconti

      ipar(ixold)  = ipar(ibasi2) + (kcol + 1 + nconti) * 3 * ncpts2

      ipar(iypre)  = ipar(ixold)  + nintmx + 1

c     The offset is different between the initial call and remeshing.
      if ((ipar(irshfg) .ne. 0) .and. (ipar(istart) .eq. 1)) then
         ipar(irpstp) = ipar(iypre) + neqpre - 1
      else
         ipar(irpstp) = ipar(iypre) + neq2 - 1
      endif

c     Check for a sufficiently large rpar floating point work array.
      if (lrp .lt. ipar(irpstp)) goto 770

c     Calculate the offsets into the integer storage array.
      ipar(iipstp) = ipivot + 3 * neq - 1

c     Check for a sufficiently large ipar integer work array.
      if (lip .lt. ipar(iipstp)) goto 780

c     Check whether it is initial call or for remeshing.
      if ((ipar(irshfg) .ne. 0) .and. (ipar(istart) .ne. 0)) goto 300

c-----------------------------------------------------------------------
c     Save the atol and rtol in rpar, the real communication
c     storage array.
      if (mflag(2) .eq. 0) then
         rpar(ipar(iatol)) = atol(1)
         rpar(ipar(irtol)) = rtol(1)
      else
         do 120 i = 1, npde
            rpar(ipar(iatol)-1+i) = atol(i)
            rpar(ipar(irtol)-1+i) = rtol(i)
  120    continue
      endif

c-----------------------------------------------------------------------
c     Perform initializations for using radau_{kcol}.
c-----------------------------------------------------------------------
c     Check whether it is initial call or for remeshing.
      if (ipar(irshfg) .eq. 0) then

c        Set the initial stepsize if applicable.
         if (mflag(6) .eq. 0) then
            if (mflag(2) .eq. 0) then
               rpar(iiniss) = max(atol(1),rtol(1))
            else
               rpar(iiniss) = zero
               do 130 i = 1, npde
                  rpar(iiniss) = max(rpar(iiniss), atol(i))
  130          continue
               do 140 i = 1, npde
                  rpar(iiniss) = max(rpar(iiniss), rtol(i))
  140          continue
            endif
            if (((tout-t0)*rpar(iiniss)) .lt. zero) 
     &         rpar(iiniss) = - rpar(iiniss)
         endif
      else
         ipar(irshin) = ipar(irshin) + 1
      endif

      call meshsq(kcol, nint, x, rpar(ipar(irwork)), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

      call colpnt(kcol, nint, ncpts1, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol1)),
     &            rpar(ipar(ixbs1)))

      icflag = 0

      call iniy(t0, npde, kcol, nint, neq1, ncpts1, mflag(5), 
     &          rpar(ipar(ixcol1)), rpar(ipar(ixbs1)),
     &          rpar(ipar(iabbk1)), rpar(ipar(ibasi1)), rpar(ipar(iy1)),
     &          ipar(ipivot), rpar(ipar(irwork)), lenin1, icflag)
 
      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

      ipar(irshfg) = 0

c-----------------------------------------------------------------------
c     Perform initializations for using radau_{kcol+1}.
c-----------------------------------------------------------------------

      call colpnt(kcol+1, nint, ncpts2, x, rpar(ipar(ih)),
     &            rpar(ipar(irwork)), rpar(ipar(ixcol2)),
     &            rpar(ipar(ixbs2)))

      icflag = 0

      call iniy(t0, npde, kcol+1, nint, neq2, ncpts2, mflag(5),
     &          rpar(ipar(ixcol2)), rpar(ipar(ixbs2)),
     &          rpar(ipar(iabbk2)), rpar(ipar(ibasi2)),
     &          rpar(ipar(iy2)), ipar(ipivot),
     &          rpar(ipar(irwork)), lenin2, icflag)
 
      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

c     Copy rpar(ipar(iy2)) to rpar(ipar(iypre)).
      call dcopy(neq2, rpar(ipar(iy2)), 1, rpar(ipar(iypre)), 1)

      do 150 i = 1, 20
         ipar(iiwork-1+i) = 0
  150 continue

      do 160 i = 1, 20
         rpar(ipar(irwork)-1+i) = zero
  160 continue

      goto 400

c-----------------------------------------------------------------------
c     This is not the first call for the problem, and integration is to
c     continue.
c-----------------------------------------------------------------------
  200 continue

c     Examine idid to determine if RADAU5 can be called again.
      if (idid .ne. 1) goto 790

      goto 400
 
c-----------------------------------------------------------------------
c     Initialization after remeshing.
c-----------------------------------------------------------------------
  300 continue

      ipar(irmesh) = ipar(irmesh) + 1

      do 310 i = 1, 20
         ipar(iiwork-1+i) = 0
  310 continue
      do 320 i = 1, 20
         rpar(ipar(irwork)-1+i) = zero
  320 continue

      rpar(ipar(irwork)-1+3) = point1

      if (nint .lt. ninold) then
         call dcopy(nintmx+1+20+neqpre, rpar(irold),
     &              1, rpar(ipar(ixold)), 1)
      else
         if (nint .gt. ninold) then
            call dcopy(nintmx+1+20+neqpre, rpar(irold),
     &                 -1, rpar(ipar(ixold)), -1)
         endif
      endif

      lenri1 = lenpd1 + kcol + 1 + nconti + (kcol + 1) * (ninpre + 1)
     &         + 2 * nconti
      lenri2 = lenpd2 + kcol + 1 + nconti + (kcol + 1) * (ninpre + 1)
     &         + 2 * nconti

      call meshsq(kcol, nint, x, rpar(ipar(irwork)+20), rpar(ipar(ih)),
     &            rpar(ipar(iexcol)), rpar(ipar(iewts)))

      call reinit(npde, kcol, kcol+1, nint, ninpre, ncpts1, neq1,
     &            neqpre, x, rpar(ipar(ixold)),
     &            rpar(ipar(iypre)), rpar(ipar(irwork)+20), lenri1,
     &            ipar(ipivot), rpar(ipar(ih)), rpar(ipar(ixbs1)),
     &            rpar(ipar(ixcol1)), rpar(ipar(ibasi1)),
     &            rpar(ipar(iy1)), rpar(ipar(iabbk1)), icflag)

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

      call reinit(npde, kcol+1, kcol+1, nint, ninpre, ncpts2, neq2,
     &            neqpre, x, rpar(ipar(ixold)),
     &            rpar(ipar(iypre)), rpar(ipar(irwork)+20), lenri2,
     &            ipar(ipivot), rpar(ipar(ih)), rpar(ipar(ixbs2)),
     &            rpar(ipar(ixcol2)), rpar(ipar(ibasi2)),
     &            rpar(ipar(iy2)), rpar(ipar(iabbk2)), icflag)

      if (icflag .ne. 0) then
         idid = -66
         goto 600
      endif

c-----------------------------------------------------------------------
c     Time integration loop for RADAU5.
c-----------------------------------------------------------------------

  400 continue

      idid = 0
      call radau5(neq, radfcn, t0, rpar(ipar(iy1)), tout, rpar(iiniss), 
     &            rtol, atol, mflag(2), radjac, radmas, solout,
     &            rpar(ipar(irwork)), lenrw, ipar(iiwork), leniw, rpar,
     &            ipar, cpar, idid)

c-----------------------------------------------------------------------
c     Check for a successful time step and decide whether to continue
c     integration or to perform a remeshing.
c-----------------------------------------------------------------------

      if (idid .le. 0) goto 600

      if (idid .eq. 2) then

c        The current step is rejected.
         if (ipar(icount) .eq. maxrsh) goto 610

c        For the first remeshing at the current step, save nintpre and
c        neqpre at the last successful step.
         if (ipar(icount) .eq. 0) then
            ninpre = nint
            neqpre = neq2
         endif

         do 410 i = 14, 20
            ipar(iradi-1+i) = ipar(iiwork-1+i)
  410    continue
         ninold = nint
         irold = ipar(ixold)

c        Update xold.
         if (ipar(icount) .eq. 0) then
            do 420 i = 1, ninpre + 1
               rpar(ipar(ixold)-1+i) = x(i)
  420       continue
         endif

         call remesh(ipar(istart), ipar(icount), nintmx,  
     &               ninold, rpar(ierrat), rpar(ipar(ierint)), 
     &               ipar(irshfg), nint, kcol, x, rpar(ipar(iework)))

         if (ipar(istart) .eq. 1) then

c           This is not the initial step.
            t0 = rpar(it0)

            ipar(istblc) = ipar(istblc) + ipar(iradi-1+17) - 1

         else

c           This is the initial step.
            t0 = torign

         endif

         goto 100

      else

c        The current step is the last step.
         goto 500

      endif

c-----------------------------------------------------------------------
c     Successful return section.
c-----------------------------------------------------------------------
  500 continue

c     Retrieve the value of mflag(1).
      mflag(1) = 1

c     Retrieve the output vector y from the rpar communication array.
      do 510 i = 1, neq1
         y(i) = rpar(ipar(iy1)-1+i)
  510 continue

c     Retrieve information on the time stepping from the ipar array.
      ipar(isteps) = ipar(istblc) + ipar(iiwork-1+17)

      return

c-----------------------------------------------------------------------
c     Unsuccessful return section.
c-----------------------------------------------------------------------
  600 continue
      write(6,9999) 'ERROR: BACOLR runtime error in time stepping.'
      write(6,9999) '       An error code and message should have'
      write(6,9999) '       been issued by RADAU5.'
      return
  610 continue
      if (ipar(istart) .eq. 1) then
         write(6,9998) 'ERROR: BACOLR has remeshed ', maxrsh,' times at'
     &                 , ' t0 =', rpar(it0)
      else
         write(6,9998) 'ERROR: BACOLR has remeshed ', maxrsh,' times at'
     &                 , ' t0 =', torign
      endif
      return

c-----------------------------------------------------------------------
c     The following section is the return point for invalid input.
c-----------------------------------------------------------------------

  710 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require:  0 <= mflag(i) <= 1, i = 1, 2, ..., 6.' 
      idid = -51
      return
c 715 continue
c     write(6,9999) 'ERROR: BACOLR input violation.'
c     write(6,9999) 'Require:  if mflag(4) = 1, ipar(8) must be set to' 
c     write(6,9999) 'be a positive integer.' 
c     idid = -52
c     return
  720 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require:  if mflag(6) = 1, tout must be in front'
      write(6,9999) 'of t0.' 
      idid = -53
      return
  725 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require:  if mflag(6) = 1, rpar(2) must be the'
      write(6,9999) 'initial stepsize, thus nonzero.' 
      idid = -54
      return
  730 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: npde > 0.'      
      idid = -55
      return
  740 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: 1 < kcol <=', mxkcol, '.'
      idid = -56
      return
  750 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: 0 < nint <=', nintmx, '.'
      idid = -57
      return
  760 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: x(1) < x(2) < ... < x(nint+1).'
      idid = -58
      return
  770 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: lrp >= ', ipar(irpstp), '.'
      idid = -59
      return
  780 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'Require: lip >= ', ipar(iipstp), '.'
      idid = -60
      return
  790 continue
      write(6,9999) 'ERROR: BACOLR input violation.'
      write(6,9999) 'IDID .ne. 1, on a continuation call of BACOLR'
      write(6,9999) 'If IDID > 1, set idid = 1 and tout (t0 < tout)'
      write(6,9999) 'If IDID < -1, the code cannot be continued due to'
      write(6,9999) '              a previous error.'
      idid = -61
      return

c-----------------------------------------------------------------------
 9998 format(a,i4,a,a,e12.5)
 9999 format(a,i8,a,i4,a,i4,a,i4,a,i4)
c-----------------------------------------------------------------------
      end
      SUBROUTINE BSPLVD ( XT, K, X, ILEFT, VNIKX, NDERIV )
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C CALCULATES THE VALUE AND THE FIRST NDERIV-1 DERIVATIVES OF ALL
C B-SPLINES WHICH DO NOT VANISH AT X.  THE ROUTINE FILLS THE TWO-
C DIMENSIONAL ARRAY VNIKX(J,IDERIV), J=IDERIV, ... ,K WITH NONZERO
C VALUES OF B-SPLINES OF ORDER K+1-IDERIV, IDERIV=NDERIV, ... ,1, BY
C REPEATED CALLS TO BSPLVN.
C
C LAST MODIFIED BY RONG WANG, DEC 13, 2006.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  BSPLVN
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   DBLE,MAX
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER K,NDERIV,ILEFT
      DOUBLE PRECISION X
      DOUBLE PRECISION XT(*),VNIKX(K,NDERIV)
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER KO,IDERIV,IDERVM,KMD,JM1,IPKMD,JLOW
      DOUBLE PRECISION A(20,20)
      DOUBLE PRECISION FKMD,DIFF,V
C-----------------------------------------------------------------------
C CONSTANT
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
C-----------------------------------------------------------------------
C LOOP INDICES
      INTEGER I,J,M,L
C-----------------------------------------------------------------------
      KO = K + 1 - NDERIV
      CALL BSPLVN(XT,KO,1,X,ILEFT,VNIKX(NDERIV,NDERIV))
      IF (NDERIV .LE. 1) GO TO 130
      IDERIV = NDERIV
      DO 20 I=2,NDERIV
        IDERVM = IDERIV-1
        DO 10 J=IDERIV,K
          VNIKX(J-1,IDERVM) = VNIKX(J,IDERIV)
   10   CONTINUE
        IDERIV = IDERVM
        CALL BSPLVN(XT,0,2,X,ILEFT,VNIKX(IDERIV,IDERIV))
   20 CONTINUE
      DO 40 I=1,K
        DO 30 J=1,K
          A(I,J) = ZERO
   30   CONTINUE
        A(I,I) = ONE
   40 CONTINUE
      KMD = K
      DO 120 M=2,NDERIV
        KMD = KMD - 1
        FKMD =  DBLE(KMD)
        I = ILEFT
        J = K
   50   CONTINUE
        JM1 = J-1
        IPKMD = I + KMD
        DIFF = XT(IPKMD) -XT(I)
        IF (JM1 .NE. 0) THEN
          IF (DIFF .NE. ZERO) THEN
            DO 60 L=1,J
              A(L,J) = (A(L,J) - A(L,J-1))/DIFF*FKMD
   60       CONTINUE
          ENDIF
          J = JM1
          I = I - 1
          GO TO 50
        ENDIF
        IF (DIFF .NE. ZERO) THEN
          A(1,1) = A(1,1)/DIFF*FKMD
        ENDIF
        DO 110 I=1,K
          V = ZERO
          JLOW = MAX(I,M)
          DO 100 J=JLOW,K
            V = A(I,J)*VNIKX(J,M) + V
  100     CONTINUE
          VNIKX(I,M) = V
  110   CONTINUE
  120 CONTINUE
  130 RETURN
      END
      SUBROUTINE BSPLVN ( XT, JHIGH, INDEX, X, ILEFT, VNIKX )
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT THE
C POINT X OF ORDER MAX(JHIGH,(J+1)(INDEX-1)) FOR THE BREAKPOINT SEQ-
C UENCE XT.  ASSUMING THAT XT(ILEFT) .LE. X .LE. XT(ILEFT+1), THE ROUT-
C INE RETURNS THE B-SPLINE VALUES IN THE ONE DIMENSIONAL ARRAY VNIKX.
C
C LAST MODIFIED BY RONG WANG, JAN 8, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  NONE
C USER ROUTINES CALLED..     NONE
C CALLED BY..                BSPLVD
C FORTRAN FUNCTIONS USED..   NONE
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      DOUBLE PRECISION XT(*),X,VNIKX(*)
      INTEGER JHIGH,INDEX,ILEFT
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER IPJ,IMJP1,JP1,JP1ML
      DOUBLE PRECISION VMPREV,VM
C-----------------------------------------------------------------------
C CONSTANT
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
      DOUBLE PRECISION DELTAM(20),DELTAP(20)
      INTEGER J
C-----------------------------------------------------------------------
C LOOP INDICE
      INTEGER L
C-----------------------------------------------------------------------
      DATA J/1/,DELTAM/20*0.D+0/,DELTAP/20*0.D+0/
      
      IF(INDEX.EQ.1) THEN
        J = 1
        VNIKX(1) = ONE
        IF (J .GE. JHIGH) GO TO 40
      ENDIF
   20 CONTINUE
      IPJ = ILEFT+J
      DELTAP(J) = XT(IPJ) - X
      IMJP1 = ILEFT-J+1
      DELTAM(J) = X - XT(IMJP1)
      VMPREV = ZERO
      JP1 = J+1
      DO 30 L=1,J
        JP1ML = JP1-L
        VM = VNIKX(L)/(DELTAP(L) + DELTAM(JP1ML))
        VNIKX(L) = VM*DELTAP(L) + VMPREV
        VMPREV = VM*DELTAM(JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      J = JP1
      IF (J .LT. JHIGH) GO TO 20
   40 RETURN
      END
      subroutine calfcn(npde, kcol, nint, ncpts, neq, xcol, fbasis,
     &                  t, y, work, fr)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radfcn. It provides a lower-level
c       interface to calculate the right side of the DAEs.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 21, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bsplines
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               T is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
c       Work storage:
        double precision        work(4*npde+2*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+2*npde*npde.
c
c       Output:
        double precision        fr(neq)
c                               fr is the vector at the right side of
c                               the DAEs.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
c
c       Indices:
        integer                 ii
        integer                 jj
        integer                 mm
        integer                 kk
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               f
c                               eval
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations.
 
      do 20 i = 1, nint
 
c        ii is the value of ileft for the current collocation point.
         ii = kcol + nconti + (i - 1) * kcol
 
         do 10 j = 1, kcol
 
c           jj is the pointer of collocation point.
            jj = (i - 1) * kcol + j + 1
 
c           mm is the pointer of fr.
            mm = (jj - 1) * npde + 1
 
c           kk is the pointer of the basis function values at
c           the current collocation point.
            kk =(jj-1)*(kcol+nconti)*3+1
 
c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ii,jj,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(kk),y)
 
c           Evaluate the function defining the PDE at the current
c           collocation point, storing the result in fr.
            call f(t, xcol(jj), work(iu), work(iux),
     &              work(iuxx), fr(mm), npde)
 
   10    continue
   20 continue

c-----------------------------------------------------------------------   
c     Calculate (fr(i), i=1, npde), which depend on the left
c     boundary point.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call bndxa(t, work(iu), work(iux), fr(1), npde)
 
c-----------------------------------------------------------------------
c     Calculate (fr(i), i=neq-npde+1, neq), which depend on the right
c     boundary point.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call bndxb(t, work(iu), work(iux), fr(neq-npde+1), npde)

      return
      end
      subroutine caljac(npde, kcol, nint, ncpts, neq, xcol, fbasis,
     &                  abdtop, abdbot, t, y, work, dfdy)

c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radjac. It provides a lower-level
c       interface to generate the Jacobian matrix at the right side
c       of the DAEs.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 22, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcolis the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*ncpts is the number of bspline
c                               coefficients (or DAEs).
c
        double precision        xcol(ncpts)
c                               xcol stores the collocation
c                               points when using kcol collocation
c                               points at each subinterval.
c
        double precision        abdtop(npde*npde*nconti)
c                               abdtop stores the top block of the ABD
c                               matrices.
c
        double precision        abdbot(npde*npde*nconti)
c                               abdbot stores the bottom block of the
c                               ABD matrices.
c
        double precision        fbasis((kcol+nconti)*3*ncpts)
c                               fbasis stores the basis function values
c                               at the collocation points. It acts like
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
c       Work storage:
        double precision        work(4*npde+5*npde*npde)
c                               work is a floating point work array
c                               of size 4*npde+5*npde*npde.
c
c       Output:
        double precision        dfdy(npde*npde*(2*nconti
     *                               +nint*kcol*(kcol+nconti)))
c                               dfdy is the ABD Jacobian matrix at the
c                               right side of the DAEs.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ifytop
c                               ifytop is the pointer into dfdy where 
c                               the top block of the ABD Jacobian is
c                               stored.
c
        integer                 ifyblk
c                               ifyblk is the pointer into dfdy where 
c                               the nint blocks in the middle of the ABD
c                               Jacobian are stored.
c
        integer                 ifybot
c                               ifybot is the pointer into dfdy where 
c                               the bottom block of the ABD Jacobian is
c                               stored.
c
        integer                 nsiztb
c                               nsiztb is the size of the top block
c                               as same as the bottom block of the ABD
c                               Jacobian.
c
        integer                 nsizbk
c                               nsizbk is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 n
c
        integer                 ii
        integer                 ij
        integer                 jj
        integer                 kk
        integer                 nn
        integer                 mm
        integer                 jk
        integer                 jk2
        integer                 jk3
        integer                 mn
        integer                 mn2
        integer                 mn3
c
c-----------------------------------------------------------------------
c       Pointers into the floating point work array work:
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idfdu
c                               work(idfdu) stores the Jacobian of f
c                               with respect to u.
c
        integer                 idfdux
c                               work(idfdux) stores the Jacobian of f
c                               with respect to u_x.
c
        integer                 idfuxx
c                               work(idfuxx) stores the Jacobian of f
c                               with respect to u_xx.
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde), That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th
c                               component of the vector b with respect
c                               to time t.
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               derivf
c                               difbxa
c                               difbxb
c                               eval
c
c-----------------------------------------------------------------------

c     Set pointers into the temporary floating point work array.
      iu     = 1
      iux    = iu     + npde
      iuxx   = iux    + npde
      idfdu  = iuxx   + npde
      idfdux = idfdu  + npde * npde
      idfuxx = idfdux + npde * npde
      idbdu  = idfuxx + npde * npde
      idbdux = idbdu  + npde * npde
      idbdt  = idbdux + npde * npde
 
c     Set the indices into dfdy which define the ABD Jacobian.
      ifytop = 1
      ifyblk = ifytop + nconti * npde * npde
      ifybot = ifyblk + nint * npde * npde * kcol * (kcol + nconti)
 
c-----------------------------------------------------------------------
c     Calculate the size of top (or bottom) block and the size of a
c     subblock in the middle.
      nsiztb = npde * npde * nconti
      nsizbk = npde * npde * kcol * (kcol + nconti)
 
c     Initialize fytop, fyblk and fybot to zero.
      do 10 i = 1, nsiztb
         dfdy(ifytop-1+i) = zero
         dfdy(ifybot-1+i) = zero
   10 continue
      do 20 i = 1, nint * nsizbk
         dfdy(ifyblk-1+i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Loop over the nint blocks of collocation equations and
c     caluculate the portion of dfdy which depends on them.
 
      do 70 i = 1, nint
 
c        ii+1 is the pointer to the first element at the i-th subblock
c        of the jacobian matrix, i = 1, nint.
         ii = ifyblk - 1 + (i - 1) * nsizbk
 
c        ij is the value of ileft for the current collocation point.
         ij = kcol + nconti + (i - 1) * kcol
 
         do 60 j = 1, kcol
 
c           jj+1 is the pointer to the first element corresponding to
c           the j-th collocation point in the i-th interval.
            jj = ii + (j - 1) * npde
 
c           mm is the index of the current collocation point.
            mm = (i - 1) * kcol + j + 1
 
c           Generate the approximate solution and its spatial
c           derivatives at the current collocation point.
            call eval(npde,kcol,ij,mm,ncpts,work(iu),work(iux),
     &                work(iuxx),fbasis(1+(mm-1)*(kcol+nconti)*3),y)
 
c           Generate dfdu, dfdux, and dfdux at the current
c           collocation point (the j-th point of the i-th
c           subinterval).
            call derivf(t, xcol(1+(i-1)*kcol+j), work(iu),
     &                  work(iux), work(iuxx), work(idfdu),
     &                  work(idfdux), work(idfuxx), npde)
 
            do 50 k = 1, kcol + nconti
 
c              kk+1 is the pointer to the first element of a npde by
c              npde submatrix, which is corresponding to the j-th
c              collocation point in the i-th interval, and the k-th
c              nonzero basis function.
               kk = jj + (k-1) * npde * npde * kcol
 
c              jk is the pointer to the k-th nonzero function at the
c              mm-th collocation point in the basis function,
c              fbasis(1).
               jk = (mm - 1) * (kcol + nconti) * 3 + k
 
c              jk2 is the pointer to the first derivative for the
c              above basis function.
               jk2 = jk + kcol + nconti
 
c              jk3 is the pointer to the second derivative for the
c              above basis function.
               jk3 = jk2 + kcol + nconti
 
               do 40 m = 1, npde
                  do 30 n = 1, npde
 
c                    nn is the pointer to the (n, m) element of the
c                    npde by npde submatrix.
                     nn = kk + (m-1)*npde*kcol + n
 
c                    mn is the pointer to the (n, m) element of dfdu.
                     mn = idfdu - 1 + (m - 1) * npde + n
 
c                    mn2 is the pointer to the (n, m) element of dfdux.
                     mn2 = mn + npde * npde
 
c                    mn3 is the pointer to the (n, m) element of dfduxx.
                     mn3 = mn2 + npde * npde
 
c                    now set up the value in pd at the place nn.
                     dfdy(nn) = work(mn) * fbasis(jk)
     &                          + work(mn2) * fbasis(jk2)
     &                          + work(mn3) * fbasis(jk3)
 
   30             continue
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Update the values at the left boundary.
      call eval(npde, kcol, kcol+2, 1, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1), y)
      call difbxa(t, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)
 
c     Update the top block of the collocation matrix dG/dY'.
      do 90 j = 1, npde
         do 80 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdtop(jj) =
     &            fbasis(2+kcol+nconti) * work(idbdux-1+mm)
            abdtop(ii) =
     &            work(idbdu-1+mm) - abdtop(jj)
   80    continue
   90 continue
 
c-----------------------------------------------------------------------
c     Update the values at the right boundary.
      call eval(npde, kcol, ncpts, ncpts, ncpts, work(iu), work(iux),
     &          work(iuxx), fbasis(1+(ncpts-1)*(kcol+nconti)*3), y)
      call difbxb(t, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)
 
c     Update the bottom block of the collocation matrix.
      do 110 j = 1, npde
         do 100 i = 1, npde
            ii = (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            abdbot(ii) =
     &            fbasis(1+kcol+kcol+nconti+(ncpts-1)*(kcol+nconti)*3)
     &            * work(idbdux-1+mm)
            abdbot(jj) =
     &            work(idbdu-1+mm) - abdbot(ii)
  100    continue
  110 continue
 
c-----------------------------------------------------------------------
c     Copy abdtop and abdbot to the corresponding parts of dfdy.

      call dcopy(nsiztb, abdtop, 1, dfdy(ifytop), 1)
      call dcopy(nsiztb, abdbot, 1, dfdy(ifybot), 1)

c-----------------------------------------------------------------------
      return
      end
      subroutine calmas(npde, kcol, nint, abdblk, am)
c-----------------------------------------------------------------------
c Purpose:
c       This subroutine is called by radmas. It provides a lower-level
c       interface to generate the mass-matrix.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 20, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcolis the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        double precision        abdblk(npde*npde*nint*kcol
     *                                 *(kcol+nconti))
c                               abdblk stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using kcol
c                               collocation points at each subinterval.
c
c       Output:
        double precision        am(npde*npde*(2*nconti
     *                             +nint*kcol*(kcol+nconti)))
c                               am is the ABD mass-matrix.
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 iamtop
c                               iamtop is the pointer into am where the
c                               top block of the ABD matrix is stored.
c
        integer                 iamblk
c                               iamblk is the pointer into am where the
c                               nint blocks in the middle of the ABD
c                               matrix are stored.
c
        integer                 iambot
c                               iambot is the pointer into am where the
c                               bottom block of the ABD matrix is
c                               stored.
c
        integer                 nels
c                               nels is the size of a subblock in
c                               the middle of ABD Jacobian.
c
c-----------------------------------------------------------------------
c       Loop indices:
        integer                 i
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)   

c     Set the pointers into the mass-matrix.  
      iamtop = 1
      iamblk = iamtop + npde*npde*nconti 
      iambot = iamblk + nint*nels

c     Initialize the top and bottom of the mass-matrix to zero. 
      do 10 i = 1, npde*npde*nconti
         am(i) = zero
   10 continue
      do 20 i = iambot, iambot+npde*npde*nconti-1
         am(i) = zero
   20 continue

c-----------------------------------------------------------------------    
c     The nint blocks at the middle of the mass-matrix are set up. 
      do 30 i = 1, nint*nels
         am(i+npde*npde*nconti) = abdblk(i)
   30 continue

      return
      end
        SUBROUTINE CCRCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C C R C M PDECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - COMPLEX(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - COMPLEX(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - COMPLEX(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C***************************************************************
C
c-----------------------------------------------------------------------
        integer n,nrwtop,novrlp,nrwblk,nclblk,nbloks,nrwbot,iflag,
     &          nrwtp1,nrowel,nrwel1,nvrlp0,i,iplus1,ipvt,j,l,incr,k,
     &          kplus1,jplus1,jminn,loop,incrj,iplusn,incrn,irwblk,
     &          ipvblk,jrwblk
c-----------------------------------------------------------------------
        DOUBLE COMPLEX TOPBLK,ARRAY,BOTBLK
        DOUBLE COMPLEX ROWPIV,ROWMLT,COLPIV,SWAP,COLMLT
        REAL ROWMAX,COLMAX,TEMPIV,ZERO,PIVMAX
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*)
        DATA ZERO/0.0/
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = ABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = ABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = MAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)/COLPIV
              TOPBLK(I,J) = COLMLT
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = MAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = ABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = ABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)/COLPIV
                 ARRAY(I,J,K) = COLMLT
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 345 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
345                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = ABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = ABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = MAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+ABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END
        SUBROUTINE CCRSLV(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B)
C
C***************************************************************
C
C  C C R S L V  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - COMPLEX(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - COMPLEX(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - COMPLEX(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - COMPLEX(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - COMPLEX(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C                    X - COMPLEX(N)
C                         THE SOLUTION VECTOR
C
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C***************************************************************
C
c-----------------------------------------------------------------------
        integer nrwtop,novrlp,nrwblk,nclblk,nbloks,nrwbot,nrwtp1,nrwbk1,
     &          nvrlp1,nrwtp0,nrwbt1,nrowel,nrwel1,nvrlp0,nblks1,nbktop,
     &          j,loop,i,incr,k,incrtp,incrj,incri,jpivot,jrwtop,ll,
     &          nrwbtl,l,l1,iplusn,incrn,ipvtn,nrwell,ipvti
c-----------------------------------------------------------------------
        DOUBLE COMPLEX TOPBLK,ARRAY,BOTBLK,B
        DOUBLE COMPLEX DOTPRD,XJ,XINCRJ,BINCRJ,SWAP
        INTEGER PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              XJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*XJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
440           CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
510        CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
        RETURN
        END
      subroutine colpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates the piecewise polynomial space breakpoint
c       sequence, and calculates the collocation point sequence.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 3, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points. 
c
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval. 
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x. 
c                               nint >= 1.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
c       Work Storage:
        double precision        work(kcol*kcol)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [a,b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i), 
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
c-----------------------------------------------------------------------
c Local Variables:  
        double precision        rho(mxkcol+1)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+1)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Generate the piecewise polynomial space breakpoint sequence.
      do 10 i = 1, kcol + nconti
         xbs(i) = x(1)
         xbs(i + ncpts) = x(nint + 1)
   10 continue
      do 30 i = 2, nint
         ii = (i - 2) * kcol + kcol + nconti
         do 20 j = 1, kcol
            xbs(ii + j) = x(i)
   20    continue
   30 continue

c-----------------------------------------------------------------------
c     Compute the Gaussian points.
      call gauleg(kcol, kcol*kcol, rho, wts, work, 2)

c     Define the collocation point sequence.
      xcol(1) = x(1)
      do 50 i = 1, nint
         ii = (i - 1) * kcol + 1
         do 40 j = 1, kcol
            xcol(ii + j) = x(i) + h(i) * rho(j)
   40    continue
   50 continue
      xcol(ncpts) = x(nint + 1)

      return
      end
        SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK
        DOUBLE PRECISION ROWMAX,ROWPIV,ROWMLT,COLMAX,COLPIV
        DOUBLE PRECISION SWAP,COLMLT,PIVMAX,ZERO,TEMPIV
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          IFLAG
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*)
        INTEGER NRWEL1,NROWEL,I,IPVT,J,IPLUS1,L,INCR,K,NRWTP1,JMINN,
     *          LOOP,JPLUS1,INCRJ,IPLUSN,IPVBLK,KPLUS1,IRWBLK,JRWBLK,
     *          NVRLP0,INCRN
        DATA ZERO/0.0D0/
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = DABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = DABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)/COLPIV
              TOPBLK(I,J) = COLMLT
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = DABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)/COLPIV
                 ARRAY(I,J,K) = COLMLT
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 345 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
345                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = DABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+DABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END
        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,JOB)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        DOUBLE PRECISION DOTPRD,BJ,XINCRJ,BINCRJ,SWAP,BI
        INTEGER NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          JOB
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)
        INTEGER NRWTP1,NRWBK1,NVRLP1,NRWBT1,NROWEL,NVRLP0,NBLKS1,
     *          NBKTOP,J,I,LOOP,INCR,INCRJ,INCRI,JPIVOT,JRWTOP,
     *          LL,L1,IPLUSN,INCRN,NRWTP0,NRWEL1,K,INCRTP,NRWBTL,
     *          IPVTN,NRWELL,IPVTI,L
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C       IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH
C       TRANSPOSE(A)*X = B.
C
        IF ( JOB .NE. 0 ) GO TO 530
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              BJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*BJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
440           CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
510        CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
C
C       RETURN FROM THE SOLUTION OF A.X = B.
        RETURN
C
C       IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:
C
  530   CONTINUE

C       FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).

        DO 540 I = 1,NRWTOP
           IPVTI = PIVOT(I)
           IF ( I .NE. IPVTI ) THEN
              SWAP = B(I)
              B(I) = B(IPVTI)
              B(IPVTI) = SWAP
           ENDIF
           BI = -B(I)
           LOOP = I+1
           DO 535 J = LOOP,NOVRLP
              B(J) = B(J) + BI*TOPBLK(I,J)
  535      CONTINUE
  540   CONTINUE

C       IN EACH BLOCK, K = 1,..,NBLOKS:

        INCR = NRWTOP
        DO 590 K = 1,NBLOKS

C          FIRST, THE FORWARD SOLUTION.

           DO 550 J = 1,NROWEL
              INCRJ = INCR + J
              DO 545 I = 1,J-1
                 B(INCRJ) = B(INCRJ) - ARRAY(I,NRWTOP+J,K)*B(INCR+I)
  545         CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(J,NRWTOP+J,K)
  550       CONTINUE

C           FORWARD MODIFICATION.

            DO 570 I = 1,NOVRLP
               INCRI = INCR + NROWEL + I
               LOOP = NRWBLK + I
               DO 560 J = 1,NROWEL
                  INCRJ = INCR + J
                  B(INCRI) = B(INCRI) - ARRAY(J,LOOP,K)*B(INCRJ)
  560          CONTINUE
  570       CONTINUE

C           NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS
C           CORRESPONDS TO THE LOOP 540 ABOVE.

            INCR = INCR + NROWEL
            DO 580 I = 1,NRWTOP
               INCRI = INCR + I
               IPVTI = PIVOT(INCRI)
               IF ( INCRI .NE. IPVTI ) THEN
                  SWAP = B(INCRI)
                  B(INCRI) = B(IPVTI)
                  B(IPVTI) = SWAP
               ENDIF
               LOOP = NROWEL + I
               BI = -B(INCRI)
               DO 575 J = I+1,NOVRLP
                  INCRJ = INCR+J
                  L = NRWBLK + J
                  B(INCRJ) = B(INCRJ) + BI*ARRAY(LOOP,L,K)
  575          CONTINUE
  580       CONTINUE
            INCR = INCR + NRWTOP
  590   CONTINUE

C       FINALLY, FINISH WITH NRWBOT SOLUTIONS:

        DO 600 J = 1,NRWBOT
           INCRJ = INCR + J
           DO 595 I = 1,J-1
              B(INCRJ) = B(INCRJ) - BOTBLK(I,J+NRWTOP)*B(INCR+I)
  595      CONTINUE
           B(INCRJ) = B(INCRJ)/BOTBLK(J,J+NRWTOP)
  600   CONTINUE


C       NOW, THE BACKWARD PASS:


C       FIRST, BACKWARD SOLUTION IN BOTBLK:

        INCRJ = INCR + NRWBOT
        DO 610 J = 1,NRWBOT-1
           INCRJ = INCRJ - 1
           DO 605 I = NRWBOT-J+1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NOVRLP-J)*B(INCRI)
  605      CONTINUE

           IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
              SWAP = B(INCRJ)
              B(INCRJ) = B(PIVOT(INCRJ))
              B(PIVOT(INCRJ)) = SWAP
           ENDIF
  610   CONTINUE

C       NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:

        DO 620 J = 1,NRWTOP
           INCRJ = INCR - J + 1
           DO 615 I = 1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NRWTP1-J)*B(INCRI)
  615      CONTINUE
  620   CONTINUE


C       NOW, IN EACH BLOCK, K = NBLOKS,..,1:
        DO 800 K = NBLOKS,1,-1

C          FIRST, THE BACKSUBSTITUIONS:

           DO 630 J = 1,NRWTOP
              INCRJ = INCR - J + 1
              LOOP = NBKTOP - J + 1
              DO 625 I = 1,J-1
                 INCRI = INCR - I + 1
                 B(INCRJ) = B(INCRJ) - ARRAY(NRWBLK-I+1,LOOP,K)*B(INCRI)
  625         CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(NRWBLK-J+1,LOOP,K)
  630      CONTINUE

C          THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:

           DO 650 J = 1,NROWEL
              INCRJ = INCR - NRWTOP -J + 1
              DO 645 I = 1,J+NRWTOP-1
                 INCRI = INCRJ + I
                 B(INCRJ) = B(INCRJ) -
     *           ARRAY(NRWBLK-NRWTOP-J+1+I,NRWBLK-J+1,K)*B(INCRI)
  645         CONTINUE
              IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(PIVOT(INCRJ))
                 B(PIVOT(INCRJ)) = SWAP
              ENDIF
  650      CONTINUE

C          NOW, THE DEFERRED OPERATIONS ON B:

           INCR = INCR - NRWBLK
           DO 660 J = 1,NRWTOP
              INCRJ = INCR + J - NRWTOP
              DO 655 I = 1,NRWBLK
                 INCRI = INCR + I
                 B(INCRJ) = B(INCRJ) - ARRAY(I,J,K)*B(INCRI)
  655        CONTINUE
  660      CONTINUE
  800   CONTINUE

C       FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:

        DO 900 J = 1,NRWTOP
           INCRJ = NRWTOP -J + 1
           DO 850 I = INCRJ+1,NRWTOP
              B(INCRJ) = B(INCRJ) - TOPBLK(I,INCRJ)*B(I)
  850      CONTINUE
           B(INCRJ) = B(INCRJ)/TOPBLK(INCRJ,INCRJ)
  900   CONTINUE
C
C       RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B

        RETURN
        END
C ***********************************************************
C
      SUBROUTINE DECOMC(N,NPDE,NINT,KCOL,FJAC,FMAS,ALPHN,BETAN,E2R,
     &                  IP2,IER)
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 23, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,IP2(N),IER
      INTEGER NPDE, NINT, KCOL  
      DOUBLE PRECISION FJAC(*),FMAS(*),ALPHN,BETAN
      DOUBLE COMPLEX E2R(*),BB
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J
      INTEGER NEQ,NSIZTB,NSIZBK,KCOLTM,ITEMP,IABDTP,IABDBK,IABDBT  
c-----------------------------------------------------------------------
c Subroutines Called:
c                               ccrcmp
c
c-----------------------------------------------------------------------

      BB = CMPLX(ALPHN,BETAN)
      ITEMP = 0
      DO 100 I = 1, 2
         KCOLTM = KCOL+I-1
         NSIZTB = NPDE*NPDE*NCONTI
         NSIZBK = NPDE*NPDE*KCOLTM*(KCOLTM+NCONTI)*NINT
 
         IABDTP = ITEMP  + 1
         IABDBK = IABDTP + NSIZTB
         IABDBT = IABDBK + NSIZBK
 
         DO 10 J = 1,NSIZTB
            E2R(IABDTP-1+J)=BB*FJAC(IABDTP-1+J)
   10    CONTINUE
 
         DO 20 J = 1, NSIZBK
            E2R(IABDBK-1+J) = -FJAC(IABDBK-1+J)+BB*FMAS(IABDBK-1+J)
   20    CONTINUE
 
         DO 30 J = 1,NSIZTB
            E2R(IABDBT-1+J)=BB*FJAC(IABDBT-1+J)
   30    CONTINUE

         NEQ = NPDE*(KCOLTM*NINT+NCONTI)
 
         CALL CCRCMP(NEQ,E2R(IABDTP),NPDE,2*NPDE,E2R(IABDBK),
     &               KCOLTM*NPDE,(KCOLTM+NCONTI)*NPDE,NINT,E2R(IABDBT),
     &               NPDE,IP2(1+(I-1)*NPDE*(KCOL*NINT+NCONTI)),IER)
         ITEMP  = NSIZTB + NSIZBK + NSIZTB
 
  100 CONTINUE
C
      RETURN
      END
C ******************************************
C     VERSION OF SEPTEMBER 18, 1995      
C ******************************************
C
      SUBROUTINE DECOMR(N,NPDE,NINT,KCOL,FJAC,FMAS,FAC1,E1,IP1,IER)
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 5, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,IP1(N),IER
      INTEGER NPDE, NINT, KCOL
      DOUBLE PRECISION FJAC(*),FMAS(*),FAC1,E1(*)
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J
      INTEGER NEQ,NSIZTB,NSIZBK,KCOLTM,ITEMP,IABDTP,IABDBK,IABDBT
c-----------------------------------------------------------------------
c Subroutines Called:
c                               crdcmp
c                               daxpy
c
c-----------------------------------------------------------------------
 
      ITEMP = 0
      DO 100 I = 1, 2
         KCOLTM = KCOL+I-1
         NSIZTB = NPDE*NPDE*NCONTI
         NSIZBK = NPDE*NPDE*KCOLTM*(KCOLTM+NCONTI)*NINT

         IABDTP = ITEMP  + 1
         IABDBK = IABDTP + NSIZTB
         IABDBT = IABDBK + NSIZBK
         DO 10 J = 1,NSIZTB
            E1(IABDTP-1+J) = ZERO
   10    CONTINUE

         CALL DAXPY(NSIZTB, FAC1, FJAC(IABDTP), 1, E1(IABDTP), 1)

         DO 20 J = 1, NSIZBK
            E1(IABDBK-1+J) = - FJAC(IABDBK-1+J)
   20    CONTINUE

         CALL DAXPY(NSIZBK, FAC1, FMAS(IABDBK), 1, E1(IABDBK), 1)

         DO 30 J = 1,NSIZTB
            E1(IABDBT-1+J) = ZERO
   30    CONTINUE
         
         CALL DAXPY(NSIZTB, FAC1, FJAC(IABDBT), 1, E1(IABDBT), 1)

         NEQ = NPDE*(KCOLTM*NINT+NCONTI)

         CALL CRDCMP(NEQ,E1(IABDTP),NPDE,2*NPDE,E1(IABDBK),KCOLTM*NPDE,
     &               (KCOLTM+NCONTI)*NPDE,NINT,E1(IABDBT),NPDE,
     &               IP1(1+(I-1)*NPDE*(KCOL*NINT+NCONTI)),IER)

         ITEMP  = NSIZTB + NSIZBK + NSIZTB
  100 CONTINUE
C
      RETURN
      END
      subroutine errest(kcol, nint, npde, neq1, neq2, npts, icount,
     &                  xsol, wts, xbs1, xbs2, y1, y2, istart, mflag2,
     &                  atol, rtol, lenwk, work, errba1, errba2, errrat,
     &                  errint, errcom, ieflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the error estimate at each subinterval
c       and for each component of PDEs, and decides whether a remeshing
c       is necessary or not.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        integer                 nintsm
        parameter              (nintsm = 15)
c                               when the current step is the first step
c                               after remeshing, we require
c                               if nint <= nintsm
c                                    errrat < saffa2
c                               else
c                                    saffa1 < errrat < saffa2.
c                               endif
c
        double precision        zero
        parameter              (zero = 0.0d0)
c
        double precision        one
        parameter              (one = 1.0d0)
c
        double precision        two
        parameter              (two = 2.0d0)
c
        double precision        saffa1
        parameter              (saffa1 = 0.1d0)
c
        double precision        saffa2
        parameter              (saffa2 = 0.4d0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq1
c                               neq1=npde*(nint*kcol+nconti) is the
c                               number of bspline coefficients (or
c                               DAEs) when using dassl_kcol.
c
        integer                 neq2
c                               neq2=neq1+npde*nint is the number of
c                               bspline coefficients (or DAEs) when
c                               using dassl_kcol+1.
c
        integer                 npts
c                               npts is the number of points in the
c                               x vector, which is equal to
c                               nint*(kcol+3).
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xsol(npts)
c                               xsol is the npts Gauss-Legend 
c                               points at which the solution are
c                               to be calculated.
c
        double precision        wts(npts)
c                               wts is the npts Gauss-Legend 
c                               weights at the corresponding xsol.
c
        double precision        xbs1((kcol+1)*nint+nconti+nconti)
c                               xbs1 is the breakpoint sequence when
c                               using dassl_kcol.
c
        double precision        xbs2((kcol+2)*nint+nconti+nconti)
c                               xbs2 is the breakpoint sequence when
c                               using dassl_kcol+1.
c
        double precision        y1(neq1)
c                               y1 is the vector of bspline 
c                               coefficients when using dassl_kcol.
c
        double precision        y2(neq2)
c                               y2 is the vector of bspline 
c                               coefficients when using dassl_kcol+1.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 mflag2
c                               mflag2 = 0, scalar atol and rtol.;
c                               mflag2 = 1, vector atol and rtol.
c
        double precision        atol(npde)
c                               atol is the absolute error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        double precision        rtol(npde)
c                               rtol is the relative error tolerance
c                               request and is a scalar quantity if
c                               mflag2 = 0.
c
        integer                 lenwk
c                               lenwk is the size of the work storage
c                               array and must satisfy:
c                               lenwk >= 2*npde*nint*(kcol+3)
c                                        +npde*nint
c
c       Work Storage:
        double precision        work(lenwk)
c                               work is a floating point work storage
c                               array of size lenwk.
c
c       output:
        double precision        errba1((kcol+nconti)*npts)
c                               errba1 is the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol.
c
        double precision        errba2((kcol+1+nconti)*npts)
c                               errba2 is the values of the nonzero
c                               basis functions at xsol when using
c                               dassl_kcol+1.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of errcom.
c
        double precision        errint(nint)
c                               errint is the error estimate at each
c                               subinterval.
c
        double precision        errcom(npde)
c                               errcom is the error estimate for
c                               each component of pdes at the whole 
c                               range, i.e. from x_a to x_b.
c
        integer                 ieflag
c                               ieflag is a status flag for remesh.
c                               ieflag = 0, indicates no need remeshing.
c                               ieflag = 1, indicates need remeshing.
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        errsum
c                               errsum is the sum of errint.
c
        double precision        errmax
c                               errmax is the maximum value of
c                               errint(i), i = 1, nint.
c
        double precision        aerr
c                               aerr is the average value of errint(i),
c                               i = 1, nint.
c
        double precision        disind
c                               disind is equal to errmax/aerr, and it
c                               indicates the error distribution over
c                               the mesh.
c
c       Pointers into the floating point work array:
        integer                 iusol1
c                               work(iusol1) stores the values at the
c                               npts points when using dassl_kcol.
c
        integer                 iusol2
c                               work(iusol2) stores the values at the
c                               npts points when using dassl_kcol+1.
c
        integer                 ierrci
c                               work(ierrci) stores the error estimate
c                               at each subinterval for each component.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j, m, ij, im, mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               errval
c
c-----------------------------------------------------------------------

c     Set the pointers into the floating point work array.
      iusol1 = 1
      iusol2 = iusol1 + npde * nint * (kcol + 3)
      ierrci = iusol2 + npde * nint * (kcol + 3)

c-----------------------------------------------------------------------
c     Generate the different values at the npts points xsol, and save 
c     in work(iusol1) and work(iusol2).
      call errval(kcol, nint, npde, neq1, kcol+3, istart, icount, xbs1,
     &            xsol, y1, errba1, work(iusol1))
      call errval(kcol+1, nint, npde, neq2, kcol+3, istart, icount,
     &            xbs2, xsol, y2, errba2, work(iusol2))
c-----------------------------------------------------------------------
c     Initialization task.
      do 10 i = 1, nint
         errint(i) = zero
   10 continue

      do 20 i = 1, npde
         errcom(i) = zero
   20 continue

      do 30 i = 1, npde * nint
         work(ierrci - 1 + i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Calculate the error estimate at each subinterval for each
c     component of PDEs.
      if (mflag2 .eq. 0) then
         do 60 m = 1, npde
            do 50 i = 1, nint
               do 40 j = 1, kcol + 3
                  ij = (i - 1) * (kcol + 3) + j
                  mm = npde * (ij - 1) + m 
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im) + ((work(iusol1-1+mm) 
     &                       - work(iusol2-1+mm)) /
     &                       (atol(1) + rtol(1)*abs(work(iusol1-1+mm))))
     &                       **2 * wts(ij)
   40          continue
   50       continue
   60    continue
      else
         do 90 m = 1, npde
            do 80 i = 1, nint
               do 70 j = 1, kcol + 3
                  ij = (i - 1) * (kcol + 3) + j
                  mm = npde * (ij - 1) + m 
                  im = ierrci - 1 + (m - 1) * nint + i
                  work(im) = work(im) + ((work(iusol1-1+mm) 
     &                       - work(iusol2-1+mm)) /
     &                       (atol(m) + rtol(m)*abs(work(iusol1-1+mm))))
     &                       **2 * wts(ij)
   70          continue
   80       continue
   90    continue
      endif

c-----------------------------------------------------------------------
c     Calculate errint and errcom.
      do 110 j = 1, npde
         do 100 i = 1, nint
            ij = ierrci - 1 + (j - 1) * nint + i
            errint(i) = errint(i) + work(ij)
            errcom(j) = errcom(j) + work(ij)
  100    continue
  110 continue

c     Take the square root and update errint and errcom.
      do 120 i = 1, nint
         errint(i) = sqrt(errint(i))
         errint(i) = errint(i) ** (one/dble((kcol+2)))
c        errint(i) = errint(i) ** (one/dble(2*(kcol+2)))
  120 continue

      do 130 i = 1, npde
         errcom(i) = sqrt(errcom(i))
  130 continue

c-----------------------------------------------------------------------
c     Decide whether remeshing is needed.
      ieflag = 0

c     update errrat.
      errrat = zero
      do 140 i = 1, npde
         if (errcom(i) .gt. errrat) then
            errrat = errcom(i)
         endif
  140 continue

c     Calculate errsum to be the sum of the errint. Find the maximum
c     errint(i) and save it in errmax.
      errsum = errint(1)
      errmax = errint(1)
      do 150 i = 2, nint
         if (errmax .lt. errint(i)) errmax = errint(i)
         errsum = errint(i) + errsum
  150 continue
 
c     Let aerr be the mean value of errint(i).
      aerr = errsum/dble(nint)
      
c     Calculate disind.
      disind = errmax/aerr

      if (disind .gt. two) then
         ieflag = 1
      else
         if ((istart .ne. 1) .or. (icount .ne. 0)) then
            if (nint .gt. nintsm) then
               if ((errrat .ge. saffa2) .or. (errrat .le. saffa1)) 
     &            ieflag = 1
            else
               if (errrat .ge. saffa2) ieflag = 1
            endif
         else
            if (errrat .ge. one) ieflag = 1
         endif
      endif

      return
      end
      subroutine errval(kcol, nint, npde, neq, nptse, istart, icount,
     &                  xbs, xsol, y, errbas, usol)
 
c-----------------------------------------------------------------------
c Purpose:
c       This routine computes the values of the (kcol+nconti) nonzero
c       bspline basis function at each Gaussian point of xsol.
c       Then determine the solution usol, which is used for error
c       estimate, at xsol.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 29, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               nint >= 1.
c
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 neq
c                               neq=npde*(kcol*nint+2) is the number of
c                               bspline coefficients.
c
        integer                 nptse
c                               nptse is the number of Gaussian points
c                               in each subinterval for the error
c                               estimate.
c
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        double precision        xbs((kcol+1)*nint+nconti+nconti)
c                               The breakpoint sequence.
c
        double precision        xsol(nptse*nint)
c                               xsol is a set of spatial points at which
c                               the solution are to be calculated for
c                               error estimate.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients.
c
c       output:
        double precision        errbas(kcol+nconti, nptse*nint)
c                               errbas is the values of the nonzero
c                               basis functions at xsol.
c
        double precision        usol(npde, nptse*nint)
c                               uval is the solution at xsol.
c
c-----------------------------------------------------------------------
c Local Variables:                                                              
        integer                 ileft
c                               breakpoint information.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 jj
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c
c-----------------------------------------------------------------------

c     check whether errbas is necessary to be calculated.
      if ((istart .eq. 1) .and. (icount .eq. 0)) goto 30

c     calculate errbas.
      do 20 i = 1, nint
         ileft = kcol + nconti + (i - 1) * kcol
         do 10 j = 1, nptse
            jj = (i - 1) * nptse + j
            call bsplvd(xbs, kcol+nconti, xsol(jj), ileft, errbas(1,jj),
     &                  1)
   10    continue
   20 continue

   30 continue

c     compute the values of usol at xsol.
      do 70 i = 1, nint
         do 60 j = 1, nptse
            jj = (i - 1) * nptse + j
            do 50 k = 1, npde
               usol(k,jj) = zero
               do 40 m = 1, kcol + nconti
                  mm = npde * (m + (i - 1) * kcol - 1) + k
                  usol(k,jj) = usol(k,jj) + y(mm) * errbas(m,jj)
   40          continue
   50       continue
   60    continue
   70 continue

      return
      end
C ***********************************************************
C
      SUBROUTINE ESTRAD(N,NPDE,KCOL,NINT,FMAS,H,DD1,DD2,DD3,FCN,
     &          NFCN,Y0,Y,X,E1,Z1,Z2,Z3,CONT,F1,F2,IP1,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 30, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        one
        parameter              (one = 1.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,NPDE,KCOL,NINT,NFCN,IP1(N),IPAR(*)
      DOUBLE PRECISION FMAS(*),H,DD1,DD2,DD3,Y0(N),Y(N),X,E1(*),Z1(N),
     &                 Z2(N),Z3(N),CONT(N),F1(N),F2(N),SCAL(N),ERR,
     &                 FAC1,RPAR(*)
      LOGICAL FIRST,REJECT
      EXTERNAL FCN
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J,K,M,L,III,KKK,MMM,NNN,II,KK,MM,NEQ1,KCOLTM,NPDTP1,
     &        NPDBK1,NPDBT1,NPDTP2,NPDBK2,NPDBT2
      DOUBLE PRECISION HEE1,HEE2,HEE3,SUM
c-----------------------------------------------------------------------
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
C
      DO 10 I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
   10 CONTINUE

      neq1 = npde*(kcol*nint+nconti)
      nnn  = npde*npde*(2*nconti+kcol*(kcol+nconti)*nint)
 
      do 80 l = 1, 2
         kcoltm = kcol + l - 1
         iii = npde*npde*kcoltm
         mmm = (l-1)*neq1
         kkk = (l-1)*nnn
         do 50 i = 1, nint
            do 40 k = 1, kcoltm
               do 30 m = 1, npde
                  sum = zero
                  ii = mmm+npde+(i-1)*npde*kcoltm+(k-1)*npde+m
                  do 20 j = 1, kcoltm + nconti
                     kk = kkk+1+(i-1)*iii*(kcoltm+nconti)+(j-1)*iii
     &                    +(k-1)*npde+npde*npde*nconti
                     mm = mmm+(i-1)*kcoltm*npde+(j-1)*npde+m
                     sum = sum + fmas(kk) * f1(mm)
   20             continue
                  f2(ii)   = sum
                  cont(ii) = sum + y0(ii)
   30          continue
   40       continue
   50    continue
         do 60 i = 1, npde
            f2(i+mmm) = zero
            cont(i+mmm) = - y0(i+mmm) * fac1
   60    continue
         do 70 i = 1, npde
            ii = i+npde*(kcoltm*nint+1)
            f2(ii+mmm) = zero
            cont(ii+mmm) = - y0(ii+mmm) * fac1 
   70    continue
   80 continue        

      npdtp1 = 1
      npdbk1 = npdtp1 + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)
 
      npdtp2 = npdbt1 + npde * npde * nconti
      npdbk2 = npdtp2 + npde * npde * nconti
      npdbt2 = npdbk2 + npde * npde * nint * (kcol + 1) * (kcol + 1 +
     &                  nconti)
 
      call crslve(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &           (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1, cont,
     &           0)
 
      call crslve(e1(npdtp2), npde, 2*npde, e1(npdbk2), (kcol+1)*npde,
     &            (kcol+1+nconti)*npde, nint, e1(npdbt2), npde,
     &            ip1(neq1+1), cont(neq1+1), 0)

      err = zero

      do 90 i = 1, n
         err = err + (cont(i)/scal(i))**2
   90 continue

      err = max(sqrt(err/n), 1.d-10)

      if (err .lt. one) return

      if (first .or. reject) then

         do 100 i = 1, n
            cont(i) = y(i) + cont(i)
  100    continue
         call fcn(n, x, cont, f1, rpar, ipar)
         nfcn = nfcn + 1
         do 110 i = 1, n
            cont(i) = f1(i) + f2(i)
  110    continue
         do 120 i = 1, npde
            cont(i) = - cont(i) * fac1
  120    continue
         do 130 i = neq1-npde+1, neq1
            cont(i) = - cont(i) * fac1 
  130    continue
         do 140 i = neq1+1, neq1+npde
            cont(i) = - cont(i) * fac1
  140    continue
         do 150 i = n-npde+1, n
            cont(i) = - cont(i) * fac1 
  150    continue

         call crslve(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &               (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1,
     &               cont, 0)
 
         call crslve(e1(npdtp2), npde, 2*npde, e1(npdbk2),
     &               (kcol+1)*npde, (kcol+1+nconti)*npde, nint,
     &               e1(npdbt2), npde, ip1(neq1+1), cont(neq1+1), 0)

         err = zero
       
         do 160 i = 1, n
            err = err + (cont(i)/scal(i))**2
  160    continue

         err = max(sqrt(err/n), 1.d-10)

      endif

      RETURN
      END

      subroutine eval(npde,kcol,ileft,icpt,ncpts,uval,uxval,uxxval,
     &                fbasis,y)

c-----------------------------------------------------------------------
c Purpose:
c       This routine evaluates u(k), ux(k), and uxx(k), k=1 to npde,
c       at the icpt-th collocation point.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Feb. 11, 2001.
c
c-----------------------------------------------------------------------
c Constants
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 icpt
c                               the index of the collocation point.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        double precision        fbasis((kcol+nconti)*3)
c                               Basis function values at the icpt-th
c                               collocation point. 
c                               fbasis(k+(j-1)*(kcol+nconti)) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti). 
c
        double precision        y(ncpts*npde)
c                               y is the vector of bspline coefficients.
c
c       Output:
        double precision        uval(npde)
c                               uval gives the approximation to 
c                               u(t,x).
c
        double precision        uxval(npde)
c                               uxval gives the approximation to 
c                               the first spatial derivative of u(t,x).
c
        double precision        uxxval(npde)
c                               uxxval gives the approximation to 
c                               the second spatial derivative of u(t,x).
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 j
        integer                 m
c-----------------------------------------------------------------------
      do 10 j = 1, npde
         uval(j)   = zero
         uxval(j)  = zero
         uxxval(j) = zero
   10 continue
      if (icpt .ne. 1 .and. icpt .ne. ncpts) then
         do 30 j = 1, npde
            do 20 m = 1, kcol + nconti
               uval(j)   = uval(j) + fbasis(m) * 
     &                     y((ileft-kcol-3+m) * npde + j)
               uxval(j)  = uxval(j) + fbasis(m+kcol+nconti) * 
     &                     y((ileft-kcol-3+m) * npde + j)
               uxxval(j) = uxxval(j) + fbasis(m+2*(kcol+nconti)) * 
     &                     y((ileft-kcol-3+m) * npde + j)
   20          continue
   30       continue
      else 
         if (icpt .eq. 1) then
            do 40 j = 1, npde
               uval(j)   = uval(j)   + fbasis(1) * y(j)
               uxval(j)  = uxval(j)  + fbasis(1+kcol+nconti) * y(j)
     &                     + fbasis(2+kcol+nconti) * y(npde + j)
               uxxval(j) = uxxval(j) + fbasis(1+2*(kcol+nconti)) * y(j)
     &                     + fbasis(2+2*(kcol+nconti)) * y(npde + j)
     &                     + fbasis(3+2*(kcol+nconti)) * y(2*npde + j)
   40       continue
         else
            do 50 j = 1, npde
               uval(j)   = uval(j)   + fbasis(kcol+nconti) 
     &                     * y((ncpts - 1) * npde + j)
               uxval(j)  = uxval(j)  + fbasis((kcol+nconti)*2) 
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*2-1)
     &                     * y((ncpts - 2) * npde + j)
               uxxval(j) = uxxval(j) + fbasis((kcol+nconti)*3) 
     &                     * y((ncpts - 1) * npde + j)
     &                     + fbasis((kcol+nconti)*3-1)
     &                     * y((ncpts - 2) * npde + j)
     &                     + fbasis((kcol+nconti)*3-2)
     &                     * y((ncpts - 3) * npde + j)
   50       continue
         endif
      endif
      return
      end

      SUBROUTINE GAULEG(N, NSQ, PTS, WTS, WORK, FLAG)
C $ID: GAULEG.F,V 1.11 1992/06/25 15:09:31 KEAST EXP $
C LAST MODIFIED BY RONG WANG, 2001/03/08
C
C     GAULEG RETURNS THE POINTS AND WEIGHTS FOR GAUSS-LEGENDRE
C     QUADRATURE OR GAUSS-LOBATTO QUADRATURE OVER THE INTERVAL 
C     [-1,1] OR [0,1].
C
C     ON INPUT:
C        
C        N     : THE NUMBER OF GAUSS-LEGENDRE POINTS.
C        NSQ   : EQUAL TO N*N, TO HANDLE STORAGE FOR EIGENVECTORS.
C        PTS   : DOUBLE PRECISION (N).
C        WTS   : DOUBLE PRECISION (N).  IF FLAG (SEE BELOW) IS 1 OR 3,
C                WTS IS USED ONLY FOR TEMPORARY WORKSPACE.
C        WORK  : DOUBLE PRECISION (NSQ), WORK SPACE FOR CALL TO IMTQL2,
C                IF WEIGHTS ARE REQUIRED. IF FLAG IS 1 OR 3, WORK IS
C                NOT REFERENCED, AND MAY BE DECLARED AS SCALAR IN THE 
C                CALLING PROGRAM.
C        FLAG  : SPECIFIES WHETHER WEIGHTS ARE ALSO REQUIRED, AND
C                WHETHER GAUSS-LEGENDRE OR LOBATTO POINTS ARE WANTED.
C                   FLAG  = 1: GAUSS-LEGENDRE POINTS ONLY OVER [-1,1];
C                         = 2: GAUSS-LEGENDRE POINTS ONLY OVER [0,1];
C                         = 3: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [-1,1];
C                         = 4: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [0,1];
C                         = 5: LOBATTO POINTS ONLY OVER [-1,1];
C                         = 6: LOBATTO POINTS ONLY OVER [0,1];
C                         = 7: LOBATTO POINTS AND WEIGHTS OVER [-1,1];
C                         = 8: LOBATTO POINTS AND WEIGHTS OVER [0,1];
C                FOR ANY OTHER VALUE, THE DEFAULT IS GAUSS-LEGENDRE
C                POINTS AND WEIGHTS OVER [-1,1].
C
C     ON OUTPUT:
C
C        FOR FLAG <> 1 OR 2:
C          PTS : PTS(I) IS THE ITH GAUSS-LEGENDRE POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C          WTS : WTS(I) IS THE ITH GAUSS-LEGENDRE WEIGHT IF FLAG <> 1.
C
C        FOR FLAG = 3 OR 4:
C          PTS : PTS(I) IS THE ITH LOBATTO POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C                CLEARLY, PTS(1) = -1.0, PTS(N) = 1.0.
C          WTS : WTS(I) IS THE ITH LOBATTO WEIGHT IF FLAG = 4.
C
C        WORK  : WORK, USED TO STORE EIGENVECTORS, UNREFERENCED IF FLAG
C                IS 1 OR 3.
C
C     SUBROUTINES USED:
C
C        IMQTL1: EISPACK ROUTINE TO COMPUTE THE EIGENVALUES OF A 
C                SYMMETRIC TRIDIAGONAL MATRIX.
C
C        IMQTL2: EISPACK ROUTINE TO COMPUTE THE EIGNEVECTORS AND 
C                EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX.
C
C     FUNCTIONS USED:
C
C        PYTHAG: EISPACK FUNCTION TO COMPUTE EUCLIDEAN NORM.
C
C     INTRINSIC FUNCTIONS USED:
C
C        SQRT, DBLE.
C
C     VERSION: JUNE 22 1992, PAT KEAST.
C
C     DECLARATIONS:
C
C        PARAMETERS:
C
      INTEGER          N, NSQ, FLAG
      DOUBLE PRECISION PTS(N), WTS(N), WORK(NSQ)
C
C        LOCAL VARIABLES:

      INTEGER          IFAIL, J, NM2
      DOUBLE PRECISION FOUR, THREE, TWO, ONE, ZERO
*     .. EXTERNAL FUNCTIONS ..
      EXTERNAL         IMTQL2 
      PARAMETER ( FOUR = 4.0D0, THREE = 3.0D0, TWO = 2.0D0, 
     *            ONE = 1.0D0, ZERO = 0.0D0 )
  
      DO 10 J = 1,N
         PTS(J) = ZERO
   10 CONTINUE
 
      DO 20 J = 1,NSQ
         WORK(J) = ZERO
   20 CONTINUE

      DO 30 J = 1,N
         WORK((J-1)*N+J) = ONE
   30 CONTINUE
 
      IF ( FLAG .LE.4 .OR. FLAG .GT. 8 ) THEN
C        GAUS-LEGENDRE POINTS AND WEIGHTS.
         DO 40 J = 1,N-1
            WTS(J+1) = DBLE(J)/SQRT(DBLE(4*J*J)-ONE)
   40    CONTINUE
 
         IF ( FLAG .EQ. 1 .OR. FLAG .EQ. 2 ) THEN
C           COMPUTE ONLY THE GAUSS-LEGENDRE POINTS OVER [-1,1].
            CALL IMTQL1(N, PTS, WTS, IFAIL )
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 2 ) THEN
               DO 45 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   45          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS OVER [-1,1].
            IFAIL = 1
C
            CALL IMTQL2(N, N, PTS, WTS, WORK, IFAIL)
    
            DO 50 J = 1,N
               WTS(J) = TWO*WORK((J-1)*N+1)*WORK((J-1)*N+1)
   50       CONTINUE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 4 ) THEN
               DO 55 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   55          CONTINUE
            ENDIF
         ENDIF
C
      ELSE
C        THE LOBATTO POINTS AND WEIGHTS.

C        FIRST, COMPUTE THE ORDER N-2 JACOBI POINTS AND/OR WEIGHTS.
         NM2 = N-2
         DO 60 J = 1,NM2-1
            WTS(J+1) = SQRT(DBLE(J*(J+2))/DBLE((2*J+1)*(2*J+3)))
   60    CONTINUE

         IF ( FLAG .EQ. 5 .OR. FLAG .EQ. 6) THEN
C           COMPUTE ONLY THE GAUSS-LOBATTO POINTS OVER [-1,1].
            CALL IMTQL1(NM2, PTS(2), WTS, IFAIL )
            PTS(1) = -ONE
            PTS(N) = ONE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 6 ) THEN
               DO 65 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   65          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS.
            IFAIL = 1
C
            CALL IMTQL2(NM2, NM2, PTS(2), WTS, WORK, IFAIL)
            PTS(1) = -ONE
            PTS(N) = ONE
    
            DO 70 J = 2,N-1
               WTS(J) = (FOUR/THREE)*WORK((J-2)*NM2+1)*WORK((J-2)*NM2+1)
     *                  /(ONE - PTS(J)*PTS(J))
   70       CONTINUE
            WTS(1) = ZERO
            DO 80 J = 2,N-1
               WTS(1) = WTS(1) - WTS(J)
   80       CONTINUE
            WTS(1) = ONE + WTS(1)/TWO
            WTS(N) = WTS(1)
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 8 ) THEN
               DO 85 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   85          CONTINUE
            ENDIF
         ENDIF
C
         RETURN
      ENDIF
*
      RETURN
      END
      subroutine iniy(t0, npde, kcol, nint, neq, ncpts, ifglin, xcol,
     &                xbs, abdblk, fbasis, y, ipivot, work, lw, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks required by 
c       inital including:
c
c               calculating the Bspline basis functions,
c               constructing abdblk of the collocation matrices and
c               determining y(t0).
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, November 8, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points. 
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
        double precision        negone
        parameter              (negone = -1.0D0)
c                                                                               
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        double precision        t0
c                               t0 is the initial time.
c   
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval. 
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x. 
c                               nint >= 1.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bsplines 
c                               coefficients (or DAEs).
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        integer                 ifglin
c                               ifglin is a flag for the boundary
c                               conditions.
c                               ifglin = 1, indicate both derichlet
c                                           boundary conditions;
c                                      = 0, else.
c
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c                               xbs(i)=x(1), i=1, kcol+nconti;
c                               xbs((i-1)*kcol+nconti+j)=x(i), 
c                                    i=2, nint;  j=1, kcol
c                               xbs(ncpts+i)=x(nint+1), i=1,kcol+nconti.
c
        integer                 lw
c                               lw is the size of the work storage 
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +2*neq+2*npde+2*npde*npde
c
c       Work Storage:
        integer                 ipivot(neq)
c                               pivoting information from the 
c                               factorization of the temporary matrix.
c
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
c       Output:
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. fbasis(k,j,i) contains the
c                               values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        double precision        y(neq)
c                               y = y(t0) is the initial vector of
c                               bspline coefficients. 
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:  
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one 
c                               collocation block of work.
c 
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the 
c                               bottom block which is required since 
c                               crdcmp overwrites the input collocation 
c                               matrix.
c
        integer                 idelta
c                               work(idelta) contains the residual which
c                               indicates how well y satisfies to the
c                               boundary condition and the initial
c                               condition at the internal collocation
c                               points.
c
        integer                 ivcol
c                               work(ivcol) contains the values of u
c                               at the internal collocation points.
c
        integer                 iu
c                               work(iu) stores the approximation to
c                               u(t,x).
c
        integer                 iux
c                               work(iux) stores the approximation to
c                               the first spatial derivative of u(t,x).
c
        integer                 iuxx
c                               work(iuxx) stores the approximation to
c                               the second spatial derivative of u(t,x).
c
        integer                 idbdu
c                               work(idbdu-1+i), i=1, npde*npde,
c                               contains dbdu(npde,npde). That is,
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        integer                 idbdux
c                               work(idbdux-1+i), i=1, npde*npde,
c                               contains dbdux(npde,npde). That is,
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        integer                 idbdt
c                               work(idbdt-1+i), i=1, npde, contains
c                               the partial derivative of the i-th              c                               component of the vector b with respect
c                               to time t.
c           
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 jj
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bndxa
c                               bndxb
c                               bsplvd
c                               difbxa
c                               difbxb
c                               eval
c                               uinit
c                               crdcmp
c                               crslve
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c                               dscal
c
c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      idelta = iabdbt + npde*npde*nconti
      ivcol  = idelta + neq
      iu     = ivcol  + neq-2*npde
      iux    = iu     + npde
      iuxx   = iux    + npde
      idbdu  = iuxx   + npde
      idbdux = idbdu  + npde*npde
      idbdt  = idbdux + npde*npde

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 20 i = 1, npde*npde*nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   20 continue
      do 30 i = 1, nint*nels
         abdblk(i) = zero
   30 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j) 
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Uinit is called to evaluate the first npde components at the
c     left boundary point, and save in y.
      call uinit(xcol(1), y(1), npde)

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 40 i = 1, npde 
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   40 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 80 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 70 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     jj is the position in the y vector where the values for the 
c     right hand side of the initial conditions, evaluated at the ii-th
c     collocation point are stored.
            jj = (ii - 1) * npde + 1

c     compute information for ii-th collocation point.
            call bsplvd(xbs,kcol+nconti,xcol(ii),ileft,fbasis(1,1,ii),3)
            call uinit(xcol(ii), y(jj), npde)

            do 60 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 50 m = 1, npde
                  mm = ll + (m -1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   50          continue
   60       continue
   70    continue
   80 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      ii = neq - npde + 1
      call uinit(xcol(ncpts), y(ii), npde)
      do 90 i = 1, npde
         ii = ((i-1)+npde)*npde + i  
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   90 continue

c-----------------------------------------------------------------------   
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)

c     Check whether both boundary conditions are derichlet boundary
c     conditions. If no, copy the values at the internal collocation
c     points to work(ivcol), which will be used for newton iterations.
      if (ifglin .eq. 0) then
         call dcopy(neq-2*npde,y(npde+1),1,work(ivcol),1)
         call dscal(neq-2*npde,negone,work(ivcol),1)
      endif

c-----------------------------------------------------------------------
c     Generate the initial vector y(t0).
c-----------------------------------------------------------------------

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999
 
c     Solve the linear system. If derichlet boundary conditions are
c     given, this gives the basis function coefficients for the initial
c     conditions, i.e. y(t0). If not, this gives the predictor of y(t0).
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,y,0)

      if (icflag .ne. 0) goto 999
 
c     Check whether both boundary conditions are derichlet boundary
c     conditions.
      if (ifglin .eq. 1) goto 999

c-----------------------------------------------------------------------
c     Newton iteration loop.

c     Calculate (work(idelta-1+i), i = npde+1, neq-npde), which depends
c     on the nint blocks in the middle of the collocation matrix A.
      call dcopy(neq-2*npde,work(ivcol),1,work(idelta+npde),1)
      do 130 i = 1, nint
         do 120 j = 1, kcol + nconti
            do 110 l = 1, kcol
               ll = 1+(i-1)*npde*npde*kcol*(kcol+nconti)
     &              +(j-1)*npde*npde*kcol+(l-1)*npde
               do 100 m = 1, npde
                  ii = idelta-1+npde+(i-1)*npde*kcol+(l-1)*npde+m
                  mm = (i-1)*kcol*npde+(j-1)*npde+m
                  work(ii) = work(ii) + abdblk(ll) * y(mm)
  100          continue
  110       continue
  120    continue
  130 continue
                                                                                
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nint*nels,abdblk,1,work(iabdbk),1)
 
c     Update the values at the left boundary.
      call eval(npde,kcol,kcol+2,1,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,1),y)
      call bndxa(t0, work(iu), work(iux), work(idelta), npde)  
      call difbxa(t0, work(iu), work(iux), work(idbdu),
     &            work(idbdux), work(idbdt), npde)
 
c     Set up the top block and save in work(iabdtp).
      do 150 j = 1, npde
         do 140 i = 1, npde
            ii = iabdtp - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i
            work(jj) = fbasis(2,2,1) * work(idbdux-1+mm)
            work(ii) = work(idbdu-1+mm) - work(jj)             
  140    continue
  150 continue  

c     Update the values at the right boundary.
      call eval(npde,kcol,ncpts,ncpts,ncpts,work(iu),work(iux),
     &          work(iuxx),fbasis(1,1,ncpts),y)
      call bndxb(t0, work(iu), work(iux), work(idelta+neq-npde), npde) 
      call difbxb(t0,work(iu),work(iux),work(idbdu),
     &            work(idbdux),work(idbdt),npde)
 
c     Set up the bottom block and save in work(iabdbt).
      do 170 j = 1, npde
         do 160 i = 1, npde
            ii = iabdbt - 1 + (j - 1) * npde + i
            jj = ii + npde * npde
            mm = (j - 1) * npde + i  
            work(ii) = fbasis(kcol+1,2,ncpts) * work(idbdux-1+mm)
            work(jj) = work(idbdu-1+mm) - work(ii)
  160    continue
  170 continue 

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) goto 999
 
c     Solve the corrector equation.  
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            work(idelta),0)

      if (icflag .ne. 0) goto 999
 
c     Now generate the corrector of y(t0).
      do 180 i = 1, neq
         y(i) = y(i) - work(idelta-1+i)
  180 continue

c-----------------------------------------------------------------------

  999 return
      end
      SUBROUTINE INTERV ( XT, LXT, X, ILEFT, MFLAG, ILO)
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C COMPUTES LARGEST ILEFT IN (1,LXT) SUCH THAT XT(ILEFT) .LE. X.  THE
C PROGRAM STARTS THE SEARCH FOR ILEFT WITH THE VALUE OF ILEFT THAT WAS
C RETURNED AT THE PREVIOUS CALL (AND WAS SAVED IN THE LOCAL VARIABLE
C ILO) TO MINIMIZE THE WORK IN THE COMMON CASE THAT THE VALUE OF X ON
C THIS CALL IS CLOSE TO THE VALUE OF X ON THE PREVIOUS CALL.  SHOULD
C THIS ASSUMPTION NOT BE VALID, THEN THE PROGRAM LOCATES ILO AND IHI
C SUCH THAT XT(ILO) .LE. X .LT. XT(IHI) AND, ONCE THEY ARE FOUND USES
C BISECTION TO FIND THE CORRECT VALUE FOR ILEFT.
C
C LAST MODIFIED BY RONG WANG, JAN 9, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  NONE
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   NONE
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER LXT,ILEFT,MFLAG
      DOUBLE PRECISION XT(LXT),X
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER ILO,IHI,ISTEP,MIDDLE
C-----------------------------------------------------------------------
      IF(MFLAG.EQ.-2) ILO = 1
      IHI = ILO + 1
      IF (IHI .LT. LXT) GO TO 20
      IF (X .GE. XT(LXT)) GO TO 110
      IF (LXT .LE. 1) GO TO 90
      ILO = LXT - 1
      GO TO 21
   20 IF (X .GE. XT(IHI)) GO TO 40
   21 IF (X .GE. XT(ILO)) GO TO 100
C-----------------------------------------------------------------------
C NOW X .LT. XT(IHI).  FIND LOWER BOUND.
C-----------------------------------------------------------------------
      ISTEP = 1
   31 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO .LE. 1) GO TO 35
      IF (X .GE. XT(ILO)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1)) GO TO 90
      GO TO 50
C-----------------------------------------------------------------------
C NOW X .GE. XT(ILO).  FIND UPPER BOUND.
C-----------------------------------------------------------------------
   40 ISTEP = 1
   41 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .GE. LXT) GO TO 45
      IF (X .LT. XT(IHI)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 41
   45 IF (X .GE. XT(LXT)) GO TO 110
      IHI = LXT
C-----------------------------------------------------------------------
C NOW XT(ILO) .LE. X .LT. XT(IHI).  NARROW THE INTERVAL.
C-----------------------------------------------------------------------
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO) GO TO 100
C-----------------------------------------------------------------------
C NOTE..  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1.
C-----------------------------------------------------------------------
      IF (X .LT. XT(MIDDLE)) GO TO 53
      ILO = MIDDLE
      GO TO 50
   53 IHI = MIDDLE
      GO TO 50
C-----------------------------------------------------------------------
C SET OUTPUT AND RETURN.
C-----------------------------------------------------------------------
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END
      subroutine meshsq(kcol, nint, x, work, h, excol, ewts)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mesh size sequence, then generates
c       the collocation points and Gaussian weights for error
c       estimate.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 5, 2001.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 mxkcol
        parameter              (mxkcol = 10)
c                               mxkcol is the maximum number of
c                               collocation points per subinterval.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval. 
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x. 
c                               nint >= 1.
c
        double precision        x(nint+1)
c                               x is the spatial mesh which divides the
c                               interval [x_a, x_b] as: x_a = x(1) <
c                               x(2) < x(3) < ... < x(nint+1) = x_b.
c
c       Work Storage:
        double precision        work((kcol+3)*(kcol+3))
c                               work is a floating point work storage
c                               array.
c
c       Output:
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        double precision        excol(nint*(kcol+3))
c                               excol is the collocation point sequence
c                               which is used for error estimate.
c
        double precision        ewts(nint*(kcol+3))
c                               ewts is the gaussian weight sequence
c                               which is used for error estimate.
c
c-----------------------------------------------------------------------
c Local Variables:  
        double precision        rho(mxkcol+3)
c                               rho stores the Gaussian points.
c
        double precision        wts(mxkcol+3)
c                               wts stores the Gaussian weights.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 ii
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               gauleg
c
c-----------------------------------------------------------------------
c     Calculate the mesh step size sequence.
      do 10 i = 1, nint
         h(i) = x(i+1)-x(i)
   10 continue

c     Compute the Gaussian points and Gaussian weights.
      call gauleg(kcol+3, (kcol+3)*(kcol+3), rho, wts,
     &            work, 4)

c     Define the collocation point sequence.
      do 30 i = 1, nint
         ii = (i - 1) * (kcol + 3)
         do 20 j = 1, kcol+3
            excol(ii + j) = x(i) + h(i) * rho(j)
            ewts(ii + j) = h(i) * wts(j)
   20    continue
   30 continue

      return
      end
      DOUBLE PRECISION FUNCTION PYTHAG(A,B) 
      DOUBLE PRECISION A,B 
C 
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW 
C 
      DOUBLE PRECISION P,Q,R,S,T 
      P = DMAX1(DABS(A),DABS(B)) 
      Q = DMIN1(DABS(A),DABS(B)) 
      IF (Q .EQ. 0.0D0) GO TO 20 
   10 CONTINUE 
         R = (Q/P)**2 
         T = 4.0D0 + R 
         IF (T .EQ. 4.0D0) GO TO 20 
         S = R/T 
         P = P + (2.0D0*S)*P 
         Q = S*Q 
      GO TO 10 
   20 PYTHAG = P 
      RETURN 
      END 
      SUBROUTINE RADAU5(N,FCN,X,Y,XEND,H,RTOL,ATOL,ITOL,JAC,MAS,SOLOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,CWORK,IDID)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
C     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
C                     M*Y'=F(X,Y).
C     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA)
C     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT.
C     CF. SECTION IV.8
C
C     AUTHORS: E. HAIRER AND G. WANNER
C              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C              CH-1211 GENEVE 24, SWITZERLAND 
C              E-MAIL:  Ernst.Hairer@math.unige.ch
C                       Gerhard.Wanner@math.unige.ch
C     
C     THIS CODE IS PART OF THE BOOK:
C         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
C         SPRINGER-VERLAG 1991, SECOND EDITION 1996.
C      
C     VERSION OF JULY 9, 1996
C        (small correction April 14, 2000)
C
C     INPUT PARAMETERS  
C     ----------------  
C     N           DIMENSION OF THE SYSTEM 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),F(N)
C                    F(1)=...   ETC.
C                 RPAR, IPAR (SEE BELOW)
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, 
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y.
C                 THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(N,X,Y,DFY,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),DFY(*)
C                    DFY(1,1)= ...
C
C     ----   MAS HAS ANALOG MEANINGS      -----
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): -
C
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
C                 MATRIX M.
c-----------------------------------------------------------------------  
c                 SUBROUTINE MAS(AM,RPAR,IPAR) 
c                 the mass-matrix AM is stored in the ABD form.
c-----------------------------------------------------------------------  
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (N,X,Y,ITOL,ATOL,RTOL,
C                                       RPAR,IPAR,IRTRN)
C                    DOUBLE PRECISION X,Y(N),CONT(LRC)
C                    ....  
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM.
C           
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE
C                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE
C                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
C                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE
C                 FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST
C                            3*lenpd+12*N+20
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
C                 IWORK(20) TO ZERO BEFORE CALLING.
C                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA.
C                 "LIWORK" MUST BE AT LEAST 3*N+20.
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH  
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. 
C
C ----------------------------------------------------------------------
C 
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK 
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
C              IT DOES NOT WORK FOR BANDED JACOBIAN, 
C              AND NOT FOR IMPLICIT SYSTEMS. 
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP.
C              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7.
C
C    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION
C              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD.
C              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED.
C              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS
C              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN
C              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.).
C              DEFAULT IS IWORK(4)=0.
C
C       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
C       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
C       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
C       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. 
C       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
C       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
C
C    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR 
C              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
C              DEFAULT IWORK(5)=N.
C
C    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0.
C
C    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0.
C
C    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY
C              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON)
C              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL
C              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1.
C              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS;
C              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES
C              OFTEN SLIGHTLY FASTER RUNS
C
C ----------
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER 
C              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO
C              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.     
C              DEFAULT 0.001D0.
C
C    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
C              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER.
C              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0)
C
C    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE
C              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A
C              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR
C              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE
C              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS
C              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD.
C              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 .
C
C    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(8) <= HNEW/HOLD <= WORK(9)
C              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0
C
C-----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS 
C     ----------------- 
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C 
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
C   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)  
C   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                      OR NUMERICALLY)
C   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS
C   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS
C   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
C   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES
C   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH
C                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS,
C                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED
C-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 1, 2003.
c
c-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS 
C *** *** *** *** *** *** *** *** *** *** *** *** ***
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,ITOL,LWORK,LIWORK,IDID,IWORK(LIWORK),IPAR(*)
      DOUBLE PRECISION X,XEND,H,Y(N),WORK(LWORK)
      DOUBLE PRECISION ATOL(*),RTOL(*),RPAR(*)
      DOUBLE COMPLEX CWORK(*)
c-----------------------------------------------------------------------
c local variables
      INTEGER NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NMAX,NIT,
     &        IEZ1,IEZ2,IEZ3,IEY0,IESCAL,IEF1,IEF2,IEF3,IECON,IEJAC,
     &        IEMAS,IEE1,IEE2R,ISTORE,IEIP1,IEIP2,IECZ2
      DOUBLE PRECISION UROUND,EXPM,QUOT,SAFE,THET,TOLST,FNEWT,QUOT1,
     &                 QUOT2,HMAX,FACL,FACR
c-----------------------------------------------------------------------
c loop indices
      INTEGER I
c-----------------------------------------------------------------------

      LOGICAL ARRET,STARTN,PRED
      EXTERNAL FCN,JAC,MAS,SOLOUT
c-----------------------------------------------------------------------
      INTEGER NCONTI,NPDE,NINT,KCOL,NSZJAC
      PARAMETER (NCONTI = 2)        
      NPDE = IPAR(1)
      KCOL = IPAR(2)
      NINT = IPAR(3)
      NSZJAC = NPDE * NPDE * (NCONTI + NINT * KCOL * (KCOL + NCONTI)
     &         + NCONTI + NCONTI + NINT * (KCOL + 1) * (KCOL + 1
     &         + NCONTI) + NCONTI)
c-----------------------------------------------------------------------
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS 
C *** *** *** *** *** *** ***
       NFCN=0
       NJAC=0
       NSTEP=0
       NACCPT=0
       NREJCT=0
       NDEC=0
       NSOL=0
       ARRET=.FALSE.
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0  
      IF (WORK(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=WORK(1)
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            ARRET=.TRUE.
         END IF
      END IF
C -------- CHECK AND CHANGE THE TOLERANCES
      EXPM=2.0D0/3.0D0
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=0.1D0*RTOL(1)**EXPM
              ATOL(1)=RTOL(1)*QUOT
          END IF
      ELSE
          DO I=1,NPDE
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              ARRET=.TRUE.
          ELSE
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=0.1D0*RTOL(I)**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END IF
          END DO
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF (IWORK(2).EQ.0) THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF (NMAX.LE.0) THEN
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS
      IF (IWORK(3).EQ.0) THEN
         NIT=7
      ELSE
         NIT=IWORK(3)
         IF (NIT.LE.0) THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS
      IF(IWORK(4).EQ.0)THEN
         STARTN=.FALSE.
      ELSE
         STARTN=.TRUE.
      END IF
C -------- PRED   STEP SIZE CONTROL
      IF(IWORK(8).LE.1)THEN
         PRED=.TRUE.
      ELSE
         PRED=.FALSE.
      END IF
C --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION
      IF (WORK(2).EQ.0.0D0) THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2)
            ARRET=.TRUE.
         END IF
      END IF
C ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
      IF (WORK(3).EQ.0.D0) THEN
         THET=0.001D0
      ELSE
         THET=WORK(3)
         IF (THET.GE.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3)
            ARRET=.TRUE.
         END IF
      END IF
C --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1.
      TOLST=RTOL(1)
      IF (WORK(4).EQ.0.D0) THEN
         FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0))
      ELSE
         FNEWT=WORK(4)
         IF (FNEWT.LE.UROUND/TOLST) THEN
            WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4)
            ARRET=.TRUE.
         END IF
      END IF
C --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
      IF (WORK(5).EQ.0.D0) THEN
         QUOT1=1.D0
      ELSE
         QUOT1=WORK(5)
      END IF
      IF (WORK(6).EQ.0.D0) THEN
         QUOT2=1.2D0
      ELSE
         QUOT2=WORK(6)
      END IF
      IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN
         WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2
         ARRET=.TRUE.
      END IF
C -------- MAXIMAL STEP SIZE
      IF (WORK(7).EQ.0.D0) THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(7)
      END IF 
C -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(8).EQ.0.D0)THEN
         FACL=5.D0
      ELSE
         FACL=1.D0/WORK(8)
      END IF
      IF(WORK(9).EQ.0.D0)THEN
         FACR=1.D0/8.0D0
      ELSE
         FACR=1.D0/WORK(9)
      END IF
      IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN
            WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9)
            ARRET=.TRUE.
         END IF
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEZ1=21
      IEZ2=IEZ1+N
      IEZ3=IEZ2+N
      IEY0=IEZ3+N
      IESCAL=IEY0+N
      IEF1=IESCAL+N
      IEF2=IEF1+N
      IEF3=IEF2+N
      IECON=IEF3+N
      IEJAC=IECON+4*N
      IEMAS=IEJAC+NSZJAC
      IEE1=IEMAS+NSZJAC
c-----------------------------------------------------------------------  
c ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN CWORK -----
      IEE2R=1
      IECZ2=IEE2R+NSZJAC
c-----------------------------------------------------------------------  
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEE1+NSZJAC-1
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
      IEIP1=21
      IEIP2=IEIP1+N
C --------- TOTAL REQUIREMENT ---------------
      ISTORE=IEIP2+N-1
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,JAC,MAS,SOLOUT,
     &   IDID,NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,STARTN,
     &   PRED,FACL,FACR, WORK(IEZ1),WORK(IEZ2),WORK(IEZ3),WORK(IEY0),
     &   WORK(IESCAL),WORK(IEF1),WORK(IEF2),WORK(IEF3),WORK(IEJAC),
     &   WORK(IEE1),CWORK(IEE2R),WORK(IEMAS),IWORK(IEIP1),IWORK(IEIP2),
     &   WORK(IECON),NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR,
     &   CWORK(IECZ2))
      IWORK(14)=NFCN
      IWORK(15)=NJAC
      IWORK(16)=NSTEP
      IWORK(17)=NACCPT
      IWORK(18)=NREJCT
      IWORK(19)=NDEC
      IWORK(20)=NSOL
C -------- RESTORE TOLERANCES
      EXPM=1.0D0/EXPM
      IF (ITOL.EQ.0) THEN
              QUOT=ATOL(1)/RTOL(1)
              RTOL(1)=(10.0D0*RTOL(1))**EXPM
              ATOL(1)=RTOL(1)*QUOT
      ELSE
          DO I=1,NPDE
              QUOT=ATOL(I)/RTOL(I)
              RTOL(I)=(10.0D0*RTOL(I))**EXPM
              ATOL(I)=RTOL(I)*QUOT
          END DO
      END IF
C ----------- RETURN -----------
      RETURN
      END
C ***********************************************************
C
      SUBROUTINE RADCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,MAS,SOLOUT,IDID,NMAX,UROUND,SAFE,THET,FNEWT,QUOT1,
     &   QUOT2,NIT,STARTN,PRED,FACL,FACR,Z1,Z2,Z3,Y0,SCAL,F1,F2,F3,
     &   FJAC,E1,E2R,FMAS,IP1,IP2,CONT,NFCN,NJAC,NSTEP,NACCPT,
     &   NREJCT,NDEC,NSOL,RPAR,IPAR,CZ2)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR RADAU5
C     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED 
C ---------------------------------------------------------- 
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 21, 2003.
c
c-----------------------------------------------------------------------
C ---------------------------------------------------------- 
C         DECLARATIONS 
C ---------------------------------------------------------- 
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,ITOL,IDID,NMAX,NIT,NFCN,NJAC,NSTEP,NACCPT,NREJCT,
     &        NDEC,NSOL,IP1(N),IP2(N),IPAR(*)
      DOUBLE PRECISION X,XEND,HMAX,H,UROUND,SAFE,THET,FNEWT,QUOT1,QUOT2,
     &                 FACL,FACR,Y(N),Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),
     &                 F1(N),F2(N),F3(N),CONT(4*N)
c-----------------------------------------------------------------------
c the dimensions of the following arrays have been changed.
      DOUBLE PRECISION E1(*), FJAC(*), FMAS(*)
      DOUBLE COMPLEX E2R(*), CZ2(N)
c-----------------------------------------------------------------------
      DOUBLE PRECISION RTOL(*),ATOL(*),RPAR(*)
      LOGICAL STARTN,PRED
      EXTERNAL FCN,JAC,MAS,SOLOUT
c-----------------------------------------------------------------------
c variables in common expression
      INTEGER NN,NN2,NN3,NN4
      DOUBLE PRECISION XSOL,HSOL,C2M1,C1M1
      COMMON /CONRA5/NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1
c-----------------------------------------------------------------------
c local variables
      INTEGER NPDE,NINT,KCOL
c-----------------------------------------------------------------------
      INTEGER NSING,IRTRN,N2,N3,IER,NEWT
      DOUBLE PRECISION SQ6,C1,C2,C1MC2,DD1,DD2,DD3,U1,ALPH,BETA,CNO,T11,
     &                 T12,T13,T21,T22,T23,T31,TI11,TI12,TI13,TI21,TI22,
     &                 TI23,TI31,TI32,TI33,POSNEG,HMAXN,HOLD,FACCON,
     &                 CFAC,HHFAC,FAC1,ALPHN,
     &                 BETAN,XPH,C1Q,C2Q,C3Q,AK1,AK2,AK3,Z1I,Z2I,Z3I,
     &                 THETA,A1,A2,A3,DYNO,DENOM,THQ,DYNOLD,THQOLD,
     &                 DYTH,QNEWT,F1I,F2I,F3I,ERR,FAC,QUOT,HNEW,FACGUS,
     &                 HACC,ERRACC,AK,ACONT3,HOPT,QT
      LOGICAL REJECT,FIRST,CALJAC
      LOGICAL LAST
c-----------------------------------------------------------------------
c loop indices
      INTEGER I,J
c-----------------------------------------------------------------------
      NPDE = IPAR(1)
      KCOL = IPAR(2)
      NINT = IPAR(3)
c-----------------------------------------------------------------------
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
C --------- DUPLIFY N FOR COMMON BLOCK CONT -----
      NN=N
      NN2=2*N
      NN3=3*N 
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
c----------------------------------------------------------------------- 
      CALL MAS(FMAS,RPAR,IPAR)
c----------------------------------------------------------------------- 
C ---------- CONSTANTS ---------
      SQ6=DSQRT(6.D0)
      C1=(4.D0-SQ6)/10.D0
      C2=(4.D0+SQ6)/10.D0
      C1M1=C1-1.D0
      C2M1=C2-1.D0
      C1MC2=C1-C2
      DD1=-(13.D0+7.D0*SQ6)/3.D0
      DD2=(-13.D0+7.D0*SQ6)/3.D0
      DD3=-1.D0/3.D0
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0
      CNO=ALPH**2+BETA**2
      U1=1.0D0/U1
      ALPH=ALPH/CNO
      BETA=BETA/CNO
      T11=9.1232394870892942792D-02
      T12=-0.14125529502095420843D0
      T13=-3.0029194105147424492D-02
      T21=0.24171793270710701896D0
      T22=0.20412935229379993199D0
      T23=0.38294211275726193779D0
      T31=0.96604818261509293619D0
      TI11=4.3255798900631553510D0
      TI12=0.33919925181580986954D0
      TI13=0.54177053993587487119D0
      TI21=-4.1787185915519047273D0
      TI22=-0.32768282076106238708D0
      TI23=0.47662355450055045196D0
      TI31=-0.50287263494578687595D0
      TI32=2.5719269498556054292D0
      TI33=-0.59603920482822492497D0
      POSNEG=SIGN(1.D0,XEND-X)
      HMAXN=MIN(ABS(HMAX),ABS(XEND-X)) 
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAXN)
      H=SIGN(H,POSNEG)
      HOLD=H
      REJECT=.FALSE.
      FIRST=.TRUE.
      LAST=.FALSE.
      IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      END IF
      HOPT=H
      FACCON=1.D0
      CFAC=SAFE*(1+2*NIT)
      NSING=0
c-----------------------------------------------------------------------
      N2=2*N
      N3=3*N
      IF (ITOL.EQ.0) THEN
          DO I=1,N
             SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
          END DO
      ELSE
          DO J=1,N/NPDE
             DO I=1,NPDE
                SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I+(J-1)*NPDE))
             END DO
          END DO
      END IF
      HHFAC=H
      CALL FCN(N,X,Y,Y0,RPAR,IPAR)
      NFCN=NFCN+1
C --- BASIC INTEGRATION STEP  
  10  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
c----------------------------------------------------------------------- 
      CALL JAC(N,X,Y,FJAC,RPAR,IPAR)
c----------------------------------------------------------------------- 
      CALJAC=.TRUE.
  20  CONTINUE
C --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS
      FAC1=U1/H
      ALPHN=ALPH/H
      BETAN=BETA/H
c----------------------------------------------------------------------- 
      CALL DECOMR(N,NPDE,NINT,KCOL,FJAC,FMAS,FAC1,E1,IP1,IER)
c----------------------------------------------------------------------- 
      IF (IER.NE.0) GOTO 78
c----------------------------------------------------------------------- 
      CALL DECOMC(N,NPDE,NINT,KCOL,FJAC,FMAS,ALPHN,BETAN,E2R,IP2,IER)
c----------------------------------------------------------------------- 
      IF (IER.NE.0) GOTO 78
      NDEC=NDEC+1
  30  CONTINUE
      NSTEP=NSTEP+1
      IF (NSTEP.GT.NMAX) GOTO 178
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND) GOTO 177
c----------------------------------------------------------------------- 
      do 31 i = 1, npde
         scal(i) = scal(i)/hhfac
   31 continue
      do 32 i = npde*(nint*kcol+1)+1, npde*(nint*kcol+2)
         scal(i) = scal(i)/hhfac
   32 continue
      do 33 i = npde*(nint*kcol+2)+1, npde*(nint*kcol+3)
         scal(i) = scal(i)/hhfac
   33 continue
      do 34 i = n-npde+1, n
         scal(i) = scal(i)/hhfac
   34 continue
c----------------------------------------------------------------------- 
      XPH=X+H
C *** *** *** *** *** *** ***
C  STARTING VALUES FOR NEWTON ITERATION
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN) THEN
         DO I=1,N
            Z1(I)=0.D0
            Z2(I)=0.D0
            Z3(I)=0.D0
            F1(I)=0.D0
            F2(I)=0.D0
            F3(I)=0.D0
         END DO
      ELSE
         C3Q=H/HOLD
         C1Q=C1*C3Q
         C2Q=C2*C3Q
         DO I=1,N
            AK1=CONT(I+N)
            AK2=CONT(I+N2)
            AK3=CONT(I+N3)
            Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3))
            Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3))
            Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3))
            Z1(I)=Z1I
            Z2(I)=Z2I
            Z3(I)=Z3I
            F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I
            F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I
            F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
            NEWT=0
            FACCON=MAX(FACCON,UROUND)**0.8D0
            THETA=ABS(THET)
  40        CONTINUE
            IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
            DO I=1,N
               CONT(I)=Y(I)+Z1(I)
            END DO
            CALL FCN(N,X+C1*H,CONT,Z1,RPAR,IPAR)
            DO I=1,N
               CONT(I)=Y(I)+Z2(I)
            END DO
            CALL FCN(N,X+C2*H,CONT,Z2,RPAR,IPAR)
            DO I=1,N
               CONT(I)=Y(I)+Z3(I)
            END DO
            CALL FCN(N,XPH,CONT,Z3,RPAR,IPAR)
            NFCN=NFCN+3
C ---     SOLVE THE LINEAR SYSTEMS
           DO I=1,N
              A1=Z1(I)
              A2=Z2(I)
              A3=Z3(I)
              Z1(I)=TI11*A1+TI12*A2+TI13*A3
              Z2(I)=TI21*A1+TI22*A2+TI23*A3
              Z3(I)=TI31*A1+TI32*A2+TI33*A3
           END DO
c----------------------------------------------------------------------- 
        CALL SLVRAD(N,NPDE,KCOL,NINT,FMAS,FAC1,ALPHN,BETAN,E1,E2R,
     &              Z1,Z2,Z3,F1,F2,F3,CONT,IP1,IP2,CZ2)
c----------------------------------------------------------------------- 
            NSOL=NSOL+1
            NEWT=NEWT+1
            DYNO=0.D0
            DO I=1,N
               DENOM=SCAL(I)
               DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2
     &          +(Z3(I)/DENOM)**2
            END DO
            DYNO=DSQRT(DYNO/N3)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
            IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
                THQ=DYNO/DYNOLD
                IF (NEWT.EQ.2) THEN
                   THETA=THQ
                ELSE
                   THETA=SQRT(THQ*THQOLD)
                END IF
                THQOLD=THQ
                IF (THETA.LT.0.99D0) THEN
                    FACCON=THETA/(1.0D0-THETA)
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
                    IF (DYTH.GE.1.0D0) THEN
                         QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
                         HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
                         H=HHFAC*H
                         REJECT=.TRUE.
                         LAST=.FALSE.
                         IF (CALJAC) GOTO 20
                         GOTO 10
                    END IF
                ELSE
                    GOTO 78
                END IF
            END IF
            DYNOLD=MAX(DYNO,UROUND)
            DO I=1,N
               F1I=F1(I)+Z1(I)
               F2I=F2(I)+Z2(I)
               F3I=F3(I)+Z3(I)
               F1(I)=F1I
               F2(I)=F2I
               F3(I)=F3I
               Z1(I)=T11*F1I+T12*F2I+T13*F3I
               Z2(I)=T21*F1I+T22*F2I+T23*F3I
               Z3(I)=T31*F1I+    F2I
            END DO
            IF (FACCON*DYNO.GT.FNEWT) GOTO 40
C --- ERROR ESTIMATION  
      CALL ESTRAD (N,NPDE,KCOL,NINT,FMAS,H,DD1,DD2,DD3,FCN,NFCN,
     &          Y0,Y,X,E1,Z1,Z2,Z3,CONT,F1,F2,IP1,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .2<=HNEW/H<=8.
      FAC=MIN(SAFE,CFAC/(NEWT+2*NIT))
      QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC))
      HNEW=H/QUOT
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LT.1.D0) THEN
C --- STEP IS ACCEPTED  
         FIRST=.FALSE.
         NACCPT=NACCPT+1
         IF (PRED) THEN
C       --- PREDICTIVE CONTROLLER OF GUSTAFSSON
            IF (NACCPT.GT.1) THEN
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE
               FACGUS=MAX(FACR,MIN(FACL,FACGUS))
               QUOT=MAX(QUOT,FACGUS)
               HNEW=H/QUOT
            END IF
            HACC=H
            ERRACC=MAX(1.0D-2,ERR)
         END IF
         HOLD=H
         X=XPH 
         DO I=1,N
            Y(I)=Y(I)+Z3(I)  
            Z2I=Z2(I)
            Z1I=Z1(I)
            CONT(I+N)=(Z2I-Z3(I))/C2M1
            AK=(Z1I-Z2I)/C1MC2
            ACONT3=Z1I/C1
            ACONT3=(AK-ACONT3)/C2
            CONT(I+N2)=(AK-CONT(I+N))/C1M1
            CONT(I+N3)=CONT(I+N2)-ACONT3
         END DO
         IF (ITOL.EQ.0) THEN
             DO I=1,N
                SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I))
             END DO
         ELSE
             DO J=1,N/NPDE
                DO I=1,NPDE
                   SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I+(J-1)*NPDE))
                END DO
             END DO
         END IF
         CALL SOLOUT(N,X,Y,ITOL,RPAR,IPAR,IRTRN)
         IF (IRTRN.NE.0) GOTO 179
         CALJAC=.FALSE.
         IF (LAST) THEN
            H=HOPT
            IDID=1
            RETURN
         END IF
         CALL FCN(N,X,Y,Y0,RPAR,IPAR)
         NFCN=NFCN+1
         HNEW=POSNEG*MIN(ABS(HNEW),HMAXN)
         HOPT=HNEW
         HOPT=MIN(H,HNEW)
         IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H)) 
         REJECT=.FALSE.
         IF ((X+HNEW/QUOT1-XEND)*POSNEG.GE.0.D0) THEN
            H=XEND-X
            LAST=.TRUE.
         ELSE
            QT=HNEW/H 
            HHFAC=H
            IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) GOTO 30
            H=HNEW 
         END IF
         HHFAC=H
         IF (THETA.LE.THET) GOTO 20
         GOTO 10
      ELSE
C --- STEP IS REJECTED  
         REJECT=.TRUE.
         LAST=.FALSE.
         IF (FIRST) THEN
             H=H*0.1D0
             HHFAC=0.1D0
         ELSE 
             HHFAC=HNEW/H
             H=HNEW
         END IF
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         IF (CALJAC) GOTO 20
         GOTO 10
      END IF
C --- UNEXPECTED STEP-REJECTION
  78  CONTINUE
      IF (IER.NE.0) THEN
          NSING=NSING+1
          IF (NSING.GE.5) GOTO 176
      END IF
      H=H*0.5D0 
      HHFAC=0.5D0
      REJECT=.TRUE.
      LAST=.FALSE.
      IF (CALJAC) GOTO 20
      GOTO 10
C --- FAIL EXIT
 176  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER
      IDID=-4
      RETURN
 177  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
      IDID=-3
      RETURN
 178  CONTINUE
      WRITE(6,979)X   
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED' 
      IDID=-2
      RETURN
C --- EXIT CAUSED BY SOLOUT
 179  CONTINUE
 979  FORMAT(' EXIT OF RADAU5 AT X=',E18.4) 
      IDID=2
      RETURN
      END
      subroutine radfcn(neq, t, y, fr, rpar, ipar) 

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the vector at the right side of the 
c       DAEs, which is required by RADAU5.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 21, 2003.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of DAEs.
c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients,
c                               including the one for radau_kcol and the
c                               one for radau_kcol+1
c
        double precision        rpar(*)
c                               rpar is the BACOLR floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLR integer work array.
c
c       Output:
        double precision        fr(neq)
c                               fr is the vector at the right side of 
c                               the DAEs.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt1
        parameter              (incpt1 =  4)
c                               ipar(incpt1) = ncpts1.
c
        integer                 ineq1
        parameter              (ineq1  =  5)
c                               ipar(ineq1) = neq1.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:    
        integer                 ixcol1
        parameter              (ixcol1 = 22)
c                               rpar(ipar(ixcol1)) stores the
c                               collocation points when using
c                               radau_kcol.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)  
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi1
        parameter              (ibasi1 = 31)
c                               rpar(ipar(ibasi1)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol.
c                               rpar(ipar(ibasi1)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 ixcol2
        parameter              (ixcol2 = 41)   
c                               rpar(ipar(ixcol2)) stores the
c                               collocation points when using
c                               radau_kcol+1.
c
        integer                 ibasi2
        parameter              (ibasi2 = 48)
c                               rpar(ipar(ibasi2)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol+1.
c                               rpar(ipar(ibasi2)) contains
c                               a three dimensional array A of size
c                               (kcol+1+nconti,3,ncpts). A(k,j,i)
c                               contains the values of the (j-1)st
c                               derivative (j=1,2,3) of the k-th
c                               non-zero basis function (k=1,...,
c                               kcol+1+nconti) at the i-th collocation
c                               point.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts1
        integer                 ncpts2
        integer                 neq1
        integer                 neq2
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               calfcn
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts1 = ipar(incpt1)
      ncpts2 = ncpts1 + nint
      neq1   = ipar(ineq1)
      neq2   = neq - neq1

c     Calculate the right side of the DAEs for radau_kcol.   
      call calfcn(npde, kcol, nint, ncpts1, neq1, rpar(ipar(ixcol1)),
     &            rpar(ipar(ibasi1)), t, y, rpar(ipar(iwkrj)), fr)

c     Calculate the right side of the DAEs for radau_kcol+1.   
      call calfcn(npde, kcol+1, nint, ncpts2, neq2, rpar(ipar(ixcol2)),
     &            rpar(ipar(ibasi2)), t, y(neq1+1), rpar(ipar(iwkrj)),
     &            fr(neq1+1))

      return
      end
      subroutine radjac(neq, t, y, dfdy, rpar, ipar)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the Jacobian matrix for the right side
c       of the DAEs, which is required by RADAU5.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, May 22, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 neq
c                               neq is the number of DAEs.
c
        double precision        t
c                               t is the current time.
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients,
c                               including the one for radau_kcol and the
c                               one for radau_kcol+1
c
        double precision        rpar(*)
c                               rpar is the BACOLR floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLR integer work array.
c
c       Output:
        double precision        dfdy(*)
c                               dfdy is the ABD Jacobian matrix for the
c                               right side of the DAE.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 incpt1
        parameter              (incpt1 =  4)
c                               ipar(incpt1) = ncpts1.
c
        integer                 ineq1
        parameter              (ineq1  =  5)
c                               ipar(ineq1) = neq1.
c
c----------------------------------------------------------------------- 
c Indirect pointers into the RPAR floating point work array:     
        integer                 ixcol1
        parameter              (ixcol1 = 22)
c                               rpar(ipar(ixcol1)) stores the
c                               collocation points when using
c                               radau_kcol.
c
        integer                 iabtp1
        parameter              (iabtp1 = 26)   
c                               rpar(ipar(iabtp1)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_kcol.
c
        integer                 iabbt1
        parameter              (iabbt1 = 28)  
c                               rpar(ipar(iabbt1)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_kcol.
c
        integer                 iwkrj
        parameter              (iwkrj  = 30)
c                               rpar(ipar(iwkrj)) stores an additional
c                               work array required by res and jac.
c
        integer                 ibasi1
        parameter              (ibasi1 = 31)
c                               rpar(ipar(ibasi1)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol.
c                               rpar(ipar(ibasi1)) contains
c                               a three dimensional array A of size
c                               (kcol+nconti,3,ncpts). A(k,j,i) contains
c                               the values of the (j-1)st derivative
c                               (j=1,2,3) of the k-th non-zero basis
c                               function (k=1,...,kcol+nconti) at the
c                               i-th collocation point.
c
        integer                 ixcol2
        parameter              (ixcol2 = 41)
c                               rpar(ipar(ixcol2)) stores the
c                               collocation points when using
c                               radau_kcol+1.
c
        integer                 iabtp2
        parameter              (iabtp2 = 45)  
c                               rpar(ipar(iabtp2)) stores the top block
c                               of the ABD collocation matrices when
c                               using radau_kcol+1.
c
        integer                 iabbt2
        parameter              (iabbt2 = 47)
c                               rpar(ipar(iabbt2)) stores the bottom
c                               block of the ABD collocation matrices
c                               when using radau_{kcol+1}.
c
        integer                 ibasi2
        parameter              (ibasi2 = 48)
c                               rpar(ipar(ibasi2)) stores the basis
c                               function values at the collocation
c                               points when using radau_kcol+1.
c                               rpar(ipar(ibasi2)) contains
c                               a three dimensional array A of size
c                               (kcol+1+nconti,3,ncpts). A(k,j,i)
c                               contains the values of the (j-1)st
c                               derivative (j=1,2,3) of the k-th
c                               non-zero basis function (k=1,...,
c                               kcol+1+nconti) at the i-th collocation
c                               point.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde
        integer                 kcol
        integer                 nint
        integer                 ncpts1
        integer                 ncpts2
        integer                 neq1
        integer                 neq2
        integer                 njcp1   
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               caljac
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      ncpts1 = ipar(incpt1)
      ncpts2 = ncpts1 + nint
      neq1   = ipar(ineq1)
      neq2   = neq - neq1
      njcp1  = npde * npde * (2 * nconti + nint * kcol *
     &         (kcol + nconti)) + 1

c     Calculate the Jacobian matrix at the right side of the DAEs for
c     radau_kcol.
      call caljac(npde, kcol, nint, ncpts1, neq1, rpar(ipar(ixcol1)),
     &            rpar(ipar(ibasi1)), rpar(ipar(iabtp1)),
     &            rpar(ipar(iabbt1)), t, y, rpar(ipar(iwkrj)), dfdy)
 
c     Calculate the Jacobian matrix at the right side of the DAEs for
c     radau_kcol+1.
      call caljac(npde, kcol+1, nint, ncpts2, neq2, rpar(ipar(ixcol2)),
     &            rpar(ipar(ibasi2)), rpar(ipar(iabtp2)),
     &            rpar(ipar(iabbt2)), t, y(neq1+1), rpar(ipar(iwkrj)),
     &            dfdy(njcp1))
 
      return
      end
      subroutine radmas(am, rpar, ipar)

c-----------------------------------------------------------------------
c Purpose:
c       This routine calculates the mass-matrix which is required by
c       the implicit Runge-Kutta DAE solver RADAU5.
c
c-----------------------------------------------------------------------  
c
c Last modified by Rong Wang, May 17, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:     
        double precision        rpar(*)
c                               rpar is the BACOLR floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLR integer work array.
c
c       Output:
        double precision        am(*)
c                               am is the ABD mass-matrix.
c
c-----------------------------------------------------------------------  
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array:
        integer                 iabbk1
        parameter              (iabbk1 = 27)  
c                               rpar(ipar(iabbk1)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               radau_kcol.
c
        integer                 iabbk2
        parameter              (iabbk2 = 46) 
c                               rpar(ipar(iabbk2)) stores the nint
c                               blocks in the middle of the ABD
c                               collocation matrices when using
c                               radau_kcol+1.
c
c-----------------------------------------------------------------------
c Local variables:
        integer                 npde 
        integer                 kcol
        integer                 nint  
        integer                 namp1   
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               calmas
c
c-----------------------------------------------------------------------

      npde   = ipar(inpde)
      kcol   = ipar(ikcol)
      nint   = ipar(inint)
      
      namp1  = npde * npde * (2 * nconti + nint * kcol *
     &         (kcol + nconti)) + 1
      
c     Calculate the mass-matrix for radau_kcol.
      call calmas(npde, kcol, nint, rpar(ipar(iabbk1)), am)

c     Calculate the mass-matrix for radau_kcol+1.
      call calmas(npde, kcol+1, nint, rpar(ipar(iabbk2)), am(namp1))

      return
      end
      subroutine reinit(npde, kcol, kold, nint, ninold, ncpts, neq,
     &                  neqold, x, xold, yold, work, lw, ipivot, h, 
     &                  xbs, xcol, fbasis, y, abdblk, icflag)

c-----------------------------------------------------------------------
c Purpose:
c       This routine performs the initialization tasks after remeshing: 
c
c               calculating the mesh step size sequence,
c               generating the piecewise polynomial space breakpoint
c               sequence,
c               calculating the collocation point sequence,
c               calculating the B-spline basis functions,
c               constructing abdblk of the collocation matrices and
c               calculating the bspline coefficients at the last step.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 15, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points. 
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       Input:
        integer                 npde
c                               npde is the number of components in
c                               the system of PDEs. npde > 0.
c
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval after
c                               remeshing. 
c
        integer                 kold
c                               kold is the number of collocation points
c                               to be used in each subinterval before
c                               remeshing. 
c
        integer                 nint
c                               nint is the number of subintervals after
c                               remeshing.
c
        integer                 ninold
c                               ninold is the number of subintervals 
c                               before remeshing.
c
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number 
c                               of collocation points.
c
        integer                 neq
c                               neq=npde*(kcol*nint+nconti) is
c                               the number of bspline 
c                               coefficients after remeshing.
c
        integer                 neqold
c                               neqold=npde*(kold*ninold+nconti) is
c                               the number of bspline 
c                               coefficients before remeshing.
c
        double precision        x(nint+1)
c                               x is the spatial mesh after remeshing.
c
        double precision        xold(ninold+1)
c                               xold is the spatial mesh before
c                               remeshing.
c
        double precision        yold(neqold)
c                               yold is the vector of bspline
c                               coefficients at the last time step.
c
        double precision        h(nint)
c                               h is the mesh step size sequence.
c
        integer                 lw
c                               lw is the size of the work storage 
c                               array and must satisfy:
c                               lw >= 2*npde*npde*nconti+
c                                     npde*npde*kcol*(kcol+nconti)*nint
c                                     +(kold+nconti)+kold*(ninold+1)
c                                     +2*nconti
c                               Since nint >= ninold/2 and kcol >=
c                               kold+1, it implies that lw >= 3*neqold.
c
c       Work Storage:
        double precision        work(lw)
c                               work is a floating point work storage
c                               array of size lw.
c
        integer                 ipivot(neq)
c                               pivoting information from the 
c                               factorization of the temporary matrix.
c
c       Output:
        double precision        xcol(ncpts)
c                               The sequence of collocation points on
c                               the interval [x_a, x_b].
c
        double precision        xbs(ncpts+kcol+nconti)
c                               The breakpoint sequence.
c
        double precision        fbasis(kcol+nconti, 3, ncpts)
c                               Basis function values at the collocation
c                               points. 
c
        double precision        y(neq)
c                               y is the vector of bspline coefficients
c                               at the last time step after remeshing. 
c
        double precision        abdblk(npde*npde*nint*kcol
     &                                 *(kcol+nconti))
c                               The nint blocks in the middle of
c                               the matrix A.
c
        integer                 icflag
c                               This is the status flag from COLROW
c                               which is called by crdcmp.
c                               icflag =  0, indicates non-singularity.
c                               icflag = -1, indicates singularity.
c                               icflag =  1, indicates invalid input.
c-----------------------------------------------------------------------
c Local Variables:  
        integer                 ileft
c                               breakpoint information.
c
        integer                 nels
c                               the number of elements in one 
c                               collocation block of work.
c 
c       Pointers into the floating point work array:
        integer                 iabdtp
c                               work(iabdtp) contains a copy of the top
c                               block which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbk
c                               work(iabdbk) contains a copy of abdblk
c                               which is required since crdcmp
c                               overwrites the input collocation matrix.
c
        integer                 iabdbt
c                               work(iabdbt) contains a copy of the 
c                               bottom block which is required since 
c                               crdcmp overwrites the input collocation 
c                               matrix.
c
        integer                 ivwork
c                               work(ivwork) is the work storage
c                               required by values.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 l
        integer                 m
        integer                 ii
        integer                 ll
        integer                 mm
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               colpnt
c                               crdcmp
c                               crslve
c                               values
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c       double precision:
c                               dcopy
c
c-----------------------------------------------------------------------

c     Generate the piecewise polynomial space breakpoint sequence,
c     and calculates the collocation point sequences.
      call colpnt(kcol, nint, ncpts, x, h, work, xcol, xbs)

c-----------------------------------------------------------------------
      nels = npde*npde*kcol*(kcol+nconti)

c     Set the pointers into the floating point work array.
      iabdtp = 1
      iabdbk = iabdtp + npde*npde*nconti
      iabdbt = iabdbk + nint*nels
      ivwork = iabdbt + npde*npde*nconti

c-----------------------------------------------------------------------
c     Call values to calculate the values at xcol and at the last 
c     time step. Then save in y.
      call values(kold, xcol, ninold, xold, npde, ncpts, 0, y, yold,
     &            work(ivwork))

c-----------------------------------------------------------------------
c     Initialize abdblk, the top block and the bottom block to zero.
      do 10 i = 1, npde * npde * nconti
         work(iabdtp-1+i) = zero
         work(iabdbt-1+i) = zero
   10 continue
      do 20 i = 1, nint*nels
         abdblk(i) = zero
   20 continue

c-----------------------------------------------------------------------
c     Bsplvd is called to compute the components of fbasis(k,i,j) 
c     associated the first collocation point. Now ileft = kcol + nconti.
      call bsplvd(xbs,kcol+nconti,xcol(1),kcol+nconti,fbasis(1,1,1),3)

c     Makeing use of the fact that only the first bspline has a nonzero
c     value at the left end point, set up the top block in work.
      do 30 i = 1, npde 
         ii = (i-1) * npde + i
         work(iabdtp-1+ii) = fbasis(1,1,1)
   30 continue

c-----------------------------------------------------------------------
c     The nint blocks at the middle of the matrix will now be set up.
      do 70 i = 1, nint

c     Make use the fact that there are kcol collocation points in each
c     subinterval to find the value of ileft.
         ileft = kcol + nconti + (i - 1) * kcol

         do 60 j = 1, kcol

c     ii is the position in xcol of the j-th collocation point of the
c     i-th subinterval.
            ii = (i-1) * kcol + 1 + j

c     compute information for ii-th collocation point.
            call bsplvd(xbs,kcol+nconti,xcol(ii),ileft,fbasis(1,1,ii),3)

            do 50 l = 1, kcol + nconti

c     generate the subblock in abdblk corresponding to the ii-th
c     collocation point.
c
               ll = (l-1)*npde*npde*kcol + (i-1)*nels + (j-1)*npde
               do 40 m = 1, npde
                  mm = ll + (m-1)*npde*kcol + m
                  abdblk(mm) = fbasis(l,1,ii)
   40          continue
   50       continue
   60    continue
   70 continue

c-----------------------------------------------------------------------
c     Now, set up the bottom block, using the fact that only the
c     last bspline basis function is non-zero at the right end point.
c     Simultaneously, set up the corresponding part of the right hand
c     side.
c
      call bsplvd(xbs,kcol+nconti,xcol(ncpts),ncpts,
     &            fbasis(1,1,ncpts),3)
      do 80 i = 1, npde
         ii = ((i-1)+npde)*npde + i  
         work(iabdbt-1+ii) = fbasis(kcol+nconti,1,ncpts)
   80 continue

c-----------------------------------------------------------------------   
c     Copy the middle of the collocation matrix into temporary storage.
      call dcopy(nels*nint,abdblk,1,work(iabdbk),1)

c-----------------------------------------------------------------------
c     Generate the vector y.

c     LU decompose the matrix.

      call crdcmp(neq,work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,ipivot,
     &            icflag)

      if (icflag .ne. 0) go to 999

c     Solve the linear system. This gives the basis function 
c     coefficients for the initial conditions, i.e. y(t0).
      call crslve(work(iabdtp),npde,2*npde,work(iabdbk),kcol*npde,
     &            (kcol+nconti)*npde,nint,work(iabdbt),npde,
     &            ipivot,y,0)
      if (icflag .ne. 0) go to 999

  999 return
      end
      subroutine remesh(istart, icount, nintmx, ninold, errrat,
     &                  errint, irshfg, nint, kcol, x, work)

c-----------------------------------------------------------------------
c Purpose:
c       This routine generates a new mesh by equidistributing the error
c       in each subinterval.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, July 11, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        double precision        point5
        parameter              (point5 = 0.5d0)
c
        double precision        one
        parameter              (one    = 1.0d0)
c
        double precision        two
        parameter              (two    = 2.0d0)
c
        double precision        saffac
        parameter              (saffac = 0.2d0)
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 istart
c                               istart is a flag to begin the code.
c                               istart = 0, it is the initial step;
c                                      = 1, it is not the initial step.
c
        integer                 icount
c                               icount is the number of remeshing times
c                               at the current step.
c
        integer                 nintmx
c                               the maximal number of subintervals that
c                               the user requires.
c
        integer                 ninold
c                               ninold is the number of subintervals
c                               before remeshing.
c
        double precision        errrat
c                               errrat is the value of the largest
c                               component of rpar(ipar(iercom)).
c
        double precision        errint(ninold)
c                               errint is the error estimate at
c                               each subintervals.
c
c       Output:
        integer                 irshfg
c                               irshfg is a flag for redefining all the
c                               pointers.
c                               irshfg = 0, initial call or continuation
c                                           calls;
c                                      = 1, remesh.
c
c       In-output: 
        integer                 kcol
c                               kcol is the number of collocation points
c                               to be used in each subinterval.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        integer                 nint
c                               nint is the number of subintervals
c                               defined by the spatial mesh x.
c                               ninmx >= nint >= 1.
c                               As input, it is the value before
c                               remeshing; as output, it is the value
c                               after remeshing.
c
        double precision        x(nintmx+1)
c                               x is the spatial mesh. As input, it is
c                               the value before remeshing; as output,
c                               it is the value after remeshing. 
c
c       Work storage:
        double precision        work(2*ninold+1)
c
c-----------------------------------------------------------------------
c Local Variables:
        double precision        aerr
        double precision        berr
c
c       Pointers into the floating point work array:
        integer                 ierror
c                               work(ierror-1+i) is the L2-norm error
c                               estimate at the first i subintervals.
c
        integer                 ixold
c                               work(ixold) contains a copy of mesh
c                               points before remeshing.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
c
c-----------------------------------------------------------------------
c Functions used:
c                               dble
c                               int
c
c-----------------------------------------------------------------------
      
c     Set the pointers into the floating point work array.
      ierror = 1
      ixold = ierror + ninold

c-----------------------------------------------------------------------
c     Update icount, irshfg and nint.
      icount = icount + 1

c     If this is the first remesh at the current step which is not the
c     initial step.
      if ((icount .eq. 1) .and. (istart .eq. 1)) then
         irshfg = 1
         goto 20
      endif

c     Update errrat.
      errrat = (errrat/saffac) ** (one/dble(kcol+2))

c     Set the upper bound and lower bound of the ratio of nint over
c     ninold.
      if (errrat .gt. two) then
         errrat = two
      else
         if (errrat .lt. point5) then
            errrat = point5
         endif
      endif

      nint = int(ninold * errrat)

c     The code does not allow nint = ninold.
      if (nint .eq. ninold) then
         nint = nint + 1
      endif

   20 continue

c-----------------------------------------------------------------------
c     Update work(ixold) to be the mesh before remeshing.
      do 30 i = 1, ninold + 1
         work(ixold-1+i) = x(i)
   30 continue

c-----------------------------------------------------------------------
c     Store work(i) to be the sum of the error at the first i
c     subintervals.
      work(ierror) = errint(1)
      do 40 i = ierror-1+2, ninold
         work(i) = errint(i) + work(i-1)
   40 continue

c     Let aerr to be the mean value of errint(i).
      aerr = work(ninold)/dble(nint)

c     Equidistribute the mesh points.
      berr = aerr
      j = 1

      do 60 i = 2, nint
   50    continue
         if (berr .gt. work(j)) then
            j = j + 1
            goto 50
         else
            if (j .eq. 1) then
               x(i) = work(ixold) + (work(ixold-1+2) - work(ixold))
     &                * berr/work(1)
            else
               x(i) = work(ixold-1+j) + (work(ixold-1+j+1) -
     &                work(ixold-1+j)) * (berr - work(j-1))/errint(j)
            endif
         endif
         berr = berr + aerr
   60 continue

      x(1) = work(ixold)
      x(nint+1) = work(ixold-1+ninold+1)

      return
      end
C ***********************************************************
C
      SUBROUTINE SLVRAD(N,NPDE,KCOL,NINT,FMAS,FAC1,ALPHN,BETAN,
     &                  E1,E2R,Z1,Z2,Z3,F1,F2,F3,CONT,IP1,IP2,CZ2)
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, April 23, 2003.
c
c-----------------------------------------------------------------------
c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c
        double precision        zero
        parameter              (zero = 0.0D0)
c
c-----------------------------------------------------------------------
c subroutine parameters
      INTEGER N,NPDE,KCOL,NINT,IP1(N),IP2(N)
      DOUBLE PRECISION FMAS(*),FAC1,ALPHN,BETAN,E1(*),
     &                 Z1(N),Z2(N),Z3(N),F1(N),F2(N),F3(N),CONT(N)
      DOUBLE COMPLEX E2R(*),CZ2(N)
c-----------------------------------------------------------------------
c local variables
      INTEGER I,J,K,M,L,III,KKK,MMM,NNN,II,KK,MM,NEQ1,KCOLTM,NPDTP1,
     &        NPDBK1,NPDBT1,NPDTP2,NPDBK2,NPDBT2
      DOUBLE PRECISION S1,S2,S3,BB

c-----------------------------------------------------------------------

      neq1 = npde*(kcol*nint+nconti)
      nnn  = npde*npde*(2*nconti+kcol*(kcol+nconti)*nint)

      do 70 l = 1, 2
         kcoltm = kcol + l - 1
         iii = npde*npde*kcoltm
         mmm = (l-1)*neq1
         kkk = (l-1)*nnn
         do 40 i = 1, nint
            do 30 k = 1, kcoltm
               do 20 m = 1, npde
                  s1 = zero
                  s2 = zero
                  s3 = zero
                  ii = mmm+npde+(i-1)*npde*kcoltm+(k-1)*npde+m
                  do 10 j = 1, kcoltm + nconti
                     kk = kkk+1+(i-1)*iii*(kcoltm+nconti)+(j-1)*iii
     &                    +(k-1)*npde+npde*npde*nconti
                     mm = mmm+(i-1)*kcoltm*npde+(j-1)*npde+m
                     bb = fmas(kk)
                     s1 = s1 - bb * f1(mm)
                     s2 = s2 - bb * f2(mm)
                     s3 = s3 - bb * f3(mm)
   10             continue
                  z1(ii)   = z1(ii) + s1 * fac1
                  z2(ii)   = z2(ii) + s2 * alphn - s3 * betan
                  cont(ii) = z3(ii) + s3 * alphn + s2 * betan
   20          continue
   30       continue
   40    continue
         do 50 i = 1, npde
            z1(i+mmm) = - z1(i+mmm) * fac1
            cont(i+mmm) = z3(i+mmm)
   50    continue
         do 60 i = 1, npde
            ii = i+npde*(kcoltm*nint+1)
            z1(ii+mmm) = - z1(ii+mmm) * fac1
            cont(ii+mmm) = z3(ii+mmm)
   60    continue
   70 continue

      do 80 i = 1, n
         cz2(i) = cmplx(z2(i),cont(i))
   80 continue

      do 90 i = 1, npde
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
   90 continue

      do 100 i = neq1-npde+1, neq1
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
  100 continue

      do 110 i = neq1+1, neq1+npde
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
  110 continue

      do 120 i = n-npde+1, n
         cz2(i) = - cmplx(alphn,betan) * cz2(i)
  120 continue

      npdtp1 = 1
      npdbk1 = npdtp1 + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)

      npdtp2 = npdbt1 + npde * npde * nconti
      npdbk2 = npdtp2 + npde * npde * nconti
      npdbt2 = npdbk2 + npde * npde * nint * (kcol + 1) * (kcol + 1 +
     &                  nconti)


      call crslve(e1(npdtp1), npde, 2*npde, e1(npdbk1), kcol*npde,
     &           (kcol+nconti)*npde, nint, e1(npdbt1), npde, ip1, z1, 0)
     
      call crslve(e1(npdtp2), npde, 2*npde, e1(npdbk2), (kcol+1)*npde,
     &            (kcol+1+nconti)*npde, nint, e1(npdbt2), npde,
     &            ip1(neq1+1), z1(neq1+1), 0)

      call ccrslv(e2r(npdtp1), npde, 2*npde, e2r(npdbk1), kcol*npde,
     &            (kcol+nconti)*npde, nint, e2r(npdbt1), npde, ip2, cz2)

      call ccrslv(e2r(npdtp2), npde, 2*npde, e2r(npdbk2), (kcol+1)*npde,
     &            (kcol+1+nconti)*npde, nint, e2r(npdbt2), npde,
     &            ip2(neq1+1), cz2(neq1+1))

      do 130 i = 1, n
         z2(i) = dble(cz2(i))
         z3(i) = dimag(cz2(i))
  130 continue

      return
      end

      subroutine solout(neq, t0, y, itol, rpar, ipar, irtrn)

c-----------------------------------------------------------------------
c Purpose:
c       This routine is called at each successful time step. It decides
c       whether a remeshing is necessary or not.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Jun 5, 2003.
c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 neq
c                               neq=npde*ncpts is the number of bspline
c                               coefficients (or DAEs).
c
        double precision        t0
c                               t0 is the current time, t0 <= tout.
c
        double precision        y(neq)
c                               y is the vector of bspline
c                               coefficients at the current time.
c
        integer                 itol
c                               itol = 0, both rtol and atol are 
c                                         scalars;
c                                    = 1, both rtol and atol are
c                                         vectors.
c
        double precision        rpar(*)
c                               rpar is the BACOLR floating point work
c                               array.
c
        integer                 ipar(*)
c                               ipar is the BACOLR integer work array.
c
c       output:                    
        integer                 irtrn
c                               irtrn is a status flag for remesh.
c                               irtrn = 0, indicates no need remeshing.
c                               irtrn = 1, indicates need remeshing.
c
c-----------------------------------------------------------------------
c Direct pointers into the IPAR integer work array:
        integer                 inpde
        parameter              (inpde  =  1)
c                               ipar(inpde) = npde
c
        integer                 ikcol
        parameter              (ikcol  =  2)
c                               ipar(ikcol) = kcol.
c
        integer                 inint
        parameter              (inint  =  3)
c                               ipar(inint) = nint.
c
        integer                 ineq1
        parameter              (ineq1  =  5)
c                               ipar(ineq1) = neq1.
c
        integer                 istalr
        parameter              (istalr = 12)
c                               ipar(istalr) is the number of accepted
c                               steps after the last successful
c                               remeshing.
c
        integer                 irshfg
        parameter              (irshfg = 14)
c                               ipar(irshfg) is a flag for redefining
c                               all the pointers.
c                               ipar(irshfg) = 0, the initial step or
c                                                 any step not needing
c                                                 remesh;
c                                            = 1, a step needing remesh.
c
        integer                 icount
        parameter              (icount = 15)
c                               ipar(icount) is the number of remeshing
c                               times at the current step.
c
        integer                 istart
        parameter              (istart = 16)
c                               ipar(istart) is a flag to begin the
c                               code.
c                               ipar(istart) = 0, the initial step;
c                                            = 1, not the initial step.
c
        integer                 iatol
        parameter              (iatol  = 32)
c                               rpar(ipar(iatol)) = atol.
c
        integer                 irtol
        parameter              (irtol  = 33)
c                               rpar(ipar(irtol)) = rtol.
c
c-----------------------------------------------------------------------
c Direct pointers into the RPAR floating point work array:
        integer                 ierrat
        parameter              (ierrat =  3)
c                               rpar(ierrat) = the value of the largest
c                               component of rpar(ipar(iercom)).
c
        integer                 it0
        parameter              (it0    =  4)
c                               rpar(it0)    = t0 at the last accepted
c                               time step.
c
c-----------------------------------------------------------------------
c Indirect pointers into the RPAR floating point work array: 
        integer                 ixbs1
        parameter              (ixbs1  = 23)
c                               rpar(ipar(ixbs1)) stores the breakpoint
c                               sequence when using radau_{kcol}.
c
        integer                 iexcol
        parameter              (iexcol = 34)
c                               rpar(ipar(iexcol)) stores the
c                               collocation points which are used for
c                               error estimate.
c
        integer                 iewts
        parameter              (iewts  = 35)
c                               rpar(ipar(iewts)) stores the gaussian
c                               weights which are used for error
c                               estimate.
c
        integer                 iebas1
        parameter              (iebas1 = 36)
c                               rpar(ipar(iebas1)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol}.
c
        integer                 iebas2
        parameter              (iebas2 = 37)
c                               rpar(ipar(iebas2)) stores the values
c                               of the nonzero basis functions at
c                               rpar(ipar(iexcol)) when using
c                               dassl_{kcol+1}.
c
        integer                 iercom
        parameter              (iercom = 38)
c                               rpar(ipar(iercom)) stores the error
c                               estimate for each component.
c
        integer                 ierint
        parameter              (ierint = 39)
c                               rpar(ipar(ierint)) stores the error
c                               estimate at each subinterval.
c
        integer                 iework
        parameter              (iework = 40)
c                               rpar(ipar(iework)) stores the floating
c                               point work array for errest.
c
        integer                 ixbs2
        parameter              (ixbs2  = 42)
c                               rpar(ipar(ixbs2)) stores the breakpoint
c                               sequence when using radau_{kcol+1}.
c
        integer                 iypre
        parameter              (iypre  = 53)
c                               rpar(ipar(iypre)) stores the values of
c                               rpar(ipar(iy2)) at the previous 6 steps.
c                               It is required for a hot restart after
c                               remeshing.
c
c-----------------------------------------------------------------------
c Local variables:
c
        integer                 neq2
c                               neq2=neq1+npde*nint is the number of
c                               bspline coefficients (or DAEs) when
c                               using radau_{kcol+1}.
c
        integer                 necpts
c                               necpts=(kcol+3)*nint is the total number
c                               of collocation points used for
c                               error estimate.
c
        integer                 lenerr
c                               lenerr is the size of the floating point
c                               work array used by ERREST.
c                               lenerr>=2*npde*necpts+npde*nint.
c
c-----------------------------------------------------------------------
c BLAS Subroutines Called:
c                               dcopy
c
c-----------------------------------------------------------------------
c Subroutines Called:
c                               errest
c
c-----------------------------------------------------------------------

      neq2 = neq - ipar(ineq1) 
      necpts = (ipar(ikcol) + 3) * ipar(inint)
      lenerr = (2 * necpts + ipar(inint)) * ipar(inpde)

      call errest(ipar(ikcol), ipar(inint), ipar(inpde), ipar(ineq1),
     &            neq2, necpts, ipar(icount), rpar(ipar(iexcol)),
     &            rpar(ipar(iewts)), rpar(ipar(ixbs1)),
     &            rpar(ipar(ixbs2)), y, y(ipar(ineq1)+1), ipar(istart),
     &            itol, rpar(ipar(iatol)), rpar(ipar(irtol)), lenerr,
     &            rpar(ipar(iework)), rpar(ipar(iebas1)),
     &            rpar(ipar(iebas2)), rpar(ierrat), rpar(ipar(ierint)),
     &            rpar(ipar(iercom)), irtrn)

      if (irtrn .ne. 0) goto 100

c     The current step is accepted.
      if (ipar(icount) .ne. 0) then
         ipar(istalr) = 1
         ipar(irshfg) = 0
      else
         ipar(istalr) = ipar(istalr) + 1
      endif

c     Update the backup information. 
      call dcopy(neq2, y(ipar(ineq1)+1), 1,  rpar(ipar(iypre)), 1)
      ipar(icount) = 0
      ipar(istart) = 1
      rpar(it0)    = t0

  100 continue

      return
      end

      subroutine values(kcol, xsol, nint, x, npde, npts, nderiv, 
     &                  usol, y, work)

c-----------------------------------------------------------------------
c Purpose:
c     This routine computes the solution u and the first nderv
c     derivatives of u at the npts points xsol. It then returns the 
c     values in the array usol.

c Constants:
        integer                 nconti
        parameter              (nconti = 2)
c                               nconti continuity conditions are imposed
c                               at the internal mesh points.
c
        double precision        zero
        parameter              (zero = 0.0D0)
c-----------------------------------------------------------------------

c
c-----------------------------------------------------------------------
c Subroutine Parameters:
c       input
        integer                 kcol
c       kcol is the number of collocation points to be used in
c       each subinterval.
c
        integer                 npts
c       npts is the number of points in the x vector.
c
        double precision        xsol(npts)
c       xsol is an arbitrary set of spatial points at which the solution
c       and the first nderv derivative values are to be calculated.
c
        integer                 nint
c       nint >= 1 is the number of subintervals defined by the spatial mesh x.
c
        double precision        x(nint+1)
c       x is the spatial mesh which divides the interval [x_a,x_b] into
c       x_a = x(1) < x(2) < x(3) < ... < x(nint+1) = x_b.
c
        integer                 npde
c       npde is the number of components in the system of PDEs. npde > 0.
c
        integer                 nderiv
c       nderiv is the number of derivatives of the solution which are
c       to be calculated.
c
        double precision        y(npde*(nint*kcol+nconti))
c       y is the vector of bspline coefficients at the final time step.
c
c       output:
        double precision        usol(npde, npts, nderiv+1)
c       usol is the solution and the spatial derivatives up to the
c       nderiv-th derivative at the given points and at the final time step.
c
c       Work Storage:
        double precision        work((kcol+nconti)*(nderiv+1)
     *                                 +kcol*(nint+1)+2*nconti)
c       work is a floating point work storage array.
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c Local Variables:
        integer                 ncpts
c                               ncpts=(kcol*nint+nconti) is the number
c                               of collocation points.
c
        integer                 ileft
c                               breakpoint information.
c
        integer                 mflag
c                               mflag is required by subroutine
c                               interv.
c
        integer                 ilo
c                               ilo is required by subroutine
c                               interv.
c       Pointers into the floating point work array:
        integer                 ixbs
c                               work(ixbs) contains the breakpoint
c                               sequence.
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
        integer                 j
        integer                 k
        integer                 m
        integer                 ii
        integer                 mj
        integer                 mm
c-----------------------------------------------------------------------
c Subroutines Called:
c                               bsplvd
c                               interv
c
c-----------------------------------------------------------------------
    
c     set up the value for ileft, mflag and ncpts.
      ileft = 0
      mflag = -2
      ncpts = nint * kcol + nconti

c     set the pointer into the floating point work array
      ixbs  = (kcol+nconti)*(nderiv+1) + 1

c     Store the piecewise polynomial space breakpoint sequence in 
c     work(ixbs).
c
      do 10 i = 1, kcol + nconti
         work(ixbs-1+i) = x(1)
         work(ixbs-1+i+ncpts) = x(nint+1)
   10 continue
      do 30 i = 2, nint
         ii = (i-2) * kcol + kcol + nconti
         do 20 k = 1, kcol
            work(ixbs-1+ii+k) = x(i)
   20 continue
   30 continue

      do 70 i = 1, npts
c
c     interv is called to compute ileft. bsplvd is called to compute
c     the values of the basis function at the required point.
         call interv(work(ixbs), ncpts, xsol(i), ileft, mflag, ilo)
         call bsplvd(work(ixbs),kcol+nconti,xsol(i),ileft,work,
     &               nderiv+1)
         ii = ileft - kcol - nconti
         do 60 j = 1, nderiv + 1
            do 50 k = 1, npde
               usol(k,i,j) = zero
               do 40 m = 1, kcol + nconti
                  mm = (m + ii - 1) * npde  
                  mj = (j - 1) * (kcol + nconti) + m
                  usol(k,i,j) = usol(k,i,j) + y(mm+k) * work(mj)
   40          continue
   50       continue
   60    continue
   70 continue
      return
      end
