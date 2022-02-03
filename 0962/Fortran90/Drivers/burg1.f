c       Problem definition for One Layer Burgers Equation.
c       PDE: u_t = eps*u_xx - u*u_x, with initial and boundary conditions
c       defined from the exact solution: 
c       u = 0.5d0 - 0.5d0 * tanh( (x-0.5d0*t-0.25d0) / (4.0d0*eps) ),
c       where eps is a problem dependent parameter.

c       This code uses loops to produce multiple independent copies of
c       the problem in order to artificially increase the computation
c       required, when npde > 1.

c-----------------------------------------------------------------------
      subroutine f(t, x, u, ux, uxx, fval, npde)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine defines the right hand side vector of the
c       npde dimensional parabolic partial differential equation
c                        ut = f(t, x, u, ux, uxx).
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        x
c                               the current spatial coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
        double precision        uxx(npde)
c                               uxx(1:npde) is the approximation of the
c                               second spatial derivative of the
c                               solution at the point (t,x).
c
c output:
        double precision        fval(npde)
c                               fval(1:npde) is the right hand side
c                               vector f(t, x, u, ux, uxx) of the pde.
c-----------------------------------------------------------------------
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i
c-----------------------------------------------------------------------
c
c     assign fval(1:npde) according to the right hand side of the pde
c     in terms of u(1:npde), ux(1:npde), uxx(1:npde).
c
      do i = 1, npde
         fval(i) = eps*uxx(i) - u(i)*ux(i)
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx, npde)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine is used to define the information about the
c       pde required to form the analytic jacobian matrix for the dae
c       or ode system. assuming the pde is of the form
c                        ut = f(t, x, u, ux, uxx)
c       this routine returns the jacobians d(f)/d(u), d(f)/d(ux), and
c       d(f)/d(uxx).
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        x
c                               the current spatial coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
        double precision        uxx(npde)
c                               uxx(1:npde) is the approximation of the
c                               second spatial derivative of the
c                               solution at the point (t,x).
c
c output:
        double precision        dfdu(npde,npde)
c                               dfdu(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dfdux(npde,npde)
c                               dfdux(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        double precision        dfduxx(npde,npde)
c                               dfduxx(i,j) is the partial derivative
c                               of the i-th component of the vector f
c                               with respect to the j-th component
c                               of the second spatial derivative of the
c                               unknown function u.
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i, j
c-----------------------------------------------------------------------
c
      
c     assign dfdu(1:npde,1:npde), dfdux(1:npde,1:npde), and
c     dfduxx(1:npde,1:npde) according to the right hand side of the pde
c     in terms of u(1:npde), ux(1:npde), uxx(1:npde).
c
      do i = 1, npde
         do j = 1, npde
            dfdu(i,j) = 0.d0
            dfdux(i,j) = 0.d0
            dfduxx(i,j) = 0.d0
         end do
         dfdu(i,i) = -ux(i)
         dfdux(i,i) = -u(i)
         dfduxx(i,i) = eps
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine bndxa(t, u, ux, bval, npde)
c-----------------------------------------------------------------------
c purpose:
c       the subroutine is used to define the boundary conditions at the
c       left spatial end point x = xa.
c                           b(t, u, ux) = 0
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,xa).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,xa).
c
c output:
        double precision        bval(npde)
c                               bval(1:npde) is the boundary contidition
c                               at the left boundary point.
c-----------------------------------------------------------------------
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i
c-----------------------------------------------------------------------
      
      do i = 1, npde
         bval(i) = u(i) - 0.5d0
     &           + 0.5d0 * tanh( (-0.5d0*t-0.25d0) / (4.0d0*eps) )
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine bndxb(t, u, ux, bval, npde)
c-----------------------------------------------------------------------
c purpose:
c       the subroutine is used to define the boundary conditions at the
c       right spatial end point x = xb.
c                           b(t, u, ux) = 0
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,xb).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,xb).
c
c output:
        double precision        bval(npde)
c                               bval(1:npde) is the boundary contidition
c                               at the right boundary point.
c-----------------------------------------------------------------------
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i
c-----------------------------------------------------------------------
      
      do i = 1, npde
c      this leads to cancellation issues (between u(i) and 0.5d))
c      which breaks the forward difference approximation of the
c      jacobian in my tests at sharp (<1d-8) tolerances.
c      bval(i)=u(i)-0.5d0+0.5d0*
c     *        tanh((1.d0-0.5d0*t-0.25d0)/(4.0d0*eps))
         bval(i) = 0.5d0 * tanh((0.75d0-0.5d0*t)/(4.0d0*eps))
     &           - 0.5d0 + u(i)
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine difbxa(t, u, ux, dbdu, dbdux, dbdt, npde)
c-----------------------------------------------------------------------
c purpose:
c       the subroutine is used to define the differentiated boundary
c       conditions at the left spatial end point x = xa. for the
c       boundary condition equation
c                              b(t, u, ux) = 0
c       the partial derivatives db/du, db/dux, and db/dt are supplied
c       by this routine.
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
c output:
        double precision        dbdu(npde,npde)
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dbdux(npde,npde)
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        double precision        dbdt(npde)
c                               dbdt(i) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to time t.
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i, j
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      
c
c     assign dbdu(1:npde,1:npde), dbdu(1:npde,1:npde), and dbdt(1:npde)
c     according to the right boundary condition equation in terms of
c     u(1:npde), ux(1:npde), uxx(1:npde).
c
      do i = 1, npde
         do j = 1, npde
            dbdu(i,j) = 0.0d0
            dbdux(i,j) = 0.0d0
         end do
         dbdu(i,i) = 1.0d0
         dbdt(i) = 0.0625d0/eps
     &           * (1.d0-(tanh(0.25d0/eps*(-0.5d0*t-0.25d0)))**2)
      end do
c      dbdt(1) = -1.d0/eps/sinh(-(0.5d0*t+0.25d0)/(4.d0*eps))**2
c
      return
      end
c-----------------------------------------------------------------------
      subroutine difbxb(t, u, ux, dbdu, dbdux, dbdt, npde)
c-----------------------------------------------------------------------
c purpose:
c       the subroutine is used to define the differentiated boundary
c       conditions at the right spatial end point 1 = xb. for the
c       boundary condition equation
c                              b(t, u, ux) = 0
c       the partial derivatives db/du, db/dux, and db/dt are supplied
c       by this routine.
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        u(npde)
c                               u(1:npde) is the approximation of the
c                               solution at the point (t,x).
c
        double precision        ux(npde)
c                               ux(1:npde) is the approximation of the
c                               spatial derivative of the solution at
c                               the point (t,x).
c
c output:
        double precision        dbdu(npde,npde)
c                               dbdu(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the unknown function u.
c
        double precision        dbdux(npde,npde)
c                               dbdux(i,j) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to the j-th component
c                               of the spatial derivative of the
c                               unknown function u.
c
        double precision        dbdt(npde)
c                               dbdt(i) is the partial derivative
c                               of the i-th component of the vector b
c                               with respect to time t.
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i, j
c-----------------------------------------------------------------------
      
c
c     assign dbdu(1:npde,1:npde), dbdu(1:npde,1:npde), and dbdt(1:npde)
c     according to the right boundary condition equation in terms of
c     u(1:npde), ux(1:npde), uxx(1:npde).
c
      do i = 1, npde
         do j = 1, npde
            dbdu(i,j) = 0.0d0
            dbdux(i,j) = 0.0d0
         end do
         dbdu(i,i) = 1.0d0
         dbdt(i) = 0.0625d0/eps
     &           *(1.d0-(tanh(0.25d0/eps*(-0.5d0*t+0.75d0)))**2)
      end do
c      dbdt(1) = -1.d0/eps/sinh(-(0.5d0*t-0.75d0)/(4.d0*eps))**2
c
      return
      end
c-----------------------------------------------------------------------
      subroutine uinit(x, u, npde)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine is used to return the npde-vector of initial
c       conditions of the unknown function at the initial time t = t0
c       at the spatial coordinate x.
c
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        double precision        x
c                               the spatial coordinate.
c
        integer                 npde
c                               the number of pdes in the system.
c
c output:
        double precision        u(npde)
c                               u(1:npde) is vector of initial values of
c                               the unknown function at t = t0 and the
c                               given value of x.
c-----------------------------------------------------------------------
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i
c-----------------------------------------------------------------------

c
c     assign u(1:npde) the initial values of u(t0,x).
c
      do i = 1, npde
         u(i) = 0.5d0 - 0.5d0 * tanh( (x-0.25d0) / (4.0d0*eps) )
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine truu(t, x, u, npde)
c-----------------------------------------------------------------------
c purpose:
c     this function provides the exact solution of the pde.
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 npde
c                               the number of pdes in the system.
c
        double precision        t
c                               the current time coordinate.
c
        double precision        x
c                               the current spatial coordinate.
c
c output:
        double precision        u(npde)
c                               u(1:npde) is the exact solution at the
c                               point (t,x).
c-----------------------------------------------------------------------
      double precision eps
      common /burger/ eps
c-----------------------------------------------------------------------
c loop indices:
        integer                 i
c-----------------------------------------------------------------------
      
      do i = 1, npde
         u(i) = 0.5d0 - 0.5d0 * tanh( (x-0.5d0*t-0.25d0) / (4.0d0*eps) )
      end do
c
      return
      end
c-----------------------------------------------------------------------
      subroutine header(nout)
c-----------------------------------------------------------------------
c purpose:
c       this subroutine writes a header describing the npde dimensional
c       parabolic partial differential equation
c                        ut = f(t, x, u, ux, uxx).
c-----------------------------------------------------------------------
c subroutine parameters:
c input:
        integer                 nout
c                               nout is the output unit number.
c-----------------------------------------------------------------------
c constants:
        double precision        t0
        parameter              (t0 = 0.0d0)
c
        double precision        xa
        parameter              (xa = 0.0d0)
c
        double precision        xb
        parameter              (xb = 1.0d0)
c-----------------------------------------------------------------------
c
      write(nout,95) 'burgers'' equation:'
      write(nout,95) 'pde:'
      write(nout,95) '   u_t = eps * u_xx - u * u_x , '
      write(nout,95) 'domain:'
      write(nout,96) '   t0 =', t0, ' < t,'
      write(nout,96) '   xa =', xa, ' <= x <= xb =', xb, ','
c
      return
   95 format(a)
   96 format(a,e13.5,a,e13.5,a,e13.5,a,e13.5,a)
      end
