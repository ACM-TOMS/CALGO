function [x, U] = LT_samples (Xspan, S, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              PDE EXAMPLE 1a
%
%     Solve, by means of ode45() MATLAB solver, ODE problems coming from the
%     application of the Laplace method to the following PDE (uniform
%     transport equation):
% 
%                 u_t (x,t) = u_x (x,t),      x>X0, t>0
%    
%     with the conditions
% 
%                 u(X0,t) = X0+t,    t>0    [boundary condition]
%                 u(x,0+) = x,       x>X0   [initial condition]
% 
%     The ODE problems are solved in the interval [Xspan(1), Xspan(end)] with
%     an error tolerance tol. Each of them has a different initial value which
%     depends on S(k).
% 
%     The parameter Nx manages the output values x where the solution will
%     be computed: if Nx > 2 then the solution has to be computed at Nx
%     equispaced values in the above interval. Otherwise the values of x are
%     chosen by the ODE solver in the same interval.
% 
%     To solve an ODE problem as U' = f(x,U) the function f(x,U) is required:
%     it is implemented by the private function odefun() below.
% 
%     Applying the LT:    U(x,s) = LT[u(x,t)]    to both the sides of the PDE,
%     we get the ODE [derivatives with respect to x]
%       U' = s*U - u(x,0+)
%       U(x0,s) = x0/s + 1/s^2    [initial condition]
% 
%     The analytical solution of the ODE is: U(x,s) = x/s + 1/s^2
%     U(x,s) is a LT fun with a pole of 2nd order at s=0 so that
%     geometrical and accuracy parameters for Talbot's method have to be
%     computed with the following informations on the singularities of LT:
%                   Nsings = 1;
%                   Ssings = 0;
%                   MULT = 2;
% 
%     The analytical solution of the above PDE is: u(x,t) = x + t
% 
%     >>>>>>>>>>>>>      PAR VERSION 2.0      Jul 5th, 2015      <<<<<<<<<<
%                             by Mariarosaria Rizzardi
%     The MATLAB Parallel Computing Toolbox is required.
%

    U0 = Xspan(1)./S + 1./S.^2; % initial conditions
    options = odeset('RelTol',tol,'Vectorized','on');

    %% SOLVE ALL THE ODEs (one for each S(k))
    %           U(h,k) = U(x(h), S(k)) = U(x,s)
    [x,U] = ode45(@(x,U) odefun(x,U,S), Xspan, U0, options);


    %% U is a MATLAB matrix (col-wise allocated in memory)
    %  The calling program is C and needs a row-wise matrix for the SUM step
    %  U is transposed here and not in the C code.
    U = U.'; % U has complex values: we use transposed, not hermitian, matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function dydt = odefun(x,U,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     ODE function f(x,U) for
%         y'=f(x,U)

    dydt = s.'.*U - x; % -u(x,0+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
