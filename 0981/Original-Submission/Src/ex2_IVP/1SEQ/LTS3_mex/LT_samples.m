function [x, U] = LT_samples (Xspan, S, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              PDE EXAMPLE 2
%
%     Solve, by means of ode45() MATLAB solver, ODE problems coming from the
%     application of the Laplace method to the following PDE:
% 
%                 u_{tx} (x,t) = exp(-x)*cos(t),      x>0, t>0
%    
%     with the conditions
% 
%                 u_x(x,0+) = 0,     x>0
%                 u(0,t)    = 0,     t>0
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
%               U' = exp(-x)/(s^2+1),        x>0
%               U(0,s) = 0
% 
%     The analytical solution of the ODE is:
%
%                         1 - exp(-x)
%               U(x,s) = -------------
%                           s^2 + 1
%
%     U(x,s) is a LT fun with simple poles at s=+/-i so that
%     geometrical and accuracy parameters for Talbot's method have to be
%     computed with the following informations on the singularities of LT:
%                   Nsings = 1;
%                   Ssings = 1i;
%                   MULT = 1;
% 
%     The analytical solution of the above PDE is: u(x,t) = sin(t)*(1-exp(-x)).
% 
%     >>>>>>>>>>>>>        VERSION 3.0     May 15th, 2015       <<<<<<<<<<
%                             by Mariarosaria Rizzardi
%

    U0 = zeros(size(S)); % initial conditions
     options = odeset('AbsTol',tol,'RelTol',tol,'Vectorized','off');
    
    % SOLVE ALL THE ODEs (one for each S(j))
    %           U(i,j) = U(x(i), S(j)) = U(x,s)
     [x,U] = ode45(@(x,U) odefun(x,U,S), Xspan, U0, options);
    %[x,U] = ode113(@(x,U) odefun(x,U,S), Xspan, U0, options);

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

	dydt = exp(-x)./((s.').^2+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
