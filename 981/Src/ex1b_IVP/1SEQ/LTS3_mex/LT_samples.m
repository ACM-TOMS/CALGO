function [x, U] = LT_samples (Xspan, S, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              PDE EXAMPLE 1b
%
%     Solve, by means of ode45() MATLAB solver, ODE problems coming from the
%     application of the Laplace method to the following PDE:
% 
%                 u_t (x,t) = u_x (x,t),      x>X0>0, t>0
%    
%    in [X0,X1]x[0,T], with the conditions
%
%           u(X0,t) = (X0+t)*sin(3*(X0+t))/6,    t>0   [boundary condition]
%           u(x,0+) = x*sin(3*x)/6,              x>X0  [initial condition]
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
%				U' = s*U - u(x,0+)
%           	U(X0,s) = [sin(3*X0)+3*X0*cos(3*X0)+s*X0*sin(3*X0)]/[6*(s^2+9)]
%               	     - [3*sin(3*X0)-s*cos(3*X0)]/(s^2+9)^2
%
%     The analytical solution of the ODE problem is:
%
%           	U(x,s) = [sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)]/[6*(s^2+9)]
%	                    - [3*sin(3*x)-s*cos(3*x)]/(s^2+9)^2
%
%     U(x,s) is a LT fun with a pole of 2nd order at s=+/-3i so that
%     geometrical and accuracy parameters for Talbot's method have to be
%     computed with the following informations on the singularities of LT:
%                   Nsings = 1;
%                   Ssings = 3i;
%                   MULT = 2;
% 
%     The analytical solution of the above PDE is: u(x,t) = (x+t)*sin(3*(x+t))/6
% 

%     >>>>>>>>>>>>>        VERSION 2.0     Jul 13th, 2013       <<<<<<<<<<
%                             by Mariarosaria Rizzardi
%

    SIN3 = sin(3*Xspan(1));	  COS3 = cos(3*Xspan(1));
	S29  = S.^2+9;
	% initial conditions
	U0 = (SIN3 + 3*Xspan(1).*COS3+S*Xspan(1).*SIN3)./(6*S29) - (3*SIN3-S*COS3)./S29.^2;
    options = odeset('AbsTol',tol,'RelTol',tol,'Vectorized','off');
    
    % SOLVE ALL THE ODEs (one for each S(j))
    %           U(i,j) = U(x(i), S(j)) = U(x,s)
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

	dydt = s.'.*U - x*sin(3*x)/6; % -u(x,0+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
