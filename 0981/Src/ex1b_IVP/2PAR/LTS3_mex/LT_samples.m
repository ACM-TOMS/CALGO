function [x, U] = LT_samples (Xspan, S, tol, thrds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              PDE EXAMPLE 1b
%
%     Solve, by means of ode45() MATLAB solver, ODE problems coming from the
%     application of the Laplace method to the following PDE (uniform
%     transport equation):
% 
%                 u_t (x,t) = u_x (x,t),      x>X0, t>0
%    
%     with the conditions
% 
%                 u(x,0+) = x*sin(3*x)/6
%                 u(X0,t) = (X0+t)*sin(3*(X0+t))/6
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
%       U(x0,s) = (sin(3*X0)+3*X0*cos(3*X0)+s*X0*sin(3*X0))/(6*(s^2+9)) - (3*sin(3*X0)-s*cos(3*X0))/(s^2+9)^2
% 
%     The analytical solution of the ODE is:
%
%                    sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)   3*sin(3*x)-s*cos(3*x)
%           U(x,s) = ---------------------------------- - ---------------------
%                                 6*(s^2+9)                     (s^2+9)^2
%
%     U(x,s) is a LT fun with a pole of 2nd order at s=0 so that
%     geometrical and accuracy parameters for Talbot's method have to be
%     computed with the following informations on the singularities of LT:
%                   Nsings = 1;
%                   Ssings = 3i;
%                   MULT = 2;
% 
%     The analytical solution of the above PDE is: u(x,t) = (x+t)*sin(3*(x+t))/6
% 
%     >>>>>>>>>>>>>      PAR VERSION 2.0      Jul 5th, 2015      <<<<<<<<<<
%                             by Mariarosaria Rizzardi
%     The MATLAB Parallel Computing Toolbox is required.
%

    % initial conditions
    U0 = (sin(3*Xspan(1))+3*Xspan(1)*cos(3*Xspan(1))+S*Xspan(1)*sin(3*Xspan(1)))./(6*(S.^2+9)) ...
       - (3*sin(3*Xspan(1))-S*cos(3*Xspan(1)))./(S.^2+9).^2;

    options = odeset('AbsTol',tol,'RelTol',tol);
    
    % SOLVE IN PARALLEL (with PCT) ALL THE ODEs (one for each S(k))
    U=zeros(numel(Xspan),numel(S));
    x=Xspan;
    parfor (k=1:numel(S), thrds) % PARALLEL ON thrds THREADS
        [~,U(:,k)] = ode45(@(x,y) odefun(x,y,S(k)), Xspan, U0(k), options);
    end

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
