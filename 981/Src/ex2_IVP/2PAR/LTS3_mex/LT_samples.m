function [x, U] = LT_samples (Xspan, S, tol, thrds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              PDE EXAMPLE 2
%
%     Solve, by means of ode45() MATLAB solver, ODE problems coming from the
%     application of the Laplace method to the following PDE problem
% 
%                 u_tx (x,t) = exp(-x)*cos(t),    x>0,    t>0
%    
%     with the conditions
% 
%                 u_x(x,0+)  = 0
%                 u(0,t)     = 0
% 
%     The ODE problems are solved in the interval [Xspan(1), Xspan(end)] with
%     an error tolerance tol. Each of them has a different initial value which
%     depends on S(k).
% 
%     To solve an ODE problem as U' = f(x,U) the function f(x,U) is required:
%     it is implemented by the private function odefun() below.
% 
%     Applying the LT:    U(x,s) = LT[u(x,t)]    to both the sides of the PDE,
%     we get the ODE [derivatives with respect to x]
%       U' = exp(-x)/(s^2+1),        x>0
%       U(0) = 0
% 
%     The analytical solution of the ODE is:
%
%                     1-exp(-x)
%           U(x,s) = -----------
%                       s^2+1
%
%     U(x,s) is a LT fun with simple poles at s=+/-1i so that
%     geometrical and accuracy parameters for Talbot's method have to be
%     computed with the following informations on the singularities of LT:
%                   Nsings = 1;
%                   Ssings = 1i;
%                   MULT = 1;
% 
%     The analytical solution of the above PDE is: u(x,t) = sin(t)*(1-exp(-x))
% 
%     >>>>>>>>>>>>>      PAR VERSION 2.0      Jul 5th, 2015      <<<<<<<<<<
%                             by Mariarosaria Rizzardi
%     The MATLAB Parallel Computing Toolbox is required.
%

    % initial conditions
    U0 = zeros(size(S));

    options = odeset('AbsTol',tol,'RelTol',tol);
    
    % SOLVE THE ODEs (one for each S(k))
    U=zeros(numel(Xspan),numel(S));
    x=Xspan';

    parfor (k=1:numel(S), thrds) % PARALLEL ON thrds THREADS
        [~,U(:,k)] = ode45(@(x,y) odefun(x,y,S(k)), x, U0(k), options);
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
%         U'=f(x,U)

	dydt = exp(-x)/(s^2+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
