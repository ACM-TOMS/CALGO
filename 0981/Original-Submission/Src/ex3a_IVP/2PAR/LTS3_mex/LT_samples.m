function [x,U] = LT_samples (Xspan, S, tol, thrds)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
%    ***        EXAMPLE 3a         ***
%	Analytical solution of PDE: u(x,t) = x*(x-1) +2*t
%
% Solve the IVP
%   U" = s*U - x*(x-1),         0 < x < L=Xspan(end)
%              u(x,0+)
%   U(0,s)   = 2/s^2 = Laplace[u(0,t)],   u(0,t)=2*t
%   U'(0,s)  = -1/s  = Laplace[u_x(0,t)], u_x(0,t)=-1
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s*U1 - x*(x-1)
%

% x: col-wise
% s: row-wise

   %options=odeset('RelTol',tol,'Vectorized','off');
    options=odeset('RelTol',tol,'Vectorized','on');

    x = Xspan';
	U = zeros(numel(x),numel(S)); % preallocate output matrix

	parfor (k=1:numel(S), thrds) % PARALLEL ON thrds THREADS
	%for k=1:numel(s)    % SEQUENTIAL FOR
         %     U(0,s(k))   U'(0,s(k))   initial conditions
         U0 = [2/S(k)^2   -1/S(k)];
         [~,sol] = ode45 (@(t,y) odefun(t,y,S(k)), x, U0, options);
        %[x,sol] = ode113(@(t,y) odefun(t,y,s(k)), x, U0, options);
         U(:,k) = sol(:,1);
	end

    %% U is a MATLAB matrix (col-wise allocated in memory)
    %  The calling program is C and needs a row-wise matrix for the SUM step
    %  U is transposed here and not in the C code.
    U = U.'; % U has complex values: we use transposed, not hermitian, matrix
end


function Uprime = odefun(x,U,s)
%%   ***   USER DEFINED INTERNAL FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% Solve the IVP
%   U" = s*U - x*(x-1),         0 < x < 1=L
%              u(x,0+)
%   U(0,s)   = 2/s^2 = Laplace[u(0,t)],  u(0,t)=2*t
%   U'(0,s)  = -1/s  = Laplace[u_x(0,t)]
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s*U1 - x*(x-1)
%

    % s: scalar value
    Uprime = [U(2);
              s*U(1) - x.*(x-1)];
end

