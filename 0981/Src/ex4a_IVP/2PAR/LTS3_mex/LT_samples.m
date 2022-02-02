function [x,U] = LT_samples (Xspan, S, tol, thrds)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
%    ***        EXAMPLE 4a         ***
%	Analytical solution of PDE: u(x,t) = (x+t)*sin(3*(x+t))/6
%
% Solve the IVP
%   U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2,    0 < x < L=Xspan(end)
%   U(0,s)  = s/(s^2+9)^2
%   U'(0,s) = s^2/(s^2+9)^2
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s^2*U1 - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
%

% x: col-wise
% s: row-wise

	options=odeset('AbsTol',tol,'RelTol',tol,'Vectorized','off');

    x = Xspan';
	U = zeros(numel(x),numel(S)); % preallocate output matrix

	parfor (k=1:numel(S), thrds) % PARALLEL ON thrds THREADS
	%for k=1:numel(s)    % SEQUENTIAL FOR
         %     U(0,s(k))   U'(0,s(k))   initial conditions
         U0 = [S(k)/(S(k)^2+9)^2   S(k)^2/(S(k)^2+9)^2];
         [~,sol] = ode45 (@(t,y) odefun(t,y,S(k)), x, U0, options);
        %[~,sol] = ode113(@(t,y) odefun(t,y,s(k)), x, U0, options);
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
%   U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2,    0 < x < L=Xspan(end)
%   U(0,s)  = s/(s^2+9)^2
%   U'(0,s) = s^2/(s^2+9)^2
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s^2*U1 - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2
%

    % s: scalar value
    Uprime = [U(2);
              s^2*U(1) - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2];
end
