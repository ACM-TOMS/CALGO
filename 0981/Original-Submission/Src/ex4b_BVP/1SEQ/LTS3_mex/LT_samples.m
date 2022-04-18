function [x,U] = LT_samples (Xspan,s,tol)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% Solve the BVP
%   U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2,      0<x<L=Xspan(end)
%   U(0,s) = s/(s^2+9)^2
%   U(L,s) = (3*L*cos(3*L) + (1+s*L)*sin(3*L))/(6*(s^2+9)) - (3*sin(3*L) - s*cos(3*L))/(s^2+9)^2
%
%   L=2*pi  ==>  U(2*pi,s) = s/(s^2+9)^2 + pi/(s^2+9)
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s^2*U1 - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2

% x: col-wise
% s: row-wise

	options = bvpset('AbsTol',tol,'RelTol',tol);

    x = Xspan;
	U = zeros(numel(x),numel(s)); % preallocate output matrix
    solinit = bvpinit(x,[1 0]); % initial guess (U1=U,U2=U')
   %solinit = bvpinit(x,[2 0]);

    %parfor k=1:numel(s) % PARALLEL FOR
     for k=1:numel(s)    % SEQUENTIAL FOR
       %solinit = bvpinit(x,[2/s(k)^2 0]);

       %sol = bvp4c(@(t,U)twoode(t,U,s(k)),@(ya,yb)twobc(ya,yb,s(k),Xspan(end)),solinit,options);
        sol = bvp5c(@(t,U)twoode(t,U,s(k)), @(ya,yb)twobc(ya,yb,s(k),Xspan(end)), solinit, options);

        % In order to avoid errors if the following warning is raised:
        %     Unable to meet the tolerance without using more than 5000 mesh points.
        y = deval(sol,x);
        U(:,k) = transpose(y(1,:)); % col-wise
    end

    %% U is a MATLAB matrix (col-wise allocated in memory)
    %  The calling program is C and needs a row-wise matrix for the SUM step
    %  U is transposed here and not in the C code.
    U = U.'; % U has complex values: we use transposed, not hermitian, matrix
end


function Uprime = twoode(x,U,s)
%%   ***   USER DEFINED INTERNAL FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
%
%% Solve the BVP
%   U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2,      x in ]0,L[
%   U(0,s)    = s/(s^2+9)^2
%   U(L,s) = (3*L*cos(3*L) + (1+s*L)*sin(3*L))/(6*(s^2+9)) - (3*sin(3*L) - s*cos(3*L))/(s^2+9)^2
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s^2*U1 - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2

    % s: scalar value
    Uprime = [U(2);
              s^2*U(1) - (s*x+1).*sin(3*x)/6 - x.*cos(3*x)/2];
end


function res = twobc(ya,yb,s,L)
%%   ***   USER DEFINED INTERNAL FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% BVP: U" = s^2*U - (s*x+1)*sin(3*x)/6 - x*cos(3*x)/2,      x in ]0,L[
%      U(0,s) = s/(s^2+9)^2
%      U(L,s) = (3*L*cos(3*L) + (1+s*L)*sin(3*L))/(6*(s^2+9)) - (3*sin(3*L) - s*cos(3*L))/(s^2+9)^2

    % s: a single value (parameter)
          % U(0,s) - ... = 0
    res = [ ya(1) - s/(s^2+9)^2;
          % U(L,s) - ... = 0
            yb(1) - (3*L*cos(3*L) + (s*L+1)*sin(3*L))/(6*(s^2+9)) + (3*sin(3*L) - s*cos(3*L))/(s^2+9)^2];
end
