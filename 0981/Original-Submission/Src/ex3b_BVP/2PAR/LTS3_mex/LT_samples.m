function [x,U] = LT_samples (Xspan,S,tol, thrds)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% Solve the BVP
%   U" = s*U - x*(x-1),        0 < x < L
%   U(0,s) = 2/s^2
%   U(L,s) = 2/s^2 + L*(L-1)/s
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s*U1 - x*(x-1)

    options = bvpset('RelTol',tol);

    x = Xspan'; % x: col-wise
    Xend = Xspan(end);
	U = zeros(numel(x),numel(S)); % preallocate output matrix
    solinit = bvpinit(x,[1 0]); % initial guess (U1=U,U2=U')

    parfor (k=1:numel(S), thrds) % PARALLEL ON thrds THREADS
        sol = bvp5c(@(t,U)twoode(t,U,S(k)),@(ya,yb)twobc(ya,yb,S(k),Xend),solinit,options);
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
% 	U" = s*U - u(x,0+)   x in ]0,1[
%      u(x,0+) = x*(x-1)
%   U(0,s)=2/s^2  <==> U(0,1)-2/s^2 = 0
%   U(1,s)=2/s^2  <==> U(1,s)-2/s^2 = 0
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s*U1 - x*(x-1)

    % s: scalar value

    %% NOT VECTORIZED ON U
    Uprime = [U(2);
              s*U(1) - x.*(x-1)];
end


function res = twobc(ya,yb,s,L)
%%   ***   USER DEFINED INTERNAL FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% BVP: U" = s*U - u(x,0+)   x in ]0,L[
%      u(x,0+) = x*(x-1)
%      U(0,s)=2/s^2           <==> U(0,1)-2/s^2 = 0
%      U(L,s)=2/s^2+L*(L-1)/s <==> U(L,s)-2/s^2-L*(L-1)/s = 0

    % s: a single value (parameter)
    res = [ ya(1) - 2/s^2;
            yb(1) - 2/s^2 - L*(L-1)/s];
end
