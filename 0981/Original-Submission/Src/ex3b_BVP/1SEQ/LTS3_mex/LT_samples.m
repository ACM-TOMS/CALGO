function [x,U] = LT_samples (Xspan,s,tol)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% Solve the BVP
%   U" = s*U - x*(x-1),        0 < x < L
%   U(0,s) = 2/s^2
%   U(L,s) = 2/s^2 + L*(L-1)/s
%
% U1 = U(x,s)
% U1' = U2             ==>   U2' = U"
% U2' = s*U1 - x*(x-1)

    options = bvpset('RelTol',tol,'Vectorized','off');

    x = Xspan; % x: row-wise, s: scalar
	U = zeros(numel(x),numel(s)); % preallocate output matrix
    solinit = bvpinit(x,[1 0]); % initial guess (U1=U,U2=U')

    %parfor k=1:numel(s) % PARALLEL FOR
     for k=1:numel(s)    % SEQUENTIAL FOR

       %sol = bvp4c(@(t,U)twoode(t,U,s(k)),@(ya,yb)twobc(ya,yb,s(k),Xspan(end)),solinit,options);
        sol = bvp5c(@(t,U)twoode(t,U,s(k)),@(ya,yb)twobc(ya,yb,s(k),Xspan(end)),solinit,options);

        % In order to avoid errors if the following warning is raised:
        %     Unable to meet the tolerance without using more than 5000 mesh points.
        y = deval(sol,x);   U(:,k) = transpose(y(1,:)); % col-wise
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

    %% NOT VECTORIZED
    % {
    Uprime = [U(2);
              s*U(1) - x.*(x-1)];
    %}

    %% VECTORIZED
    %{
    Uprime = [U(2,:);
              s*U(1,:) - x.*(x-1)];
    %}
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
