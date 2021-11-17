function [Y,T]=wrapping(c,problem,h,N)
% Function <<wrapping>> tests a linear multistep method.
% Usage:
% wrapping(c,problem,h,N)
% INPUT:
%   c       : the coeffs of the method in symbolic form [a,b]
%   problem : a handle to the function f, such that y'(t)=f(t,y(t))
%           : The detailed communication protocol with this function
%          	: is the following:
%           : Call problem(0,0,1) to display a message
%           : Call problem(h,k,2) to calculate the first k-values with
%           :       step size -h
%           : Call problem(Y,n_omega,3) to calculate a set of
%           :       n_omega frequencies given the current value Y
%           : Call problem(Y,t,4) to calculate the Jacobian given
%           :       the current value Y and the current time -t
%           : Call problem(t,Y) to calculate the derivative at time -t
%           :       given Y=Y(t)
%   h       : the integration step
%   N       : the iterations number
% OUTPUT:
%   Y       : the values of the calculated function. The value
%             at position -i-, corresponds to t=(i-1)*h
%   T       : a vector with the times at which Y is calculated
%--------------------------------------------------------------------------
% Display comment about problem
problem(0,0,1);
% Initialize
k=max(size(c))/2;   % The number of steps
y0=problem(h,k-1,2);% The k-1 initial values
[tt,s]=size(y0);    % s is the dimension of the ode
% Create arrays
Y=zeros(N,s);       % Here the calculations are stored
Yt=zeros(k,s);      % Store the last (k-1) derivatives
Y(1:k-1,:)=y0(1:k-1,:);
for j=1:k-1
    Yt(j,:)=problem((j-1)*h,Y(j,:));
end
% Find symbolic variables in coeffs
sv=findsym(c);      % Get the symbolic variables
q='h';              % One of them is h
j=0;
if ~isempty(sv)
    [v,remain]=strtok(sv,', ');
    while ~isempty(v)
        if ~strcmp(v,'h')
            j=j+1;
            q=strvcat(q,v);
        end
        [v,remain]=strtok(remain,', ');
    end
end
omegas=j;           % This is the number of frequencies
tic;
for j=k:N
    % Evaluate coeffs. This is done inside the loop in order
    % to take into account the new frequency
    om=problem(Y(j-1,:),omegas,3);  % Get the frequencies
    r=subs(c,q(1),h);               % and substitute
    for l=1:omegas
        r=subs(r,q(l+1,:),om(l));
    end
    % Finally, get the a's and b's
    a=double(r(1:k));
    b=double(r(k+1:2*k));
    % Calculate the contribution of the past values
    S=zeros(1,s);
    for l=1:k-1
        S=S+h*b(l+1)*Yt(k-l,:)-a(l+1)*Y(j-l,:);
    end
    if (b(1)==0)
        Y(j,:)=S(:)/a(1);
    else
        % Use Newton-Raphson to solve 
        Y(j,:)=Y(j-1,:);    % Start with the previous value
        err=ones(s,1);      % Initialize error (just for entering the loop)
        watch_=0;           % This is a watchdog for breaking the loop
        while (abs(err)>1e-15)&(watch_<20)
            J_1=problem(Y(j,:),j*h,4);
            J=a(1)*eye(s)-h*b(1)*J_1;
            J=inv(J);
            F=a(1)*Y(j,:)-h*b(1)*problem(j*h,Y(j,:))-S;
            err=J*F';
            Y(j,:)=Y(j,:)-err';
            watch_=watch_+1;
        end
    end
    % Shift previous values
    for l=1:k-2
        Yt(l,:)=Yt(l+1,:);
    end
    Yt(k-1,:)=problem(j*h,Y(j,:));
    if mod(j,1000)==0
        rest=(N-j)/1000;
        f=toc;
        v=sprintf('Remaining time~%f sec',rest*f);
        disp(v);
        tic;
    end
end
% Fill the array T with times at which the function has been evaluated
T=zeros(N,1);
for j=1:N
    T(j)=(j-1)*h;
end