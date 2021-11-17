function yt=f_01(x,y,mode)
% Function <<f_01>> gives the functionality for the cubic oscillator.
% Usage:
% yt=f_01(x,y,mode)
% INPUT:
%           : Call problem(0,0,1) to display a message
%           : Call problem(h,k,2) to calculate the first k-values with
%           :       step size -h
%           : Call problem(Y,n_omega,3) to calculate a set of
%           :       n_omega frequencies given the current value Y
%           : Call problem(Y,t,4) to calculate the Jacobian given
%           :       the current value Y and the current time -t
%           : Call problem(t,Y) to calculate the derivative at time -t
%           :       given Y=Y(t)
% OUTPUT:
%   yt      : depending on mode
%--------------------------------------------------------------------------
if nargin==3
    switch mode
        case 1  % Display a message
            disp('Cubic Oscillator:y''(t)=0.001*y(t)^3-y(t)');
        case 2  % Calculate initial values
            yt=zeros(y,2);
            yt(1,1)=1;
            yt(1,2)=0;
            for j=2:y
                [t,val,dval,ia,ir]=rkn86(inline('0.001*yy^3-yy','yy'),(j-1)*x,j*x,yt(j-1,1),yt(j-1,2),1e-12);
                yt(j,1)=val(ia+1);
                yt(j,2)=dval(ia+1);
            end
        case 3  % Calculate frequencies
            yt=sqrt(1-0.75*0.001);
            for j=2:y
                yt=[yt;3*yt(1)];
            end
        case 4  % Calculate Jacobian
            yt=[0,1;0.003*x(1)^2-1,0];
    end
    return;
end
% Calculate derivative
yt=[y(2),0.001*y(1)^3-y(1)];
