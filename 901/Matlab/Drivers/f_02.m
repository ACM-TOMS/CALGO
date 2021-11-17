function yt=f_02(x,y,mode)
% Function <<f_02>> gives the functionality for the 2-body problem.
% Usage:
% yt=f_02(x,y,mode)
% INPUT:
%           : Call problem(0,0,1) to display a message
%           : Call problem(h,k,2) to calculate the first k-values with
%           :       step size -h
%           : Call problem(Y,n_omega,3) to calculate a set of
%           :       n_omega frequencies given the current value Y
%           : Call problem(Y,t,4) to calculate the Jacobian given
%           :       the current value Y and the current time -t
%           : Call f_02(ecc,0,5) to set the eccentricity to the 
%                   value ecc
%           : Call problem(t,Y) to calculate the derivative at time -t
%           :       given Y=Y(t)
% OUTPUT:
%   yt      : depending on mode
%--------------------------------------------------------------------------
persistent ecc;
if isempty(ecc)
    ecc=0.0;
end
if nargin==3
    switch mode
        case 1  % Display a message
            v=sprintf('2-body problem:r''=-r/|r|^3. Eccentricity=%f',ecc);
            disp(v);
        case 2  % Calculate initial values
            yt=zeros(y,4);
            yt(1,1)=1-ecc;
            yt(1,2)=0;
            yt(1,3)=0;
            yt(1,4)=sqrt((1+ecc)/(1-ecc));
            for j=2:y
                u=(j-2)*x;
                t=(j-1)*x;
                err=1;
                reps=0;
                while (abs(err)>1e-15)&(reps<20)
                    err=-(u-ecc*sin(u)-t)/(1-ecc*cos(u));
                    u=u+err;
                    reps=reps+1;
                end
                yt(j,1:2)=[cos(u)-ecc,sqrt(1-ecc^2)*sin(u)];
                yt(j,3:4)=[-sin(u)/(1-ecc*cos(u)),sqrt(1-ecc^2)*cos(u)/(1-ecc*cos(u))];
            end
        case 3  % Calculate frequencies
            %yt=1/(1-ecc^2-ecc*x(1));
            %yt=(x(1)^2+x(2)^2)^(-3/4);
            yt=(x(1)^2+x(2)^2)^(-3/2);
            for j=2:y
                yt=[yt;j*yt(1)];
            end
        case 4  % Calculate Jacobian
            Aux1=(x(1)^2+x(2)^2)^(-5/2);
            Aux2=(x(1)^2+x(2)^2)^(-3/2);
            Aux3=3*x(1)*x(2)*Aux1;
            yt=[0,0,1,0;0,0,0,1;3*x(1)^2*Aux1-Aux2,Aux3,0,0;Aux3,3*x(2)^2*Aux1-Aux2,0,0];
        case 5 % Set the eccentricity
            ecc=x;
    end
    return;
end
% Calculate derivative
U=(y(1)^2+y(2)^2)^(3/2);
yt=[y(3),y(4),-y(1)/U,-y(2)/U];
