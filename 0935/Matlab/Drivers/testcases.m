% Test cases 1-23 for IIPBF outlined in Ratnanather et. al
% Author: Aug 30, 2012 Jung Hun Kim
% Actual tests used in paper are cases 1-18; cases 19-23 optional

function [results, absolute, relative, nevals, required] = testcases(n,varargin)
tic;format long e

switch n
    case 1
        % f(x) = 1
        % B = J_0(x) J_1(3*x/2)
        rho = 1; tau = 1.5;
        a = 0; b = 1;
        type = 'JJ';
        answer = 2/3;
        fx = @(x)1;
    case 2
        % f(x) = power(x,-4);
        % B = J_0(x) J_5(2*x)
        rho = 1; tau = 2;
        a = 0; b = 5;
        type = 'JJ';
        answer = 27/4096;
        fx = @(x)power(x,-4);
    case 3
        % f(x) = x/(1+x^2);
        % B = J_0(x) J_20(1.1*x)
        rho = 1; tau = 1.1;
        a = 0; b = 20;
        type = 'JJ';
        answer = -6.050747903049e-3;
        fx = @(x)x ./(1+power(x,2));
    case 4
        % f(x) = 1/x
        % B = J_0(x) J_1(x)
        rho = 1;tau = 1;
        a = 0;b = 1;
        type = 'JJ';
        answer = 2/pi;
        fx = @(x)1./x;
    case 5
        % f(x) = power(x,-2)
        % B = J_1(x) J_1(x)
        rho = 1;tau = 1;
        a = 1;b = 1;
        type = 'JJ';
        answer = 4/(3*pi);
        fx = @(x)1./power(x,2);
    case 6
        % f(x) = x * K_0(u*x);
        % B = J_0(rho*x) J_0(tau*x)
        u = varargin{1};
        rho = varargin{2};tau = varargin{3};
        a = 0; b = 0;
        type = 'JJ';
        temp1 = power(power(u,2) + power(rho,2) + power(tau,2),2);
        temp2 = 4*power(rho,2)*power(tau,2);
        answer = power(temp1-temp2,-.5);
        fx = @(x) x.*besselk(0,u.*x);
    case 7
        % f(x) = x^2 * K_1(u*x)
        % B = J_1(rho*x) J_1(tau*x)
        u = varargin{1};
        rho = varargin{2};tau = varargin{3};
        a = 1;b = 1;
        type = 'JJ';
        temp1 = power(power(u,2) + power(rho,2) + power(tau,2),2);
        temp2 = 4*power(rho,2)*power(tau,2);
        answer= 4*u*rho*tau/power(temp1-temp2,3/2);
        fx = @(x) power(x,2).*besselk(1,u.*x);
    case 8
        % f(x) = exp(-2*u*x);
        % B =  J_0(x) Y_0(x)
        u = varargin{1};
        rho = 1;tau = 1;
        a = 0;b = 0;
        type = 'JY';
        temp1 = @(y)(-1./sqrt(1+(u.^2).*power(cos(y),2)))./(pi);
        answer = quadgk(temp1,0,pi/2);
        fx =@(x) exp(-2.*u.*x);
    case 9
        % f(x) = x * exp(-(x^2) / u);
        % B = J_2(x) Y_2(x)
        u = varargin{1};
        rho = 1;tau = 1;
        a = 2;b = 2;
        type = 'JY';
        answer=4/(u*pi) - 2/pi - (u*besselk(2,.5*u))/(2*pi*exp(.5*u));
        fx =@(x) x.*exp(-power(x,2)./u);
    case 10
        % f(x) = x^3 * exp(-(x^2) / u);
        % B = J_2(x) Y_2(x)
        u = varargin{1};
        rho = 1;tau = 1;
        a = 2;b = 2;
        type = 'JY';
        answer = -4./(pi) + ...
            (power(u,2).*(2+u).*besselk(0,u./2))./(4.*pi.*exp(u./2)) + ...
            (u.*(8+4.*u+power(u,2)).*besselk(1,u/2))./(4.*pi.*exp(u./2));
        fx =@(x) power(x,3).*exp(-power(x,2)./u);
        
    case 11
        % f(x) = exp(-u*x);
        % B = Y_0(rho*x) Y_0(tau*x)
        % rho > tau
        u=varargin{1};
        rho=varargin{2};tau=varargin{3};
        a = 0; b = 0;
        type='YY';
        answer= quadgk(case11(u,rho,tau),rho,Inf);
        fx = @(x)exp(-u*x);
    case 12
        % f(x) = exp(-u*x);
        % B = Y_0(rho*x) Y_0(tau*x)
        % rho > tau
        u=varargin{1};
        rho=varargin{2};tau=varargin{3};
        a = 0; b = 0;
        type='JJ';
        temp1=1/(sqrt(u^2 + (rho+tau)^2));
        temp=2/pi*temp1;
        answer= temp*ellipk(2*sqrt(rho*tau)*temp1);
        fx = @(x)exp(-u*x);
    case 13
        % f(x) = exp(-u*x);
        % B = J_0(rho*x) Y_0(tau*x)
        u=varargin{1};
        rho=varargin{2};tau=varargin{3};
        a = 0; b = 0;
        type='JY';
%         answer= quadgk(case13(u,rho,tau),0,rho);
        answer= quadgk(case13(u,rho,tau),tau,Inf);
        fx = @(x)exp(-u*x);
    case 14
        % f(x) = 1;
        % B = J_0(rho*x) J_0(tau*x)
        rho=varargin{1};tau=varargin{2};
        a = 0;b = 0;
        type='JJ';
        answer = (2/(pi.*rho))*ellipke(power(tau/rho,2));
        fx = @(x)1;
    case 15
        % f(x) = 1;
        % B = Y_0(rho*x) Y_0(tau*x)
        rho=varargin{1};tau=varargin{2};
        a = 0; b = 0;
        type = 'YY';
        answer = (2/(pi.*rho))*ellipke(power(tau/rho,2));
        fx = @(x)1;
    case 16
        % f(x) = 1;
        % B = J_0(rho*x) Y_0(tau*x)
        rho=varargin{1};tau=varargin{2};
        a = 0; b = 0;
        type = 'JY';
        if rho>=tau
            answer=- (2/(pi.*rho))*ellipke(1-(tau^2/rho^2));
        else
            disp('rho<tau');
            answer =(2/(pi.*rho))*ellipke(1-(rho^2/tau^2));
        end
        fx = @(x)1;
    case 17
        % f(x) = x^-.5 * exp(-u*x);
        % B = J_0(rho*x) Y_0(tau*x)
        u = varargin{1};
        rho=varargin{2};tau=varargin{3};
        a = 1; b = .5;
        type = 'JJ';
        rho_tau2=power(rho+tau,2);rho_tau=power(rho-tau,2); u2=power(u,2);
        l_1 = .5*(sqrt(rho_tau2+u2)-sqrt(rho_tau+u2));
        answer = (sqrt(2)/(sqrt(tau*pi)*rho)) * (tau - sqrt(power(tau,2) - power(l_1,2)));
        fx = @(x)power(x,-.5).*exp(-u.*x);
    case 18
        % f(x) = 1;
        % B = J_0(rho*x) Y_0(tau*x)
        u = varargin{1};
        rho=varargin{2};tau=varargin{3};
        a = 2; b = 1.5;
        type = 'JJ';
        rho_tau2=power(rho+tau,2);rho_tau=power(rho-tau,2); u2=power(u,2);
        l_1 = .5*(sqrt(rho_tau2+u2)-sqrt(rho_tau+u2));
        temp1= (2*power(tau,1.5))/(sqrt(2*pi)*power(rho,2));
        tau2=power(tau,2);l_12=power(l_1,2);
        answer=temp1 * ((2/3) - (power(tau2-l_12,.5)/tau) + (power(tau2-l_12,1.5)/(3*tau2*tau)) ); 
        fx = @(x)power(x,-.5).*exp(-u.*x);
    case 19
        % f(x) = x*e^{-2ux);
        % B = J_0(rho*x) J_1(rho*x)
        u=varargin{1};
        rho=varargin{2};tau=rho;
        if (rho ~= tau)
            disp('rho==tau');exit;
        end
        k=rho/sqrt(rho^2+u^2);
        [K,E]=ellipke(k^2);
        a=0;b=1;
        type= 'JJ';
        answer=(K-E)/(2*pi*rho*sqrt(rho^2+u^2));
        fx = @(x)x.*exp(-2*u*x);    
    case 20
        % f(x) =e^{-2ux);
        % B = J_0(rho*x) J_0(rho*x)
        u=varargin{1};
        rho=varargin{2};tau=rho;
        if (rho ~= tau)
            error('rho has to equal tau');
        end
        k=rho/sqrt(rho^2+u^2);
        [K,E]=ellipke(k^2);
        a=0;b=0;
        type= 'JJ';
        answer=K/(pi*sqrt(rho^2+u^2));
        fx = @(x)exp(-2*u*x);
        
    case 21
        % f(x) =e^{-2ux);
        % B = J_1(rho*x) J_1(rho*x)
        u=varargin{1};
        rho=varargin{2};tau=rho;
        if (rho ~= tau)
            disp('rho has to equal tau');
        end
        k=rho/sqrt(rho^2+u^2);
        [K,E]=ellipke(k^2);
        a=1;b=1;
        type= 'JJ';
        answer=((2*u^2+rho^2)*K-2*(u^2+rho^2)*E)/(pi*rho^2*sqrt(rho^2+u^2));
        fx = @(x)exp(-2*u*x);
    case 22
        % f(x) = x^2 exp(-2ux);
        % B = J_1(x) J_1(x)
        u=varargin{1};
        rho=1;tau=1;%unused
        
        a=1;b=a;
        type= 'JJ';
        answer=3/(4*pi)*quadgk(@(x)(cos(x)).^2./(u^2+(cos(x)).^2).^(2.5),0,pi/2);
        fx = @(x)x.^2.*exp(-2*u*x);
    case 23
        % f(x) = x^(b-a+1)/(u^2+x^2);
        % B = J_a(rho*x) J_a(tau*x)
        u=varargin{3};
        rho=varargin{4};
        tau=varargin{5};
        a=varargin{1};
        b=varargin{2};
        type= 'JJ';
        answer=u^(b-a)*besseli(a,rho*u)*besselk(a,tau*u);
        fx = @(x)x.^(b-a+1)./(u.^2+x.^2);
    otherwise
        error('Cases 1-23 only');
end;

absolute = zeros(11,1);
relative = zeros(11,1);
results = zeros(11,1);
nevals = zeros(11,1);
required = zeros(11,1);

for i = 1:11
    
    abserr = (1.0e-3)/power(10,i);
    relerr = (1.0e-11);
    
    [result,err, neval]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type);
    results(i) = result;
    required(i) = abserr;
    relative(i) = err;
    absolute(i) = abs(answer-result);
    nevals(i) = neval;
    
end;

toc
