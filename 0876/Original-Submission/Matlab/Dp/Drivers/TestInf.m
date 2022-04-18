function TestInf(num_problem,lambda,AbsTol,RelTol,section)
% Test problems 1,...,14 for solving an integral equation of the form
%                        Inf
%  lambda*x(s) - Integral    K(s,t)*x(t) dt = RHS(s)
%                         0
% The kernel K(s,t) is assumed to be moderately smooth on all of 
% [0,Inf)x[0,Inf) and a compact operator on a suitable space.
%
% The default value of lambda is 1.
%
% The default tolerances are AbsTol = 1e-4 and RelTol = 1e-2. 
%
% The parameter section affects the display of results and the assessment
% of error. The continuous extension will be evaluated at 150 equally 
% spaced points in [0,section].  The default value of section is 10.
%
% As an example, use     TestInf(10,2,[],[],15)
%
% TestInf calls the function intmi.

behavior = 0; 

if nargin < 2 || isempty(lambda)
    lambda = 1;
end    
if nargin < 3 || isempty(AbsTol)
    AbsTol = 1e-4;
end
if nargin < 4 || isempty(RelTol)
    RelTol = 1e-2;
end
if nargin < 5
    section = 10;
end

while (num_problem < 1) || (num_problem > 14)
    num_problem = input('Choose a problem in the range 1,...,14: ');
end

if (num_problem == 2) && (lambda ~= 1)
    warning('TestInf:ResetLambda','lambda reset to 1 for Problem 2')
    lambda = 1;
end

% Compute and plot the solution:
[soln,errest,cond] = Fie(lambda,0,Inf,behavior,@kernel,@RHS,AbsTol,RelTol);
t = soln.s; sol = soln.x; 
nfinal = length(t);

% Interpolate the solution for assessing the error and a smooth graph.
tint = linspace(0,section,150);
xint = ntrpFie(soln,tint);
plot(tint,xint)
title(['Infinite interval problem ',num2str(num_problem),...
       ':   Final n = ',num2str(nfinal)])
xlabel(['a = 0,   b = \infty,   \lambda = ',...
        num2str(lambda),',   AbsTol = ',num2str(AbsTol),...
        ',   RelTol = ',num2str(RelTol)])

% Study the error and condition of the integral equation.
true_error = norm(true_soln(t) - sol,inf);
max_error = norm(true_soln(tint)-xint,inf);
disp(' ')
disp(['Condition number = ',num2str(cond)])
disp(' ')
fprintf('Approximate bound on error at nodes       = %6.1e\n',errest)   
fprintf('Actual error at nodes                     = %6.1e\n',true_error)
fprintf('\nActual error at 150 equally spaced points = %6.1e\n',max_error)

%====Nested functions=======================================================
function kst = kernel(s,t)
    switch num_problem
        case {1,3,4,5,12}
            % This is from paper of I. Sloan
            kst = 1./(1 + s.^2 + t.^2);
        case 2
            % Artificial problem with degenerate kernel.
            kst = exp(1 - t);
        case {6,8,10} 
            kst = 1./(1 + s.^2 + t.^2).^2;
        case {7,11} 
            kst = 1./(1 + s.^2 + t.^2).^3;
        case {9,13,14}
            kst = exp(-.2*(s+t).^2);
    end
end % kernel

function ts = true_soln(t)
    switch num_problem
        case 1
            ts = 1./(1 + t.^2);
        case 2
            ft = 1./((t-3).^2 + 1/2) + 1./((t-9).^2 + 1/3);
            ts = ft + (exp(1)/(1-exp(1)))*3.484412728591260e-001;
        case {3,10,13}
            ts = cos(t)./(1 + t.^2);
        case {4,11}
            ts = cos(t)./(1 + t.^2).^2;
        case {5,6,7,9}
            ts = t./(t+1);
        case {8,12,14}
            ts = cos(t).*exp(-.2*t);
    end
end % true_soln

function rs = RHS(s)
% The variable intgrl is the integral over [a,b] with
% respect to t of the integrand   kernel(s,t)*true_soln(t)
    switch num_problem
        case 1
            rs = zeros(size(s));
            rs(s == 0) = 1 - pi/4;
            ndx = find(s ~= 0);
            rs(ndx) = 1./(1 + s(ndx).^2) - ...
                    (pi./(2*s(ndx).^2)).*(1 - 1./sqrt(1 + s(ndx).^2));
        case 2
            rs = 1./((s-3).^2 + 1/2) + 1./((s-9).^2 + 1/3);
        case {3,4,5,12}
            num_int_nodes = 300;
            [tx,w] = intmi(0,1,num_int_nodes);
            x = (1-tx)./tx; 
            n = length(s);
            rs = zeros(size(s));
            for i=1:n
                fx = true_soln(x)./(1+s(i)^2+x.^2);
                rs(i) = sum(fx.*w./tx.^2);
            end
            rs = lambda*true_soln(s) - rs;
        case {6,8,10}
            num_int_nodes = 300;
            [tx,w] = intmi(0,1,num_int_nodes);
            x = (1-tx)./tx; 
            n = length(s);
            rs = zeros(size(s));
            for i=1:n
                fx = true_soln(x)./(1+s(i)^2+x.^2).^2;
                rs(i) = sum(fx.*w./tx.^2);
            end
            rs = lambda*true_soln(s) - rs;            
        case {7,11}
            num_int_nodes = 300;
            [tx,w] = intmi(0,1,num_int_nodes);
            x = (1-tx)./tx; 
            n = length(s);
            rs = zeros(size(s));
            for i=1:n
                fx = true_soln(x)./(1+s(i)^2+x.^2).^3;
                rs(i) = sum(fx.*w./tx.^2);
            end
            rs = lambda*true_soln(s) - rs;        
        case {9,13,14}
            num_int_nodes = 300;
            [tx,w] = intmi(0,1,num_int_nodes);
            x = (1-tx)./tx; 
            n = length(s);
            rs = zeros(size(s));
            for i=1:n
                fx = true_soln(x).*exp(-0.2*(s(i)+x).^2);
                rs(i) = sum(fx.*w./tx.^2);
            end
            rs = lambda*true_soln(s) - rs;        
    end
end % RHS

%=========================================================================
end % TestInf
