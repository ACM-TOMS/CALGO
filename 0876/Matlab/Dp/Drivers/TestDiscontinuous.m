function TestDiscontinuous(num_problem,lambda,a,b,AbsTol,RelTol,c)
% Test problems 1,...,5 for solving an integral equation of the form
%                        b
%  lambda*x(s) - Integral  K(s,t)*x(t) dt = RHS(s)
%                        a
% as s varies over [a,b].  The kernel K(s,t) is assumed to be moderately
% smooth on all of [a,b]x[a,b} except across the diagonal s = t where
% there is a discontinuity in a low-order derivative, which corresponds
% to BEHAVIOR = 2.  This behavior is typical of a Green's function.
% 
% All the test problems require [a,b] to be [0,1] except #3 which allows
% an interval of the form [0,b].  Problems #4 and #5 require lambda = 1. 
% Problem #1 requires the input of a parameter c. If c is not specified in 
% the call list for this problem, it will be requested by the program. 
%
% The default tolerances are AbsTol = 1e-6 and RelTol = 1e-3.
%
% As an example, use     TestDiscontinuous(3,1,0,5)

behavior = 2;

if nargin < 5 || isempty(AbsTol)
    AbsTol = 1e-6;
end
if nargin < 6 || isempty(RelTol)
    RelTol = 1e-3;
end
while (num_problem < 1) || (num_problem > 5)
    num_problem = input('Choose a problem in the range 1,...,5: ');
end
if num_problem == 1 && nargin < 7
    c = input('Supply problem parameter c > 0: ');
end

% Test the input data:
if num_problem == 3
    if a ~= 0
       warning('TestDiscontinuous:ResetInterval','a reset to 0.')                    
       a = 0;
    end
elseif ~(a == 0 && b == 1)
    warning('TestDiscontinuous:ResetInterval','[a,b] reset to [0,1].')                    
    a = 0; b = 1;
end
if (num_problem == 4 || num_problem == 5) && lambda ~= 1
    warning('TestDiscontinuous:ResetLambda','lambda reset to 1.')                        
    lambda = 1;
end

% Compute and plot the solution:
[soln,errest,cond] = Fie(lambda,a,b,behavior,@kernel,@RHS,AbsTol,RelTol);
t = soln.s; sol = soln.x; 
nfinal = length(t) - 1;

% Interpolate the solution for assessing the error and if necessary,
% use to get a smooth graph.
tint = linspace(a,b,150);
xint = ntrpFie(soln,tint);
if nfinal < 150
    plot(tint,xint)
else
    plot(t,sol)
end
title(['Discontinuous kernel problem ',num2str(num_problem),...
       ':   Final n = ',num2str(nfinal)])
if num_problem == 1
    xlabel(['a = ',num2str(a),',   b = ',num2str(b),',   \lambda = ',...
        num2str(lambda),',   AbsTol = ',num2str(AbsTol),...
        ',   RelTol = ',num2str(RelTol),',   c = ',num2str(c)])
else
    xlabel(['a = ',num2str(a),',   b = ',num2str(b),',   \lambda = ',...
        num2str(lambda),',   AbsTol = ',num2str(AbsTol),...
        ',   RelTol = ',num2str(RelTol)])
end

% Study the error and condition of the integral equation.
true_error = norm(true_soln(t) - sol,inf);
max_error = norm(true_soln(tint)-xint,inf);
disp(' ')
disp(['Condition number = ',num2str(cond)])
disp(' ')
fprintf('Approximate bound on error at nodes       = %6.1e\n',errest)   
fprintf('Actual error at nodes                     = %6.1e\n',true_error)
fprintf('\nActual error at 150 equally spaced points = %6.1e\n',max_error)
   
%===Nested functions=======================================================
function kst = kernel(s,t)
    switch num_problem
        case {1,2,3}
            % This problem has a Green's kernel for the standard BVP for
            % the differential equation x''(s) = f(s), x(0) = x(b) = 0.
            % NOTE Problems 1 and 2 have b = 1.
            [m,n] = size(s);
            kst = zeros(m,n);
            for i=1:m
                for j=1:n
                    if s(i,j) < t(i,j)
                        kst(i,j) = -s(i,j)*(b-t(i,j));
                    else
                        kst(i,j) = -t(i,j)*(b-s(i,j));
                    end
                end
            end
        case 4
            % Example on p. 178 of D.F. Mayers, Chap. 14 in L. Fox, 
            % Numerical Solution of Ordinary and Partial Differential 
            % Equations, Pergamon Press, London, 1962.
            kst = abs(s - t);     
        case 5
            % Example of documentation for NAG code D05AAF.
            [m,n] = size(s);
            kst = zeros(m,n);
            for i=1:m
                for j=1:n
                    if s(i,j) < t(i,j)
                        kst(i,j) = s(i,j)*(1-t(i,j));
                    else
                        kst(i,j) = t(i,j)*(1-s(i,j));
                    end
                end
            end
    end
end % kernel

function ts = true_soln(t)
    switch num_problem
        case 1
            ts = (t.^c).*(1 - t);
        case {2,3}
            ts = exp(-t).*sin(t);
        case 4
            ts = t;
        case 5
            ts = sin(pi*t);
    end
end % true_soln

function rs = RHS(s)
% The variable intgrl is the integral over [a,b] with
% respect to t of the integrand   kernel(s,t)*true_soln(t)
    switch num_problem
        case 1
            z1 = (s.^(c+2))/(c+2);
            z2 = 1/(c+1) - s/(c+3);
            z3 = 2*s/((c+1)*(c+2)*(c+3));
            intgrl = z1.*z2 - z3;
        case 2
            intgrl = (s*(1-cos(1)*exp(-1)) - 1 + cos(s).*exp(-s))/2;
        case 3
            intgrl = (s*(1-cos(b)*exp(-b)) - b + b*cos(s).*exp(-s))/2;
    end
    switch num_problem
        case 4
            rs = -(2*s.^3 - 9*s + 2)/6;
        case 5
            rs = (1 - 1/pi^2)*sin(pi*s);
        otherwise
            rs = lambda*true_soln(s) - intgrl;
    end
end % RHS
%==========================================================================

end % TestDiscontinuous
