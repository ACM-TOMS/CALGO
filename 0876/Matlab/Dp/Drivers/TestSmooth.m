function TestSmooth(num_problem,lambda,a,b,AbsTol,RelTol,c)
% Test problems 1,...,8 for solving an integral equation of the form
%                        b
%  lambda*x(s) - Integral  K(s,t)*x(t) dt = RHS(s)
%                        a
% as s varies over [a,b].  The kernel K(s,t) is assumed to be moderately
% smooth on all of [a,b]x[a,b}, which corresponds to BEHAVIOR = 1.
%
% Some of the test problems place restrictions on [a,b]. Some require the 
% input of a parameter c. If c is not specified in the call list for one 
% of these problems, it will be requested by the program.  As the parameter 
% lambda approaches an element in the spectrum of the integral operator, 
% the problem becomes increasingly ill-conditioned. Accordingly, some of 
% the problems place restrictions on lambda. 
%
% The default tolerances are AbsTol = 1e-6 and RelTol = 1e-3.
%
% As an example, use     TestSmooth(8,1,-pi,pi)

behavior = 1;

if nargin < 5 || isempty(AbsTol)
    AbsTol = 1e-6;
end
if nargin < 6 || isempty(RelTol)
    RelTol = 1e-3;
end
while (num_problem < 1) || (num_problem > 8)
    num_problem = input('Choose a problem in the range 1,...,8: ');
end

if nargin < 7 || isempty(c)
    c = -1;
end
switch num_problem
    case {4,6,7}
        while ~(c > 0)
            c = input('Supply problem parameter c > 0: ');
        end
    case 5
        % The eigenvalues of this integral operator are c^j for
        % j = 0,1,... and -c^j for j = 1,2,... These values should be
        % avoided when choosing the value of lambda.    
        while ~((0 < c) && (c < 1))
            c = input('Supply c such that 0 < c < 1: ');
        end
end

% Test the input data:
switch num_problem
    case 1
        % The value of lambda can be anything other than 0 and b-a.  The 
        % equation should be ill-conditioned as lambda approaches either
        % of these forbidden choices.
        if (lambda == 0) || (lambda == b-a)
            error('For this problem, lambda must not be 0 or b-a.')
        end
    case 3
        if a <= -1/3
            error('For this problem, a must be larger than -1/3.')
        end  
    case {4,5,6}
        if ~(a == 0 && b == 1)
            warning('TestSmooth:ResetInterval','[a,b] reset to [0,1].')                        
            a = 0; b = 1;
        end
    case 7
        if ~(a == 0 && b == 1)
            warning('TestSmooth:ResetInterval','[a,b] reset to [0,1].')                                    
            a = 0; b = 1;        
        end     
        if (lambda == 0) || (lambda == 0.125)
            error('For this problem, lambda must not be 0 or 0.125.')
        end
    case 8
        if ~(a == -pi && b == pi)
            warning('TestSmooth:ResetInterval','[a,b] reset to [-pi,pi]].')                                    
            a = -pi; b = pi;
        end
        if lambda ~= 1
            warning('TestSmooth:ResetLambda','lambda reset to 1.')                        
            lambda = 1;
        end
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
title(['Smooth kernel problem ',num2str(num_problem),...
       ':   Final n = ',num2str(nfinal)])
if num_problem < 4 || num_problem == 8
    xlabel(['a = ',num2str(a),',   b = ',num2str(b),',   \lambda = ',...
        num2str(lambda),',   AbsTol = ',num2str(AbsTol),...
        ',   RelTol = ',num2str(RelTol)])
else
    xlabel(['a = ',num2str(a),',   b = ',num2str(b),',   \lambda = ',...
        num2str(lambda),',   AbsTol = ',num2str(AbsTol),...
        ',   RelTol = ',num2str(RelTol),',   c = ',num2str(c)])
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

%====Nested functions=======================================================
function kst = kernel(s,t)
    switch num_problem
        case 1
            % This is a trivial test case.  It should give an exact answer
            % with only rounding error in the solution.  
            kst = ones(size(s));
        case 2
            % For this problem, the Simpson quadrature error is constant
            % with respect to the choice of s.
            kst = s.*t.^2;
        case 3
            kst = 1 ./ (1+2*s+t);
        case 4
            % This is a standard case with the Runge kernel.  It is very
            % peaked along s = t when c is a small number.
            kst = c ./ (c*c+(s-t).^2);
        case 5
            % This problem has a periodic kernel with a periodic solution.
            kst = (1-c^2) ./ (1+c^2-2*c*cos(2*pi*(s+t)));
        case 6
            % This problem has a Green's kernel for the standard BVP for
            % the differential equation x''(s) = f(s), x(0) = x(b) = 0.
            [m,n] = size(s);
            kst = zeros(m,n);
            for i=1:m
                for j=1:n
                    if s(i,j) < t(i,j)
                        kst(i,j) = -s(i,j)*(1-t(i,j));
                    else
                        kst(i,j) = -t(i,j)*(1-s(i,j));
                    end
                end
            end
        case 7
            kst = (s.^3).*t.^4;
        case 8
            % Example 6 of Chap. 4 in S.G. Mikhlin and K.L. Smolitskiy, 
            % Approximate Methods for Solution of Differential and 
            % Integral Equations, Elsevier, London, 1967. This is the
            % plane interior Dirichlet problem for an ellipse. 
            kst = -(0.3/pi)./(1 - 0.64*cos( (s + t)/2 ).^2);
     end
end % kernel

function ts = true_soln(t)
    switch num_problem
        case 1
            ts = ones(size(t));
        case 2
            ts = t.^2;
        case 3
            ts = t;
        case 4
            ts = 0.06+t.*(-0.8+t);
        case 5
            ts = cos(6*pi*t);
        case 6
            ts = (t.^c).*(1-t);
        case 7
            ts = t.^c;
        case 8
            ts = 8.5 + (128/17)*cos(2*t);
    end
end % true_soln

function rs = RHS(s)
% The variable intgrl is the integral over [a,b] with
% respect to t of the integrand   kernel(s,t)*true_soln(t)
    switch num_problem
        case 1
            intgrl = (b-a)*ones(size(s));
        case 2
            intgrl = s*(b^5 - a^5)/5;
        case 3
            intgrl = (b-a) - (1+2*s).*log((1+2*s+b)./(1+2*s+a));
        case 4
            z1 = log((c^2+(1-s).^2)./(c^2+s.^2));
            z2 = atan((1-s)/c) + atan(s/c);
            intgrl = c + c*(s-0.4).*z1 + (0.06-c^2-0.8*s+s.^2).*z2;
        case 5
            intgrl = (c^3)*cos(6*pi*s);
        case 6
            z1 = (s.^(c+2))/(c+2);
            z2 = 1/(c+1) - s/(c+3);
            z3 = 2*s/((c+1)*(c+2)*(c+3));
            intgrl = z1.*z2 - z3;
        case 7
            intgrl = (s.^3)/(c + 5);
    end
    if num_problem == 8
        rs = 25 - 16*sin(s).^2;
    else
        rs = lambda*true_soln(s) - intgrl;
    end
end % RHS
%==========================================================================

end % TestSmooth