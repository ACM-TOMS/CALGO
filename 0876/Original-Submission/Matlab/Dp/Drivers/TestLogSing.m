function TestLogSing(num_problem,lambda,a,b,AbsTol,RelTol)
% Test problems 1,...,5 for solving an integral equation of the form
%                        b
%  lambda*x(s) - Integral  L(s,t)*log|s-t| dt = RHS(s)
%                        a
% as s varies over [a,b].  L(s,t) is assumed to be moderately smooth
% smooth on all of [a,b]x[a,b].  These problems have BEHAVIOR = 3.
%
% The default tolerances are AbsTol = 1e-6 and RelTol = 1e-3.
%
% As an example, use     TestLogSing(5,1,0,2)
%
% TestLogSing calls the function intmi.

behavior = 3;

if nargin < 5 || isempty(AbsTol)
    AbsTol = 1e-6;
end
if nargin < 6 || isempty(RelTol)
    RelTol = 1e-3;
end
while (num_problem < 1) || (num_problem > 5)
    num_problem = input('Choose a problem in the range 1,...,5: ');
end

% Test the input data:
if num_problem == 5
    % The interval is [0,b].
    if a ~= 0
        error('For this problem, the interval must be [0,b].')
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
title(['Logarithmic singularity problem ',num2str(num_problem),...
       ':   Final n = ',num2str(nfinal)])
xlabel(['a = ',num2str(a),',   b = ',num2str(b),',   \lambda = ',...
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

%===Nested functions=======================================================
function kst = kernel(s,t)
    if num_problem == 4
        kst = t;
    else
        kst = ones(size(s));
    end
end % kernel

function ts = true_soln(t)
    switch num_problem
        case 1
            ts = ones(size(t));
        case {2,4}
            ts = t.^2;
        case 3
            ts = t.^3;
        case 5
            ts = zeros(size(t));
            m = length(t);
            for i = 1:m
                if t(i) == 0
                    ts(i) = 0;
                else
                    ts(i) = t(i)*log(t(i));
                end
            end
    end
end % true_soln

function rs = RHS(s)
% The variable intgrl is the integral over [a,b] with
% respect to t of the integrand kernel(s,t)*true_soln(t).
    intgrl = zeros(size(s));
    switch num_problem
        case 1
            % The true solution is x(t) = 1. The error should be zero.
            m = length(s);
            L = b - a;
            for i = 1:m
                if s(i) == a || s(i) == b
                    intgrl(i) = L*log(L) - L;
                else
                    intgrl(i) = (b-s(i))*log(b-s(i)) + ...
                                (s(i)-a)*log(s(i)-a) - L;
                end
            end
        case 2
            % The true solution is x(t) = t^2.
            m = length(s);
            L = b - a; 
            for i = 1:m
                if s(i) == a
                    lnL = log(L);
                    temp1 = L*lnL - L;
                    temp2 = (L^2)*lnL/2 - (L^2)/4;
                    temp3 = (L^3)*lnL/3 - (L^3)/9;
                    intgrl(i) = (a^2)*temp1 + 2*a*temp2 + temp3;
                elseif s(i) == b
                    lnL = log(L);
                    temp1 = L*lnL - L;
                    temp2 = -(L^2)*lnL/2 + (L^2)/4;
                    temp3 = (L^3)*lnL/3 - (L^3)/9;
                    intgrl(i) = (b^2)*temp1 + 2*b*temp2 + temp3;
                else
                    bs = b-s(i); sa = s(i)-a;
                    lnbs = log(bs); lnsa = log(sa);
                    temp1 = bs*lnbs + sa*lnsa - L;
                    temp2 = ((bs^2)*lnbs - (sa^2)*lnsa)/2 - ((bs^2)-(sa^2))/4;
                    temp3 = ((bs^3)*lnbs + (sa^3)*lnsa)/3 - ((bs^3)+(sa^3))/9;
                    intgrl(i) = (s(i)^2)*temp1 + 2*s(i)*temp2 + temp3;
                end
            end
        case {3,4}
           % The true solution is x(t) = t^3 and t^2, respectively.
            m = length(s);
            L = b - a; 
            for i = 1:m
                if s(i) == a
                    lnL = log(L);
                    temp1 = L*lnL - L;
                    temp2 = (L^2)*lnL/2 - (L^2)/4;
                    temp3 = (L^3)*lnL/3 - (L^3)/9;
                    temp4 = (L^4)*lnL/4 - (L^4)/16;
                    intgrl(i) = (a^3)*temp1 + 3*(a^2)*temp2 + 3*a*temp3 + temp4;
                elseif s(i) == b
                    lnL = log(L);
                    temp1 = L*lnL - L;
                    temp2 = -(L^2)*lnL/2 + (L^2)/4;
                    temp3 = (L^3)*lnL/3 - (L^3)/9;
                    temp4 = -(L^4)*lnL/4 + (L^4)/16;
                    intgrl(i) = (b^3)*temp1 + 3*(b^2)*temp2 + 3*b*temp3 + temp4;
                else
                    bs = b-s(i); sa = s(i)-a;
                    lnbs = log(bs); lnsa = log(sa);
                    temp1 = bs*lnbs + sa*lnsa - L;
                    temp2 = ((bs^2)*lnbs - (sa^2)*lnsa)/2 - ((bs^2)-(sa^2))/4;
                    temp3 = ((bs^3)*lnbs + (sa^3)*lnsa)/3 - ((bs^3)+(sa^3))/9;
                    temp4 = ((bs^4)*lnbs - (sa^4)*lnsa)/4 - ((bs^4)-(sa^4))/16;
                    intgrl(i) = (s(i)^3)*temp1 + 3*(s(i)^2)*temp2 + 3*s(i)*temp3 + temp4;  
                end
            end
        case 5
            % Use INTMI to generate the nodes and weights for
            % numerical integration of these singular integrals.
            num_nodes = 100;
            [nodes_norm,weights_norm] = intmi(0,1,num_nodes);
            m = length(s);
            for i = 1:m
                si = s(i);
                % Integrate over [a,si].
                x = a + (si-a)*nodes_norm;
                w = (si-a)*weights_norm;
                tempsum = 0;
                for j = 1:num_nodes
                    if ((x(j) ~= si) && (x(j) ~= 0))
                        tempsum = tempsum + w(j)*x(j)*log(x(j))*log(si-x(j));
                    end
                end
                % Integrate over [si,b].
                x = si + (b-si)*nodes_norm;
                w = (b-si)*weights_norm;
                for j = 1:num_nodes
                    if ((x(j) ~= si) && (x(j) ~= 0))
                        tempsum = tempsum + w(j)*x(j)*log(x(j))*log(x(j)-si);
                    end
                end
                intgrl(i) = tempsum;
            end
    end        
    rs = lambda*true_soln(s) - intgrl;
end % RHS
%==========================================================================

end % TestLogSing
