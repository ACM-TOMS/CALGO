function TestAlgSing(num_problem,lambda,a,b,alpha,AbsTol,RelTol)
% Test problems 1,...,4 for solving an integral equation of the form
%                        b
%  lambda*x(s) - Integral  L(s,t)/|s-t|^alpha dt = RHS(s)
%                        a
% as s varies over [a,b].  L(s,t) is assumed to be moderately smooth
% smooth on all of [a,b]x[a,b].  These problems have BEHAVIOR = ALPHA
% with 0 < ALPHA < 1.
%
% The default tolerances are AbsTol = 1e-6 and RelTol = 1e-3.
%
% As an example, use     TestAlgSing(3,5,0,1,0.5)
%
% TestAlgSing calls the functions alg_integrals and intmi.

behavior = alpha;
if (alpha < 0) || (1 < alpha)
    error('Must have 0 < ALPHA < 1.')
end

if nargin < 6
    AbsTol = 1e-6;
end
if nargin < 7
    RelTol = 1e-3;
end
while (num_problem < 1) || (num_problem > 4)
    num_problem = input('Choose a problem in the range 1,...,4: ');
end

% Test the input data:
switch num_problem
    case {1,2}
        % Use [a,b] = [-1,1].
        if ~(a == -1 && b == 1)
            warning('TestAlgSing:ResetInterval','[a,b] reset to [-1,1].')
            a = -1; b = 1;
        end
    case 3
        % Use [a,b] = [0,1] and alpha = 0.5.
        if ~(a == 0 && b == 1)
            warning('TestAlgSing:ResetInterval','[a,b] reset to [0,1].')            
            a = 0; b = 1;
        end
        if alpha ~= 0.5
            warning('TestAlgSing:ResetAlpha','alpha reset to 0.5.')            
            alpha = 0.5;
            behavior = alpha;
        end
    case 4
        % Use [a,b] = [0,pi/2] and alpha = 0.5.
        if ~(a == 0 && b == pi/2)
            warning('TestAlgSing:ResetInterval','[a,b] reset to [0,pi/2].')
            a = 0; b = pi/2;
        end
        if alpha ~= 0.5
            warning('TestAlgSing:ResetAlpha','alpha reset to 0.5.')            
            alpha = 0.5;
            behavior = alpha;
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
title(['Algebraic singularity problem ',num2str(num_problem),...
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
if num_problem == 3 || num_problem == 4
    disp('EVALUATING RHS IS QUITE EXPENSIVE FOR THIS EQUATION.')
    disp(' ')
end
disp(['Condition number = ',num2str(cond)])
disp(' ')
fprintf('Approximate bound on error at nodes       = %6.1e\n',errest)   
fprintf('Actual error at nodes                     = %6.1e\n',true_error)
fprintf('\nActual error at 150 equally spaced points = %6.1e\n',max_error)

%===Nested functions=======================================================
function kst = kernel(s,t)
    kst = ones(size(s));
    if num_problem == 4
        temp = s - t;
        ndx = find(temp ~= 0);         
        kst(ndx) = sqrt( abs( temp(ndx) ./ sin(temp(ndx)) ) );    
    end
end % kernel

function ts = true_soln(t)
    switch num_problem
        case 1
            ts = ones(size(t));
        case 2
            ts = t.^2;
        case {3,4}
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
            for i = 1:m
                c = s(i);
                intgrl(i) = ((1-c)^(1-alpha) + (1+c)^(1-alpha))/(1-alpha);
            end
        case 2
            % The true solution is x(t) = t^2. The error should be zero.
            m = length(s);
            for i = 1:m
                c = s(i);
                rs = alg_integrals(c,alpha);
                intgrl(i) = rs(3);
            end   
        case {3,4}
            % The true solution is x(s) = s*log(s). Use the intmi routine to 
            % generate the nodes and weights for numerical integration of these 
            % singular integrals.
            num_nodes = 100;
            [nodes_norm,weights_norm] = intmi(0,1,num_nodes);
            m = length(s);
            for i = 1:m
                si = (s(i)/2)^0.25;
                % Integrate over the interval 0 <= t <= s(i)/2.
                x = si*nodes_norm;
                w = si*weights_norm;
                tempsum = 0;
                for j = 1:num_nodes
                    if ((x(j) ~= si) && (x(j) ~= 0))
                        tempsum = tempsum + w(j)*(x(j)^7)*log(x(j))*...
                            kernel(s(i),x(j)^4)/sqrt(s(i)-x(j)^4);
                    end
                end
                tempsum1 = 16*tempsum;
                si = sqrt(s(i)/2);
                % Integrate over the interval s(i)/2 <= t <= s(i).
                x = si*nodes_norm;
                w = si*weights_norm;
                tempsum = 0;
                for j = 1:num_nodes
                    if ((x(j) ~= si) && (x(j) ~= 0))
                        tempsum = tempsum + w(j)*(s(i)-x(j)^2)*...
                            kernel(s(i),s(i)-x(j)^2)*log(s(i)-x(j)^2);
                    end
                end
                tempsum2 = 2*tempsum;
                si = sqrt(b-s(i));
                % Integrate over the interval s(i) <= t <= b.
                x = si*nodes_norm;
                w = si*weights_norm;
                tempsum = 0;
                for j = 1:num_nodes
                    if ((x(j) ~= si) && (x(j) ~= 0))
                        tempsum = tempsum + w(j)*(s(i)+x(j)^2)*...
                            kernel(s(i),s(i)+x(j)^2)*log(s(i)+x(j)^2);
                    end
                end
                tempsum3 = 2*tempsum;
                intgrl(i) = tempsum1 + tempsum2 + tempsum3;
            end  
    end
    rs = lambda*true_soln(s) - intgrl;
end % RHS
%==========================================================================

end % TestAlgSing
