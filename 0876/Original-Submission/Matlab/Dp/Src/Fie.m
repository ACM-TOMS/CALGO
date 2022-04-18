function [sol,errest,cond] = Fie(lambda,a,b,behavior,kernel,RHS,AbsTol,RelTol)
% Solve the Fredholm integral equation
%                         b
%   lambda*x(s) - integral   K(s,t)x(t) dt = RHS(s)
%                         a
%
% INPUT
%
% LAMBDA non-zero scalar used to define the integral equation. 
%
% A < B are either both finite or A = 0 and B = Inf.
%
% BEHAVIOR The kernel K(s,t) is expected to be moderately smooth on all 
%          of R = [a,b]x[a,b] except possibly across tbe diagonal s = t. 
%          For finite A and B, BEHAVIOR specifies how the kernel behaves 
%          across the diagonal s = t. Several kinds of behavior are allowed:
%
%    1     K(s,t) is smooth. 
%          
%    2     K(s,t) is discontinuous in a low-order derivative.
%
%    3     K(s,t) = L(s,t)*log|s-t|.
%
%  ALPHA   K(s,t) = L(s,t)/|s-t|^alpha for 0 < alpha < 1. At present, we 
%          recommend choosing alpha <= 0.75.  
%
%          When [A,B] = [0,Inf), K(s,t) must be smooth across the diagonal
%          so the input value of BEHAVIOR is ignored.
%
% KERNEL Handle of a function for evaluating the kernel of the integral 
% equation.  If BEHAVIOR is 3 or ALPHA, the function is to evaluate L(s,t)
% and otherwise, K(s,t). In either case the function must accept matrices 
% s and t defining a grid.
%
% RHS Handle of a function for evaluating the right hand side of the 
% integral equation.  For a column vector s, it must return a column
% vector of the same size.
%
% OPTIONAL INPUT
%
% AbsTol, RelTol: Absolute and relative error tolerances. Default
% values are AbsTol = 1e-6, RelTol = 1e-3. The solver attempts to 
% find an approximate solution z for which 
%
%          norm(x - z,inf) <= atol + rtol*norm(z,inf)
%
% Here atol = max(AbsTol,0) and rtol = max(RelTol,100*eps).
%
%
% OUTPUT
%
% SOL is a solution structure containing information about the problem,
% values of an approximate solution on a mesh, an error estimate, and 
% quantities needed to evaluate an interpolant anywhere in [a,b]. 
%
% Often the mesh is sufficiently fine that values on the mesh give a 
% satisfactory graph of the solution.  These values are in the fields
%
% SOL.S  Nodes at which the approximate solution is computed.
%
% SOL.X  Approximate solution x(s) at corresponding entries of SOL.S.
%
% If the graph is not smooth or values are required at specific points,
% the approximate solution is evaluated with a call of the form  
%
%                  xint = ntrpFie(sol,sint). 
%
% The natural interpolant is evaluated at all the entries of a vector 
% SINT and returned as corresponding entries of XINT.  This function
% calls the functions kernel and RHS.  
%
% ERREST  An estimate of the maximum error at nodes.  If the number 
% of subdivisions of [a,b] exceeds the maximum allowed, FIE warns that
% the value of ERREST corresponding to the solution returned does not
% pass the error test.
%
% OPTIONAL OUTPUT
%
% COND  The integral equation is approximated by a system of linear
% equations Az = b. COND is an estimate of the condition of A in the 
% 1-norm.  b consists of samples of RHS(s), so a large value of COND 
% indicates that small changes in RHS(s) can lead to large changes in 
% x(s), i.e., the integral equation is ill-conditioned.
%
%
% Fie calls the functions alg_integrals and log_integrals.

if lambda == 0
    error('LAMBDA must be non-zero.')
end

if a >= b
    error('Must have A < B.')
end
if isinf(b) && a ~= 0
    error('When B = Inf, A must be 0.')
end

if nargin < 7 || isempty(AbsTol)
    atol = 1e-6;
else
    atol = max(AbsTol,0);
end
if nargin < 8 || isempty(RelTol)
    rtol = 1e-3;
else
    rtol = max(RelTol,100*eps);
end

if isinf(b)
    behavior = 0;
elseif (0 < behavior) && (behavior < 1)
    alpha = behavior;
    behavior = 4;
elseif (behavior ~= 1) && (behavior ~= 2) && (behavior ~= 3)
    error('BEHAVIOR must be 1,2,3 or ALPHA with 0 < ALPHA < 1.')
end 
    

% The main loop, with the number of subdivisions to be doubled
% repeatedly until a solution with sufficient accuracy is obtained.

% Create and initialize solution structure sol.
n = 2^3;
[told,solnold] = mod_simp(n);
if behavior == 0
    sol.s = (1 - told) ./ told;
else
    sol.s = told;
end
sol.x = solnold;
sol.lambda = lambda;
sol.kernel = kernel;
sol.behavior = behavior;
if behavior == 4
    sol.alpha = alpha;
end
sol.RHS = RHS;

max_m = 9;
for m = 4:max_m
    n = 2^m;
    [t,soln] = mod_simp(n);    % Solve with n subdivisions of [a,b].
    
    % Calculate the norm of the difference between the current numerical
    % solution and the preceding one. In several cases we must interpolate
    % to get approximations on the same mesh.
    switch behavior
        case 0
            t_inf = (1 - t) ./ t;
            delta = norm(soln - ntrpFie(sol,t_inf),inf);        
        case {1,2}
            delta = norm(soln(1:2:n+1) - solnold,inf);
        case {3,4}
            delta = norm(soln - ntrpFie(sol,t),inf);
    end

    % Update solution structure.
    solnold = soln;
    if behavior == 0
        sol.s = (1 - t) ./ t;
    else
        sol.s = t;
    end
    sol.x = soln;

    if m > 4
        % Estimate contraction rate and use to predict error.  
        % Be conservative and do not use if not contracting.
        if behavior == 0
            rate = 0.5;
        else
            rate = min(0.5,max(delta/deltaold,0.0625));
        end
        errest = (rate/(1 - rate))*delta;  
        if errest < (atol + rtol*norm(soln,inf))
            break;
        end
    end    
    deltaold = delta;
    if m == max_m
        warning('Fie:Failure','Failed to pass the error test.')
    end
end

% Inexpensive approximation of condition of matrix in 1-norm.
if nargout == 3
    cond = condest(kermat);
end


%===Nested functions=======================================================

function [t,soln] = mod_simp(n)
% Solve with n subdivisions of [a,b] using a quadrature rule.  Simpson's 
% rule is used for smooth kernels.  Product rules based on quadratic
% interpolation (like Simpson's rule) are used for kernels that are not
% smooth.  In addition, the mesh is graded at end points for integrable 
% singularities. Infinite intervals are transformed to [0,1] and a Gauss 
% formula is used to avoid evaluating at an end point.
    switch behavior
        case 0
            % Composite 2-point Gauss quadrature on [0,1].
            % Assumes n is even.
            idx = 1:2:n-1;
            t = zeros(n,1);
            t(idx)   = idx' - 1/sqrt(3);
            t(idx+1) = idx' + 1/sqrt(3);
            t = t/n;
        case {1,2}
            h = (b - a)/n;
            t = linspace(a,b,n+1)';
            % Weights for Simpson's rule.  Assumes n is even.
            wt = ones(1,n+1); wt(3:2:n-1) = 2; wt(2:2:n) = 4;
            wt = (h/3)*wt;
        case 3  % Use a graded mesh optimal for the max norm.
            t = gmesh(a,b,n,3)';
        case 4  % Use a graded mesh optimal for the L2 norm.       
            t = gmesh(a,b,n,3/(1.5 - alpha))';
    end
    
    [S,T] = meshgrid(t,t);    
    if behavior == 0
        rhs = RHS((1 - t)./t);
        Kvals = ( kernel((1 - S)./S,(1 - T)./T) ./ T.^2)'; 
        Msize = n;
    else
        rhs = RHS(t);
        Kvals = kernel(S,T)';
        Msize = n + 1;
    end
    
    kermat = zeros(Msize);
    switch behavior
        case 0  % Infinite interval, smooth across s = t.
            kermat = Kvals/n;    
        case 1  % Smooth across s = t.
            for col = 1:n+1
                kermat(:,col) = Kvals(:,col)*wt(col);
            end
        case 2  % Discontinous in a low-order derivative across s = t.
            for col = 1:n+1
                kermat(:,col) = Kvals(:,col)*wt(col);
            end
            tL = kernel(t,t-h/2);
            tR = kernel(t,t+h/2); 
            for i = 2:2:n
                kermat(i,i-1) = (h/6)*(Kvals(i,i-1) + 1.5*tL(i) - 0.5*tR(i));
                if i > 2
                    kermat(i,i-1) = kermat(i,i-1) + (h/3)*Kvals(i,i-1); 
                end
                kermat(i,i) = (h/6)*(3*tL(i) + 2*Kvals(i,i) + 3*tR(i));
                kermat(i,i+1) = (h/6)*(Kvals(i,i+1) + 1.5*tR(i) - 0.5*tL(i));
                if i < n
                    kermat(i,i+1) = kermat(i,i+1) + (h/3)*Kvals(i,i+1);   
                end
            end 
        case 3  % Behaves like log|s-t| near s = t.
            for j = 2:2:n
                h = t(j) - t(j-1);
                hlnh3 = h*log(h)/3;
                for i = 1:n+1
                    beta = (t(i) - t(j))/h;
                    ints = log_integrals(beta);
                    I1 =   hlnh3 + h*(ints(3) - ints(2))/2;
                    I2 = 4*hlnh3 + h*(ints(1) - ints(3));
                    I3 =   hlnh3 + h*(ints(3) + ints(2))/2;
                    kermat(i,j-1) = kermat(i,j-1) + I1;
                    kermat(i,j) = I2;
                    kermat(i,j+1) = I3;
                end
            end
            kermat = Kvals.*kermat;            
        case 4  % Behaves like 1/|s-t|^alpha near s = t.  
            for j = 2:2:n
                h = t(j) - t(j-1);
                h_power = 0.5*h^(1-alpha);
                for i = 1:n+1
                    beta = (t(i) - t(j))/h;
                    ints = alg_integrals(beta,alpha);
                    I1 =   h_power*(ints(3) - ints(2));
                    I2 = 2*h_power*(ints(1) - ints(3));
                    I3 =   h_power*(ints(3) + ints(2));
                    kermat(i,j-1) = kermat(i,j-1) + I1;
                    kermat(i,j) = I2;
                    kermat(i,j+1) = I3;
                end
            end
            kermat = Kvals.*kermat;
    end 
    
    for j = 1:Msize
        kermat(j,j) = kermat(j,j) - lambda;
    end
    soln = kermat \ (-rhs);
    
end % mod_simp

%== Nested function used for behavior = 3,4 ===============================
function t = gmesh(a,b,n,q)
% n must be an integer evenly divisible by 4.
    nhalf = n/2;
    temp = 0.5*(b - a)*((2/nhalf)*(0:(nhalf/2))).^q;
    temp2(1:2:nhalf+1) = temp;
    temp2(2:2:nhalf) = 0.5*(temp(1:nhalf/2)+temp(2:nhalf/2+1));
    t = [ (a + temp2), (b - temp2(end-1:-1:1)) ];
end % gmesh  

%==========================================================================

end % Fie
