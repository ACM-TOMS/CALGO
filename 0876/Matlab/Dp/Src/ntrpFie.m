function xint = ntrpFie(sol,sint)
% Evaluate the Nystroem interpolant for a solution structure sol at the
% points SINT.  XINT has the same shape as SINT.  ntrpFie must have 
% access to the functions sol.kernel and sol.RHS used to compute sol. 
%
% ntrpFie calls the functions alg_integrals and log_integrals.

    iscol = (size(sint,2) == 1);
    if iscol, sint = sint'; end
    t = sol.s; x = sol.x;            % local variables
    behavior = sol.behavior;
    if behavior == 0
        n = length(t);     % The number of quadrature node points.
        rhs = sol.RHS(sint(:))';        
        sint_fin = 1 ./ (1 + sint);
        t_fin = 1 ./ (1 + t);        
        [S,T] = meshgrid(sint_fin,t_fin);
        S_inf = (1 - S) ./ S;
        T_inf = (1 - T) ./ T;
        kermat = ( sol.kernel(S_inf,T_inf) ./ T.^2 )'/n;             
    else
        rhs = sol.RHS(sint(:))';
        [S,T] = meshgrid(sint,t);
        kermat = sol.kernel(S,T)';
        n = length(t) - 1;   % The number of subdivisions.
    end
    
    switch sol.behavior
    % Solve with n subdivisions of [a,b] using a quadrature rule.  Simpson's 
    % rule is used for smooth kernels.  Product rules based on quadratic
    % interpolation (like Simpson's rule) are used for kernels that are not
    % smooth.  In addition, the mesh is graded at end points for integrable 
    % singularities. Infinite intervals are transformed to [0,1] and a Gauss 
    % formula is used to avoid evaluating at an end point.        
        case 0
            % Infinite interval, smooth across s = t.
            xint = (rhs(:) + kermat*x)' / sol.lambda;        
        case 1
            % Smooth across s = t.
            h = t(2) - t(1);      
            % Weights for Simpson's rule.  This assumes n is even.
            wt = ones(1,n+1); wt(3:2:n-1) = 2; wt(2:2:n) = 4;
            wtxsoln = (h/3)*wt .* x';
            xint = (rhs(:) + kermat*wtxsoln(:))' / sol.lambda; 
        case 2
            % Discontinous in a low-order derivative across s = t.
            h = t(2) - t(1);      
            % Weights for Simpson's rule.  This assumes n is even.
            wt = ones(1,n+1); wt(3:2:n-1) = 2; wt(2:2:n) = 4;
            wtxsoln = (h/3)*wt .* x';
            xint = zeros(size(sint));   
            D1 = diff(x); D2 = diff(D1);
            for i = 1:length(sint)
                tau = sint(i);
                j = floor((tau - t(1))/(2*h));
                if j == n/2, j = j-1; end
                j = 2*j+2;
                % tau is in [t(j-1),t(j+1)]. Interpolate on
                % this subinterval with a quadratic.
                submesh = [(tau + t(j-1))/2, tau, (tau + t(j+1))/2];
                K_vals = sol.kernel(tau*ones(1,3),submesh);
                
                % Simpson's rule on [t(j-1),tau].
                mu = (tau - t(j-1))/h; nu = mu/2;
                SR_left = x(j-1)*kermat(i,j-1) ...
                   + 4*quadratic(nu,x(j-1),D1(j-1),D2(j-1))*K_vals(1) ...
                   + quadratic(mu,x(j-1),D1(j-1),D2(j-1))*K_vals(2);
                SR_left = SR_left*(tau - t(j-1))/6;
                
                % Simpson's rule on [tau,t(j+1)].
                nu = 1 + mu/2;
                SR_right = quadratic(mu,x(j-1),D1(j-1),D2(j-1))*K_vals(2) ...
                   + 4*quadratic(nu,x(j-1),D1(j-1),D2(j-1))*K_vals(3) ...
                   + x(j+1)*kermat(i,j+1);
                SR_right = SR_right*(t(j+1) - tau)/6;
                
                % Simpson's rule on the whole interval.
                SR = dot(kermat(i,1:j-2),wtxsoln(1:j-2)) + SR_left ...
                    + SR_right + dot(kermat(i,j+2:n+1),wtxsoln(j+2:n+1));
                if j == 2
                    SR = SR + (h/3)*kermat(i,j+1)*x(j+1);
                elseif j == n
                    SR = SR + (h/3)*kermat(i,j-1)*x(j-1);
                else
                    SR = SR + (h/3)*kermat(i,j-1)*x(j-1)...
                            + (h/3)*kermat(i,j+1)*x(j+1);
                end
                
                xint(i) = (rhs(i) + SR)/sol.lambda;
            end  
        case 3  
            % Behaves like log|s-t| near s = t.
            xint = zeros(size(sint)); 
            for j = 2:2:n
                h = t(j) - t(j-1);
                hlnh3 = h*log(h)/3;
                for i=1:length(sint)
                    tau = sint(i);
                    beta = (tau-t(j))/h;
                    ints = log_integrals(beta);
                    I1 =   hlnh3 + h*(ints(3) - ints(2))/2;
                    I2 = 4*hlnh3 + h*(ints(1) - ints(3));
                    I3 =   hlnh3 + h*(ints(3) + ints(2))/2;
                    xint(i) = xint(i) + I1*x(j-1)*kermat(i,j-1)...
                              + I2*x(j)*kermat(i,j) + I3*x(j+1)*kermat(i,j+1);
                end
            end
            xint = (rhs + xint)/sol.lambda;
        case 4  
            % Behaves like 1/|s-t|^alpha near s = t.
            xint = zeros(size(sint)); 
            alpha = sol.alpha;            
            for j = 2:2:n
                h = t(j) - t(j-1);
                h_power = 0.5*h^(1-alpha);
                for i = 1:length(sint)
                    tau = sint(i);
                    beta = (tau - t(j))/h;
                    ints = alg_integrals(beta,alpha);
                    I1 =   h_power*(ints(3) - ints(2));
                    I2 = 2*h_power*(ints(1) - ints(3));
                    I3 =   h_power*(ints(3) + ints(2));
                    xint(i) = xint(i) + I1*x(j-1)*kermat(i,j-1)...
                           + I2*x(j)*kermat(i,j) + I3*x(j+1)*kermat(i,j+1);
                end
            end
            xint = (rhs + xint)/sol.lambda;                    
    end
    
    if iscol, xint = xint'; end              
    
%== Nested function used for behavior = 2 =================================
function q = quadratic(mu,d0,d1,d2)
    q = d0 + mu*(d1 + (mu - 1)*d2/2);
end % quadratic
%==========================================================================

end % ntrpFie
