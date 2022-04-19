function U = LT_samples (NXYval,XY,NOPTS,S,tol,thrds)
%%   ***   USER DEFINED FUNCTION   ***   %%%%%%%%%%%%%%%%%%%%
% Solve for each point on the Talbot contour the BVP:
%
%       \Delta U - s*U = - u(x,y,0),    (x,y) in [0,1]x[0,1]
%                      = x*(1-x) + y*(1-y) initial condition
%
%       U(0,y,s) = 4/s^2 + y*(y-1)/s
%       U(1,y,s) = 4/s^2 + y*(y-1)/s
%       U(x,0,s) = 4/s^2 + x*(x-1)/s
%       U(x,1,s) = 4/s^2 + x*(x-1)/s
%
% which is a Helmholtz equation with Dirichlet boundary conditions.
%
% INPUT PARAMETERS
% XY    matrix with 2 columns XY = [X Y] where the two arrays contain the
%       cartesian coordinates of vectorized internal mesh points.
%       X and Y are arrays of length (NXYval-2)^2.
%
% OUTPUT PARAMETERS
% U     U is a matrix of size (NXYval-2)^2 x NOPTS, where
%       NOPTS = numel(S).
%       Each column contains the solution U(X,Y) solution at internal mesh points.
%       Each column is computed by solving a linear system.

    %% (opposite) initial condition
    u0 = @(x,y) x.*(1-x) + y.*(1-y); % = -u(x,y,0) = - [x*(x-1) + y*(y-1)]

    %% boundary conditions
    Uy = @(y,s) 4./s.^2 + y.*(y-1)./s; % U(0,y)=U(1,y)
    Ux = @(x,s) 4./s.^2 + x.*(x-1)./s; % U(x,0)=U(x,1)


    %% INTERNAL MESH POINTS: P(i,j),  i,j = 2,...,N-1
    n = NXYval-2;
    h = 1/(NXYval-1); % mesh step
    X = XY(:,1);
    Y = XY(:,2);

    %% BUILD THE LINEAR SYSTEM A*v = b for each S(K)
    U = zeros(numel(X),NOPTS); % preallocate

    parfor (K=1:NOPTS, thrds) % PARALLEL ON thrds THREADS
        s=S(K);
        d = -(4+s*h^2); % diagonal element
        A = sparse(eye(numel(X)))*d;
        b = h^2*u0(X,Y); % constant terms of the linear system -h^2*u(x,y,0)
        %% ultra-internal knots
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        P(i-1,j)=U(k-n)                        %
        %                                |                              %
        %                                |                              %
        %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k) ---- P(i,j+1)=U(k+1)     %
        %                                |                              %
        %                                |                              %
        %                        P(i+1,j)=U(k+n)                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        for j=3:n       % mesh point index x(j)
            for i=3:n   % mesh point index y(i)
                k=(j-2)*n+(i-2)+1; % index of unknown and of row in A
                A(k,k-1)=1; % U(k-1)
                A(k,k+1)=1; % U(k+1)
                A(k,k-n)=1; % U(k-n)
                A(k,k+n)=1; % U(k+n)
            end
        end
        %%  1) LEFT AND UPPER SIDES
        i=2;
        for j=3:n
            % 1.1) left side (without corner points)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 P(i-1,j)=U(k-n)                               %
            %                         |                                     %
            %                         |                                     %
            %                   P(i,j)=U(k) ---- P(i,j+1)=U(k+1)            %
            %                         |                                     %
            %                         |                                     %
            %                 P(i+1,j)=U(k+n)                               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            k=(j-2)*n+(i-2)+1; % index of unknown and of row in A
            A(k,k+1)=1; % U(k+1)
            A(k,k-n)=1; % U(k-n)
            A(k,k+n)=1; % U(k+n)
            b(k) = b(k) - Uy(Y(k),s); % constant term
            %  1.2) upper side (without corner points)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k) ---- P(i,j+1)=U(k+1)     %
            %                                |                              %
            %                                |                              %
            %                        P(i+1,j)=U(k+n)                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            k=(i-2)*n+(j-2)+1; % swapped i and j
            A(k,k-1)=1; % U(k-1)
            A(k,k+1)=1; % U(k+1)
            A(k,k+n)=1; % U(k+n)
            b(k) = b(k) - Ux(X(k),s); % constant term
        end
        %%  2) LOWER AND RIGHT SIDES
        j=n+1;
        for i=3:n
            % 2.1) lower side (without corner points)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                        P(i-1,j)=U(k-n)                        %
            %                                |                              %
            %                                |                              %
            %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k) ---- P(i,j+1)=U(k+1)     %
            %                                                               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            k=(j-2)*n+(i-2)+1;
            A(k,k-1)=1; % U(k-1)
            A(k,k+1)=1; % U(k+1)
            A(k,k-n)=1; % U(k-n)
            b(k) = b(k) - Ux(X(k),s); % constant term
            %  2.2) right side (without corner points)
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                        P(i-1,j)=U(k-n)                        %
            %                                |                              %
            %                                |                              %
            %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k)                          %
            %                                |                              %
            %                                |                              %
            %                        P(i+1,j)=U(k+n)                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            k=(i-2)*n+(j-2)+1; % swapped i and j
            A(k,k-1)=1; % U(k-1)
            A(k,k-n)=1; % U(k-n)
            A(k,k+n)=1; % U(k+n)
            b(k) = b(k) - Uy(Y(k),s);
        end
        %%  3) CORNER POINTS
        b(1)         = b(1)         - Uy(Y(1),s) - Ux(X(1),s);
        b(n)         = b(n)         - Uy(Y(n),s) - Ux(X(n),s);
        b((n-1)*n+1) = b((n-1)*n+1) - Uy(Y((n-1)*n+1),s) - Ux(X((n-1)*n+1),s);
        b(end)       = b(end)       - Uy(Y(end),s) - Ux(X(end),s);
        %%  3.1) up and left vertex
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                   P(i,j)=U(k) ---- P(i+1,j)=U(k+1)            %
        %                         |                                     %
        %                         |                                     %
        %                 P(i+1,j)=U(k+n)                               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        i=2; j=2; k=(j-2)*n+(i-2)+1;
            A(k,k+1)=1; % U(k+1)
            A(k,k+n)=1; % U(k+n)
        %%  3.2) down and left vertex
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 P(i-1,j)=U(k-n)                               %
        %                         |                                     %
        %                         |                                     %
        %                   P(i,j)=U(k) ---- P(i+1,j)=U(k+1)            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        i=2; j=n+1; k=(j-2)*n+(i-2)+1;
            A(k,k+1)=1; % U(k+1)
            A(k,k-n)=1; % U(k-n)
        %%  3.3) up and right vertex
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k)                          %
        %                                |                              %
        %                                |                              %
        %                        P(i+1,j)=U(k+n)                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        i=n+1; j=2; k=(j-2)*n+(i-2)+1;
            A(k,k-1)=1; % U(k-1)
            A(k,k+n)=1; % U(k+n)
        %%  3.4) down and right vertex
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        P(i-1,j)=U(k-n)                        %
        %                                |                              %
        %                                |                              %
        %     P(i,j-1)=U(k-1) ---- P(i,j)=U(k)                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        i=n+1; j=n+1; k=(j-2)*n+(i-2)+1;
            A(k,k-1)=1; % U(k-1)
            A(k,k-n)=1; % U(k-n)
    
        %% SOLVE THE SYSTEM for S(K)
        U(:,K) = A\b;
    end % parfor
    U = U.'; % transpose matrix (row-wise)
end

