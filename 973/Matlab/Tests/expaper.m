function [ N, X , T , an , ae, c, nk, re] = expaper( ex )

%EXPAPER Numerical examples considered in the paper
%   Input: 
%       EX = integer between 1 and 8, corresponding to Example 7.EX
%   Output:
%       N = vector containing the number of nodes and weights in the
%           quadrature formulae
%       X(1,:) = vector containing the exact (maximal) relative error as a 
%           function of N for RFEJER with one input argument
%       X(2,:) = vector containing the exact (maximal) relative error as a 
%           function of N for RFEJER with two input arguments
%   (the following outputs are not defined for EX == 2 and EX == 7)
%       T(1,:) = vector containing the average time (of 100 measurements) 
%           to construct the N-point quadrature formula
%       T(2,:) = vector containing the average time (of 100 measurements) 
%           to construct the k-point quadrature formulae for k=c:i:N
%       T(3,:) = vector containing the average time (of 100 measurements)
%           for the semi-automatic integrator
%       [AN,AE,C,Y] = RFEJER(SGL,FUN [,ARRAY]) where SGL and FUN are
%           respectively the sequence of poles and the integrand from
%           Example 7.EX
%       NK = length(Y(1,:))
%       RE = exact relative error on the obtained approximation(s) AN
%
%   Remarks:
%       * For EX == 1 a choice should be made between OMEG == 1.1 and     
%         OMEG == 1.001 by removing/adding `%' on line 36
%       * Possible messages of the form "warning: division by zero" can be
%         avoided by using the alternatives from lines 325-353 instead on 
%         lines 38, 52, 76, and 110-113

switch ex
    case 1, % Example 7.1
        omeg = 1.1; 
        %omeg = 1.001;     % remove/add `%' at the beginning of this line
        fun = @(x)(pi.*x./omeg)./sin(pi.*x./omeg);
        % fun = @(x) sinf(x,omeg); % use this instead of the previous line
                                   % to avoid a warning message of the form
                                   % "warning: division by zero"
        sgl1=zeros(1,100);
        sgl1(1:2:end) = omeg*[1:1:50];
        sgl1(2:2:end) = -sgl1(1:2:end);
        if (omeg == 1.1),
            exact = 4.4677736463877657892361233985034;
        else % omeg = 1.001
            exact = 12.929256850002296208119060780724; 
        end
    case 2, % Example 7.2
        omeg = 1.1;
        fun = @(x)(pi.*x./omeg)./sin(pi.*x./omeg);
        % fun = @(x) sinf(x,omeg); % use this instead of the previous line
                                   % to avoid warning messages of the form
                                   % "warning: division by zero"
        sgl1=zeros(5,100); 
        sgl1(1,1:2:end) = omeg*[1:1:50];    % optimal set
        sgl1(1,2:2:end) = -sgl1(1,1:2:end);
        
        del = 1e-10; % small perturbations
        sgl1(2,:) = (1+del)*sgl1(1,:);      % set 1
        sgl1(3,:) = (1-del)*sgl1(1,:);      % set 2

        rho = omeg;
        sgl1(4,1:2:end) = rho;              % set 3
        sgl1(4,2:2:end) = -rho;
        
        rho = 1.3;
        sgl1(5,1:2:end) = rho;              % set 4
        sgl1(5,2:2:end) = -rho;        
        
        sgl1(6,:) = inf;                    % set 5
        exact = 4.4677736463877657892361233985034;
    case 3, % Example 7.3
        omeg = 0.1;
        fun = @(x)(pi.*x./omeg)./sinh(pi.*x./omeg);
        % fun = @(x) sihf(x,omeg); % use this instead of the previous line
                                   % to avoid warning messages of the form
                                   % "warning: division by zero"
        sgl1=zeros(1,100);
        sgl1(1:2:end) = -sqrt(-1)*omeg*[1:1:50];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = 0.15707963267939592622470892489958;
    case 4, % Example 7.4
        omeg = 1.1;
        fun = @(x) sin(1./(x-omeg));
        sgl1 = omeg*ones(1,100);
        exact = -1.19245706732219214078770039225;
    case 5, % Example 7.5
        omeg = 0.1;
        fun = @(x) 1./(exp(pi.*x./omeg)+1);
        sgl1 = zeros(1,100);
        sgl1(1:2:end) = sqrt(-1)*omeg*[1:2:100];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = 1;
    case 6, % Example 7.6
        fun = @(x) 1./sqrt((x+3).*(x+2));
        sgl1 = -2.5*ones(1,100);
        exact = 0.87116861981054736678487404465635;
    case 7, % Example 7.7
        fun = @(x) exp(x)./(25*x.^2+1);
        sgl1 = [.2*sqrt(-1)*[1,-1],inf*ones(1,98)];
        sgl2 = inf*ones(1,100);
        exact = 0.579408643011693305029559997358;
    otherwise % Example 7.8
        omeg = [2, 1.1, 1.01, 1.001];
        fun = @(x) [(pi.*x./omeg(1))./sin(pi.*x./omeg(1)), ...
            (pi.*x./omeg(2))./sin(pi.*x./omeg(2)); ...
            (pi.*x./omeg(3))./sin(pi.*x./omeg(3)), ...
            (pi.*x./omeg(4))./sin(pi.*x./omeg(4))];
        % fun = @(x) [sinf(x,omeg(1)), ... % use this instead of the 
        %     sinf(x,omeg(2)); ...         % previous lines to avoid
        %     sinf(x,omeg(3)), ...         % warning messages of the form
        %     sinf(x,omeg(4))];            % "warning: division by zero"       
        sgl1 = zeros(1,102);
        sgl1(1:6:end) = -omeg(4)*[1:1:17];
        sgl1(3:6:end) = -omeg(3)*[1:1:17];
        sgl1(5:6:end) = -omeg(2)*[1:1:17];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = [2.33248723224655024110707565174 , ...
            4.4677736463877657892361233985034; ...
                 8.43018458047084205897126420479  , ...
                 12.929256850002296208119060780724];
end

N = [2:1:41];
n = length(N);
T = zeros(3,n);
X = zeros(2,n);
 
if (ex < 7) && ~(ex == 2), % Examples 7.1, 7.3, 7.4, 7.5, and 7.6
    
    % exact relative errors and execution times 
    
    for M = 1:1:n,
        disp(sprintf('Execution time measurements countdown: %1.3f', n-M+1));
        sgl = sgl1(1:N(M)-1);
        for K=1:100,
            tic;
            [x,Aerr] = rfejer(sgl);
            A=fun(x(1,:))*x(2,:)';
            t1(K) = toc;
            tic;
            [an1,ae1,c1,y1] = rfejer(sgl,fun);
            t2(K) = toc;
        end
        T(1,M) = mean(t1);
        T(3,M) = mean(t2);
        X(1,M) = abs(exact-A)/abs(exact);
        X(2,M) = abs(exact-an1)/abs(exact);
        
        if M == 1,
            T(2,1) = T(1,1);
        elseif (ex == 7), % c = 3, i = 1
            T(2,M) = sum(T(1,2:M));
        elseif (ex == 3) || (ex == 5), % c = 3, i = 2
            T(2,M) = sum(T(1,2:2:M));
        else
            T(2,M) = sum(T(1,1:M));
        end
    end 
    index = find(X==0);
    X(index)=eps/10;
    
    % convergence / deterioration
    
    [an,ae,c,y] = rfejer(sgl1,fun);
    nk = length(y(1,:));
    re = abs(exact-an)/abs(exact);

elseif (ex == 8), % Example 7.8
    
    % exact relative errors and execution times 
    
    for M = 1:1:n,
        disp(sprintf('Execution time measurements countdown: %1.3f', n-M+1));
        sgl = sgl1(1:N(M)-1);
        for K=1:100,
            tic;
            [x,Aerr] = rfejer(sgl);
            F = fun(x(1,:));
            m = length(x(1,:));
            for j = 1:2,
                A(:,j) = F(:,(j-1)*m+1:j*m)*x(2,:)';
            end
            t1(K) = toc;
            tic;
            [an1,ae1,c1,y1] = rfejer(sgl,fun,'Array',true);
            t2(K) = toc;
        end
        T(1,M) = mean(t1);
        T(2,M) = sum(T(1,1:M));
        T(3,M) = mean(t2);
        Xa = abs(exact-A)./abs(exact);
        X(1,M) = max(Xa(:));
        Xa = abs(exact-an1)./abs(exact);
        X(2,M) = max(Xa(:));
    end 
    index = find(X==0);
    X(index)=eps/10;
    
    % convergence / deterioration
    
    [an,ae,c,y] = rfejer(sgl1,fun,'Array',true);
    nk = length(y(1,:));
    re = abs(exact-an)./abs(exact);
    
elseif (ex == 2), % Example 7.2
    
    X = inf*ones(5,n);
    T = NaN; an = NaN; ae = NaN; c = NaN; nk = NaN; re = NaN; 
    
    % n-point rational Fejer

    for M = 1:1:n,
        % optimal set
        sgl = sgl1(1,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(1,M) = abs(exact-A)/abs(exact);
        
        % set 1
        sgl = sgl1(2,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(2,M) = abs(exact-A)/abs(exact);
        
        % set 2
        sgl = sgl1(3,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(3,M) = abs(exact-A)/abs(exact);
        
        % set 3
        sgl = sgl1(4,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(4,M) = abs(exact-A)/abs(exact);        
        
        % set 4 
        sgl = sgl1(5,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(5,M) = abs(exact-A)/abs(exact);       
        
        % set 5 
        sgl = sgl1(6,1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(6,M) = abs(exact-A)/abs(exact);  
            
    end 
    index = find(X==0);
    X(index)=eps/10;
    X1 = X(1,1:19); X2 = X(2,1:19); X3 = X(3,1:19);
    N1 = N(1:19);
    X4 = X(4,:); X5 = X(5,:); X6 = X(6,:);
    figure;
    semilogy(N1,X1,'-',N1,X2,'--',N1,X3,'-.',N,X4,'-+',N,X5,'--o',N,X6,'-o');
    xlabel('n');
    ylabel('relative error');
    hleg1 = legend('opt','Delta = 1e-10','Delta = - 1e-10','rho = omega',...
        'rho = 1.3', 'infinity');
else % Example 7.7
    X = inf*ones(2,n);
    T = NaN; an = NaN; ae = NaN; c = NaN; nk = NaN; re = NaN; 
    
    % n-point rational Fejer

    for M = 1:1:n,
        % rational
        sgl = sgl1(1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(1,M) = abs(exact-A)/abs(exact);
        
        % polynomial
        sgl = sgl2(1:N(M)-1);
        [x,Aerr] = rfejer(sgl);
        A = fun(x(1,:))*x(2,:)';
        X(2,M) = abs(exact-A)/abs(exact);
    end 
    index = find(X==0);
    X(index)=eps/10;
    X1 = X(1,2:end); N1=N(2:end);
    X2 = X(2,:); 
    figure;
    semilogy(N1,X1,'-',N,X2,'--');
    xlabel('n');
    ylabel('relative error');
    hleg1 = legend('rational','polynomial');
end

if ~(ex==2) && ~(ex==7),
    
    if ~(ex==3) && ~(ex==5),
        p = [1:1:n];
    else
        p = [2:2:n];
    end
        
    figure;
    semilogy(N(p),X(1,p),'-',N(p),X(2,p),'--');
    xlabel('n');
    ylabel('Relative error');
    title('Relative error against number of points');
    
    figure;
    loglog(T(3,p),X(2,p));
    xlabel('t');
    ylabel('Relative error');
    title('Relative error against execution time');  
    
    T1 = T(1,p); T2 = T(2,p); T3 = T(3,p); N1=N(p);
    figure;
    semilogy(N1,T3,'-',N1,T2,'--',N1,T1,'-.');
    xlabel('n');
    ylabel('t');
    hleg1 = legend('t_{n}^{(sa)}','t_{c:i:n}','t_n');
    title('Execution time against number of points');

end

end

%--------------------------------------------------------------------------

function [ f ] = sinf( x , w )
%SINF @(x)(pi.*x./omeg)./sin(pi.*x./omeg)

warning('off'); % division by zero is allowed in the next line
f = (pi.*x./w)./sin(pi.*x./w);
warning('on');
index = find(x==0);
if ~isempty(index),
    f(index)=1;
end

end

%--------------------------------------------------------------------------

function [ f ] = sihf( x , w )
%SINF @(x)(pi.*x./omeg)./sinh(pi.*x./omeg)

warning('off'); % division by zero is allowed in the next line
f = (pi.*x./w)./sinh(pi.*x./w);
warning('on');
index = find(x==0);
if ~isempty(index),
    f(index)=1;
end

end
