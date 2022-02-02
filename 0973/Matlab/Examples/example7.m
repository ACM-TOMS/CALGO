function example7(param,intquad)

%EXAMPLE7 Integrands that are more expensive to evaluate
%   This function requires Walter Gautschi's subroutines: 
%       r_jacobi, r_mod, stieltjes, lanczos, mcdis, quadrat, gauss, 
%       and gauss_rational 
%   that can be found on:
%       https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%   Input:
%       PARAM = positive integer corresponding to a matrix of size 2Nx2N, 
%               with N = 25 * 2^PARAM, in the definition of the integrand
%       INTQUAD = 1 for comparison with integral.m
%               = 2 for comparison with quadv.m
%   WARNING: The execution time of EXAMPLE7 grows rapidly for increasing
%            PARAM.


N = 25*2^param;
Na = 18;

disp('Integrands that are more expensive to evaluate:')

disp('We consider the integral int_{-1}^{1}f(x)dx, for a vector-valued')
disp('function of the form:')
disp('         f(x)=(A-x*I)^{-1}B,')
disp('where')
disp('         A is a random antisymmetric matrix of size 2Nx2N,') 
disp('         I is the 2Nx2N identity matrix,') 
disp('         B is a random column vector of length 2N, and')
disp(sprintf('         N=%1.0f,',N))
disp('which has complex conjugate singularities on the imaginary axis.')

disp('Figure 1 shows the ''exact'' maximal relative error against')
disp('execution time of MATLAB''s built-in automatic integrator integral.m')
disp('(respectively quadv.m) compared with respectively')
disp('1. the (2n+1)-point rational Fejer without error estimates,')
disp('2. the (2n+1)-point rational Fejer with error estimates,')
disp('3. the semi-automatic ( (3:2:2n+1)-points ) rational Fejer,')
disp('4. the (n+1)-point rational Gauss-Legendre, and')
disp('5. the (2:1:n+1)-points rational Gauss-Legendre,')
disp(sprintf('for n=1,...,%1.0f.',Na))
disp('Figure 2 shows the exact maximal relative errors against n together')
disp('with the estimated maximal relative error for the n-point')
disp('rational Fejer.')
disp('Figure 3 shows the imaginary part of the singularities of f(x) in')
disp('the upper half-plane.')
disp('The exact value of the integral is obtained from the MATLAB''s') 
disp('built-in automatic integrators themselves with highest possible')
disp('precision.')

% Integrand

A = randn(2*N);
A = A-A';
B = randn(2*N,1);
fun = @(x) (A-x*eye(2*N))\B;

% Computing the exact value

disp('Busy to compute the exact value')

if intquad==1,
    Q = integral(fun,-1,1,'RelTol',1e-15,'AbsTol',0,'ArrayValued',true);
else
    Q = quadv(fun,-1,1,1e-15);
end

% Initialisation

v = 12;
P1 = 10.^(-[1:v]);
Ei = zeros(1,v); 
Ti = Ei;

E = zeros(4,Na);
Es = zeros(1,Na);
T = E;

c = Na+v+1;

% Iterations for rational Fejer and rational Gauss-Legendre 

for j = 1:Na,
    
    disp(sprintf('Execution time measurements countdown: %1.0f', c-j));
        
    % 1. semi-automatic
    for k=1:10,
        tic;
        sgl1 = eig(A);
        sgl1 = sgl1(1:2:end);
        sgl1 = sort(sgl1,'ascend');
        [X1,e1,r1] = rfejer(sgl1(1:j),fun,'Array',true);
        t1(k) = toc;
    end
    T(1,j) = mean(t1);
    E(1,j) = max(max(abs(X1-Q)./abs(Q)));
    
    % 2. n-point rational Fejer with error estimates
    for k=1:10,
        tic;
        sgl2 = eig(A);
        sgl2 = sgl2(1:2:end);   
        sgl2 = sort(sgl2,'ascend');
        [X2,e2,r2] = rfejer(sgl2(1:j),fun,'Array',true,'Nmin',Na+1);
        t2(k) = toc;
    end
    T(2,j) = mean(t2);
    E(2,j) = max(max(abs(X2-Q)./abs(Q)));
    Es(j) = max(max(e2));
    
    % 3. n-point rational Fejer without error estimates
    for k=1:10,
        tic;
        sgl3 = eig(A);
        sgl3 = sgl3(1:2:end);
        sgl3 = sort(sgl3,'ascend');
        [XL,e3,r3] = rfejer(sgl3(1:j));
    
        X3 = zeros(size(fun(0)));
        for K = 1:length(XL), 
            X3 = X3+fun(XL(1,K))*XL(2,K);
        end
        t3(k) = toc;
    end
    T(3,j) = mean(t3);
    E(3,j) = max(max(abs(X3-Q)./abs(Q)));
    
    % 4. n-point rational Gauss-Legendre
    for k=1:10,
        tic;
        sgl4 = eig(A);
        sgl4 = sort(sgl4,'ascend');
        % warning('off', 'Octave:possible-matlab-short-circuit-operator');
        XLG = quadWGR2( sgl4(1:2*j),j+1,1000 );
        % warning('on', 'Octave:possible-matlab-short-circuit-operator');
    
        X4 = zeros(size(fun(0)));
        for KG = 1:length(XLG), 
            X4 = X4+fun(XLG(1,KG))*XLG(2,KG);
        end
        t4(k) = toc;
    end
    T(4,j) = mean(t4);
    E(4,j) = max(max(abs(X4-Q)./abs(Q)));
     
end

c = c-Na;

% Iterations for integral.m (respectively quadv.m)

if intquad==1,
    for j = 1:v,
        
        disp(sprintf('Execution time measurements countdown: %1.0f', c-j));
        
        p1 = P1(j);

    % 5. MATLAB's built-in automatic integrator integral.m
        for k = 1:10,
            tic;
            q = integral(fun,-1,1,'RelTol',p1,'AbsTol',0,'ArrayValued',true);
            t4(k) = toc;
        end
        Ti(j) = mean(t4);
        Ei(j) = max(max(abs(q-Q)./abs(Q)));
    end
else
    for j = 1:v,
        
        disp(sprintf('Execution time measurements countdown: %1.0f', c-j));
                
        p1 = P1(j);

    % 6. MATLAB's built-in automatic integrator quadv.m
        for k = 1:10,
            tic;
            q = quadv(fun,-1,1,p1);
            t4(k) = toc;
        end
        Ti(j) = mean(t4);
        Ei(j) = max(max(abs(q-Q)./abs(Q)));
    end
end

% Display results

disp('Figure 1 shows the ''exact'' maximal relative error against')
disp('execution time of MATLAB''s built-in automatic integrator integral.m')
disp('(respectively quadv.m) compared with respectively')
disp('1. the (2n+1)-point rational Fejer without error estimates,')
disp('2. the (2n+1)-point rational Fejer with error estimates,')
disp('3. the semi-automatic ( (3:2:2n+1)-points ) rational Fejer,')
disp('4. the (n+1)-point rational Gauss-Legendre, and')
disp('5. the (2:1:n+1)-points rational Gauss-Legendre,')
disp(sprintf('for n=1,...,%1.0f.',Na))
disp('Figure 2 shows the exact maximal relative errors against n together')
disp('with the estimated maximal relative error for the (n+1)-point')
disp('rational Fejer.')
disp('Figure 3 shows the imaginary part of the singularities of f(x) in')
disp('the upper half-plane.')
disp('The exact value of the integral is obtained from the MATLAB''s') 
disp('built-in automatic integrators themselves with highest possible')
disp('precision.')

for k=1:length(T(4,:)),
    Ts(k) = sum(T(4,1:k));
end

if intquad==1,
    figure
    loglog(T(3,:),E(3,:),'--',T(2,:),E(2,:),'-.',T(1,:),E(1,:),':',...
        T(4,:),E(4,:),'*-',Ts,E(4,:),'o-',Ti,Ei,'-');
    xlabel('time');
    ylabel('relative error');
    hleg1 = legend('RatFejNoEstim','RatFejWithEstim','SemiAutom',...
        'nRatGaussLeg','(2:1:n)RatGaussLeg','integral.m', ...
        'Location','NorthEastOutside');
    title('Maximal relative error against execution time');
else
    figure
    loglog(T(3,:),E(3,:),'--',T(2,:),E(2,:),'-.',T(1,:),E(1,:),':',...
        T(4,:),E(4,:),'*-',Ts,E(4,:),'o-',Ti,Ei,'-');
    xlabel('time');
    ylabel('relative error');
    hleg1 = legend('RatFejNoEstim','RatFejWithEstim','SemiAutom',...
        '(n+1)RatGaussLeg','(2:1:n+1)RatGaussLeg','quadv.m', ...
        'Location','NorthEastOutside');
    title('Maximal relative error against execution time');
end

figure
semilogy([3:2:2*Na+1],E(2,:),'--',[2:1:Na+1],E(4,:),'*-',...
    [3:2:2*Na+1],Es,'-.');
xlabel('n');
ylabel('relative error');
hleg1 = legend('Exact RatFej','Exact RatGaussLeg','Estimate RatFej');
title('Exact relative error versus estimated relative error')

sglr = eig(A);
sglr = sglr(1:2:end);
sglr = sort(sglr,'ascend');
figure
semilogy(imag(sglr),'*');
xlabel('k');
ylabel('imag(pole_k)');
title('Imaginary part of the singularities')

end

%--------------------------------------------------------------------------

function xw = quadWGR2( sgl , N, Nm )

%QUADWG Rational Gauss quadrature (Walter Gautschi) 
%   N : number of points in the rational Gauss quadrature formula
%   sgl : sequence of poles
%   xw(1,:) the nodes 
%   xw(2,:) the weights

clear mc mp iq idelta irout AB Z eps0 Nmax ab0 M
global mc mp iq idelta irout AB Z eps0 Nmax ab0 M

% constants
mc = 1; mp = 0; iq = 1; idelta = 2; irout = 1;  %irout = 0;
AB=[-1 1]; eps0 = 100 *eps; Nmax = Nm;

Z = ZmatrixR2(-1./sgl);
M = size(Z,1); 
ab0 = r_jacobi(Nmax);
[abmod,Ncap,kount]=r_mod(N,ab0);
xw = gauss_rational(N,abmod);
xw=real(xw)';

end


%--------------------------------------------------------------------------

function [ Zm ] = ZmatrixR2( a )

n= length(a);
k=1;
     while k <= n
          d = find(a == a(k));
          Zm(k,1) = a(k);
          Zm(k,2) = length(d);
          a(d(2:end))=[];
          m=single(Zm(k,2));
          n=n-m +1 ;
          k=k+1;
     end
end



