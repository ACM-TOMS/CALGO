function [ N, X , T , E] = expaperWG( ex )

%EXPAPERWG Numerical examples considered in the paper
%   This function requires Walter Gautschi's subroutines: 
%       r_jacobi, r_mod, stieltjes, lanczos, mcdis, quadrat, gauss, 
%       and gauss_rational 
%   that can be found on:
%       https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%   Input: 
%       EX = integer between 1 and 8 (with exclusion of 2 and 7),
%            corresponding to Example 7.EX
%   Output:
%       N = vector containing the number of nodes and weights in the
%           rational Gauss quadrature formulae
%       X = vector containing the exact (maximal) relative error as a 
%           function of N
%       T = vector containing the average time (of 10 measurements)
%           to construct the N-point rational Gauss quadrature formula (WG)
%       E = vector containing estimates for the maximale error on the
%           computed weights as a function of N
%
%   Remark:
%       * For EX == 1 a choice should be made between OMEG == 1.1 and     
%         OMEG == 1.001 by removing/adding `%' on line 29 
    
switch ex
    case 1, % Example 7.1
        omeg = 1.1; 
        %omeg = 1.001;     % remove/add `%' at the beginning of this line 
        fun = @(x)(pi.*x./omeg)./sin(pi.*x./omeg);
        sgl1=zeros(1,100);
        sgl1(1:2:end) = omeg(1)*[1:1:50];
        sgl1(2:2:end) = -sgl1(1:2:end);
        if (omeg == 1.1),
            exact = 4.4677736463877657892361233985034;
            Nmin = [51,45,43,49,51,61,57,65,73,81,67,73,79,85,91,97,69,73,77,81,85];
        else
            exact = 12.929256850002296208119060780724;
            Nmin = [443,565,1000,489,532,445,967,1000,523,421,727,409,1000,...
         421,451,1000,715,469,419,1000,1000];  % omeg = 1.001
        end
    case 2, % Example 7.2
        X=NaN; T=NaN; E=NaN;
        disp('We did not use WG in Example 7.2');
    case 3, % Example 7.3
        omeg = 0.1;
        fun = @(x)(pi.*x./omeg)./sinh(pi.*x./omeg);
        sgl1=zeros(1,100);
        sgl1(1:2:end) = -sqrt(-1)*omeg*[1:1:50];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = 0.15707963267939592622470892489958;          
        Nmin = [213,199,233,211,205,211,209,271,221,287,265,235,253,271,...
            225,239,253,229,241,253,265];
    case 4, % Example 7.4
        omeg = 1.1;
        fun = @(x) sin(1./(x-omeg));
        sgl1 = omeg*ones(1,100);
        exact = -1.19245706732219214078770039225;       
        Nmin = [51,53,67,169,1000*ones(1,17)];
    case 5, % Example 7.5
        omeg = 0.1;
        fun = @(x) 1./(exp(pi.*x./omeg)+1);
        sgl1 = zeros(1,100);
        sgl1(1:2:end) = sqrt(-1)*omeg*[1:2:100];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = 1;
        Nmin = [213,199,233,211,205,211,209,235,221,243,265,235,253,211,...
            225,239,217,229,241,253,265];
    case 6, % Example 7.6
        fun = @(x) 1./sqrt((x+3).*(x+2));
        sgl1 = -2.5*ones(1,100);
        exact = 0.87116861981054736678487404465635;
        Nmin = [13,17,25,33,31,37,43,49,37,41,45,49,53,57,61,65,69,73,...
            77,81,85];
	%Nmin = [13,17,25,33,31,37,43,49,1000*ones(1,13)];
    case 7, % Example 7.7
        X=NaN; T=NaN; E=NaN;
        disp('We did not use WG in Example 7.7');
    otherwise % Example 7.8
        omeg = [2, 1.1, 1.01, 1.001];
        fun = @(x) [(pi.*x./omeg(1))./sin(pi.*x./omeg(1)), ...
            (pi.*x./omeg(2))./sin(pi.*x./omeg(2)); ...
            (pi.*x./omeg(3))./sin(pi.*x./omeg(3)), ...
            (pi.*x./omeg(4))./sin(pi.*x./omeg(4))];
        sgl1 = zeros(1,102);
        sgl1(1:6:end) = -omeg(4)*[1:1:17];
        sgl1(3:6:end) = -omeg(3)*[1:1:17];
        sgl1(5:6:end) = -omeg(2)*[1:1:17];
        sgl1(2:2:end) = -sgl1(1:2:end);
        exact = [2.33248723224655024110707565174 , ...
            4.4677736463877657892361233985034; ...
                 8.43018458047084205897126420479  , ...
                 12.929256850002296208119060780724];
        Nmin = [571,565,1000,489,1411,1000*ones(1,5),1167,409,...
            1000*ones(1,6),419,1000,1000]; 
end

N = [1:1:21];
n = length(N);
T = zeros(1,n);
 
if (ex < 7) && ~(ex == 2), % Examples 7.1, 7.3, 7.4, 7.5, and 7.6
    
    % n-point rational Gauss

    for M = 1:1:n,
        disp(sprintf('Execution time measurements countdown: %1.3f', n-M+1));
        sgl = sgl1(1:2*N(M)-1);
        m = N(M); Nm = Nmin(M);
        if (ex==3) || (ex==5), % (n+1)-point rational Gauss
            m = m+1;
            sgl= [sgl,sgl1(2*N(M))];
        end
        for K=1:10,
            tic;
            % warning('off', 'Octave:possible-matlab-short-circuit-operator');
            x = quadWG( sgl,m,Nm );
            % warning('on', 'Octave:possible-matlab-short-circuit-operator');
            A = fun(x(1,:))*x(2,:)';
            t1(K) = toc;
        end
        T(M) = mean(t1);
        X(M) = abs(exact-A)/abs(exact);
        test = errW([sgl(1:m-1),inf],x);
        E(M) = max(abs(test));
    end 
    index = find(X==0);
    X(index)=eps/10;
    
elseif (ex == 8), % Example 7.8
    
    % n-point rational Gauss
    
    for M = 1:1:n,
        disp(sprintf('Execution time measurements countdown: %1.3f', n-M+1));
        sgl = sgl1(1:2*N(M)-1);
        m = N(M); Nm = Nmin(M);
        for K=1:10,
            tic;
            % warning('off', 'Octave:possible-matlab-short-circuit-operator');
            x = quadWG( sgl,m,Nm);
            % warning('on', 'Octave:possible-matlab-short-circuit-operator');
            F = fun(x(1,:));
            m = length(x(1,:));
            for j = 1:2,
                A(:,j) = F(:,(j-1)*m+1:j*m)*x(2,:)';
            end
            t1(K) = toc;
        end
        T(M) = mean(t1);
        Xa = abs(exact-A)./abs(exact);
        X(M) = max(Xa(:));
        test = errW([sgl(1:m-1),inf],x);
        E(M) = max(abs(test));
    end 
    index = find(X==0);
    X(index)=eps/10;
    
end

if ~(ex==2) && ~(ex==7),
    if (ex==3) || (ex==5),
        N = N+1;
    end
    for k=1:length(T),
        Ts(k) = sum(T(1:k));
    end
    
    figure;
    semilogy(N,X);
    xlabel('n');
    ylabel('Relative error');
    title('Relative error against number of points');
    
    figure;
    loglog(Ts,X);
    xlabel('t');
    ylabel('Relative error');
    title('Relative error against execution time');
    
    figure;
    semilogy(N,Ts,'-',N,T,'--');
    xlabel('n');
    ylabel('t');
    hleg1 = legend('t_{1:1:n}','t_n');
    title('Execution time against number of points');

end

end

%--------------------------------------------------------------------------

function xw = quadWG( sgl , N, Nm )

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

Z = Zmatrix(-1./sgl);
M = size(Z,1); 
ab0 = r_jacobi(Nmax);
[abmod,Ncap,kount]=r_mod(N,ab0);
xw = gauss_rational(N,abmod);
xw=real(xw)';

end


%--------------------------------------------------------------------------

function [ Zm ] = Zmatrix( a )

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
