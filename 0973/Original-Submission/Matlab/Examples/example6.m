% Example 6: Semi-automatic rational Fejer and adaptive methods

clear; 
close; 
clc; 

disp('Semi-automatic and Adaptive')
disp('One multiple pole A');
disp('int_{-1}^{1} sin( 1/(x-A) ) dx');

w = 1.1;
fun = @(x) sin(1./(x-w)); % on the interval [-1,1]
exact = -1.19245706732219214078770039225;

disp(sprintf('\nPole is %1.3f',w));

M = 10;
N = 15;
NP = 100;

for K = 1:M,
    
    disp(sprintf('Execution time measurements countdown: %1.3f', M-K+1));

% fixed number of points

    tic;
    [y,e]=rfejer(w*ones(1,N-1)); 
    x = fun(y(1,:))*y(2,:)';
    t(K) = toc;

% semi-automatic

    tic;
    [xsa,esa,csa,ysa]=rfejer(w*ones(1,NP),fun); 
    tsa(K) = toc;
    
% optimal fixed number of points

    Nopt = length(ysa(1,:));
    tic;
    [yopt,eopt]=rfejer(w*ones(1,Nopt-1)); 
    xopt = fun(yopt(1,:))*yopt(2,:)';
    topt(K) = toc;
    
% semi-automatic with limited number of iterations

    Ninit = Nopt-2;
    tic;
    [xsaI,esaI,csaI,ysaI]=rfejer(w*ones(1,NP),fun,'Nmin',Ninit); 
    tsaI(K) = toc;
    
% adaptive 1

    tic;
    ival = [-1,.5];
    [fout,sglout]=transf(fun,w,ival);
    [xap1(1),eap1(1),cap1(1)]=rfejer(sglout*ones(1,N-1),fout,'Nmin',N);
    ival = [.5,.875];
    [fout,sglout]=transf(fun,w,ival);
    [xap1(2),eap1(2),cap1(2)]=rfejer(sglout*ones(1,N-1),fout,'Nmin',N);
    ival = [.875,.96875];
    [fout,sglout]=transf(fun,w,ival);
    [xap1(3),eap1(3),cap1(3)]=rfejer(sglout*ones(1,N-1),fout,'Nmin',N);
    ival = [.96875,1];
    [fout,sglout]=transf(fun,w,ival);
    [xap1(4),eap1(4),cap1(4)]=rfejer(sglout*ones(1,N-1),fout,'Nmin',N);
    xa1 = sum(xap1);
    ta1(K) = toc;
    
% adaptive 2

    tic;
    k = 7;
    p = nextpoint(1,w,k*w);
    ival = [p,1]; NN = 1;
    while (p > -1) && (N < 100),
        NN = NN+1;
        p2 = nextpoint(p,w,k*w);
        if (p2 > -1),
            ival = [ival;[p2,p]];
        else
            ival = [ival;[-1,p]];
        end
        p = p2;
    end
    [fout,sglout]=transf(fun,w,ival);
    w2 = k*w;
    [xap,eap,cap]=rfejer(w2*ones(1,N-1),fout,'Nmin',N);
    xa = sum(xap(:));
    ta(K) = toc;
    
% quadgk

    tic;
    [Q,E] = quadgk(fun,-1,1,'RelTol',100*eps,'AbsTol',0);
    tQ(K)=toc;
    
end

% optimal fixed number of points
topt = mean(topt);
erropt = abs(xopt-exact)/abs(exact);
Lopt=length(yopt(1,:));

disp(sprintf('\n N-point rational Fejer with N = %1.3f', Lopt));
disp(sprintf('Exact relative error: %1.16e', erropt));
disp(sprintf('Execution time: %1.16e',topt));

% semi-automatic
tsa = mean(tsa);
errsa = abs(xsa-exact)/abs(exact);
Lsa=length(ysa(1,:));

disp(sprintf('\n Semi-automatic rational Fejer:'));
disp(sprintf('Number of poles: %1.3f', NP));
disp(sprintf('Initial iteration: %1.3f', 1));
disp(sprintf('Number of iterations: %1.3f', Lsa));
disp(sprintf('Estimated relative error: %1.16e', esa));
disp(sprintf('Exact relative error: %1.16e', errsa));
disp(sprintf('Execution time: %1.16e',tsa));

% semi-automatic with limited number of iterations
tsaI = mean(tsaI);
errsaI = abs(xsaI-exact)/abs(exact);
LsaI=length(ysa(1,:));

disp(sprintf('\n Semi-automatic rational Fejer with limited number of iterations:'));
disp(sprintf('Number of poles: %1.3f', NP));
disp(sprintf('Initial iteration: %1.3f', Ninit));
disp(sprintf('Number of iterations: %1.3f', LsaI-Ninit+1));
disp(sprintf('Estimated relative error: %1.16e', esaI));
disp(sprintf('Exact relative error: %1.16e', errsaI));
disp(sprintf('Execution time: %1.16e',tsaI));

% fixed number of points
t = mean(t);
err = abs(x-exact)/abs(exact);
L=length(y(1,:));

disp(sprintf('\n N-point rational Fejer with N = %1.3f', N));
disp(sprintf('Exact relative error: %1.16e', err));
disp(sprintf('Execution time: %1.16e',t));

% adaptive 1
ta1 = mean(ta1);
ea1 = sum(abs(xap1).*eap1)/abs(xa1);
erra1 = abs(xa1-exact)/abs(exact);

disp(sprintf('\n Adaptive rational Fejer (1):'));
disp(sprintf('Number of subintervals: %1.3f',4));
disp(sprintf('Number of different quadrature formulae: %1.3f',4));
disp(sprintf('Number of points in each subinterval: %1.3f',N));
disp(sprintf('Estimated relative error: %1.16e', ea1));
disp(sprintf('Exact relative error: %1.16e', erra1));
disp(sprintf('Execution time: %1.16e',ta1));

% adaptive 2
ta = mean(ta);
ea = sum(abs(xap).*eap)/abs(xa);
erra = abs(xa-exact)/abs(exact);
La = length(ival(:,1));

disp(sprintf('\n Adaptive rational Fejer (2):'));
disp(sprintf('Number of subintervals: %1.3f', La));
disp(sprintf('Number of different quadrature formulae: %1.3f',1));
disp(sprintf('Number of points in each subinterval: %1.3f',N));
disp(sprintf('Estimated relative error: %1.16e', ea));
disp(sprintf('Exact relative error: %1.16e', erra));
disp(sprintf('Execution time: %1.16e',ta));

% quadgk
tQ = mean(tQ);
ERRQ = abs(Q-exact)/abs(exact);
EQ = abs(E)/abs(Q);

disp(sprintf('\n Matlab''s built-in QUADGK:'));
disp(sprintf('Minimal number of subintervals: %1.3f',10));
disp(sprintf('Number of points in each subinterval: %1.3f',15));
disp(sprintf('Estimated relative error: %1.16e', EQ));
disp(sprintf('Exact relative error: %1.16e', ERRQ));
disp(sprintf('Execution time: %1.16e',tQ)); 
