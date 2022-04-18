%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% This is a demo driver for the More-Sorensen Sequential
%% (MSS) method.
%% 
%% For details, see
%%
%% "MSS: MATLAB Software for L-BFGS Trust-Region 
%%       Subproblems for Large-Scale Optimization"
%%
%% Authors: Jennifer B. Erway and Roummel F. Marcia
%%
%% Paper available at 
%%
%% http://www.wfu.edu/~erwayjb/publications.html
%%
%% To run: type >> demo
%% Change problem size by changing n in demo.m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

global S Y RHO;

%set parameters
M     =    1;   %number of L-BFGS updates (must be greater 0)
n     = 1000;   %dimension of problem

%form the matrices S and Y from quasi-Newton method
fprintf('\nGenerating random problem...');
S   = randn(n,M);
Y   = randn(n,M);
RHO = [];

%ensure B is positive definite
for i=1:M
  if S(:,i)'*Y(:,i)<0
     S(:,i)=-S(:,i);
  end
  RHO = [RHO , 1/(Y(:,i)'*S(:,i))];
end

%set gamma
if M>0
  gamma_qN = (S(:,M)'*Y(:,M)) / (Y(:,M)'*Y(:,M));
else
  gamma_qN = 1;
end

%generate right-hand side and delta
g     = randn(n,1);
delta = rand(1,1);
fprintf(' done.\n');

%set tolerance
tol    = norm(g)*1e-5;  %1e-2
tolmin = min( tol, sqrt(eps) );  

fprintf('Parameters:\n');
fprintf('  n       =  %8d\n', n);
fprintf('  delta   =  %8.2e \n', delta); 
fprintf('  gamma   =  %8.2e \n', gamma_qN);
fprintf('  ||g||   =  %8.2e \n', norm(g));
fprintf('  tau     =  %8.2e \n', tolmin);

%ms
if n<10000
  fprintf('\nRunning More-Sorensen method...\n');
  tic
  [s1,sigma1] = ms(g,delta,gamma_qN,tolmin);
  fprintf('Total time: %8.2e\n', toc);
end

%mss.m
fprintf('\nRunning More-Sorensen Sequential (MSS) method...\n');
tic
[s2,sigma2] = mss(g,delta,gamma_qN,tolmin);
fprintf('Total time: %8.2e\n', toc);


fprintf('\n');



