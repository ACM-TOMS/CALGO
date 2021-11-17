% Example 2: Numerical approximation of an integral on the interval [-1,1]  
%            where the integrand has one multiple singularity

clear; 
close; 
clc; 

disp('One multiple pole A');
disp('int_{-1}^{1} ( 1/2 + (x-A)^(-100) ) dx');

fun = @(x,a) 1/2 + 1./(x-a).^100;
A = [1.1:.1:2]; % ten different cases

for j = 1:length(A),
    disp(sprintf('\nPole is %1.3f',A(j)));
    n = 32;
    disp(sprintf('The maximal number of iterations is %1.3f',n+1));
    fx = @(x) fun(x,A(j));
    sgl = A(j)*ones(1,n);
	[NumInt,Err] = rfejer(sgl,fx);
    disp(sprintf('Computed value: %1.16e',NumInt));
    disp(sprintf('Estimated relative error: %1.16e',Err));
    Exact = 1-1/99*(1/(1-A(j))^99+1/(1+A(j))^99);
    ErrExact = abs(1-NumInt/Exact);
    disp(sprintf('Exact value: %1.16e',Exact));
    disp(sprintf('Exact relative error: %1.16e',ErrExact));
end
