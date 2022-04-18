function [f,variables,x,options]=cyclic_cubic(n)
% Generate the cyclic cubic system of n x n :
%    x^3_i - x_{i+1} x_{i+2} = 0, for i = 1, 2,..., n, 
%    x_{n+1} = x_1 and x_{n+2} = x_2 at zero (0,...,0).
% The calling sequence is:
%     [f,variables,x,options]=cyclic_cubic(n)
% where n is the size of system, f, variables, zeros and options can be used 
% by multiplicity directly to compute the multiplicity structure for the 
% cyclic cubic system. See Test_cyclic_cubic for details.

variables=cell(1,n+2);
for j=1:n+2
    variables{j}=['x',num2str(j)];
end
f=cell(1,n+2);
f{1}=['x',num2str(n+1),'-x1'];
f{2}=['x',num2str(n+2),'-x2'];
for j=1:n
    f{j+2}=['x',num2str(j),'^3-x',num2str(j+1),'*x',num2str(j+2)];
end
x=zeros(1,n+2)+zeros(1,n+2)*1i;
options=optset('EqnType','Poly');
