function [f,variables,zeros,options]=KSS(n)
% Generate the KSS system of n x n :
%     x_i^2+\sum_{j=1}^{n}x_j-2x_i-(n-1) = 0, 
%     for  i=1,2,...,n at the zero (1,...,1).
% The calling sequence is:
%     [f,variables,x,options]=KSS(n)
% where n is the size of system, f, variables, zeros and options can be used 
% by multiplicity directly to compute the multiplicity structure for the KSS 
% system. See Test_KSS for details.

variables=cell(1,n);
f=cell(1,n);
for j=1:n
    variables{j}=['x',num2str(j)];
    f{j}=['x',num2str(j),'^2-x',num2str(j),'-',num2str(n-1)];
    for k=1:n
        if k~=j
            f{j}=[f{j} '+x',num2str(k)];
        end
    end
end
zeros=ones(1,n);
options=optset('EqnType','Poly');
