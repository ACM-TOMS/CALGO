function [f,variables,x]=ten_fold(n)
% Generate the ten fold system of n x n :
%     x_1 +...+ x_{i}  =  0, for i = 1, 2,..., n-2,
%     (x_1 + ...  + x_{n-1})^2  =  0, 
%     (x_1 + \cdots + x_{n})^5  =  0 at (0,...,0).
% The calling sequence is:
%     [f,variables,x,options]=ten_fold(n)
% where n is the size of system, f, variables, zeros and options can be used 
% by multiplicity directly to compute the multiplicity structure for the ten 
% fold system. See Test_ten_fold for details.

variables=cell(1,n);
f=cell(1,n);
for j=1:n
    variables{j}=['x',num2str(j)];
end
for j=1:n-2
    if j==1
        f{j}='x1';
    else
        f{j}=[f{j-1},'+x',num2str(j)];
    end
end
f{n-1}='(';
f{n}='(';
for j=1:n-1
    if j==1
        f{n-1}=[f{n-1},variables{j}];
        f{n}=[f{n},variables{j}];
    else
        f{n-1}=[f{n-1},'+',variables{j}];
        f{n}=[f{n},'+',variables{j}];
    end
end
f{n-1}=[f{n-1},')^2'];
f{n}=[f{n},'+x',num2str(n),')^5'];
x=zeros(1,n)+zeros(1,n)*1i;
