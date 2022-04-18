function X=colex(n,k,a,b)
% 
% Generate all integer compositions of n with k parts, each between 0 and b
% 
% n: integer whose compositions are to be generated
% k: number of parts
% a: minimum value of parts
% b: maximum value of parts
%
% The algorithm is directly taken from 
% Vincent Vajnovszki, Generating permutations with a given major index, http://arxiv.org/abs/1302.6558
%
% Matlab implementation:
% Steffen Eger, steffen.eger@yahoo.com, 11/5/2013
% 
%

X=[];
if k*a>n || k*b<n || n<0 || k<0 
  return
end
[minimum,im] = generateMin(n-k*a,k,b-a);
c = zeros(k,1);
% not impossible
if im==0
	% first, compute all compositions with parts p such that 0<=p<=b-a
	X=genColex(n-k*a,k,k,b-a,c,minimum,[]);
	% then add a such that a<=p<=b
	X = X+a;
end
end

% this implements Algorithm 2 in the paper
% b is the upper bound
function X=genColex(n,r,k,b,c,minimum,X)
if n==0
  % leave this if you just want to have it printed out (might be memory-saving)
  % disp(c')
  % leave this if want to store results in matrix X
  X(end+1,:)=c;
else
  if c(r)==b
    r = r-1;
  end
  l = minimum(n);
  for i=l:r
    if i==l
	e = n-(l-1)*b;
    else
	e = 1;
    end
    c(i) = c(i)+e;
    X=genColex(n-e,i,k,b,c,minimum,X);
    c(i) = c(i)-e;
  end
end
end

function [m,impossible]=generateMin(n,k,b)
m=zeros(n,1);
impossible=0;
for i=1:n
   q = find(i,k,b);
   if q==-1
	impossible=1;
	break
   end
   m(i) = q;
end
end

function t=find(n,k,b)
t=-1;
for s=1:k
  if s*b>=n
	t=s;
	break;
  end
end
end