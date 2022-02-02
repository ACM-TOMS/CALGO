function [c, d, hvt, sigma] = findcd(sigma, hvt, d)
 
n = size(sigma,1);
 
%permute rows to get HVT on the diagonal
sigma = sigma(hvt,:); 
 
if nargin==2
   d = max(sigma);
end;
 
%diagonal of sigma
diag_sigma = diag(sigma);
%initial c
c = d' - diag_sigma;
d_old = d;
 
%sigma may contain -oo
%add 1 to sigma so we do not have 0's
S = sigma+1;
% replace -inf's by 0's
i = find(S==-inf); 
S(i) = 0;
% i,j are indices of elements > -inf
% s is the elements S(i,j)+1
[i,j,s] = find(S);
 
while 1
   %add c to S
   C = sparse(i,j,s+c(i));
   %take the max of each column and subtract 1 to account for the +1
   d = max(C)-1; 
   
   if d == d_old
       break;
   end
   c = d' - diag_sigma;
   d_old = d;
end
 
c(hvt) = c; % permute c back
end