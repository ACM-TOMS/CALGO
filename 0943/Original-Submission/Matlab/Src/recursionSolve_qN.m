function [out] = recursionSolve_qN(a,b,b0,sigma,z)
%
%
% This function directly solves systems of the form (B+\sigma I)x=z, 
% where B is an L-BFGS quasi-Newton matrix.  
%
% Based on the paper:
% "Limited-Memory BFGS Systems with Diagonal Updates" 
% by Jennifer B. Erway and Roummel F. Marcia
% Linear Algebra and its Applications
% Volume 437, Issue 1 (2012) pp 333-344.
%
% Copyright (2011): Jennifer B. Erway and Roummel F. Marcia
%
% A link to the paper and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% 
% This code is distributed under the terms of the GNU General Public License
% 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------


% Inputs: a(:,k) and b(:,k) make up the kth L-bfgs update and
% B_0 = b0*I is the initial L-BFGS matrix, i.e., 
% B_k = B_0 - sum_{i=1}^k a(:,i)a(:,i)^T + sum_{i=1}^{k} b(:,i)b(:,i)^T

% initialization
k0 = size(a,2);
k  = 2*k0;
inv_c0 = 1/( b0 + sigma );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursion formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invCz  = inv_c0*z;

for j = 1:k;
    if (mod(j,2)==1)
        c  = a(:,(j+1)/2);
    else
        c  = b(:,j/2);
    end
    p(:,j) = inv_c0*c;
    for t = 1:j-1
        p(:,j) = p(:,j) + (-1)^(t+1)*v(t)*(p(:,t)'*c)*p(:,t);
    end
    v(j) = 1/(1 + ((-1)^j)*c'*p(:,j));
    invCz = invCz + (-1)^(j+1)*v(j)*(p(:,j)'*z)*p(:,j);

end

out = invCz;



