function [objPoly,ineqPolySys,lbd,ubd] = Rosenbrock(nDim,s);
%
% Rosenbrock function, which is
% described in "Newton-Type Minimization via the Lanczos method",
% SIAM J.Numer.Anal., 21, p.770-788.
%
%   1 + \sum_{i=1}^nDim ( 100(x_i-x_{i-1}^2)^2 + (1-x_i)^2 ). 
%
% <Input> 
% nDim: the dimension of the function.
% s 
%   = -1 for the constraint x_1 <= 0,
%   = 0  for no constraint,
%   = 1  for the constraint x_1 >= 0.
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 4;
objPoly.noTerms  = 5*(nDim-1)+1;
objPoly.supports = sparse(objPoly.noTerms,objPoly.dimVar);
objPoly.coef     = sparse(objPoly.noTerms,1);

for i=2:nDim
  objPoly.supports(5*(i-2)+1,i)   = 2;
  objPoly.supports(5*(i-2)+2,i)   = 1;
  objPoly.supports(5*(i-2)+2,i-1) = 2;
  objPoly.supports(5*(i-2)+3,i-1) = 4;
  objPoly.supports(5*(i-2)+4,i)   = 1;
  objPoly.coef(5*(i-2)+1)   = 101;
  objPoly.coef(5*(i-2)+2)   = -200;
  objPoly.coef(5*(i-2)+3)   = 100;
  objPoly.coef(5*(i-2)+4)   = -2;
  objPoly.coef(5*(i-2)+5)   = 1;
end
objPoly.coef(5*(nDim-1)+1) = 1;

%%%%%%%%%%%%%
% Important %
%%%%%%%%%%%%%
% objPoly constructed above has multiple terms with a common support, 
% so we need to simplify objPoly by combining two terms with a common 
% support.
objPoly = simplifyPolynomial(objPoly);
%

ineqPolySys = [];

lbd = -1.0e+10*ones(1,nDim);
ubd = 1.0e+10*ones(1,nDim);
if nargin == 2
    if s > 0
        lbd(1,1) = 0;
    elseif s < 0 
        ubd(1,1) = 0;
    end
end
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/Rosenbrock.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
