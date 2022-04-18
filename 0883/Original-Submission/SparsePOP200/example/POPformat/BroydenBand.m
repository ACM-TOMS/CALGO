function [objPoly,ineqPolySys,lbd,ubd] = BroydenBand(nDim);

% BroydenBand(nDim);
% 
% The Broyden Banded function, which is
% described in "Testing Unconstrained Optimization Software",
% J.J.More et.al, ACM Trans. Math. Soft., 7, p.17-41
%
% <Input> 
% nDim: the dimension of the function
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

%%default
u = 5;
ell = 1;
%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nDim
  poly{i}.typeCone = 1;
  poly{i}.sizeCone = 1;
  poly{i}.dimVar   = nDim;
  poly{i}.degree   = 3;
  
  I = zeros(1,nDim);
  I(max(1,i-u):min(nDim,i+ell)) = 1;
  idx2 = find(I); 
  I(i) = 0;
  idx = find(I);
  poly{i}.noTerms  = 2*length(idx)+2;
  
  poly{i}.supports =zeros(poly{i}.noTerms,poly{i}.dimVar);
  poly{i}.supports(1,i) = 1;
  poly{i}.supports(2,i) = 3;
  if ~isempty(idx) 
    for j=1:length(idx)
      poly{i}.supports(2*j+1,idx(j)) = 1;
      poly{i}.supports(2*j+2,idx(j)) = 2;
    end
    poly{i}.coef     = [2;5;-ones(2*length(idx),1)];
  else
    poly{i}.coef     = [2;5];
  end
end

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 6;
objPoly.noTerms  = 0;
objPoly.supports = [];

objPoly.coef     = [];
for i=1:nDim
  temp = multiplyPolynomials(poly{i}, poly{i});
  objPoly.noTerms = objPoly.noTerms + temp.noTerms + poly{i}.noTerms+1;
  objPoly.supports= [objPoly.supports;temp.supports;poly{i}.supports;zeros(1,nDim)];
  objPoly.coef    = [objPoly.coef;temp.coef;2*poly{i}.coef;1];
end
objPoly = simplifyPolynomial(objPoly);
ineqPolySys = [];
lbd = -1.0e+10*ones(1,nDim);
ubd  = 1.0e+10*ones(1,nDim);

%writePolynomials(1,objPoly);

return;

function polyOut = multiplyPolynomials(poly1,poly2, order)
%
% Multiplication of two real valued polynomials 
%       polynomial.nDim --- a positive integer 
%       polynomial.degree --- a positive integer 
%       polynomial.noTerms --- a positive integer 
%       polynomial.supports --- a noTerms \times nDim matrix
%       polynomial.coef --- a noTerms \times 1 vector

%Now, this function can deal with only typeCone = 1 case.
%2005 Jan 9

if nargin < 3
  order = 'grevlex';
end

if poly1.dimVar ~= poly2.dimVar
  error('!!! poly1.nDim ~= poly2.nDim !!!');
end

polyIn.typeCone = poly1.typeCone;
polyIn.sizeCone = poly1.sizeCone;
polyIn.dimVar = poly1.dimVar;
polyIn.noTerms = poly1.noTerms * poly2.noTerms; 
polyIn.supports = sparse(polyIn.noTerms,polyIn.dimVar);
UsedVariable = find(any([poly1.supports;poly2.supports],1));

idx1 = repmat((1:poly1.noTerms),poly2.noTerms,1);
idx1 = idx1(:);
FirSup = poly1.supports(idx1,UsedVariable);
idx2 = (1:poly2.noTerms)';
idx2 = repmat(idx2,poly1.noTerms,1);
SecSup = poly2.supports(idx2,UsedVariable);
polyIn.supports(:,UsedVariable) = FirSup + SecSup;

tempCoef = poly2.coef * poly1.coef';
polyIn.coef = tempCoef(:);
    
degree = sum(polyIn.supports(:,UsedVariable),2);
polyIn.degree = max(degree);
polyOut = simplifyPolynomial(polyIn,'grevlex');
clear degree polyIn UsedVariable tempCoef idx1 idx2 FirSup SecSup poly1 poly2
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/BroydenBand.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
