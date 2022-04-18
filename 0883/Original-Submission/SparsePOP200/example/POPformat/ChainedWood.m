function [objPoly,ineqPolySys,lbd,ubd] = ChainedWood(nDim)
%
% ChainedWood(nDim);
% 
% The chained wood function, which is
% described in "Testing a class of methods for solving minimization
% problems with simple bounds on the variables", Math.Comp., 50,
% p.399-430
%
% <Input> 
% nDim: the dimension of the function
% Remark: nDim needs to be a multiple of 4. 
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

if mod(nDim,4) ~=0
	error('## nDim should be a multiple of 4 ##');
end

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = ceil(nDim/2-1);
poly = cell(N,1);
for j=1:N
  poly{j}.typeCone = 1;
  poly{j}.sizeCone = 1;
  poly{j}.dimVar   = nDim;
  poly{j}.degree   = 4;
  poly{j}.noTerms  = 21;
  poly{j}.supports = sparse(poly{j}.noTerms,poly{j}.dimVar);
  i = 2*j-1;
  poly{j}.supports(1,i+1)   = 2;
  poly{j}.supports(2,i)     = 2;
  poly{j}.supports(2,i+1)   = 1;
  poly{j}.supports(3,i)     = 4;
  poly{j}.supports(5,i)     = 1;
  poly{j}.supports(6,i)     = 2;
  poly{j}.supports(7,i+3)   = 2;
  poly{j}.supports(8,i+3)   = 1;
  poly{j}.supports(8,i+2)   = 2;
  poly{j}.supports(9,i+2)   = 4;
  poly{j}.supports(11,i+2)  = 1;
  poly{j}.supports(12,i+2)  = 2;
  poly{j}.supports(13,i+1)  = 2;
  poly{j}.supports(14,i+3)  = 2;
  poly{j}.supports(16,i+1)  = 1;
  poly{j}.supports(16,i+3)  = 1;
  poly{j}.supports(17,i+1)  = 1;
  poly{j}.supports(18,i+3)  = 1;
  poly{j}.supports(19,i+1)  = 2;
  poly{j}.supports(20,i+1)  = 1;
  poly{j}.supports(20,i+3)  = 1;
  poly{j}.supports(21,i+3)  = 2;
  
  poly{j}.coef = [100;-200;100;1;-2;1;90;-180;90;1;-2;1;10;10;40;20;-40;-40;0.1;-0.2;0.1];
  
end

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 4;
objPoly.noTerms  = 1;
objPoly.supports = sparse(1,nDim);
objPoly.coef     = 1;
for i=1:N
  objPoly.noTerms = objPoly.noTerms + poly{i}.noTerms;
  objPoly.supports= [objPoly.supports;poly{i}.supports];
  objPoly.coef    = [objPoly.coef;poly{i}.coef];
end
objPoly = simplifyPolynomial(objPoly);
ineqPolySys = [];
lbd = -1.0e+10*ones(1,nDim);
ubd  = 1.0e+10*ones(1,nDim);

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/ChainedWood.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
