function [objPoly,ineqPolySys,low,up] = optControl2(N);

% optControl2(N);
% 
% The Discrite-Time optimal control problem, which is
% described in "An efficient trust region method for unconstrained
% discrete-time optimal control problems", Comp.Opt.Appl,4,p.47-66
%
% <Input> 
% N: 2N is the number of variable x and y.
%
% <Output>
% [objPoly,ineqPolySys,low,up]

N_x = N-1;
N_y = N-1;
nDim = 2*N-2;

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 2;
objPoly.noTerms  = 2*N-2;
objPoly.supports = sparse(objPoly.noTerms,objPoly.dimVar);

objPoly.supports = 2*speye(objPoly.noTerms);
objPoly.coef = repmat(1/N,objPoly.noTerms,1);
objPoly = simplifyPolynomial(objPoly);

%Constraints--Inequalities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-1
  ineqPolySys{i}.typeCone = -1;
  ineqPolySys{i}.sizeCone = 1;
  ineqPolySys{i}.dimVar   = nDim;
  ineqPolySys{i}.degree   = 2;
  ineqPolySys{i}.noTerms  = 4;
end
ineqPolySys{1}.noTerms  = 3;

for i=1:N-1
  ineqPolySys{i}.supports = sparse(ineqPolySys{i}.noTerms,ineqPolySys{i}.dimVar);
end
ineqPolySys{1}.degree = 1;
ineqPolySys{1}.supports(1,1) =1;
ineqPolySys{1}.supports(2,N_x+1) =1;
ineqPolySys{1}.coef =[1/N;1;-(1+1/N)];
ineqPolySys{1} = simplifyPolynomial(ineqPolySys{1});

for i=2:N-1
  ineqPolySys{i}.supports(1,i)=1;
  ineqPolySys{i}.supports(2,N_x+i)=1;
  ineqPolySys{i}.supports(3,N_x+i-1)=1;
  ineqPolySys{i}.supports(4,N_x+i-1)=2;
  ineqPolySys{i}.coef = [1/N;1;-1;-1/N];
  ineqPolySys{i} = simplifyPolynomial(ineqPolySys{i});
end

low = -1.0e+10*ones(1,nDim);
up  = 1.0e+10*ones(1,nDim);
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/optControl2.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
