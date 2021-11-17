function [objPoly,I01,Icomp] = QmatToSlackBoxQOP(Qmat0,lambda)
%QMATTOSLACKBOXQOP Add slack variables and convert it to BBCPOP
% Usage: [objPoly,I01,Icomp] = QmatToSlackBinQOP(Qmat0,lambda);
%
% Suppose that we have a binary QOP:
% min. x'Qx s.t. x\in[0,1]^n            (1)
%
% QMATTOSLACKBOXQOP first reformulate (1) by adding slack variable u 
% min. x'Qx s.t. x+u=e,  x,u\in[0,1]^n
%
% and then convert it to BBCPOP by applying Lagrangian relaxation
% min. x'Qx + lambda\|x+u-e\|^2  
% s.t. x,u\in[0,1]^n                    (2)
%
% It is reported that a DNN relaxation of (2) is stronger than that of (1).
% See Section 6.2 and 7.3 of
% 
% S. Kim, M. Kojima, and K.-C. Toh.
% Doubly Nonnegative Relaxations for Quadratic and Polynomial Optimization
% Problems with Binary and Box Constraints,
% Research Rport B-483, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, July 2016. 
% 
	if nargin==1
		lambda = 1.0e4; 	
	end
	nDim = size(Qmat0,1);
	normQmat0=norm(Qmat0,inf); 
	lambda = normQmat0*lambda;
	% Introducing slack variables x + u - e = 0;
	Qmat = [[Qmat0,sparse(nDim,nDim)];sparse(nDim,2*nDim)];
	objPoly = quad2Poly(Qmat,[],[]);
	% objPoly.sizeCone=1;
	% objPoly.typeCone=1; 
	for p=1:nDim
		poly.supports=sparse([2,3],[p,nDim+p],[1,1],3,2*nDim,3);
		poly.coef=[-1;1;1];
		poly= multiplyPoly(poly,poly);
		poly.coef=lambda*poly.coef;
		objPoly=addPoly(objPoly,poly); 
	end
	objPoly.sizeCone=1;
	objPoly.typeCone=1; 
	objPoly.dimVar=2*nDim;
	objPoly.noTerms=size(objPoly.coef,1);
    objPoly.UbdIX=nDim+1;
	
  	Icomp=[]; 
	I01=false(1,2*nDim);
end
