function [objPoly,ineqPolySys,lbd,ubd] = optControl(M,Nx,Ny,mu);

% optControl(M,Nx,Ny,mu);
% 
% Thhe Discrite-Time optimal control problem, which is
% described in "An efficient trust region method for unconstrained
% discrete-time optimal control problems", Comp.Opt.Appl,4,p.47-66
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example is problem 1 
% in "An efficient trust region method for unconstrained
% disrete-time optimal control problems" T.F. Coleman et.al
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% <Input> 
% M:  the number of variable vectors
% Nx: the number of variables x
% Ny: the number of variables y
% mu: nonnegative scalar parameter
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

if mu < 0 
	error('## mu is nonnegative ##');
end

nDim = (Nx+Ny)*(M-1);
%%
%% Order of variable
%% y_{2,1},y_{2,2},...,y_{2,Ny},...,y_{M,Ny},x_{1,1},,,,
%% x_{1,Nx},...,x_{M-1,Nx}
%%
%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 4;
objPoly.noTerms  = 5*nDim+1;
objPoly.supports = zeros(objPoly.noTerms,objPoly.dimVar);

for i=1:(M-1)*Ny
  objPoly.supports(5*(i-1)+1,i) = 4;
  objPoly.supports(5*(i-1)+2,i) = 3;
  objPoly.supports(5*(i-1)+3,i) = 2;
  objPoly.supports(5*(i-1)+4,i) = 1;
  objPoly.coef(5*(i-1)+1,1)     = 1;
  objPoly.coef(5*(i-1)+2,1)     = 1;
  objPoly.coef(5*(i-1)+3,1)     = 3/8;
  objPoly.coef(5*(i-1)+4,1)     = 1/16;
  objPoly.coef(5*(i-1)+5,1)     = 1/256;
end
for i=(M-1)*Ny+1:nDim
  objPoly.supports(5*(i-1)+1,i) = 4;
  objPoly.supports(5*(i-1)+2,i) = 3;
  objPoly.supports(5*(i-1)+3,i) = 2;
  objPoly.supports(5*(i-1)+4,i) = 1;
  objPoly.coef(5*(i-1)+1,1)     = 1;
  objPoly.coef(5*(i-1)+2,1)     = 2;
  objPoly.coef(5*(i-1)+3,1)     = 3/2;
  objPoly.coef(5*(i-1)+4,1)     = 1/2;
  objPoly.coef(5*(i-1)+5,1)     = 1/16; 
end
objPoly.coef(5*nDim+1,1) = Ny/256;
objPoly = simplifyPolynomial(objPoly);

%Constraints--Inequalities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(M-1)*Ny
  ineqPolySys{i}.typeCone = -1;
  ineqPolySys{i}.sizeCone = 1;
  ineqPolySys{i}.dimVar   = nDim;
  if abs(mu) > 1.0e-5 
    if i < Ny+1
      ineqPolySys{i}.degree   = 1;
    else
      ineqPolySys{i}.degree   = 2;
    end
    ineqPolySys{i}.noTerms  = 1+Ny+Nx+Nx*Ny;
  else
    ineqPolySys{i}.degree = 1;
    ineqPolySys{i}.noTerms  = 1+Ny+Nx;
  end
end
for i=1:Ny
  ineqPolySys{i}.noTerms  = 1+Nx;
end

for i=1:(M-1)*Ny
  ineqPolySys{i}.supports = sparse(ineqPolySys{i}.noTerms,ineqPolySys{i}.dimVar);
end
%%make A,B,C
A = 0.5*speye(Ny);
A = A + diag(0.25*ones(Ny-1,1),1);
A = A + diag(-0.25*ones(Ny-1,1),-1);
for i=1:Ny
  for j=1:Nx
  B(i,j) = i-j; 
  C(i,j) = i+j;
  end
end
B = B/(Nx+Ny);
C = C*mu/(Nx+Ny);

%%y_{2,j}
for j=1:Ny
  ineqPolySys{j}.supports(1,j) = 1;
  for k=1:Nx
    ineqPolySys{j}.supports(k+1,(M-1)*Ny+k) = 1;
  end
  ineqPolySys{j}.coef =[1;-B(j,:)'];
  I = find(abs(ineqPolySys{j}.coef) > 1.0e-5);
  ineqPolySys{j}.coef = ineqPolySys{j}.coef(I,:);
  ineqPolySys{j}.supports = ineqPolySys{j}.supports(I,:);
  ineqPolySys{i}.noTerms = length(I);
  ineqPolySys{j}= simplifyPolynomial(ineqPolySys{j});
end


%%y_{i+1,j}

for i=1:M-2
  for j=1:Ny
    coef = [];
    ncount = 1;
    ineqPolySys{Ny*i+j}.supports(ncount,Ny*i+j) = 1;
    coef = [1];
    for k=1:Ny
      ncount = ncount+1;
      ineqPolySys{Ny*i+j}.supports(ncount,Ny*(i-1)+k) = 1;
    end
    coef = [coef;-A(j,:)'];
    for k=1:Nx
      ncount = ncount+1;
      ineqPolySys{Ny*i+j}.supports(ncount,Ny*(M-1)+Nx*i+k) = ...
	  1;  
    end
    coef = [coef;-B(j,:)'];
    if abs(mu) > 1.0e-5
      for k=1:Ny
	for ell=1:Nx
	  ncount = ncount+1;
	  ineqPolySys{Ny*i+j}.supports(ncount,Ny*(M-1)+Nx*i+ell) ...
	      = 1;
	  ineqPolySys{Ny*i+j}.supports(ncount,Ny*(i-1)+k) = ...
	      1;
	end
	coef = [coef;-C(k,:)'];  
      end
    end
    ineqPolySys{Ny*i+j}.coef = coef;
    I = find(abs(ineqPolySys{Ny*i+j}.coef) > 1.0e-5);
    ineqPolySys{Ny*i+j}.coef = ineqPolySys{Ny*i+j}.coef(I,:);
    ineqPolySys{Ny*i+j}.supports = ineqPolySys{Ny*i+j}.supports(I,:);
    ineqPolySys{Ny*i+j}.noTerms = length(I);
    ineqPolySys{Ny*i+j}= simplifyPolynomial(ineqPolySys{Ny*i+j});
  end
end
lbd = -1.0e+10*ones(1,nDim);
ubd  = 1.0e+10*ones(1,nDim);
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/optControl.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
