%
%   File:              tcheigs_lstrs_gateway.m
%
%   Description:       Computes two smallest eigenpairs
%                      of bordered matrix using ARPACK
%                      Implicitly Restarted Lanczos Method 
%                      implemented in the Matlab routine EIGS 
%
%                      This version implements a Tchebyshev Polynomial
%                      Filter (T): the eigenvalues of T(B) (and not of B)
%                      are computed. The eigenvalues of B are recovered
%                      by means of the Rayleigh quotient.
%
%   Authors:           Marielba Rojas
%                        mr@imm.dtu.dk
%                      Sandra Santos
%                      Danny Sorensen
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              15 March 2007
%
%   Functions called:  ceil, diag, length, eigs_lstrs, threshold
%
%   Call: [nconv,lambda1,y1,lambda2,y2,v1] =    ...
%         tcheigs_lstrs_gateway(Balpha,eopts);
%

function [nconv,lambda1,y1,lambda2,y2,v1] =     ...
          tcheigs_lstrs_gateway(Balpha,eopts)

   k     = eopts.k;
   eopts = rmfield(eopts,'k');

   firsteopts   = eopts;
   firsteopts.p = 4;
 
   m = Balpha.dim;
    
   [W,D,flag,nconv,smallritz] = ...
    eigs_lstrs(@matvec,m,1,'LA',firsteopts,Balpha);

   D = diag(D);               
   b = D(1);
   a = threshold(smallritz,b,m);

   eopts.maxit = ceil(eopts.maxit/3);
   maxdegree   = 10;

   [W,D,flag,nconv,smallritz,v1] = ...
    eigs_lstrs(@tchmatvec,m,k,'SA',eopts,Balpha,a,b,maxdegree);

   D = diag(D);

%  Computes eigenvalues via Rayleigh quotient
   for i = 1:nconv
       D(i) = W(:,i)'*matvec(W(:,i),Balpha);
   end

   lambda1 = D(1);
   y1      = W(:,1);
   lambda2 = D(2);
   y2      = W(:,2);


%
%   Call: lb = threshold(a,b,n);
%
%   Functions called:  acosh, cosh
%

function [lb] = threshold(a,b,n)

     delta = 10.0;
     ub = b;        % largest eigenvalue
     par = cosh(acosh(delta)/n);
     lb = (2*a+(par-1)*ub) / (par+1);

