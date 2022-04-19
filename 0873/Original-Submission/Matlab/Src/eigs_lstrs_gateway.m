%
%   File:              eigs_lstrs_gateway.m
%
%   Description:       Computes two smallest eigenpairs
%                      of bordered matrix using ARPACK
%                      Implicitly Restarted Lanczos Method 
%                      implemented in the Matlab routine EIGS 
%
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
%   Functions called:  diag, length, eigs_lstrs, rmfield
%
%   Call: [nconv,lambda1,y1,lambda2,y2,v1] =    ...
%         eigs_lstrs_gateway(Balpha,eopts);
%

function [nconv,lambda1,y1,lambda2,y2,v1] =     ...
          eigs_lstrs_gateway(Balpha,eopts)

   k     = eopts.k;
   eopts = rmfield(eopts,'k');
 
   m     = Balpha.dim;

   [W,D,flag,nconv,smallritz,v1] = ...
    eigs_lstrs(@matvec,m,k,'sa',eopts,Balpha);

   D = diag(D);               

   lambda1 = D(1);
   y1      = W(:,1);
   lambda2 = D(2);
   y2      = W(:,2);
