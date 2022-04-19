%
%   File:              eig_gateway.m
%
%   Description:       Computes two smallest eigenpairs
%                      of bordered matrix 
%                         | alpha  g' |
%                         |   g    Hg |
%                      using routine EIG from Matlab
%
%                      Assumes  H = A'*A  for
%                      Constrained Least Squares Problems
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
%   Functions called:  eig, diag, sort
%
%   Call: [nconv,lambda1,y1,lambda2,y2,v1] = eig_gateway(Balpha,...);
%

function [nconv,lambda1,y1,lambda2,y2,v1] = ...
          eig_gateway(Balpha,varargin)

   B       = [Balpha.alpha Balpha.g'; Balpha.g Balpha.H];
   [W,D]   = eig(full(B));
   D       = diag(D);               
   [D,ind] = sort(D);        
   W       = W(:,ind);

   nconv   = 2;
   lambda1 = D(1);
   y1      = W(:,1);
   lambda2 = D(2);
   y2      = W(:,2);

   v1      = y1;
