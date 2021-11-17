%
%   File:              matvec.m
%
%   Description:       Matrix-Vector Product with
%                      Bordered Matrix 
%           
%                          | alpha  g' |
%                          |   g    H  |
%           
%                      It also computes Hx if alpha = 0 and g = 0 
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
%   Functions called:  any, feval, isa, length
%
%   Call: w = matvec(v,Balpha);
%

function [w] = matvec(v,Balpha)

%
% Declare and update global variable that counts matrix-vector products
% Routines that use it: lstrs_method, matvec, output
%
global mvp_lstrs;
mvp_lstrs = mvp_lstrs + 1;

     if (Balpha.bord)
        n  = length(v);
        nu = v(1);
        v  = v(2:n);
     end

     if (isa(Balpha.H,'double'))
        w = Balpha.H*v;
     else     
        w = feval(Balpha.H,v,Balpha.Hpar);
     end

     if (Balpha.bord)
        w = [Balpha.alpha * nu + Balpha.g'*v; nu*Balpha.g + w];
     end

