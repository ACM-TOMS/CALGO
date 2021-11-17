%
%   File:              quadratic.m
%
%   Description:       computes the value of the quadratic
%                      psi(x) = 1/2 x'Hx + g'x
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
%   Functions called:  matvec
%
%   Call: psi = quadratic(Balpha,x);
%

function [psi] = quadratic(Balpha,x)

      Balpha.bord = 0;
      psi         = 0.5*x'*matvec(x,Balpha) + Balpha.g'*x;
