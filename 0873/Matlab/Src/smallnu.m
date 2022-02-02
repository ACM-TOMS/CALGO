%
%   File:              smallnu.m
%
%   Description:       Checks for small first component of eigenvector
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
%   Functions called:  none
%
%   Call: small = smallnu(absnu,noru,norg,tol)
%

function [small] = smallnu(absnu,noru,norg,tol)

         small = absnu * norg <= tol * noru;
         
