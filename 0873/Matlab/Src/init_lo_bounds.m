%
%   File:              init_lo_bounds.m
%
%   Description:       Computes initial alphaL (lower bound for alpha) and 
%                      deltaL (lower bound for delta1)
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
%   Functions called:  disp, sprintf
%
%   Call: alphaL = init_lo_bounds(epair,norg,Delta,message_level);
%

function [alphaL] = init_lo_bounds(epair,norg,Delta,message_level)

      deltaL = epair.lambda;

      alphaL = deltaL - norg/Delta;

     if (message_level == 2)
        disp('--------------')
        disp('init_lo_bounds')
        disp('--------------')
        disp(sprintf('\nalphaL = %6.4e',alphaL))
        disp(sprintf('deltaL = %6.4e\n',deltaL))
     end
