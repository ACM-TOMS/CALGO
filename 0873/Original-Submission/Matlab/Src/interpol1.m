%
%   File:              interpol1.m
%
%   Description:       1-point interpolation scheme
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
%   Call: alpha = interpol1(ipoint,Delta,alpha);
%

function [alpha] = interpol1(ipoint,Delta,alpha)

      lambda = ipoint.lambda;
      norx   = ipoint.norx;

      alpha  = alpha + ((alpha - lambda)/norx) * ...
                       ((Delta - norx)/Delta) * (Delta + 1/norx);
