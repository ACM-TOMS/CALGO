%
%   File:              upd_param0.m
%
%   Description:       Computes alpha by one-point interpolation
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
%   Functions called:  interpol1, safe_alpha1
%
%   Call: alpha = upd_param0(ipoint,Delta,alpha,alphaL,alphaU,deltaU,message_level);
%

function [alpha] = upd_param0(ipoint,Delta,alpha,alphaL,alphaU,deltaU,message_level)

    alpha  = interpol1(ipoint,Delta,alpha);
    alpha  = safe_alpha1(alpha,alphaL,alphaU,deltaU,ipoint,message_level);

