%
%   File:              upd_paramk.m
%
%   Description:       Updates and safeguards alpha
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
%   Functions called:  interpol2, safe_alphak
%
%   Call: alpha = upd_paramk(ipointa,ipointb,Delta,alphaL,alphaU,deltaU,message_level);
%

function [alpha] = upd_paramk(ipointa,ipointb,Delta,alphaL,alphaU,deltaU,message_level)

   alpha  = interpol2(ipointa,ipointb,deltaU,Delta,alphaL,alphaU,message_level);
   alpha  = safe_alphak(alpha,alphaL,alphaU,deltaU,ipointa,ipointb,message_level);

