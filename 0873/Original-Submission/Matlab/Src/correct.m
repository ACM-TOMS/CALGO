%
%   File:              correct.m
%
%   Description:       Computes correction vector in the Hard Case
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
%   Functions called:  norm, sqrt
%
%   Call: x = correct(z,x,norx,Delta);
%


function [x] = correct(z,x,norx,Delta)

    z     = z/norm(z);
    tau   =  - z'*x;
    tau1  = (Delta - norx)*(Delta + norx);
    tau   = - tau1/(tau + sign(tau)*sqrt(tau^2+tau1));

    x     = x + tau*z;
