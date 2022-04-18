%
%   File:              uutmatvec.m
%
%   Description:       (I-2uu') D (I-2uu') times vector
%                      D is a diagonal matrix, u is a unit vector
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
%   Call: w = uutmatvec(v,uutmatvecpar);
%

function [w] = uutmatvec(v,uutmatvecpar)

     d = uutmatvecpar.d;
     u = uutmatvecpar.u;

     w = v - 2 * (u'*v) * u;
     w = d .* w;
     w = w - 2 * (u'*w) * u;

