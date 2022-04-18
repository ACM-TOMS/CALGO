%
%   File:         atamv.m
%
%   Description:  A'*A times a vector
%
%   Authors:      Marielba Rojas
%                    mr@imm.dtu.dk
%                 Sandra Santos
%                 Danny Sorensen
%
%   Version:      1.2
%
%   System:       MATLAB 6.0 or higher
%
%   Date:         15 March 2007
%
%   Call: w = atamv(v,atamvpar);
%

function [w] = atamv(v,atamvpar)

w = atamvpar.A*v;
w = (w'*atamvpar.A)';
