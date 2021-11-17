function x = m1( y, M1 );
%
% Function for (left-hand side factor) of preconditioner.
% Called by idrs_ex4.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

x = M1.L\ y;
return
