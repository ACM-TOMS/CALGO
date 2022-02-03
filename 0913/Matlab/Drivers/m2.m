function x = m2( y, M2 );
%
% Function for (right-hand side factor) of preconditioner.
% Called by idrs_ex4.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

x = M2.U\y;
return
