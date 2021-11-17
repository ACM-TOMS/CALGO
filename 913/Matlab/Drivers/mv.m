function x = mv( y, A );
%
% Function for (preconditioned) matrix-vector multiplication.
% Symmetric Gauss-Seidel preconditioner + Eisenstat's trick.
% Called in idrs_ex2.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

t = A.U\y;
x = t + (A.L\( y - t ));
return
