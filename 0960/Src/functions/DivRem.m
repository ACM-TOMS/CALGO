function [q,r] = DivRem(a,b)
%DIVREM [q,r]=DivRem(a,b) is an error free transformation for the division such
% that q=fl(a/b) and a=b*q+r
%DIVREM Error free transformation of a/b into a=b*q+r with q=fl(a/b)
%
%   [q,r] = DivRem(a,b)
%
%On return, a=b*q+r and q=fl(a/b) provided no underflow occurs.
%Input a,b may be vectors or matrices as well, in single or double precision.
%
%Follows N. Louvet: Algorithmes compensés en arithm\'etique flottante: pr\'ecision, validation, performances. Ph.D. Dissertation, 2007.

    q=a./b;
    [x,y]=TwoProduct(q,b);
    r=(a-x)-y;
end

