function [meqn, mvar] = getMissing(sadata)
%[meqn, mvar] = getMissing(sadata)
%returns the indices of equations (in vector meqn) and indices of variables
%(in vector mvar) that are missing in the DAE definition. 
%If no equations is missing, meqn=[], and if no variable is missing, mvar=[].
%sadata is an object returned by daeSA.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@illPosed1,n,G,L,c);
%       [meqn,mvar] = getMissing(sadata);
%
%produces
%
%       meqn =
%            3
%       mvar =
%            []
% 
%That is, equation 3 in illPosed1.m is not evaluated.
%
%See also daeSA, isSWP, showStruct.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

[meqn, mvar] = DAESAgetMissing(sadata);

end