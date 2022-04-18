function [c, d, cl, dl] = getOffsets(sadata)
%getOffsets returns the offsets of a DAE.
%[c,d] = getOffsets(sadata)
%[c,d,cl,dl] = getOffsets(sadata)
%sadata is an object returned by daeSA.
%
%If the DAE is structurally well-posed (SWP), c and d are row vectors 
%containing global equation and variable offsets, respectively; and
%cl and dl are row vectors containing local equation and variable offsets,
%respectively. If the DAE is not SWP, all output vectors are [].
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       [c,d,cl,dl] = getOffsets(sadata)
%
%produces
%
%c =
%     4     4     6     0     0     2
%d =
%     6     6     4     2     3     0
%cl =
%     0     0     2     0     0     0
%dl =
%     2     2     0     0     3     0
%
%See also daeSA, isSWP, showStruct.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce

if ~any(nargout==[2 4])
    error('Number of output arguments must be 2 or 4.');
end
[c, d, cl, dl] = DAESAgetOffsets(sadata);
end
