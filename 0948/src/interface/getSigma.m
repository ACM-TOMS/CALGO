function sigma = getSigma(sadata)
%getSigma returns the signature matrix of a DAE in dense form.
%sigma = getSigma(sadata)
%sadata is an object returned by daeSA.
%
%Example:
%       n = 6; G = 9.8; L = 1.0; c = 0.1;
%       sadata = daeSA(@modified2pendula,n,G,L,c);
%       sigma = getSigma(sadata)
%
%outputs
%
%sigma =
%     2  -Inf     0  -Inf  -Inf  -Inf
%  -Inf     2     0  -Inf  -Inf  -Inf
%     0     0  -Inf  -Inf  -Inf  -Inf
%  -Inf  -Inf  -Inf     2  -Inf     0
%  -Inf  -Inf  -Inf  -Inf     3     0
%  -Inf  -Inf     2     0     0  -Inf
%
% See also daeSA, showStruct
%
% Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
%
sigma = DAESAgetSigma(sadata);
end