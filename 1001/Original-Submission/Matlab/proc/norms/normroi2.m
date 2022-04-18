%% normroi2
% Computes the norm of a vector on the region of interest as normroi but
% always with q = 2.
%% Syntax
%
%   n = normroi2(x,seti)
%
%% See Also
%
% * <normroi.html>
%
%% Code
function n = normroi2(x,seti)
% As normroi but always with Q = 2 (!) (and not from Q = seti.qNorm).
% Is e.g. used in testDerivative.
% It is also used in stepsize choice by Barzilai-Borwein (minShrink.m) --
% not available in public version.

n = sqrt(abs(innerroi(x,x,seti)));
end
