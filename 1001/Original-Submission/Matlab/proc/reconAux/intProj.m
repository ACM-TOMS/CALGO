%% intProj
% Projects each value of a vector or matrix inside a interval.
%
%% Syntax
%
%   res = intProj(a,b,t)
%
%% Description
% |res = intProj(a,b,t)| projects each value of a vector or matrix |t| 
% into the interval $[a,b]$.
%
%% Example
%
% *Code*
%
%   t = [-3.1, -2, 1.5, 5];
%   res = intProj(-2,4,t);
%   A = [t; res]
%
% *Result*
%
%   A = 
%
%   -3.1000   -2.0000    1.5000    5.0000
%   -2.0000   -2.0000    1.5000    4.0000
%
% * First row: original values.
% * Second row: projected values into interval [-2,4].
%
%% Input Arguments
%
% * a and b     :   real numbers to define the interval [a,b]
% * t           :   vector or matrix (of any size)
%
%% Output Arguments
% * res         : result of the projection of each value of t into [a,b].
%
%% More About
%
% This function is used in pda.m.
%
%% See Also
%
% * <pda.html>
%
function res = intProj(a,b,t) % intervall projection
%SS_{a,b}(t) = a if t_i \leq a; t_i if t_i \in (a,b); b if t_i \geq b
%res = @(a,b,t) a.*(t<=a) + t.*and(a<t,t<b) + b.*(b<=t);
res = a.*(t<=a) + t.*and(a<t,t<b) + b.*(b<=t);
%SSabAlternative = @(t) t./max(1,abs(t)); % only if a = b = 1
end

