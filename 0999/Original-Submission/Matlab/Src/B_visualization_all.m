function B_visualization_all(P, n, varargin)

% Visualization of all B-splines
%
% INPUT
%   P     : B-spline patch
%   n     : number of evaluation points (optional)
%   specs : pass any number of plot specifications (optional)

if nargin < 2
   n = 100;
end
[a, b] = B_domain(P);
xx = linspace(a, b, n);
M = B_evaluation_all(P, xx);
plot(xx, M, varargin{:});

end
