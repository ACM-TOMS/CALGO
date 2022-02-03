function B_visualization_spline(P, cc, n, varargin)

% Visualization of a spline
%
% INPUT
%   P     : B-spline patch
%   cc    : vector of coefficients
%   n     : number of evaluation points (optional)
%   specs : pass any number of plot specifications (optional)

if nargin < 3
   n = 100;
end
[a, b] = B_domain(P);
xx = linspace(a, b, n);
ss = B_evaluation_spline(P, cc, xx);
plot(xx, ss, varargin{:});

end
