function MDB_visualization_spline(MP, H, cc, n, varargin)

% Visualization of a multi-degree spline
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%   cc    : vector of coefficients
%   n     : number of evaluation points (optional)
%   specs : pass any number of plot specifications (optional)

if nargin < 4
   n = 100;
end
[a, b] = MDB_domain(MP);
xx = linspace(a, b, n);
ss = MDB_evaluation_spline(MP, H, cc, xx);
plot(xx, ss, varargin{:});

end
