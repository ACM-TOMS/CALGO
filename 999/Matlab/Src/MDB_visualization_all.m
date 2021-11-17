function MDB_visualization_all(MP, H, n, varargin)

% Visualization of all MDB-splines
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%   n     : number of evaluation points (optional)
%   specs : pass any number of plot specifications (optional)

if nargin < 3
   n = 100;
end
[a, b] = MDB_domain(MP);
xx = linspace(a, b, n);
M = MDB_evaluation_all(MP, H, xx);
plot(xx, M, varargin{:});

end
