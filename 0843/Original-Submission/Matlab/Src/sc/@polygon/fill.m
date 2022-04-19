function h = fill(p,varargin)

% Plot a polygon in filled style (see FILL).

v = vertex(p);
h = fill(real(v),imag(v),varargin{:});
