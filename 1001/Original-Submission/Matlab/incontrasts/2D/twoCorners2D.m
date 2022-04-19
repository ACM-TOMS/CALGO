%% Two corners in 2D
% For reference see <incontrastsRef.html>.

function q = twoCorners2D(X1,X2,varargin)

seti = varargin{end};
% R = seti.rCD/2;

q = corner2D(X1,X2,seti);

q = q+rot90(q,2);

q = double(q);

end
