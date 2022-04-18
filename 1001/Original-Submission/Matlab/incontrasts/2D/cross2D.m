%% Cross in 2D
% For reference see <incontrastsRef.html>.

function q = cross2D(X1,X2,varargin)

seti = varargin{end};
R = seti.rCD/2;

w = 1/8*R; % width
l = 7/8*R; % length
w = w/2;
l = l/2;

q1 = (-w < X1) & (X1 < w) & (-l < X2) & (X2 < l); % vertical line
q2 = (-w < X2) & (X2 < w) & (-l < X1) & (X1 < l); % horizontal line
q = q1 | q2; % to avoid q1+q2 (because some elements would have contrast 2...)

q = double(q);

end
