%% Rectangle in 2D
% For reference see <incontrastsRef.html>.

function q = rectangle2D(X1,X2,varargin)

qValue = 0.5;

w = 0.03; % width
h = 0.06; % height

right = 0.01; % move it right
up    = -0.02; % move it up (or down with a negative number)

q = (-w/2+right <= X1) & (X1 <= w/2+right) & (-h/2+up <= X2) & (X2 <= h/2+up);
q = qValue*q;
q = double(q);

end
