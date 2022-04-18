%% Triangle in 2D
% For reference see <incontrastsRef.html>.

function q = triangle2D(X1,X2,varargin)

seti = varargin{end};
R = seti.rCD/2;

% three lines
f1 = @(t) -3*t-R;
f2 = @(t) -1/2*t+R/4;
f3 = @(t) t-1/2*R;

q1 = (f1(X1) <= X2);
q2 = (f2(X1) > X2);
q3 = (f3(X1) <= X2);

q = q1 & q2 & q3;

q = 0.8*q;

q = double(q);

end
