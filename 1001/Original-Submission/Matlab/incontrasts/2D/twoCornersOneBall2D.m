%% Two corners and on ball in 2D
% For reference see <incontrastsRef.html>.

function q = twoCornersOneBall2D(X1,X2,varargin)

seti = varargin{end};
R = seti.rCD/2;

q = twoCorners2D(X1,X2,seti);

% q = 0.8*q + 1.0*((X2-9/16*R).^2 + (X1+9/16*R).^2 <= (1/16*R)^2);

q = 0.8*q + 1.0*((X2-15/32*R).^2 + (X1+15/32*R).^2 <= (1/16*R)^2);

q = double(q);

end
