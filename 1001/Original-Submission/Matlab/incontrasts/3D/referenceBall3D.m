%% Ball in 3D
% For reference see <incontrastsRef.html>.

function q = referenceBall3D(X1,X2,X3, varargin)
% seti.qBall: contrast value q of the ball
% seti.rBall: radius

seti = varargin{end};

seti = checkfield(seti,'qBall',0.8);
seti = checkfield(seti,'rBall',0.015);

qBall = seti.qBall;
rBall = seti.rBall;

r = sqrt(X1.^2+X2.^2+X3.^2);
q = qBall*( r <= rBall );

q = double(q);

end
