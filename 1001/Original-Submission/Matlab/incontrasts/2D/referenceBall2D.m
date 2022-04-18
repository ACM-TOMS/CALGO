%% Ball in 2D
% For reference see <incontrastsRef.html>.

function q = referenceBall2D(X1,X2, varargin)
% seti.qBall: contrast value q of the ball
% seti.rBall: radius

seti = varargin{end};

seti = checkfield(seti,'qBall',0.8);
seti = checkfield(seti,'rBall',0.015);

qBall = seti.qBall;
rBall = seti.rBall;

r = hypot(X1, X2);
q = qBall*( r < rBall );
q = double(q);

end
