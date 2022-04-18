%% Ball in 2D
% For reference see <incontrastsRef.html>.

function q = referenceBallSmooth2D(X1,X2, varargin)
% seti.qBall: contrast value q of the ball
% seti.rBall: radius

seti = varargin{end};

seti = checkfield(seti,'qBall',100);
seti = checkfield(seti,'rBall',0.05);

qBall = seti.qBall;
rBall = seti.rBall;

r = sqrt(X1.^2+X2.^2);
q = rBall-r;
q(q<=0) = 0;
q(q>0) = qBall*q(q>0).^2;
q = double(q);

end
