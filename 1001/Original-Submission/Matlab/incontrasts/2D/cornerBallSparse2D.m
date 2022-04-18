%% Non-constant corner, ball and sparse object in 2D
% For reference see <incontrastsRef.html>.

function q = cornerBallSparse2D(X1,X2,varargin)
% two corners and one ball with non constant contrast

seti = varargin{end};
R = seti.rCD/2;

n = sqrt(length(X1)); % n = nCD or nROI

% -- corner: start
qCorner = corner2D(X1,X2,seti); % size 1 x seti.nROI^2 or 1 x seti.nCD^2
qCorner = 0.8*rot90(qCorner,2); % corner right top (0.8 as in twoCornersOneBall2D.m)

% new layer...
lspace = linspace(0,1,n); % values between 0 and 1
L = repmat(lspace,[n 1]); % Layer: left high, right low...
% L is a matrix... but we need a vector...
L = reshape(L,[1 n^2]);
L = 1.5.*L;

qCorner = L.*qCorner;
%-- corner: end

% -- ball: start
qBall = 1.0*((X2-15/32*R).^2 + (X1+15/32*R).^2 <= (1/16*R)^2);
% -- ball: end

% -- sparse: start
qLeft = corner2D(X1,X2,seti);
y = ZeroOne(linspace(1,18,n)); % first version
%y = ZeroOne(linspace(1,128,n)); % second version (finer)
[~, Y2] = meshgrid(y,y);
qSparse = reshape(Y2.*reshape(qLeft,[n n]),[1 n^2]);
% -- sparse: end

% -- all 3 objects together --
q = qCorner + qBall + qSparse;
q = double(q);
end

function y = ZeroOne(x)
% input  x is a vector
% output y is a vector
f = zeros(1,length(x));
y = zeros(1,length(x));
for i = 1:length(x)
    f(i) = floor(x(i));
    if mod(f(i),2) == 0
        % even
        y(i) = 1;
    else
        % odd
        y(i) = 0;
    end
end

end
