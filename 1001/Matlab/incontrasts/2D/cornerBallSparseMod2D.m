%% Non-constant corner, ball and sparse object in 2D
% For reference see <incontrastsRef.html>.

function q = cornerBallSparseMod2D(X1,X2,varargin)
% rotate 90 degree of: two corners and one ball with non constant contrast

seti = varargin{end};
q = cornerBallSparse2D(X1,X2,seti);
q = rot90(q,2); % rotates 2*90 degrees counterclockwise
q = double(q);
end

