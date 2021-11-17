%% Cross in 3D
% For reference see <incontrastsRef.html>.

function q = cross3D(X1,X2,X3, varargin)
% 3D cross-shaped object

seti = varargin{end};
R = seti.rCD/2;

p0 = -3/8;
p1 = -1/8;
l = 7/8; % length of bars
w = 1/8; % width of bars

% first horizontal bar of cross
q1 = (p0*R <= X1) & (X1 < (p0+l)*R) & (p1*R <= X2) & (X2 < (p1+w)*R) & (p1*R <= X3) & (X3 < (p1+w)*R);
% second horizontal bar of cross
q2 = (p1*R <= X1) & (X1 < (p1+w)*R) & (p0*R <= X2) & (X2 < (p0+l)*R) & (p1*R <= X3) & (X3 < (p1+w)*R);
% vertical bar of cross
q3 = (p1*R <= X1) & (X1 < (p1+w)*R) & (p1*R <= X2) & (X2 < (p1+w)*R) & (p0*R <= X3) & (X3 < (p0+l)*R);
% center block
q4 = (p1*R <= X1) & (X1 < (p1+w)*R) & (p1*R <= X2) & (X2 < (p1+w)*R) & (p1*R <= X3) & (X3 < (p1+w)*R);

q = 0.8*q1 + 0.6*q2 + 1.0*q3 - (0.8+0.6)*q4; % without center block...

q = double(q);
end

