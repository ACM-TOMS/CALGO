%% Corner in 2D
% For reference see <incontrastsRef.html>.

function q = corner2D(X1,X2,varargin)

seti = varargin{end};
R = seti.rCD/2;

% p1 = -4/8; % top and left
% l = 9/8; % length
% w = 1/8; % width

% p1 = -5/8;
% l = 7/8;
% w = 1/8;

p1 = -4/8;
l = 5/8;
w = 1/16;

p2 = p1+l; % bottom and right
p3 = p1+w;

q1 = (p1*R <= X1) & (X1 < p2*R) & (p1*R <= X2) & (X2 < p3*R); % top line
q2 = (p1*R <= X1) & (X1 < p3*R) & (p3*R <= X2) & (X2 < p2*R); % left line
q = q1 | q2;

q = double(q);

end