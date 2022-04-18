%% Corner in 3D (tripod)
% For reference see <incontrastsRef.html>.
%
function q = corner3D(X1,X2,X3, varargin)
% 3D angle-shaped object (tripod)

seti = varargin{end};
R = seti.rCD/2;

p1 = -4/8;
l = 7/8;
w = 1/8;

p2 = p1+l; % bottom and right
p3 = p1+w;

% first horizontal arc of angle
q1 = 1.0*((p1*R <= X1) & (X1 < p2*R) & (p1*R <= X2) & (X2 < p3*R) & (p1*R <= X3) & (X3 < p3*R)); 
% second horizontal arc of angle
q2 = 0.8*((p1*R <= X1) & (X1 < p3*R) & ((p1+w)*R <= X2) & (X2 < p2*R) & (p1*R <= X3) & (X3 < p3*R)); 
% vertical arc of angle
q3 = 0.6*((p1*R <= X1) & (X1 < p3*R) & (p1*R <= X2) & (X2 < p3*R) & ((p1+w)*R <= X3) & (X3 < p2*R)); 
q = q1 + q2 + q3; % interesting is the middle of all arcs...

q = double(q);
end

