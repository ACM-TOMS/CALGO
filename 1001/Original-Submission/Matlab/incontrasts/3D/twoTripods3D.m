%% Two tripods in 3D
% For reference see <incontrastsRef.html>.
%
function q = twoTripods3D(X1,X2,X3, varargin)
% 3D cube-like object

seti = varargin{end};
R = seti.rCD/2;

p = -4/8; % start position
l = 7/8; % length
w = 1/8; % width

pl = p+l;
pw = p+w;
s = l-w; % shift
% ps = p+s;

%%
% *First Tripod*
q1 = arm(R,X1,X2,X3,[p,pl,p,pw,p,pw]); % front bottom
q2 = arm(R,X1,X2,X3,[p,pw,pw,pl,p,pw]);
q3 = arm(R,X1,X2,X3,[p,pw,p,pw,pw,pl]);

% *Second Tripod*
q4 = arm(R,X1,X2,X3,[p,pl,p+s,pw+s,p+s,pw+s]); % back top (+s in x2 and +s in x3)
q5 = arm(R,X1,X2,X3,[p+s,pw+s,pw,pl,p+s,pw+s]);
q6 = arm(R,X1,X2,X3,[p+s,pw+s,p+s,pw+s,p,pl-w]);

q = q1 + q2 + q3 + q4 + q5 + q6;
q = double(q);
end

function res = arm(R,X1,X2,X3,lim)
% lim : vector of size 1 x 6 with limits for X1, X2 and X3
lim = R*lim;
res = ((lim(1) <= X1) & (X1 < lim(2)) & (lim(3) <= X2) & (X2 < lim(4)) & (lim(5) <= X3) & (X3 < lim(6))); 
end
