%% Cube-like object in 3D
% For reference see <incontrastsRef.html>.
%
% *Update 20180912*
%
% * Problem: The arms overlap each other.
% * Solution: use q = logical(q) as a workaround.
%
function q = cubeLike3D(X1,X2,X3, varargin)
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

% % all arms with length l in x1 (without short notation)
% q1 = arm(R,X1,X2,X3,[p,p+l,p,p+w,p,p+w]); % front bottom
% q2 = arm(R,X1,X2,X3,[p,p+l,p,p+w,p+l-w,p+w+l-w]); % front top (+l-w in x3)
% q3 = arm(R,X1,X2,X3,[p,p+l,p+l-w,p+w+l-w,p,p+w]); % back bottom (+l-w in x2)
% q4 = arm(R,X1,X2,X3,[p,p+l,p+l-w,p+w+l-w,p+l-w,p+w+l-w]); % back top (+l-w in x2 and +l-w in x3)

% all arms with length l in x1
q1 = arm(R,X1,X2,X3,[p,pl,p,pw,p,pw]); % front bottom
q2 = arm(R,X1,X2,X3,[p,pl,p,pw,p+s,pw+s]); % front top (+s in x3)
q3 = arm(R,X1,X2,X3,[p,pl,p+s,pw+s,p,pw]); % back bottom (+s in x2)
q4 = arm(R,X1,X2,X3,[p,pl,p+s,pw+s,p+s,pw+s]); % back top (+s in x2 and +s in x3)

% all arms with length l in x2 (limits for x1 and x2 were interchanged compared to q1 to q4)
q5 = arm(R,X1,X2,X3,[p,pw,p,pl,p,pw]);
q6 = arm(R,X1,X2,X3,[p,pw,p,pl,p+s,pw+s]);
q7 = arm(R,X1,X2,X3,[p+s,pw+s,p,pl,p,pw]);
q8 = arm(R,X1,X2,X3,[p+s,pw+s,p,pl,p+s,pw+s]);

% all arms with length l in x3 (limits for x1 and x3 were interchanged compared to q1 to q4)
q9  = arm(R,X1,X2,X3,[p,pw,p,pw,p,pl]);
q10 = arm(R,X1,X2,X3,[p+s,pw+s,p,pw,p,pl]);
q11 = arm(R,X1,X2,X3,[p,pw,p+s,pw+s,p,pl]);
q12 = arm(R,X1,X2,X3,[p+s,pw+s,p+s,pw+s,p,pl]);

q = q1 + q2 + q3 + q4 + q5 + q6 + q7 + q8 + q9 + q10 + q11 + q12;
q = logical(q);
q = double(q);
end

function res = arm(R,X1,X2,X3,lim)
% lim : vector of size 1 x 6 with limits for X1, X2 and X3
lim = R*lim;
res = ((lim(1) <= X1) & (X1 < lim(2)) & (lim(3) <= X2) & (X2 < lim(4)) & (lim(5) <= X3) & (X3 < lim(6))); 
end
