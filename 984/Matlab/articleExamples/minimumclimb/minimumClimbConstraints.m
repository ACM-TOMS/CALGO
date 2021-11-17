function c = minimumClimbConstraints(z,probinfo)
% Constraint wrapper for minimum climb problem
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

D = probinfo.LGR.Dmatrix;

Y = z(probinfo.map.state); 
U = z(probinfo.map.control);
YNp1 = z(probinfo.map.stateNp1).';
tf = z(probinfo.map.tf);

% Build X = [Y U]
X = [Y, U];
auxdata = probinfo.auxdata;

% Get Dynamics
F = minimumClimbDynamics(X,auxdata);

% Defect Constraints
C = D*[Y;YNp1] - (tf/2)*F;

% Unroll
c = C(:);
end
