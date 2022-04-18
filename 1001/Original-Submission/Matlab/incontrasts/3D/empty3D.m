%% Empty contrast in 3D
% For reference see <incontrastsRef.html>.

function q = empty3D(X1,X2,X3,varargin)
% empty contrast, no obstacle (scattering object)

% q = 0; % simple, but then anisotropic case is detected...

q = 0.*X1;
q = double(q);

end

