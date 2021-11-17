%% Empty contrast in 2D
% For reference see <incontrastsRef.html>.

function q = empty2D(X1,X2,varargin)
% empty contrast, no obstacle (scattering object)

% q = 0; % simple, but then anisotropic case is detected...

q = 0.*X1;
q = double(q);

end
