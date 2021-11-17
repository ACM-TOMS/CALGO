function [N,T]=Mesh3d(p)
% MESH3d generates simple tetrahedra meshes 
%   [N,T]=Mesh3d(p); generates an initial coarse tetrahedra mesh for
%   each value of p. The result is a list of tetrahedra T which points
%   into the list of nodes N containing x and y coordinates.  The
%   tetrahedra contains in entries 5 to 8 the neighboring tetrahedra
%   indices, and as a guard # of tetrahedra + 1 if there is no
%   neighbor.

if p==1                    % a cube
  N=[0 1 0 0 1 1 0 1
     0 0 1 0 0 1 1 1
     0 0 0 1 1 0 1 1];
  T=[1 2 3 4 6 5 6 6
     2 6 3 8 6 6 5 6
     2 8 4 5 5 6 6 6
     3 4 8 7 5 6 6 6
     2 3 4 8 1 4 3 2];
elseif p==2                      % just one tetrahedra as default
  N=[0 1 0 0
     0 0 1 0
     0 0 0 1];
  T=[1 2 3 4 2 2 2 2];
end
