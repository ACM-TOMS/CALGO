function [N,T]=Mesh2d(p)
% MESH2D generates simple 2d triangular meshes 
%   [N,T]=Mesh2d(p); generates an initial coarse triangular mesh for
%   each value of p. The result is a list of triangles T which points
%   into the list of nodes N containing x and y coordinates.  The
%   triangle contains in entries 4 to 6 the neighboring triangle
%   indices, and as a guard # of triangles + 1 if there is no
%   neighbor.

if p==1
  N=[0 1 0 1
     0 0 1 1];
  T=[1 2 4 3 3 2
     1 4 3 1 3 3];
elseif p==2
  N=[0 1 0 1 0.5
     0 0 1 1  1];
  T=[1 2 5 4 2 3
     2 4 5 4 4 1
     1 5 3 1 4 4];  
end
