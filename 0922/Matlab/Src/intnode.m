function [Inode,Bnode]=intnode(mesh)
%------------- Find interior and boundary nodes ---------------------------
% find boundary points using boundary edge
Bnode = sort(unique([mesh.e(1,:) mesh.e(2,:)]))';
% find interior points using boundary points
Inode = setdiff((1:length(mesh.p)),Bnode)';
