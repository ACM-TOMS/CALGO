function PlotTetrahedron(X,c)
% PLOTTETRAHEDRON plots a tetrahedron in 3d
%   PlotTetrahedron(X) plots the tetrahedron given by the four
%   points in the matrix X, stored column-wise, with the color c

line(X(1,:),X(2,:),X(3,:),'Color',c)
line(X(1,[4 1]),X(2,[4 1]),X(3,[4 1]),'Color',c)
line(X(1,[1 3]),X(2,[1 3]),X(3,[1 3]),'Color',c)
line(X(1,[2 4]),X(2,[2 4]),X(3,[2 4]),'Color',c)

