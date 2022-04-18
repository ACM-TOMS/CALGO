function v=PyramidVolume(P,c)
% PYRAMIDVOLUME computes the volume of a pyramid
%   v=PyramidVolume(P,c) computes for the pyramid in 3d with base
%   polygone given by the nodes P stored columnwise and the top c
%   its volume v.

a=0;                           % area
for i=3:size(P,2)
  a=a+norm(cross(P(:,i)-P(:,1),P(:,i-1)-P(:,1)))/2;
end;
no=cross(P(:,3)-P(:,1),P(:,2)-P(:,1));
no=no/norm(no);
h=(c-P(:,1))'*no;
v=abs(a*h/3);
