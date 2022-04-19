function [k,P]=InsertPoint(p,P);
% INSERTPOINT inserts a new intersection point
%   [k,P]=InsertPoint(p,P); inserts an intersection point p into
%   the list of intersection points P. If the point p is already
%   contained in P, the point is not inserted, and k is the index
%   of the point found in the list. Otherwise, the new point is
%   inserted, and k is its index in the list.

ep=10*eps;                         % tolerance for identical nodes
k=1;
while k<=size(P,2) & norm(p-P(:,k))>ep
  k=k+1;
end;
if k>size(P,2)                     % new node is not contained in P
  P(:,k)=p;
end;
