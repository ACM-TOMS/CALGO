function b=NewFace(H,id)
% NEWFACE checks if index set is contained already in H
%   b=NewFace(H,id) checks if a permutation of the vector id is
%   contained in a row of the matrix H.

n=size(H,1);
m=length(id);
A=zeros(n,m);
A(1:n,1:size(H,2))=H;
id=sort(id);
i=0;
b=1==1;                    % true
while b & i<n
  i=i+1;
  b=norm(sort(A(i,1:m))-id)>0;
end;
