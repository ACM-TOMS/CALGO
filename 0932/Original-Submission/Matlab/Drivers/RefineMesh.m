function [Nr,Tr]=RefineMesh(N,T);
% REFINEMESH refines the mesh by a factor of 4
%   [Nr,Tr]=RefineMesh(N,T); refines the mesh given by the nodes N and
%   the triangles T by cutting each triangle into four smaller ones. 
%   The neighboring information in the triangles is also updated,
%   and the boundary guard is still number of trianlges plus one,
%   as in NewMesh.

Nr=N;                            % new node list starts with old one
nn=size(N,2);
Tr=[];                           % triangles start from scratch
nt=0;                            % number of new triangles
ep=10*eps;                       % to check identical nodes
for j=1:size(T,1),               % quadrisect all old triangles
  i=T(j,1:3); n=N(:,i);          % old nodes of current triangle
  n(:,4)=(n(:,1)+n(:,2))/2;      % 3 new nodes to be used
  n(:,5)=(n(:,1)+n(:,3))/2;
  n(:,6)=(n(:,2)+n(:,3))/2;
  for k=4:6,                     % insert new nodes in Nr if necessary
    l=find(Nr(1,:)==n(1,k));
    m=find(Nr(2,l)==n(2,k));
    if isempty(m),               % node not found, hence insert
      nn=nn+1;
      Nr(:,nn)=n(:,k);
      i(k)=nn;
    else
      i(k)=l(m);                 % node was found
    end;
  end; 
  Tr(nt+1,:)=[i(1) i(4) i(5) -1 nt+2 -1]; 
  Tr(nt+2,:)=[i(5) i(4) i(6) nt+1 nt+3 nt+4];  
  Tr(nt+3,:)=[i(6) i(4) i(2) nt+2 -1 -1];  
  Tr(nt+4,:)=[i(6) i(3) i(5) -1 -1 nt+2];  
  nt=nt+4;
end;
for j=1:size(T,1),               % updating neighboring information
  i=T(j,1:3);          
  ne=T(j,4:6);
  if ne(1)==size(T,1)+1          % real boundary
    Tr(4*(j-1)+1,4)=size(Tr,1)+1; Tr(4*(j-1)+3,5)=size(Tr,1)+1;
  else  
    id=find(T(ne(1),1:3)==i(2));
    if id==1
      Tr(4*(j-1)+1,4)=4*(ne(1)-1)+3; Tr(4*(j-1)+3,5)=4*(ne(1)-1)+1;
    elseif id==2
      Tr(4*(j-1)+1,4)=4*(ne(1)-1)+4; Tr(4*(j-1)+3,5)=4*(ne(1)-1)+3;
    else
      Tr(4*(j-1)+1,4)=4*(ne(1)-1)+1; Tr(4*(j-1)+3,5)=4*(ne(1)-1)+4;
    end;
  end;
  if ne(2)==size(T,1)+1          % real boundary
    Tr(4*(j-1)+3,6)=size(Tr,1)+1; Tr(4*(j-1)+4,4)=size(Tr,1)+1;
  else    
    id=find(T(ne(2),1:3)==i(3));
    if id==1
      Tr(4*(j-1)+3,6)=4*(ne(2)-1)+3; Tr(4*(j-1)+4,4)=4*(ne(2)-1)+1;
    elseif id==2
      Tr(4*(j-1)+3,6)=4*(ne(2)-1)+4; Tr(4*(j-1)+4,4)=4*(ne(2)-1)+3;
    else
      Tr(4*(j-1)+3,6)=4*(ne(2)-1)+1; Tr(4*(j-1)+4,4)=4*(ne(2)-1)+4;
    end;
  end;
  if ne(3)==size(T,1)+1          % real boundary
    Tr(4*(j-1)+4,5)=size(Tr,1)+1; Tr(4*(j-1)+1,6)=size(Tr,1)+1;
  else
    id=find(T(ne(3),1:3)==i(1));
    if id==1
      Tr(4*(j-1)+4,5)=4*(ne(3)-1)+3; Tr(4*(j-1)+1,6)=4*(ne(3)-1)+1;
    elseif id==2
      Tr(4*(j-1)+4,5)=4*(ne(3)-1)+4; Tr(4*(j-1)+1,6)=4*(ne(3)-1)+3;
    else
      Tr(4*(j-1)+4,5)=4*(ne(3)-1)+1; Tr(4*(j-1)+1,6)=4*(ne(3)-1)+4;
    end;
  end;
end;
