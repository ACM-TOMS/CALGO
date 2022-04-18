function [Nr,Tr]=RefineMesh3d(N,T);
% REFINEMESH3D refines the mesh by a factor of 8
%   [Nr,Tr]=RefineMesh3d(N,T); refines the mesh given by the nodes N and
%   the tetrahedra T by cutting each tetrahedra into 8 smaller ones. 
%   The neighboring information in the tetrahedra is also updated,
%   and the boundary guard is still the number of tetrahedra plus one,
%   as in NewMesh3d.

Nr=N;                            % new node list starts with old one
nn=size(N,2);
Tr=[];                           % tetrahedra start from scratch
nt=0;                            % number of new tetrahedra
ep=10*eps;                       % to check identical nodes
for j=1:size(T,1),               % divide all old tetrahedra
  i=T(j,1:4); n=N(:,i);          % old nodes of current tetrahedra
  n(:,5)=(n(:,1)+n(:,2))/2;      % 6 new nodes to be used
  n(:,6)=(n(:,1)+n(:,3))/2;
  n(:,7)=(n(:,2)+n(:,3))/2;
  n(:,8)=(n(:,1)+n(:,4))/2;
  n(:,9)=(n(:,2)+n(:,4))/2;
  n(:,10)=(n(:,3)+n(:,4))/2;
  for k=5:10,                    % insert new nodes in Nr if necessary
    [i(k),Nr]=InsertPoint(n(:,k),Nr);
  end;                           % new tetrahedra with known neighborinfo
  Tr(nt+1,:)=[i(1) i(5) i(6)  i(8)  -1  nt+5  -1   -1 ]; % corner tetrahedra 
  Tr(nt+2,:)=[i(5) i(2) i(7)  i(9)  -1   -1  nt+6  -1 ];  
  Tr(nt+3,:)=[i(6) i(7) i(3)  i(10) -1   -1   -1  nt+7];  
  Tr(nt+4,:)=[i(8) i(9) i(10) i(4) nt+8  -1   -1   -1 ];  
  Tr(nt+5,:)=[i(5) i(7) i(6)  i(8)  -1  nt+7 nt+1 nt+6]; % interior tetrahedra 
  Tr(nt+6,:)=[i(5) i(9) i(7)  i(8) nt+2 nt+8 nt+5  -1 ];  
  Tr(nt+7,:)=[i(6) i(7) i(10) i(8) nt+3 nt+8  -1  nt+5];  
  Tr(nt+8,:)=[i(8) i(9) i(7) i(10) nt+6  -1  nt+7 nt+4];  
  nt=nt+8;
end;
nn=[1 2 3 5                      % new neighbor indices
    2 3 4 8
    3 4 1 7
    4 1 2 6];
for j=1:size(T,1),               % updating remaining neighborinfo
  ne=T(j,5:8);
  if ne(1)==size(T,1)+1          % real boundary
    Tr(8*(j-1)+nn(1,:),5)=size(Tr,1)+1;
  else  
    fn=FindFace(T(j,1:3),T(ne(1),1:4));
    Tr(8*(j-1)+nn(1,:),5)=8*(ne(1)-1)+fn;
  end;
  if ne(2)==size(T,1)+1          % real boundary
    Tr(8*(j-1)+nn(2,:),6)=size(Tr,1)+1;
  else    
    fn=FindFace(T(j,2:4),T(ne(2),1:4));
    Tr(8*(j-1)+nn(2,:),6)=8*(ne(2)-1)+fn;
  end;
  if ne(3)==size(T,1)+1          % real boundary
    Tr(8*(j-1)+nn(3,:),7)=size(Tr,1)+1;
  else
    fn=FindFace(T(j,[3 4 1]),T(ne(3),1:4));
    Tr(8*(j-1)+nn(3,:),7)=8*(ne(3)-1)+fn;
  end;
  if ne(4)==size(T,1)+1          % real boundary
    Tr(8*(j-1)+nn(4,:),8)=size(Tr,1)+1;
  else
    fn=FindFace(T(j,[4 1 2]),T(ne(4),1:4));
    Tr(8*(j-1)+nn(4,:),8)=8*(ne(4)-1)+fn;
  end;
end;

function f=FindFace(X,Y)
% FINDFACE find face in tetrahedra 
%   f=FindFace(X,Y) finds face X in thetrahedra Y, where the face and
%   tetrahedra are given by node indices.

f(1)=find(Y==X(1));
f(2)=find(Y==X(2));
f(3)=find(Y==X(3));
switch sum(f)
 case 6                 % face 1 2 3
  f(4)=5;
 case 9                 % face 2 3 4
  f(4)=8;
 case 8                 % face 3 4 1
  f(4)=7; 
 case 7                % face 4 1 2
  f(4)=6;
end; 
f=f';

