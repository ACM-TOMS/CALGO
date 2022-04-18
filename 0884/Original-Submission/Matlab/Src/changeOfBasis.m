function [C,B,b,Th]=changeOfBasis(p)

% CHANGEOFBASIS Change of basis
% 
% [C,B,b,Th]=CHANGEOFBASIS(P) Computes matrices C, B, Th 
% and vector b used in the change of coordinates between the 
% reference triangle and the triangle given by P. 
%
% P is a 2 x 3 matrix with the coordinates of the vertices 
% of the triangle. 

if ~isequal(size(p),[2 3])
    error('The third argument must be a 2x3 matrix')
end

v=[p(:,2)- p(:,1), p(:,3)-p(:,1), p(:,3)-p(:,2)];
[B,b,Th]=afftrans(p);
sides=diag([norm(v(:,1)),norm(v(:,2)),norm(v(:,3))]);
aux=sides^(-2)*[0             1;...
               -1             0;...
             -1/sqrt(2) -1/sqrt(2)]*B';
R=[0 -1; 1 0];
g=dot(aux',R*v);  
f=dot(aux',v);
D=blkdiag(eye(3), B',B',B',Th,Th,Th,[diag(g) diag(f)]);
E=zeros(24,21);
E(1:21,:)=blkdiag(eye(18),sides);
w=[v(1,:).^2; 2*v(1,:).*v(2,:); v(2,:).^2]';
E(22:24,1:3)= 15/8*[-1  1 0;  -1  0 1; 0 -1 1];
E(22:24,4:9)=-7/16 *[v(:,1)'  v(:,1)'  0          0;...
                     v(:,2)'   0       0         v(:,2)';...
                       0       0      v(:,3)'    v(:,3)'];
E(22:24,10:18)=1/32*[-w(1,:)  w(1,:)  0      0          0;...
                     -w(2,:)   0      0      0        w(2,:);...
                       0       0      0   -w(3,:)     w(3,:)];
                   
C=D*E;  
return

function [B,b,Th]=afftrans(p)
% AFFTRANS affine transformation 
%
% [B,b,TH]=AFFTRANS(P) Computes the matrix B 2x2 and the 
% 2 x 1 vector b of the affine transformation mapping the 
% reference triangle onto the triangle with vertices P. 
% TH is a 3 x 3 matrix with the change of variables
% for the second derivatives 
%
% P is a 2 x 3 matrix with the coordinates of the vertices 
% of the triangle. 

B=[p(:,2)-p(:,1), p(:,3)-p(:,1)];
b=p(:,1);
Th=kron(B',B');
Th(3,:)=[];
Th(:,2)=Th(:,2)+Th(:,3); Th(:,3)=[];

return