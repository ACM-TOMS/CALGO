function [Smat, Mmat, Mnmat] = assemble(mesh)
%
% Construct stifness and mass matrices
%
[weight,place]=quad_setup2;  % get 3 point quad rule on the triangle
nq=length(weight);

yloc=zeros(3,nq);
gyloc=zeros(2,3,nq);

for r=1:nq
  [yloc(:,r),gyloc(:,:,r)]=phihat(place(:,r));
end
[dum,nt]=size(mesh.t);
[dum,nv]=size(mesh.p);


ii=zeros(9*nt,1);
jj=zeros(9*nt,1);
zzM=zeros(9*nt,1);
zzS=zeros(9*nt,1);
zzMn = zeros(9*nt,1);
for it=1:nt
    % The coordinates of the vertices of triangle 'it'
    v=mesh.p(:,mesh.t(1:3,it));
    % The index number for triangle 'it'
    index=mesh.t(4,it);
    
    B=[v(:,2)-v(:,1),v(:,3)-v(:,1)];
    detB=abs(det(B));
    
    b=v(:,1);
    n=zeros(nq,1);
	[L,U,P]=lu(B');
	for r=1:nq
	    gphi(:,:,r)=U\(L\(P*gyloc(:,:,r)));
	end
    for j=1:nq
      x=B*place(:,j)+b;
      n(j)=Rindex(x,index);
    end

    % Construct local matrices
    % Stiffness matrix
    Mloc=zeros(3,3);
    for r=1:nq
        Mloc=Mloc+((gphi(:,:,r)'*gphi(:,:,r)))*weight(r);
    end
    Mloc=Mloc*detB/2;
    indices=[mesh.t(1:3,it)'];
    is=9*(it-1)+(1:3);
    for j=1:3
        ii(is+(j-1)*3)=indices;
        jj(is+(j-1)*3)=indices(j);
        zzS(is+(j-1)*3)=Mloc(j,:);
    end
    % Mass matrix
    Mloc=zeros(3,3);
    for r=1:nq
        Mloc=Mloc+((yloc(:,r)*yloc(:,r)'))*weight(r);
    end
    Mloc=Mloc*detB/2;
    for j=1:3
        zzM(is+(j-1)*3)=Mloc(j,:);
    end
    % Weighted maass matrix
    Mloc=zeros(3,3);
    for r=1:nq
        Mloc=Mloc+((yloc(:,r)*yloc(:,r)')*n(r))*weight(r);
    end
    Mloc=Mloc*detB/2;
    for j=1:3
        zzMn(is+(j-1)*3)=Mloc(j,:);
    end
end
Mmat=sparse(ii,jj,zzM,nv,nv);
Smat=sparse(ii,jj,zzS,nv,nv);
Mnmat=sparse(ii,jj,zzMn,nv,nv);

function [weight,place]=quad_setup2()
%
% 3 point quadrature on the reference element(2nd order exactly)
% 
w1=1/3; 
weight=[w1,w1,w1];
place(:,1)=[0;1/2];
place(:,2)=[1/2;0];
place(:,3)=[1/2;1/2];

function [y,grady]=phihat(xhat)
%
% computes all basis function values and all gradients at the point
% xhat on the reference element
%
y=zeros(3,1);
grady=zeros(2,3);

y(1) = 1 - xhat(1) - xhat(2);
y(2) = xhat(1);
y(3) = xhat(2);

grady(:,1) = [-1; -1];
grady(:,2) = [1; 0];
grady(:,3) = [0; 1];
