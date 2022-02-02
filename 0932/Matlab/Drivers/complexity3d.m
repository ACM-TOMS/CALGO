% 3D EXAMPLE3D compares in computing time for three dimensional projections

nr=3;                           % number of mesh refinements
np=5;                           % number of projection calculations

vtb=zeros(1,nr); vtn=zeros(1,nr); s=zeros(1,nr);
for k=1:nr
  [N1,T1]=Mesh3d(1); [N2,T2]=Mesh3d(1); N2=(N2/1.5);
  for l=1:k
    [N1,T1]=RefineMesh3d(N1,T1); [N2,T2]=RefineMesh3d(N2,T2);
  end
  tb=0;tn=0;
  for i=1:np
    N2(:,1:end)=N2(:,1:end)+(0.02*rand(size(N2(:,1:end))));
    t0=clock; InterfaceMatrix3d(N1,T1,N2,T2); tn=tn+etime(clock,t0);
    t0=clock; InterfaceMatrix3dBruteForce(N1,T1,N2,T2); tb=tb+etime(clock,t0);
  end;
  vtn(k)=tn; vtb(k)=tb; s(k)=size(T1,1);
end;
c1=vtb(end)/s(end)^2; c2=vtn(end)/s(end); % align theoretical and numerical curves
loglog(s,vtb,'--o',s,vtn,'--*',s,c1*s.^2,'-o',s,c2*s,'-*');
legend('Brute Force','New Method','O(n^2)','O(n)',2);
xlabel('number of triangles'); ylabel('time for 20 projections');
