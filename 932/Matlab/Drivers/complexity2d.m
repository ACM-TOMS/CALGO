% 2D COMPLEXITY compares the computing time for 2d projections

nr=3;                           % number of mesh refinements
np=20;                          % number of projection calculations
vtb=zeros(1,nr); vtn=zeros(1,nr); s=zeros(1,nr);
for k=1:nr
  [N1,T1]=Mesh2d(1); [N2,T2]=Mesh2d(2);
  for l=1:k
    [N1,T1]=RefineMesh(N1,T1); [N2,T2]=RefineMesh(N2,T2);
  end
  tb=0;tn=0;
  for i=1:np
    N1(:,5:end)=N1(:,5:end)+(0.1*rand(size(N1(:,5:end)))-0.05);
    N2(:,5:end)=N2(:,5:end)+(0.1*rand(size(N2(:,5:end)))-0.05);  
    t0=clock; M=InterfaceMatrix(N1,T1,N2,T2); tn=tn+etime(clock,t0);
    t0=clock; Mb=InterfaceMatrixBruteForce(N1,T1,N2,T2); tb=tb+etime(clock,t0);
  end
  vtn(k)=tn; vtb(k)=tb; s(k)=size(T1,1);
end
c1=vtb(end)/s(end)^2; c2=vtn(end)/s(end); % align theoretical and numerical curves
loglog(s,vtb,'--o',s,vtn,'--*',s,c1*s.^2,'-o',s,c2*s,'-*');
legend('Brute Force','New Method','O(n^2)','O(n)',2);
xlabel('number of triangles'); ylabel('time for 20 projections');
