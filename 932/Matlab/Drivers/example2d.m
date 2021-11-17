% 2D EXAMPLE: test the advancing front program for two dimensional meshes

[N1,T1]=Mesh2d(1); [N1,T1]=RefineMesh(N1,T1); [N1,T1]=RefineMesh(N1,T1);
[N2,T2]=Mesh2d(2); [N2,T2]=RefineMesh(N2,T2); [N2,T2]=RefineMesh(N2,T2);

fig=figure(1); clf; set(fig,'DoubleBuffer','on');
PlotMesh(N1,T1,'b'); PlotMesh(N2,T2,'r');
xlabel('x'); ylabel('y'); title('nonconforming grids in 2d');

M=InterfaceMatrix(N1,T1,N2,T2);