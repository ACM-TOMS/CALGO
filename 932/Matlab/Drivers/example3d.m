% 3D EXAMPLE test the advancing front program for three dimensional meshes

[N1,T1]=Mesh3d(1); [N2,T2]=Mesh3d(1); N2=(N2/1.5);

fig=figure(1); clf; set(fig,'DoubleBuffer','on');
PlotMesh3d(N1,T1,'b'); PlotMesh3d(N2,T2,'r');
xlabel('x'); ylabel('y'); zlabel('z'); title('nonconforming grids in 3d');

InterfaceMatrix3d(N1,T1,N2,T2);

