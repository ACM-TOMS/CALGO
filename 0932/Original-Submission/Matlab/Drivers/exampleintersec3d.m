% EXAMPLEINTERSEC3D test the intersection program for three dimensional meshes

X=[0   1   1/2  1/2
   0   0   1    1/2
   0   0   0    1]
Y=[0   1   1/2  1/2
   1   1  -1/2  1/2
   1/2 1/2 1/2 -1/2]

[P,nc,H,v]=Intersection3d(X,Y);
