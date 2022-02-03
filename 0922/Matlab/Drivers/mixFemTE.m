clear all
close all
% load meshfile 
mfile = input('Input name of the meshfile --> ');
load(mfile);

n=input('Input the supremum of index of refraction --> ');

noe=input('Input the number of transmission eigenvalues --> ');

% constructing the mass and stiffness matrix, also find the boundary nodes.
disp(['Construct the stiffness and mass matrices ... '])
[S, M, Mn] = assemble(mesh);
disp(['find interior and boundary nodes ... '])
[Inode,Bnode]=intnode(mesh);

% [S, M, Inode, Bnode] = constructurematrix(mesh);
% compute the first eigenvalue for -\Delta
disp(['Compute the first Dirichlet eigenvalue ... '])
lambda=DirichletEig(S(Inode,Inode), M(Inode,Inode));

disp(['Compute the transmission eigenvalues ... '])
% constructure the matrix for A and B.
[A,B]=MixMethod(S, M, Mn, Inode, Bnode);

% compute the first noe transmission eigenvalue.
k=sptarnite(A,B,lambda/n,noe)
