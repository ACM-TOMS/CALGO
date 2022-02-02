function [A,B]=MixMethod(S, Ma, Mn, Inode, Bnode)
% construct the matrix for the transmission eigenvalue problem
[sizeI dum]=size(Inode); [sizeB dum]=size(Bnode); sizeT=sizeI+sizeB;
A = sparse(sizeI+sizeT, sizeI+sizeT); B = sparse(sizeI+sizeT, sizeI+sizeT);
A(1:sizeI, 1:sizeT) = S(Inode,:); 
A(sizeI+1:sizeI+sizeT, sizeT+1:sizeI+sizeT) = S(:,Inode);
A(sizeI+1:sizeI+sizeT,1:sizeT) = Mn-Ma;
B(1:sizeI, 1:sizeT) = Ma(Inode,:); 
B(sizeI+1:sizeI+sizeT, sizeT+1:sizeI+sizeT) = Mn(:,Inode);


