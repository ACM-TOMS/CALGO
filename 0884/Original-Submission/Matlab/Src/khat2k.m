function [zx,zy]= khat2k(x,y,B,b)

% KHAT2K computes the image by an affine transformacion 
%
% [ZX,ZY]= KHATK2(X,Y,B,b) returns in Z1,Z2 the images of 
% (X,Y) by the affine transformation in R^2 given by the 
% 2x2 matrix B and the vector 2x1 b.
%
% X,Y must be vectors of the same length. ZX,ZY are vector of the 
% same length as X and Y.

zx=zeros(size(x));
zy=zeros(size(x));

ch=B*[x(:)'; y(:)']+b *ones(1,length(x));
zx(:)=ch(1,:);
zy(:)=ch(2,:);

return