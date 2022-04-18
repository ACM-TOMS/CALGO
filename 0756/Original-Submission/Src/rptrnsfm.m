function z = rptrnsfm(y,corners)
%RPTRNSFM (not intended for calling directly by the user)
%       Transform optimization vars to prevertices for rectangle
%       parameter problem. 
%
%       Written by Toby Driscoll.  Last updated 5/23/95.

n = length(y)+3;
z = zeros(n,1);
z(corners(1)-1:-1:1) = cumsum(-exp(y(corners(1)-1:-1:1)));
z(corners(1)+1:corners(2)-1) = cumsum(exp(y(corners(1):corners(2)-2)));
z(corners(4)+1:n) = cumsum(-exp(y(corners(4)-2:n-3)));
z(corners(4)-1:-1:corners(3)+1) = cumsum(...
    exp(y(corners(4)-3:-1:corners(3)-1)));
xr = z([corners(2)-1,corners(3)+1]);
z(corners(2)) = mean(xr)+sqrt(diff(xr/2)^2+exp(2*y(corners(2)-1)));
z(corners(3)) = z(corners(2));
z(corners(2)+1:corners(3)-1) = z(corners(2)) + cumsum(...
    exp(y(corners(2):corners(3)-2)));
z(corners(3):n) = i + z(corners(3):n);

