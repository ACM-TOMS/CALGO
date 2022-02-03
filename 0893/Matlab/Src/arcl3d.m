function [t,ier] = arcl3d(x,y,z)
% arcl3d:  Computes cumulative arc lengths along a space curve
%
% USAGE:  [t,ier] = arcl3d(x,y,z);
%
%   Given an ordered sequence of N points (X,Y,Z) defining a
% polygonal curve in 3-space, this function computes the
% sequence T of cumulative arc lengths along the curve:
% T(1) = 0 and, for 2 <= K <= N, T(K) is the sum of
% Euclidean distances between (X(I-1),Y(I-1),Z(I-1)) and
% (X(I),Y(I),Z(I)) for I = 2 to K.  A closed curve corre-
% sponds to X(1) = X(N), Y(1) = Y(N), and Z(1) = Z(N).  More
% generally, duplicate points are permitted but must not be
% adjacent.  Thus, T contains a strictly increasing sequence
% of values which may be used as parameters for fitting a
% smooth curve to the sequence of points.  
%
% On input:
%
%       X,Y,Z = Vectors of length N containing the coordi-
%               nates of the points.
%
% On output:
%
%       T = Vector of size(X) containing cumulative arc 
%           lengths defined above.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = I if X(I) = X(I+1), Y(I) = Y(I+1), and 
%                     Z(I) = Z(I+1) for some I in the range
%                     1 to N-1.
%
% Modules required by ARCL3D:  None
%
%***********************************************************

% Set ds to the vector of arc lengths, and compute t.

ds = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
m = size(x);
if (m(1) == 1)
   t = [0 cumsum(ds)];
else
   t = [0; cumsum(ds)];
end
if (nargout > 1) 

% Test for a zero arc length.

   if (all(ds))
      ier = 0;
   else
      ier = find(~ds,1);  % ier = index of first zero
   end
end
return;

end  % arcl3d
