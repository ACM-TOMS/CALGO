function i = intrvl(t,x)
% intrvl:  Finds interval containing a point
%
% USAGE:  i = intrvl(t,x);
%
%   This function returns the index of the left end of an
% interval (defined by an increasing sequence X) that
% contains the value T.  The method consists of first test-
% ing the interval returned by a previous call, if any, and
% then using a binary search if necessary.
%
% On input:
%
%       T = Point to be located.
%
%       X = Vector of length N assumed (without a test) to
%           contain a strictly increasing sequence of
%           values.  N >= 2.
%
% On output:
%
%       I = Index defined as follows:
%
%           I = 1    if T < X(2) or N <= 2,
%           I = N-1  if  T >= X(N-1), and
%           X(I) <= T < X(I+1) otherwise.
%
% Modules required by INTRVL:  None
%
%***********************************************************

persistent il

n = length(x);

if (isempty(il)), il = 1; end
if (il >= 1  &&  il < n)
   if (x(il) <= t  &&  t < x(il+1)) 
      i = il;
      return;
   end
end

% Initialize low and high indexes.

il = 1;
ih = n;

% Binary search:

while (true)
   if (ih <= il+1)
      i = il;
      return;
   end
   k = (il+ih)/2.0;
   if (t < x(k)) 
      ih = k;
   else
      il = k;
   end
end  % while

end  % intrvl
