% MAKEMISSING  Create missing value test cases
%
%   X = MAKEMISSING(X,i) creates a testcase for missing value functions by
%   assigning a subset of X the value nan. The chosen subset is determined by i
%   as follows:
%                       i=-1  X(1,1) only is missing
%                         -2  X(end,end) only is missing
%                         -3  first four values of X(:) missing
%                         -4  last four values of X(:) missing
%                         -5 or 'miss-5a' 5% missing, all in first quarter
%                         -6 or 'miss-5b' same as i=5
%                         -7 or 'miss-25' 25% missing, single runs at start
%                         >=0 i% of the values (chosen at random) missing (i<90)
%
%   (with r even, the miss-25 option makes half of the series have the first
%   half of the values missing, with r odd (r+1)/2 series have a little less
%   than half of the values at the beginning missing)

function X = makemissing(X, i)
  % Drop some values from X (make them missing by replacing with nan), with a
  % pattern according to i.
  switch i
    case 0, % nothing missing
    case -1, X(1) = nan;
    case -2, X(end) = nan;
    case -3, X(1:4) = nan;
    case -4, X(end-3:end) = nan;
    case {-5,'miss-5a'}
      L = length(X(:));
      k = ceil(5*L/100);
      X(randset(0.25*L, k)) = nan;
    case {-6,'miss-5b'}
      X = makemissing(X,5);
    case {-7,'miss-25'}
      [r,n] = size(X);
      R = ceil(r/2);
      k = round(n*0.25*r/R);
      X(1:R, 1:k) = nan;
    otherwise
      ascertain(i<=90)
      L = length(X(:));
      k = ceil(i*L/100);
      X(randset(L, k)) = nan;
  end
end

function S = randset(n, k)
  % S set to a k-element random set of numbers in {1,2,...,n}
  S = unique(floor(n*rand(2*k,1))+1);
  while length(S) < k, S = union(S, floor(n*rand(k,1))+1); end
  J = randperm(length(S(:)));
  S = S(J(1:k));
end
