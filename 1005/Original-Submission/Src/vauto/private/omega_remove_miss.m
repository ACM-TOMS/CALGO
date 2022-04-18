%OMEGA_REMOVE_MISS  Remove rows and columns of missing values from Omega
%
%  [Suo, Olowo] = OMEGA_REMOVE_MISS(Su, Olow, miss) removes rows and columns
%  given by the vector miss(:) from Omega, keeping only rows and columns
%  corresponding to observed values. Omega is stored in two parts, Su and Olow
%  as returned by omega_build, and the updated Omega is similarly returned in
%  two parts, a full upper left partition Suo and a sparse block-band lower part
%  Olowo (which will in general have variable sized blocks). The r×n logical
%  matrix miss should be true in missing locations, thus specifiying the
%  structure of Olowo.
%
%  [Suo,Olowo,Suod,Olowod] = OMEGA_REMOVE_MISS(Su,Olow,miss,Sud,Olowd) also
%  updates the derivative of Omega, stored in the three dimensional arrays Sud
%  and Olowd as calculated by omega_build_deriv.

function [Suo, Olowo, varargout] = omega_remove_miss(Su, Olow, miss, Sud, Olowd)
  [r, n] = size(miss);
  ko = find_missing_info(miss);
  DIFF = nargout > 2;
  h = size(Su,1)/r;
  q = size(Olow,2)/r - 1;
  observed = ~miss;
  obs1 = observed(:,1:h);
  obs2 = observed(:,h+1:end);
  Suo = Su(obs1, obs1);  % copy Su-partit. (also OK because of how Matlab works)
  if DIFF
    Suod = Sud(obs1, obs1, :);
    Olowod = zeros(size(Olowd));
  end
  % REMOVE COLUMNS FROM LOWER PART(S):
  ro = diff(ko);
  J = 1 : r;
  K = h+1-q : h+1;
  Olowo = zeros(size(Olow));
  for t = 1:n-h
    obsk = observed(:,K);
    m = sum(ro(K));
    Olowo(J, 1:m) = Olow(J, obsk(:));
    if DIFF, Olowod(J, 1:m, :) = Olowd(J, obsk(:), :); end
    J = J+r;
    K = K+1;
  end
  Olowo = Olowo(obs2, :);  % remove rows
  if DIFF
    Olowod = Olowod(obs2, :, :);
    [varargout{1:2}] = deal(Suod, Olowod);
  end
end
