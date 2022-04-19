% IS_STATIONARY  Check VARMA model for stationarity
%   IS_STATIONARY(A, Sig) returns true if the VARMA-process x(t) = A1·x(t-1) +
%   ... + Ap·x(t-p) + y(t) is stationary, false otherwise. A is the r × r·p
%   matrix [A1...Ap]. The time series y(t) is given by a moving average process,
%   y(t) = B1·y(t-1) + ... + Bq·y(t-q) + eps(t), where eps(t) is N(0,Sig).
%   Whether x is stationary or not does not depend on B, so it need not be a
%   parameter.
%
%   IS_STATIONARY(A, Sig, PLU) uses PLU from a previous call to vyw_factorize.
%
%   The test is based on solving the modified Yule-Walker equations for the AR
%   process x(t)=A1·x(t-1)+...+Ap·x(t-p)+eps(t) with eps(t) ~ N(0,I), but the
%   solution to these equations gives a postive definite Su matrix (see
%   omega_build) if and only if all the roots of the polynomial fi(b) = I - A1·b
%   - A2·b^2 - ... - Ap·B^p are outside the unit circle, which characterizes a
%   stationary VARMA process.

function stat = is_stationary(A, Sig, PLU)
  [R,pp] = chol(Sig);
  if pp>0, stat=false; return, end;
  if isempty(A), stat=true; return, end
  r = size(A,1);
  p = length(A)/r;
  if nargin==2
    PLU = vyw_factorize(A);
    if ~isempty(PLU) && ~isempty(PLU{1}) && PLU{1}(1) == 0
      stat=false; 
      return
    end
  end
  I = eye(r);
  S = vyw_solve(A, PLU, {I});
  Su = omega_build(S, {I}, {I}, p, p);
  [R, pp] = chol(Su');
  stat = pp==0;
end
